/**
 * @file dist_array_c.c
 * @brief Distributed data structure of multiple global arrays
 *
 * @copyright Copyright  (C)  2014  Thomas Jahns <jahns@dkrz.de>
 *                                  Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_extents.h"
#include "core/ppm_extents_mp.h"
#include "core/ppm_xfuncs.h"

#include "ppm/dist_array.h"

enum {
  exposed_flag = 1,
  transfer_mode_mask = 1 << 1,
};

struct PPM_global_array_desc_
{
  struct PPM_extent rect[PPM_dma_max_rank];
  unsigned a_rank;
};

struct PPM_dm_array_cache_entry
{
  void *base;
  MPI_Aint *offset;
  unsigned long access_stamp;
  MPI_Datatype composite_dt;
  int rank, composite_count;
};

struct PPM_dist_mult_array
{
  unsigned num_sub_arrays;
  unsigned max_a_rank;
  MPI_Comm comm, comm_self_dup;
  int comm_size, comm_rank;
  /* points to array of size [num_sub_arrays] holding extent of
   * each element datatype in char */
  MPI_Aint *dt_extents;
  /* points to array of size [num_sub_arrays] where dt_is_basic[i] iff
   * element_dt[i] is one of the pre-defined data types */
  bool *dt_is_basic;
  MPI_Win win;
  int flags;
  enum PPM_dma_sync_mode sync_mode;
  unsigned long access_stamp, valid_stamp;
  MPI_Aint max_win_size;
  struct PPM_global_array_desc_ *sub_arrays_global_desc;
  MPI_Datatype *element_dt;
  /* size of cache (always at least 1 (because of local data in cache[0]) */
  size_t cache_size;
  /* points to array of cache entries for remote data (cache[0]
   * always refers to local data) */
  struct PPM_dm_array_cache_entry *cache;
  /* actually a multidimensional array of desc[comm_size][num_sub_arrays][max_a_rank] */
  struct PPM_extent local_chunks[];
};

static size_t
init_local_mem(struct PPM_dist_mult_array *dm_array,
               enum PPM_dma_sync_mode sync_mode);
static void compute_max_win_size(struct PPM_dist_mult_array *dm_array);
static struct PPM_dist_mult_array *
dm_array_alloc(size_t cache_size,
               size_t num_sub_arrays,
               size_t max_sub_array_rank,
               int comm_size);
static bool
is_basic_type(MPI_Datatype dt);

static inline size_t
adjust_cache_size(size_t cache_size, enum PPM_dma_sync_mode sync_mode,
                  size_t max_sub_array_rank, int comm_size)
{
  if (sync_mode == PPM_dma_sync_mode_passive_target)
  {
    /* establish number of cache entries to create */
    if (cache_size == 0)
    {
      cache_size = (size_t)ceil(sqrt((double)comm_size));
      cache_size
        = cache_size < max_sub_array_rank ? max_sub_array_rank : cache_size;
    }
    /* implicitely limits cache size such that increment below cannot
     * overflow */
    cache_size = cache_size < (size_t)comm_size
      ? cache_size : (size_t)comm_size - 1;
  }
  else /*    sync_mode == PPM_dma_sync_mode_active_target
        * || sync_mode == PPM_dma_sync_mode_local_only */
    cache_size = 0;
  ++cache_size;
  return cache_size;
}

struct PPM_dist_mult_array *
PPM_dist_mult_array_new(size_t num_sub_arrays,
                        const struct PPM_global_array_desc *sub_arrays,
                        const struct PPM_extent *local_chunk,
                        MPI_Comm comm,
                        size_t cache_size,
                        enum PPM_dma_sync_mode sync_mode)
{
  /* find maximal rank of any sub array */
  size_t max_sub_array_rank = 0;
  /* this limit is MPI-imposed */
  assert(num_sub_arrays <= INT_MAX);
  assert(sync_mode == PPM_dma_sync_mode_passive_target
         || sync_mode == PPM_dma_sync_mode_active_target
         || sync_mode == PPM_dma_sync_mode_local_only);
  for (size_t i = 0; i < num_sub_arrays; ++i)
    if (max_sub_array_rank < sub_arrays[i].a_rank && sub_arrays[i].a_rank > 0)
      max_sub_array_rank = sub_arrays[i].a_rank;
  assert(max_sub_array_rank < PPM_dma_max_rank);
  /* get data about communication partners */
  int comm_size, comm_rank;
  xmpi(MPI_Comm_size(comm, &comm_size));
  xmpi(MPI_Comm_rank(comm, &comm_rank));
  cache_size
    = adjust_cache_size(cache_size, sync_mode, max_sub_array_rank, comm_size);
  struct PPM_dist_mult_array *dm_array
    = dm_array_alloc(cache_size, num_sub_arrays, max_sub_array_rank, comm_size);
  dm_array->num_sub_arrays = (unsigned)num_sub_arrays;
  dm_array->max_a_rank = (unsigned)max_sub_array_rank;
  dm_array->sync_mode = sync_mode;
  xmpi(MPI_Comm_dup(comm, &dm_array->comm));
  xmpi(MPI_Comm_dup(MPI_COMM_SELF, &dm_array->comm_self_dup));
  dm_array->comm_size = comm_size;
  dm_array->comm_rank = comm_rank;
  int defaultFlags = 0;
  {
    const char *ppm_tm_default_str
      = getenv("PPM_DIST_MULT_ARRAY_TRANSFER_MODE");
    if (ppm_tm_default_str && !strcmp(ppm_tm_default_str, "bytes"))
      defaultFlags |= PPM_dma_transfer_mode_bytes << 1;
  }
  dm_array->flags = defaultFlags;
  MPI_Aint *restrict dt_extents = dm_array->dt_extents;
  {
    MPI_Aint lb;
#if !defined __INTEL_COMPILER && !defined __PGI
    struct PPM_extent (*local_chunks)[max_sub_array_rank]
      = ((struct PPM_extent (*)[max_sub_array_rank])
         (dm_array->local_chunks
          + (size_t)comm_rank * num_sub_arrays * max_sub_array_rank));
#else
    struct PPM_extent *local_chunks = dm_array->local_chunks
      + (size_t)comm_rank * num_sub_arrays * max_sub_array_rank;
#endif
    MPI_Datatype *restrict element_dt = dm_array->element_dt;
    bool *restrict dt_is_basic = dm_array->dt_is_basic;
    for (size_t i = 0; i < num_sub_arrays; ++i)
    {
      bool is_basic = is_basic_type(sub_arrays[i].element_dt);
      dt_is_basic[i] = is_basic;
      if (!is_basic)
        xmpi(MPI_Type_dup(sub_arrays[i].element_dt, element_dt + i));
      else
        element_dt[i] = sub_arrays[i].element_dt;
      xmpi(MPI_Type_get_extent(element_dt[i], &lb, dt_extents + i));
      size_t a_rank = (size_t)sub_arrays[i].a_rank;
      dm_array->sub_arrays_global_desc[i].a_rank = (unsigned)a_rank;
      for (size_t j = 0; j < a_rank; ++j)
      {
#if !defined __INTEL_COMPILER && !defined __PGI
        local_chunks[i][j] = local_chunk[i * max_sub_array_rank + j];
#else
        local_chunks[i * max_sub_array_rank + j]
          = local_chunk[i * max_sub_array_rank + j];
#endif
        dm_array->sub_arrays_global_desc[i].rect[j]
          = sub_arrays[i].rect[j];
      }
    }
    /* it should be possible to use MPI_DATATYPE_NULL for sendtype, but
     * OpenMPI 1.4.x complains :-( */
    xmpi(MPI_Allgather(MPI_IN_PLACE, 0, PPM_extent_mp,
                       dm_array->local_chunks,
                       (int)(max_sub_array_rank * num_sub_arrays),
                       PPM_extent_mp, comm));
  }
  init_local_mem(dm_array, sync_mode);
  compute_max_win_size(dm_array);
  if (sync_mode == PPM_dma_sync_mode_passive_target)
  {
    xmpi(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, comm_rank,
                      MPI_MODE_NOCHECK, dm_array->win));
  }
  return dm_array;
}

struct PPM_dist_mult_array *
PPM_dist_mult_array_copy(struct PPM_dist_mult_array *dm_array)
{
  struct PPM_dist_mult_array *dm_array_new
    = dm_array_alloc(dm_array->cache_size, dm_array->num_sub_arrays,
                     dm_array->max_a_rank, dm_array->comm_size);
  xmpi(MPI_Comm_dup(dm_array->comm, &dm_array_new->comm));
  xmpi(MPI_Comm_dup(MPI_COMM_SELF, &dm_array_new->comm_self_dup));
  size_t num_sub_arrays
    = dm_array_new->num_sub_arrays = dm_array->num_sub_arrays;
  size_t max_a_rank = dm_array_new->max_a_rank = dm_array->max_a_rank;
  dm_array_new->cache_size = dm_array->cache_size;
  enum PPM_dma_sync_mode sync_mode
    = dm_array_new->sync_mode = dm_array->sync_mode;
  size_t comm_size = (size_t)(dm_array_new->comm_size = dm_array->comm_size),
    comm_rank = (size_t)(dm_array_new->comm_rank = dm_array->comm_rank);
  memcpy(dm_array_new->local_chunks, dm_array->local_chunks,
         sizeof (*dm_array_new->local_chunks) * comm_size * num_sub_arrays
         * max_a_rank);
  memcpy(dm_array_new->dt_extents, dm_array->dt_extents,
         sizeof (*dm_array_new->dt_extents) * num_sub_arrays);
  bool *dt_is_basic = dm_array_new->dt_is_basic;
  memcpy(dt_is_basic, dm_array->dt_is_basic,
         sizeof (*dm_array_new->dt_is_basic) * num_sub_arrays);
  MPI_Datatype *restrict element_dt = dm_array->element_dt,
    *restrict element_dt_new = dm_array_new->element_dt;
  for (size_t i = 0; i < num_sub_arrays; ++i)
    if (!dt_is_basic[i])
      xmpi(MPI_Type_dup(element_dt[i], element_dt_new + i));
    else
      element_dt_new[i] = element_dt[i];
  memcpy(dm_array_new->sub_arrays_global_desc, dm_array->sub_arrays_global_desc,
         sizeof (*dm_array_new->sub_arrays_global_desc) * num_sub_arrays);
  dm_array_new->flags = 0;
  size_t local_win_size = init_local_mem(dm_array_new, sync_mode);
  dm_array_new->max_win_size = dm_array->max_win_size;
  if (sync_mode == PPM_dma_sync_mode_passive_target)
  {
    xmpi(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, (int)comm_rank,
                      MPI_MODE_NOCHECK, dm_array_new->win));
  }
  memcpy(dm_array_new->cache[0].base, dm_array->cache[0].base,
         local_win_size);
  enum {
    fence_assert = MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE,
  };
  if (sync_mode == PPM_dma_sync_mode_active_target)
    xmpi(MPI_Win_fence(fence_assert, dm_array_new->win));
  return dm_array_new;
}

/* round size to next multiple of factor */
static inline size_t
roundUpToMultiple(size_t size, size_t factor)
{
  return (size + factor - 1)/factor * factor;
}

static void
resize_cache(struct PPM_dist_mult_array *dm_array,
             size_t cache_size);

static struct PPM_dist_mult_array *
dm_array_alloc(size_t cache_size,
               size_t num_sub_arrays,
               size_t max_sub_array_rank,
               int comm_size)
{
  struct PPM_dist_mult_array *dm_array;
  size_t local_chunks_size
    = ((sizeof (dm_array->local_chunks[0]) * num_sub_arrays * max_sub_array_rank
        * (size_t)comm_size)),
    header_size = (offsetof(struct PPM_dist_mult_array, local_chunks[0])
                   + local_chunks_size),
    dt_extents_size = sizeof (dm_array->dt_extents[0]) * num_sub_arrays,
    element_dt_size = sizeof (dm_array->element_dt[0]) * num_sub_arrays,
    sub_arrays_global_desc_size = (sizeof (*dm_array->sub_arrays_global_desc)
                                   * num_sub_arrays),
    dt_is_basic_size = sizeof (*dm_array->dt_is_basic) * num_sub_arrays,
    dt_extents_ofs = roundUpToMultiple(header_size,
                                       __alignof(MPI_Aint)),
    element_dt_ofs = roundUpToMultiple(dt_extents_ofs + dt_extents_size,
                                       __alignof(MPI_Datatype)),
    sub_arrays_global_desc_ofs
    = roundUpToMultiple(element_dt_ofs + element_dt_size,
                        __alignof(struct PPM_global_array_desc_)),
    dt_is_basic_ofs = roundUpToMultiple(sub_arrays_global_desc_ofs
                                        + sub_arrays_global_desc_size,
                                        __alignof(bool)),
    dm_array_alloc_size = dt_is_basic_ofs + dt_is_basic_size;
  dm_array = xmalloc(dm_array_alloc_size);
  dm_array->dt_extents
    = (MPI_Aint *)((unsigned char *)dm_array + dt_extents_ofs);
  dm_array->dt_is_basic
    = (bool *)((unsigned char *)dm_array + dt_is_basic_ofs);
  dm_array->element_dt
    = (MPI_Datatype *)((unsigned char *)dm_array + element_dt_ofs);
  dm_array->sub_arrays_global_desc
    = (struct PPM_global_array_desc_ *)
  ((unsigned char *)dm_array + sub_arrays_global_desc_ofs);
  dm_array->cache_size = 0;
  dm_array->cache = NULL;
  dm_array->num_sub_arrays = (unsigned)num_sub_arrays;
  resize_cache(dm_array, cache_size);
  dm_array->access_stamp = 1;
  dm_array->valid_stamp = 1;
  return dm_array;
}

static void
delete_cache_entry(struct PPM_dm_array_cache_entry *entry);

static void
resize_cache(struct PPM_dist_mult_array *dm_array,
             size_t cache_size)
{
  size_t cache_size_prev = dm_array->cache_size;
  struct PPM_dm_array_cache_entry *cache = dm_array->cache;
  MPI_Aint *offsets_prev,
    *offsets = offsets_prev = cache ? cache[0].offset : (MPI_Aint *)NULL;
  if (cache_size < cache_size_prev)
    for (size_t i = cache_size; i < cache_size_prev; ++i)
      delete_cache_entry(cache + i);
  dm_array->cache = cache = xrealloc(cache, cache_size * sizeof(*cache));
  size_t num_sub_arrays = dm_array->num_sub_arrays;
  offsets = xrealloc(offsets, cache_size * num_sub_arrays * sizeof(*offsets));
  /* update offsets if changed */
  if (offsets != offsets_prev)
  {
    size_t cache_rewrite_size = cache_size_prev <= cache_size ?
      cache_size_prev : cache_size;
    for (size_t i = 0; i < cache_rewrite_size; ++i)
      cache[i].offset = offsets + i * num_sub_arrays;
  }
  /* initialize newly created entries */
  for (size_t i = cache_size_prev; i < cache_size; ++i)
  {
    cache[i].base = NULL;
    cache[i].offset = offsets + i * num_sub_arrays;
    cache[i].access_stamp = 0;
    cache[i].rank = -1;
    cache[i].composite_dt = MPI_DATATYPE_NULL;
    cache[i].composite_count = 1;
  }
  dm_array->cache_size = cache_size;
}


static void
create_mpi_win(struct PPM_dist_mult_array *dm_array,
               size_t local_win_size,
               enum PPM_dma_sync_mode sync_mode);

static size_t
get_local_win_size(struct PPM_dist_mult_array *dm_array);

static size_t
init_local_mem(struct PPM_dist_mult_array *dm_array,
               enum PPM_dma_sync_mode sync_mode)
{
  size_t local_win_size = get_local_win_size(dm_array);
  xmpi(MPI_Alloc_mem((MPI_Aint)local_win_size, MPI_INFO_NULL,
                     &dm_array->cache[0].base));
  create_mpi_win(dm_array, local_win_size, sync_mode);
  dm_array->flags &= ~exposed_flag;
  dm_array->cache[0].rank = dm_array->comm_rank;
  return local_win_size;
}

static void
create_mpi_win(struct PPM_dist_mult_array *dm_array,
               size_t local_win_size,
               enum PPM_dma_sync_mode sync_mode)
{
  if (sync_mode == PPM_dma_sync_mode_passive_target
      || sync_mode == PPM_dma_sync_mode_active_target)
  {
    MPI_Info win_create_info;
    xmpi(MPI_Info_create(&win_create_info));
    xmpi(MPI_Info_set(win_create_info, "no_locks",
                      sync_mode == PPM_dma_sync_mode_passive_target
                      ? "false" : "true"));
    xmpi(MPI_Win_create(dm_array->cache[0].base, (MPI_Aint)local_win_size, 1,
                        win_create_info, dm_array->comm, &dm_array->win));
    xmpi(MPI_Info_free(&win_create_info));
  }
}

static void *
PPM_align_address_for_dt(const void *addr, MPI_Datatype dt_prev,
                         MPI_Datatype dtype2align_for,
                         MPI_Aint dtype2align_for_extent,
                         bool dtype2align_for_is_basic)
{
  intptr_t iaddr = (intptr_t)addr, oaddr;
  if (dt_prev == MPI_DATATYPE_NULL || dt_prev == dtype2align_for)
    oaddr = iaddr;
  else if (dtype2align_for_is_basic)
    /* align basic type to multiple of its extent */
    oaddr = ((iaddr + dtype2align_for_extent - 1)
             / dtype2align_for_extent) * dtype2align_for_extent;
  else
  {
#if PPM_MAXIMUM_ALIGNMENT & (PPM_MAXIMUM_ALIGNMENT - 1) == 0
    oaddr = (iaddr + PPM_MAXIMUM_ALIGNMENT - 1) & -PPM_MAXIMUM_ALIGNMENT;
#else
    oaddr
      = ((iaddr + PPM_MAXIMUM_ALIGNMENT - 1) / PPM_MAXIMUM_ALIGNMENT)
      * PPM_MAXIMUM_ALIGNMENT;
#endif
  }
  return (void *)oaddr;
}



static void
compute_max_win_size(struct PPM_dist_mult_array *dm_array)
{
  size_t max_remote_win_size = 0, comm_size = (size_t)dm_array->comm_size,
    num_sub_arrays = dm_array->num_sub_arrays,
    max_sub_array_rank = dm_array->max_a_rank;
  const struct PPM_global_array_desc_ *sub_arrays
    = dm_array->sub_arrays_global_desc;
  const struct PPM_extent (*chunk)[num_sub_arrays][max_sub_array_rank]
    = (const struct PPM_extent (*)[num_sub_arrays][max_sub_array_rank])
    dm_array->local_chunks;
  const MPI_Datatype *restrict element_dt = dm_array->element_dt;
  const MPI_Aint *restrict dt_extents = dm_array->dt_extents;
  const bool *restrict dt_is_basic = dm_array->dt_is_basic;
  for (size_t rank = 0; rank < comm_size; ++rank)
  {
    size_t remote_win_size = 0;
    MPI_Datatype prev_dt = MPI_DATATYPE_NULL;
    for (size_t i = 0; i < num_sub_arrays; ++i)
    {
      remote_win_size = (size_t)PPM_align_address_for_dt(
        (void *)(remote_win_size
                 + ((size_t)dt_extents[i]
                    * (size_t)PPM_extents_size(sub_arrays[i].a_rank,
                                               chunk[rank][i]))),
        prev_dt, element_dt[i], dt_extents[i], dt_is_basic[i]);
      prev_dt = element_dt[i];
    }
    max_remote_win_size = max_remote_win_size >= remote_win_size
      ? max_remote_win_size : remote_win_size;
  }
  dm_array->max_win_size = (MPI_Aint)max_remote_win_size;
}

static size_t
subarray_offset(size_t sub_array_idx,
                size_t max_sub_array_rank,
                const struct PPM_global_array_desc_ *sub_arrays,
                const struct PPM_extent (*chunk)[max_sub_array_rank],
                const MPI_Datatype *element_dt,
                const MPI_Aint *dt_extents,
                const bool *dt_is_basic)
{
  size_t ofs = 0;
  MPI_Datatype dt_prev = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < sub_array_idx; ++i)
  {
    ofs = (size_t)PPM_align_address_for_dt(
      (void *)(ofs
               + ((size_t)dt_extents[i]
                  * (size_t)PPM_extents_size(sub_arrays[i].a_rank, chunk[i]))),
      dt_prev, element_dt[i], dt_extents[i], dt_is_basic[i]);
    dt_prev = element_dt[i];
  }
  return ofs;
}


/**
 * computes offset each sub-array is stored at in cache memory
 * and returns total size of corresponding window area
 */
static size_t
compute_cache_addr(size_t num_sub_arrays,
                   size_t max_sub_array_rank,
                   struct PPM_dm_array_cache_entry *cache_entry,
                   const struct PPM_global_array_desc_ *sub_arrays,
                   const struct PPM_extent (*chunk)[max_sub_array_rank],
                   const MPI_Datatype *element_dt,
                   const MPI_Aint *dt_extents,
                   const bool *dt_is_basic)
{
  size_t win_size = 0;
  MPI_Datatype dt_prev = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < num_sub_arrays; ++i)
  {
    cache_entry->offset[i] = (MPI_Aint)win_size;
    win_size = (size_t)PPM_align_address_for_dt(
      (void *)(win_size
               + ((size_t)dt_extents[i]
                  * (size_t)PPM_extents_size(sub_arrays[i].a_rank, chunk[i]))),
      dt_prev, element_dt[i], dt_extents[i], dt_is_basic[i]);
    dt_prev = element_dt[i];
  }
  return win_size;
}

static size_t
get_local_win_size(struct PPM_dist_mult_array *dm_array)
{
  int comm_rank = dm_array->comm_rank;
  size_t max_sub_array_rank = dm_array->max_a_rank;
  size_t num_sub_arrays = dm_array->num_sub_arrays;
  const struct PPM_extent (*local_chunks)[max_sub_array_rank]
    = (const struct PPM_extent (*)[max_sub_array_rank])
    (dm_array->local_chunks
     + max_sub_array_rank * (size_t)comm_rank * num_sub_arrays);
  size_t local_win_size
    = compute_cache_addr(num_sub_arrays, max_sub_array_rank,
                         dm_array->cache + 0,
                         dm_array->sub_arrays_global_desc,
                         local_chunks,
                         dm_array->element_dt,
                         dm_array->dt_extents,
                         dm_array->dt_is_basic);
  return local_win_size;
}


/**
 * determine if an MPI datatype is a basic datatype
 */
static bool
is_basic_type(MPI_Datatype dt)
{
  bool is_basic_type_;
  int num_int, num_addr, num_dt, combiner;
  xmpi(MPI_Type_get_envelope(dt, &num_int, &num_addr, &num_dt, &combiner));
  if (combiner == MPI_COMBINER_NAMED)
    is_basic_type_ = true;
  else if (combiner == MPI_COMBINER_DUP)
  {
    MPI_Aint cont_a[1];
    MPI_Datatype cont_dt[1];
    int cont_i[1];
    xmpi(MPI_Type_get_contents(dt, num_int, num_addr, num_dt,
                               cont_i, cont_a, cont_dt));
    is_basic_type_ = is_basic_type(cont_dt[0]);
    xmpi(MPI_Type_get_envelope(cont_dt[0], &num_int, &num_addr, &num_dt,
                               &combiner));

    if (combiner != MPI_COMBINER_NAMED)
      xmpi(MPI_Type_free(cont_dt));
  }
  else
    is_basic_type_ = false;
  return is_basic_type_;
}

void
PPM_dist_mult_array_expose(struct PPM_dist_mult_array *dm_array);
void
PPM_dist_mult_array_unexpose(struct PPM_dist_mult_array *dm_array);

void
PPM_dist_mult_array_delete(struct PPM_dist_mult_array *dm_array)
{
  if (dm_array->flags & exposed_flag)
    PPM_dist_mult_array_unexpose(dm_array);
  int comm_rank = dm_array->comm_rank;
  enum PPM_dma_sync_mode sync_mode = dm_array->sync_mode;
  if (sync_mode == PPM_dma_sync_mode_passive_target)
    xmpi(MPI_Win_unlock(comm_rank, dm_array->win));
  if (sync_mode != PPM_dma_sync_mode_local_only)
    xmpi(MPI_Win_free(&dm_array->win));
  {
    size_t cache_size = dm_array->cache_size;
    for (size_t i = 0; i < cache_size; ++i)
      delete_cache_entry(dm_array->cache + i);
  }
  size_t num_sub_arrays = dm_array->num_sub_arrays;
  MPI_Datatype *restrict element_dt = dm_array->element_dt;
  const bool *restrict dt_is_basic = dm_array->dt_is_basic;
  for (size_t i = 0; i < num_sub_arrays; ++i)
    if (!dt_is_basic[i])
      xmpi(MPI_Type_free(element_dt + i));
  xmpi(MPI_Comm_free(&dm_array->comm));
  xmpi(MPI_Comm_free(&dm_array->comm_self_dup));
  free(dm_array->cache[0].offset);
  free(dm_array->cache);
  free(dm_array);
}

static void
free_cache_entry_dt(struct PPM_dm_array_cache_entry *entry)
{
  MPI_Datatype composite_dt = entry->composite_dt;
  if (entry->rank >= 0 && composite_dt != MPI_DATATYPE_NULL
      && composite_dt != MPI_BYTE)
    xmpi(MPI_Type_free(&entry->composite_dt));
}

static void
delete_cache_entry(struct PPM_dm_array_cache_entry *entry)
{
  free_cache_entry_dt(entry);
  if (entry->base)
    xmpi(MPI_Free_mem(entry->base));
}

/* retrieve TYPE(c_ptr) corresponding to a cached sub-array chunk */
static void *
dist_mult_array_cache_chunk_ptr_c(struct PPM_dist_mult_array *dm_array,
                                  size_t sub_array_idx,
                                  size_t cache_idx)
{
  return (void *)((unsigned char *)dm_array->cache[cache_idx].base
                  + dm_array->cache[cache_idx].offset[sub_array_idx]);
}

void *
PPM_dist_mult_array_local_ptr(struct PPM_dist_mult_array *dm_array,
                              size_t sub_array_idx)
{
  return dist_mult_array_cache_chunk_ptr_c(dm_array, sub_array_idx, 0);
}

void
PPM_dist_mult_array_expose(struct PPM_dist_mult_array *dm_array)
{
  enum {
    fence_assert = MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE,
  };
  assert(!(dm_array->flags & exposed_flag));
  if (dm_array->sync_mode == PPM_dma_sync_mode_passive_target)
  {
    xmpi(MPI_Win_unlock(dm_array->comm_rank, dm_array->win));
    xmpi(MPI_Barrier(dm_array->comm));
  }
  else if (dm_array->sync_mode == PPM_dma_sync_mode_active_target)
    xmpi(MPI_Win_fence(fence_assert, dm_array->win));
  dm_array->flags |= exposed_flag;
  dm_array->valid_stamp = dm_array->access_stamp + 1;
}

void
PPM_dist_mult_array_unexpose(struct PPM_dist_mult_array *dm_array)
{
  enum {
    fence_assert = MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED,
  };
  assert(dm_array->flags & exposed_flag);
  if (dm_array->sync_mode == PPM_dma_sync_mode_passive_target)
  {
    xmpi(MPI_Barrier(dm_array->comm));
    xmpi(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, dm_array->comm_rank,
                      MPI_MODE_NOCHECK, dm_array->win));
  }
  else if (dm_array->sync_mode == PPM_dma_sync_mode_active_target)
    xmpi(MPI_Win_fence(fence_assert, dm_array->win));
  dm_array->flags &= ~exposed_flag;
}

void
PPM_dist_mult_array_rma_sync(struct PPM_dist_mult_array *dm_array)
{
  enum {
    fence_assert = MPI_MODE_NOSTORE | MPI_MODE_NOPUT,
  };
  assert((dm_array->flags & exposed_flag)
         && dm_array->sync_mode == PPM_dma_sync_mode_active_target);
  /* dm_array%sync_mode == sync_mode_active_target */
  xmpi(MPI_Win_fence(fence_assert, dm_array->win));
}

int
PPM_dist_mult_array_get_sync_mode(struct PPM_dist_mult_array *dm_array)
{
  return dm_array->sync_mode;
}

void
PPM_dist_mult_array_set_sync_mode(struct PPM_dist_mult_array *dm_array,
                                  enum PPM_dma_sync_mode sync_mode,
                                  size_t cache_size)
{
  enum {
    fence_assert = MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED,
  };
  enum PPM_dma_sync_mode prev_sync_mode = dm_array->sync_mode;
  if (dm_array->flags & exposed_flag)
    PPM_dist_mult_array_unexpose(dm_array);
  if (prev_sync_mode != sync_mode)
  {
    if (prev_sync_mode != PPM_dma_sync_mode_local_only)
    {
      if (prev_sync_mode == PPM_dma_sync_mode_passive_target)
        xmpi(MPI_Win_unlock(dm_array->comm_rank, dm_array->win));
      xmpi(MPI_Win_free(&dm_array->win));
    }
    cache_size
      = adjust_cache_size(cache_size, sync_mode, dm_array->max_a_rank,
                          dm_array->comm_size);
    resize_cache(dm_array, cache_size);
    create_mpi_win(dm_array, get_local_win_size(dm_array), sync_mode);
    if (sync_mode == PPM_dma_sync_mode_passive_target)
      xmpi(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, dm_array->comm_rank,
                        MPI_MODE_NOCHECK, dm_array->win));
    dm_array->sync_mode = sync_mode;
  }
  else if (sync_mode == PPM_dma_sync_mode_passive_target)
  {
    cache_size
      = adjust_cache_size(cache_size, sync_mode, dm_array->max_a_rank,
                          dm_array->comm_size);
    if (cache_size != dm_array->cache_size)
      resize_cache(dm_array, cache_size);
  }
}


/**
 * perform passive mode RMA get operation
 */
static void cache_get(struct PPM_dist_mult_array *dm_array,
                      size_t cache_entry, int rank)
{
  struct PPM_dm_array_cache_entry *cache = dm_array->cache + cache_entry;
  void *baseptr = cache->base;
  MPI_Win win = dm_array->win;
  MPI_Datatype composite_dt = cache->composite_dt;
  int composite_count = cache->composite_count;
  xmpi(MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, win));
  xmpi(MPI_Get(baseptr, composite_count, composite_dt,
               rank, 0, composite_count, composite_dt, win));
  xmpi(MPI_Win_unlock(rank, win));
}

/**
 * map remote memory of rank to cache entry and return cache entry index
 */
static size_t
dist_mult_array_cache_rank(struct PPM_dist_mult_array *dm_array, int rank)
{
  if (rank == dm_array->comm_rank)
    return 0;
  size_t num_sub_arrays = dm_array->num_sub_arrays;
  assert(dm_array->flags & exposed_flag);
  unsigned long stamp = dm_array->access_stamp + 1, oldest_stamp = stamp;
  size_t oldest = SIZE_MAX;
  size_t cache_size = dm_array->cache_size;
  for (size_t i = 1; i < cache_size; ++i)
    if (dm_array->cache[i].rank == rank)
    {
      if (dm_array->cache[i].access_stamp < dm_array->valid_stamp)
        cache_get(dm_array, i, rank);
      dm_array->cache[i].access_stamp = stamp;
      dm_array->access_stamp = stamp;
      return i;
    }
    else if (dm_array->cache[i].access_stamp < oldest_stamp)
    {
      oldest_stamp = dm_array->cache[i].access_stamp;
      oldest = i;
    }
  if (dm_array->cache[oldest].rank >= 0)
    /* free cache resources here */
    free_cache_entry_dt(dm_array->cache + oldest);
  /* (re-)populate free cache entry */
  if (dm_array->cache[oldest].base == NULL)
    xmpi(MPI_Alloc_mem(dm_array->max_win_size, MPI_INFO_NULL,
                       &dm_array->cache[oldest].base));
  /* create datatype for MPI_Get */
  int current_transfer_mode = (dm_array->flags & transfer_mode_mask) >> 1;
  size_t max_a_rank = dm_array->max_a_rank;
  const struct PPM_extent (*local_chunks)[max_a_rank]
    = (const struct PPM_extent (*)[max_a_rank])
    (dm_array->local_chunks + num_sub_arrays * max_a_rank * (size_t)rank);
  size_t remote_win_size
    = compute_cache_addr(num_sub_arrays, max_a_rank, dm_array->cache + oldest,
                         dm_array->sub_arrays_global_desc,
                         local_chunks, dm_array->element_dt,
                         dm_array->dt_extents, dm_array->dt_is_basic);
  switch (current_transfer_mode)
  {
  case PPM_dma_transfer_mode_struct:
    {
      enum { autoComponentCount = 16 };
      int auto_component_counts[autoComponentCount],
        *component_count
        = num_sub_arrays <= autoComponentCount
        ? auto_component_counts
        : xmalloc(num_sub_arrays * sizeof (*component_count));
      struct PPM_global_array_desc_ *sub_arrays_global_desc
        = dm_array->sub_arrays_global_desc;
      for (size_t i = 0; i < num_sub_arrays; ++i)
      {
        component_count[i]
          = PPM_extents_size(sub_arrays_global_desc[i].a_rank, local_chunks[i]);
      }
      xmpi(MPI_Type_create_struct((int)num_sub_arrays, component_count,
                                  dm_array->cache[oldest].offset,
                                  dm_array->element_dt,
                                  &dm_array->cache[oldest].composite_dt));
      dm_array->cache[oldest].composite_count = 1;
      xmpi(MPI_Type_commit(&dm_array->cache[oldest].composite_dt));
      if (num_sub_arrays > autoComponentCount)
        free(component_count);
    }
    break;
  case PPM_dma_transfer_mode_bytes:
    dm_array->cache[oldest].composite_count = (int)remote_win_size;
    dm_array->cache[oldest].composite_dt = MPI_BYTE;
    break;
  }
  cache_get(dm_array, oldest, rank);
  dm_array->cache[oldest].rank = rank;
  dm_array->cache[oldest].access_stamp = stamp;
  dm_array->access_stamp = stamp;
  return oldest;
}

static int
dist_mult_array_coord2rank(struct PPM_dist_mult_array *dm_array,
                           size_t sub_array, const int32_t coord[])
{
  size_t num_sub_arrays = dm_array->num_sub_arrays,
    max_sub_array_rank = dm_array->max_a_rank;
  unsigned comm_rank = (unsigned)dm_array->comm_rank,
    comm_size = (unsigned)dm_array->comm_size;
  unsigned max_dist = ((unsigned)comm_size + 1)/2;
  struct PPM_extent (*local_chunks)[num_sub_arrays][max_sub_array_rank]
    = (struct PPM_extent (*)[num_sub_arrays][max_sub_array_rank])
    dm_array->local_chunks;
  struct PPM_global_array_desc_ *sub_array_global_desc
    = dm_array->sub_arrays_global_desc + sub_array;
  size_t sub_array_rank = sub_array_global_desc->a_rank;
  for (size_t dist = 0; dist <= max_dist; ++dist)
  {
    size_t cand_rank = (comm_rank + dist) % comm_size,
      cand_rank_other = (comm_rank + comm_size - dist) % comm_size;
    if (PPM_coord_is_contained_in_extents(
          sub_array_rank, coord, local_chunks[cand_rank][sub_array]))
      return (int)cand_rank;
    else if (PPM_coord_is_contained_in_extents(
               sub_array_rank, coord, local_chunks[cand_rank_other][sub_array]))
      return (int)cand_rank_other;
  }
  PPM_abort(dm_array->comm, "invalid coordinate", __FILE__, __LINE__);
}


/**
 * establish remote MPI rank holding data of sub-array at coordinate
 * and fetch data if not already in-cache
 */
static size_t
dist_mult_array_get_cache_idx(struct PPM_dist_mult_array *dm_array,
                              size_t sub_array,
                              int32_t coord[])
{
#ifndef NDEBUG
  size_t num_sub_arrays = dm_array->num_sub_arrays;
#endif
  assert(sub_array < num_sub_arrays);
  assert(PPM_coord_is_contained_in_extents(
           dm_array->sub_arrays_global_desc[sub_array].a_rank,
           coord, dm_array->sub_arrays_global_desc[sub_array].rect));
  return dist_mult_array_cache_rank(
    dm_array, dist_mult_array_coord2rank(dm_array, sub_array, coord));
}

static void
dist_mult_array_get_deferred(struct PPM_dist_mult_array *dm_array,
                             size_t sub_array, int32_t coord[], void *v_out)
{
  size_t num_sub_arrays = dm_array->num_sub_arrays,
    max_sub_array_rank = dm_array->max_a_rank;
  int src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord);
  struct PPM_global_array_desc_ *restrict sub_arrays_global_desc
    = dm_array->sub_arrays_global_desc;
  const struct PPM_extent (*chunk)[num_sub_arrays][max_sub_array_rank]
    = (const struct PPM_extent (*)[num_sub_arrays][max_sub_array_rank])
    dm_array->local_chunks;
  size_t ofs
    = (src_comm_rank == dm_array->comm_rank)
    ? (size_t)dm_array->cache[0].offset[sub_array]
    : subarray_offset(sub_array, max_sub_array_rank, sub_arrays_global_desc,
                      chunk[src_comm_rank], dm_array->element_dt,
                      dm_array->dt_extents, dm_array->dt_is_basic);
  unsigned a_rank = sub_arrays_global_desc[sub_array].a_rank;
  int32_t coord_base[PPM_dma_max_rank];
  for (size_t i = 0; i < a_rank; ++i)
    coord_base[i] = chunk[src_comm_rank][sub_array][i].first;
  MPI_Aint byte_offset = (MPI_Aint)ofs,
    ofs_factor = dm_array->dt_extents[sub_array];
  for (size_t i = 0; i < a_rank; ++i)
  {
    byte_offset += (MPI_Aint)(coord[i] - coord_base[i]) * ofs_factor;
    ofs_factor
      *= chunk[src_comm_rank][sub_array][i].size;
  }
  MPI_Datatype dt = dm_array->element_dt[sub_array];
  if (src_comm_rank == dm_array->comm_rank)
    xmpi(MPI_Sendrecv((unsigned char *)dm_array->cache[0].base + byte_offset,
                      1, dt, 0, 0, v_out, 1, dt, 0, 0,
                      MPI_COMM_SELF, MPI_STATUS_IGNORE));
  else
    xmpi(MPI_Get(v_out, 1, dt,
                 src_comm_rank, byte_offset, 1, dt, dm_array->win));
}



void
PPM_dist_mult_array_get(struct PPM_dist_mult_array *dm_array,
                        size_t sub_array,
                        int32_t coord[],
                        void *v_out)
{
  if (dm_array->sync_mode == PPM_dma_sync_mode_passive_target
      || dm_array->sync_mode == PPM_dma_sync_mode_local_only)
  {
    size_t cache_idx;
    if (dm_array->sync_mode == PPM_dma_sync_mode_passive_target)
      cache_idx = dist_mult_array_get_cache_idx(dm_array, sub_array, coord);
    else /* dm_array->sync_mode == PPM_dma_sync_mode_local_only */
    {
#ifndef NDEBUG
      size_t num_sub_arrays = dm_array->num_sub_arrays,
        max_sub_array_rank = dm_array->max_a_rank;
      assert(sub_array < num_sub_arrays);
      struct PPM_extent (*local_chunks)[num_sub_arrays][max_sub_array_rank]
        = (struct PPM_extent (*)[num_sub_arrays][max_sub_array_rank])
        dm_array->local_chunks;
      size_t sub_array_rank
        = dm_array->sub_arrays_global_desc[sub_array].a_rank;
      size_t comm_rank = (size_t)dm_array->comm_rank;
#endif
      assert(PPM_coord_is_contained_in_extents(
               sub_array_rank, coord, local_chunks[comm_rank][sub_array]));
      cache_idx = 0;
    }
    size_t sub_array_rank = dm_array->sub_arrays_global_desc[sub_array].a_rank,
      top_sub_array_rank = sub_array_rank - 1;
    struct PPM_dm_array_cache_entry *cache_entry = dm_array->cache + cache_idx;
    size_t rank = (size_t)cache_entry->rank,
      num_sub_arrays = dm_array->num_sub_arrays,
      max_sub_array_rank = dm_array->max_a_rank;
    struct PPM_extent *chunk = dm_array->local_chunks
      + rank * num_sub_arrays * max_sub_array_rank
      + sub_array * max_sub_array_rank;
    size_t dt_extent = (size_t)dm_array->dt_extents[sub_array];
    size_t ofs = (size_t)cache_entry->offset[sub_array]
      + ((size_t)(coord[top_sub_array_rank] - chunk[top_sub_array_rank].first)
         * dt_extent);
    size_t accum = 1;
    for (size_t i = top_sub_array_rank; i > 0; --i)
    {
      ofs += (size_t)(coord[i-1] - chunk[i-1].first) * accum * dt_extent;
      accum *= (size_t)chunk[i-1].size;
    }
    void *v_in = (unsigned char *)cache_entry->base + ofs;
    xmpi(MPI_Sendrecv(v_in, 1, dm_array->element_dt[sub_array], 0, 0,
                      v_out, 1, dm_array->element_dt[sub_array], 0, 0,
                      dm_array->comm_self_dup, MPI_STATUS_IGNORE));
  }
  else /* dm_array->sync_mode == PPM_dma_sync_mode_active_target */
    dist_mult_array_get_deferred(dm_array, sub_array, coord, v_out);
}

unsigned
PPM_dist_mult_array_a_rank(struct PPM_dist_mult_array *dm_array,
                           size_t sub_array_idx)
{
  return dm_array->sub_arrays_global_desc[sub_array_idx].a_rank;
}


void
PPM_dist_mult_array_rank_rect(struct PPM_dist_mult_array *dm_array,
                              size_t sub_array_idx,
                              int rank,
                              struct PPM_extent *rect)
{
  size_t rank_ = (rank < 0) ? (size_t)dm_array->comm_rank : (size_t)rank;
  size_t max_a_rank = dm_array->max_a_rank,
    a_rank = dm_array->sub_arrays_global_desc[sub_array_idx].a_rank,
    num_sub_arrays = dm_array->num_sub_arrays;
  struct PPM_extent *restrict rank_chunks
    = dm_array->local_chunks + (rank_ * num_sub_arrays
                                + sub_array_idx) * max_a_rank;
  for (size_t i = 0; i < a_rank; ++i)
    rect[i] = rank_chunks[i];
}

MPI_Comm
PPM_dist_mult_array_comm(struct PPM_dist_mult_array *dm_array)
{
  return dm_array->comm;
}

void
PPM_dist_mult_array_set_transfer_mode(struct PPM_dist_mult_array *dm_array,
                                      int mode)
{
  int flags = dm_array->flags,
    current_transfer_mode = (flags & transfer_mode_mask) >> 1;
  if ((mode != PPM_dma_transfer_mode_struct
       && mode != PPM_dma_transfer_mode_bytes)
      || (mode == PPM_dma_transfer_mode_bytes
          && dm_array->sync_mode != PPM_dma_sync_mode_passive_target))
    PPM_abort(dm_array->comm, "invalid transfer mode requested",
              __FILE__, __LINE__);
  if (current_transfer_mode != mode)
    dm_array->flags = ((flags & ~transfer_mode_mask) | (mode << 1));
}

int
PPM_dist_mult_array_get_transfer_mode(struct PPM_dist_mult_array *dm_array)
{
  return (dm_array->flags & transfer_mode_mask) >> 1;
}


/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
