/**
 * @file xt_idxsection.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include <config.h>
#endif

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt_arithmetic_util.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxsection.h"
#include "xt_idxsection_internal.h"
#include "xt/xt_mpi.h"
#include "xt/mergesort.h"
#include "xt/xt_sort.h"
#include "xt/quicksort.h"
#include "xt_idxlist_unpack.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "instr.h"

static void
idxsection_delete(Xt_idxlist data);

static size_t
idxsection_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxsection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm);

static Xt_idxlist
idxsection_copy(Xt_idxlist idxlist);

static void
idxsection_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static const Xt_int *
idxsection_get_indices_const(Xt_idxlist idxlist);

static void
idxsection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             int * num_stripes);

static int
idxsection_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index);

static int
idxsection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position);

static int
idxsection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                     int * position, int offset);
static size_t
idxsection_get_positions_of_indices(Xt_idxlist body_idxlist,
                                    Xt_int const *selection_idx,
                                    size_t num_selection, int *positions,
                                    int single_match_only);

static Xt_int
idxsection_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxsection_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxsection_vtable = {
              .delete                       = idxsection_delete,
              .get_pack_size                = idxsection_get_pack_size,
              .pack                         = idxsection_pack,
              .copy                         = idxsection_copy,
              .get_indices                  = idxsection_get_indices,
              .get_indices_const            = idxsection_get_indices_const,
              .get_index_stripes            = idxsection_get_index_stripes,
              .get_index_at_position        = idxsection_get_index_at_position,
              .get_indices_at_positions     = NULL,
              .get_position_of_index        = idxsection_get_position_of_index,
              .get_positions_of_indices     = idxsection_get_positions_of_indices,
              .get_position_of_index_off    = idxsection_get_position_of_index_off,
              .get_positions_of_indices_off = NULL,
              .get_min_index                = idxsection_get_min_index,
              .get_max_index                = idxsection_get_max_index,
              .get_bounding_box             = NULL,
              .idxlist_pack_code            = SECTION,
};

/* descriptor for per-dimension extent and stride */
struct dim_desc
{
  Xt_int global_size, global_stride, local_start, local_stride;
  int local_size;
};

static MPI_Datatype dim_desc_dt;


typedef struct Xt_idxsection_ *Xt_idxsection;

struct Xt_idxsection_ {

  struct Xt_idxlist_ parent;

  Xt_int *index_array_cache;

  Xt_int global_start_index;
  Xt_int local_start_index;

  Xt_int min_index_cache;
  Xt_int max_index_cache;
  int ndim;
  struct dim_desc dims[];
};

static int
idxsection_get_num_indices(Xt_idxsection section);

#if __GNUC__ >= 11
__attribute__ ((access (none, 1)))
int MPI_Get_address(XT_MPI_SEND_BUF_CONST void *location, MPI_Aint *address);
#endif

void
xt_idxsection_initialize(void)
{
  struct dim_desc dim_desc;

  MPI_Aint base_address, local_size_address;

  MPI_Get_address(&dim_desc, &base_address);
  MPI_Get_address(&dim_desc.local_size, &local_size_address);

  enum { num_dt_components = 2 };
  int block_lengths[num_dt_components] = { 4, 1 };
  MPI_Aint displacements[num_dt_components]
    = {0, local_size_address - base_address };
  MPI_Datatype types[num_dt_components]
    = { Xt_int_dt, MPI_INT },
    dim_desc_dt_unaligned;
  xt_mpi_call(MPI_Type_create_struct(num_dt_components,
                                     block_lengths, displacements, types,
                                     &dim_desc_dt_unaligned), Xt_default_comm);
  xt_mpi_call(MPI_Type_create_resized(dim_desc_dt_unaligned, 0,
                                      (MPI_Aint)sizeof(dim_desc),
                                      &dim_desc_dt), Xt_default_comm);
  xt_mpi_call(MPI_Type_free(&dim_desc_dt_unaligned), Xt_default_comm);
  xt_mpi_call(MPI_Type_commit(&dim_desc_dt), Xt_default_comm);
}

void
xt_idxsection_finalize(void)
{
  xt_mpi_call(MPI_Type_free(&dim_desc_dt), Xt_default_comm);
}

Xt_idxlist xt_idxsection_new(Xt_int start, int num_dimensions,
                             const Xt_int global_size[num_dimensions],
                             const int local_size[num_dimensions],
                             const Xt_int local_start[num_dimensions]) {

  INSTR_DEF(instr,"xt_idxsection_new")
  INSTR_START(instr);
  // ensure that yaxt is initialized
  assert(xt_initialized());

  Xt_idxsection idxsection = NULL;
  int num_indices;
  if (num_dimensions > 0) {
    idxsection = xmalloc(sizeof (*idxsection)
                         + (size_t)num_dimensions *
                         sizeof (idxsection->dims[0]));

    idxsection->global_start_index = start;
    idxsection->ndim               = num_dimensions;

    idxsection->index_array_cache  = NULL;

    for (int i = 0; i < num_dimensions; ++i) {
      idxsection->dims[i].global_size = global_size[i];
      idxsection->dims[i].local_size  = local_size[i];
      idxsection->dims[i].local_start = local_start[i];
    }
    num_indices = idxsection_get_num_indices(idxsection);
  } else {
    num_indices = 0;
  }
  if (num_indices == 0) {
    free(idxsection);
    return xt_idxempty_new();
  }
  Xt_idxlist_init(&idxsection->parent, &idxsection_vtable, num_indices);
  idxsection->local_start_index = start;
  idxsection->dims[num_dimensions - 1].global_stride =
    (Xt_int)(Xt_isign(global_size[num_dimensions - 1]) *
             isign(local_size[num_dimensions - 1]));
  idxsection->dims[num_dimensions - 1].local_stride = 1;

  // compute local and global stride
  // (local stride is always positive, global stride can be negative)
  for (int i = num_dimensions - 2; i >= 0; --i) {
    idxsection->dims[i].global_stride
      = (Xt_int)(idxsection->dims[i+1].global_stride * global_size[i + 1]
                 * Xt_isign((Xt_int)(idxsection->dims[i+1].global_stride
                                     * global_size[i + 1]))
                 * Xt_isign((Xt_int)(global_size[i]
                                     * (Xt_int)isign(local_size[i]))));
    idxsection->dims[i].local_stride
      = (Xt_int)(idxsection->dims[i + 1].local_stride * local_size[i + 1]
                 * Xt_isign((Xt_int)(idxsection->dims[i + 1].local_stride
                                     * (Xt_int)local_size[i + 1])));
  }

  // compute the local start index
  // depends on global size and sign of local and global size
  for (int i = num_dimensions - 1; i >= 0; --i) {
    if (global_size[i] > 0)
      idxsection->local_start_index
        = (Xt_int)(idxsection->local_start_index
                   + (XT_INT_ABS(idxsection->dims[i].global_stride)
                      * local_start[i]));
    else
      idxsection->local_start_index
        = (Xt_int)(idxsection->local_start_index
                   - (XT_INT_ABS(idxsection->dims[i].global_stride)
                      * (global_size[i] + local_start[i] + 1)));
    if (local_size[i] < 0)
      idxsection->local_start_index
        = (Xt_int)(idxsection->local_start_index
                   - (XT_INT_ABS(idxsection->dims[i].global_stride)
                      * (local_size[i] + 1)
                      * Xt_isign(global_size[i])));
  }

  // due the possibility of negative local and global sizes, the minimum and
  // maximum can be in any corner of the n-dimensional section
  idxsection->min_index_cache = idxsection->local_start_index;
  idxsection->max_index_cache = idxsection->local_start_index;
  for (int i = 0; i < num_dimensions; ++i) {

    // if either local and global size are negative
    if ((global_size[i] < 0) ^ (local_size[i] < 0))
      idxsection->min_index_cache
        = (Xt_int)(idxsection->min_index_cache
                   + (Xt_int)(idxsection->dims[i].global_stride *
                              (Xt_int)(abs(local_size[i]) - 1)));
    else // if local and global size are both positive or negative
      idxsection->max_index_cache
        = (Xt_int)(idxsection->max_index_cache
                   + (Xt_int)(idxsection->dims[i].global_stride *
                              (Xt_int)(abs(local_size[i]) - 1)));
  }

  INSTR_STOP(instr);

  return (Xt_idxlist)idxsection;
}

static void
idxsection_delete(Xt_idxlist data) {

  if (data == NULL) return;

  Xt_idxsection section = (Xt_idxsection)data;

  free(section->index_array_cache);
  free(section);
}

static size_t
idxsection_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

  Xt_idxsection section = (Xt_idxsection)data;

  int size_header, size_dim_descs, size_xt_int;

  xt_mpi_call(MPI_Pack_size(2, MPI_INT, comm, &size_header), comm);
  xt_mpi_call(MPI_Pack_size(2, Xt_int_dt, comm, &size_xt_int), comm);
  xt_mpi_call(MPI_Pack_size(section->ndim, dim_desc_dt,
                            comm, &size_dim_descs), comm);

  return (size_t)size_header + (size_t)size_dim_descs
    + (size_t)size_xt_int;
}

static void
idxsection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm) {

  INSTR_DEF(instr,"idxsection_pack")
  INSTR_START(instr);

  assert(data);
  Xt_idxsection section = (Xt_idxsection)data;
  int header[2] = { SECTION, section->ndim };
  Xt_int starts[2]
    = { section->global_start_index, section->local_start_index };
  xt_mpi_call(MPI_Pack(header, 2, MPI_INT, buffer,
                       buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(starts, 2, Xt_int_dt, buffer,
                       buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(section->dims, section->ndim, dim_desc_dt,
                       buffer, buffer_size, position, comm), comm);
  INSTR_STOP(instr);
}

Xt_idxlist xt_idxsection_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm) {

  INSTR_DEF(instr,"xt_idxsection_unpack")
  INSTR_START(instr);

  int ndim;
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position, &ndim, 1, MPI_INT,
                         comm), comm);
  Xt_idxsection section
    = xmalloc(sizeof (*section) + (size_t)ndim * sizeof(section->dims[0]));
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &section->global_start_index, 1, Xt_int_dt, comm),
              comm);
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &section->local_start_index, 1, Xt_int_dt, comm),
              comm);
  assert(ndim > 0);
  section->index_array_cache = NULL;
  section->ndim = ndim;

  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         section->dims, ndim, dim_desc_dt, comm),comm);

  // due to the possibility of negative local and global sizes, the minimum and
  // maximum can be in any corner of the n-dimensional section
  section->min_index_cache = section->local_start_index;
  section->max_index_cache = section->local_start_index;
  for (int i = 0; i < ndim; ++i) {

    // if either local or global size is negative
    if ((section->dims[i].global_size < 0) ^ (section->dims[i].local_size < 0))
      section->min_index_cache =
        (Xt_int)(section->min_index_cache + section->dims[i].global_stride
                 * (abs(section->dims[i].local_size) - 1));
    else // if local and global size are both positive or negative
      section->max_index_cache
        = (Xt_int)(section->max_index_cache
                   + (Xt_int)(section->dims[i].global_stride
                              * (Xt_int)(abs(section->dims[i].local_size)
                                         - 1)));
  }
  int num_indices = idxsection_get_num_indices(section);
  Xt_idxlist_init(&section->parent, &idxsection_vtable, num_indices);
  INSTR_STOP(instr);
  return (Xt_idxlist)section;
}

Xt_idxlist
xt_idxsection_get_intersection_with_other_idxlist(Xt_idxlist src_idxsection,
                                                  Xt_idxlist dst_idxlist,
                                                  Xt_config XT_UNUSED(config))
{
  // intersection between an idxsection and a general idxlist:
  //
  // performance picture:
  //  - src_idxsection is treated as too big for elemental transforms/access
  //  - dst_idxlist is considered to be small enough (subdomain like) for elemental usage

  INSTR_DEF(instr,"idxsection_get_intersection_with_other_idxlist")
  INSTR_START(instr);

  size_t num_dst_idx = (size_t)(xt_idxlist_get_num_indices(dst_idxlist));

  Xt_int const* dst_idx = xt_idxlist_get_indices_const(dst_idxlist);
  int single_match_only = 0;

  Xt_int const * sorted_dst_idx;
  Xt_int * temp_dst_idx = NULL;

  for (size_t i = 1; i < num_dst_idx; ++i)
    if (dst_idx[i] < dst_idx[i-1])
      goto unsorted;
  sorted_dst_idx = dst_idx;
  goto get_pos;
unsorted:
  temp_dst_idx = xmalloc(num_dst_idx * sizeof(*temp_dst_idx));
  memcpy(temp_dst_idx, dst_idx, num_dst_idx * sizeof(*temp_dst_idx));

  xt_mergesort_index(temp_dst_idx, (int)num_dst_idx, NULL, 0);
  sorted_dst_idx = temp_dst_idx;
get_pos:;
  int *pos = xmalloc(num_dst_idx * sizeof(*pos));
  size_t num_unmatched = idxsection_get_positions_of_indices(
    src_idxsection, sorted_dst_idx, num_dst_idx, pos,
    single_match_only);
  size_t num_inter_idx = num_dst_idx - num_unmatched;
  Xt_idxlist result;
  if (num_inter_idx != 0) {
    Xt_int *intersection = xmalloc(num_inter_idx * sizeof(*intersection));

    for(size_t i = 0, j = 0; i < num_dst_idx && j < num_inter_idx; i++) {
      intersection[j] = sorted_dst_idx[i];
      j += (pos[i] >= 0);
    }

    result = xt_idxvec_new(intersection, (int)num_inter_idx);
    free(intersection);
  } else {
    result = xt_idxempty_new();
  }

  free(temp_dst_idx);
  free(pos);

  INSTR_STOP(instr);
  return result;
  // return xt_idxvec_pinned_new(intersection, num_inter_idx);
}

Xt_idxlist
xt_idxsection_get_intersection(Xt_idxlist idxlist_src,
                               Xt_idxlist idxlist_dst,
                               Xt_config config) {
  INSTR_DEF(instr,"idxsection_get_intersection.part")

  // both lists are index section:

  Xt_idxsection idxsection_src = (Xt_idxsection)idxlist_src,
    idxsection_dst = (Xt_idxsection)idxlist_dst;

  if (idxsection_src->ndim != idxsection_dst->ndim ||
      idxsection_src->global_start_index != idxsection_dst->global_start_index)
    return xt_default_isect(idxlist_src, idxlist_dst, config);

  int i;

  // the size of first global dimension is irrelevant,
  // the others have to be identically
  for (i = 1; i < idxsection_src->ndim; ++i)
    if (XT_INT_ABS(idxsection_src->dims[i].global_size)
        != XT_INT_ABS(idxsection_dst->dims[i].global_size))
      return xt_idxsection_get_intersection_with_other_idxlist(
        idxlist_src, idxlist_dst, config);

  Xt_int *local_start, *global_size;
  int *local_size;

  INSTR_START(instr);

  // dimension information for the intersection
  local_start = xmalloc((size_t)idxsection_src->ndim * sizeof(*local_start));
  local_size = xmalloc((size_t)idxsection_src->ndim * sizeof(*local_size));
  global_size = xmalloc((size_t)idxsection_src->ndim * sizeof(*global_size));

  // indices in an intersection have to be sorted in ascending order. therefore,
  // local and global sizes of the intersection have to be positive

  for (i = 0; i < idxsection_src->ndim; ++i) {

    Xt_int src_start, src_end, dst_start, dst_end, local_end;

    // the start value is the minmum position in the current dimension (with positive
    // size)
    // in case the global size of src or dst is negative the start value has be be
    // adjusted accordingly

    if (idxsection_src->dims[i].global_size >= 0)
      src_start = idxsection_src->dims[i].local_start;
    else
      src_start = (Xt_int)(-idxsection_src->dims[i].global_size
                           - abs(idxsection_src->dims[i].local_size)
                           - idxsection_src->dims[i].local_start);

    if (idxsection_dst->dims[i].global_size >= 0)
      dst_start = idxsection_dst->dims[i].local_start;
    else
      dst_start = (Xt_int)(-idxsection_dst->dims[i].global_size
                           - abs(idxsection_dst->dims[i].local_size)
                           - idxsection_dst->dims[i].local_start);

    src_end = (Xt_int)(src_start
                       + (Xt_int)abs(idxsection_src->dims[i].local_size));
    dst_end = (Xt_int)(dst_start
                       + (Xt_int)abs(idxsection_dst->dims[i].local_size));

    local_start[i] = (src_start > dst_start)?src_start:dst_start;
    local_end = (src_end > dst_end)?dst_end:src_end;

    if (local_end <= local_start[i]) {
      free(global_size);
      free(local_size);
      free(local_start);
      INSTR_STOP(instr);
      return xt_idxempty_new();
    }

    local_size[i] = (int)(local_end - local_start[i]);
    global_size[i] = XT_INT_ABS(idxsection_src->dims[i].global_size);
  }

  Xt_idxlist intersection
    = xt_idxsection_new(idxsection_src->global_start_index,
                        idxsection_src->ndim, global_size,
                        local_size, local_start);

  free(global_size);
  free(local_size);
  free(local_start);

  INSTR_STOP(instr);
  return intersection;
}

static Xt_idxlist
idxsection_copy(Xt_idxlist idxlist) {

  Xt_idxsection src = (Xt_idxsection)idxlist;

  int num_dimensions = src->ndim;

  Xt_idxsection idxsection = xmalloc(sizeof (*idxsection)
                                     + (size_t)num_dimensions
                                     * sizeof (idxsection->dims[0]));
  *idxsection = *src;
  idxsection->index_array_cache  = NULL;

  memcpy(idxsection->dims, src->dims, (size_t)num_dimensions *
                                      sizeof (src->dims[0]));

  return (Xt_idxlist)idxsection;
}

static int
idxsection_get_num_indices(Xt_idxsection section) {

  int i;
  long long size = 1;

  for (i = 0; i < section->ndim; ++i)
    size *= abs(section->dims[i].local_size);
  assert(size <= INT_MAX);

  return (int)size;
}


static int
idxsection_get_indices_any(Xt_int start_index, Xt_int *indices,
                           int ndim, struct dim_desc dims[ndim])
{

  int abs_local_size = abs(dims[0].local_size);

  if (ndim == 1)
  {
    if (dims[0].global_stride > 0)
      for (int i = 0; i < abs_local_size; ++i)
        indices[i] = (Xt_int)(start_index + i);
    else
      for (int i = 0; i < abs_local_size; ++i)
        indices[i] = (Xt_int)(start_index - i);
    return abs_local_size;
  }
  else
  {
    int indices_written = 0, overflow = 0;
    assert(ndim > 1);
    for (int dim_ofs = 0; dim_ofs < abs_local_size; ++dim_ofs)
    {
      int indices_written_temp
        = idxsection_get_indices_any(
          (Xt_int)(start_index
                   + dim_ofs * dims[0].global_stride),
          indices + indices_written,
          ndim - 1, dims + 1);
      overflow |= (indices_written_temp > INT_MAX - indices_written);
      indices_written += indices_written_temp;
    }
    assert(!overflow);
    return indices_written;
  }
}

static void
idxsection_create_index_array_cache(Xt_idxsection section)
{
  size_t num_indices = (size_t)section->parent.num_indices;
  Xt_int *indices = section->index_array_cache
    = xmalloc(num_indices * sizeof(*(section->index_array_cache)));
  idxsection_get_indices_any(section->local_start_index, indices,
                             section->ndim, section->dims);
}


static void
idxsection_get_indices(Xt_idxlist idxlist, Xt_int *indices) {
  INSTR_DEF(instr,"idxsection_get_indices")
  INSTR_START(instr);
  Xt_idxsection section = (Xt_idxsection)idxlist;

  int num_indices = idxlist->num_indices;

  if (num_indices > 0) {
    // if the indices are already computed
    if (section->index_array_cache != NULL);
    else
      idxsection_create_index_array_cache(section);
    memcpy(indices, section->index_array_cache,
           (size_t)num_indices * sizeof(*indices));
  }
  INSTR_STOP(instr);
}

static Xt_int const*
idxsection_get_indices_const(Xt_idxlist idxlist) {

  Xt_idxsection idxsection = (Xt_idxsection)idxlist;
  if (idxsection->index_array_cache == NULL)
    idxsection_create_index_array_cache(idxsection);
  return idxsection->index_array_cache;
}


static void
idxsection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe **stripes,
                             int * num_stripes) {

  INSTR_DEF(instr,"idxsection_get_index_stripes.part")

  Xt_idxsection section = (Xt_idxsection)idxlist;

  int ndim = section->ndim;
  struct dim_desc *restrict dims = section->dims;

  size_t nstripes = dims[ndim-1].local_size != 0;

  for (int i = 0; i < ndim-1; ++i)
    nstripes *= (size_t)abs(dims[i].local_size);

  if (nstripes == 0) {
    *num_stripes = (int)nstripes;
    return;
  }

  INSTR_START(instr);

  struct Xt_stripe *restrict p = *stripes;
  if ((size_t)*num_stripes < nstripes)
    p = xrealloc(p, nstripes * sizeof(**stripes));

  enum { curr_local_position_auto_size=16 };
  Xt_int curr_local_position_auto[curr_local_position_auto_size];
  Xt_int *restrict curr_local_position;
  if (ndim-2 <= curr_local_position_auto_size) {
    curr_local_position = curr_local_position_auto;
    for (int i = 0; i < ndim-1; ++i)
      curr_local_position[i] = 0;
  } else
    curr_local_position
      = xcalloc((size_t)(ndim-2), sizeof(*curr_local_position));

  for (size_t i = 0; i < nstripes; ++i) {

    p[i].start  = section->local_start_index;
    p[i].nstrides = abs(dims[ndim-1].local_size);
    p[i].stride = 1;

    for (int j = 0; j < ndim - 1; ++j)
      p[i].start = (Xt_int)(p[i].start
                            + curr_local_position[j]
                            * dims[j].global_stride);

    for (int j = ndim - 2; j >= 0; --j)
      if (curr_local_position[j] < abs(dims[j].local_size) - 1) {
        curr_local_position[j]++;
        break;
      } else
        curr_local_position[j] = 0;
  }
  *stripes = p;
  *num_stripes = (int)nstripes;
  if (curr_local_position != curr_local_position_auto)
    free(curr_local_position);

  INSTR_STOP(instr);
}

static int
idxsection_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index) {

  Xt_idxsection section = (Xt_idxsection)idxlist;

  if (position < 0) return 1;

  Xt_int temp_index;

  temp_index = section->local_start_index;

  int dim;
  Xt_int curr_local_position;
  long long pos = (long long)position;

  for (dim = 0; dim < section->ndim; ++dim) {

    curr_local_position = (Xt_int)(pos / (long long)section->dims[dim].local_stride);

    if (curr_local_position >= abs(section->dims[dim].local_size))
      return 1;

    temp_index = (Xt_int)(temp_index
                          + curr_local_position
                          * section->dims[dim].global_stride);
    /* FIXME: assert(section->dims[dim].local_stride is in [-INT_MAX,INT_MAX]) */
    pos %= (long long)section->dims[dim].local_stride;
  }

  *index = temp_index;

  return 0;
}

static int
idxsection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position) {

  INSTR_DEF(instr,"idxsection_get_position_of_index.part")

  Xt_idxsection section = (Xt_idxsection)idxlist;
  *position = -1;

  if (index < section->min_index_cache || index > section->max_index_cache)
    return 1;

  int retval = 1;

  INSTR_START(instr);

  // normalise index (global start of indices at 0)
  index = (Xt_int)(index - section->global_start_index);

  int i;
  int temp_position = 0;

  for (i = 0; i < section->ndim; ++i) {

    Xt_int abs_global_stride
      = XT_INT_ABS(section->dims[i].global_stride);

    Xt_int curr_global_position
      = (Xt_int)(index / abs_global_stride);

    // in case the global size is negative, we have to adjust the global position,
    // because the ordering of indices in this dimension is inverted
    if (section->dims[i].global_size < 0)
      curr_global_position
        = (Xt_int)(-section->dims[i].global_size - curr_global_position - 1);

    index = (Xt_int)(index % abs_global_stride);

    if (curr_global_position < section->dims[i].local_start)
      goto fun_exit;

    Xt_int curr_local_position
      = (Xt_int)(curr_global_position - section->dims[i].local_start);

    // same adjustment for local position as for the global one before
    if (section->dims[i].local_size < 0)
      curr_local_position
        = (Xt_int)(-section->dims[i].local_size - curr_local_position - 1);

    if (curr_local_position >= abs(section->dims[i].local_size))
      goto fun_exit;

    temp_position += (int)(curr_local_position * section->dims[i].local_stride);
  }

  *position = temp_position;

  retval = 0;

 fun_exit: ;
  INSTR_STOP(instr);
  return retval;
}

static size_t
idxsection_get_positions_of_indices_v1(Xt_idxlist body_idxlist,
                                       const Xt_int selection_idx[],
                                       size_t num_selection,
                                       int positions[],
                                       int single_match_only) {

  INSTR_DEF(instr,"idxsection_get_positions_of_indices_v1.part")
  if (num_selection == 1)
    return (size_t)(idxsection_get_position_of_index(body_idxlist,
                                                     *selection_idx,
                                                     positions));

  size_t num_unmatched = 0;

  if (!single_match_only) {
    // this is the easy case, we don't care about multiple uses of the same position
    for (size_t i = 0; i < num_selection; ++i)
      num_unmatched
        += (size_t)(idxsection_get_position_of_index(body_idxlist,
                                                     selection_idx[i],
                                                     &positions[i]));
    return num_unmatched;
  }

  INSTR_START(instr);

  for (size_t i = 1; i < num_selection; ++i)
    if (selection_idx[i] < selection_idx[i-1])
      goto unsorted;

  // indices are sorted
  {
    // we need an index that is different from the current one
    Xt_int prev_index = (Xt_int)(selection_idx[0] - 1);

    for (size_t i = 0; i < num_selection; i++) {

      Xt_int curr_index = selection_idx[i];

      if (prev_index != curr_index) {

        num_unmatched
          += (size_t)(idxsection_get_position_of_index(body_idxlist, curr_index,
                                                       positions + i));
        prev_index = curr_index;

      } else {
        // for an idxsection there is a unique map from indices to positions,
        // we got the same index again, so there is no match left:
        positions[i] = -1;
        num_unmatched++;
      }
    }
  }
  goto end;
  // indices are not sorted
unsorted:
  {
    // the remaining (single_match_only) case follows:
    idxpos_type *v = xmalloc(num_selection * sizeof(*v) );
    for (size_t i = 0; i < num_selection; i++) {
      v[i].idx = selection_idx[i];
      v[i].pos = (int)i;
    }
    xt_mergesort_idxpos(v, num_selection);
    Xt_int last_jx = (Xt_int)(v[0].idx - 1); // any index that does not equal v[0].idx will do
    for (size_t i = 0; i < num_selection; i++) {
      int j = v[i].pos;
      Xt_int jx = v[i].idx;
      if (jx != last_jx) {
        num_unmatched
          += (size_t)(idxsection_get_position_of_index(body_idxlist, jx,
                                                       &positions[j]));
        last_jx = jx;
      } else {
        // for an idxsection there is a unique map from indices to positions,
        // we got the same index again, so there is no match left:
        positions[j] = -1;
        num_unmatched++;
      }
    }
    free(v);
  }
end:
  INSTR_STOP(instr);
  return num_unmatched;
}

static size_t
idxsection_get_positions_of_indices_v2(Xt_idxlist body_idxlist,
                                       const Xt_int selection_idx[],
                                       size_t num_selection, int positions[],
                                       int single_match_only) {

  INSTR_DEF(instr,"idxsection_get_positions_of_indices_v2.part")

  if (num_selection == 1)
    return (size_t)(idxsection_get_position_of_index(body_idxlist,
                                                     *selection_idx, positions));

  INSTR_START(instr);

  Xt_int * temp_selection_idx = NULL;
  const Xt_int *restrict sorted_selection_idx;
  int * selection_pos = NULL;

  for (size_t i = 1; i < num_selection; ++i)
    if (selection_idx[i] < selection_idx[i-1])
      goto unsorted;

  sorted_selection_idx = selection_idx;
  goto sorted;
unsorted:
  // the indices are not sorted
  temp_selection_idx
    = xmalloc(num_selection * sizeof(*temp_selection_idx));
  memcpy(temp_selection_idx, selection_idx,
         num_selection * sizeof(*temp_selection_idx));
  selection_pos = xmalloc(num_selection * sizeof(*selection_pos));

  xt_assign_id_map_int(num_selection, selection_pos, 0);
  xt_quicksort_xt_int_permutation(temp_selection_idx, num_selection,
                                  selection_pos);
  sorted_selection_idx = temp_selection_idx;

sorted: ;
  const Xt_int *body_indices = idxsection_get_indices_const(body_idxlist);
  size_t num_body_indices = (size_t)xt_idxlist_get_num_indices(body_idxlist);

  // Xt_int last_idx = sorted_selection_idx[0] - 1;
  // 
  // for (i = 0, j = 0; i < num_selection && j < num_body_indices; ++i) {
  // 
  //   while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;
  // 
  //   if (j >= num_body_indices) break;
  // 
  //   if (!single_match_only)
  //     positions[(selection_pos == NULL)?i:selection_pos[i]] =
  //       (body_indices[j] == sorted_selection_idx[i])?j:-1;
  //   else
  //     positions[selection_pos[i]] =
  //       ((last_idx == sorted_selection_idx[i]) ||
  //       (body_indices[j] != sorted_selection_idx[i]))?-1:j;
  // }

  // the following loops are an unrolled version of the one above
  size_t i = 0;
  if (!single_match_only) {

    if (selection_pos == NULL) {
      for (size_t j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[i] = (body_indices[j] == sorted_selection_idx[i])?(int)j:-1;
      }
    } else {
      for (size_t j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[selection_pos[i]] = (body_indices[j] == sorted_selection_idx[i])?(int)j:-1;
      }
    }
  } else {

    Xt_int last_idx = (Xt_int)(sorted_selection_idx[0] - 1);

    if (selection_pos == NULL) {
      for (size_t j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[i] = ((last_idx == sorted_selection_idx[i]) ||
                        (body_indices[j] != sorted_selection_idx[i]))?-1:(int)j;

        last_idx = sorted_selection_idx[i];
      }
    } else {
      for (size_t j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[selection_pos[i]] = ((last_idx == sorted_selection_idx[i]) ||
                                       (body_indices[j] != sorted_selection_idx[i]))?-1:(int)j;

        last_idx = sorted_selection_idx[i];
      }
    }
  }

  // process indices that were not handled by the loop above
  if (selection_pos == NULL)
    for (; i < num_selection; ++i)
      positions[i] = -1;
  else
    for (; i < num_selection; ++i)
      positions[selection_pos[i]] = -1;

  free(temp_selection_idx);
  free(selection_pos);

  size_t num_unmatched = 0;

  // count the number of unmachted indices
  for (size_t j = 0; j < num_selection; ++j)
    num_unmatched += positions[j] == -1;

  INSTR_STOP(instr);
  return num_unmatched;
}

static size_t
idxsection_get_positions_of_indices_recursive(Xt_int index_offset,
                                              int position_offset,
                                              const Xt_int indices[],
                                              size_t num_indices,
                                              int positions[],
                                              int ndim,
                                              struct dim_desc dims[ndim])
{
  size_t num_processed = 0;

  Xt_int abs_global_size = XT_INT_ABS(dims[0].global_size);
  int abs_local_size = abs(dims[0].local_size);
  Xt_int abs_global_stride = XT_INT_ABS(dims[0].global_stride);
  Xt_int abs_local_stride = XT_INT_ABS(dims[0].local_stride);

  if (ndim == 1)
  {

    Xt_int curr_position;

    Xt_int tmp_local_start = dims[0].local_start;

    // we want to work on ascending indices in the lowest dimension -> have to
    // adjust in case of negative global size
    if (dims[0].global_size < 0)
      tmp_local_start = (Xt_int)(-dims[0].global_size - tmp_local_start -
                                 abs_local_size);

    Xt_int min_index = (Xt_int)(index_offset + tmp_local_start * abs_global_stride);

    // set all indices that are smaller than the minimum to "not found"
    while ((num_processed < num_indices)
           && (indices[num_processed] < min_index))
      positions[num_processed++] = -1;

    // if either the local or the global dimension is negative
    if ((dims[0].global_stride < 0) ^ (dims[0].local_stride < 0)) {

      // for as long as we are in the range local section of the current
      // global dimension
      while ((num_processed < num_indices) &&
             ((curr_position = (Xt_int)(indices[num_processed] - min_index)) <
             abs_local_size)) {

        positions[num_processed++] = position_offset
          + (int)(abs_local_size - curr_position - 1);
      }
    } else { // if the local and global dimension are both negative or positive

      // for as long as we are in the range local section of the current
      // global dimension
      while ((num_processed < num_indices) &&
             ((curr_position = (Xt_int)(indices[num_processed] - min_index)) <
             abs_local_size)) {

        positions[num_processed++] = position_offset + (int)curr_position;
      }
    }

    // for all remaining indices that are in the current global dimension but not
    // within the local section
    while ((num_processed < num_indices) &&
           (indices[num_processed] < index_offset + abs_global_size))
      positions[num_processed++] = -1;

  } else {

    assert(ndim > 1);

    Xt_int tmp_local_start = dims[0].local_start;

    // we want to work on ascending indices in the lowest dimension -> have to
    // adjust in case of negative global size
    if (dims[0].global_size < 0)
      tmp_local_start = (Xt_int)(-dims[0].global_size - tmp_local_start -
                                 abs_local_size);

    Xt_int min_index
      = (Xt_int)(index_offset + tmp_local_start * abs_global_stride);

    // set all indices that are smaller than the minimum to "not found"
    while ((num_processed < num_indices)
           && (indices[num_processed] < min_index))
      positions[num_processed++] = -1;

    // while there are indices that have not yet been processed
    while (num_processed < num_indices) {

      Xt_int curr_global_position, curr_local_position;

      // compute global position of the smallest index that has not yet been processed
      curr_global_position = (Xt_int)((indices[num_processed] - index_offset) /
                                      abs_global_stride);

      // if the position is outside of the range of the current dimension
      if (curr_global_position >= tmp_local_start + abs_local_size)
        break;

      // if either the local or the global dimension is negative
      if ((dims[0].global_size < 0) ^ (dims[0].local_size < 0))

        curr_local_position = (Xt_int)(abs_local_size - curr_global_position
                                       + tmp_local_start - 1);
      else // if the local and global dimension are both negative or positive
        curr_local_position = (Xt_int)(curr_global_position - tmp_local_start);

      Xt_int curr_index_offset
        = (Xt_int)(index_offset + curr_global_position * abs_global_stride);
      /* FIXME: no guarantee curr_local_position * abs_local_stride
       * <= INT_MAX */
      int position_offset_ = (int)(curr_local_position * abs_local_stride);

      num_processed += idxsection_get_positions_of_indices_recursive(
        curr_index_offset, position_offset_, indices + num_processed,
        num_indices - num_processed, positions + num_processed, ndim-1,
        dims + 1);
    }
  }

  return num_processed;
}

static size_t
idxsection_get_positions_of_indices_v3(Xt_idxlist body_idxlist,
                                       const Xt_int *restrict selection_idx,
                                       size_t num_selection,
                                       int *restrict positions,
                                       int single_match_only) {

  INSTR_DEF(instr,"idxsection_get_positions_of_indices_v3.part")
  INSTR_DEF(instr2,"idxsection_get_positions_of_indices_recursive")

  Xt_idxsection section = (Xt_idxsection)body_idxlist;

  if (num_selection == 1)
    return (size_t)(idxsection_get_position_of_index(body_idxlist,
                                                     *selection_idx,
                                                     positions));

  INSTR_START(instr);

  const Xt_int * restrict sorted_selection_idx;
  Xt_int *temp_selection_idx = NULL;
  int *sorted_positions;
  int *selection_pos = NULL;

  for (size_t i = 1; i < num_selection; ++i)
    if (selection_idx[i] < selection_idx[i-1])
      goto unsorted_selection;

  sorted_selection_idx = selection_idx;
  sorted_positions = positions;
  goto sorted_selection;
  // if the selection is not sorted
unsorted_selection:
  temp_selection_idx
    = xmalloc(num_selection * sizeof(*temp_selection_idx));
  {
    size_t num_sp_alloc = num_selection;
#if defined _CRAYC && _RELEASE_MAJOR < 9
    num_sp_alloc = (num_sp_alloc + _MAXVL_32 - 1) & ~(_MAXVL_32 - 1);
#endif
    size_t total_alloc = num_sp_alloc + num_selection;
    sorted_positions
      = xmalloc(total_alloc * sizeof(*sorted_positions));
    selection_pos = sorted_positions + num_sp_alloc;
  }
  memcpy(temp_selection_idx, selection_idx,
         num_selection * sizeof(*temp_selection_idx));

  xt_assign_id_map_int(num_selection, selection_pos, 0);
  xt_quicksort_xt_int_permutation(temp_selection_idx, num_selection,
                                  selection_pos);
  sorted_selection_idx = temp_selection_idx;
sorted_selection:

  INSTR_START(instr2);

  size_t num_processed
    = idxsection_get_positions_of_indices_recursive(
      section->global_start_index,
      0, sorted_selection_idx, num_selection,
      sorted_positions, section->ndim,
      section->dims);

  INSTR_STOP(instr2);

  // set remaining index positions to -1
  for (size_t i = num_processed; i < num_selection; ++i)
    sorted_positions[i] = -1;

  // apply single match only rule
  if (single_match_only)
    for (size_t i = 1; i < num_processed; ++i)
      if (sorted_selection_idx[i] == sorted_selection_idx[i-1])
        sorted_positions[i] = -1;

  // convert positions if unsorted
  if (sorted_selection_idx != selection_idx) {

    for (size_t i = 0; i < num_selection; ++i)
      positions[i] = sorted_positions[selection_pos[i]];

    free(sorted_positions);
    free(temp_selection_idx);
  }

  // count the number of unmached indices
  size_t num_unmatched = num_selection - num_processed;

  for (size_t i = 0; i < num_processed; ++i)
    num_unmatched += positions[i] == -1;

  INSTR_STOP(instr);

  return num_unmatched;
}

static size_t
idxsection_get_positions_of_indices(Xt_idxlist body_idxlist,
                                    const Xt_int *selection_idx,
                                    size_t num_selection, int * positions,
                                    int single_match_only) {

  INSTR_DEF(instr,"idxsection_get_positions_of_indices")
  Xt_idxsection section = (Xt_idxsection)body_idxlist;
  size_t retval = 0;


  INSTR_START(instr);

  // if any dimension of the body index list is negative we have to use the
  // v3 version, because the other version cannot handle negative sizes
  for (int i = 0; i < section->ndim; ++i)
    if ((section->dims[i].local_size < 0) || (section->dims[i].global_size < 0)) {
      retval = idxsection_get_positions_of_indices_v3(body_idxlist, selection_idx,
                                                      num_selection, positions,
                                                      single_match_only);
      goto fun_exit;
    }

  {
    size_t num_section_indices
      = (size_t)(xt_idxlist_get_num_indices(body_idxlist));

    /*
     * if the indices are already cached or (if the caching would not
     * consume too much memory and the number of selection indices are
     * sufficient to justify the use of cached indices)
     */
    if ((section->index_array_cache != NULL) ||
        (((size_t)num_section_indices * sizeof(Xt_int)
          <= (size_t)128 * 1024U * 1024U)
         && (num_section_indices <= 1000 * num_selection))) {
      retval = idxsection_get_positions_of_indices_v2(body_idxlist, selection_idx,
                                                      num_selection, positions,
                                                      single_match_only);
      goto fun_exit;
    }
    else {
      retval = idxsection_get_positions_of_indices_v1(body_idxlist, selection_idx,
                                                      num_selection, positions,
                                                      single_match_only);
      goto fun_exit;
    }
  }

 fun_exit: ;
  INSTR_STOP(instr);
  return retval;
}

static int
idxsection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset) {

  int temp_position;
  // we make use of the uniqueness of the index-to-position relation:
  if (idxsection_get_position_of_index(idxlist, index, &temp_position))
    return 1;

  if (temp_position < offset)
    return 1;

  *position = temp_position;

  return 0;
}

static Xt_int
idxsection_get_min_index(Xt_idxlist idxlist) {

  Xt_idxsection section = (Xt_idxsection)idxlist;
  return section->min_index_cache;
}

static Xt_int
idxsection_get_max_index(Xt_idxlist idxlist) {

  Xt_idxsection section = (Xt_idxsection)idxlist;
  return section->max_index_cache;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
