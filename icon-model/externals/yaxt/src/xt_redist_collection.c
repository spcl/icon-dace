/**
 * @file xt_redist_collection.c
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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_collection.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt/xt_request.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_config_internal.h"

enum { DEFFAULT_DATATYPE_CACHE_SIZE=16 };

static void
redist_collection_delete(Xt_redist redist);

static Xt_redist
redist_collection_copy(Xt_redist redist);

static void
redist_collection_s_exchange(Xt_redist redist, int num_src_arrays,
                             const void **src_data, void **dst_data);

static void
redist_collection_a_exchange(Xt_redist redist, int num_src_arrays,
                             const void **src_data, void **dst_data,
                             Xt_request *request);

static void
redist_collection_s_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data);

static void
redist_collection_a_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data,
                              Xt_request *request);

static int redist_collection_get_num_msg(Xt_redist redist,
                                         enum xt_msg_direction direction);

static MPI_Datatype
redist_collection_get_MPI_Datatype(Xt_redist redist, int rank,
                                   enum xt_msg_direction direction, bool do_dup);

static int
redist_collection_get_msg_ranks(Xt_redist redist,
                                enum xt_msg_direction direction,
                                int *restrict *ranks);

static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_collection_vtable = {
  .copy                  = redist_collection_copy,
  .delete                = redist_collection_delete,
  .s_exchange            = redist_collection_s_exchange,
  .a_exchange            = redist_collection_a_exchange,
  .s_exchange1           = redist_collection_s_exchange1,
  .a_exchange1           = redist_collection_a_exchange1,
  .get_num_msg           = redist_collection_get_num_msg,
  .get_msg_MPI_Datatype  = redist_collection_get_MPI_Datatype,
  .get_msg_ranks         = redist_collection_get_msg_ranks,
  .get_MPI_Comm          = redist_collection_get_MPI_Comm
};

struct exchanger_cache
{
  size_t token;
  /* array structure: [cache_size][2][num_redists] */
  MPI_Aint *displacements;
  Xt_exchanger * exchangers;
};

typedef struct Xt_redist_collection_ *Xt_redist_collection;

struct Xt_redist_collection_ {

  const struct xt_redist_vtable *vtable;

  struct exchanger_cache cache;

  unsigned num_redists, cache_size;

  unsigned nmsg[2];
  int *msg_ranks;

  struct Xt_config_ config;
  MPI_Comm comm;
  int tag_offset;

  MPI_Datatype all_component_dt[];
};

static void align_component_dt(unsigned num_redists, unsigned nmsgs,
                               const Xt_redist *redists,
                               int *restrict in_ranks[num_redists],
                               const size_t num_ranks[num_redists],
                               int *out_ranks,
                               MPI_Datatype *component_dt,
                               enum xt_msg_direction direction)
{
  size_t rank_pos[num_redists];
  for (size_t j = 0; j < num_redists; ++j)
    rank_pos[j] = 0;
  if (nmsgs) {
    /* find ranks and corresponding component datatypes */
    for (size_t i = 0; i < nmsgs; ++i) {
      int min_rank = INT_MAX;
      for (size_t j = 0; j < num_redists; ++j)
        if (rank_pos[j] < num_ranks[j] && in_ranks[j][rank_pos[j]] < min_rank)
          min_rank = in_ranks[j][rank_pos[j]];

      for (size_t j = 0; j < num_redists; ++j)
        component_dt[i * num_redists + j] =
          (rank_pos[j] < num_ranks[j] && in_ranks[j][rank_pos[j]] == min_rank)
          ? xt_redist_get_MPI_Datatype(redists[j], min_rank, direction, true)
          : MPI_DATATYPE_NULL;

      out_ranks[i] = min_rank;
      for (size_t j = 0; j < num_redists; ++j)
        rank_pos[j]
          += (rank_pos[j] < num_ranks[j] && in_ranks[j][rank_pos[j]] == min_rank);
    }
  }
  free(in_ranks[0]);
}

/* not yet used cache entries are marked with -1 as first displacement,
 * which becomes 0 later on through use */
static inline void
init_cache(struct exchanger_cache *cache, size_t cache_size,
           unsigned num_redists)
{
  cache->exchangers = xcalloc(cache_size, sizeof(*(cache->exchangers)));
  size_t num_displ = cache_size * num_redists;
  MPI_Aint *restrict q = cache->displacements
    = xmalloc(2 * num_displ * sizeof (*q));
  for (size_t i = 0; i < 2 * num_displ; i += num_redists)
    q[i] = (MPI_Aint)-1;
  cache->token = 0;
}

static inline void
destruct_cache(struct exchanger_cache *cache,
               size_t cache_size)
{
  for (size_t i = 0; i < cache_size; ++i)
    if (cache->exchangers[i] != NULL)
      xt_exchanger_delete(cache->exchangers[i]);
  free(cache->exchangers);
  free(cache->displacements);
}

Xt_redist xt_redist_collection_new(Xt_redist *redists, int num_redists,
                                   int cache_size, MPI_Comm comm)
{
  return xt_redist_collection_custom_new(redists, num_redists, cache_size,
                                         comm, (Xt_config)&xt_default_config);
}

static Xt_redist_collection
alloc_redist_coll(size_t num_redists, size_t nmsg_send, size_t nmsg_recv,
                  Xt_config config, MPI_Comm comm)
{
  size_t nmsg = nmsg_recv + nmsg_send,
    header_size = sizeof (struct Xt_redist_collection_),
    all_component_dt_size = sizeof (MPI_Datatype) * num_redists * nmsg,
    ranks_size = nmsg * sizeof (int),
    team_data_size,
    own_size = header_size + all_component_dt_size + ranks_size;
  Xt_exchanger_new exchanger_new = (Xt_exchanger_new)0;
  if (!config->exchanger_team_share) {
    exchanger_new
      = xt_config_get_exchange_new_by_comm(config, comm);
    team_data_size
      = xt_exchanger_new_team_get_share_size(exchanger_new);
    /* round up own size to align beginning of team_data */
    own_size = ((own_size + team_share_align - 1)
                / team_share_align) * team_share_align;
  } else
    team_data_size = 0;
  Xt_redist_collection redist_coll
    = xmalloc(own_size + team_data_size);
  redist_coll->vtable = &redist_collection_vtable;
  redist_coll->num_redists = (unsigned)num_redists;
  redist_coll->nmsg[SEND] = (unsigned)nmsg_send;
  redist_coll->nmsg[RECV] = (unsigned)nmsg_recv;
  redist_coll->msg_ranks
    = (int *)(redist_coll->all_component_dt + nmsg * num_redists);
  redist_coll->config = *config;
  /* all mpi datatypes are created for the exclusive use by the
   * underlying exchanger, no need for dup'ing them */
  redist_coll->config.flags |= exch_no_dt_dup;
  if (!config->exchanger_team_share) {
    redist_coll->config.exchanger_team_share
      = (unsigned char *)redist_coll + own_size;
    xt_exchanger_new_team_share_default_init(
      exchanger_new, redist_coll->config.exchanger_team_share);
  }
  return redist_coll;
}

Xt_redist xt_redist_collection_custom_new(Xt_redist *redists, int num_redists,
                                          int cache_size, MPI_Comm comm,
                                          Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized());
#ifndef NDEBUG
  if (num_redists > 0) ; else
    Xt_abort(comm, "ERROR: invalid number of redists passed",
             "xt_redist_collection_custom_new", __LINE__);
  if (cache_size >= -1) ; else
    Xt_abort(comm, "ERROR: invalid cache size in xt_redist_collection_new",
             __FILE__, __LINE__);
#endif
  unsigned num_redists_ = (unsigned)num_redists;
  size_t num_ranks[2][num_redists_];
  int *restrict ranks[2][num_redists_];
  unsigned nmsg_send = xt_redist_agg_msg_count(num_redists_, SEND, redists,
                                               num_ranks[SEND], ranks[SEND]),
    nmsg_recv = xt_redist_agg_msg_count(num_redists_, RECV, redists,
                                        num_ranks[RECV], ranks[RECV]);
  int tag_offset;
  MPI_Comm new_comm = xt_mpi_comm_smart_dup(comm, &tag_offset);
  Xt_redist_collection redist_coll
    = alloc_redist_coll(num_redists_, nmsg_send, nmsg_recv, config, new_comm);
  redist_coll->cache_size
    = (cache_size == -1)
    ?((unsigned)DEFFAULT_DATATYPE_CACHE_SIZE):(unsigned)cache_size;

  redist_coll->comm = new_comm;
  redist_coll->tag_offset = tag_offset;

  xt_redist_check_comms(redists, num_redists, comm);

  MPI_Datatype *all_component_dt = redist_coll->all_component_dt;
  align_component_dt(num_redists_, nmsg_send, redists,
                     ranks[SEND], num_ranks[SEND], redist_coll->msg_ranks,
                     all_component_dt, SEND);
  align_component_dt(num_redists_, nmsg_recv, redists,
                     ranks[RECV], num_ranks[RECV],
                     redist_coll->msg_ranks+nmsg_send,
                     all_component_dt + nmsg_send * num_redists_, RECV);
  init_cache(&redist_coll->cache, redist_coll->cache_size, num_redists_);

  return (Xt_redist)redist_coll;
}


static void
create_all_dt_for_dir(
  unsigned num_messages, size_t num_redists,
  const int ranks[num_messages],
  const MPI_Datatype *component_dt,
  const MPI_Aint displacements[num_redists],
  struct Xt_redist_msg redist_msgs[num_messages],
  MPI_Comm comm)
{
  int block_lengths[num_redists];

  for (size_t i = 0; i < num_redists; ++i)
    block_lengths[i] = 1;
  for (size_t i = 0; i < num_messages; ++i) {
    redist_msgs[i].datatype
      = xt_create_compound_datatype(num_redists, displacements,
                                    component_dt + i * num_redists,
                                    block_lengths, comm);
    redist_msgs[i].rank = ranks[i];
  }
}

static void
compute_displ(const void *const *data, unsigned num_redists,
              MPI_Aint displacements[num_redists])
{
  if (num_redists) {
    MPI_Aint base_addr, offset;
    base_addr = (MPI_Aint)(intptr_t)(const void *)data[0];
    displacements[0] = 0;
    for (size_t i = 1; i < num_redists; ++i) {
      offset = (MPI_Aint)(intptr_t)(const void *)data[i];
      displacements[i] = offset - base_addr;
    }
  }
}

static size_t
lookup_cache_index(unsigned num_redists,
                   const MPI_Aint displacements[2][num_redists],
                   const MPI_Aint (*cached_displacements)[2][num_redists],
                   size_t cache_size)
{
  for (size_t i = 0; i < cache_size &&
       cached_displacements[i][SEND][0] == (MPI_Aint)0 &&
       cached_displacements[i][RECV][0] == (MPI_Aint)0; ++i) {
    bool mismatch = false;
    for (size_t d = 0; d < 2; ++d)
      for (size_t j = 0; j < num_redists; ++j)
        mismatch |= (displacements[d][j] != cached_displacements[i][d][j]);
    if (!mismatch) return i;
  }
  return cache_size;
}

enum { redist_msg_stack_alloc_lim = 16 };

static Xt_exchanger
create_exchanger(struct Xt_redist_collection_ *redist_coll,
                 size_t num_redists,
                 MPI_Aint displacements[2][num_redists],
                 struct Xt_redist_msg *msgs)
{
  unsigned nmsg_send = redist_coll->nmsg[SEND],
    nmsg_recv = redist_coll->nmsg[RECV];
  const MPI_Datatype *all_component_dt = redist_coll->all_component_dt;
  MPI_Comm comm = redist_coll->comm;
  create_all_dt_for_dir(nmsg_send, num_redists,
                        redist_coll->msg_ranks, all_component_dt,
                        displacements[SEND], msgs, comm);
  create_all_dt_for_dir(nmsg_recv, num_redists,
                        redist_coll->msg_ranks+nmsg_send,
                        all_component_dt + nmsg_send * num_redists,
                        displacements[RECV],
                        msgs + nmsg_send, comm);
  Xt_exchanger_new exchanger_new
    = xt_config_get_exchange_new_by_comm(&redist_coll->config, comm);
  return exchanger_new((int)nmsg_send, (int)nmsg_recv,
                       msgs, msgs + nmsg_send, comm, redist_coll->tag_offset,
                       &redist_coll->config);
}


static Xt_exchanger
get_exchanger(struct Xt_redist_collection_ *redist_coll,
              const void *const * src_data, const void *const * dst_data,
              unsigned num_redists)
{
  MPI_Aint displacements[2][num_redists];
  unsigned nmsg_send = redist_coll->nmsg[SEND],
    nmsg_recv = redist_coll->nmsg[RECV];
  compute_displ(src_data, num_redists, displacements[SEND]);
  compute_displ(dst_data, num_redists, displacements[RECV]);

  Xt_exchanger exchanger;

  struct exchanger_cache *restrict cache = &redist_coll->cache;
  size_t cache_size = redist_coll->cache_size,
    cache_index = 0;
  if (cache_size > 0)
  {
    cache_index
      = lookup_cache_index(num_redists,
                           (const MPI_Aint (*)[num_redists])displacements,
                           (const MPI_Aint (*)[2][num_redists])
                           cache->displacements,
                           cache_size);

    if (cache_index != cache_size)
      return cache->exchangers[cache_index];

    cache_index = cache->token;
    cache->token = (cache_index + 1) % cache_size;
    memcpy(cache->displacements + cache_index * 2 * num_redists,
           displacements, sizeof (displacements));

    if (cache->exchangers[cache_index] != NULL)
      xt_exchanger_delete(cache->exchangers[cache_index]);
  }
  size_t nmsg = (size_t)nmsg_send + nmsg_recv;
  struct Xt_redist_msg *restrict p = xmalloc(nmsg * sizeof (*p));

  exchanger =
    create_exchanger(redist_coll, num_redists, displacements, p);
  free(p);
  if (cache_size > 0)
    cache->exchangers[cache_index] = exchanger;
  return exchanger;
}

static inline Xt_redist_collection
xrc(void *redist)
{
  return (Xt_redist_collection)redist;
}

static void
redist_collection_s_exchange(Xt_redist redist, int num_arrays,
                             const void **src_data, void **dst_data) {

  Xt_redist_collection redist_coll = xrc(redist);

  if (num_arrays != (int)redist_coll->num_redists)
    Xt_abort(redist_coll->comm, "ERROR: wrong number of arrays in "
             "redist_collection_s_exchange", __FILE__, __LINE__);


  Xt_exchanger exchanger = get_exchanger(redist_coll,
                                         src_data,
                                         (const void *const (*))dst_data,
                                         redist_coll->num_redists);

  xt_exchanger_s_exchange(exchanger, src_data[0], dst_data[0]);

  if (redist_coll->cache_size == 0)
    xt_exchanger_delete(exchanger);
}

static void
redist_collection_a_exchange(Xt_redist redist, int num_arrays,
                             const void **src_data, void **dst_data,
                             Xt_request *request) {

  Xt_redist_collection redist_coll = xrc(redist);

  if (num_arrays != (int)redist_coll->num_redists)
    Xt_abort(redist_coll->comm, "ERROR: wrong number of arrays in "
             "redist_collection_a_exchange", __FILE__, __LINE__);


  Xt_exchanger exchanger = get_exchanger(redist_coll,
                                         src_data,
                                         (const void *const (*))dst_data,
                                         redist_coll->num_redists);

  xt_exchanger_a_exchange(exchanger, src_data[0], dst_data[0], request);

  if (redist_coll->cache_size == 0)
    xt_exchanger_delete(exchanger);

}

static void
copy_component_dt(size_t num_component_dt,
                  const MPI_Datatype *component_dt_orig,
                  MPI_Datatype *component_dt_copy,
                  MPI_Comm comm)
{
  for (size_t i = 0; i < num_component_dt; ++i)
  {
    MPI_Datatype orig_dt = component_dt_orig[i];
    if (orig_dt != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_dup(orig_dt, component_dt_copy + i), comm);
    else
      component_dt_copy[i] = orig_dt;
  }
}

static Xt_redist
redist_collection_custom_copy(Xt_redist redist, Xt_config config);

static Xt_redist
redist_collection_copy(Xt_redist redist)
{
  Xt_redist_collection redist_coll = xrc(redist);
  unsigned nmsg_send = redist_coll->nmsg[SEND],
    nmsg_recv = redist_coll->nmsg[RECV];
  size_t nmsg = (size_t)nmsg_recv + nmsg_send;
  struct Xt_config_ config = redist_coll->config;
  /* if redist comes with an integrated team share, cause creation of
   * a new one */
  void *team_share = config.exchanger_team_share;
  if ((unsigned char *)team_share > (unsigned char *)redist_coll
      && (unsigned char *)team_share <
      (unsigned char *)(redist_coll->msg_ranks + nmsg) + team_share_align)
    config.exchanger_team_share = NULL;
  return redist_collection_custom_copy(redist, &config);
}

static Xt_redist
redist_collection_custom_copy(Xt_redist redist, Xt_config config)
{
  Xt_redist_collection redist_coll = xrc(redist);
  unsigned num_redists = redist_coll->num_redists,
    nmsg_send = redist_coll->nmsg[SEND],
    nmsg_recv = redist_coll->nmsg[RECV];
  size_t nmsg = (size_t)nmsg_recv + nmsg_send;
  int tag_offset;
  MPI_Comm copy_comm
    = xt_mpi_comm_smart_dup(redist_coll->comm, &tag_offset);
  Xt_redist_collection redist_copy
    = alloc_redist_coll(num_redists, nmsg_send, nmsg_recv, config, copy_comm);
  redist_copy->comm = copy_comm;
  redist_copy->tag_offset = tag_offset;

  memcpy(redist_copy->msg_ranks, redist_coll->msg_ranks,
         sizeof (*redist_copy->msg_ranks) * nmsg);
  copy_component_dt(num_redists * nmsg,
                    redist_coll->all_component_dt,
                    redist_copy->all_component_dt, copy_comm);
  unsigned cache_size = redist_coll->cache_size;
  redist_copy->cache_size = cache_size;
  init_cache(&redist_copy->cache, cache_size, num_redists);
  return (Xt_redist)redist_copy;
}

static void
free_component_dt(size_t num_dt, MPI_Datatype *all_component_dt, MPI_Comm comm)
{
  for (size_t i = 0; i < num_dt; ++i)
    if (all_component_dt[i] != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_free(all_component_dt + i), comm);
}

static void
redist_collection_delete(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  unsigned num_redists = redist_coll->num_redists;
  size_t nmsg = (size_t)redist_coll->nmsg[RECV] + redist_coll->nmsg[SEND];
  free_component_dt(nmsg * num_redists, redist_coll->all_component_dt,
                    redist_coll->comm);

  destruct_cache(&redist_coll->cache, redist_coll->cache_size);

  void *team_share = redist_coll->config.exchanger_team_share;
  if ((unsigned char *)team_share > (unsigned char *)redist_coll
      && (unsigned char *)team_share
      < (unsigned char *)(redist_coll->msg_ranks + nmsg) + team_share_align) {
    Xt_exchanger_new exchanger_new
      = xt_config_get_exchange_new_by_comm(&redist_coll->config,
                                           redist_coll->comm);
    xt_exchanger_new_team_share_destroy(exchanger_new, team_share);
  }
  xt_mpi_comm_smart_dedup(&(redist_coll->comm), redist_coll->tag_offset);

  free(redist_coll);
}

static int redist_collection_get_num_msg(Xt_redist redist,
                                         enum xt_msg_direction direction)
{
  return (int)(xrc(redist)->nmsg[direction]);
}

static MPI_Datatype
redist_collection_get_MPI_Datatype(Xt_redist redist, int XT_UNUSED(rank),
                                   enum xt_msg_direction XT_UNUSED(direction),
                                   bool XT_UNUSED(do_dup))
{
  Xt_redist_collection redist_coll = xrc(redist);

  Xt_abort(redist_coll->comm, "ERROR: datatype retrieval is not"
           " supported for this xt_redist type (Xt_redist_collection)",
           __FILE__, __LINE__);

  return MPI_DATATYPE_NULL;
}

static void
redist_collection_s_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data)
{

  Xt_redist_collection redist_coll = xrc(redist);
  if (redist_coll->num_redists == 1)
    redist_collection_s_exchange(redist, 1, &src_data, &dst_data);
  else
    Xt_abort(redist_coll->comm, "ERROR: s_exchange1 is not implemented for"
             " this xt_redist type (Xt_redist_collection)", __FILE__, __LINE__);
}

static void
redist_collection_a_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data,
                              Xt_request *request)
{

  Xt_redist_collection redist_coll = xrc(redist);
  if (redist_coll->num_redists == 1)
    redist_collection_a_exchange(redist, 1, &src_data, &dst_data, request);
  else
    Xt_abort(redist_coll->comm, "ERROR: a_exchange1 is not implemented for"
             " this xt_redist type (Xt_redist_collection)", __FILE__, __LINE__);
}

static int
redist_collection_get_msg_ranks(Xt_redist redist,
                                enum xt_msg_direction direction,
                                int *restrict *ranks)
{
  Xt_redist_collection redist_coll = xrc(redist);
  unsigned nmsg_direction = redist_coll->nmsg[direction],
    nmsg_send = redist_coll->nmsg[SEND];
  size_t nmsg = (size_t)nmsg_direction + redist_coll->nmsg[!direction];
  int *ranks_orig
    = (int *)(redist_coll->all_component_dt + nmsg * redist_coll->num_redists)
    + (((unsigned)direction-1) & nmsg_send);
  if (nmsg_direction > 0) {
    int *ranks_ = *ranks;
    if (!ranks_)
      ranks_ = *ranks = xmalloc(nmsg_direction * sizeof (*ranks_));
    memcpy(ranks_, ranks_orig, nmsg_direction * sizeof (*ranks));
  }
  return (int)nmsg_direction;
}


static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  return redist_coll->comm;
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
