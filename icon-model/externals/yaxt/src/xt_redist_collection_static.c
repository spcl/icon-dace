/**
 * @file xt_redist_collection_static.c
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

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_collection_static.h"
#include "xt/xt_redist_single_array_base.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"
#include "xt_config_internal.h"


static void
generate_msg_infos(size_t nmsg, size_t num_redists,
                   struct Xt_redist_msg *msgs,
                   const MPI_Aint displacements[num_redists],
                   const Xt_redist redists[num_redists],
                   const size_t num_ranks[num_redists],
                   const int *restrict ranks[num_redists],
                   MPI_Comm comm,
                   enum xt_msg_direction direction)
{
  if (nmsg) {
    size_t rank_pos[num_redists];
    for (size_t j = 0; j < num_redists; ++j)
      rank_pos[j] = 0;
    MPI_Datatype datatypes[num_redists];
    int block_lengths[num_redists];
    for (size_t i = 0; i < num_redists; ++i)
      block_lengths[i] = 1;
    for (size_t i = 0; i < nmsg; ++i) {
      int min_rank = INT_MAX;
      for (size_t j = 0; j < num_redists; ++j)
        if (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] < min_rank)
          min_rank = ranks[j][rank_pos[j]];

      for (size_t j = 0; j < num_redists; ++j)
        datatypes[j] =
          (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] == min_rank)
          ? xt_redist_get_MPI_Datatype(redists[j], min_rank, direction, false)
          : MPI_DATATYPE_NULL;

      msgs[i].rank = min_rank;
      msgs[i].datatype
        = xt_create_compound_datatype(num_redists, displacements, datatypes,
                                      block_lengths, comm);
      for (size_t j = 0; j < num_redists; ++j) {
        rank_pos[j]
          += (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] == min_rank);
      }
    }
  }
}

Xt_redist
xt_redist_collection_static_new(Xt_redist *redists, int num_redists,
                                const MPI_Aint src_displacements[num_redists],
                                const MPI_Aint dst_displacements[num_redists],
                                MPI_Comm comm)
{
  return xt_redist_collection_static_custom_new(
    redists, num_redists, src_displacements, dst_displacements, comm,
    (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_collection_static_custom_new(
  Xt_redist * redists, int num_redists,
  const MPI_Aint src_displacements[num_redists],
  const MPI_Aint dst_displacements[num_redists],
  MPI_Comm comm, Xt_config config) {
  // ensure that yaxt is initialized
  assert(xt_initialized());

  int tag_offset;
  MPI_Comm new_comm = xt_mpi_comm_smart_dup(comm, &tag_offset);

  xt_redist_check_comms(redists, num_redists, comm);

  size_t num_redists_ = num_redists >= 0 ? (size_t)num_redists : 0;
  int *restrict ranks[2][num_redists_];
  size_t num_ranks[2][num_redists_];
  size_t nmsg[2];
  /* get lists of ranks to send/receive message to/from */
  for (size_t i = 0; i < 2; ++i)
    nmsg[i] = xt_redist_agg_msg_count(num_redists_, (enum xt_msg_direction)i,
                                      redists, num_ranks[i], ranks[i]);
  size_t nmsg_sum = nmsg[SEND] + nmsg[RECV];
  struct Xt_redist_msg *msgs = xmalloc(sizeof (*msgs) * nmsg_sum);
  for (size_t i = 0; i < 2; ++i) {
    size_t ofs = i == 0 ? 0 : nmsg[SEND];
    const MPI_Aint *disp = i == 0 ? src_displacements : dst_displacements;
    generate_msg_infos(nmsg[i], num_redists_, msgs+ofs, disp, redists,
                       num_ranks[i], (const int *restrict (*))ranks[i],
                       new_comm, (enum xt_msg_direction)i);
    free(ranks[i][0]);
  }

  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;

  Xt_redist redist_collection =
    xt_redist_single_array_base_custom_new(
      (int)nmsg[SEND], (int)nmsg[RECV], msgs, msgs+nmsg[SEND], new_comm,
      &config_);

  free(msgs);
  xt_mpi_comm_smart_dedup(&new_comm, tag_offset);
  return redist_collection;
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
