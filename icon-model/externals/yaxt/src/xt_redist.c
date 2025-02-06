/**
 * @file xt_redist.c
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

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>

#include "core/core.h"
#include "xt/xt_core.h"
#include "xt/xt_redist.h"
#include "xt/xt_mpi.h"
#include "xt/xt_request.h"
#include "xt/xt_sort.h"
#include "core/ppm_xfuncs.h"
#include "xt_redist_internal.h"

Xt_redist xt_redist_copy(Xt_redist redist) {

  return redist->vtable->copy(redist);
}

void xt_redist_delete(Xt_redist redist) {

  redist->vtable->delete(redist);
}

void xt_redist_s_exchange(Xt_redist redist, int num_arrays,
                          const void **src_data, void **dst_data) {

  redist->vtable->s_exchange(redist, num_arrays, src_data, dst_data);
}

void xt_redist_a_exchange(Xt_redist redist, int num_arrays,
                          const void **src_data, void **dst_data,
                          Xt_request *request) {

  redist->vtable->a_exchange(redist, num_arrays, src_data, dst_data, request);
}

void xt_redist_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data) {

  redist->vtable->s_exchange1(redist, src_data, dst_data);
}

void xt_redist_a_exchange1(Xt_redist redist, const void *src_data,
                           void *dst_data, Xt_request *request) {

  redist->vtable->a_exchange1(redist, src_data, dst_data, request);
}

int xt_redist_get_num_send_msg(Xt_redist redist) {

  return redist->vtable->get_num_msg(redist, SEND);
}

int xt_redist_get_num_recv_msg(Xt_redist redist) {

  return redist->vtable->get_num_msg(redist, RECV);
}

MPI_Datatype xt_redist_get_send_MPI_Datatype(Xt_redist redist, int rank) {

  return redist->vtable->get_msg_MPI_Datatype(redist, rank, SEND, true);
}

MPI_Datatype xt_redist_get_recv_MPI_Datatype(Xt_redist redist, int rank) {

  return redist->vtable->get_msg_MPI_Datatype(redist, rank, RECV, true);
}

MPI_Datatype xt_redist_get_MPI_Datatype(Xt_redist redist, int rank,
                                        enum xt_msg_direction direction,
                                        bool do_dup)
{
  return redist->vtable->get_msg_MPI_Datatype(redist, rank, direction, do_dup);
}

MPI_Comm xt_redist_get_MPI_Comm(Xt_redist redist) {

  return redist->vtable->get_MPI_Comm(redist);
}

int xt_redist_get_msg_ranks(Xt_redist redist, enum xt_msg_direction direction,
                            int *restrict *ranks)
{
  return redist->vtable->get_msg_ranks(redist, direction, ranks);
}


void
xt_redist_check_comms(Xt_redist *redists, int num_redists, MPI_Comm comm) {
  int result;

  for (int i = 0; i < num_redists; ++i) {

    if (redists[i] == NULL)
      Xt_abort(comm, "ERROR: invalid redist; cannot build "
               "redistribution collection\n", __FILE__, __LINE__);

    xt_mpi_call(MPI_Comm_compare(xt_redist_get_MPI_Comm(redists[i]),
                                 comm, &result), comm);

    if ((result != MPI_IDENT) && (result != MPI_CONGRUENT))
      Xt_abort(comm, "ERROR: MPI communicators do not match; cannot build "
               "redistribution collection\n", __FILE__, __LINE__);
  }
}

static size_t
xt_ranks_uniq_count(size_t num_rank_sets,
                    const size_t *restrict num_ranks,
                    const int *const ranks[num_rank_sets])
{
  size_t rank_pos[num_rank_sets];
  for (size_t j = 0; j < num_rank_sets; ++j)
    rank_pos[j] = 0;
  bool ranks_left;
  size_t num_messages = 0;
  do {
    int min_rank = INT_MAX;
    /* find minimal rank in list, guaranteed to be smaller than comm_size */
    for (size_t j = 0; j < num_rank_sets; ++j)
      if (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] < min_rank)
        min_rank = ranks[j][rank_pos[j]];
    ranks_left = false;
    /* increment list index for all redists matching minimal rank and
     * see if any ranks are left */
    for (size_t j = 0; j < num_rank_sets; ++j) {
      rank_pos[j]
        += (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] == min_rank);
      ranks_left |= (rank_pos[j] < num_ranks[j]);
    }
    ++num_messages;
  } while (ranks_left);
  return num_messages;
}

unsigned
xt_redist_agg_msg_count(size_t num_redists, enum xt_msg_direction direction,
                        const Xt_redist redists[num_redists],
                        size_t num_ranks[num_redists],
                        int *restrict ranks[num_redists])
{
  bool ranks_left = false;
  /* get lists of ranks to send/receive message to/from */
  size_t num_ranks_total = 0;
  for (size_t j = 0; j < num_redists; ++j) {
    size_t redist_num_ranks
      = (size_t)(redists[j]->vtable->get_num_msg(redists[j], direction));
    num_ranks[j] = redist_num_ranks;
    num_ranks_total += redist_num_ranks;
  }
  if (num_ranks_total) {
    int *ranks_buf = xmalloc(num_ranks_total * sizeof (*ranks_buf));
    size_t ofs = 0;
    for (size_t j = 0; j < num_redists; ++j) {
      ranks[j] = ranks_buf + ofs;
      size_t nranks
        = (size_t)xt_redist_get_msg_ranks(redists[j], direction, ranks + j);
      ranks_left |= (nranks > 0);
      /* sort list */
      xt_sort_int(ranks[j], nranks);
      ofs += nranks;
    }
  } else
    for (size_t j = 0; j < num_redists; ++j)
      ranks[j] = NULL;
  /* count number of different ranks to send/receive message to/from */
  size_t num_messages = ranks_left
    ? xt_ranks_uniq_count(num_redists, num_ranks, (const int *const *)ranks)
    : 0;
  return (unsigned)num_messages;
}

MPI_Datatype
xt_create_compound_datatype(size_t count,
                            const MPI_Aint displacements[count],
                            const MPI_Datatype datatypes[count],
                            const int block_lengths[count],
                            MPI_Comm comm)
{
  size_t num_datatypes = 0;
  /* allocate more than max_auto_dt datatype items from heap */
  enum { max_auto_dt = 8 };
  for (size_t i = 0; i < count; ++i)
    num_datatypes += (datatypes[i] != MPI_DATATYPE_NULL);
  MPI_Datatype *datatypes_, dt_auto[max_auto_dt];
  MPI_Aint *displacements_, disp_auto[max_auto_dt];
  int *block_lengths_, bl_auto[max_auto_dt];

  if (num_datatypes != count) {
    if (num_datatypes > max_auto_dt) {
      size_t buf_size = num_datatypes * sizeof(*datatypes_)
        + num_datatypes * sizeof(*displacements_)
        + num_datatypes * sizeof(*block_lengths_);
      displacements_ = xmalloc(buf_size);
      datatypes_ = (MPI_Datatype *)(displacements_ + num_datatypes);
      block_lengths_ = (int *)(datatypes_ + num_datatypes);
    } else {
      datatypes_ = dt_auto;
      displacements_ = disp_auto;
      block_lengths_ = bl_auto;
    }
    num_datatypes = 0;

    for (size_t i = 0; i < count; ++i) {
      if (datatypes[i] != MPI_DATATYPE_NULL) {

        datatypes_[num_datatypes] = datatypes[i];
        displacements_[num_datatypes] = displacements[i];
        block_lengths_[num_datatypes] = block_lengths[i];
        ++num_datatypes;
      }
    }
  } else {
    datatypes_ = (MPI_Datatype *)datatypes;
    displacements_ = (MPI_Aint *)displacements;
    block_lengths_ = (int *)block_lengths;
  }
  MPI_Datatype datatype;
  if (num_datatypes > 1)
    xt_mpi_call(MPI_Type_create_struct((int)num_datatypes, block_lengths_,
                                       displacements_, datatypes_, &datatype),
                comm);
  else if (displacements_[0] == 0)
    xt_mpi_call(MPI_Type_dup(datatypes_[0], &datatype), comm);
  else
    xt_mpi_call(MPI_Type_create_hindexed(1, (int [1]){1}, displacements_,
                                         datatypes_[0], &datatype), comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  if (num_datatypes != count && num_datatypes > max_auto_dt)
    free(displacements_);

  return datatype;
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
