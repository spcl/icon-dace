/**
 * @file xt_redist_repeat.c
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

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_repeat.h"
#include "xt/xt_redist_single_array_base.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"
#include "xt_config_internal.h"

static void
generate_msg_infos(struct Xt_redist_msg *restrict msgs,
                   MPI_Aint extent, const int *displacements, Xt_redist redist,
                   int num_repetitions, MPI_Comm comm,
                   enum xt_msg_direction direction)
{
  int *restrict ranks = NULL;
  size_t num_ranks
    = (size_t)xt_redist_get_msg_ranks(redist, direction, &ranks);
  for (size_t i = 0; i < num_ranks; ++i) {
    MPI_Datatype datatype
      = xt_redist_get_MPI_Datatype(redist, ranks[i], direction, false);
    MPI_Aint curr_lb, curr_extent;
    MPI_Datatype datatype_with_extent;

    // adjust extent of datatype to match the displacements
    xt_mpi_call(MPI_Type_get_extent(datatype, &curr_lb, &curr_extent), comm);
    xt_mpi_call(MPI_Type_create_resized(datatype, curr_lb, extent,
                                        &datatype_with_extent), comm);

    msgs[i].rank = ranks[i];
    msgs[i].datatype
      = xt_mpi_generate_datatype(displacements, num_repetitions,
                                 datatype_with_extent, comm);
    MPI_Type_free(&datatype_with_extent);
  }
  free(ranks);
}

Xt_redist
xt_redist_repeat_asym_new(Xt_redist redist, MPI_Aint src_extent,
                          MPI_Aint dst_extent, int num_repetitions,
                          const int src_displacements[num_repetitions],
                          const int dst_displacements[num_repetitions])
{
  return xt_redist_repeat_asym_custom_new(
    redist, src_extent, dst_extent, num_repetitions,
    src_displacements, dst_displacements, (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_repeat_asym_custom_new(Xt_redist redist, MPI_Aint src_extent,
                                 MPI_Aint dst_extent, int num_repetitions,
                                 const int src_displacements[num_repetitions],
                                 const int dst_displacements[num_repetitions],
                                 Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized());

  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_redist_get_MPI_Comm(redist), &tag_offset);

  if (num_repetitions < 1)
    Xt_abort(comm, "ERROR: invalid number of repetitions (Xt_redist_repeat)\n",
             __FILE__, __LINE__);

  int nmsg[2];
  for (int i = 0; i < 2; ++i)
    nmsg[i] = redist->vtable->get_num_msg(redist, (enum xt_msg_direction)i);
  size_t num_messages = (size_t)nmsg[SEND] + (size_t)nmsg[RECV];
  struct Xt_redist_msg *msgs = xmalloc(sizeof (*msgs) * num_messages);

  generate_msg_infos(msgs, src_extent, src_displacements,
                     redist, num_repetitions, comm, SEND);

  generate_msg_infos(msgs+nmsg[0], dst_extent, dst_displacements,
                     redist, num_repetitions, comm, RECV);

  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;

  Xt_redist result = xt_redist_single_array_base_custom_new(
    nmsg[SEND], nmsg[RECV], msgs, msgs+nmsg[SEND], comm, &config_);
  free(msgs);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

Xt_redist xt_redist_repeat_new(Xt_redist redist, MPI_Aint src_extent,
                               MPI_Aint dst_extent, int num_repetitions,
                               const int displacements[num_repetitions]) {
  return xt_redist_repeat_asym_custom_new(
    redist, src_extent, dst_extent, num_repetitions, displacements,
    displacements, (Xt_config)&xt_default_config);
}

Xt_redist xt_redist_repeat_custom_new(Xt_redist redist, MPI_Aint src_extent,
                                      MPI_Aint dst_extent, int num_repetitions,
                                      const int displacements[num_repetitions],
                                      Xt_config config)
{
  return xt_redist_repeat_asym_custom_new(
    redist, src_extent, dst_extent, num_repetitions, displacements,
    displacements, config);
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
