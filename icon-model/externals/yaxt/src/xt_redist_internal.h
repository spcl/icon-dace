/**
 * @file xt_redist_internal.h
 * @brief redistribution of data, non-public declarations
 *
 * contains declaration the redistribution data structure, which
 * is derived from one or more xt_xmaps
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
#ifndef XT_REDIST_INTERNAL_H
#define XT_REDIST_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <stdlib.h>
#include <mpi.h>

#include "core/ppm_visibility.h"
#include "xt/xt_redist.h"
#include "xt/xt_request.h"

enum xt_msg_direction {SEND, RECV};

struct xt_redist_vtable {

  Xt_redist (*copy)(Xt_redist);
  void (*delete)(Xt_redist);
  void (*s_exchange)(Xt_redist, int, const void **, void **);
  void (*a_exchange)(Xt_redist, int, const void **, void **, Xt_request *);
  void (*s_exchange1)(Xt_redist, const void *, void *);
  void (*a_exchange1)(Xt_redist, const void *, void *, Xt_request *);
  MPI_Datatype (*get_msg_MPI_Datatype)(Xt_redist, int, enum xt_msg_direction,
                                       bool need_dup);
  int (*get_num_msg)(Xt_redist, enum xt_msg_direction);
  int (*get_msg_ranks)(Xt_redist, enum xt_msg_direction, int *restrict *);
  MPI_Comm (*get_MPI_Comm)(Xt_redist);
};

struct Xt_redist_ {
  const struct xt_redist_vtable *vtable;
};

PPM_DSO_INTERNAL void
xt_redist_msgs_strided_copy(size_t n,
                            const struct Xt_redist_msg *restrict src,
                            size_t src_stride,
                            struct Xt_redist_msg *restrict dst,
                            size_t dst_stride,
                            MPI_Comm comm, bool dt_dup);

PPM_DSO_INTERNAL void
xt_redist_msgs_strided_destruct(size_t n, struct Xt_redist_msg *msgs,
                                MPI_Comm comm, size_t ofs_stride);

static inline void
xt_redist_msgs_free(size_t n, struct Xt_redist_msg *msgs, MPI_Comm comm)
{
  xt_redist_msgs_strided_destruct(n, msgs, comm, sizeof (*msgs));
  free(msgs);
}

/**
 * Checks whether a number of redists are based on the same communicator. This
 * is a requirement in case these redists are to be combined into a redist
 * collection.
 * @param[in] redists     redistribution objects to be checked
 * @param[in] num_redists number of redistribution objects is redists
 * @param[in] comm        reference communicator
 * @remark In case this check fails, it will abort the program.
 */
PPM_DSO_INTERNAL void
xt_redist_check_comms(Xt_redist *redists, int num_redists, MPI_Comm comm);

/**
 * Gets the ranks of all processes that receive data from/send data to the local
 * process in the exchanges defined by the redist.
 * @param[in]  redist    redistribution object
 * @param[in]  direction specifices whether ranks for the outgoing or incoming
 *                       messages are requested
 * @param[in,out] ranks  ranks for all outgoing/incoming messages
 * @return number of outgoing/incoming message
 * @remark the user needs to ensure that array ranks is big enough to hold all
 *         ranks, each element of ranks must be either a pointer to a valid
 *         output array or NULL, in which case it will be allocated
 */
PPM_DSO_INTERNAL int
xt_redist_get_msg_ranks(Xt_redist redist, enum xt_msg_direction direction,
                        int *restrict *ranks);

/**
 * Gets a MPI derived datatype that encodes all data sent/received in a
 * specified message.
 * @param[in] redist    redistribution object
 * @param[in] rank      MPI rank of the communicator partner
 * @param[in] direction specifices whether the datatype for an outgoing or
 *                      incoming message is requested
 * @param[in] do_dup    if true only return MPI_Datatype_dup of
 *                      internally stored datatype
 * @return Datatype for the specified message. The return value is
 *         MPI_DATATYPE_NULL, if no data for the specified message.
 */
PPM_DSO_INTERNAL MPI_Datatype
xt_redist_get_MPI_Datatype(Xt_redist redist, int rank,
                           enum xt_msg_direction direction, bool do_dup);

/**
 * Generates a new MPI derived datatype from a number of MPI derived datatypes.
 * @param[in] count         number of datatypes
 * @param[in] displacements byte displacement of each block
 * @param[in] datatypes     type of elements in each block
 * @param[in] block_lengths number of elements in each block
 * @param[in] comm          communicator
 * @return new MPI derived datatype
 */
PPM_DSO_INTERNAL MPI_Datatype
xt_create_compound_datatype(size_t count,
                            const MPI_Aint displacements[count],
                            const MPI_Datatype datatypes[count],
                            const int block_lengths[count],
                            MPI_Comm comm);

/**
 * Determines number of processes that receive data from/send data to the local
 * processes by any redists provided to this routine.
 * @param[in]  num_redists number of redistribution objects
 * @param[in]  direction   specifices whether the of outgoing or incoming
 *                         message is to be determined
 * @param[in]  redists     redistribution objects
 * @param[out] num_ranks   number of incoming/outgoing messages per redist
 * @param[out] ranks       ranks of communicator partners for each redist
 * @return Number of processes that receive data from/send data to the
 *         local process.
 * @remark the user needs to ensure that rank arrays are big enough to hold all
 *         ranks
 */
PPM_DSO_INTERNAL unsigned
xt_redist_agg_msg_count(size_t num_redists, enum xt_msg_direction direction,
                        const Xt_redist redists[num_redists],
                        size_t num_ranks[num_redists],
                        int *restrict ranks[num_redists]);

#endif
/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
