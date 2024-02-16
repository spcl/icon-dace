/**
 * @file xt_exchanger_irecv_isend.c
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

#include "core/ppm_xfuncs.h"
#include "xt_config_internal.h"
#include "xt/xt_mpi.h"
#include "xt/xt_request_msgs.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_simple_base.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static MPI_Request *
redist_msgs_to_req(const void *src_data, void *dst_data,
                   int nsend, int nrecv,
                   const struct Xt_redist_msg *send_msgs,
                   const struct Xt_redist_msg *recv_msgs,
                   int tag_offset, MPI_Comm comm)
{
  MPI_Request *requests
    = xmalloc((size_t)(nrecv + nsend) * sizeof (*requests));

  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Irecv(dst_data, 1, recv_msgs[i].datatype,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+i), comm);

  for (int i = 0; i < nsend; ++i)
    xt_mpi_call(MPI_Isend(CAST_MPI_SEND_BUF(src_data), 1, send_msgs[i].datatype,
                          send_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+nrecv+i), comm);
  return requests;
}

static void
xt_exchanger_irecv_isend_s_exchange(const void *src_data, void *dst_data,
                                    int nsend, int nrecv,
                                    const struct Xt_redist_msg *send_msgs,
                                    const struct Xt_redist_msg *recv_msgs,
                                    int tag_offset, MPI_Comm comm) {

  MPI_Request *requests = redist_msgs_to_req(src_data, dst_data, nsend, nrecv,
                                             send_msgs, recv_msgs, tag_offset,
                                             comm);
  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE), comm);
  free(requests);
}

static void
xt_exchanger_irecv_isend_a_exchange(const void *src_data, void *dst_data,
                                    int nsend, int nrecv,
                                    const struct Xt_redist_msg * send_msgs,
                                    const struct Xt_redist_msg * recv_msgs,
                                    int tag_offset, MPI_Comm comm,
                                    Xt_request *request) {

  MPI_Request *tmp_requests
    = redist_msgs_to_req(src_data, dst_data, nsend, nrecv,
                         send_msgs, recv_msgs, tag_offset, comm);
  Xt_request requests = xt_request_msgs_new(nrecv + nsend, tmp_requests, comm);

  free(tmp_requests);

  *request = requests;
}

Xt_exchanger
xt_exchanger_irecv_isend_new(int nsend, int nrecv,
                             const struct Xt_redist_msg *send_msgs,
                             const struct Xt_redist_msg *recv_msgs,
                             MPI_Comm comm, int tag_offset,
                             Xt_config config) {

  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */
  return xt_exchanger_simple_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                      comm, tag_offset,
                                      xt_exchanger_irecv_isend_s_exchange,
                                      xt_exchanger_irecv_isend_a_exchange,
                                      config);
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
