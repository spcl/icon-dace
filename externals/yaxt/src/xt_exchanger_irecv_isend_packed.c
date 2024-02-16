/**
 * @file xt_exchanger_irecv_isend_packed.c
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
#include "xt/xt_mpi.h"
#include "xt/xt_request_msgs_packed.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger_irecv_isend_packed.h"
#include "xt_exchanger_simple_base.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static void
xt_exchanger_irecv_isend_packed_s_exchange(
  const void *src_data, void *dst_data,
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs, const struct Xt_redist_msg *recv_msgs,
  int tag_offset, MPI_Comm comm)
{
  enum { AUTO_ALLOC_SIZE = 32, };
  MPI_Request *requests, requests_auto[AUTO_ALLOC_SIZE];
  int *buffer_sizes, buffer_sizes_auto[AUTO_ALLOC_SIZE];

  size_t num_tx = (size_t)nrecv + (size_t)nsend;
  if (num_tx <= AUTO_ALLOC_SIZE) {
    requests = requests_auto;
    buffer_sizes = buffer_sizes_auto;
  } else {
    requests = xmalloc(num_tx * sizeof (*requests));
    buffer_sizes = xmalloc(num_tx * sizeof (*buffer_sizes));
  }

  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Pack_size(1, recv_msgs[i].datatype, comm, buffer_sizes+i),
                comm);
  for (int i = 0; i < nsend; ++i)
    xt_mpi_call(MPI_Pack_size(1, send_msgs[i].datatype, comm,
                              buffer_sizes+nrecv+i), comm);
  size_t buffer_size = 0;
  for (size_t i = 0; i < num_tx; ++i)
    buffer_size += (size_t)buffer_sizes[i];

  unsigned char *buffer = xmalloc(buffer_size);

  size_t ofs = 0;
  for (int i = 0; i < nrecv; ++i) {
    int recv_size = buffer_sizes[i];
    xt_mpi_call(MPI_Irecv(buffer + ofs, recv_size, MPI_PACKED,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+i), comm);
    ofs += (size_t)recv_size;
  }

  for (int i = 0; i < nsend; ++i) {
    int position = 0;
    xt_mpi_call(MPI_Pack(CAST_MPI_SEND_BUF(src_data), 1, send_msgs[i].datatype,
                         buffer + ofs, buffer_sizes[nrecv+i], &position,
                         comm), comm);
    xt_mpi_call(MPI_Isend(buffer + ofs, position, MPI_PACKED,
                          send_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+nrecv+i), comm);
    ofs += (size_t)position;
  }

  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE), comm);

  ofs = 0;
  for (int i = 0; i < nrecv; ++i) {
    int position = 0, recv_size = buffer_sizes[i];
    xt_mpi_call(MPI_Unpack(buffer + ofs, recv_size, &position, dst_data,
                           1, recv_msgs[i].datatype, comm), comm);
    ofs += (size_t)recv_size;
  }

  free(buffer);
  if (num_tx > AUTO_ALLOC_SIZE) {
    free(buffer_sizes);
    free(requests);
  }
}

static void
xt_exchanger_irecv_isend_packed_a_exchange(const void *src_data, void *dst_data,
                                           int nsend, int nrecv,
                                           const struct Xt_redist_msg * send_msgs,
                                           const struct Xt_redist_msg * recv_msgs,
                                           int tag_offset, MPI_Comm comm,
                                           Xt_request *request) {

  MPI_Request * tmp_requests =
    xmalloc((size_t)(nrecv + nsend) * sizeof (*tmp_requests));
  void ** buffers =
    xmalloc((size_t)(nrecv + nsend) * sizeof (*buffers));

  int buffer_size;
  for (int i = 0; i < nrecv; ++i) {
    xt_mpi_call(MPI_Pack_size(1, recv_msgs[i].datatype, comm, &buffer_size),
                comm);
    buffers[i] = xmalloc((size_t)buffer_size);
    xt_mpi_call(MPI_Irecv(buffers[i], buffer_size, MPI_PACKED,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          tmp_requests+i), comm);
  }

  for (int i = 0; i < nsend; ++i) {
    int position = 0;
    xt_mpi_call(MPI_Pack_size(1, send_msgs[i].datatype, comm, &buffer_size),
                comm);
    buffers[nrecv + i] = xmalloc((size_t)buffer_size);
    xt_mpi_call(MPI_Pack((void*)src_data, 1, send_msgs[i].datatype,
                         buffers[nrecv + i], buffer_size, &position,
                         comm), comm);
    xt_mpi_call(MPI_Isend(buffers[nrecv + i], buffer_size, MPI_PACKED,
                          send_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          tmp_requests+nrecv+i), comm);
  }

  MPI_Datatype * datatypes = xmalloc((size_t)nrecv * sizeof (*datatypes));
  for (int i = 0; i < nrecv; ++i) datatypes[i] = recv_msgs[i].datatype;

  Xt_request requests =
    xt_request_msgs_packed_new(nrecv + nsend, tmp_requests, comm, nrecv, nsend,
                               datatypes, buffers, buffers + nrecv, dst_data);

  free(datatypes);
  free(buffers);
  free(tmp_requests);

  *request = requests;
}

Xt_exchanger
xt_exchanger_irecv_isend_packed_new(int nsend, int nrecv,
                                    const struct Xt_redist_msg *send_msgs,
                                    const struct Xt_redist_msg *recv_msgs,
                                    MPI_Comm comm, int tag_offset,
                                    Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */
  return xt_exchanger_simple_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                      comm, tag_offset,
                                      xt_exchanger_irecv_isend_packed_s_exchange,
                                      xt_exchanger_irecv_isend_packed_a_exchange,
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
