/**
 * @file xt_exchanger_irecv_isend_ddt_packed.c
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
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
#include "xt_request_msgs_ddt_packed.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger_irecv_isend_ddt_packed.h"
#include "xt_exchanger_simple_base.h"
#include "xt_ddt_internal.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static void
xt_exchanger_irecv_isend_ddt_packed_s_exchange(
  const void *src_data, void *dst_data,
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs, const struct Xt_redist_msg *recv_msgs,
  int tag_offset, MPI_Comm comm) {

  XT_GPU_INSTR_PUSH(xt_exchanger_irecv_isend_ddt_packed_s_exchange);

  enum { AUTO_ALLOC_SIZE = 32, };
  MPI_Request *requests, requests_auto[AUTO_ALLOC_SIZE];
  size_t *buffer_sizes, buffer_sizes_auto[AUTO_ALLOC_SIZE];
  Xt_ddt *ddts, ddts_auto[AUTO_ALLOC_SIZE];

  size_t num_tx = (size_t)nrecv + (size_t)nsend;
  if (num_tx <= AUTO_ALLOC_SIZE) {
    requests = requests_auto;
    buffer_sizes = buffer_sizes_auto;
    ddts = ddts_auto;
  } else {
    requests = xmalloc(num_tx * sizeof (*requests));
    buffer_sizes = xmalloc(num_tx * sizeof (*buffer_sizes));
    ddts = xmalloc(num_tx * sizeof(*ddts));
  }

  size_t recv_buffer_size = 0;
  size_t send_buffer_size = 0;
  for (int i = 0; i < nrecv; ++i)
    recv_buffer_size +=
      ((buffer_sizes[i] =
          xt_ddt_get_pack_size_internal(
            ((ddts[i] = xt_ddt_from_mpi_ddt(recv_msgs[i].datatype))))));
  for (int i = 0; i < nsend; ++i)
    send_buffer_size +=
      ((buffer_sizes[nrecv+i] =
          xt_ddt_get_pack_size_internal(
            ((ddts[nrecv+i] = xt_ddt_from_mpi_ddt(send_msgs[i].datatype))))));

  enum xt_memtype src_data_memtype = xt_gpu_get_memtype(src_data);
  enum xt_memtype dst_data_memtype = xt_gpu_get_memtype(dst_data);
  unsigned char * send_buffer =
    xt_gpu_malloc(send_buffer_size, src_data_memtype);
  unsigned char * recv_buffer =
    xt_gpu_malloc(recv_buffer_size, dst_data_memtype);

  size_t ofs = 0;
  for (int i = 0; i < nrecv; ++i) {
    int recv_size = (int)buffer_sizes[i];
    xt_mpi_call(MPI_Irecv(recv_buffer + ofs, recv_size, MPI_BYTE,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+i), comm);
    ofs += (size_t)recv_size;
  }

  ofs = 0;
  for (int i = 0; i < nsend; ++i) {
    size_t send_size = buffer_sizes[nrecv+i];
    xt_ddt_pack_internal(
      ddts[nrecv+i], CAST_MPI_SEND_BUF(src_data), send_buffer + ofs,
      src_data_memtype);
    xt_mpi_call(MPI_Isend(send_buffer + ofs, (int)send_size, MPI_BYTE,
                          send_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+nrecv+i), comm);
    ofs += send_size;
  }

  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE), comm);

  ofs = 0;
  for (int i = 0; i < nrecv; ++i) {
    size_t recv_size = buffer_sizes[i];
    xt_ddt_unpack_internal(
      ddts[i], recv_buffer + ofs, dst_data, dst_data_memtype);
    ofs += recv_size;
  }

  xt_gpu_free(recv_buffer, dst_data_memtype);
  xt_gpu_free(send_buffer, src_data_memtype);
  if (num_tx > AUTO_ALLOC_SIZE) {
    free(ddts);
    free(buffer_sizes);
    free(requests);
  }
  XT_GPU_INSTR_POP; // xt_exchanger_irecv_isend_ddt_packed_s_exchange
}

static void
xt_exchanger_irecv_isend_ddt_packed_a_exchange(const void *src_data, void *dst_data,
                                               int nsend, int nrecv,
                                               const struct Xt_redist_msg * send_msgs,
                                               const struct Xt_redist_msg * recv_msgs,
                                               int tag_offset, MPI_Comm comm,
                                               Xt_request *request) {

  XT_GPU_INSTR_PUSH(xt_exchanger_irecv_isend_ddt_packed_a_exchange);

  MPI_Request * tmp_requests =
    xmalloc((size_t)(nrecv + nsend) * sizeof (*tmp_requests));
  void ** buffers =
    xmalloc((size_t)(nrecv + nsend) * sizeof (*buffers));

  enum xt_memtype src_data_memtype = xt_gpu_get_memtype(src_data);
  enum xt_memtype dst_data_memtype = xt_gpu_get_memtype(dst_data);

  Xt_ddt * recv_ddts = xmalloc((size_t)nrecv * sizeof(*recv_ddts));
  for (int i = 0; i < nrecv; ++i) {
    recv_ddts[i] = xt_ddt_from_mpi_ddt(recv_msgs[i].datatype);
    size_t buffer_size = xt_ddt_get_pack_size_internal(recv_ddts[i]);
    buffers[i] = xt_gpu_malloc(buffer_size, dst_data_memtype);
    xt_mpi_call(MPI_Irecv(buffers[i], (int)buffer_size, MPI_BYTE,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          tmp_requests+i), comm);
  }

  for (int i = 0; i < nsend; ++i) {
    Xt_ddt send_ddt = xt_ddt_from_mpi_ddt(send_msgs[i].datatype);
    size_t buffer_size = xt_ddt_get_pack_size_internal(send_ddt);
    buffers[nrecv + i] = xt_gpu_malloc(buffer_size, src_data_memtype);
//! \todo merge all packing kernels into single kernel call -> less overhead,
//!       but not overlapping of packing and sending
    xt_ddt_pack_internal(
      send_ddt, src_data, buffers[nrecv + i], src_data_memtype);
    xt_mpi_call(MPI_Isend(buffers[nrecv + i], (int)buffer_size, MPI_BYTE,
                          send_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          tmp_requests+nrecv+i), comm);
  }

  Xt_request requests =
    xt_request_msgs_ddt_packed_new(
      nrecv + nsend, tmp_requests, comm, nrecv, nsend,
      recv_ddts, buffers, buffers + nrecv, dst_data,
      src_data_memtype, dst_data_memtype);

  free(recv_ddts);
  free(buffers);
  free(tmp_requests);

  *request = requests;

  XT_GPU_INSTR_POP; // xt_exchanger_irecv_isend_ddt_packed_a_exchange
}

Xt_exchanger
xt_exchanger_irecv_isend_ddt_packed_new(int nsend, int nrecv,
                                        const struct Xt_redist_msg *send_msgs,
                                        const struct Xt_redist_msg *recv_msgs,
                                        MPI_Comm comm, int tag_offset,
                                        Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */
  return
    xt_exchanger_simple_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                 comm, tag_offset,
                                 xt_exchanger_irecv_isend_ddt_packed_s_exchange,
                                 xt_exchanger_irecv_isend_ddt_packed_a_exchange,
                                 (xt_simple_create_omp_share_func)0,
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
