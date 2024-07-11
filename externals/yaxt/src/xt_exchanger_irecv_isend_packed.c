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

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_request_msgs_ebuf.h"
#include "xt_request_msgs_ebuf_internal.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger_irecv_isend_packed.h"
#include "xt_exchanger_simple_base.h"

/* unfortunately GCC 11 to 13 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static size_t
get_buffer_offsets(size_t *restrict buf_ofs,
                   int nsend, int nrecv,
                   const struct Xt_redist_msg *send_msgs,
                   const struct Xt_redist_msg *recv_msgs,
                   MPI_Comm comm)
{
  int buf_size;
  size_t accum = 0;
  for (int i = 0; i < nrecv; ++i) {
    buf_ofs[i] = accum;
    xt_mpi_call(MPI_Pack_size(1, recv_msgs[i].datatype, comm, &buf_size),
                comm);
    accum += (size_t)buf_size;
  }
  for (int i = 0; i < nsend; ++i) {
    buf_ofs[nrecv+i] = accum;
    xt_mpi_call(MPI_Pack_size(1, send_msgs[i].datatype, comm,
                              &buf_size), comm);
    accum += (size_t)buf_size;
  }
  buf_ofs[nsend+nrecv] = accum;
  return accum;
}

static size_t
get_buffer_size(int nsend, int nrecv,
                const struct Xt_redist_msg *send_msgs,
                const struct Xt_redist_msg *recv_msgs,
                MPI_Comm comm)
{
  int buf_size;
  size_t accum = 0;
  for (int i = 0; i < nrecv; ++i) {
    xt_mpi_call(MPI_Pack_size(1, recv_msgs[i].datatype, comm, &buf_size),
                comm);
    accum += (size_t)buf_size;
  }
  for (int i = 0; i < nsend; ++i) {
    xt_mpi_call(MPI_Pack_size(1, send_msgs[i].datatype, comm,
                              &buf_size), comm);
    accum += (size_t)buf_size;
  }

  return accum;
}

static void
start_packed_transfer(unsigned char *buffer,
                      int send_start,
                      const size_t *buf_ofs,
                      const void *src_data, int nsend, int nrecv,
                      const struct Xt_redist_msg *send_msgs,
                      const struct Xt_redist_msg *recv_msgs,
                      MPI_Comm comm, int tag_offset,
                      MPI_Request *requests)
{
  for (int i = 0; i < nrecv; ++i) {
    int recv_size = (int)(buf_ofs[i+1] - buf_ofs[i]);
    xt_mpi_call(MPI_Irecv(buffer + buf_ofs[i], recv_size, MPI_PACKED,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+i), comm);
  }

  for (int i = 0; i < nsend; ++i) {
    int position = 0;
    int buf_size = (int)(buf_ofs[send_start+i+1] - buf_ofs[send_start+i]);
    xt_mpi_call(MPI_Pack(CAST_MPI_SEND_BUF(src_data), 1, send_msgs[i].datatype,
                         buffer + buf_ofs[send_start+i], buf_size, &position,
                         comm), comm);
    xt_mpi_call(MPI_Isend(buffer + buf_ofs[send_start+i], position, MPI_PACKED,
                          send_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          requests+nrecv+i), comm);
  }
}

enum { AUTO_ALLOC_SIZE = 32, };

static void
xt_exchanger_irecv_isend_packed_s_exchange(
  const void *src_data, void *dst_data,
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs, const struct Xt_redist_msg *recv_msgs,
  int tag_offset, MPI_Comm comm)
{
  MPI_Request *requests, requests_auto[AUTO_ALLOC_SIZE];
  size_t *buf_ofs, buf_ofs_auto[AUTO_ALLOC_SIZE+1];

  size_t num_tx = (size_t)nrecv + (size_t)nsend;
  if (num_tx <= AUTO_ALLOC_SIZE) {
    requests = requests_auto;
    buf_ofs = buf_ofs_auto;
  } else {
    requests = xmalloc((num_tx+(num_tx&1)) * sizeof (*requests) + (num_tx+1) * sizeof (*buf_ofs));
    buf_ofs = (void *)(requests + (num_tx+(num_tx&1)));
  }

  size_t buffer_size = get_buffer_offsets(buf_ofs, nsend, nrecv,
                                          send_msgs, recv_msgs, comm);
  unsigned char *buffer = xmalloc(buffer_size);

  start_packed_transfer(buffer, nrecv, buf_ofs,
                        src_data, nsend, nrecv,
                        send_msgs, recv_msgs,
                        comm, tag_offset,
                        requests);

  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE), comm);

  for (int i = 0; i < nrecv; ++i) {
    int position = 0, recv_size = (int)(buf_ofs[i+1]-buf_ofs[i]);
    xt_mpi_call(MPI_Unpack(buffer + buf_ofs[i], recv_size, &position, dst_data,
                           1, recv_msgs[i].datatype, comm), comm);
  }

  free(buffer);
  if (num_tx > AUTO_ALLOC_SIZE)
    free(requests);
}

#ifdef _OPENMP
static void
xt_exchanger_irecv_isend_packed_s_exchange_omp(
  const void *src_data, void *dst_data,
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs, const struct Xt_redist_msg *recv_msgs,
  int tag_offset, MPI_Comm comm)
{
  MPI_Request *requests, requests_auto[AUTO_ALLOC_SIZE];
  size_t *buf_ofs, buf_ofs_auto[AUTO_ALLOC_SIZE+1];

  size_t num_tx = (size_t)nrecv + (size_t)nsend;
  if (num_tx <= AUTO_ALLOC_SIZE) {
    requests = requests_auto;
    buf_ofs = buf_ofs_auto;
  } else {
    requests = xmalloc((num_tx+(num_tx&1)) * sizeof (*requests) + (num_tx+1) * sizeof (*buf_ofs));
    buf_ofs = (size_t *)(requests + (num_tx+(num_tx&1)));
  }

  size_t buffer_size = get_buffer_offsets(buf_ofs, nsend, nrecv,
                                          send_msgs, recv_msgs, comm);
  unsigned char *buffer = xmalloc(buffer_size);

#pragma omp parallel
  {
    int num_threads = omp_get_num_threads(),
      tid = omp_get_thread_num();
    int start_send = (nsend * tid) / num_threads,
      nsend_ = (nsend * (tid+1)) / num_threads - start_send,
      start_recv = (nrecv * tid) / num_threads,
      end_recv = (nrecv * (tid+1)) / num_threads,
      nrecv_ = end_recv - start_recv,
      nreq = nrecv_+nsend_,
      start_req = start_send+start_recv;
    start_packed_transfer(buffer, start_send+nrecv-start_recv,
                          buf_ofs+start_recv,
                          src_data, nsend_, nrecv_,
                          send_msgs+start_send, recv_msgs+start_recv,
                          comm, tag_offset,
                          requests+start_req);

    xt_mpi_call(MPI_Waitall(nreq, requests+start_req, MPI_STATUSES_IGNORE),
                comm);
    for (int i = start_recv; i < end_recv; ++i) {
      int position = 0, recv_size = (int)(buf_ofs[i+1]-buf_ofs[i]);
      xt_mpi_call(MPI_Unpack(buffer + buf_ofs[i], recv_size, &position,
                             dst_data, 1, recv_msgs[i].datatype, comm), comm);
    }
  }

  free(buffer);
  if (num_tx > AUTO_ALLOC_SIZE)
    free(requests);
}
#endif

struct inventory_header
{
  int nsend, nrecv;
  void *dst_data;
};

/*
 * layout of inventory buffer created by
 * xt_exchanger_irecv_isend_packed_a_exchange:
 * struct inventory_header header;
 * size_t buf_ofs[nrecv+1];
 * MPI_Datatype dt[nrecv];
 * unsigned char packed_data_buf[sum_of_pack_and_unpack_sizes];
 */

static void
finalize_packed_a_exchange(Xt_request request, void *buf)
{
  struct inventory_header *header = buf;
  int nsend = header->nsend,
    nrecv = header->nrecv;
  (void)nsend;
  MPI_Comm comm = xt_request_msgs_ebuf_get_comm(request);
  MPI_Datatype *datatypes = (MPI_Datatype *)
    (void *)((unsigned char *)(header+1) + sizeof (size_t) * ((size_t)nrecv+1));
  void *dst_data = header->dst_data;
  size_t *buf_ofs = (void *)(header+1);
  for (int i = 0; i < nrecv; ++i) {
    int position = 0, buffer_size = (int)(buf_ofs[i+1]-buf_ofs[i]);
    xt_mpi_call(MPI_Unpack((unsigned char *)buf + buf_ofs[i], buffer_size,
                           &position, dst_data,
                           1, datatypes[i], comm), comm);
  }
  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Type_free(datatypes+i), comm);
}

static void
xt_exchanger_irecv_isend_packed_a_exchange(const void *src_data, void *dst_data,
                                           int nsend, int nrecv,
                                           const struct Xt_redist_msg * send_msgs,
                                           const struct Xt_redist_msg * recv_msgs,
                                           int tag_offset, MPI_Comm comm,
                                           Xt_request *request)
{
  size_t *buf_ofs, buf_ofs_auto[AUTO_ALLOC_SIZE+1];

  size_t num_tx = (size_t)nrecv + (size_t)nsend;
  if (num_tx <= AUTO_ALLOC_SIZE)
    buf_ofs = buf_ofs_auto;
  else
    buf_ofs = xmalloc((num_tx+1) * sizeof (*buf_ofs));

  size_t buffer_size = get_buffer_offsets(buf_ofs, nsend, nrecv,
                                          send_msgs, recv_msgs, comm);
  size_t inventory_size
    = sizeof (struct inventory_header)
    + sizeof (size_t) * ((size_t)nrecv+1)
    + sizeof (MPI_Datatype) * (size_t)nrecv;
  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_NONE);
  Xt_request requests
    = xt_request_msgs_ebuf_alloc(nrecv + nsend, comm,
                                 inventory_size + buffer_size, &conf);
  xt_request_msgs_ebuf_set_finalizer(requests, finalize_packed_a_exchange);

  MPI_Request *tmp_requests
    = xt_request_msgs_ebuf_get_req_ptr(requests);

  unsigned char *buffer
    = xt_request_msgs_ebuf_get_extra_buf(requests);

  start_packed_transfer(buffer+inventory_size, nrecv, buf_ofs,
                        src_data, nsend, nrecv,
                        send_msgs, recv_msgs,
                        comm, tag_offset,
                        tmp_requests);

  size_t *buf_ofs_ = (void *)(buffer + sizeof (struct inventory_header));
  for (int i = 0; i <= nrecv; ++i) {
    buf_ofs_[i] = buf_ofs[i]+inventory_size;
  }

  struct inventory_header *header = (void *)buffer;
  header->nsend = nsend;
  header->nrecv = nrecv;
  header->dst_data = dst_data;

  MPI_Datatype *datatypes = (void *)(buffer + sizeof (struct inventory_header)
                                     + sizeof (size_t) * ((size_t)nrecv+1));
  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Type_dup(recv_msgs[i].datatype, datatypes+i), comm);
  *request = requests;

  if (num_tx > AUTO_ALLOC_SIZE)
    free(buf_ofs);
}

#ifdef _OPENMP
static void
finalize_packed_a_exchange_mt(Xt_request request, void *buf)
{
  struct inventory_header *header = buf;
  int nsend = header->nsend,
    nrecv = header->nrecv;
  (void)nsend;
  MPI_Comm comm = xt_request_msgs_ebuf_get_comm(request);
  MPI_Datatype *datatypes = (MPI_Datatype *)
    ((unsigned char *)(header+1) + sizeof (size_t) * ((size_t)nrecv+1));
  void *dst_data = header->dst_data;
  size_t *buf_ofs = (void *)(header+1);
  int num_threads = omp_get_num_threads(),
    tid = omp_get_thread_num();
  int start_recv = (nrecv * tid) / num_threads,
    end_recv = (nrecv * (tid+1)) / num_threads;
  for (int i = start_recv; i < end_recv; ++i) {
    int position = 0, buffer_size = (int)(buf_ofs[i+1]-buf_ofs[i]);
    xt_mpi_call(MPI_Unpack((unsigned char *)buf + buf_ofs[i], buffer_size,
                           &position, dst_data,
                           1, datatypes[i], comm), comm);
  }
  for (int i = start_recv; i < end_recv; ++i)
    xt_mpi_call(MPI_Type_free(datatypes+i), comm);
}

static void
xt_exchanger_irecv_isend_packed_a_exchange_omp(const void *src_data, void *dst_data,
                                               int nsend, int nrecv,
                                               const struct Xt_redist_msg * send_msgs,
                                               const struct Xt_redist_msg * recv_msgs,
                                               int tag_offset, MPI_Comm comm,
                                               Xt_request *request)
{
  size_t *buf_ofs, buf_ofs_auto[AUTO_ALLOC_SIZE+1];

  size_t num_tx = (size_t)nrecv + (size_t)nsend;
  if (num_tx <= AUTO_ALLOC_SIZE)
    buf_ofs = buf_ofs_auto;
  else
    buf_ofs = xmalloc((num_tx+1) * sizeof (*buf_ofs));

  size_t buffer_size = get_buffer_offsets(buf_ofs, nsend, nrecv,
                                          send_msgs, recv_msgs, comm);
  size_t inventory_size
    = sizeof (struct inventory_header)
    + sizeof (size_t) * ((size_t)nrecv+1)
    + sizeof (MPI_Datatype) * (size_t)nrecv;
  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_OPENMP);
  Xt_request requests
    = xt_request_msgs_ebuf_alloc(nrecv + nsend, comm,
                                 inventory_size + buffer_size, &conf);
  xt_request_msgs_ebuf_set_finalizer(
    requests, finalize_packed_a_exchange_mt);
#pragma omp parallel firstprivate(requests, inventory_size)
  {
    MPI_Request *prequests
      = xt_request_msgs_ebuf_get_req_ptr(requests);

    unsigned char *buffer
      = xt_request_msgs_ebuf_get_extra_buf(requests);

    int num_threads = omp_get_num_threads(),
      tid = omp_get_thread_num();
    int start_send = (nsend * tid) / num_threads,
      nsend_ = (nsend * (tid+1)) / num_threads - start_send,
      start_recv = (nrecv * tid) / num_threads,
      end_recv = (nrecv * (tid+1)) / num_threads,
      nrecv_ = end_recv - start_recv,
      start_req = start_send+start_recv;

    start_packed_transfer(buffer+inventory_size,
                          start_send+nrecv-start_recv,
                          buf_ofs+start_recv,
                          src_data, nsend_, nrecv_,
                          send_msgs+start_send, recv_msgs+start_recv,
                          comm, tag_offset,
                          prequests+start_req);

    size_t *buf_ofs_ = (void *)(buffer + sizeof (struct inventory_header));
    for (int i = start_recv; i < end_recv+(end_recv==nrecv); ++i) {
      buf_ofs_[i] = buf_ofs[i]+inventory_size;
    }
#pragma omp master
    {
      struct inventory_header *header = (void *)buffer;
      header->nsend = nsend;
      header->nrecv = nrecv;
      header->dst_data = dst_data;
    }
    MPI_Datatype *datatypes = (void *)(buffer + sizeof (struct inventory_header)
                                       + sizeof (size_t) * ((size_t)nrecv+1));
    for (int i = start_recv; i < end_recv; ++i)
      xt_mpi_call(MPI_Type_dup(recv_msgs[i].datatype, datatypes+i), comm);
  }
  *request = requests;

  if (num_tx > AUTO_ALLOC_SIZE)
    free(buf_ofs);
}
#endif

static Xt_exchanger_omp_share
xt_exchanger_irecv_isend_packed_create_omp_share(
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs,
  const struct Xt_redist_msg *recv_msgs,
  MPI_Comm comm)
{
  size_t buf_size = get_buffer_size(nsend, nrecv, send_msgs, recv_msgs, comm);
  size_t inventory_size
    = sizeof (struct inventory_header)
    + sizeof (size_t) * ((size_t)nrecv+1)
    + sizeof (MPI_Datatype) * (size_t)nrecv;
  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_OPENMP);
  Xt_request shared_req
    = xt_request_msgs_ebuf_alloc(nsend+nrecv, comm,
                                 inventory_size + buf_size, &conf);
  return (Xt_exchanger_omp_share)shared_req;
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
  static const xt_simple_s_exchange_func
    s_exch_by_mthread_mode[] = {
    xt_exchanger_irecv_isend_packed_s_exchange,
#ifdef _OPENMP
    xt_exchanger_irecv_isend_packed_s_exchange_omp,
#else
    (xt_simple_s_exchange_func)0,
#endif
  };
  static const xt_simple_a_exchange_func
    a_exch_by_mthread_mode[] = {
    xt_exchanger_irecv_isend_packed_a_exchange,
#ifdef _OPENMP
    xt_exchanger_irecv_isend_packed_a_exchange_omp,
#else
    (xt_simple_a_exchange_func)0,
#endif
  };
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  return xt_exchanger_simple_base_new(
    nsend, nrecv, send_msgs, recv_msgs,
    comm, tag_offset,
    s_exch_by_mthread_mode[mthread_mode],
    a_exch_by_mthread_mode[mthread_mode],
    xt_exchanger_irecv_isend_packed_create_omp_share,
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
