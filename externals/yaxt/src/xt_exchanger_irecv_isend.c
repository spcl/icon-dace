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

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/ppm_xfuncs.h"
#include "xt_config_internal.h"
#include "xt/xt_mpi.h"
#include "xt/xt_request_msgs.h"
#include "xt_request_msgs_internal.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_simple_base.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static void
redist_msgs_to_req(const void *src_data, void *dst_data,
                   int nsend, int nrecv,
                   const struct Xt_redist_msg *send_msgs,
                   const struct Xt_redist_msg *recv_msgs,
                   int tag_offset, MPI_Comm comm,
                   MPI_Request *requests)
{
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
}

static void
xt_exchanger_irecv_isend_s_exchange(const void *src_data, void *dst_data,
                                    int nsend, int nrecv,
                                    const struct Xt_redist_msg *send_msgs,
                                    const struct Xt_redist_msg *recv_msgs,
                                    int tag_offset, MPI_Comm comm) {
  MPI_Request *requests
    = xmalloc((size_t)(nrecv + nsend) * sizeof (*requests));

  redist_msgs_to_req(src_data, dst_data, nsend, nrecv,
                     send_msgs, recv_msgs, tag_offset,
                     comm, requests);
  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE), comm);
  free(requests);
}

#ifdef _OPENMP
static void
xt_exchanger_irecv_isend_a_exchange_mt(const void *src_data, void *dst_data,
                                       int nsend, int nrecv,
                                       const struct Xt_redist_msg *send_msgs,
                                       const struct Xt_redist_msg *recv_msgs,
                                       int tag_offset, MPI_Comm comm,
                                       Xt_exchanger_omp_share shared_req) {
  MPI_Request *requests
    = xt_request_msgs_get_req_ptr((Xt_request)shared_req);
  int num_threads = omp_get_num_threads(),
    tid = omp_get_thread_num();
  int start_send = (nsend * tid) / num_threads,
    nsend_ = (nsend * (tid+1)) / num_threads - start_send,
    start_recv = (nrecv * tid) / num_threads,
    nrecv_ = (nrecv * (tid+1)) / num_threads - start_recv;

  redist_msgs_to_req(src_data, dst_data,
                     nsend_, nrecv_,
                     send_msgs+start_send, recv_msgs+start_recv, tag_offset,
                     comm, requests+start_recv+start_send);
}

static void
xt_exchanger_irecv_isend_s_exchange_mt(const void *src_data, void *dst_data,
                                       int nsend, int nrecv,
                                       const struct Xt_redist_msg *send_msgs,
                                       const struct Xt_redist_msg *recv_msgs,
                                       int tag_offset, MPI_Comm comm,
                                       Xt_exchanger_omp_share shared_req)
{
  int num_threads = omp_get_num_threads(),
    tid = omp_get_thread_num();
  xt_exchanger_irecv_isend_a_exchange_mt(src_data, dst_data,
                                         nsend, nrecv,
                                         send_msgs,
                                         recv_msgs,
                                         tag_offset, comm,
                                         shared_req);
  MPI_Request *requests
    = xt_request_msgs_get_req_ptr((Xt_request)shared_req);
  int start_send = (nsend * tid) / num_threads,
    nsend_ = (nsend * (tid+1)) / num_threads - start_send,
    start_recv = (nrecv * tid) / num_threads,
    nrecv_ = (nrecv * (tid+1)) / num_threads - start_recv,
    start_req = start_send+start_recv,
    nreq_ = nrecv_+nsend_;
  xt_mpi_call(MPI_Waitall(nreq_, requests+start_req,
                          MPI_STATUSES_IGNORE), comm);
}

static void
xt_exchanger_irecv_isend_s_exchange_omp(const void *src_data, void *dst_data,
                                        int nsend, int nrecv,
                                        const struct Xt_redist_msg *send_msgs,
                                        const struct Xt_redist_msg *recv_msgs,
                                        int tag_offset, MPI_Comm comm)
{
  Xt_exchanger_omp_share shared_req
    = (Xt_exchanger_omp_share)xt_request_msgs_alloc(nsend+nrecv, comm,
                                                    &xt_default_config);
#pragma omp parallel
  xt_exchanger_irecv_isend_s_exchange_mt(src_data, dst_data,
                                         nsend, nrecv,
                                         send_msgs,
                                         recv_msgs,
                                         tag_offset, comm,
                                         shared_req);
  free(shared_req);
}


#endif

static void
xt_exchanger_irecv_isend_a_exchange(const void *src_data, void *dst_data,
                                    int nsend, int nrecv,
                                    const struct Xt_redist_msg * send_msgs,
                                    const struct Xt_redist_msg * recv_msgs,
                                    int tag_offset, MPI_Comm comm,
                                    Xt_request *request) {

  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_NONE);
  Xt_request requests_ = *request
    = xt_request_msgs_alloc(nrecv + nsend, comm, &conf);

  MPI_Request *requests
    = xt_request_msgs_get_req_ptr(requests_);
  redist_msgs_to_req(src_data, dst_data, nsend, nrecv,
                     send_msgs, recv_msgs, tag_offset, comm, requests);
}

#ifdef _OPENMP
static void
xt_exchanger_irecv_isend_a_exchange_omp(const void *src_data, void *dst_data,
                                        int nsend, int nrecv,
                                        const struct Xt_redist_msg * send_msgs,
                                        const struct Xt_redist_msg * recv_msgs,
                                        int tag_offset, MPI_Comm comm,
                                        Xt_request *request) {

  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_OPENMP);
  Xt_exchanger_omp_share shared_req
    = (Xt_exchanger_omp_share)(
      *request = xt_request_msgs_alloc(nrecv + nsend, comm, &conf));
#pragma omp parallel
  xt_exchanger_irecv_isend_a_exchange_mt(src_data, dst_data,
                                         nsend, nrecv,
                                         send_msgs, recv_msgs,
                                         tag_offset, comm,
                                         shared_req);
}
#endif

static Xt_exchanger_omp_share
xt_exchanger_irecv_isend_create_omp_share(
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs,
  const struct Xt_redist_msg *recv_msgs,
  MPI_Comm comm)
{
  (void)send_msgs;
  (void)recv_msgs;
  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_OPENMP);
  int nmsg = nsend + nrecv;
  return (Xt_exchanger_omp_share)xt_request_msgs_alloc(nmsg, comm, &conf);
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
  static const xt_simple_s_exchange_func
    s_exch_by_mthread_mode[] = {
    xt_exchanger_irecv_isend_s_exchange,
#ifdef _OPENMP
    xt_exchanger_irecv_isend_s_exchange_omp,
#else
    (xt_simple_s_exchange_func)0,
#endif
  };
  static const xt_simple_a_exchange_func
    a_exch_by_mthread_mode[] = {
    xt_exchanger_irecv_isend_a_exchange,
#ifdef _OPENMP
    xt_exchanger_irecv_isend_a_exchange_omp,
#else
    (xt_simple_a_exchange_func)0,
#endif
  };
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  return xt_exchanger_simple_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                      comm, tag_offset,
                                      s_exch_by_mthread_mode[mthread_mode],
                                      a_exch_by_mthread_mode[mthread_mode],
                                      xt_exchanger_irecv_isend_create_omp_share,
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
