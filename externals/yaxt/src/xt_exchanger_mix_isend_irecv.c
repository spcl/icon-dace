/**
 * @file xt_exchanger_mix_isend_irecv.c
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
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_config_internal.h"
#include "xt/xt_mpi.h"
#include "xt/xt_request_msgs.h"
#include "xt_request_msgs_internal.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_mix_isend_irecv.h"

/* unfortunately GCC 11 to 13 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static Xt_exchanger
xt_exchanger_mix_isend_irecv_copy(Xt_exchanger exchanger,
                                  MPI_Comm newComm, int new_tag_offset);
static void xt_exchanger_mix_isend_irecv_delete(Xt_exchanger exchanger);
static void
xt_exchanger_mix_isend_irecv_s_exchange(Xt_exchanger exchanger,
                                        const void *src_data,
                                        void *dst_data);
#ifdef _OPENMP
static void
xt_exchanger_mix_isend_irecv_s_exchange_omp(Xt_exchanger exchanger,
                                            const void *src_data,
                                            void *dst_data);
#endif
static void
xt_exchanger_mix_isend_irecv_a_exchange(
  Xt_exchanger exchanger, const void * src_data, void * dst_data,
  Xt_request *request);
#ifdef _OPENMP
static void
xt_exchanger_mix_isend_irecv_a_exchange_omp(
  Xt_exchanger exchanger, const void * src_data, void * dst_data,
  Xt_request *request);
#endif
static Xt_exchanger_omp_share
xt_exchanger_mix_isend_irecv_create_omp_share(Xt_exchanger exchanger);

static int
xt_exchanger_mix_isend_irecv_get_msg_ranks(Xt_exchanger exchanger,
                                           enum xt_msg_direction direction,
                                           int *restrict *ranks);

static MPI_Datatype
xt_exchanger_mix_isend_irecv_get_MPI_Datatype(Xt_exchanger exchanger,
                                              int rank,
                                              enum xt_msg_direction direction,
                                              bool do_dup);

const struct xt_exchanger_vtable
xt_exchanger_mix_isend_irecv_vtable = {
  .copy = xt_exchanger_mix_isend_irecv_copy,
  .delete = xt_exchanger_mix_isend_irecv_delete,
  .s_exchange = xt_exchanger_mix_isend_irecv_s_exchange,
  .a_exchange = xt_exchanger_mix_isend_irecv_a_exchange,
  .get_msg_ranks = xt_exchanger_mix_isend_irecv_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_mix_isend_irecv_get_MPI_Datatype,
  .create_omp_share = xt_exchanger_mix_isend_irecv_create_omp_share,
}
#ifdef _OPENMP
  ,
xt_exchanger_mix_isend_irecv_auto_omp_vtable = {
  .copy = xt_exchanger_mix_isend_irecv_copy,
  .delete = xt_exchanger_mix_isend_irecv_delete,
  .s_exchange = xt_exchanger_mix_isend_irecv_s_exchange_omp,
  .a_exchange = xt_exchanger_mix_isend_irecv_a_exchange_omp,
  .get_msg_ranks = xt_exchanger_mix_isend_irecv_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_mix_isend_irecv_get_MPI_Datatype,
  .create_omp_share = xt_exchanger_mix_isend_irecv_create_omp_share,
}
#endif
;

typedef struct Xt_exchanger_mix_isend_irecv_ * Xt_exchanger_mix_isend_irecv;

#define MSG_DIR(msg) ((msg).rank < 0)

struct Xt_exchanger_mix_isend_irecv_ {

  const struct xt_exchanger_vtable * vtable;

  int n, tag_offset;
  MPI_Comm comm;
  struct Xt_redist_msg msgs[];
};

static Xt_exchanger_mix_isend_irecv
xt_exchanger_mix_isend_irecv_alloc(size_t nmsg,
                                   Xt_config config)
{
  (void)config;
  Xt_exchanger_mix_isend_irecv exchanger;
  size_t header_size = sizeof (*exchanger),
    body_size = sizeof (struct Xt_redist_msg) * nmsg;
  exchanger = xmalloc(header_size + body_size);
  exchanger->n = (int)nmsg;
#ifdef _OPENMP
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  if (mthread_mode == XT_MT_OPENMP)
    exchanger->vtable = &xt_exchanger_mix_isend_irecv_auto_omp_vtable;
  else
#endif
    exchanger->vtable = &xt_exchanger_mix_isend_irecv_vtable;
  return exchanger;
}

static inline int
adjusted_rank(int r, int comm_rank, int comm_size)
{
  int r_ = r & INT_MAX;
  return r_ + (r_ <= comm_rank ? comm_size : 0);
}

#define XT_SORTFUNC_DECL static
#define SORT_TYPE struct Xt_redist_msg
#define SORT_TYPE_SUFFIX redist_msg_mix
#define SORT_TYPE_CMP_LT(u,v,i,j)                       \
  (adjusted_rank((u).rank, comm_rank, comm_size)        \
   < adjusted_rank((v).rank, comm_rank, comm_size))
#define SORT_TYPE_CMP_LE(u,v,i,j)                       \
  (adjusted_rank((u).rank, comm_rank, comm_size)        \
   <= adjusted_rank((v).rank, comm_rank, comm_size))
#define SORT_TYPE_CMP_EQ(u,v,i,j) \
  (((u).rank & INT_MAX) == ((v).rank & INT_MAX))
#define XT_SORT_EXTRA_ARGS_DECL , int comm_rank, int comm_size
#define XT_SORT_EXTRA_ARGS_PASS , comm_rank, comm_size
#define XT_SORT_VECSWAP_EXTRA_ARGS_DECL
#define XT_SORT_VECSWAP_EXTRA_ARGS_PASS

#include "xt_quicksort_base.h"


Xt_exchanger
xt_exchanger_mix_isend_irecv_new(int nsend, int nrecv,
                                 const struct Xt_redist_msg *send_msgs,
                                 const struct Xt_redist_msg *recv_msgs,
                                 MPI_Comm comm, int tag_offset,
                                 Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */

  assert((nsend >= 0) & (nrecv >= 0));
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  Xt_exchanger_mix_isend_irecv exchanger
    = xt_exchanger_mix_isend_irecv_alloc(nmsg, config);
  exchanger->comm = comm;
  exchanger->tag_offset = tag_offset;
  struct Xt_redist_msg *restrict msgs = exchanger->msgs;
  bool dt_dup = !(config->flags & exch_no_dt_dup);
  xt_redist_msgs_strided_copy((size_t)nsend, send_msgs, sizeof (send_msgs[0]),
                              msgs, sizeof (msgs[0]), comm, dt_dup);
  xt_redist_msgs_strided_copy((size_t)nrecv, recv_msgs, sizeof (recv_msgs[0]),
                              msgs+nsend, sizeof (msgs[0]), comm,
                              dt_dup);
  for (size_t i = 0; i < (size_t)nrecv; ++i)
    msgs[i + (size_t)nsend].rank |= ~INT_MAX;

  {
    int comm_size, comm_rank, is_inter;
    xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
    xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
    int (*get_comm_size)(MPI_Comm, int *)
      = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
    xt_mpi_call(get_comm_size(comm, &comm_size), comm);
    xt_quicksort_redist_msg_mix(msgs, nmsg, comm_rank, comm_size);
  }

  for (size_t i = 1; i < nmsg; ++i) {

    if ((msgs[i-1].rank & INT_MAX) == (msgs[i].rank & INT_MAX)
        && MSG_DIR(msgs[i]) == SEND) {

      struct Xt_redist_msg temp = msgs[i-1];
      msgs[i-1] = msgs[i];
      msgs[i] = temp;
      i++;
    }
  }

  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_mix_isend_irecv_copy(Xt_exchanger exchanger,
                                  MPI_Comm new_comm, int new_tag_offset)
{
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t nmsg = (size_t)exchanger_msr->n;
  /* fixme: needs to use custom config */
  Xt_exchanger_mix_isend_irecv exchanger_copy
    = xt_exchanger_mix_isend_irecv_alloc(nmsg, &xt_default_config);
  exchanger_copy->comm = new_comm;
  exchanger_copy->tag_offset = new_tag_offset;
  exchanger_copy->vtable = exchanger_msr->vtable;
  struct Xt_redist_msg *restrict new_msgs = exchanger_copy->msgs,
    *restrict orig_msgs = exchanger_msr->msgs;
  xt_redist_msgs_strided_copy(nmsg, orig_msgs, sizeof (*orig_msgs),
                              new_msgs, sizeof (*new_msgs),
                              new_comm, true);
  return (Xt_exchanger)exchanger_copy;
}

enum { max_on_stack_req = 16 };

static void xt_exchanger_mix_isend_irecv_delete(Xt_exchanger exchanger) {

  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;

  size_t nmsg = (size_t)exchanger_msr->n;
  struct Xt_redist_msg *restrict msgs = exchanger_msr->msgs;

  xt_redist_msgs_strided_destruct(nmsg, msgs, exchanger_msr->comm,
                                  sizeof (*msgs));
  free(exchanger_msr);
}

static inline void
redist_msgs_to_req(size_t nmsg,
                   const struct Xt_redist_msg *restrict msgs,
                   const void *src_data, void *dst_data,
                   MPI_Request *requests,
                   MPI_Comm comm, int tag_offset)
{
  for (size_t i = 0; i < nmsg; ++i) {
    typedef int (*ifp)(void *buf, int count, MPI_Datatype datatype, int dest,
                       int tag, MPI_Comm comm, MPI_Request *request);
    ifp op = MSG_DIR(msgs[i]) == SEND ? (ifp)MPI_Isend : (ifp)MPI_Irecv;
    void *data = MSG_DIR(msgs[i]) == SEND ? (void *)src_data : dst_data;
    xt_mpi_call(op(data, 1, msgs[i].datatype, msgs[i].rank & INT_MAX,
                   tag_offset + xt_mpi_tag_exchange_msg,
                   comm, requests+i), comm);
  }
}


static void xt_exchanger_mix_isend_irecv_s_exchange(Xt_exchanger exchanger,
                                                    const void * src_data,
                                                    void * dst_data) {

  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;

  if (exchanger_msr->n > 0) {
    size_t nmsg = (size_t)exchanger_msr->n;
    MPI_Request req_buf[max_on_stack_req];
    MPI_Request *requests
      = nmsg <= max_on_stack_req
      ? req_buf : xmalloc(nmsg * sizeof (*requests));
    redist_msgs_to_req(nmsg, exchanger_msr->msgs,
                       src_data, dst_data, requests,
                       exchanger_msr->comm, exchanger_msr->tag_offset);
    xt_mpi_call(MPI_Waitall((int)nmsg, requests, MPI_STATUSES_IGNORE),
                exchanger_msr->comm);
    if (requests != req_buf)
      free(requests);
  }
}

static void xt_exchanger_mix_isend_irecv_a_exchange(
  Xt_exchanger exchanger, const void * src_data, void * dst_data,
  Xt_request *request) {

  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;

  if (exchanger_msr->n > 0) {
    size_t nmsg = (size_t)exchanger_msr->n;
    struct Xt_config_ conf = xt_default_config;
    xt_config_set_redist_mthread_mode(&conf, XT_MT_NONE);
    Xt_request requests
      = xt_request_msgs_alloc((int)nmsg, exchanger_msr->comm, &conf);
    MPI_Request *requests_
      = xt_request_msgs_get_req_ptr(requests);
    redist_msgs_to_req(nmsg, exchanger_msr->msgs,
                       src_data, dst_data, requests_,
                       exchanger_msr->comm, exchanger_msr->tag_offset);
    *request = requests;
  } else
    *request = XT_REQUEST_NULL;
}

#ifdef _OPENMP
static void
xt_exchanger_mix_isend_irecv_a_exchange_mt(Xt_exchanger exchanger,
                                           const void * src_data,
                                           void * dst_data,
                                           Xt_exchanger_omp_share shared_req)
{
  MPI_Request *requests
    = xt_request_msgs_get_req_ptr((Xt_request)shared_req);
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t num_threads = (size_t)omp_get_num_threads(),
    tid = (size_t)omp_get_thread_num();
  size_t nmsg = (size_t)exchanger_msr->n,
    start = (nmsg * tid) / num_threads,
    nmsg_ = (nmsg * (tid+1)) / num_threads - start;
  redist_msgs_to_req(nmsg_, exchanger_msr->msgs+start,
                     src_data, dst_data, requests+start,
                     exchanger_msr->comm, exchanger_msr->tag_offset);
}

static void
xt_exchanger_mix_isend_irecv_a_exchange_omp(Xt_exchanger exchanger,
                                            const void *src_data,
                                            void *dst_data,
                                            Xt_request *request)
{
  Xt_exchanger_omp_share shared_req
    = xt_exchanger_mix_isend_irecv_create_omp_share(exchanger);
#pragma omp parallel
  xt_exchanger_mix_isend_irecv_a_exchange_mt(exchanger, src_data, dst_data,
                                             shared_req);
  *request = (Xt_request)shared_req;
}

static void
xt_exchanger_mix_isend_irecv_s_exchange_mt(Xt_exchanger exchanger,
                                           const void * src_data,
                                           void * dst_data,
                                           Xt_exchanger_omp_share shared_req)
{
  MPI_Request *requests
    = xt_request_msgs_get_req_ptr((Xt_request)shared_req);
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t num_threads = (size_t)omp_get_num_threads(),
    tid = (size_t)omp_get_thread_num();
  size_t nmsg = (size_t)exchanger_msr->n,
    start = (nmsg * tid) / num_threads,
    nmsg_ = (nmsg * (tid+1)) / num_threads - start;
  redist_msgs_to_req(nmsg_, exchanger_msr->msgs+start,
                     src_data, dst_data, requests+start,
                     exchanger_msr->comm, exchanger_msr->tag_offset);
  xt_mpi_call(MPI_Waitall((int)nmsg_, requests+start,
                          MPI_STATUSES_IGNORE), exchanger_msr->comm);
}

static void
xt_exchanger_mix_isend_irecv_s_exchange_omp(Xt_exchanger exchanger,
                                            const void * src_data,
                                            void * dst_data)
{
  Xt_exchanger_omp_share shared_req
    = xt_exchanger_mix_isend_irecv_create_omp_share(exchanger);
#pragma omp parallel
  xt_exchanger_mix_isend_irecv_s_exchange_mt(exchanger, src_data, dst_data,
                                             shared_req);
  free(shared_req);
}
#endif

static Xt_exchanger_omp_share
xt_exchanger_mix_isend_irecv_create_omp_share(Xt_exchanger exchanger)
{
  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_OPENMP);
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  return (Xt_exchanger_omp_share)xt_request_msgs_alloc(exchanger_msr->n,
                                                       exchanger_msr->comm,
                                                       &conf);
}

static int
xt_exchanger_mix_isend_irecv_get_msg_ranks(Xt_exchanger exchanger,
                                           enum xt_msg_direction direction,
                                           int *restrict *ranks)
{
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t nmsg = 0, nmsg_all = (size_t)exchanger_msr->n;
  const struct Xt_redist_msg *restrict msgs = exchanger_msr->msgs;
  for (size_t i = 0; i < nmsg_all; ++i)
    nmsg += MSG_DIR(msgs[i]) == direction;
  int *restrict ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc(nmsg * sizeof (*ranks_));
  for (size_t i = 0, j = (size_t)-1; i < nmsg_all; ++i)
    if (MSG_DIR(msgs[i]) == direction)
      ranks_[++j] = msgs[i].rank & INT_MAX;
  return (int)nmsg;
}

static MPI_Datatype
xt_exchanger_mix_isend_irecv_get_MPI_Datatype(Xt_exchanger exchanger,
                                              int rank,
                                              enum xt_msg_direction direction,
                                              bool do_dup)
{
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t nmsg = (size_t)exchanger_msr->n;
  struct Xt_redist_msg *restrict msgs = exchanger_msr->msgs;
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  if (direction == RECV) rank |= ~INT_MAX;
  for (size_t i = 0; i < nmsg; ++i)
    if (msgs[i].rank == rank) {
      if (do_dup)
        xt_mpi_call(MPI_Type_dup(msgs[i].datatype, &datatype_copy),
                    exchanger_msr->comm);
      else
        datatype_copy = msgs[i].datatype;
      break;
    }
  return datatype_copy;
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
