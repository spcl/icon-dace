/**
 * @file xt_request_msgs_ebuf.c
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
#include <stdlib.h>
#include <string.h>

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
#include "xt_request_internal.h"

static void xt_request_msgs_ebuf_wait(Xt_request request);
static int xt_request_msgs_ebuf_test(Xt_request request);

#ifdef _OPENMP
static void
xt_request_msgs_ebuf_wait_omp(Xt_request request);
static int
xt_request_msgs_ebuf_test_omp(Xt_request request);
#endif

static const struct Xt_request_vtable request_msgs_ebuf_vtable = {

  .wait = xt_request_msgs_ebuf_wait,
  .test = xt_request_msgs_ebuf_test,
}
#ifdef _OPENMP
  , request_msgs_ebuf_auto_omp_vtable = {
  .wait = xt_request_msgs_ebuf_wait_omp,
  .test = xt_request_msgs_ebuf_test_omp,
};
#endif
;

typedef struct Xt_request_msgs_ebuf_ * Xt_request_msgs_ebuf;

struct Xt_request_msgs_ebuf_ {
  const struct Xt_request_vtable *vtable;
  MPI_Comm comm;
  int n_requests;
  void (*finalizer)(Xt_request, void *);
  size_t extra_buf_ofs;
  MPI_Request requests[];
};

#define no_finalizer ((void (*)(Xt_request, void *))0)

Xt_request
xt_request_msgs_ebuf_alloc(int n_requests,
                           MPI_Comm comm,
                           size_t extra_buf_size,
                           Xt_config config)
{
  size_t hdr_size = sizeof(struct Xt_request_msgs_ebuf_),
    bufr_size = (size_t)n_requests * sizeof(MPI_Request)
    + (size_t)n_requests * sizeof(int);
  /* round allocation up to next multiple of sizeof (void *) to place
   * beginning of user requested extra buffer */
  size_t ofs = ((hdr_size + bufr_size + sizeof (void *) - 1)
                   & -sizeof(void *));
  Xt_request_msgs_ebuf request = xmalloc(ofs + extra_buf_size);
#ifndef _OPENMP
  (void)config;
#else
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  if (mthread_mode == XT_MT_OPENMP)
    request->vtable = &request_msgs_ebuf_auto_omp_vtable;
  else
#endif
  request->vtable = &request_msgs_ebuf_vtable;
  request->comm = comm;
  request->n_requests = n_requests;
  request->finalizer = no_finalizer;
  request->extra_buf_ofs = ofs;
  return (Xt_request)request;
}

static inline void *
extra_buf(Xt_request_msgs_ebuf request_msgs_ebuf)
{
  size_t ofs = request_msgs_ebuf->extra_buf_ofs;
  return (unsigned char *)request_msgs_ebuf + ofs;
}

Xt_request
xt_request_msgs_ebuf_new(int n_requests,
                         size_t extra_buf_size,
                         Xt_fill_ebuf_requests init,
                         void *data,
                         Xt_request_msgs_ebuf_finalizer finalize,
                         MPI_Comm comm) {
  return xt_request_msgs_ebuf_custom_new(n_requests, extra_buf_size, init, data,
                                         finalize, comm, &xt_default_config);
}

Xt_request
xt_request_msgs_ebuf_custom_new(int n_requests,
                                size_t extra_buf_size,
                                Xt_fill_ebuf_requests init,
                                void *data,
                                Xt_request_msgs_ebuf_finalizer finalize,
                                MPI_Comm comm,
                                Xt_config config)
{
  assert(n_requests >= 0);
  Xt_request_msgs_ebuf request
    = (Xt_request_msgs_ebuf)xt_request_msgs_ebuf_alloc(
      n_requests, comm, extra_buf_size, config);
  xt_request_msgs_ebuf_set_finalizer((Xt_request)request, finalize);
  init((Xt_request)request, request->requests, extra_buf(request),
       data, comm);
  return (Xt_request)request;
}

static void
finish_requests(Xt_request_msgs_ebuf request_msgs_ebuf)
{
  if (request_msgs_ebuf->finalizer != no_finalizer)
    request_msgs_ebuf->finalizer((Xt_request)request_msgs_ebuf,
                                   extra_buf(request_msgs_ebuf));
  free(request_msgs_ebuf);
}

static void xt_request_msgs_ebuf_wait(Xt_request request) {

  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  xt_mpi_call(MPI_Waitall(request_msgs_ebuf->n_requests,
                          request_msgs_ebuf->requests, MPI_STATUSES_IGNORE),
              request_msgs_ebuf->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif
  finish_requests(request_msgs_ebuf);
}

static int xt_request_msgs_ebuf_test(Xt_request request) {

  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;

  int flag =
    xt_mpi_test_some(&(request_msgs_ebuf->n_requests),
                     request_msgs_ebuf->requests,
                     (int *)(request_msgs_ebuf->requests
                             + request_msgs_ebuf->n_requests),
                     request_msgs_ebuf->comm);

  if (flag)
    finish_requests(request_msgs_ebuf);

  return flag;
}

MPI_Request *
xt_request_msgs_ebuf_get_req_ptr(Xt_request request)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  return request_msgs_ebuf->requests;
}

MPI_Comm
xt_request_msgs_ebuf_get_comm(Xt_request request)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  return request_msgs_ebuf->comm;
}

void
xt_request_msgs_ebuf_set_finalizer(Xt_request request,
                                   Xt_request_msgs_ebuf_finalizer finalizer)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  request_msgs_ebuf->finalizer = finalizer;
}

void *
xt_request_msgs_ebuf_get_extra_buf(Xt_request request)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  return extra_buf(request_msgs_ebuf);
}


#ifdef _OPENMP
static void
finish_requests_mt(Xt_request_msgs_ebuf request_msgs_ebuf)
{
  if (request_msgs_ebuf->finalizer != no_finalizer)
    request_msgs_ebuf->finalizer((Xt_request)request_msgs_ebuf,
                                   extra_buf(request_msgs_ebuf));
#pragma omp barrier
#pragma omp master
  free(request_msgs_ebuf);
}

static void
xt_request_msgs_ebuf_wait_omp(Xt_request request)
{
#pragma omp parallel
  {
    Xt_request_msgs_ebuf request_msgs = (Xt_request_msgs_ebuf)request;
    int num_threads = omp_get_num_threads(),
      tid = omp_get_thread_num();
    int nreq = request_msgs->n_requests,
      start_req = (nreq * tid) / num_threads,
      nreq_ = (nreq * (tid+1)) / num_threads - start_req;

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
    xt_mpi_call(MPI_Waitall(nreq_, request_msgs->requests+start_req,
                            MPI_STATUSES_IGNORE), request_msgs->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif
#pragma omp barrier
    finish_requests_mt(request_msgs);
  }
}

static int
xt_request_msgs_ebuf_test_omp(Xt_request request)
{
  bool flag = true;
#pragma omp parallel
  {
    Xt_request_msgs_ebuf request_msgs = (Xt_request_msgs_ebuf)request;

    int *ops_completed_buffer
      = (int *)(request_msgs->requests + request_msgs->n_requests);
    bool flag_
      = xt_mpi_test_some_mt(&request_msgs->n_requests,
                            request_msgs->requests,
                            ops_completed_buffer,
                            request_msgs->comm);
#pragma omp master
    flag = flag_;
#pragma omp barrier
    if (flag_)
      finish_requests_mt(request_msgs);
  }

  return flag;
}
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
