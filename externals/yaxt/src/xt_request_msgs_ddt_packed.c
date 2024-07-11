/**
 * @file xt_request_msgs_ddt_packed.c
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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"

#include "xt_request_msgs_ddt_packed.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_request_internal.h"
#include "xt_ddt_internal.h"
#include "xt_request_internal.h"
#include "xt_request_msgs_ebuf_internal.h"

struct Xt_request_msgs_ddt_packed {
  int n_packed, n_tmp_buffers;
  void *dst_data;
  Xt_ddt *ddts;
  enum xt_memtype packed_memtype;
  enum xt_memtype tmp_memtype;
  void *buffers[];
};

static void
xt_request_msgs_packed_ddt_finalize(Xt_request request, void *ebuf)
{
  (void)request;
  //! \todo merge all unpacking kernels into single kernel call -> less overhead
  struct Xt_request_msgs_ddt_packed *request_msgs_ddt_packed = ebuf;
  int n_packed = request_msgs_ddt_packed->n_packed,
    n_tmp_buffers = request_msgs_ddt_packed->n_tmp_buffers;
  void **buffers = request_msgs_ddt_packed->buffers;
  Xt_ddt *ddts = request_msgs_ddt_packed->ddts;
  void *dst_data = request_msgs_ddt_packed->dst_data;
  enum xt_memtype pmemtype = request_msgs_ddt_packed->packed_memtype;
  for (int i = 0; i < n_packed; ++i) {
    xt_ddt_unpack_internal(ddts[i], buffers[i], dst_data, pmemtype);
    xt_gpu_free(buffers[i], pmemtype);
    xt_ddt_delete(ddts[i]);
  }
  for (int i = n_packed; i < n_packed + n_tmp_buffers; ++i)
    xt_gpu_free(buffers[i], request_msgs_ddt_packed->tmp_memtype);
  free(ddts);
}


Xt_request
xt_request_msgs_ddt_packed_new(int n_requests,
                               const MPI_Request requests[n_requests],
                               MPI_Comm comm, int n_packed,
                               int n_tmp_buffers,
                               Xt_ddt packed_ddts[n_packed],
                               void *packed_data[n_packed],
                               void *tmp_buffers[n_tmp_buffers],
                               void * dst_data,
                               enum xt_memtype packed_memtype,
                               enum xt_memtype tmp_memtype) {

  assert(n_requests >= 0 && n_packed >= 0 && n_tmp_buffers >= 0);
  size_t hdr_size = sizeof(struct Xt_request_msgs_ddt_packed),
    bufp_size
    = ((size_t)n_packed + (size_t)n_tmp_buffers) * sizeof(void *);

  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_NONE);
  Xt_request request
    = xt_request_msgs_ebuf_alloc(n_requests, comm, hdr_size+bufp_size, &conf);
  MPI_Request *requests_
    = xt_request_msgs_ebuf_get_req_ptr(request);
  memcpy(requests_, requests,
         (size_t)n_requests * sizeof(*requests_));
  struct Xt_request_msgs_ddt_packed *ebuf
    = xt_request_msgs_ebuf_get_extra_buf(request);
  ebuf->n_packed = n_packed;
  ebuf->n_tmp_buffers = n_tmp_buffers;
  Xt_ddt *ddts = ebuf->ddts = xmalloc((size_t)n_packed * sizeof(*ddts));
  for (int i = 0; i < n_packed; ++i) {
    ddts[i] = packed_ddts[i];
    xt_ddt_inc_ref_count(ddts[i]);
  }
  ebuf->packed_memtype = packed_memtype;
  ebuf->tmp_memtype = tmp_memtype;
  memcpy(ebuf->buffers, packed_data,
         (size_t)n_packed * sizeof(*ebuf->buffers));
  memcpy(ebuf->buffers + n_packed, tmp_buffers,
         (size_t)n_tmp_buffers * sizeof(*ebuf->buffers));
  ebuf->dst_data = dst_data;
  xt_request_msgs_ebuf_set_finalizer(request,
                                     xt_request_msgs_packed_ddt_finalize);
  return request;
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
