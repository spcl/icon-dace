/**
 * @file xt_request_msgs_ebuf_internal.h
 * @brief internal interfaces for xt_request_msgs_ebuf
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
#ifndef XT_REQUEST_MSGS_EBUF_INTERNAL_H
#define XT_REQUEST_MSGS_EBUF_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_request_msgs_ebuf.h"
#include "xt/xt_config.h"
#include "core/ppm_visibility.h"

/**
 * @param n_requests number of requests that can be stored
 * @param comm communicator to use for error handling
 * @param config custom configuration object handle
 * @param extra_buf_size size of buffer to allocate additionally
 * @return Xt_request_msgs object ready for \a n MPI requests
 */
PPM_DSO_INTERNAL Xt_request
xt_request_msgs_ebuf_alloc(int n_requests,
                           MPI_Comm comm,
                           size_t extra_buf_size,
                           Xt_config config);

/**
 * @param request xt_request_msgs_ebuf object
 * @return Pointer to first element of MPI_Request array for requests
 * stored in request object
 */
PPM_DSO_INTERNAL MPI_Request *
xt_request_msgs_ebuf_get_req_ptr(Xt_request request);

PPM_DSO_INTERNAL void
xt_request_msgs_ebuf_set_finalizer(Xt_request request,
                                   Xt_request_msgs_ebuf_finalizer finalize);

PPM_DSO_INTERNAL void *
xt_request_msgs_ebuf_get_extra_buf(Xt_request request);

PPM_DSO_INTERNAL MPI_Comm
xt_request_msgs_ebuf_get_comm(Xt_request request);


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
