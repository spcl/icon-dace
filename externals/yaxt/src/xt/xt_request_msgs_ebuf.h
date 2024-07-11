/**
 * @file xt_request_msgs_ebuf.h
 * @brief functions to create collection of request handles augmented
 *        with user-defined buffer
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
#ifndef XT_REQUEST_MSGS_EBUF_H
#define XT_REQUEST_MSGS_EBUF_H

#include <mpi.h>

#include "xt/xt_request.h"
#include "xt/xt_config.h"

/**
 * function called after the xt_request_msgs_ebuf object has been
 * created to actually create the requests
 * @param requests_msgs handle of requests just created
 * @param requests points to requests array which must be initialized
 * in this function
 * @param ebuf points to memory of extra_buf_size chars, aligned like
 * a data pointer
 */
typedef void (*Xt_fill_ebuf_requests)(Xt_request requests_msgs,
                                      MPI_Request *requests,
                                      void *ebuf,
                                      void *data,
                                      MPI_Comm comm);

/**
 * signature of function to be called after all requests of the
 * xt_request_msgs_ebuf object have finished
 */
typedef void (*Xt_request_msgs_ebuf_finalizer)(Xt_request request_msgs,
                                               void *ebuf);

/**
 * constructor for message request handle associated with extra buffer
 * @param[in] n_requests     number of requests needed
 * @param[in] extra_buf_size amount of extra bytes to allocate and
 *                           make available to \a init
 * @param[in] init           function to call before returning which
 *                           will then initialize the requests array
 *                           completely and the extra buffer as
 *                           needed
 * @param[in] data           pointer to caller-defined additional information
 *                           needed by and passed to init
 * @param[in] finalize       pointer to function to call when deleting
 *                           the request handle object (because all
 *                           requests are finished)
 * @param[in] comm           MPI communicator
 * @return request handle associated with extra buffer
 */
Xt_request
xt_request_msgs_ebuf_new(int n_requests,
                         size_t extra_buf_size,
                         Xt_fill_ebuf_requests init,
                         void *data,
                         Xt_request_msgs_ebuf_finalizer finalize,
                         MPI_Comm comm);

/**
 * customizable constructor for message request handle
 * @param[in] n_requests    number of entries in requests array
 * @param[in] extra_buf_size amount of extra bytes to allocate and
 *                           make available to \a init
 * @param[in] init           function to call before returning which
 *                           will then initialize the requests array
 *                           completely and the extra buffer as
 *                           needed
 * @param[in] data           pointer to caller-defined additional information
 *                           needed by and passed to init
 * @param[in] finalize       pointer to function to call when deleting
 *                           the request handle object (because all
 *                           requests are finished)
 * @param[in] comm           MPI communicator
 * @param[in] config         configuration object handle
 * @return request handle associated with extra buffer
 */
Xt_request
xt_request_msgs_ebuf_custom_new(int n_requests,
                                size_t extra_buf_size,
                                Xt_fill_ebuf_requests init,
                                void *data,
                                Xt_request_msgs_ebuf_finalizer finalize,
                                MPI_Comm comm,
                                Xt_config config);

#endif // XT_REQUEST_MSGS_EBUF_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
