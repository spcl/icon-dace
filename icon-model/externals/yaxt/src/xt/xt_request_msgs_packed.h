/**
 * @file xt_request_msgs_packed.h
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
#ifndef XT_REQUEST_MSGS_PACKED_H
#define XT_REQUEST_MSGS_PACKED_H

#include <mpi.h>

#include "xt/xt_request.h"

/**
 * constructor for message request handle
 * @param[in] n_requests    number of entries in requests array
 * @param[in] requests      array containg MPI requests
 * @param[in] comm          MPI communicator
 * @param[in] n_packed      number of entries in datatypes and packed_data
 * @param[in] n_tmp_buffers number of entries in tmp_buffers
 * @param[in] datatypes     array of datatypes to be used for unpacking
 * @param[in] packed_data   array of buffers containing packed data
 * @param[in] tmp_buffers   array of buffers that need to be freed after the
 *                          completion of exchange
 * @param[in] unpacked_data target buffer for unpacking
 * @remark ownership of the MPI requests is passed to the Xt_request object,
 *         however the caller remains the owner of the requests array
 * @remark ownership of the MPI datatypes and the array datatypes remain with
 *         the caller
 * @remark ownership of the buffers packed_data[0..n_packed-1] and
 *         tmp_buffers[0..n_tmp_buffers-1] is passed to the Xt_request object,
 *         however the caller remain the owner of the packed_data and
 *         tmp_buffers array
 */
Xt_request xt_request_msgs_packed_new(int n_requests,
                                      const MPI_Request requests[n_requests],
                                      MPI_Comm comm, int n_packed,
                                      int n_tmp_buffers,
                                      const MPI_Datatype datatypes[n_packed],
                                      void *packed_data[n_packed],
                                      void *tmp_buffers[n_tmp_buffers],
                                      void *unpacked_data);

#endif // XT_REQUEST_MSGS_PACKED_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
