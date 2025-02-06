/**
 * @file xt_redist.h
 * @brief redistribution of data
 *
 * contains declaration the redistribution data structure, which
 * is derived from one or more xt_xmaps
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

#ifndef XT_REDIST_H
#define XT_REDIST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_core.h"
#include "xt/xt_request.h"

struct Xt_redist_msg {

  int rank;
#if SIZEOF_MPI_DATATYPE == 2 * SIZEOF_INT
  int padding;
#endif
  MPI_Datatype datatype;
};

/**
 * \example test_redist_common.c
 */

/**
 * redist copy constructor
 *
 * @param[in,out] redist redistribution structure
 */
Xt_redist xt_redist_copy(Xt_redist redist);

/**
 * destructor
 *
 * @param[in,out] redist redistribution structure
 */
void xt_redist_delete(Xt_redist redist);

/**
 * synchronous redistribution of data
 *
 * @param[in]     redist     redistribution structure
 * @param[in]     num_arrays number of base addresses in src_data and dst_data
 * @param[in]     src_data   array containing the addresses of the first
 *                           elements of the input data
 * @param[in,out] dst_data   array containing the addresses of the first
 *                           elements of the output data
 *
 * @remark The above implies that NULL or any other invalid pointer
 * must not be used in either @a src_data or @a dst_data.
 * @par See Also
 *   \ref correct_addresses
 */
void xt_redist_s_exchange(Xt_redist redist, int num_arrays,
                          const void **src_data, void **dst_data);

/**
 * asynchronous redistribution of data
 *
 * @param[in]     redist     redistribution structure
 * @param[in]     num_arrays number of base addresses in src_data and dst_data
 * @param[in]     src_data   array containing the addresses of the first
 *                           elements of the input data
 * @param[in,out] dst_data   array containing the addresses of the first
 *                           elements of the output data
 * @param[out] request pointer to a request object that can be used to complete
 *          an asynchronous exchange
 *
 * @remark The above implies that NULL or any other invalid pointer
 * must not be used in either @a src_data or @a dst_data.
 * @par See Also
 *   \ref correct_addresses
 */
void xt_redist_a_exchange(Xt_redist redist, int num_arrays,
                          const void **src_data, void **dst_data,
                          Xt_request *request);

/**
 * synchronous redistribution of data - single array case
 *
 * @param[in]     redist   redistribution structure
 * @param[in]     src_data address of the first element of the input data
 * @param[in,out] dst_data address of the first element of the output data
 *
 * @remark The above implies that NULL or any other invalid pointer
 * must not be used in either @a src_data or @a dst_data.
 */
void xt_redist_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data);

/**
 * asynchronous redistribution of data - single array case
 *
 * @param[in]     redist   redistribution structure
 * @param[in]     src_data address of the first element of the input data
 * @param[in,out] dst_data address of the first element of the output data
 * @param[out] request pointer to a request object that can be used to complete
 *          an asynchronous exchange
 *
 * @remark The above implies that NULL or any other invalid pointer
 * must not be used in either @a src_data or @a dst_data.
 */
void xt_redist_a_exchange1(Xt_redist redist, const void *src_data,
                           void *dst_data, Xt_request *request);

/**
 * gets the number of messages send from the local process in an exchange
 * operation
 *
 * @param[in] redist redistribution structure
 * @return number of messages sent in the exchange operation
 */
int xt_redist_get_num_send_msg(Xt_redist redist);

/**
 * gets the number of messages received by the local process in an exchange
 * operation
 *
 * @param[in] redist redistribution structure
 * @return number of messages received in the exchange operation
 */
int xt_redist_get_num_recv_msg(Xt_redist redist);

/**
 * gets a copy of the MPI_Datatype used for the data of the send operation with
 * the given rank
 *
 * @param[in] redist redistribution structure
 * @param[in] rank   MPI rank
 * \return MPI_Datatype for the data of the send operation with the given rank
 *
 * \remarks returns MPI_DATATYPE_NULL if there is no send operation with the given rank
 */
MPI_Datatype xt_redist_get_send_MPI_Datatype(Xt_redist redist, int rank);

/**
 * gets a copy of the MPI_Datatype used for the data of the recv operation with
 * the given rank
 *
 * @param[in] redist redistribution structure
 * @param[in] rank   MPI rank
 * \return MPI_Datatype for the data of the recv operation with the given rank
 *
 * \remarks returns MPI_DATATYPE_NULL if there is no recv operation with the given rank
 */
MPI_Datatype xt_redist_get_recv_MPI_Datatype(Xt_redist redist, int rank);

/**
 * returns a MPI communicator, which the redistribution is based on
 * @param[in] redist redistribution structure
 * \return MPI communicator, which the redistribution is based on
 */
MPI_Comm xt_redist_get_MPI_Comm(Xt_redist redist);

#endif // XT_REDIST_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
