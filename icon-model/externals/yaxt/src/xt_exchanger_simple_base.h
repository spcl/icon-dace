/**
 * @file xt_exchanger_simple_base.h
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
#ifndef XT_EXCHANGER_SIMPLE_BASE_H
#define XT_EXCHANGER_SIMPLE_BASE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_core.h"
#include "xt/xt_config.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"

typedef void (*xt_simple_s_exchange_func)(
  const void *src_data, void *dst_data, int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs, const struct Xt_redist_msg *recv_msgs,
  int tag_offset, MPI_Comm comm);

typedef void (*xt_simple_a_exchange_func)(
  const void *src_data, void *dst_data, int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs, const struct Xt_redist_msg *recv_msgs,
  int tag_offset, MPI_Comm comm, Xt_request *request);

typedef Xt_exchanger_omp_share (*xt_simple_create_omp_share_func)(
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs, const struct Xt_redist_msg *recv_msgs,
  MPI_Comm comm);

/**
 * constructor for an exchanger using asynchronous send and recv
 * @param[in] nsend      number of send messages
 * @param[in] nrecv      number of receive messages
 * @param[in] send_msgs  array with send messages
 * @param[in] recv_msgs  array with receive messages
 * @param[in] comm       MPI communicator that is to be used for the
 *                       communication
 * @param[in] tag_offset tag
 * @param[in] s_func     function pointer used for synchronous exchanges
 * @param[in] a_func     function pointer used for asynchronous exchanges
 * @param[in] create_omp_share_func function pointer used to create
 * shared state for multi-threaded calls
 * @param[in] config     optional customization parameters
 * @remark tag_offset + xt_mpi_tag_exchange_msg must not
 *         be used on @a comm by any other part of the program during the
 *         lifetime of the created exchanger object
 */
PPM_DSO_INTERNAL Xt_exchanger
xt_exchanger_simple_base_new(
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs,
  const struct Xt_redist_msg *recv_msgs,
  MPI_Comm comm, int tag_offset,
  xt_simple_s_exchange_func s_func,
  xt_simple_a_exchange_func a_func,
  xt_simple_create_omp_share_func create_omp_share_func,
  Xt_config config);

PPM_DSO_INTERNAL extern const struct xt_exchanger_vtable
xt_exchanger_simple_base_vtable;


#endif // XT_EXCHANGER_SIMPLE_BASE_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
