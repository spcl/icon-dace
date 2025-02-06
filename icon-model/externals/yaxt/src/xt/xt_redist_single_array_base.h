/**
 * @file xt_redist_single_array_base.h
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

#ifndef XT_REDIST_SINGLE_ARRAY_BASE_H
#define XT_REDIST_SINGLE_ARRAY_BASE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_config.h"
#include "xt/xt_xmap.h"
#include "xt/xt_redist.h"

/** \example test_redist_single_array_base.c
 *  */
/** \example test_redist_single_array_base_f.f90
 *  */
/** \example test_redist_single_array_base_parallel.c
 *  */
/** \example test_redist_single_array_base_parallel_f.f90
 *  */

/**
 * constructor for a redistribution using point to point
 * communication for the exchange using default configuration
 * @param[in] nsend     number of send messages
 * @param[in] nrecv     number of receive messages
 * @param[in] send_msgs array with send messages
 * @param[in] recv_msgs array with receive messages
 * @param[in] comm      MPI communicator that is to be used for the
 *                      communication
 */
Xt_redist
xt_redist_single_array_base_new(int nsend, int nrecv,
                                const struct Xt_redist_msg send_msgs[],
                                const struct Xt_redist_msg recv_msgs[],
                                MPI_Comm comm);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange with custom configuration
 * @param[in] nsend     number of send messages
 * @param[in] nrecv     number of receive messages
 * @param[in] send_msgs array with send messages
 * @param[in] recv_msgs array with receive messages
 * @param[in] comm      MPI communicator that is to be used for the
 *                      communication
 * @param[in] config    configuration object for custom settings
 */
Xt_redist
xt_redist_single_array_base_custom_new(int nsend, int nrecv,
                                       const struct Xt_redist_msg send_msgs[],
                                       const struct Xt_redist_msg recv_msgs[],
                                       MPI_Comm comm, Xt_config config);

#endif // XT_REDIST_SINGLE_ARRAY_BASE_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
