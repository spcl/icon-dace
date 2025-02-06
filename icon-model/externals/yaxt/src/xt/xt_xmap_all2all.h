/**
 * @file xt_xmap_all2all.h
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

#ifndef XT_XMAP_ALL2ALL_H
#define XT_XMAP_ALL2ALL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_core.h"

/** \example test_xmap_all2all.c
 */
/** \example test_xmap_all2all_f.f90
 */
/** \example test_xmap_all2all_parallel.c
 */
/** \example test_xmap_all2all_parallel_f.f90
 */
/** @example test_xmap_all2all_fail.c
 */
/** @example test_xmap_all2all_fail_f.f90
 */

/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm \n
 * it uses all-to-all communication for the initialisation
 *
 * @param[in] src_idxlist source index list
 * @param[in] dst_idxlist destination index list
 * @param[in] comm        MPI communicator that contains all processes
 *                        that take part in the exchange (xt_xmap_all2all_new
 *                        will make its own copy of comm)
 */
Xt_xmap xt_xmap_all2all_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                            MPI_Comm comm);

/**
 * constructor for an exchange map with custom parameters\n
 * this operation is collective over all processes in comm \n
 * it uses all-to-all communication for the initialisation
 *
 * @param[in] src_idxlist source index list
 * @param[in] dst_idxlist destination index list
 * @param[in] comm        MPI communicator that contains all processes
 *                        that take part in the exchange (xt_xmap_all2all_new
 *                        will make its own copy of comm)
 * @param[in] config      custom parameters
 */
Xt_xmap
xt_xmap_all2all_custom_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                           MPI_Comm comm, Xt_config config);

#endif // XT_XMAP_ALL2ALL_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
