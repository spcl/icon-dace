/**
 * @file xt_redist_collection.h
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

#ifndef XT_REDIST_COLLECTION_H
#define XT_REDIST_COLLECTION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_config.h"
#include "xt/xt_redist.h"

/** \example test_redist_collection.c
 */
/** \example test_redist_collection_f.f90
 */
/** \example test_redist_collection_parallel.c
 */
/** \example test_redist_collection_parallel_f.f90
 */

/**
 * constructor for a redistribution collection that is comprised
 * of multiple other redistributions with default settings
 * @param[in] redists     redistributions
 * @param[in] num_redists number of redistributions
 * @param[in] cache_size  number of cached datatypes required for
 *                        different input array displacements \n
 *                        n = -1 => default cache size \n
 *                        n =  0 => no cache \n
 *                        n >  0 => cache size = n
 * @param[in] comm        MPI communicator
 *
 * @remarks all redistributions need to be based on the same
 *          MPI communicator
 * @par See Also
 *   \ref correct_addresses
 */
Xt_redist xt_redist_collection_new(Xt_redist * redists, int num_redists,
                                   int cache_size, MPI_Comm comm);

/**
 * constructor for a redistribution collection that is comprised
 * of multiple other redistributions with custom settings
 * @param[in] redists     redistributions
 * @param[in] num_redists number of redistributions
 * @param[in] cache_size  number of cached datatypes required for
 *                        different input array displacements \n
 *                        n = -1 => default cache size \n
 *                        n =  0 => no cache \n
 *                        n >  0 => cache size = n
 * @param[in] comm        MPI communicator
 * @param[in] config      custom configuration object
 *
 * @remarks all redistributions need to be based on the same
 *          MPI communicator
 * @par See Also
 *   \ref correct_addresses
 */
Xt_redist xt_redist_collection_custom_new(Xt_redist * redists, int num_redists,
                                          int cache_size, MPI_Comm comm,
                                          Xt_config config);

#endif // XT_REDIST_COLLECTION_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
