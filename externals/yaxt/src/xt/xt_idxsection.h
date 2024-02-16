/**
 * @file xt_idxsection.h
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

#ifndef XT_IDXSECTION_H
#define XT_IDXSECTION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_idxlist.h"

/** \example test_idxsection.c
 */
/** \example test_idxsection_f.f90
 */

/**
 * generates an index list that is comprised of a section of a set of
 * indices that are arranged in an n-dimensional cartesian coordinate
 * system. The linear index correspondence is computed such that the right-most
 * dimension is the one where indices increase fastest (i.e. according
 * to C convention).

 * @param[in] start          lowest index of the global array (typically 0 or 1)
 * @param[in] num_dimensions number of dimensions
 * @param[in] global_size    global size of each dimension
 * @param[in] local_size     size of the local section in each dimension
 * @param[in] local_start    vector with the lowest position in each dimension
 *                           of the local window within the global index space
 *
 * \remarks Negative values for global_size and local_size are allowed. Negative
 *          signs do not change the set of selected indices, only their ordering
 *          (see \ref idxsection_docu_neg_size).
 */
Xt_idxlist xt_idxsection_new(Xt_int start, int num_dimensions,
                             const Xt_int global_size[num_dimensions],
                             const int local_size[num_dimensions],
                             const Xt_int local_start[num_dimensions]);

Xt_idxlist xt_idxsection_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm);

#endif // XT_IDXSECTION_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
