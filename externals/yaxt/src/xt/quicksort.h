#ifndef XT_QUICKSORT_H
#define XT_QUICKSORT_H

/**
 * @file quicksort.h
 * @brief quicksort declaration
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>

#include "xt/sort_common.h"

/** quicksort changing structured values
  *
  * @param[in,out]  a            array of ints to be sorted
  * @param[in]      n            length of data
 **/
void xt_quicksort_int(int a[], size_t n);

/** quicksort changing structured values
  *
  * @param[in,out]  v            array to be sorted
  * @param[in]      n            length of data
 **/
void xt_quicksort_idxpos(idxpos_type v[], size_t n);

/** quicksort changing values and indices
  *
  * @param[in,out]  a            array to be sorted
  * @param[in]      n            number of elements in a and idx
  * @param[in,out]  idx          old index of sorted returned a
  * @param[in]      reset_index  override given idx by identity idx
  */
void xt_quicksort_index(Xt_int *restrict a, int n, int *restrict idx,
                        int reset_index);

/** quicksort changing values and indices
  *
  * @param[in,out]  a            array to be sorted
  * @param[in]      n            number of elements in a and idx
  * @param[in,out]  permutation  contents permuted exactly as a and
  *                              used to resolve ordering if two
  *                              elements of \a a have the same value
  */
void
xt_quicksort_xt_int_permutation(Xt_int *restrict a, size_t n,
                                int *restrict permutation);

/** quicksort sorting values and indices
  *
  * @param[in,out]  a            array of n ints to be sorted
  * @param[in]      n            length of data
  * @param[in,out]  permutation  contents permuted exactly as a and
  *                              used to resolve ordering if two
  *                              elements of \a a have the same value
  */
void
xt_quicksort_int_permutation(int *restrict a, size_t n, int *restrict permutation);

#endif  /* XT_QUICKSORT_H */

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
