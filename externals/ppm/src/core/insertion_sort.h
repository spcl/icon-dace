/**
 * @file insertion_sort.h
 * @brief perform insertion sort on array
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
 *
 */

#ifndef PPM_INSERTION_SORT_H
#define PPM_INSERTION_SORT_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdlib.h>
#include <core/fptr_api.h>

/**
 * @brief perform insertion sort on array a
 * @param a array to sort
 * @param n number of elements in a
 * @param es size of element of a in bytes
 * @param data optional extra information to use in comparison
 * @param cmp comparison function
 */
void
PPM_insertion_sort(void *a, size_t n, size_t es, void *data,
                   PPM_CompareWithData cmp);


/**
 * @brief perform insertion of element a[n] into partially sorted
 * array a, i.e where elements 0 up to n - 1 are already in sort order
 *
 * This is a special case of PPM_insertion_sort_once where inversion
 * is equal to n-1.
 *
 * @param a array to insert into
 * @param n number of elements in a
 * @param es size of element of a in bytes
 * @param data optional extra information to use in comparison
 * @param cmp comparison function
 */
void
PPM_sorted_insertion(void *a, size_t n, size_t es, void *data,
                     PPM_CompareWithData cmp);

/**
 * @brief Perform sort on element a[inversion] only.
 *
 * It is assumed, that a[0..inversion-1,inversion+1..n-1] forms an
 * already sorted array, into which the misplaced element at index inversion
 * needs to be integrated.
 * @param a array to insert into
 * @param n number of elements in a
 * @param es size of element of a in bytes
 * @param inversion position of element not in sort order
 * @param data optional extra information to use in comparison
 * @param cmp comparison function
 */
void
PPM_insertion_sort_once(void *a, size_t n, size_t es, size_t inversion,
                        void *data, PPM_CompareWithData cmp);
/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */

#endif
