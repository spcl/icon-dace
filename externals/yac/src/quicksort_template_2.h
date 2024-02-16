/**
 * @file quicksort_template_2.h
 *
 * @copyright Copyright  (C)  2021 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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

#define MAX_STACK_SIZE (64)

/**
  * Non-recursive stack version of Quicksort
  *
  *   ... from N. Wirth's Pascal Book, 'Algorithms + Data Structures = Programms'.
  *       by Alan Miller ( 19 Jul 1995 )
  *
  *  taken from:
  * - http://www.nag.com/nagware/examples.asp
  * - http://www.nag.com/nagware/Examples/nur.f90
  *
  *  see also:
  * - http://en.wikipedia.org/wiki/Quicksort
  *
  * @param[in,out]  a data to be sorted
  * @param[in]      n length of data
  * @param[out]     b arrays sorting along with a
  * @param[out]     c arrays sorting along with a
  *
 **/
void YAC_QSORT_NAME(
  YAC_QSORT_ARRAY_TYPE_A * a, size_t n,
  YAC_QSORT_ARRAY_TYPE_B * b, YAC_QSORT_ARRAY_TYPE_C * c) {

  if (n < 2) return;

  /*  Initialization */

  struct {
    size_t l, r;
  } stack[MAX_STACK_SIZE];

  size_t stack_size = 1;
  stack[0].l = 0;
  stack[0].r = n-1;

  /*  Start sorting
            
  ... keep taking the top request from the stack until stack_size = 0.

  */

  while  ( stack_size > 0 ) {

    stack_size--;
    size_t l = stack[stack_size].l;
    size_t r = stack[stack_size].r;

    /* ... keep splitting a[l], ... ,a[r] until l>= r. */

    while ( l < r ) {

      size_t i = l;
      size_t j = r;
      YAC_QSORT_ARRAY_TYPE_A x = a[(l+r+1) / 2];

      do {

        /* Search from lower end */
        while(a[i] < x) i++;

        /* Search from upper end */
        while(x < a[j]) j--;

        if (i > j) break;

        /* Swap positions i & j */

        YAC_QSORT_ARRAY_TYPE_A temp_a = a[i];
        a[i] = a[j];
        a[j] = temp_a;

        if (b != NULL) {
           YAC_QSORT_ARRAY_TYPE_B temp_b = b[i];
           b[i] = b[j];
           b[j] = temp_b;
        }
        if (c != NULL) {
           YAC_QSORT_ARRAY_TYPE_C temp_c = c[i];
           c[i] = c[j];
           c[j] = temp_c;
        }

        i++;
        j--;

      } while (i < j);

      if ( j-l >= r-i ) {

        if ( l < j ) {

          stack[stack_size].l = l;
          stack[stack_size].r = j;
          stack_size++;
        }

        l = i;

      } else {

        if ( i < r ) {

          stack[stack_size].l = i;
          stack[stack_size].r = r;
          stack_size++;
        }
        r = j;
      }
    } /* ( l < r ) */
  } /* ( stack_size /= 0 ) */
}
#undef MAX_STACK_SIZE
