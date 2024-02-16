/**
 * @file quicksort.c
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
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

#include <stdlib.h>
#include <string.h>

#include "utils.h"

#define YAC_QSORT_ARRAY_TYPE int
#define YAC_QSORT_INDEX_TYPE int
#define YAC_QSORT_NAME yac_quicksort_index
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE yac_int
#define YAC_QSORT_INDEX_TYPE size_t
#define YAC_QSORT_NAME yac_quicksort_index_yac_int_size_t
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE yac_int
#define YAC_QSORT_INDEX_TYPE yac_int
#define YAC_QSORT_NAME yac_quicksort_index_yac_int_yac_int
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE yac_int
#define YAC_QSORT_INDEX_TYPE uint64_t
#define YAC_QSORT_NAME yac_quicksort_index_yac_int_uint64_t
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE yac_int
#define YAC_QSORT_INDEX_TYPE int
#define YAC_QSORT_NAME yac_quicksort_index_yac_int_int
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE size_t
#define YAC_QSORT_INDEX_TYPE size_t
#define YAC_QSORT_NAME yac_quicksort_index_size_t_size_t
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE uint64_t
#define YAC_QSORT_INDEX_TYPE size_t
#define YAC_QSORT_NAME yac_quicksort_index_uint64_t_size_t
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE int
#define YAC_QSORT_INDEX_TYPE size_t
#define YAC_QSORT_NAME yac_quicksort_index_int_size_t
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE size_t
#define YAC_QSORT_INDEX_TYPE yac_int
#define YAC_QSORT_NAME yac_quicksort_index_size_t_yac_int
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE size_t
#define YAC_QSORT_INDEX_TYPE int
#define YAC_QSORT_NAME yac_quicksort_index_size_t_int
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE size_t
#define YAC_QSORT_INDEX_TYPE void*
#define YAC_QSORT_NAME yac_quicksort_index_size_t_void_p
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE int
#define YAC_QSORT_INDEX_TYPE yac_int
#define YAC_QSORT_NAME yac_quicksort_index_int_yac_int
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE int
#define YAC_QSORT_INDEX_TYPE double
#define YAC_QSORT_NAME yac_quicksort_index_int_double
#include "quicksort_template.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_INDEX_TYPE
#undef YAC_QSORT_ARRAY_TYPE

#define YAC_QSORT_ARRAY_TYPE_A size_t
#define YAC_QSORT_ARRAY_TYPE_B size_t
#define YAC_QSORT_ARRAY_TYPE_C double
#define YAC_QSORT_NAME yac_quicksort_index_size_t_size_t_double
#include "quicksort_template_2.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_ARRAY_TYPE_A
#undef YAC_QSORT_ARRAY_TYPE_B
#undef YAC_QSORT_ARRAY_TYPE_C

#define YAC_QSORT_ARRAY_TYPE_A yac_int
#define YAC_QSORT_ARRAY_TYPE_B yac_int
#define YAC_QSORT_ARRAY_TYPE_C double
#define YAC_QSORT_NAME yac_quicksort_index_yac_int_yac_int_double
#include "quicksort_template_2.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_ARRAY_TYPE_A
#undef YAC_QSORT_ARRAY_TYPE_B
#undef YAC_QSORT_ARRAY_TYPE_C

#define YAC_QSORT_ARRAY_TYPE_A yac_int
#define YAC_QSORT_ARRAY_TYPE_B yac_int
#define YAC_QSORT_ARRAY_TYPE_C size_t
#define YAC_QSORT_NAME yac_quicksort_index_yac_int_yac_int_size_t
#include "quicksort_template_2.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_ARRAY_TYPE_A
#undef YAC_QSORT_ARRAY_TYPE_B
#undef YAC_QSORT_ARRAY_TYPE_C

#define YAC_QSORT_ARRAY_TYPE_A int
#define YAC_QSORT_ARRAY_TYPE_B size_t
#define YAC_QSORT_ARRAY_TYPE_C size_t
#define YAC_QSORT_NAME yac_quicksort_index_int_size_t_size_t
#include "quicksort_template_2.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_ARRAY_TYPE_A
#undef YAC_QSORT_ARRAY_TYPE_B
#undef YAC_QSORT_ARRAY_TYPE_C

#define YAC_QSORT_ARRAY_TYPE_A int
#define YAC_QSORT_ARRAY_TYPE_B size_t
#define YAC_QSORT_ARRAY_TYPE_C yac_int
#define YAC_QSORT_NAME yac_quicksort_index_int_size_t_yac_int
#include "quicksort_template_2.h"
#undef YAC_QSORT_NAME
#undef YAC_QSORT_ARRAY_TYPE_A
#undef YAC_QSORT_ARRAY_TYPE_B
#undef YAC_QSORT_ARRAY_TYPE_C

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
  * @param[in,out]  a_           data to be sorted
  * @param[in]      count        length of data
  * @param[in]      size         size of single data element
  * @param[in]      compare      compare routine for data elements
  * @param[out]     idx          old index of sorted returned a
  *
 **/
void yac_qsort_index(
  void * a_, size_t count, size_t size,
  int (*compare)(void const *, void const *), size_t * idx) {

  if (count < 2) return;

  /*  Initialization */

  struct {
    size_t l, r;
  } stack[MAX_STACK_SIZE];

  size_t stack_size = 1;
  stack[0].l = 0;
  stack[0].r = count-1;

  void * x = xmalloc(2 * size);
  void * temp = (unsigned char*)x + size;
  unsigned char * a = a_;

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
      memcpy(x, &a[((l+r+1) / 2)*size], size);

      do {

        /* Search from lower end */
        while(compare(&a[i*size], x) < 0) i++;

        /* Search from upper end */
        while(compare(x, &a[j*size]) < 0) j--;

        if (i > j) break;

        /* Swap positions i & j */

        memcpy(temp, &a[i*size], size);
        memcpy(&a[i*size], &a[j*size], size);
        memcpy(&a[j*size], temp, size);

        if (idx != NULL) {
           size_t temp_idx = idx[i];
           idx[i] = idx[j];
           idx[j] = temp_idx;
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

  free(x);
}
#undef MAX_STACK_SIZE
