/**
 * @file test_quicksort.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include <stdio.h>
#include "utils.h"
#include "tests.h"


int main (void) {

  { // descending order
#define LEN 128
    int a[LEN];
    int idx[LEN];


    for (int i = 0; i < LEN; ++i) {
      a[i]   = LEN-i-1;
      idx[i] = i;
    }

#ifdef VERBOSE
    for (int i = 0; i < LEN; ++i) {
      printf ("Unsorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
    }
#endif

    yac_quicksort_index ( a, (size_t)LEN, idx );

    for (int i = 0; i < LEN; ++i) {

      if ((i != a[i]) || (idx[i] != LEN-i-1)) INC_ERR;

#ifdef VERBOSE
      printf ("Sorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
#endif
    }
#undef LEN
  }

  { // ascending order
#define LEN 128
    int a[LEN];
    int idx[LEN];

    for (int i = 0; i < LEN; ++i) {
      a[i] = i;
      idx[i] = i;
    }

#ifdef VERBOSE
    for (int i = 0; i < LEN; ++i) {
      printf ("Unsorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
    }
#endif

    yac_quicksort_index ( a, (size_t)LEN, idx );

    for (int i = 0; i < LEN; ++i) {

      if ((i != a[i]) || (idx[i] != i)) INC_ERR;

#ifdef VERBOSE
      printf ("Sorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
#endif
    }
#undef LEN
  }

  return TEST_EXIT_CODE;
}
