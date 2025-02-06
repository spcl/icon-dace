/**
 * @file test_sort.c
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

#include <assert.h>

//#define VERBOSE

#include <yaxt.h>

#include "tests.h"

static void
test1(void (*anysort_index)(Xt_int * a, int n, int * idx, int reset_index)) {
  Xt_int ivec[]        = {7,5,3,1,2,2,3,3};
  static const Xt_int sorted_ivec[] = {1,2,2,3,3,3,5,7};
  int pvec[]           = {1,2,3,4,5,6,7,8};
  static const int sorted_pvec[]    = {3,4,5,2,6,7,1,0};

  size_t tsize = sizeof(ivec) / sizeof(ivec[0]);
  assert(tsize == sizeof(pvec) / sizeof(pvec[0]));

  anysort_index(ivec, (int)tsize, pvec, 1);
  bool mismatch_values = false, mismatch_positions = false;
  for(size_t i=0; i < tsize; i++) {
    mismatch_values |= ( ivec[i] != sorted_ivec[i] );
    mismatch_positions |= ( pvec[i] != sorted_pvec[i] );
  }
  if ( mismatch_values )
    PUT_ERR("wrong sorting values\n");
  if ( mismatch_positions )
      PUT_ERR("wrong sorting positions\n");
}

static void
test2(void (*sort_int)(int *a, size_t n))
{
  int ivec[] = {7,5,3,1,2,2,3,3};
  static const int sorted_ivec[] = {1,2,2,3,3,3,5,7};
  size_t ivec_size = sizeof (ivec) / sizeof (ivec[0]);
  sort_int(ivec, ivec_size);
  bool mismatch = false;
  for(size_t i=0; i < ivec_size; i++)
    mismatch |= ( ivec[i] != sorted_ivec[i] );
  if (mismatch)
    PUT_ERR("wrong sorting values\n");
  size_t n = 40000;
  int *large_ivec = malloc(n * sizeof (*large_ivec));
  for (size_t i = 0; i < n; ++i)
    large_ivec[i] = (int)(n - i);
  sort_int(large_ivec, n);
  mismatch = false;
  for (size_t i = 0; i < n; ++i)
    mismatch |= (large_ivec[i] != (int)(i + 1));
  if (mismatch)
    PUT_ERR("wrong sorting values\n");
  size_t ofs = 50;
  for (size_t i = 0; i < n; ++i)
    large_ivec[i] = (int)(((i + ofs)%n) / 2);
  sort_int(large_ivec, n);
  mismatch = false;
  for (size_t i = 0; i < n; ++i)
    mismatch |= (large_ivec[i] != (int)(i/2));
  if (mismatch)
    PUT_ERR("wrong sorting values\n");
  free(large_ivec);
}


int main(void) {

  test1(xt_quicksort_index);
  test1(xt_mergesort_index);
  test1(xt_sort_index);

  test2(xt_sort_int);

  return TEST_EXIT_CODE;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
