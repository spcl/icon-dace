/**
 * @file test_mergesort.c
 *
 * @copyright Copyright  (C)  2014 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
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

struct test_struct {

  double dummy;
  int idx;
  double dummy_;
};

static int compare_test_struct(const void * a, const void * b) {

  int idx_a = ((struct test_struct *)a)->idx;
  int idx_b = ((struct test_struct *)b)->idx;

  return ((idx_a) > (idx_b)) - ((idx_a) < (idx_b));
}

int main (void) {

  {
    struct test_struct a[128];
    struct test_struct ref[128];
    size_t len = sizeof(a)/sizeof(a[0]);

    for (size_t i = 0; i < len; ++i) {
      a[i].dummy = -1;
      a[i].dummy_ = -1;
      ref[i].dummy = -1;
      ref[i].dummy_ = -1;
    }

    for (size_t i = 0; i < len; ++i) {
      a[i].idx = i;
      ref[i].idx = i;
    }

    qsort(ref, len, sizeof(ref[0]), compare_test_struct);
    yac_mergesort(a, len, sizeof(ref[0]), compare_test_struct);

    for (size_t i = 0; i < len; ++i) if (a[i].idx != ref[i].idx) INC_ERR;
  }

  {
    struct test_struct a[128];
    struct test_struct ref[128];

    size_t len = sizeof(a)/sizeof(a[0]);
    for (size_t i = 0; i < len; ++i) {
      a[i].dummy = -1;
      a[i].dummy_ = -1;
      ref[i].dummy = -1;
      ref[i].dummy_ = -1;
    }

    for (size_t i = 0; i < len; ++i) {
      a[i].idx = len - i;
      ref[i].idx = len - i;
    }

    qsort(ref, len, sizeof(ref[0]), compare_test_struct);
    yac_mergesort(a, len, sizeof(ref[0]), compare_test_struct);

    for (size_t i = 0; i < len; ++i) if (a[i].idx != ref[i].idx) INC_ERR;
  }

  {
    struct test_struct a[128];
    struct test_struct ref[128];

    size_t len = sizeof(a)/sizeof(a[0]);
    for (size_t i = 0; i < len; ++i) {
      a[i].dummy = -1;
      a[i].dummy_ = -1;
      ref[i].dummy = -1;
      ref[i].dummy_ = -1;
    }

    for (size_t i = 0; i < len; ++i) {
      a[i].idx = i;
      ref[i].idx = i;
    }
    for (size_t i = 32; i < 64; ++i) {
      a[i].idx = 63 - i;
      ref[i].idx = 63 - i;
    }
    for (size_t i = 96; i < 128; ++i) {
      a[i].idx = 127 - i;
      ref[i].idx = 127 - i;
    }

    qsort(ref, len, sizeof(ref[0]), compare_test_struct);
    yac_mergesort(a, len, sizeof(ref[0]), compare_test_struct);

    for (size_t i = 0; i < len; ++i) if (a[i].idx != ref[i].idx) INC_ERR;
  }

  {
    struct test_struct a[128];
    struct test_struct ref[128];

    size_t len = sizeof(a)/sizeof(a[0]);
    for (size_t i = 0; i < len; ++i) {
      a[i].dummy = -1;
      a[i].dummy_ = -1;
      ref[i].dummy = -1;
      ref[i].dummy_ = -1;
    }

    for (size_t i = 0; i < len; ++i) {
      a[i].idx = 127 - i;
      ref[i].idx = 127 - i;
    }
    for (size_t i = 32; i < 64; ++i) {
      a[i].idx = i;
      ref[i].idx = i;
    }
    for (size_t i = 96; i < 128; ++i) {
      a[i].idx = i;
      ref[i].idx = i;
    }

    qsort(ref, len, sizeof(ref[0]), compare_test_struct);
    yac_mergesort(a, len, sizeof(ref[0]), compare_test_struct);

    for (size_t i = 0; i < len; ++i) if (a[i].idx != ref[i].idx) INC_ERR;
  }

  for (int i = 0; i < 10; ++i) {

    struct test_struct a[128];
    struct test_struct ref[128];

    size_t len = sizeof(a)/sizeof(a[0]);
    for (size_t j = 0; j < len; ++j) {
      a[j].dummy = -1;
      a[j].dummy_ = -1;
      ref[j].dummy = -1;
      ref[j].dummy_ = -1;
    }

    for (size_t j = 0; j < len; ++j) ref[j].idx = a[j].idx = rand();

    qsort(ref, len, sizeof(ref[0]), compare_test_struct);
    yac_mergesort(a, len, sizeof(ref[0]), compare_test_struct);

    for (size_t j = 0; j < len; ++j) if (a[j].idx != ref[j].idx) INC_ERR;
  }

   return TEST_EXIT_CODE;
}

