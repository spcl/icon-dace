// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>

#include "tests.h"
#include "utils_core.h"

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

