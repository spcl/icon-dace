// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include "utils_core.h"
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

    for (int i = 0; i < LEN; ++i) {
      printf ("Unsorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
    }

    yac_quicksort_index ( a, (size_t)LEN, idx );

    for (int i = 0; i < LEN; ++i) {

      if ((i != a[i]) || (idx[i] != LEN-i-1)) INC_ERR;

      printf ("Sorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
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

    for (int i = 0; i < LEN; ++i) {
      printf ("Unsorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
    }

    yac_quicksort_index ( a, (size_t)LEN, idx );

    for (int i = 0; i < LEN; ++i) {

      if ((i != a[i]) || (idx[i] != i)) INC_ERR;

      printf ("Sorted list is i=%i  a[i]=%i idx[i]=%i\n ", i, a[i], idx[i]);
    }
#undef LEN
  }

  return TEST_EXIT_CODE;
}
