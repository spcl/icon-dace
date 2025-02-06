// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

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
  * @param[in,out]  a            data to be sorted
  * @param[in]      n            length of data
  * @param[out]     idx          old index of sorted returned a
  *
 **/
void YAC_QSORT_NAME(
  YAC_QSORT_ARRAY_TYPE * a, size_t n, YAC_QSORT_INDEX_TYPE * idx) {

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
      YAC_QSORT_ARRAY_TYPE x = a[(l+r+1) / 2];

      do {

        /* Search from lower end */
        while(a[i] < x) i++;

        /* Search from upper end */
        while(x < a[j]) j--;

        if (i > j) break;

        /* Swap positions i & j */

        YAC_QSORT_ARRAY_TYPE temp_a = a[i];
        a[i] = a[j];
        a[j] = temp_a;

        if (idx != NULL) {
           YAC_QSORT_INDEX_TYPE temp_idx = idx[i];
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
}
#undef MAX_STACK_SIZE
