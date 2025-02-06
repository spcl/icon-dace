// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <string.h>

#include "ensure_array_size.h"
#include "utils_core.h"

struct run {

  size_t count;
  int ordering;
};

static void determine_runs(char * base, size_t num, size_t size,
                           int (*compar)(const void *,const void *),
                           struct run ** runs, size_t * num_runs) {

  size_t runs_array_size = 0;
  *runs = NULL;
  *num_runs = 1;

  ENSURE_ARRAY_SIZE(*runs, runs_array_size, *num_runs);
  (*runs)->count = 1;
  (*runs)->ordering = -1;

  char * curr_element = base + size;
  int curr_order = -1;
  size_t * curr_count = &((*runs)->count);

  for (size_t i = 1; i < num; ++i) {

    int order = compar(curr_element - size, curr_element);

    if ((order == 0) || !((order < 0) ^ (curr_order < 0))) {

      ++*curr_count;

    } else {

      ENSURE_ARRAY_SIZE(*runs, runs_array_size, *num_runs + 1);

      curr_order = order;
      curr_count = &((*runs)[*num_runs].count);
      (*runs)[*num_runs].count = 1;
      (*runs)[*num_runs].ordering = order;
      ++*num_runs;
    }

    curr_element += size;
  }
}

static void merge(char * base_a, size_t num_a, int a_ascending,
                  char * base_b, size_t num_b, int b_ascending,
                  size_t size,
                  int (*compar)(const void *,const void *),
                  char * target) {

  char * curr_a = base_a + ((a_ascending)?0:((num_a - 1) * size));
  char * curr_b = base_b + ((b_ascending)?0:((num_b - 1) * size));
  int inc_a = (a_ascending)?(size):(-size);
  int inc_b = (b_ascending)?(size):(-size);
  char * end_a = curr_a + inc_a * num_a;
  char * end_b = curr_b + inc_b * num_b;

  while (curr_a != end_a && curr_b != end_b) {

    char * from;

    if (compar(curr_a, curr_b) <= 0) {
      from = curr_a;
      curr_a += inc_a;
    } else {
      from = curr_b;
      curr_b += inc_b;
    }

    memcpy(target, from, size);
    target += size;
  }

  size_t diff = (size_t)(end_a - curr_a);
  if (diff) {
    memcpy(target, curr_a, diff);
    target += diff;
  }

  diff = (size_t)(end_b - curr_b);
  if (diff) {
    memcpy(target, curr_b, diff);
    target += diff;
  }
}

static void merge_runs(
  char * base, size_t size, int (*compar)(const void *,const void *),
  struct run ** runs, size_t * num_runs, char * buffer) {

  size_t num_merges = *num_runs >> 1;

  struct run * curr_runs = *runs;

  for (size_t i = 0; i < num_merges; ++i) {

    merge(base, curr_runs[0].count, curr_runs[0].ordering < 0,
          base + curr_runs[0].count * size, curr_runs[1].count,
          curr_runs[1].ordering < 0, size, compar, buffer);

    base += (curr_runs[0].count + curr_runs[1].count) * size;
    buffer += (curr_runs[0].count + curr_runs[1].count) * size;

    (*runs)[i].count = curr_runs[0].count + curr_runs[1].count;
    (*runs)[i].ordering = -1;

    curr_runs += 2;
  }

  if ((*num_runs)  & 1) {
  
    (*runs)[num_merges] = (*runs)[*num_runs-1];
    memcpy(buffer, base, (*runs)[num_merges].count * size);
  }

  *num_runs = num_merges + ((*num_runs) & 1);
}

/**
  * Natural Merge sort
  *
 **/
void yac_mergesort(void * base, size_t num, size_t size,
               int (*compar)(const void *,const void *)) {

  if (size == 0 || num == 0) return;

  struct run * runs;
  size_t num_runs;

  char * buffer = xmalloc(num * size);

  determine_runs((char*)base, num, size, compar, &runs, &num_runs);

  while (num_runs > 1) {
    merge_runs((char*)base, size, compar, &runs, &num_runs, buffer);
    merge_runs(buffer, size, compar, &runs, &num_runs, (char*)base);
  };

  free(runs);
  free(buffer);
}

