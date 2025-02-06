// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <stdlib.h>

#include "ensure_array_size.h"
#include "utils_core.h"

void
yac_realloc_array(void **array, size_t elem_size, size_t *curr_array_size,
                  size_t requested_size)
{
  const size_t array_inc_size = (1024 + elem_size - 1)/ elem_size;
  *curr_array_size = array_inc_size
    * ((requested_size + array_inc_size) / array_inc_size);
  *array = xrealloc(*array, *curr_array_size * elem_size);
  if ((*array == NULL) && ((*curr_array_size * elem_size) != 0)) {
    // GCOVR_EXCL_START
    fprintf(stderr, "error in realloc\n");
    exit(EXIT_FAILURE);
    // GCOVR_EXCL_STOP
  }
}

