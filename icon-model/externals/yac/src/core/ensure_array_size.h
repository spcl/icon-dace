// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ENSURE_ARRAY_SIZE_H
#define ENSURE_ARRAY_SIZE_H

#include <stdlib.h>

void
yac_realloc_array(void **array, size_t elem_size, size_t *curr_array_size,
                  size_t requested_size);

#define ENSURE_ARRAY_SIZE(arrayp, curr_array_size, req_size)            \
  do {                                                                  \
    if ((size_t)(req_size) > (size_t)(curr_array_size))                 \
    {                                                                   \
      size_t casize = (curr_array_size);                                \
                                                                        \
      yac_realloc_array((void **)&(arrayp), sizeof(*(arrayp)), &casize, \
                        (req_size));                                    \
      (curr_array_size) = casize;                                       \
    }                                                                   \
  }                                                                     \
  while(0)

#define ENSURE_BYTE_ARRAY_SIZE(arrayp, curr_array_size, req_size) \
  do {                                                            \
    if ((size_t)(req_size) > (size_t)(curr_array_size))           \
    {                                                             \
       size_t casize = (curr_array_size);                         \
                                                                  \
       yac_realloc_array(&(arrayp), 1, &casize, (req_size));      \
                                                                  \
       (curr_array_size) = casize;                                \
    }                                                             \
  }                                                               \
  while(0)

#endif

