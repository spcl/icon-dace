// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_FILE_H
#define INTERP_METHOD_FILE_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

/** \example test_interp_method_file_parallel.c
 * A test for the parallel file interpolation method.
 */

#define YAC_WEIGHT_FILE_VERSION_STRING "yac weight file 1.0"

#define YAC_INTERP_FILE_WEIGHT_FILE_NAME_DEFAULT (NULL)

struct interp_method * yac_interp_method_file_new(
  char const * weight_file_name);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_FILE_H
