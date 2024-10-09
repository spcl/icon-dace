// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_NCC_H
#define INTERP_METHOD_NCC_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

/** \example test_interp_method_ncc_parallel.c
 * A test for the parallel nearest corner cells interpolation method.
 */

enum yac_interp_ncc_weight_type {
  YAC_INTERP_NCC_AVG   = 0, //!< average of n source points
  YAC_INTERP_NCC_DIST  = 1, //!< distance weighted average of n source points
};

#define YAC_INTERP_NCC_WEIGHT_TYPE_DEFAULT (0)
#define YAC_INTERP_NCC_PARTIAL_COVERAGE_DEFAULT (0)

struct interp_method * yac_interp_method_ncc_new(
  enum yac_interp_ncc_weight_type weight_type,
  int partial_coverage);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_NCC_H
