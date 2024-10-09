// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_AVG_H
#define INTERP_METHOD_AVG_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

/** \example test_interp_method_avg_parallel.c
 * A test for the parallel average interpolation method.
 */

enum yac_interp_avg_weight_type {
  YAC_INTERP_AVG_ARITHMETIC  = 0,  // simple average
  YAC_INTERP_AVG_DIST = 1, // distance weighted
  YAC_INTERP_AVG_BARY = 2, // barycentric coordinate based
};

#define YAC_INTERP_AVG_WEIGHT_TYPE_DEFAULT (0)
#define YAC_INTERP_AVG_PARTIAL_COVERAGE_DEFAULT (0)

struct interp_method * yac_interp_method_avg_new(
  enum yac_interp_avg_weight_type weight_type,
  int partial_coverage);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_AVG_H
