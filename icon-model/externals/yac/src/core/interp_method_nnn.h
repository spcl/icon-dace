// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_NNN_H
#define INTERP_METHOD_NNN_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

/** \example test_interp_method_nnn_parallel.c
 * A test for the parallel nearest neighbour interpolation method.
 */

/** \example test_interp_method_rbf_parallel.c
 * A test for the parallel radial basis function interpolation method.
 */

enum yac_interp_nnn_weight_type {
  YAC_INTERP_NNN_AVG   = 0, //!< average of n source points
  YAC_INTERP_NNN_DIST  = 1, //!< distance weighted average of n source points
  YAC_INTERP_NNN_GAUSS = 2, //!< distance with Gauss weights of n source points
  YAC_INTERP_NNN_RBF   = 3, //!< radial basis functions
  YAC_INTERP_NNN_ZERO  = 4, //!< all weights are set to zero
};

#define YAC_INTERP_NNN_WEIGHTED_DEFAULT (0)
#define YAC_INTERP_NNN_N_DEFAULT (1)
#define YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT (0.0)
#define YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT (0.1)

#define YAC_INTERP_RBF_N_DEFAULT (9)
#define YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT (0.0)
#define YAC_INTERP_RBF_SCALE_DEFAULT (1.487973e+01)
#define YAC_INTERP_RBF_KERNEL_DEFAULT (0)

struct yac_nnn_config {
  enum yac_interp_nnn_weight_type type;
  size_t n;
  double max_search_distance;
  union {
    double rbf_scale;
    double gauss_scale;
  } data;
};

struct interp_method * yac_interp_method_nnn_new(struct yac_nnn_config config);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_NNN_H
