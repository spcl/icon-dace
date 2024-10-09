// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tests.h"
#include "utils_common.h"

int main (int argc, char** argv) {

  if (argc != 2) PUT_ERR("wrong number of arguments");

  enum {
    SRC_COUNT = 10,
    TGT_COUNT = 5,
  };
  double weights[TGT_COUNT][SRC_COUNT] =
    {{0.0,0.1,0.0,0.3,0.0,0.0,0.2,0.0,0.4,0.0},
     {0.0,0.0,0.0,0.0,0.5,0.5,0.0,0.0,0.0,0.0},
     {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1},
     {0.2,0.0,0.2,0.0,0.2,0.0,0.2,0.0,0.2,0.0},
     {0.0,0.3,0.0,0.0,0.3,0.0,0.0,0.4,0.0,0.0}};

  // generate source field
  double src_field[SRC_COUNT];
  for (size_t i = 0; i < SRC_COUNT; ++i) src_field[i] = (double)(i+1);

  // compute reference target field
  double ref_tgt_field[TGT_COUNT];
  for (size_t i = 0; i < TGT_COUNT; ++i) {
    double tgt_value = 0.0;
    for (size_t j = 0; j < SRC_COUNT; ++j)
      tgt_value += src_field[j] * weights[i][j];
    ref_tgt_field[i] = tgt_value;
  }

#ifdef YAC_OPENMP_ENABLED
  int ref_num_threads = atoi(argv[1]);
#else
  UNUSED(argv);
#endif

  size_t prefix_num_src_per_tgt[TGT_COUNT + 1];
  size_t src_idx[SRC_COUNT * TGT_COUNT];
  double compact_weights[SRC_COUNT * TGT_COUNT];

  // generate sparse weight matrix
  prefix_num_src_per_tgt[0] = 0;
  for (size_t i = 0, k = 0; i < TGT_COUNT; ++i) {
    for (size_t j = 0; j < SRC_COUNT; ++j) {
      if (weights[i][j] != 0.0) {
        src_idx[k] = j;
        compact_weights[k] = weights[i][j];
        ++k;
      }
    }
    prefix_num_src_per_tgt[i+1] = k;
  }

  // compute target field using sparse weight matrix and OpenMP (if available)
  double tgt_field[TGT_COUNT];
  YAC_OMP_PARALLEL
  {
    YAC_OMP_FOR
    for (size_t i = 0; i < TGT_COUNT; ++i) {
#ifdef YAC_OPENMP_ENABLED
      if (omp_get_num_threads() != ref_num_threads)
        PUT_ERR("wrong number of threads");
#endif
      double tgt_value = 0.0;
      size_t const j_bound = prefix_num_src_per_tgt[i+1];
      for (size_t j = prefix_num_src_per_tgt[i]; j < j_bound; ++j) {
        tgt_value += src_field[src_idx[j]] * compact_weights[j];
      }
      tgt_field[i] = tgt_value;
    }
  }

  // check results
  for (size_t i = 0; i < TGT_COUNT; ++i)
    if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1.0e-6)
      PUT_ERR("wrong tgt_field value");

  return TEST_EXIT_CODE;
}
