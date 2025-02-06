// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <math.h>
#include "generate_reg2d.h"

void yac_generate_reg2d_decomp(
  int num_points[2], int total_num_procs, int * num_procs) {

  int flag = num_points[0] > num_points[1];
  float ratio = (float)(num_points[flag^1]) / (float)(num_points[flag]);

  int prev_num_procs[2];

  prev_num_procs[!flag] = total_num_procs;
  prev_num_procs[flag] = 1;

  float prev_proc_ratio = (float)total_num_procs;

  for (int i = total_num_procs-1; i > 0; --i) {

    if (total_num_procs%i == 0) {

      float curr_proc_ratio = (float)(i) / (float)(total_num_procs/i);

      double temp_a = fabs(ratio - curr_proc_ratio);
      double temp_b = fabs(ratio - prev_proc_ratio);
      if (temp_a < temp_b) {
        prev_num_procs[!flag] = i;
        prev_num_procs[flag] = total_num_procs/i;
        prev_proc_ratio = curr_proc_ratio;
      }
    }
  }

  num_procs[0] = prev_num_procs[0];
  num_procs[1] = prev_num_procs[1];
}

