// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DIST_GRID_UTILS_H
#define DIST_GRID_UTILS_H

#include "basic_grid.h"

/** \example test_dist_grid_utils.c
 * A test for some basic tools used to generate distributed grids.
 */

struct yac_basic_grid_data yac_generate_basic_grid_data_reg2d(
  double const * coordinates_x, double const * coordinates_y,
  size_t const num_cells[2], size_t const local_start[2],
  size_t const local_count[2], int with_halo);

struct yac_basic_grid * yac_generate_basic_grid_reg2d(
  char const * name, double const * coordinates_x,
  double const * coordinates_y, size_t const num_cells[2],
  size_t const local_start[2], size_t const local_count[2], int with_halo);

#endif // DIST_GRID_UTILS_H
