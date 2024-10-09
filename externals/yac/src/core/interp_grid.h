// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_GRID_H
#define INTERP_GRID_H

#include "dist_grid.h"

// YAC PUBLIC HEADER START

struct yac_interp_grid;

/** \example test_interp_grid_parallel.c
 * A test for parallel grid interpolation.
 */

/**
 * generate a interpolation grid
 * @param[in] grid_pair      distributed grid pair
 * @param[in] src_grid_name  name of the source grid
 * @param[in] tgt_grid_name  name of the target grid
 * @param[in] num_src_fields number of source fields
 * @param[in] src_fields     specifies the source fields
 * @param[in] tgt_field      specifies the target field
 * @return interpolation grid
 */
struct yac_interp_grid * yac_interp_grid_new(
  struct yac_dist_grid_pair * grid_pair,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t num_src_fields, struct yac_interp_field const * src_fields,
  struct yac_interp_field const tgt_field);

/**
 * deletes all memory associated with the provided interpolation grid
 * @param[inout] interp_grid interpolation grid
 */
void yac_interp_grid_delete(struct yac_interp_grid * interp_grid);

// YAC PUBLIC HEADER STOP

#endif // INTERP_GRID_H
