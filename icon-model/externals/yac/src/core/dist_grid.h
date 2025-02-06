// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DIST_GRID_H
#define DIST_GRID_H

#include <mpi.h>

#include "basic_grid.h"

// YAC PUBLIC HEADER START

/** \example test_dist_grid_pair_parallel.c
 * A test for distributed grid pairs.
 */

struct yac_dist_grid_pair;

/**
 * generate a distributed grid pair
 * @param[in] grid_a      first basic grid
 * @param[in] grid_b      second basic grid
 * @param[in] comm        communicator containing the ranks of all processes
 *                        participating in this call
 * @return distributed grid pair
 * @remark this routine is collective for all processes in comm
 */
struct yac_dist_grid_pair * yac_dist_grid_pair_new(
  struct yac_basic_grid * grid_a, struct yac_basic_grid * grid_b,
  MPI_Comm comm);

/**
 * delete all memory associated to the provided grid pair
 * @param[inout] grid_pair
 */
void yac_dist_grid_pair_delete(struct yac_dist_grid_pair * grid_pair);

// YAC PUBLIC HEADER STOP

#endif // DIST_GRID_H
