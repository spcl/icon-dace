// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <mpi.h>

#include "utils_common.h"
#include "grid_cell.h"
#include "basic_grid_data.h"

/** \example test_common.c
 * A test for test-internal functions - not relevant for the external user.
 */

#define TO_POINTER(a) (to_pointer(&a, sizeof(a)))

struct yac_grid_cell generate_cell_deg(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners);
struct yac_grid_cell generate_cell_rad(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners);
struct yac_grid_cell generate_cell_3d(
  yac_coordinate_pointer coords, enum yac_edge_type * edge_type,
  size_t num_corners);
int intersect(enum yac_edge_type edge_type_a,
              double lon_a, double lat_a, double lon_b, double lat_b,
              enum yac_edge_type edge_type_b,
              double lon_c, double lat_c, double lon_d, double lat_d,
              double * intersection);
void * to_pointer(void * data, size_t size_data);
int double_are_equal(double a, double b);
int double_are_unequal(double a, double b);
void set_even_io_rank_list(MPI_Comm comm);
void clear_yac_io_env();
void check_basic_grid_data(
  struct yac_basic_grid_data grid_a, struct yac_basic_grid_data grid_b,
  char const * grid_name);

#endif // TEST_COMMON_H

