// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <limits.h>
#include <mpi.h>
#include <yaxt.h>

#include "tests.h"
#include "dist_grid_utils.h"
#include "basic_grid.h"
#include "interp_weights_internal.h"
#include "duplicate_stencils.h"

static struct yac_basic_grid_data
  generate_dummy_grid_data(
    size_t local_num_cells, size_t global_num_cells, size_t num_cells_offset);

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_size != 3) {
    PUT_ERR("This test requires 3 processes\n");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  // generate dummy grid
  struct yac_basic_grid_data grid_data =
    generate_dummy_grid_data(5, 15, comm_rank * 5);
  struct yac_basic_grid * grid = yac_basic_grid_new("grid", grid_data);

  // generate weights
  struct yac_interp_weights * weights =
    yac_interp_weights_new(
      MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

  struct remote_point tgts_data[5];
  struct remote_points tgts = {.data = tgts_data};
  struct remote_point srcs[5];

  yac_int tgts_global_ids[3][5] = {{0,1,2,5,6},{7,8,9},{-1}};
  int tgts_ranks[3][5] = {{0,0,0,1,1},{1,1,1},{-1}};
  uint64_t tgts_orig_poses[3][5] = {{0,1,2,0,1},{2,3,4},{UINT64_MAX}};
  int srcs_ranks[3][5] = {{1,1,1,1,2},{2,2,2},{-1}};
  uint64_t srcs_orig_poses[3][5] = {{0,1,2,3,0},{1,2,3},{UINT64_MAX}};
  size_t num_tgts[3] = {5, 3, 0};

  // generate initial weights
  tgts.count = num_tgts[comm_rank];
  for (size_t i = 0; i < tgts.count; ++i) {
    tgts.data[i].global_id = tgts_global_ids[comm_rank][i];
    tgts.data[i].data.count = 1;
    tgts.data[i].data.data.single.rank = tgts_ranks[comm_rank][i];
    tgts.data[i].data.data.single.orig_pos = tgts_orig_poses[comm_rank][i];
    srcs[i].data.count = 1;
    srcs[i].data.data.single.rank = srcs_ranks[comm_rank][i];
    srcs[i].data.data.single.orig_pos = srcs_orig_poses[comm_rank][i];
  }
  yac_interp_weights_add_direct(weights, &tgts, srcs);

  // duplicate some stencils
  yac_int tgt_orig_global_id[3][4] = {{8,9},{-1},{0,1,8,9}};
  size_t tgt_duplicated_idx[3][4] = {{3,4},{SIZE_MAX},{0,1,3,4}};
  size_t nbr_duplicated[3] = {2,0,4};
  yac_duplicate_stencils(
    weights, grid, tgt_orig_global_id[comm_rank],
    tgt_duplicated_idx[comm_rank], nbr_duplicated[comm_rank], YAC_LOC_CELL);

  // get interpolation from the weights
  struct yac_interpolation * interpolation =
    yac_interp_weights_get_interpolation(
      weights, YAC_MAPPING_ON_SRC, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

  // execute interpolation
  double src_field_data[3][4] =
    {{-1, -1, -1, -1}, {0, 1, 2, 3}, {4, 5, 6, 7}};
  double ** src_fields[1] = {(double*[1]){src_field_data[comm_rank]}};
  double * tgt_field[1] = {(double[5]){-1,-1,-1,-1,-1}};
  yac_interpolation_execute(interpolation, src_fields, tgt_field);

  // check results
  double ref_tgt_field[3][5] =
    {{0,1,2,6,7}, {3,4,5,6,7}, {0,1,-1,6,7}};
  for (size_t i = 0; i < 5; ++i)
    if (ref_tgt_field[comm_rank][i] != tgt_field[0][i])
      PUT_ERR("ERROR in yac_duplicate_stencils");

  // cleanup
  yac_interpolation_delete(interpolation);
  yac_interp_weights_delete(weights);
  yac_basic_grid_delete(grid);

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static struct yac_basic_grid_data
  generate_dummy_grid_data(
    size_t local_num_cells, size_t global_num_cells, size_t num_cells_offset) {

  double coordinates_x[global_num_cells + 1];
  double coordinates_y[2] = {0.0, 1.0};
  size_t num_cells[2] = {global_num_cells, 1};
  size_t local_start[2] = {num_cells_offset, 0};
  size_t local_count[2] = {local_num_cells, 1};

  for (size_t i = 0; i <= global_num_cells; ++i) coordinates_x[i] = (double)i;

  return
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells, local_start, local_count, 0);
}
