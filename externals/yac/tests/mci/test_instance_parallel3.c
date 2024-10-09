// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <mpi.h>
#include <yaxt.h>

#include "tests.h"
#include "test_common.h"
#include "instance.h"
#include "yac.h"
#include "yac_mpi.h"
#include "event.h"

char * str_logical[2] = {"true", "false"};

int main (void) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size != 9) {

    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }

  struct yac_instance * instance =
    yac_instance_new(MPI_COMM_WORLD);

  // generate coupling configuration
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_instance_def_datetime(
    instance, "2008-03-09T16:05:07", "2008-03-10T16:05:07");
  struct yac_interp_stack_config * interp_stack_config =
    yac_interp_stack_config_new();
  yac_interp_stack_config_add_fixed(interp_stack_config, -1.0);
  char const * coupling_timestep = yac_time_to_ISO("60", C_SECOND);
  int mapping_side[3][3] = {{-1, 0, 0},
                            { 0,-1, 1},
                            { 1, 1,-1}};
  for (size_t src_comp_idx = 0; src_comp_idx < 3; ++src_comp_idx) {
    char src_comp_name[8], src_grid_name[8];
    sprintf(src_comp_name, "comp_%zu", src_comp_idx + 1);
    sprintf(src_grid_name, "grid_%zu", src_comp_idx + 1);
    for (size_t tgt_comp_idx = 0; tgt_comp_idx < 3; ++tgt_comp_idx) {
      char tgt_comp_name[8], tgt_grid_name[8];
      sprintf(tgt_comp_name, "comp_%zu", tgt_comp_idx + 1);
      sprintf(tgt_grid_name, "grid_%zu", tgt_comp_idx + 1);
      if (src_comp_idx == tgt_comp_idx) continue;
      char field_name[32];
      sprintf(field_name, "field_from_%zu_to_%zu",
              src_comp_idx + 1, tgt_comp_idx + 1);
      yac_instance_def_couple(
        instance, src_comp_name, src_grid_name, field_name,
        tgt_comp_name, tgt_grid_name, field_name, coupling_timestep,
        YAC_REDUCTION_TIME_ACCUMULATE, interp_stack_config, 0, 0, NULL,
        mapping_side[src_comp_idx][tgt_comp_idx], 1.0, 0.0, 0, NULL, NULL);
    }
  }
  yac_interp_stack_config_delete(interp_stack_config);

  char * component_names[3] = {"comp_1", "comp_2", "comp_3"};
  char const * local_component_names[2];
  int comp_ranks[2];
  int comp_idx[3][2] = {{0,1},{0,2},{1,2}};

  for (unsigned i = 0; i < 2; ++i)
    local_component_names[i] = component_names[comp_idx[rank/3][i]];

  yac_instance_def_components(instance, local_component_names, 2);

  for (unsigned i = 0; i < 2; ++i) {
    MPI_Comm component_comm =
      yac_instance_get_comps_comm(
        instance, &(local_component_names[i]), 1);
    MPI_Comm_rank(component_comm, comp_ranks + i);
    MPI_Comm_free(&component_comm);
  }

  struct yac_basic_grid * grid[2];
  size_t num_vertices[2] = {2,2};
  int cyclic[2] = {0,0};
  double coordinates_x[2][2];
  double coordinates_y[2] = {0,1};
  yac_int global_cell_ids[2][1];
  yac_int global_corner_ids[2][4];
  yac_int global_edge_ids[2][4];
  int cell_core_mask = {1};
  int corner_core_mask[4] = {1,1,1,1};
  int edge_core_mask[4] = {1,1,1,1};
  int ref_global_edge_ids[6][4] = {{0,3,13,1},{2,5,14,3},{4,7,15,5},
                                    {6,9,16,7},{8,11,17,9},{10,12,18,11}};

  for (unsigned i = 0; i < 2; ++i) {
    char grid_name[1024];
    sprintf(grid_name, "grid_%d", comp_idx[rank/3][i] + 1);
    coordinates_x[i][0] = comp_ranks[i];
    coordinates_x[i][1] = comp_ranks[i] + 1;
    global_cell_ids[i][0] = comp_ranks[i];
    global_corner_ids[i][0] = comp_ranks[i] + 0;
    global_corner_ids[i][1] = comp_ranks[i] + 1;
    global_corner_ids[i][2] = comp_ranks[i] + 7;
    global_corner_ids[i][3] = comp_ranks[i] + 8;
    global_edge_ids[i][0] = ref_global_edge_ids[comp_ranks[i]][0];
    global_edge_ids[i][1] = ref_global_edge_ids[comp_ranks[i]][1];
    global_edge_ids[i][2] = ref_global_edge_ids[comp_ranks[i]][2];
    global_edge_ids[i][3] = ref_global_edge_ids[comp_ranks[i]][3];
    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg_2d_deg(
        num_vertices, cyclic, coordinates_x[i], coordinates_y);
    grid_data.cell_ids = TO_POINTER(global_cell_ids[i]);
    grid_data.vertex_ids = TO_POINTER(global_corner_ids[i]);
    grid_data.edge_ids = TO_POINTER(global_edge_ids[i]);
    grid_data.core_cell_mask = TO_POINTER(cell_core_mask);
    grid_data.core_vertex_mask = TO_POINTER(corner_core_mask);
    grid_data.core_edge_mask = TO_POINTER(edge_core_mask);
    grid[i] = yac_basic_grid_new(grid_name, grid_data);
    yac_basic_grid_add_coordinates_nocpy(grid[i], YAC_LOC_CORNER, NULL);
  }

  struct yac_interp_field interp_fields[2];

  for (unsigned i = 0; i < 2; ++i) {
    interp_fields[i].location = YAC_LOC_CORNER;
    interp_fields[i].coordinates_idx = 0;
    interp_fields[i].masks_idx = SIZE_MAX;
  }

  char field_name[2][2][2][1024];

  for (unsigned i = 0; i < 2; ++i) {
    int curr_comp_idx = comp_idx[rank/3][i];
    for (unsigned j = 0; j < 2; ++j) {
      int remote_comp_idx = (3 + curr_comp_idx + j + 1) % 3;
      sprintf(field_name[i][j][0], "field_from_%d_to_%d",
              curr_comp_idx + 1, remote_comp_idx + 1);
      sprintf(field_name[i][j][1], "field_from_%d_to_%d",
              remote_comp_idx + 1, curr_comp_idx + 1);
    }
  }

  char const * timestep = yac_time_to_ISO("10", C_SECOND);
  for (unsigned i = 0; i < 2; ++i)
    for (unsigned j = 0; j < 2; ++j)
      for (unsigned k = 0; k < 2; ++k)
        yac_instance_add_field(
          instance, field_name[i][j][k], local_component_names[i],
          grid[i], &interp_fields[i], 1, 1, timestep);

  yac_instance_setup(instance, grid, sizeof(grid) / sizeof(grid[0]));

  yac_instance_delete(instance);

  for (unsigned i = 0; i < 2; ++i) yac_basic_grid_delete(grid[i]);

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
