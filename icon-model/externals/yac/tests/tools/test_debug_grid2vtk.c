// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <unistd.h> 
#include <stdlib.h>
#include <libgen.h>
#include <string.h>
#include <math.h>

#include "tests.h"
#include "yac_core.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

#define TO_POINTER(a) (to_pointer(&a, sizeof(a)))

static void * to_pointer(void * data, size_t size_data) {

  void * ret_value = malloc(size_data);
  memcpy(ret_value, data, size_data);
  return ret_value;
}

static void LLtoXYZ(double lon, double lat, double p_out[]) {

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

int main (int argc, char* argv[]) {

  if (argc != 2) {
    PUT_ERR("ERROR: wrong number of arguments");
    return TEST_EXIT_CODE;
  }

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int comm_size, comm_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if ((comm_size != 2) && (comm_rank == 0)) {
    fputs("wrong number of processes (has to be four)\n", stderr);
    exit(EXIT_FAILURE);
  }

  char const * debug_grid2vtk_exec = argv[1];

  { // test using artifical grid

    char const * grid_filename = "test_debug_grid.nc";
    char const * vtk_filename = "test_debug_grid.vtk";
    char const * grid_name = "dummy_grid";

    // write dummy grid file
    enum {NUM_VERTICES = 9, NUM_CELLS = 2, NUM_PROCS = 2};
    int num_cells[NUM_PROCS] = {2,1};
    int num_vertices_per_cell[NUM_PROCS][NUM_CELLS] = {{3,3},{4}};
    double x_vertices[NUM_VERTICES] = {0,1,2, 0,1,2, 0,1,2};
    double y_vertices[NUM_VERTICES] = {0,0,0, 1,1,1, 2,2,2};
    int cell_to_vertex[NUM_PROCS][6] = {{1,2,5, 1,5,4}, {3,4,7,6}};
    int core_cell_mask[NUM_PROCS][NUM_CELLS] = {{1,1}, {0}};
    double x_cells[NUM_PROCS][NUM_CELLS] = {{1.75, 1.25},{0.5}};
    double y_cells[NUM_PROCS][NUM_CELLS] = {{0.25, 0.75},{1.5}};
    yac_int cell_ids[NUM_PROCS][NUM_CELLS] = {{0,1}, {2}};
    double cell_centers[NUM_CELLS][3];
    for (int i = 0; i < num_cells[comm_rank]; ++i)
      LLtoXYZ(
        x_cells[comm_rank][i]*YAC_RAD, y_cells[comm_rank][i]*YAC_RAD,
        cell_centers[i]);


    for (int with_gid = 0; with_gid < 2; ++with_gid) {
      for (int with_cmk = 0; with_cmk < 2; ++with_cmk) {
        for (int with_center = 0; with_center < 2; ++with_center) {

          struct yac_basic_grid * dummy_grid =
            yac_basic_grid_unstruct_deg_new(
              grid_name, NUM_VERTICES, num_cells[comm_rank],
              num_vertices_per_cell[comm_rank],
              x_vertices, y_vertices, cell_to_vertex[comm_rank]);

          struct yac_basic_grid_data * dummy_grid_data =
            yac_basic_grid_get_data(dummy_grid);

          if (with_gid)
            dummy_grid_data->cell_ids =
              TO_POINTER(cell_ids[comm_rank]);
          if (with_cmk)
            dummy_grid_data->core_cell_mask =
              TO_POINTER(core_cell_mask[comm_rank]);
          if (with_center)
            yac_basic_grid_add_coordinates(
              dummy_grid, YAC_LOC_CELL, cell_centers, num_cells[comm_rank]);

          yac_basic_grid_to_file_parallel(
            dummy_grid, grid_filename, MPI_COMM_WORLD);

          yac_basic_grid_delete(dummy_grid);

          if (comm_rank == 0) {
            char cmd[1024];
            sprintf(
              cmd, "%s -f %s -n %s -o %s",
              debug_grid2vtk_exec, grid_filename, grid_name, vtk_filename);
            if (system(cmd)) PUT_ERR("failed to execute");
          }

          // ensure that all processes finished reading the file
          MPI_Barrier(MPI_COMM_WORLD);
          if (comm_rank == 0) {
            unlink(grid_filename);
            unlink(vtk_filename);
          }
        }
      }
    }
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
