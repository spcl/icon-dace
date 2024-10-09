// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>

#include "tests.h"
#include "test_common.h"
#include "basic_grid.h"
#include "io_utils.h"
#include "geometry.h"
#include "grid_file_common.h"

#include <mpi.h>
#include <netcdf.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

static void delete_file(char const * filename);

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  enum {NUM_PROCS = 3};
  if (comm_size != NUM_PROCS) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  {
    struct {
      size_t nbr_vertices[NUM_PROCS][2];
      double lon_vertices[NUM_PROCS][6];
      double lat_vertices[NUM_PROCS][6];
      yac_int cell_global_ids[NUM_PROCS][25];
      int core_cell_mask[NUM_PROCS][25];
      int has_vertex_data;
      yac_int vertex_global_ids[NUM_PROCS][36];
      int core_vertex_mask[NUM_PROCS][36];
      int has_edge_data;
      yac_int edge_global_ids[NUM_PROCS][60];
      int core_edge_mask[NUM_PROCS][60];
    } decomp[] =
      {{.nbr_vertices = {{2,5},{5,5},{6,2}},
        .lon_vertices = {{0,1},{1,2,3,4,5},{0,1,2,3,4,5}},
        .lat_vertices = {{0,1,2,3,4},{0,1,2,3,4},{4,5}},
        .cell_global_ids = {{0,5,10,15},
                            {1,2,3,4, 6,7,8,9, 11,12,13,14, 16,17,18,19},
                            {20,21,22,23,24}},
        .core_cell_mask = {{1,0,1,0},
                           {0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0},
                           {1,0,1,0,1}},
        .vertex_global_ids = {{0,1,6,7,12,13,18,19,24,25},
                              {1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17,
                               19,20,21,22,23, 25,26,27,28,29},
                              {24,25,26,27,28,29, 30,31,32,33,34,35}},
        .core_vertex_mask = {{1,0,
                              0,1,
                              1,0,
                              0,1,
                              1,0},
                             {0,1,0,1,0,
                              1,0,1,0,1,
                              0,1,0,1,0,
                              1,0,1,0,1,
                              0,1,0,1,0},
                             {1,0,1,0,1,0,
                              0,1,0,1,0,1}},
        .edge_global_ids = {{0,1,3,
                             11,12,14,
                             22,23,25,
                             33,34,36,
                             44},
                            {2,3,4,5,6,7,8,9,10,
                             13,14,15,16,17,18,19,20,21,
                             24,25,26,27,28,29,30,31,32,
                             35,36,37,38,39,40,41,42,43,
                             46,48,50,52},
                            {44,45,46,47,48,49,50,51,52,53,54,
                             55,56,57,58,59}},
        .core_edge_mask = {{1,0,0,
                            1,0,0,
                            1,0,0,
                            1,0,0,
                            1},
                           {1,0,1,0,1,0,1,0,0,
                            1,0,1,0,1,0,1,0,0,
                            1,0,1,0,1,0,1,0,0,
                            1,0,1,0,1,0,1,0,0,
                            1,1,1,1},
                           {1,0,1,0,1,0,1,0,1,0,0,
                            1,1,1,1,1}},
        .has_vertex_data = 1,
        .has_edge_data = 1},
       {.nbr_vertices = {{6,6},{0,0},{0,0}},
        .lon_vertices = {{0,1,2,3,4,5},{-1},{-1}},
        .lat_vertices = {{0,1,2,3,4,5},{-1},{-1}},
        .cell_global_ids = {{0,1,2,3,4,
                             5,6,7,8,9,
                             10,11,12,13,14,
                             15,16,17,18,19,
                             20,21,22,23,24}, {-1}, {-1}},
        .core_cell_mask = {{1,0,1,0,1,
                            0,1,0,1,0,
                            1,0,1,0,1,
                            0,1,0,1,0,
                            1,0,1,0,1}, {-1}, {-1}},
        .core_vertex_mask = {{1,0,1,0,1,0,
                              0,1,0,1,0,1,
                              1,0,1,0,1,0,
                              0,1,0,1,0,1,
                              1,0,1,0,1,0,
                              0,1,0,1,0,1},
                             {-1},{-1}},
        .vertex_global_ids = {{0,1,2,3,4,5,
                               6,7,8,9,10,11,
                               12,13,14,15,16,17,
                               18,19,20,21,22,23,
                               24,25,26,27,28,29,
                               30,31,32,33,34,35},
                              {-1},{-1}},
        .has_vertex_data = 1,
        .has_edge_data = 0},
       {.nbr_vertices = {{3,3},{6,6},{0,0}},
        .lon_vertices = {{0,1,2},{0,1,2,3,4,5},{-1}},
        .lat_vertices = {{0,1,2},{0,1,2,3,4,5},{-1}},
        .cell_global_ids = {{0,1,
                             5,6},
                            {0,1,2,3,4,
                             5,6,7,8,9,
                             10,11,12,13,14,
                             15,16,17,18,19,
                             20,21,22,23,24}, {-1}},
        .core_cell_mask = {{1,0,
                            0,1},
                           {1,0,1,0,1,
                            0,1,0,1,0,
                            1,0,1,0,1,
                            0,1,0,1,0,
                            1,0,1,0,1}, {-1}},
        .has_vertex_data = 0,
        .has_edge_data = 0}};
    enum {NUM_DECOMP = sizeof(decomp) / sizeof(decomp[0])};
    char const * grid_filename = "test_basic_grid_parallel.nc";

    if (yac_file_exists(grid_filename)) delete_file(grid_filename);

    // for all grid decompositions
    for (size_t decomp_idx = 0; decomp_idx < NUM_DECOMP; ++decomp_idx) {
      for (int with_cell_centers = 0; with_cell_centers <= 1; ++with_cell_centers) {
        for (int with_core_cell_mask = 0; with_core_cell_mask <= 1; ++with_core_cell_mask) {
          for (int with_cell_global_ids = 0; with_cell_global_ids <= 1; ++with_cell_global_ids) {
            for (int with_core_vertex_mask = 0; with_core_vertex_mask <= 1; ++with_core_vertex_mask) {
              for (int with_vertex_global_ids = 0; with_vertex_global_ids <= 1; ++with_vertex_global_ids) {
                for (int with_core_edge_mask = 0; with_core_edge_mask <= 1; ++with_core_edge_mask) {
                  for (int with_edge_global_ids = 0; with_edge_global_ids <= 1; ++with_edge_global_ids) {

                    if ((!decomp[decomp_idx].has_vertex_data &&
                         with_core_vertex_mask) ||
                        (!decomp[decomp_idx].has_vertex_data &&
                         with_vertex_global_ids) ||
                        (!decomp[decomp_idx].has_edge_data &&
                         with_core_edge_mask) ||
                        (!decomp[decomp_idx].has_edge_data &&
                         with_edge_global_ids)) continue;

                    // set up distributed grid
                    int cyclic[2] = {0,0};

                    struct yac_basic_grid * grid;
                    char const * grid_name = "test_grid";
                    if (decomp[decomp_idx].nbr_vertices[comm_rank][0] *
                        decomp[decomp_idx].nbr_vertices[comm_rank][1] > 0) {
                      struct yac_basic_grid_data grid_data =
                        yac_generate_basic_grid_data_reg_2d_deg(
                          decomp[decomp_idx].nbr_vertices[comm_rank], cyclic,
                          decomp[decomp_idx].lon_vertices[comm_rank],
                          decomp[decomp_idx].lat_vertices[comm_rank]);
                      grid = yac_basic_grid_new(grid_name, grid_data);
                    } else {
                      grid = yac_basic_grid_empty_new(grid_name);
                    }

                    // compute cell centers
                    struct yac_basic_grid_data * basic_grid_data =
                      yac_basic_grid_get_data(grid);
                    yac_coordinate_pointer cell_center_coords =
                      malloc(basic_grid_data->num_cells * sizeof(*cell_center_coords));
                    for (size_t i = 0; i < basic_grid_data->num_cells; ++i) {
                      cell_center_coords[i][0] = 0.0;
                      cell_center_coords[i][1] = 0.0;
                      cell_center_coords[i][2] = 0.0;
                      for (int j = 0; j < basic_grid_data->num_vertices_per_cell[i]; ++j)
                        for (size_t k = 0; k < 3; ++k)
                          cell_center_coords[i][k] +=
                            basic_grid_data->vertex_coordinates[
                              basic_grid_data->cell_to_vertex[
                                basic_grid_data->cell_to_vertex_offsets[i] + j]][k];
                      normalise_vector(cell_center_coords[i]);
                    }
                    if (with_cell_centers && (basic_grid_data->num_cells > 0))
                      yac_basic_grid_add_coordinates(
                        grid, YAC_LOC_CELL, cell_center_coords, basic_grid_data->num_cells);
                    free(cell_center_coords);
                    if (with_cell_global_ids && (basic_grid_data->num_cells > 0))
                      basic_grid_data->cell_ids =
                        to_pointer(
                          decomp[decomp_idx].cell_global_ids[comm_rank],
                          basic_grid_data->num_cells * sizeof(decomp[0].cell_global_ids[0]));
                    if (with_core_cell_mask && (basic_grid_data->num_cells > 0))
                      basic_grid_data->core_cell_mask =
                        to_pointer(
                          decomp[decomp_idx].core_cell_mask[comm_rank],
                          basic_grid_data->num_cells * sizeof(decomp[0].core_cell_mask[0]));
                    if (with_vertex_global_ids && (basic_grid_data->num_vertices > 0))
                      basic_grid_data->vertex_ids =
                        to_pointer(
                          decomp[decomp_idx].vertex_global_ids[comm_rank],
                          basic_grid_data->num_vertices * sizeof(decomp[0].vertex_global_ids[0]));
                    if (with_core_vertex_mask && (basic_grid_data->num_vertices > 0))
                      basic_grid_data->core_vertex_mask =
                        to_pointer(
                          decomp[decomp_idx].core_vertex_mask[comm_rank],
                          basic_grid_data->num_vertices * sizeof(decomp[0].core_vertex_mask[0]));
                    if (with_edge_global_ids && (basic_grid_data->num_edges > 0))
                      basic_grid_data->edge_ids =
                        to_pointer(
                          decomp[decomp_idx].edge_global_ids[comm_rank],
                          basic_grid_data->num_edges * sizeof(decomp[0].edge_global_ids[0]));
                    if (with_core_edge_mask && (basic_grid_data->num_edges > 0))
                      basic_grid_data->core_edge_mask =
                        to_pointer(
                          decomp[decomp_idx].core_edge_mask[comm_rank],
                          basic_grid_data->num_edges * sizeof(decomp[0].core_edge_mask[0][0]));

                    // set up io configuration
                    clear_yac_io_env();
                    setenv("YAC_IO_RANK_LIST", "0,2", 1);
                    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "3", 1);

                    // write distributed grid to file
                    yac_basic_grid_to_file_parallel(grid, grid_filename, MPI_COMM_WORLD);

                    yac_basic_grid_delete(grid);

                    // check generated grid file

                    if (comm_rank == 0) {
                      enum {NUM_CELLS = 25, NUM_CRN = 4};
                      double ref_cla[NUM_CELLS][NUM_CRN] =
                        {{0,0,1,1},{0,0,1,1},{0,0,1,1},{0,0,1,1},{0,0,1,1},
                        {1,1,2,2},{1,1,2,2},{1,1,2,2},{1,1,2,2},{1,1,2,2},
                        {2,2,3,3},{2,2,3,3},{2,2,3,3},{2,2,3,3},{2,2,3,3},
                        {3,3,4,4},{3,3,4,4},{3,3,4,4},{3,3,4,4},{3,3,4,4},
                        {4,4,5,5},{4,4,5,5},{4,4,5,5},{4,4,5,5},{4,4,5,5}};
                      double ref_clo[NUM_CELLS][NUM_CRN] =
                        {{0,1,1,0},{1,2,2,1},{2,3,3,2},{3,4,4,3},{4,5,5,4},
                        {0,1,1,0},{1,2,2,1},{2,3,3,2},{3,4,4,3},{4,5,5,4},
                        {0,1,1,0},{1,2,2,1},{2,3,3,2},{3,4,4,3},{4,5,5,4},
                        {0,1,1,0},{1,2,2,1},{2,3,3,2},{3,4,4,3},{4,5,5,4},
                        {0,1,1,0},{1,2,2,1},{2,3,3,2},{3,4,4,3},{4,5,5,4}};
                      double ref_lat[NUM_CELLS] =
                        {0.5,0.5,0.5,0.5,0.5,
                        1.5,1.5,1.5,1.5,1.5,
                        2.5,2.5,2.5,2.5,2.5,
                        3.5,3.5,3.5,3.5,3.5,
                        4.5,4.5,4.5,4.5,4.5};
                      double ref_lon[NUM_CELLS] =
                        {0.5,1.5,2.5,3.5,4.5,
                        0.5,1.5,2.5,3.5,4.5,
                        0.5,1.5,2.5,3.5,4.5,
                        0.5,1.5,2.5,3.5,4.5,
                        0.5,1.5,2.5,3.5,4.5};
                      int ref_cell_global_ids[NUM_CELLS] =
                        {0,1,2,3,4,
                        5,6,7,8,9,
                        10,11,12,13,14,
                        15,16,17,18,19,
                        20,21,22,23,24};
                      int ref_core_cell_mask[NUM_CELLS] =
                        {1,0,1,0,1,
                        0,1,0,1,0,
                        1,0,1,0,1,
                        0,1,0,1,0,
                        1,0,1,0,1};
                      int ref_vertex_global_ids[NUM_CELLS][NUM_CRN] =
                        {{0,1,7,6},{1,2,8,7},{2,3,9,8},{3,4,10,9},{4,5,11,10},
                         {6,7,13,12},{7,8,14,13},{8,9,15,14},{9,10,16,15},{10,11,17,16},
                         {12,13,19,18},{13,14,20,19},{14,15,21,20},{15,16,22,21},{16,17,23,22},
                         {18,19,25,24},{19,20,26,25},{20,21,27,26},{21,22,28,27},{22,23,29,28},
                         {24,25,31,30},{25,26,32,31},{26,27,33,32},{27,28,34,33},{28,29,35,34}};
                      int ref_core_vertex_mask[NUM_CELLS][NUM_CRN] =
                        {{1,0,1,0},{0,1,0,1},{1,0,1,0},{0,1,0,1},{1,0,1,0},
                         {0,1,0,1},{1,0,1,0},{0,1,0,1},{1,0,1,0},{0,1,0,1},
                         {1,0,1,0},{0,1,0,1},{1,0,1,0},{0,1,0,1},{1,0,1,0},
                         {0,1,0,1},{1,0,1,0},{0,1,0,1},{1,0,1,0},{0,1,0,1},
                         {1,0,1,0},{0,1,0,1},{1,0,1,0},{0,1,0,1},{1,0,1,0}};
                      int ref_edge_global_ids[NUM_CELLS][NUM_CRN] =
                        {{0,3,11,1},{2,5,13,3},{4,7,15,5},{6,9,17,7},{8,10,19,9},
                         {11,14,22,12},{13,16,24,14},{15,18,26,16},{17,20,28,18},{19,21,30,20},
                         {22,25,33,23},{24,27,35,25},{26,29,37,27},{28,31,39,29},{30,32,41,31},
                         {33,36,44,34},{35,38,46,36},{37,40,48,38},{39,42,50,40},{41,43,52,42},
                         {44,47,55,45},{46,49,56,47},{48,51,57,49},{50,53,58,51},{52,54,59,53}};
                      int ref_core_edge_mask[NUM_CELLS][NUM_CRN] =
                        {{1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},
                         {1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},
                         {1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},
                         {1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},
                         {1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0}};

                      check_grid_file(
                        grid_filename, grid_name, NUM_CELLS, NUM_CRN,
                        &ref_cla[0][0], &ref_clo[0][0],
                        with_cell_centers?ref_lat:NULL,
                        with_cell_centers?ref_lon:NULL,
                        with_cell_global_ids?ref_cell_global_ids:NULL,
                        with_core_cell_mask?ref_core_cell_mask:NULL,
                        with_vertex_global_ids?&ref_vertex_global_ids[0][0]:NULL,
                        with_core_vertex_mask?&ref_core_vertex_mask[0][0]:NULL,
                        with_edge_global_ids?&ref_edge_global_ids[0][0]:NULL,
                        with_core_edge_mask?&ref_core_edge_mask[0][0]:NULL);
                    } // comm_rank == 0
                    MPI_Barrier(MPI_COMM_WORLD);
                  } // with edge global ids
                } // with core edge mask
              } // with vertex global ids
            } // with core vertex mask
            delete_file(grid_filename);
          } // with cell global ids
        } // with core cell mask
      } // with cell centers
    } // decomp
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void delete_file(char const * filename) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) unlink(filename);
  MPI_Barrier(MPI_COMM_WORLD);
}
