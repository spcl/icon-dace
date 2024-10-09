// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "tests.h"
#include "test_common.h"
#include "read_icon_grid.h"
#include "geometry.h"
#include "grid2vtk.h"
#include "io_utils.h"

#include <netcdf.h>

#include <mpi.h>

/* the grid file contains the following grid
 * (numbers == global cell/vertex indices, XX == masked cell)
 *
 *00-00-01-02-02-05-03-08-04   ------------------------
 * \ 00 /\ 11 /\ 12 /\ 13 /    \    /\ XX /\    /\    /
 * 01 03 04 06 07 09 10 11      \  /  \  /  \  /  \  /
 *   \/ 01 \/ 10 \/ 14 \/        \/    \/ XX \/    \/
 *   05-12-06-14-07-17-08         ------------------
 *    \ 02 /\ 09 /\ 15 /          \    /\ XX /\    /
 *    13 15 16 18 19 20            \  /  \  /  \  /
 *      \/ 03 \/ 08 \/              \/    \/ XX \/
 *      09-21-10-23-11               ------------
 *       \ 04 /\ 07 /                \    /\ XX /
 *       22 24 25 26                  \  /  \  /
 *         \/ 05 \/                    \/    \/
 *         12-27-13                     ------
 *          \ 06 /                      \    /
 *          28 29                        \  /
              \/                          \/
 *            14
 */

// global grid data
double vlon[15] = {-2,-1,0,1,2, -1.5,-0.5,0.5,1.5, -1,0,1, -0.5,0.5, 0};
double vlat[15] = {1.5,1.5,1.5,1.5,1.5, 0.5,0.5,0.5,0.5, -0.5,-0.5,-0.5,
                   -1.5,-1.5, -2.5};
double clon[16] = {1,1,0,0,-1,-1,-2,-1,0,0,1,1,1,1,1,0};
double clat[16] = {-1.5,-1,-1,-0.5,-0.5,0,0,0.5,0.5,0,0,-0.5,0.5,1,1,1.5};
int mask[16] = {0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0};

int vertex_of_cell[3][16] = {{1,2,6,7,10,11,13,11,8,7,3,2,3,4,4,8},
                             {2,6,7,10,11,13,14,12,11,8,7,3,4,5,8,9},
                             {6,7,10,11,13,14,15,14,12,11,8,7,8,9,9,12}};
int edge_of_cell[3][16] =   {{1,4,13,16,22,25,28,24,19,15,7,3,6,9,10,18},
                             {4,13,16,22,25,28,30,27,24,19,15,7,10,12,18,21},
                             {2,5,14,17,23,26,29,26,20,17,8,5,8,11,11,20}};
int cells_of_vertex[6][15] = {{1, 1,11,13,14,1, 2, 9,14,3, 4, 8,5,6,7},
                              {0, 2,12,14, 0,2, 3,10,15,4, 5, 9,6,7,0},
                              {0,12,13,15, 0,3, 4,11,16,5, 6,16,7,8,0},
                              {0, 0, 0, 0, 0,0,10,13, 0,0, 8, 0,0,0,0},
                              {0, 0, 0, 0, 0,0,11,15, 0,0, 9, 0,0,0,0},
                                {0, 0, 0, 0, 0,0,12,16, 0,0,10, 0,0,0,0}};
int vertex_of_edge[2][30] = {{1,1,2,2,2,3,3,3,4,4,4,5, 6,6,7,7,7,8,8,8,9, 10,10,11,11,11,12, 13,13,14},
                             {2,6,3,6,7,4,7,8,5,8,9,9, 7,10,8,10,11,9,11,12,12, 11,13,12,13,14,14, 14,15,15}};

static void write_test_grid_file(char const * file_name);

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_size, comm_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if ((comm_size != 4) && (comm_rank == 0)) {
    fputs("wrong number of processes (has to be four)\n", stderr);
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < sizeof(vlon) / sizeof(vlon[0]); ++i) vlon[i] *= YAC_RAD;
  for (size_t i = 0; i < sizeof(vlat) / sizeof(vlat[0]); ++i) vlat[i] *= YAC_RAD;
  for (size_t i = 0; i < sizeof(clon) / sizeof(clon[0]); ++i) clon[i] *= YAC_RAD;
  for (size_t i = 0; i < sizeof(clat) / sizeof(clat[0]); ++i) clat[i] *= YAC_RAD;

  char const filename[] = "test_read_icon_parallel_grid.nc";

  if (comm_rank == 0) write_test_grid_file(filename);

  // ensure that the grid file exists
  MPI_Barrier(MPI_COMM_WORLD);

  for (int n = -1; n < comm_size + 2; ++n) {

    int nbr_vertices, nbr_cells;
    int * num_vertices_per_cell;
    int * cell_to_vertex;
    int * global_cell_ids;
    int * cell_owner;
    int * global_vertex_ids;
    int * vertex_owner;
    double * x_vertices, * y_vertices;
    double * x_cells, * y_cells;
    int * cell_mask;

    yac_read_icon_grid_information_parallel(
      filename, MPI_COMM_WORLD, &nbr_vertices, &nbr_cells, &num_vertices_per_cell,
      &cell_to_vertex, &global_cell_ids, &cell_owner, &global_vertex_ids,
      &vertex_owner, &x_vertices, &y_vertices, &x_cells, &y_cells, &cell_mask);

    int ref_nbr_cells[4] = {11, 9, 14, 9};
    int ref_nbr_vertices[4] = {11, 10, 13, 10};
    int ref_cell_to_vertex[14][3];
    int ref_global_cell_ids[4][14] = {{0,1,2,3,4,5,7,8,9,10,11},
                                      {4,5,6,7,2,3,8,9,15},
                                      {8,9,10,11,0,1,2,3,4,5,7,12,14,15},
                                      {12,13,14,15,7,8,9,10,11}};
    int ref_cell_owner[14];
    int ref_global_vertex_ids[4][13] = {{0,1,2,5,6,7,9,10,11,12,13},
                                        {5,6,7,8,9,10,11,12,13,14},
                                        {0,1,2,3,5,6,7,8,9,10,11,12,13},
                                        {1,2,3,4,6,7,8,10,11,13}};
    int ref_vertex_owner[13];

    {
      int all_cell_to_vertex[16][3] =
        {{0,1,5}, {1,5,6}, {5,6,9}, {6,9,10}, {9,10,12}, {10,12,13}, {12,13,14},
         {10,11,13}, {7,10,11}, {6,7,10}, {2,6,7}, {1,2,6}, {2,3,7}, {3,4,8},
         {3,7,8}, {7,8,11}};
      for (int i = 0; i < ref_nbr_cells[comm_rank]; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < ref_nbr_vertices[comm_rank]; ++k) {
            if (all_cell_to_vertex[ref_global_cell_ids[comm_rank][i]][j] ==
                ref_global_vertex_ids[comm_rank][k]) {
              ref_cell_to_vertex[i][j] = k;
              break;
            }
          }
        }
      }
      int global_cell_owner[16] = {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3};
      for (int i = 0; i < ref_nbr_cells[comm_rank]; ++i) {
        ref_cell_owner[i] =
          global_cell_owner[ref_global_cell_ids[comm_rank][i]];
        if (ref_cell_owner[i] == comm_rank) ref_cell_owner[i] = -1;
      }
      int global_vertex_owner[15] = {0,2,3,3,3,0,2,3,3,1,2,3,1,1,1};
      for (int i = 0; i < ref_nbr_vertices[comm_rank]; ++i) {
        ref_vertex_owner[i] =
          global_vertex_owner[ref_global_vertex_ids[comm_rank][i]];
        if (ref_vertex_owner[i] == comm_rank) ref_vertex_owner[i] = -1;
      }
    }

    if (nbr_cells != ref_nbr_cells[comm_rank])
      PUT_ERR("wrong number of vertices\n");
    if (nbr_vertices != ref_nbr_vertices[comm_rank])
      PUT_ERR("wrong number of vertices\n");
    for (int i = 0; i < nbr_cells; ++i)
      if (num_vertices_per_cell[i] != 3)
        PUT_ERR("wrong number of vertices per cell\n");
    for (int i = 0; i < nbr_cells; ++i) {
      for (int j = 0; j < 3; ++j)
        if (cell_to_vertex[3 * i + j] != ref_cell_to_vertex[i][j])
          PUT_ERR("error in cell_to_vertex\n");
      if (global_cell_ids[i] != ref_global_cell_ids[comm_rank][i])
        PUT_ERR("wrong global cell id\n");
      if (cell_owner[i] != ref_cell_owner[i])
        PUT_ERR("wrong cell owner\n");
      if (double_are_unequal(x_cells[i],
                             clon[ref_global_cell_ids[comm_rank][i]]))
        PUT_ERR("wrong cell longitude\n");
      if (double_are_unequal(y_cells[i],
                             clat[ref_global_cell_ids[comm_rank][i]]))
        PUT_ERR("wrong cell latitude\n");
      if (cell_mask[i] != mask[ref_global_cell_ids[comm_rank][i]])
        PUT_ERR("wrong cell mask\n");
    }
    for (int i = 0; i < nbr_vertices; ++i) {
      if (global_vertex_ids[i] != ref_global_vertex_ids[comm_rank][i])
        PUT_ERR("wrong global vertex id\n");
      if (vertex_owner[i] != ref_vertex_owner[i]) PUT_ERR("wrong vertex owner\n");
      if (double_are_unequal(x_vertices[i],
                             vlon[ref_global_vertex_ids[comm_rank][i]]))
        PUT_ERR("wrong vertex longitude\n");
      if (double_are_unequal(y_vertices[i],
                            vlat[ref_global_vertex_ids[comm_rank][i]]))
        PUT_ERR("wrong vertex latitude\n");
    }

    free(num_vertices_per_cell);
    free(cell_to_vertex);
    free(cell_owner);
    free(vertex_owner);
    free(global_cell_ids);
    free(global_vertex_ids);
    free(x_vertices);
    free(y_vertices);
    free(x_cells);
    free(y_cells);
    free(cell_mask);
  }

  {
    set_even_io_rank_list(MPI_COMM_WORLD);
    struct yac_basic_grid_data grid_data =
      yac_read_icon_basic_grid_data_parallel(filename, MPI_COMM_WORLD);
    yac_basic_grid_data_free(grid_data);
    clear_yac_io_env();
  }

  // ensure that all processes finished reading the file
  MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank == 0) unlink(filename);

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void write_test_grid_file(char const * file_name) {

  int ncid;

  // create file
  yac_nc_create(file_name, NC_CLOBBER, &ncid);

  int dim_cell_id, dim_vertex_id, dim_edge_id, dim_nv_id, dim_ne_id, dim_nc_id;

  // define dimensions
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "cell", 16, &dim_cell_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "vertex", 15, &dim_vertex_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "edge", 30, &dim_edge_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "nv", 3, &dim_nv_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "ne", 6, &dim_ne_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "nc", 2, &dim_nc_id));


  int var_vlon_id, var_vlat_id, var_clon_id, var_clat_id, var_mask_id,
      var_v2c_id, var_c2v_id, var_c2e_id, var_v2e_id;

  // define variables
  YAC_HANDLE_ERROR(nc_def_var(ncid, "vlon", NC_DOUBLE, 1, &dim_vertex_id,
                          &var_vlon_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "vlat", NC_DOUBLE, 1, &dim_vertex_id,
                          &var_vlat_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "clon", NC_DOUBLE, 1, &dim_cell_id,
                          &var_clon_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "clat", NC_DOUBLE, 1, &dim_cell_id,
                          &var_clat_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "cell_sea_land_mask", NC_INT, 1, &dim_cell_id,
                          &var_mask_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "vertex_of_cell", NC_INT, 2, (int[]){dim_nv_id, dim_cell_id},
      &var_c2v_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "edge_of_cell", NC_INT, 2, (int[]){dim_nv_id, dim_cell_id},
      &var_c2e_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "cells_of_vertex", NC_INT, 2, (int[]){dim_ne_id, dim_vertex_id},
      &var_v2c_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "edge_vertices", NC_INT, 2, (int[]){dim_nc_id, dim_edge_id},
      &var_v2e_id));

  // end definition
  YAC_HANDLE_ERROR(nc_enddef(ncid));

  // write grid data
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_vlon_id, vlon));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_vlat_id, vlat));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_clon_id, clon));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_clat_id, clat));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_mask_id, mask));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_c2v_id, &(vertex_of_cell[0][0])));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_c2e_id, &(edge_of_cell[0][0])));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_v2c_id, &(cells_of_vertex[0][0])));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_v2e_id, &(vertex_of_edge[0][0])));

  YAC_HANDLE_ERROR(nc_close(ncid));
}
