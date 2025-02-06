// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid.h"

// YAC PUBLIC HEADER START

/** \example test_read_mpiom.c
 * This contains examples for yac_read_mpiom_basic_grid_data.
 */

/**
 * reads in an mpiom grid netcdf file and generates a struct
 * %yac_basic_grid_data from it
 * @param[in] filename name of the mpiom grid netcdf file
 * @returns %yac_basic_grid_data structure that contains the mpiom grid
 */
struct yac_basic_grid_data yac_read_mpiom_basic_grid_data(
  char const * filename);

/**
 * reads in an mpiom grid netcdf file and generates a struct
 * %yac_basic_grid_data from it
 * @param[in] filename name of the mpiom grid netcdf file
 * @param[in] gridname name of the grid
 * @returns %yac_basic_grid structure that contains the mpiom grid
 */
struct yac_basic_grid * yac_read_mpiom_basic_grid(
  char const * filename, char const * gridname);

/**
 * reads in an mpiom grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 * @param[in]  filename              name of the mpiom grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] cellmask              mask value for cells
 *
 */
void yac_read_mpiom_grid_information(const char * filename, int * num_vertices,
                                     int * num_cells, int ** num_vertices_per_cell,
                                     int ** cell_to_vertex, double ** x_vertices,
                                     double ** y_vertices, double ** x_cells,
                                     double ** y_cells, int ** cellmask);

/**
 * reads in an partition mpiom grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 * @param[in]  filename              name of the mpiom grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] global_cell_id        global cell IDs
 * @param[out] cell_mask             mask value for cells
 * @param[out] cell_core_mask        cell core mask
 * @param[out] global_corner_id      global corner IDs
 * @param[out] corner_core_mask      corner core mask
 * @param[out] rank                  local MPI rank
 * @param[out] size                  number of MPI ranks
 *
 */
void yac_read_part_mpiom_grid_information(
  const char * filename, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, double ** x_vertices,
  double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size);

// YAC PUBLIC HEADER STOP
