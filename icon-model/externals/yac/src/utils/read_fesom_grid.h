// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid.h"

// YAC PUBLIC HEADER START

/** \example test_read_fesom.c
 * This contains examples for yac_read_fesom_basic_grid_data.
 */

/**
 * reads in an fesom grid netcdf file and generates a struct
 * %yac_basic_grid_data from it
 * @param[in] filename name of the fesom grid netcdf file
 * @returns %yac_basic_grid_data structure that contains the fesom grid
 */
struct yac_basic_grid_data yac_read_fesom_basic_grid_data(
  char const * filename);

/**
 * reads in an fesom grid netcdf file and generates a struct
 * %yac_basic_grid from it
 * @param[in] filename name of the fesom grid netcdf file
 * @param[in] gridname name of the grid
 * @returns %yac_basic_grid structure that contains the fesom grid
 */
struct yac_basic_grid * yac_read_fesom_basic_grid(
  char const * filename, char const * gridname);

/**
 * reads in an fesom grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 * @param[in]  filename              name of the fesom grid netcdf file
 * @param[out] nbr_vertices          number of vertices in the grid
 * @param[out] nbr_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 *
 */
void yac_read_fesom_grid_information(const char * filename, int * nbr_vertices,
                                     int * nbr_cells, int ** num_vertices_per_cell,
                                     int ** cell_to_vertex, double ** x_vertices,
                                     double ** y_vertices, double ** x_cells,
                                     double ** y_cells);

// YAC PUBLIC HEADER STOP
