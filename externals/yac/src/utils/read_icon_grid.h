// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid.h"

// YAC PUBLIC HEADER START

#include <mpi.h>

/** \example test_read_icon.c
 * This contains examples for yac_read_icon_basic_grid_data.
 */

/** \example test_read_icon_parallel.c
 * This contains examples for parallel use of yac_read_icon_basic_grid_data.
 */

/**
 * reads in an icon grid netcdf file and generates a struct
 * %yac_basic_grid_data from it
 * @param[in] filename name of the icon grid netcdf file
 * @returns %basic_grid_data structure that contains the icon grid
 */
struct yac_basic_grid_data yac_read_icon_basic_grid_data(
  char const * filename);

/**
 * reads in an icon grid netcdf file and generates a struct
 * %basic_grid from it
 * @param[in] filename name of the icon grid netcdf file
 * @param[in] gridname name of the grid
 * @returns %yac_basic_grid structure that contains the icon grid
 */
struct yac_basic_grid * yac_read_icon_basic_grid(
  char const * filename, char const * gridname);

/**
 * reads in an icon grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 *
 * @param[in]  filename              name of the icon grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] cell_mask             mask for cells
 */
void yac_read_icon_grid_information(const char * filename, int * num_vertices,
                                    int * num_cells, int ** num_vertices_per_cell,
                                    int ** cell_to_vertex, double ** x_vertices,
                                    double ** y_vertices, double ** x_cells,
                                    double ** y_cells, int ** cell_mask);

/**
 * reads in an icon grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 *
 * @param[in]  filename              name of the icon grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] global_cell_id        global cell IDs
 * @param[out] cell_mask             mask for cells
 * @param[out] cell_core_mask        cell core mask
 * @param[out] global_corner_id      global corner IDs
 * @param[out] corner_core_mask      corner core mask
 * @param[out] rank                  local MPI rank
 * @param[out] size                  number of MPI ranks
 */
void yac_read_part_icon_grid_information(
  const char * filename, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size);

/**
 * reads in an icon grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface. The reading is done
 * in parallel and a basic domain decomposition is applied.
 *
 * @param[in]  filename              name of the icon grid netcdf file
 * @param[in]  comm                  MPI communicator containing all proceses
 *                                   that will get a part of the grid
 * @param[out] num_vertices          number of vertices in the local part of the
 *                                   grid
 * @param[out] num_cells             number of cells in the local part of the
 *                                   grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] global_cell_ids       global ids of local cells (core and halo
 *                                   cells)
 * @param[out] cell_owner            owner of each cell (locally owned cells
 *                                   are marked with -1)
 * @param[out] global_vertex_ids     global ids of local vertices (core and halo
 *                                   vertices)
 * @param[out] vertex_owner          owner of each vertex (locally owned
 *                                   vertices are marked with -1)
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] cell_mask             mask for cells
 */

void yac_read_icon_grid_information_parallel(
  const char * filename, MPI_Comm comm, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, int ** global_cell_ids,
  int ** cell_owner, int ** global_vertex_ids, int ** vertex_owner,
  double ** x_vertices, double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** cell_mask);

struct yac_basic_grid_data yac_read_icon_basic_grid_data_parallel(
  const char * filename, MPI_Comm comm);

struct yac_basic_grid * yac_read_icon_basic_grid_parallel(
  char const * filename, char const * gridname, MPI_Comm comm);

/**
 * destroys remaining icon grid data
 *
 * @param[out] cell_mask             mask for cells
 * @param[out] global_cell_id        global cell IDs
 * @param[out] cell_core_mask        cell core mask
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] global_corner_id      global corner IDs
 * @param[out] corner_core_mask      corner core mask
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 */

void yac_delete_icon_grid_data( int ** cell_mask,
                                int ** global_cell_id,
                                int ** cell_core_mask,
                                int ** num_vertices_per_cell,
                                int ** global_corner_id,
                                int ** corner_core_mask,
                                int ** cell_to_vertex,
                                double ** x_cells,
                                double ** y_cells,
                                double ** x_vertices,
                                double ** y_vertices);

// YAC PUBLIC HEADER STOP
