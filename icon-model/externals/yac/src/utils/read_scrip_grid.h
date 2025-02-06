// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid.h"

// YAC PUBLIC HEADER START

/** \example test_read_scrip.c
 * This contains examples for yac_read_scrip_basic_grid
 */

/**
 * reads in an grid in SCRIP format
 * @param[in] grid_filename    name of the SCRIP grid netcdf file
 * @param[in] mask_filename    name of the SCRIP mask netcdf file
 * @param[in] grid_name        name of the grid in the file
 * @param[in] valid_mask_value value that marks cells as valid
 * @param[in]  use_ll_edges    if possible represent all edges using
 *                             lon/lat circles
 * @returns %yac_basic_grid_data structure that contains the grid
 */
struct yac_basic_grid_data yac_read_scrip_basic_grid_data(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value, int use_ll_edges);

/**
 * reads in grid data from a SCRIP formated file
 * @param[in]  grid_filename        name of the SCRIP grid netcdf file
 * @param[in]  mask_filename        name of the SCRIP mask netcdf file
 * @param[in]  grid_name            name of the grid in the file
 * @param[in]  valid_mask_value     value that marks cells as valid
 * @param[in]  name                 name of the grid
 * @param[in]  use_ll_edges         if possible represent all edges using
 *                                  lon/lat circles
 * @param[out] cell_coord_idx       index at which cell centers are registerd
 *                                  in the basic grid
 * @param[out] duplicated_cell_idx  indices of all duplicated cells
 * @param[out] orig_cell_global_ids global ids of the original cells
 * @param[out] nbr_duplicated_cells number of duplicated cells
 * @return basic grid
 * @remark This routine will allocate the arrays duplicated_cell_idx and
 *         orig_cell_global_ids. The user is responsible for freeing them.
 * @remark NULL is a valid argument for cell_coord_idx, duplicated_cell_idx,
 *         orig_cell_global_ids, and nbr_duplicated_cells
 */
struct yac_basic_grid * yac_read_scrip_basic_grid(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value, char const * name,
  int use_ll_edges, size_t * cell_coord_idx,
  size_t ** duplicated_cell_idx, yac_int ** orig_cell_global_ids,
  size_t * nbr_duplicated_cells);

/**
 * reads in grid data from a SCRIP formated file in parallel and applies a
 * IO decomposition to it
 * @param[in]  grid_filename        name of the SCRIP grid netcdf file
 * @param[in]  mask_filename        name of the SCRIP mask netcdf file
 * @param[in]  comm                 MPI communicator containing all proceses
 *                                  that will get a part of the grid
 * @param[in]  grid_name            name of the grid in the file
 * @param[in]  valid_mask_value     value that marks cells as valid
 * @param[in]  name                 name of the grid
 * @param[in]  use_ll_edges         if possible represent all edges using
 *                                  lon/lat circles
 * @param[out] cell_coord_idx       index at which cell centers are registerd
 *                                  in the basic grid
 * @param[out] duplicated_cell_idx  indices of all duplicated cells
 * @param[out] orig_cell_global_ids global ids of the original cells
 * @param[out] nbr_duplicated_cells number of duplicated cells
 * @return basic grid
 * @remark This routine will allocate the arrays duplicated_cell_idx and
 *         orig_cell_global_ids. The user is responsible for freeing them.
 * @remark NULL is a valid argument for cell_coord_idx, duplicated_cell_idx,
 *         orig_cell_global_ids, and nbr_duplicated_cells
 */
struct yac_basic_grid * yac_read_scrip_basic_grid_parallel(
  char const * grid_filename, char const * mask_filename,
  MPI_Comm comm, char const * grid_name, int valid_mask_value,
  char const * name, int use_ll_edges, size_t * cell_coord_idx,
  size_t ** duplicated_cell_idx, yac_int ** orig_cell_global_ids,
  size_t * nbr_duplicated_cells);

/**
 * reads in an grid in SCRIP format
 * @param[in]  grid_filename         name of the SCRIP grid netcdf file
 * @param[in]  mask_filename         name of the SCRIP mask netcdf file
 * @param[in]  grid_name             name of the grid in the file
 * @param[in]  valid_mask_value      value that marks cells as valid
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] x_vertices            longitude coordinates of the vertices
 * @param[out] y_vertices            latitude coordinates of the vertices
 * @param[out] x_cells               longitude coordinates of cell points
 * @param[out] y_cells               latitude coordinates of cell points
 * @param[out] cell_to_vertex        vertices indices per cell
 * @param[out] cell_core_mask        cell core mask
 * @param[out] duplicated_cell_idx  indices of all duplicated cells
 * @param[out] orig_cell_idx        indices of the original cells
 * @param[out] nbr_duplicated_cells number of duplicated cells
 * @remark NULL is a valid argument for duplicated_cell_idx,
 *         orig_cell_idx, and nbr_duplicated_cells
 */
void yac_read_scrip_grid_information(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value,
  size_t * num_vertices, size_t * num_cells, int ** num_vertices_per_cell,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells,
  int ** cell_to_vertex, int ** cell_core_mask, size_t ** duplicated_cell_idx,
  size_t ** orig_cell_idx, size_t * nbr_duplicated_cells);

// YAC PUBLIC HEADER STOP
