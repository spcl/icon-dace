// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef GRID2VTK_H
#define GRID2VTK_H

#include "basic_grid.h"
#include "grid_cell.h"

// YAC PUBLIC HEADER START

/** \example test_grid2vtk.c
 * This contains some examples on how to use the
 * \ref yac_write_basic_grid_data_to_file and
 * \ref yac_write_grid_cells_to_file routine.
 */


/**
 * writes a \ref yac_basic_grid_data "basic grid data" to file
 * @param[in] grid basic grid data
 * @param[in] name file name (".vtk" will be added to this)
 * @remark as a reference, the equator will be included in the file
 */
void yac_write_basic_grid_data_to_file(
  struct yac_basic_grid_data * grid, char const * name);

/**
 * writes a \ref yac_basic_grid "basic grid" to file
 * @param[in] grid basic grid
 * @param[in] name file name (".vtk" will be added to this)
 * @remark as a reference, the equator will be included in the file
 */
void yac_write_basic_grid_to_file(
  struct yac_basic_grid * grid, char const * name);

/**
 * writes a list of cells into a vtk file, which can be visualised by paraview
 * @param[in] cells               list of cells
 * @param[in] num_cells           number of entries in cells
 * @param[in] name                file name (".vtk" will be added to this)
 * @param[in] num_points_per_edge each cell edge will be approximated by N
 *                                straight lines in 3D space, where N is
 *                                (num_points_per_edge - 1)
 * @remark the edge type will be taken into account
 */
void yac_write_grid_cells_to_file(
  struct yac_grid_cell * cells, size_t num_cells, char * name,
  size_t num_points_per_edge);

// YAC PUBLIC HEADER STOP

#endif // GRID2VTK_H

