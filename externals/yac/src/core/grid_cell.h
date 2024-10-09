// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>

#ifndef GRID_CELL_H
#define GRID_CELL_H

// YAC PUBLIC HEADER START

enum yac_edge_type {
   YAC_GREAT_CIRCLE_EDGE = 0, //!< great circle
   YAC_LAT_CIRCLE_EDGE   = 1, //!< latitude circle
   YAC_LON_CIRCLE_EDGE   = 2, //!< longitude circle
};

struct yac_grid_cell {
   double (*coordinates_xyz)[3];
   enum yac_edge_type * edge_type;
   size_t num_corners;
   size_t array_size; //!< size in elements of the arrays: coordinates_x,
                      //!< coordinates_y, edge_type and 1/3 of coordinates_xyz
};

enum yac_cell_type {
  YAC_LON_LAT_CELL,
  YAC_LAT_CELL,
  YAC_GREAT_CIRCLE_CELL,
  YAC_MIXED_CELL
};

/**
 * initiates a grid_cell object
 * before the first being used a grid_cell object has to be initialised
 * @param[in] cell object to be initialised
 * @see free_grid_cell
 * @see get_grid_cell
 */
void yac_init_grid_cell(struct yac_grid_cell * cell);

/**
 * copies a given grid cell
 * @param[in]  in_cell  cell to be copied
 * @param[out] out_cell copied cell
 * @remarks out_cell needs to be a cell that has previously been
 *          initialised or a cell that already contains valid data
 */
void yac_copy_grid_cell(struct yac_grid_cell in_cell, struct yac_grid_cell * out_cell);

/**
 * frees all memory associated with a grid_cell object and reinitialised
 * the cell
 * @param[in,out] cell
 */
void yac_free_grid_cell(struct yac_grid_cell * cell);

// YAC PUBLIC HEADER STOP

#ifdef YAC_DEBUG_GRID_CELL
/**
 * prints out info about a grid_cell object and reinitialised
 * the cell
 * @param[in] stream
 * @param[in] cell
 * @param[in] name
 */
void yac_print_grid_cell(FILE * stream, struct yac_grid_cell cell, char * name);
#endif

#endif // GRID_CELL_H

