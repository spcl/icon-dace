// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid.h"

// YAC PUBLIC HEADER START

/** \example test_read_cube_csv.c
 * This contains examples for yac_read_cube_csv_grid.
 */

/**
 * reads in an cube sphere grid csv file and generates a struct %grid from it
 * @param[in] filename name of the cube sphere grid csv file
 * @returns %grid structure that contains the cube sphere grid
 */
struct yac_basic_grid_data yac_read_cube_csv_grid(char * filename);

/**
 * reads in an cube sphere grid csv file and return the grid information in
 * a format that is supported by the YAC user interface.
 * @param[in]  filename              name of the cube sphere grid csv file
 * @param[out] nbr_vertices          number of vertices in the grid
 * @param[out] nbr_cells             number of cells in the grid
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 */
void yac_read_cube_csv_grid_information(const char * filename, int * nbr_vertices,
                                        int * nbr_cells, int ** cell_to_vertex,
                                        double ** x_vertices, double ** y_vertices);

// YAC PUBLIC HEADER STOP
