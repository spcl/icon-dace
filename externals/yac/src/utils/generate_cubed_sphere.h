// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef GENERATE_CUBED_SPHERE_H
#define GENERATE_CUBED_SPHERE_H

#include "basic_grid.h"

// YAC PUBLIC HEADER START

/** \example test_generate_cubed_sphere.c
 * A test for the generation of cubed sphere grids used in several toy models.
 */


void yac_generate_cubed_sphere_grid_information(
  unsigned n, unsigned * num_cells, unsigned * num_vertices,
  double ** x_vertices, double ** y_vertices, double ** z_vertices,
  unsigned ** vertices_of_cell, unsigned ** face_id);

void yac_generate_part_cube_grid_information(
  unsigned n, unsigned * nbr_vertices, unsigned * nbr_cells,
  unsigned ** num_vertices_per_cell, unsigned ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size);

struct yac_basic_grid_data yac_generate_cubed_sphere_grid(unsigned n);

struct yac_basic_grid * yac_generate_cubed_sphere_basic_grid(
  char const * name, size_t n);

void yac_write_cubed_sphere_grid(unsigned n, char const * filename);

// YAC PUBLIC HEADER STOP

#endif // GENERATE_CUBED_SPHERE_H

