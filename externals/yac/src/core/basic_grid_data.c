// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <string.h>

#include "utils_core.h"
#include "area.h"
#include "basic_grid_data.h"

void yac_basic_grid_data_compute_cell_areas(
  struct yac_basic_grid_data grid, double * cell_areas) {

  // determine maximum number of vertices per cell
  int max_num_vertices_per_cell = 0;
  for (size_t i = 0; i < grid.num_cells; ++i)
    if (grid.num_vertices_per_cell[i] > max_num_vertices_per_cell)
      max_num_vertices_per_cell = grid.num_vertices_per_cell[i];

  struct yac_grid_cell cell =
    {.coordinates_xyz =
        xmalloc(
          (size_t)max_num_vertices_per_cell * sizeof(*(cell.coordinates_xyz))),
    .edge_type =
       xmalloc(
         (size_t)max_num_vertices_per_cell * sizeof(*(cell.edge_type))),
    .num_corners = (size_t)max_num_vertices_per_cell,
    .array_size = (size_t)max_num_vertices_per_cell};

  // for all cells in the grid
  for (size_t cell_idx = 0; cell_idx < grid.num_cells; ++cell_idx) {

    // generate grid cell
    size_t num_vertices = grid.num_vertices_per_cell[cell_idx];
    size_t * cell_to_vertex =
      grid.cell_to_vertex + grid.cell_to_vertex_offsets[cell_idx];
    size_t * cell_to_edge =
      grid.cell_to_edge + grid.cell_to_edge_offsets[cell_idx];
    for (size_t i = 0; i < num_vertices; ++i) {

      size_t vertex_idx = cell_to_vertex[i];
      memcpy(
        cell.coordinates_xyz[i], grid.vertex_coordinates[vertex_idx],
        3 * sizeof(cell.coordinates_xyz[i][0]));
      size_t edge_idx = cell_to_edge[i];
      cell.edge_type[i] = grid.edge_type[edge_idx];
    }
    cell.num_corners = num_vertices;

    cell_areas[cell_idx] = yac_huiliers_area(cell);
  }

  free(cell.edge_type);
  free(cell.coordinates_xyz);
}

void yac_basic_grid_data_free(struct yac_basic_grid_data grid) {

  free(grid.vertex_coordinates);
  free(grid.cell_ids);
  free(grid.vertex_ids);
  free(grid.edge_ids);
  free(grid.core_cell_mask);
  free(grid.core_vertex_mask);
  free(grid.core_edge_mask);
  free(grid.num_vertices_per_cell);
  free(grid.num_cells_per_vertex);
  free(grid.cell_to_vertex);
  free(grid.cell_to_vertex_offsets);
  free(grid.cell_to_edge);
  if (grid.cell_to_vertex_offsets != grid.cell_to_edge_offsets)
    free(grid.cell_to_edge_offsets);
  free(grid.vertex_to_cell);
  free(grid.vertex_to_cell_offsets);
  free(grid.edge_to_vertex);
  free(grid.edge_type);
}
