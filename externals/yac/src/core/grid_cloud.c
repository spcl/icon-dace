// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "basic_grid_data.h"
#include "geometry.h"
#include "utils_common.h"

static struct yac_basic_grid_data yac_generate_basic_grid_data_cloud_(
  size_t nbr_points, double *x_points, double *y_points,
  void (*LLtoXYZ_ptr)(double, double, double[])) {

  yac_coordinate_pointer vertex_coordinates =
    xmalloc(nbr_points * sizeof(*vertex_coordinates));
  for (size_t i = 0; i < nbr_points; ++i)
    LLtoXYZ_ptr(
      x_points[i], y_points[i], &(vertex_coordinates[i][0]));

  struct yac_basic_grid_data grid;
  grid.vertex_coordinates = vertex_coordinates;
  grid.cell_ids = NULL;
  grid.vertex_ids = NULL;
  grid.edge_ids = NULL;
  grid.num_cells = 0;
  grid.num_vertices = nbr_points;
  grid.num_edges = 0;
  grid.core_cell_mask = NULL;
  grid.core_vertex_mask = NULL;
  grid.core_edge_mask = NULL;
  grid.num_vertices_per_cell = NULL;
  grid.num_cells_per_vertex =
    xcalloc(nbr_points, sizeof(*grid.num_cells_per_vertex));
  grid.cell_to_vertex = NULL;
  grid.cell_to_vertex_offsets = NULL;
  grid.cell_to_edge = NULL;
  grid.cell_to_edge_offsets = NULL;
  grid.vertex_to_cell = NULL;
  grid.vertex_to_cell_offsets =
    xcalloc(nbr_points, sizeof(*grid.vertex_to_cell_offsets));
  grid.edge_to_vertex = NULL;
  grid.edge_type = NULL;
  grid.num_total_cells = 0;
  grid.num_total_vertices = nbr_points;
  grid.num_total_edges = 0;
  return grid;
}

struct yac_basic_grid_data yac_generate_basic_grid_data_cloud(
  size_t nbr_points, double *x_points, double *y_points) {

  return
    yac_generate_basic_grid_data_cloud_(
      nbr_points, x_points, y_points, LLtoXYZ);
}

struct yac_basic_grid_data yac_generate_basic_grid_data_cloud_deg(
  size_t nbr_points, double *x_points, double *y_points) {

  return
    yac_generate_basic_grid_data_cloud_(
      nbr_points, x_points, y_points, LLtoXYZ_deg);
}
