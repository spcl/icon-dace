// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid_data.h"
#include "grid_reg2d_common.h"
#include "geometry.h"
#include "utils_common.h"

static struct yac_basic_grid_data yac_generate_basic_grid_data_curve_2d_(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices,
  void (*LLtoXYZ_ptr)(double, double, double[])) {

  YAC_ASSERT(
    !cyclic[1],
    "ERROR(yac_generate_basic_grid_data_curve_2d): "
    "cyclic[1] != 0 not yet supported")

  size_t num_cells_2d[2] =
    {nbr_vertices[0] - (cyclic[0]?0:1), nbr_vertices[1] - (cyclic[1]?0:1)};
  size_t num_vertices_2d[2] = {num_cells_2d[0] + 1, num_cells_2d[1] + 1};
  size_t num_vertices = num_vertices_2d[0] * num_vertices_2d[1];
  size_t num_edges =
    (num_cells_2d[0] + 1) * num_cells_2d[1] +
    num_cells_2d[0] * (num_cells_2d[1] + 1);

  yac_coordinate_pointer vertex_coordinates =
    xmalloc(num_vertices * sizeof(*vertex_coordinates));
  for (size_t i = 0; i < num_vertices; ++i)
    LLtoXYZ_ptr(lon_vertices[i], lat_vertices[i],
                &(vertex_coordinates[i][0]));

  enum yac_edge_type * edge_type = xmalloc(num_edges * sizeof(*edge_type));
  for (size_t i = 0; i < num_edges; ++i) edge_type[i] = YAC_GREAT_CIRCLE_EDGE;

  struct yac_basic_grid_data grid =
    yac_generate_basic_grid_data_reg2d_common(nbr_vertices, cyclic);
  grid.vertex_coordinates = vertex_coordinates;
  grid.edge_type = edge_type;
  return grid;
}

struct yac_basic_grid_data yac_generate_basic_grid_data_curve_2d(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_generate_basic_grid_data_curve_2d_(
      nbr_vertices, cyclic, lon_vertices, lat_vertices, LLtoXYZ);
}

struct yac_basic_grid_data yac_generate_basic_grid_data_curve_2d_deg(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_generate_basic_grid_data_curve_2d_(
      nbr_vertices, cyclic, lon_vertices, lat_vertices, LLtoXYZ_deg);
}
