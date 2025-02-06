// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "dist_grid_utils.h"
#include "ppm/ppm_xfuncs.h"
#include "geometry.h"

struct temp_corner {
  yac_int id;
  size_t idx;
  int mask;
  double coord[3];
};

struct temp_edge {
  yac_int id;
  size_t idx;
  int mask;
  enum yac_edge_type type;
  yac_int vertex_idx[2];
};

struct temp_cell_data {
  size_t corner_idx, edge_idx;
};

struct temp_cell {
  yac_int id;
  int mask;
  size_t num_corners;
  struct temp_cell_data * data;
};

static struct yac_basic_grid_data generate_local_grid_data(
  size_t num_global_cells, struct temp_cell * global_cells,
  size_t num_global_vertices, struct temp_corner * global_vertices,
  size_t num_global_edges, struct temp_edge * global_edges,
  size_t num_local_cells, size_t * local_cell_indices, int with_halo) {

  for (size_t i = 0; i < num_global_cells; ++i) global_cells[i].mask = 2;
  for (size_t i = 0; i < num_global_vertices; ++i) global_vertices[i].mask = 2;
  for (size_t i = 0; i < num_global_edges; ++i) global_edges[i].mask = 2;

  for (size_t i = 0; i < num_local_cells; ++i) {
    struct temp_cell * curr_cell = global_cells + local_cell_indices[i];
    curr_cell->mask = 1;
    for (size_t j = 0; j < curr_cell->num_corners; ++j) {
      global_vertices[curr_cell->data[j].corner_idx].mask = 1;
      global_edges[curr_cell->data[j].edge_idx].mask = 1;
    }
  }

  if (with_halo) {
    for (size_t i = 0; i < num_global_cells; ++i) {
      if (global_cells[i].mask == 1) continue;
      int mask_flag = 0;
      for (size_t j = 0; j < global_cells[i].num_corners; ++j)
        if (global_vertices[global_cells[i].data[j].corner_idx].mask == 1)
          mask_flag = 1;
      if (mask_flag) {
        global_cells[i].mask = 0;
        for (size_t j = 0; j < global_cells[i].num_corners; ++j) {
          if (global_vertices[global_cells[i].data[j].corner_idx].mask == 2)
            global_vertices[global_cells[i].data[j].corner_idx].mask = 0;
          if (global_edges[global_cells[i].data[j].edge_idx].mask == 2)
            global_edges[global_cells[i].data[j].edge_idx].mask = 0;
        }
      }
    }
  }

  size_t num_cells = 0;
  size_t num_vertices = 0;
  size_t num_edges = 0;
  size_t num_vertices_per_cell_sum = 0;

  for (size_t i = 0; i < num_global_cells; ++i) {
    if (global_cells[i].mask != 2) {
      ++num_cells;
      num_vertices_per_cell_sum += global_cells[i].num_corners;
    }
  }
  for (size_t i = 0; i < num_global_vertices; ++i)
    if (global_vertices[i].mask != 2) ++num_vertices;
  for (size_t i = 0; i < num_global_edges; ++i)
    if (global_edges[i].mask != 2) ++num_edges;

  yac_coordinate_pointer vertex_coordinates =
    xmalloc(num_vertices * sizeof(*vertex_coordinates));
  yac_int * cell_ids = xmalloc(num_cells * sizeof(*cell_ids));
  yac_int * vertex_ids = xmalloc(num_vertices * sizeof(*vertex_ids));
  yac_int * edge_ids = xmalloc(num_edges * sizeof(*edge_ids));
  int * core_cell_mask = xmalloc(num_cells * sizeof(*core_cell_mask));
  int * core_vertex_mask = xmalloc(num_vertices * sizeof(*core_vertex_mask));
  int * core_edge_mask = xmalloc(num_edges * sizeof(*core_edge_mask));
  int * num_vertices_per_cell = xmalloc(num_cells * sizeof(*num_vertices_per_cell));
  int * num_cells_per_vertex = xcalloc(num_vertices, sizeof(*num_cells_per_vertex));
  size_t * cell_to_vertex = xmalloc(num_vertices_per_cell_sum * sizeof(*cell_to_vertex));
  size_t * cell_to_vertex_offsets = xmalloc(num_cells * sizeof(*cell_to_vertex_offsets));
  size_t * cell_to_edge = xmalloc(num_vertices_per_cell_sum * sizeof(*cell_to_edge));;
  size_t * cell_to_edge_offsets = cell_to_vertex_offsets;
  size_t * vertex_to_cell = xmalloc(num_vertices_per_cell_sum * sizeof(*vertex_to_cell));
  size_t * vertex_to_cell_offsets = xmalloc((num_vertices+1) * sizeof(*vertex_to_cell_offsets));
  yac_size_t_2_pointer edge_to_vertex = xmalloc(num_edges * sizeof(*edge_to_vertex));
  enum yac_edge_type * edge_type = xmalloc(num_edges * sizeof(*edge_type));

  for (size_t i = 0, j = 0; i < num_global_vertices; ++i) {
    if (global_vertices[i].mask != 2) {
      global_vertices[i].idx = j;
      for (int k = 0; k < 3; ++k)
        vertex_coordinates[j][k] = global_vertices[i].coord[k];
      vertex_ids[j] = global_vertices[i].id;
      core_vertex_mask[j] = global_vertices[i].mask;
      ++j;
    }
  }
  for (size_t i = 0, j = 0; i < num_global_edges; ++i) {
    if (global_edges[i].mask != 2) {
      global_edges[i].idx = j;
      edge_ids[j] = global_edges[i].id;
      core_edge_mask[j] = global_edges[i].mask;
      edge_type[j] = global_edges[i].type;
      edge_to_vertex[j][0] = global_vertices[global_edges[i].vertex_idx[0]].idx;
      edge_to_vertex[j][1] = global_vertices[global_edges[i].vertex_idx[1]].idx;
      ++j;
    }
  }
  for (size_t i = 0, j = 0, l = 0; i < num_global_cells; ++i) {
    if (global_cells[i].mask != 2) {
      cell_ids[j] = global_cells[i].id;
      core_cell_mask[j] = global_cells[i].mask;
      size_t num_corners = global_cells[i].num_corners;
      num_vertices_per_cell[j] = num_corners;
      for (size_t k = 0; k < num_corners; ++k, ++l) {
        cell_to_vertex[l] =
          global_vertices[global_cells[i].data[k].corner_idx].idx;
        cell_to_edge[l] =
          global_edges[global_cells[i].data[k].edge_idx].idx;
      }
      ++j;
    }
  }

  for (size_t i = 0, k = 0; i < num_cells; ++i)
    for (int j = 0; j < num_vertices_per_cell[i]; ++j, ++k)
      num_cells_per_vertex[cell_to_vertex[k]]++;

  for (size_t i = 0, accu = 0; i < num_cells; ++i) {
    cell_to_vertex_offsets[i] = accu;
    accu += num_vertices_per_cell[i];
  }

  vertex_to_cell_offsets[0] = 0;
  for (size_t i = 0, accu = 0; i < num_vertices; ++i) {
    vertex_to_cell_offsets[i+1] = accu;
    accu += num_cells_per_vertex[i];
  }
  for (size_t i = 0, k = 0; i < num_cells; ++i)
    for (int j = 0; j < num_vertices_per_cell[i]; ++j, ++k)
      vertex_to_cell[vertex_to_cell_offsets[cell_to_vertex[k]+1]++] = i;

  struct yac_basic_grid_data grid_data;

  grid_data.vertex_coordinates = vertex_coordinates;
  grid_data.cell_ids = cell_ids;
  grid_data.vertex_ids = vertex_ids;
  grid_data.edge_ids = edge_ids;
  grid_data.num_cells = num_cells;
  grid_data.num_vertices = num_vertices;
  grid_data.num_edges = num_edges;
  grid_data.core_cell_mask = core_cell_mask;
  grid_data.core_vertex_mask = core_vertex_mask;
  grid_data.core_edge_mask = core_edge_mask;
  grid_data.num_vertices_per_cell = num_vertices_per_cell;
  grid_data.num_cells_per_vertex = num_cells_per_vertex;
  grid_data.cell_to_vertex = cell_to_vertex;
  grid_data.cell_to_vertex_offsets = cell_to_vertex_offsets;
  grid_data.cell_to_edge = cell_to_edge;
  grid_data.cell_to_edge_offsets = cell_to_edge_offsets;
  grid_data.vertex_to_cell = vertex_to_cell;
  grid_data.vertex_to_cell_offsets = vertex_to_cell_offsets;
  grid_data.edge_to_vertex = edge_to_vertex;
  grid_data.edge_type = edge_type;
  grid_data.num_total_cells = num_cells;
  grid_data.num_total_vertices = num_vertices;
  grid_data.num_total_edges = num_edges;

  return grid_data;
}

struct yac_basic_grid_data yac_generate_basic_grid_data_reg2d(
  double const * global_coords_x, double const * global_coords_y,
  size_t const num_global_cells_[2], size_t const local_start[2],
  size_t const local_count[2], int with_halo) {

  size_t num_global_cells = num_global_cells_[0] * num_global_cells_[1];
  size_t num_global_vertices =
    (num_global_cells_[0] + 1) * (num_global_cells_[1] + 1);
  size_t num_global_edges =
    (num_global_cells_[0] + 1) * num_global_cells_[1] +
    (num_global_cells_[1] + 1) * num_global_cells_[0];

  struct temp_cell_data * cell_data =
    xmalloc(4 * num_global_cells * sizeof(*cell_data));
  struct temp_cell * global_cells =
    xmalloc(num_global_cells * sizeof(*global_cells));
  struct temp_corner * global_corners =
    xmalloc(num_global_vertices * sizeof(*global_corners));
  struct temp_edge * global_edges =
    xmalloc(num_global_edges * sizeof(*global_edges));

  for (size_t i = 0; i < num_global_cells; ++i) {

    struct temp_cell * curr_cell = global_cells + i;
    curr_cell->id = i;
    curr_cell->num_corners = 4;
    curr_cell->data = cell_data + 4 * i;

    size_t y_index = i / num_global_cells_[0];
    size_t x_index = i - y_index * num_global_cells_[0];

    curr_cell->data[0].corner_idx =  y_index      * (num_global_cells_[0] + 1) + x_index;
    curr_cell->data[1].corner_idx =  y_index      * (num_global_cells_[0] + 1) + x_index + 1;
    curr_cell->data[2].corner_idx = (y_index + 1) * (num_global_cells_[0] + 1) + x_index + 1;
    curr_cell->data[3].corner_idx = (y_index + 1) * (num_global_cells_[0] + 1) + x_index;

    curr_cell->data[0].edge_idx =  (2 * num_global_cells_[0] + 1) * y_index + 2 * x_index;
    curr_cell->data[1].edge_idx = curr_cell->data[0].edge_idx + 3 - (x_index+1 == num_global_cells_[0]);
    if (y_index+1 == num_global_cells_[1]) 
      curr_cell->data[2].edge_idx = (2 * num_global_cells_[0] + 1) * (y_index + 1) + x_index;
    else
      curr_cell->data[2].edge_idx = curr_cell->data[0].edge_idx + 2 * num_global_cells_[0] + 1;
    curr_cell->data[3].edge_idx = curr_cell->data[0].edge_idx + 1;
  }
  for (size_t i = 0; i < num_global_vertices; ++i) {
    global_corners[i].id = i;
    size_t y_index = i / (num_global_cells_[0] + 1);
    LLtoXYZ(global_coords_x[i - y_index * (num_global_cells_[0] + 1)],
            global_coords_y[y_index], &(global_corners[i].coord[0]));
  }
  for (size_t i = 0; i < num_global_edges; ++i) {

    global_edges[i].id = i;

    size_t y_index = i / (2*num_global_cells_[0]+1);
    size_t x_index = i - y_index * (2*num_global_cells_[0]+1);
    if (y_index == num_global_cells_[1])  {
      global_edges[i].type = YAC_LAT_CIRCLE_EDGE;
      x_index *= 2;
    } else if (x_index == 2*num_global_cells_[0]) {
      global_edges[i].type = YAC_LON_CIRCLE_EDGE;
      ++x_index;
    } else global_edges[i].type = (x_index&1)?YAC_LON_CIRCLE_EDGE:YAC_LAT_CIRCLE_EDGE;
    global_edges[i].vertex_idx[0] =
      x_index / 2 + y_index * (num_global_cells_[0] + 1);
    global_edges[i].vertex_idx[1] =
    global_edges[i].vertex_idx[0] + 1 + ((x_index & 1)?num_global_cells_[0]:0);
  }

  size_t num_local_cells = local_count[0] * local_count[1];
  size_t * local_cell_indices =
    xmalloc(num_local_cells * sizeof(*local_cell_indices));
  for (size_t i = 0, k = 0; i < local_count[1]; ++i)
    for (size_t j = 0; j < local_count[0]; ++j, ++k)
      local_cell_indices[k] =
        j + local_start[0] + (i + local_start[1]) * num_global_cells_[0];

  struct yac_basic_grid_data grid_data =
    generate_local_grid_data(
      num_global_cells, global_cells,
      num_global_vertices, global_corners,
      num_global_edges, global_edges,
      num_local_cells, local_cell_indices, with_halo);

  free(local_cell_indices);
  free(global_edges);
  free(global_corners);
  free(global_cells);
  free(cell_data);

  return grid_data;
}

struct yac_basic_grid * yac_generate_basic_grid_reg2d(
  char const * name, double const * coordinates_x,
  double const * coordinates_y, size_t const num_cells[2],
  size_t const local_start[2], size_t const local_count[2], int with_halo) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start, local_count, with_halo));
}
