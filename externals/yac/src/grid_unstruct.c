/**
 * @file grid_unstruct.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <string.h>

#include "grid.h"
#include "geometry.h"

struct temp_edge {
  yac_int vertex[2];
  size_t cell_to_edge_idx;
};

static int compare_temp_edges(void const * a, void const * b) {

  struct temp_edge const * edge_a = (struct temp_edge const *)a;
  struct temp_edge const * edge_b = (struct temp_edge const *)b;

  if (edge_a->vertex[0] != edge_b->vertex[0])
    return (edge_a->vertex[0] > edge_b->vertex[0])?1:-1;
  return (edge_a->vertex[1] > edge_b->vertex[1]) -
         (edge_a->vertex[1] < edge_b->vertex[1]);
}

static struct basic_grid_data yac_generate_basic_grid_data_unstruct_(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell_,
  double *x_vertices, double *y_vertices, int *cell_to_vertex_,
  void (*LLtoXYZ_ptr)(double, double, double[])) {

  coordinate_pointer vertex_coordinates =
    xmalloc(nbr_vertices * sizeof(*vertex_coordinates));
  for (size_t i = 0; i < nbr_vertices; ++i)
    LLtoXYZ_ptr(
      x_vertices[i], y_vertices[i], &(vertex_coordinates[i][0]));

  int * num_vertices_per_cell =
    xmalloc(nbr_cells * sizeof(*num_vertices_per_cell));
  memcpy(num_vertices_per_cell, num_vertices_per_cell_,
         nbr_cells * sizeof(*num_vertices_per_cell));

  size_t total_num_cell_corners = 0;
  size_t * cell_to_vertex_offsets =
    xmalloc(nbr_cells * sizeof(*cell_to_vertex_offsets));
  for (size_t i = 0; i < nbr_cells; ++i) {
    cell_to_vertex_offsets[i] = total_num_cell_corners;
    total_num_cell_corners += (size_t)(num_vertices_per_cell[i]);
  }

  for (size_t i = 0; i < total_num_cell_corners; ++i){
    YAC_ASSERT_F(cell_to_vertex_[i] >= 0,
      "ERROR(yac_generate_basic_grid_data_unstruct_): "
      "Index %zu in cell_to_vertex (%d) is negative.",
      i, cell_to_vertex_[i])
      YAC_ASSERT_F(cell_to_vertex_[i] < nbr_vertices,
      "ERROR(yac_generate_basic_grid_data_unstruct_): "
      "Index %zu in cell_to_vertex (%d) is not smaller than nbr_vertices (%zu).",
        i, cell_to_vertex_[i], nbr_vertices)
  }

  size_t * cell_to_vertex =
    xmalloc(total_num_cell_corners * sizeof(*cell_to_vertex));
  for (size_t i = 0; i < total_num_cell_corners; ++i)
    cell_to_vertex[i] = (size_t)(cell_to_vertex_[i]);

  int * num_cells_per_vertex =
    xcalloc(nbr_vertices, sizeof(*num_cells_per_vertex));
  for (size_t i = 0; i < total_num_cell_corners; ++i)
    num_cells_per_vertex[cell_to_vertex[i]]++;

  size_t * vertex_to_cell =
    xmalloc(total_num_cell_corners * sizeof(*vertex_to_cell));
  size_t * vertex_to_cell_offsets =
    xmalloc((nbr_vertices + 1) * sizeof(*vertex_to_cell_offsets));
  vertex_to_cell_offsets[0] = 0;
  for (size_t i = 0, accu = 0; i < nbr_vertices; ++i) {
    vertex_to_cell_offsets[i+1] = accu;
    accu += num_cells_per_vertex[i];
  }
  for (size_t i = 0, k = 0; i < nbr_cells; ++i) {
    size_t curr_num_vertices = num_vertices_per_cell[i];
    for (size_t j = 0; j < curr_num_vertices; ++j, ++k)
      vertex_to_cell[vertex_to_cell_offsets[cell_to_vertex[k]+1]++] = i;
  }

  struct temp_edge * temp_edges =
    xmalloc(total_num_cell_corners * sizeof(*temp_edges));
  for (size_t i = 0, offset = 0, k = 0; i < nbr_cells; ++i) {
    size_t * curr_cell_to_vertex = cell_to_vertex + offset;
    size_t curr_num_edges = num_vertices_per_cell[i];
    offset += curr_num_edges;
    for (size_t j = 0; j < curr_num_edges; ++j, ++k) {
      int order =
        curr_cell_to_vertex[j] > curr_cell_to_vertex[(j+1)%curr_num_edges];
      temp_edges[k].vertex[order] = curr_cell_to_vertex[j];
      temp_edges[k].vertex[order^1] = curr_cell_to_vertex[(j+1)%curr_num_edges];
      temp_edges[k].cell_to_edge_idx = k;
    }
  }
  qsort(temp_edges, total_num_cell_corners,
        sizeof(*temp_edges), compare_temp_edges);

  size_t * cell_to_edge =
    xmalloc(total_num_cell_corners * sizeof(*cell_to_edge));
  size_t nbr_edges = 0;
  size_t_2_pointer edge_to_vertex = (size_t_2_pointer)temp_edges;
  for (size_t i = 0, prev_indices[2] = {SIZE_MAX, SIZE_MAX};
       i < total_num_cell_corners; ++i) {

    size_t curr_cell_to_edge_idx = temp_edges[i].cell_to_edge_idx;
    if ((prev_indices[0] != temp_edges[i].vertex[0]) ||
        (prev_indices[1] != temp_edges[i].vertex[1])) {

      prev_indices[0] = temp_edges[i].vertex[0];
      prev_indices[1] = temp_edges[i].vertex[1];
      edge_to_vertex[nbr_edges][0] = prev_indices[0];
      edge_to_vertex[nbr_edges][1] = prev_indices[1];
      ++nbr_edges;
    }

    cell_to_edge[curr_cell_to_edge_idx] = nbr_edges - 1;
  }
  edge_to_vertex =
    xrealloc(edge_to_vertex, nbr_edges * sizeof(*edge_to_vertex));

  enum yac_edge_type * edge_type = xmalloc(nbr_edges * sizeof(*edge_type));
  for (size_t i = 0; i < nbr_edges; ++i) edge_type[i] = GREAT_CIRCLE_EDGE;

  struct basic_grid_data grid;
  grid.vertex_coordinates = vertex_coordinates;
  grid.cell_ids = NULL;
  grid.vertex_ids = NULL;
  grid.edge_ids = NULL;
  grid.num_cells = nbr_cells;
  grid.num_vertices = nbr_vertices;
  grid.num_edges = nbr_edges;
  grid.core_cell_mask = NULL;
  grid.core_vertex_mask = NULL;
  grid.core_edge_mask = NULL;
  grid.num_vertices_per_cell = num_vertices_per_cell;
  grid.num_cells_per_vertex = num_cells_per_vertex;
  grid.cell_to_vertex = cell_to_vertex;
  grid.cell_to_vertex_offsets = cell_to_vertex_offsets;
  grid.cell_to_edge = cell_to_edge;
  grid.cell_to_edge_offsets = cell_to_vertex_offsets;
  grid.vertex_to_cell = vertex_to_cell;
  grid.vertex_to_cell_offsets = vertex_to_cell_offsets;
  grid.edge_to_vertex = edge_to_vertex;
  grid.edge_type = edge_type;
  grid.num_total_cells = nbr_cells;
  grid.num_total_vertices = nbr_vertices;
  grid.num_total_edges = nbr_edges;
  return grid;
}

struct basic_grid_data yac_generate_basic_grid_data_unstruct(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell_,
  double *x_vertices, double *y_vertices, int *cell_to_vertex_) {

  return
    yac_generate_basic_grid_data_unstruct_(
      nbr_vertices, nbr_cells, num_vertices_per_cell_,
      x_vertices, y_vertices, cell_to_vertex_, LLtoXYZ);
}

struct basic_grid_data yac_generate_basic_grid_data_unstruct_deg(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell_,
  double *x_vertices, double *y_vertices, int *cell_to_vertex_) {

  return
    yac_generate_basic_grid_data_unstruct_(
      nbr_vertices, nbr_cells, num_vertices_per_cell_,
      x_vertices, y_vertices, cell_to_vertex_, LLtoXYZ_deg);
}

struct basic_grid_data yac_generate_basic_grid_data_unstruct_ll_(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell_,
  double *x_vertices, double *y_vertices, int *cell_to_vertex_,
  void (*LLtoXYZ_ptr)(double, double, double[])) {

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_unstruct_(
      nbr_vertices, nbr_cells, num_vertices_per_cell_,
      x_vertices, y_vertices, cell_to_vertex_, LLtoXYZ_ptr);

  for (size_t i = 0; i < grid_data.num_edges; ++i) {
    int is_lon_edge =
      (fabs(x_vertices[grid_data.edge_to_vertex[i][0]] -
            x_vertices[grid_data.edge_to_vertex[i][1]]) < yac_angle_tol) ||
      ((fabs(fabs(y_vertices[grid_data.edge_to_vertex[i][0]]) -
             M_PI_2) < yac_angle_tol) ^
       (fabs(fabs(y_vertices[grid_data.edge_to_vertex[i][1]]) -
             M_PI_2) < yac_angle_tol));
    int is_lat_edge =
      fabs(y_vertices[grid_data.edge_to_vertex[i][0]] -
           y_vertices[grid_data.edge_to_vertex[i][1]]) < yac_angle_tol;
    YAC_ASSERT_F(
      is_lon_edge ^ is_lat_edge,
      "ERROR(yac_generate_basic_grid_data_unstruct_ll): "
      "edge is neither lon nor lat ((%lf,%lf),(%lf,%lf))",
      x_vertices[grid_data.edge_to_vertex[i][0]],
      y_vertices[grid_data.edge_to_vertex[i][0]],
      x_vertices[grid_data.edge_to_vertex[i][1]],
      y_vertices[grid_data.edge_to_vertex[i][1]]);
    grid_data.edge_type[i] = (is_lon_edge)?LON_CIRCLE_EDGE:LAT_CIRCLE_EDGE;
  }

  return grid_data;
}

struct basic_grid_data yac_generate_basic_grid_data_unstruct_ll(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell_,
  double *x_vertices, double *y_vertices, int *cell_to_vertex_) {

  return
    yac_generate_basic_grid_data_unstruct_ll_(
      nbr_vertices, nbr_cells, num_vertices_per_cell_,
      x_vertices, y_vertices, cell_to_vertex_, LLtoXYZ);
}

struct basic_grid_data yac_generate_basic_grid_data_unstruct_ll_deg(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell_,
  double *x_vertices, double *y_vertices, int *cell_to_vertex_) {

  return
    yac_generate_basic_grid_data_unstruct_ll_(
      nbr_vertices, nbr_cells, num_vertices_per_cell_,
      x_vertices, y_vertices, cell_to_vertex_, LLtoXYZ_deg);
}
