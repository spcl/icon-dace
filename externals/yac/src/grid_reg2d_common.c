/**
 * @file grid_reg2d_common.c
 *
 * @copyright Copyright  (C)  2017 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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

#include "grid_reg2d_common.h"

struct basic_grid_data yac_generate_basic_grid_data_reg2d_common(
  size_t nbr_vertices[2], int cyclic[2]) {

  YAC_ASSERT(
    !cyclic[1],
    "ERROR(yac_generate_basic_grid_data_reg2d_common): "
    "grids that are cyclic in the y direction are not yet supported")

  size_t num_cells_2d[2] =
    {nbr_vertices[0] - ((cyclic[0])?0:1), nbr_vertices[1] - 1};
  size_t num_vertices_2d[2] = {num_cells_2d[0] + ((cyclic[0])?0:1),
                               num_cells_2d[1] + 1};
  size_t num_cells = num_cells_2d[0] * num_cells_2d[1];
  size_t num_vertices = num_vertices_2d[0] * num_vertices_2d[1];
  size_t num_edges =
    (num_cells_2d[1] + 1) *                 num_cells_2d[0] +
    (num_cells_2d[0] + ((cyclic[0])?0:1)) * num_cells_2d[1];

  int * num_vertices_per_cell =
    xmalloc(num_cells * sizeof(*num_vertices_per_cell));
  for (size_t i = 0; i < num_cells; ++i) num_vertices_per_cell[i] = 4;

  int * num_cells_per_vertex =
    xmalloc(num_vertices * sizeof(*num_cells_per_vertex));
  for (size_t i = 0; i < num_vertices; ++i) num_cells_per_vertex[i] = 4;
  if (!cyclic[0]) {
    for (size_t i = 0, j = num_vertices - num_vertices_2d[0];
         i < num_vertices_2d[0]; ++i, ++j) {
      num_cells_per_vertex[i] = 2;
      num_cells_per_vertex[j] = 2;
    }
    for (size_t i = 0, j = num_vertices_2d[0] - 1;
         i < num_vertices; i+=num_vertices_2d[0], j+=num_vertices_2d[0]) {
      num_cells_per_vertex[i] = 2;
      num_cells_per_vertex[j] = 2;
    }
    num_cells_per_vertex[0] = 1;
    num_cells_per_vertex[num_vertices_2d[0] - 1] = 1;
    num_cells_per_vertex[num_vertices - num_vertices_2d[0]] = 1;
    num_cells_per_vertex[num_vertices - 1] = 1;
  } else {
    for (size_t i = 0, j = num_vertices - num_vertices_2d[0];
         i < num_vertices_2d[0]; ++i, ++j) {
      num_cells_per_vertex[i] = 2;
      num_cells_per_vertex[j] = 2;
    }
  }

  size_t * cell_to_vertex = xmalloc(4 * num_cells * sizeof(*cell_to_vertex));
  size_t * cell_to_vertex_offsets =
    xmalloc(num_cells * sizeof(*cell_to_vertex_offsets));
  if (!cyclic[0]) {
    for (size_t i = 0, k = 0; i < num_cells_2d[1]; ++i) {
      size_t temp_cell_to_vertex[4] =
        {i * num_vertices_2d[0],
        i * num_vertices_2d[0] + 1,
        (i + 1) * num_vertices_2d[0] + 1,
        (i + 1) * num_vertices_2d[0]};
      for (size_t j = 0; j < num_cells_2d[0]; ++j)
        for (size_t l = 0; l < 4; ++l, ++k)
          cell_to_vertex[k] = temp_cell_to_vertex[l] + j;
    }
  } else {
    size_t temp_cell_to_vertex[4] =
      {0, 1, num_vertices_2d[0] + 1, num_vertices_2d[0]};
    for (size_t i = 0, k = 0, m = 0; i < num_cells_2d[1]; ++i)
      for (size_t j = 0; j < num_cells_2d[0]; ++j, ++m)
        for (int l = 0; l < 4; ++l, ++k)
          cell_to_vertex[k] = temp_cell_to_vertex[l] + m;
    for (size_t i = 0; i < num_cells_2d[1]; ++i) {
      cell_to_vertex[4*((i+1)*num_cells_2d[0]-1)+1] = i * num_cells_2d[0];
      cell_to_vertex[4*((i+1)*num_cells_2d[0]-1)+2] = (i + 1) * num_cells_2d[0];
    }
  }

  for (size_t i = 0; i < num_cells; ++i) cell_to_vertex_offsets[i] = 4 * i;

  size_t * cell_to_edge = xmalloc(4 * num_cells * sizeof(*cell_to_edge));
  if (!cyclic[0]) {
    for (size_t i = 0, k = 0; i < num_cells_2d[1]; ++i) {
      size_t edge_id_offset = (2 * num_cells_2d[0] + 1);
      size_t edge_id = i * edge_id_offset;
      for (size_t j = 0; j < num_cells_2d[0]; ++j, k += 4, edge_id += 2) {

        cell_to_edge[k+0] = edge_id;
        cell_to_edge[k+1] = edge_id + 3;
        cell_to_edge[k+2] = edge_id + edge_id_offset;
        cell_to_edge[k+3] = edge_id + 1;
      }
    }
    for (size_t i = num_edges - num_cells_2d[0],
                j = 4 * num_cells_2d[0] * (num_cells_2d[1] - 1) + 2;
                i < num_edges; ++i, j += 4) cell_to_edge[j] = i;
    for (size_t i = 0, j = 4 * (num_cells_2d[0] - 1) + 1;
                i < num_cells_2d[1]; ++i, j += 4 * num_cells_2d[0])
      cell_to_edge[j]--;
  } else {
    size_t edge_id_offset = 2 * num_cells_2d[0];
    for (size_t i = 0, k = 0; i < num_cells_2d[1]; ++i) {
      size_t edge_id = i * edge_id_offset + 1;
      for (size_t j = 0; j < num_cells_2d[0]; ++j, k += 4, edge_id += 2) {

        cell_to_edge[k+0] = edge_id;
        cell_to_edge[k+1] = edge_id + 3;
        cell_to_edge[k+2] = edge_id + edge_id_offset;
        cell_to_edge[k+3] = edge_id + 1;
      }
    }
    for (size_t i = 0; i < num_cells_2d[1]; ++i) {
      cell_to_edge[4*i*num_cells_2d[0]+0]--;
      cell_to_edge[4*i*num_cells_2d[0]+2]--;
      cell_to_edge[4*((i+1)*num_cells_2d[0]-2)+1]--;
    }
    for (size_t i = 0, * curr_cell_to_edge; i < num_cells_2d[1]; ++i) {
      curr_cell_to_edge = &cell_to_edge[4*((i+1)*num_cells_2d[0]-1)];
      size_t edge_id = 1 + i * edge_id_offset;
      curr_cell_to_edge[0] = edge_id;
      curr_cell_to_edge[1] = edge_id + 1;
      curr_cell_to_edge[2] = edge_id + edge_id_offset;
      curr_cell_to_edge[3] = edge_id + edge_id_offset - 2;
    }
    for (size_t i = 0,
         * curr_edge = cell_to_edge + 4 * (num_cells - num_cells_2d[0] + 1) + 2,
         edge_id = 2 * (num_cells_2d[0] * num_cells_2d[1] + 1);
         i < num_cells_2d[0] - 2; ++i, curr_edge += 4, ++edge_id) {
      *curr_edge = edge_id;
    }
  }

  size_t vertex_to_cell_size;
  size_t * vertex_to_cell;
  if (!cyclic[0]) {
    vertex_to_cell_size =
      4 * 1 + // 4 corners with 1 cell per vertex
      2 * 2 * (num_cells_2d[0] - 1) + // lower and upper edge with 2 cells per vertex
      2 * 2 * (num_cells_2d[1] - 1) + // 2 side edges with 2 cells per vertex
      4 * (num_cells_2d[0] - 1) * (num_cells_2d[1] - 1); // core vertices
    vertex_to_cell =
      xmalloc(vertex_to_cell_size * sizeof(*vertex_to_cell));
    size_t offset = 0;
    vertex_to_cell[offset++] = 0;
    for (size_t i = 0; i < num_cells_2d[0] - 1; ++i) {
      vertex_to_cell[offset++] = i;
      vertex_to_cell[offset++] = i + 1;
    }
    vertex_to_cell[offset++] = num_cells_2d[0] - 1;
    for (size_t i = 0, cell_idx = 0; i < num_cells_2d[1] - 1; ++i) {

      vertex_to_cell[offset++] = cell_idx;
      vertex_to_cell[offset++] = cell_idx + num_cells_2d[0];
      for (size_t j = 0; j < num_cells_2d[0] - 1; ++j, ++cell_idx) {
        vertex_to_cell[offset++] = cell_idx;
        vertex_to_cell[offset++] = cell_idx + 1;
        vertex_to_cell[offset++] = cell_idx + num_cells_2d[0];
        vertex_to_cell[offset++] = cell_idx + num_cells_2d[0] + 1;
      }
      vertex_to_cell[offset++] = cell_idx;
      vertex_to_cell[offset++] = cell_idx + num_cells_2d[0];
      ++cell_idx;
    }
    vertex_to_cell[offset++] = num_cells - num_cells_2d[0];
    for (size_t i = num_cells - num_cells_2d[0]; i < num_cells - 1; ++i) {
      vertex_to_cell[offset++] = i;
      vertex_to_cell[offset++] = i + 1;
    }
    vertex_to_cell[offset++] = num_cells - 1;
  } else {
    vertex_to_cell_size = 4 * (num_vertices - num_vertices_2d[0]);
    vertex_to_cell =
      xmalloc(vertex_to_cell_size * sizeof(*vertex_to_cell));

    size_t offset = 0;
    for (size_t i = 0; i < num_cells_2d[0]; ++i) {
      vertex_to_cell[offset++] = i - 1;
      vertex_to_cell[offset++] = i;
    }
    vertex_to_cell[0] = num_cells_2d[0] - 1;

    size_t vertex_idx = 0;
    for (size_t i = 0; i < num_cells_2d[1] - 1; ++i) {
      for (size_t j = 0; j < num_cells_2d[0]; ++j, ++vertex_idx) {
        vertex_to_cell[offset++] = vertex_idx - 1;
        vertex_to_cell[offset++] = vertex_idx;
        vertex_to_cell[offset++] = vertex_idx + num_cells_2d[0];
        vertex_to_cell[offset++] = vertex_idx + num_cells_2d[0] - 1;
      }
      vertex_to_cell[offset - 4 * num_cells_2d[0]] = (vertex_idx - 1);
      vertex_to_cell[offset - 4 * num_cells_2d[0] + 3] =
        vertex_idx - 1 + num_cells_2d[0];
    }

    for (size_t i = 0; i < num_cells_2d[0]; ++i, ++vertex_idx) {
      vertex_to_cell[offset++] = vertex_idx - 1;
      vertex_to_cell[offset++] = vertex_idx;
    }
    vertex_to_cell[offset - 2 * num_cells_2d[0]] = (vertex_idx - 1);
  }

  size_t * vertex_to_cell_offsets =
    xmalloc(num_vertices * sizeof(*vertex_to_cell_offsets));
  for (size_t i = 0, accu = 0; i < num_vertices; ++i) {
    vertex_to_cell_offsets[i] = accu;
    accu += num_cells_per_vertex[i];
  }

  size_t_2_pointer edge_to_vertex =
    xmalloc(num_edges * sizeof(*edge_to_vertex));
  if (!cyclic[0]) {
    for (size_t i = 0, k = 0, vertex_idx = 0; i < num_cells_2d[1]; ++i) {
      for (size_t j = 0; j < num_cells_2d[0]; ++j, k += 2, ++vertex_idx) {
        edge_to_vertex[k][0] = vertex_idx;
        edge_to_vertex[k][1] = vertex_idx + 1;
        edge_to_vertex[k+1][0] = vertex_idx;
        edge_to_vertex[k+1][1] = vertex_idx + num_vertices_2d[0];
      }
      edge_to_vertex[k][0] = vertex_idx + num_vertices_2d[0];
      edge_to_vertex[k][1] = vertex_idx;
      ++vertex_idx;
      ++k;
    }
    for (size_t i = num_edges - num_cells_2d[0],
         vertex_idx = num_vertices - num_vertices_2d[0]; i < num_edges;
         ++i, ++vertex_idx) {
      edge_to_vertex[i][0] = vertex_idx;
      edge_to_vertex[i][1] = vertex_idx + 1;
    }
  } else {
    for (size_t i = 1, vertex_idx = 0;
         i < 2 * num_cells_2d[0] * num_cells_2d[1]; i += 2, ++vertex_idx) {
      edge_to_vertex[i][0] = vertex_idx;
      edge_to_vertex[i][1] = vertex_idx + 1;
      edge_to_vertex[i+1][0] = vertex_idx;
      edge_to_vertex[i+1][1] = vertex_idx + num_vertices_2d[0];
    }
    for (size_t i = 0; i <= num_cells_2d[1]; ++i) {
      edge_to_vertex[2*i*num_cells_2d[0]][0] = i * num_cells_2d[0];
      edge_to_vertex[2*i*num_cells_2d[0]][1] = i * num_cells_2d[0] + 1;
      edge_to_vertex[2*i*num_cells_2d[0]+1][0] = i * num_cells_2d[0];
      edge_to_vertex[2*i*num_cells_2d[0]+1][1] = (i + 1) * num_cells_2d[0]-1;
    }
    for (size_t i = 1; i <= num_cells_2d[1]; ++i) {
      edge_to_vertex[2*i*num_cells_2d[0]-1][0] = i * num_cells_2d[0] - 1;
      edge_to_vertex[2*i*num_cells_2d[0]-1][1] = (i + 1) * num_cells_2d[0] - 1;
    }
    for (size_t i = 2; i < num_cells_2d[0]; ++i) {
      edge_to_vertex[2*num_cells_2d[0]*num_cells_2d[1]+i][0] =
        num_cells_2d[0]*num_cells_2d[1]+i-1;
      edge_to_vertex[2*num_cells_2d[0]*num_cells_2d[1]+i][1] =
        num_cells_2d[0]*num_cells_2d[1]+i;
    }
  }

  struct basic_grid_data grid;
  grid.vertex_coordinates = NULL;
  grid.cell_ids = NULL;
  grid.vertex_ids = NULL;
  grid.edge_ids = NULL;
  grid.num_cells = num_cells;
  grid.num_vertices = num_vertices;
  grid.num_edges = num_edges;
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
  grid.edge_type = NULL;
  grid.num_total_cells = num_cells;
  grid.num_total_vertices = num_vertices;
  grid.num_total_edges = num_edges;
  return grid;
}
