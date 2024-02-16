/**
 * @file grid_reg2d.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
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
#include "grid.h"
#include "grid_reg2d_common.h"
#include "geometry.h"

static struct basic_grid_data yac_generate_basic_grid_data_reg_2d_(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices,
  void (*LLtoXYZ_ptr)(double, double, double[])) {

  YAC_ASSERT(
    !cyclic[1],
    "ERROR(yac_generate_basic_grid_data_reg_2d): "
    "cyclic[1] != 0 not yet supported")

  size_t num_cells_2d[2] =
    {nbr_vertices[0] - (cyclic[0]?0:1), nbr_vertices[1] - (cyclic[1]?0:1)};
  size_t num_vertices_2d[2] = {num_cells_2d[0] + (cyclic[0]?0:1), num_cells_2d[1] + 1};
  size_t num_vertices = num_vertices_2d[0] * num_vertices_2d[1];
  size_t num_edges =
    (num_cells_2d[0] + (cyclic[0]?0:1)) * num_cells_2d[1] +
    num_cells_2d[0] * (num_cells_2d[1] + 1);

  coordinate_pointer vertex_coordinates =
    xmalloc(num_vertices * sizeof(*vertex_coordinates));
  for (size_t i = 0, k = 0; i < num_vertices_2d[1]; ++i) {
    for (size_t j = 0; j < num_vertices_2d[0]; ++j, ++k) {
      LLtoXYZ_ptr(lon_vertices[j], lat_vertices[i],
                  &(vertex_coordinates[k][0]));
    }
  }

  enum yac_edge_type * edge_type = xmalloc(num_edges * sizeof(*edge_type));
  if (!cyclic[0]) {
    for (size_t i = 0, k = 0; i < num_cells_2d[1]; ++i) {
      for (size_t j = 0; j < num_vertices_2d[0]; ++j, k += 2) {
        edge_type[k] = LAT_CIRCLE_EDGE;
        edge_type[k+1] = LON_CIRCLE_EDGE;
      }
      --k;
      edge_type[k-1] = LON_CIRCLE_EDGE;
    }
  } else {
    for (size_t i = 1; i < num_edges; ++i)
      edge_type[i] = (i&1)?LAT_CIRCLE_EDGE:LON_CIRCLE_EDGE;
    for (size_t i = 0; i < num_cells_2d[1]; ++i) {
      edge_type[2*i*num_cells_2d[0]] = LAT_CIRCLE_EDGE;
      edge_type[2*(i+1)*num_cells_2d[0]-1] = LON_CIRCLE_EDGE;
    }
  }
  for (size_t i = num_edges - num_cells_2d[0]; i < num_edges; ++i)
    edge_type[i] = LAT_CIRCLE_EDGE;

  struct basic_grid_data grid =
    yac_generate_basic_grid_data_reg2d_common(nbr_vertices, cyclic);
  grid.vertex_coordinates = vertex_coordinates;
  grid.edge_type = edge_type;
  return grid;
}

struct basic_grid_data yac_generate_basic_grid_data_reg_2d(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_generate_basic_grid_data_reg_2d_(
      nbr_vertices, cyclic, lon_vertices, lat_vertices, LLtoXYZ);
}

struct basic_grid_data yac_generate_basic_grid_data_reg_2d_deg(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_generate_basic_grid_data_reg_2d_(
      nbr_vertices, cyclic, lon_vertices, lat_vertices, LLtoXYZ_deg);
}
