/**
 * @file grid2vtk.h
 *
 * @copyright Copyright  (C)  2018 DKRZ, MPI-M
 *
 * @version 1.0
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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


#include <stdio.h>
#include <string.h>
#include <float.h>

#include "grid2vtk.h"
#include "vtk_output.h"
#include "geometry.h"
#include "ensure_array_size.h"
#include "utils.h"

static void
get_bounding_circle_point_data(struct bounding_circle bnd_circle,
                               size_t num_points, double * points) {

  double n_x[3] = {1, 0, 0};
  double n_y[3] = {0, 1, 0};
  double angle_x = get_vector_angle(n_x, bnd_circle.base_vector);
  double angle_y = get_vector_angle(n_y, bnd_circle.base_vector);
  double * n = (angle_x > angle_y)?n_x:n_y;

  double temp_vector[3];
  double start_vector[3];

  crossproduct_d(bnd_circle.base_vector, n, temp_vector);
  rotate_vector2(
    temp_vector, bnd_circle.inc_angle, bnd_circle.base_vector, start_vector);

  double angle = 2.0 * M_PI / (num_points + 1);

  for (size_t i = 0; i < num_points; ++i)
    rotate_vector(
      bnd_circle.base_vector, angle * (double)i, start_vector, points + 3 * i);
}

static void
get_equator_point_data(size_t num_points, double * points) {

  // the equator is approximated by a polygon with num_points corners

  struct bounding_circle equator;
  equator.base_vector[0] = 0;
  equator.base_vector[1] = 0;
  equator.base_vector[2] = 1;
  equator.inc_angle = SIN_COS_M_PI_2;
  equator.sq_crd = DBL_MAX;

  get_bounding_circle_point_data(equator, num_points, points);
}

void yac_write_basic_grid_to_file(
  struct basic_grid_data * grid, char * name) {

  size_t num_grid_corners = grid->num_vertices;
  size_t num_grid_cells = grid->num_cells;

  size_t const num_points_per_great_circle = 71;

  //-------------------------------------
  // generate point data for the vtk file
  // it will contain data for all corners of the grid and the points for the equator
  //-------------------------------------
  coordinate_pointer points =
    xmalloc((num_grid_corners + num_points_per_great_circle) *sizeof(*points));

  // get the point data of the grid
  memcpy(points, grid->vertex_coordinates, num_grid_corners * sizeof(*points));
  // get the point data of the equator
  get_equator_point_data(
    num_points_per_great_circle, &(points[num_grid_corners][0]));

  //------------------------------------
  // generate cell data for the vtk file
  //------------------------------------

  // has entries for all cells of the grid and one for the equator
  unsigned * num_points_per_cell =
    xmalloc((grid->num_cells + 1) * sizeof(*num_points_per_cell));

  // get the number of points per grid cell
  for (size_t i = 0; i < grid->num_cells; ++i)
    num_points_per_cell[i] = (unsigned)(grid->num_vertices_per_cell[i]);

  // get the number of corners for the equator
  num_points_per_cell[grid->num_cells] = (unsigned)num_points_per_great_circle;

  // count the total number of cell points
  size_t total_num_cell_points = 0;
  for (size_t i = 0; i < grid->num_cells; ++i)
     total_num_cell_points += num_points_per_cell[i];

  unsigned * cell_data =
    xmalloc((total_num_cell_points + num_points_per_great_circle) *
            sizeof(*cell_data));

  // get the cell data of the grid
  for (size_t i = 0; i < total_num_cell_points; ++i)
    cell_data[i] = (unsigned)(grid->cell_to_vertex[i]);

  // write the cell data for the great circle
  for (size_t i = 0; i < num_points_per_great_circle; ++i)
    cell_data[i + total_num_cell_points] = (unsigned)(num_grid_corners + i);

  //------------------------------------
  // generate scalar data
  // we write a cell type that differentiates grid cells and the equator
  //------------------------------------
  unsigned * cell_type = xmalloc((num_grid_cells + 1) * sizeof(*cell_type));
  int * core_cell_mask = xmalloc((num_grid_cells + 1) * sizeof(*core_cell_mask));

  for (size_t i = 0; i < num_grid_cells; ++i) cell_type[i] = 0;
  cell_type[num_grid_cells] = 1;

  if (grid->core_cell_mask != NULL) {
    memcpy(core_cell_mask, grid->core_cell_mask,
           num_grid_cells * sizeof(*core_cell_mask));
  } else {
    for (size_t i = 0; i < num_grid_cells; ++i) core_cell_mask[i] = 1;
  }
  core_cell_mask[num_grid_cells] = 0;

  //------------------------------------
  // generate the file name
  //------------------------------------

  char * filename = xcalloc((strlen(name) + 5), sizeof(*filename));
  strcpy(filename, name);
  strcat(filename, ".vtk");

  //------------------------------------
  // generate the actual vtk file
  //------------------------------------
  VTK_FILE * file;

  file = vtk_open(filename, "grid data");
  vtk_write_point_data(file, &(points[0][0]),
                       num_grid_corners + num_points_per_great_circle);
  vtk_write_cell_data(file, cell_data, num_points_per_cell, num_grid_cells + 1);
  vtk_write_cell_scalars_uint(file, cell_type, num_grid_cells + 1, "cell_type");
  vtk_write_cell_scalars_int(file, core_cell_mask, num_grid_cells + 1, "core_cell_mask");

  vtk_close(file);

  //------------------------------------
  // some final cleanup
  //------------------------------------

  free(filename);
  free(core_cell_mask);
  free(cell_type);
  free(cell_data);
  free(num_points_per_cell);
  free(points);
}

static void get_edge_points(
  double * a, double * b, enum yac_edge_type edge_type, double (**points)[3],
  size_t * points_array_size, size_t * num_points, size_t num_points_per_edge) {

  double edge_length = get_vector_angle(a, b);

  if (edge_length <= yac_angle_tol) {

    ENSURE_ARRAY_SIZE(
      *points, *points_array_size, *num_points + 1);

    (*points)[*num_points][0] = a[0];
    (*points)[*num_points][1] = a[1];
    (*points)[*num_points][2] = a[2];

    *num_points += 1;

    return;
  }

  ENSURE_ARRAY_SIZE(
    *points, *points_array_size, *num_points + num_points_per_edge);

  if (edge_type == LAT_CIRCLE_EDGE) {

    double a_lon, b_lon, lat;
    XYZtoLL(a, &a_lon, &lat);
    XYZtoLL(b, &b_lon, &lat);

    double d_lon_angle = get_angle(b_lon, a_lon) / (double)num_points_per_edge;

    for (size_t i = 0, offset = *num_points; i < num_points_per_edge;
         ++i, ++offset)
      LLtoXYZ(a_lon + d_lon_angle * (double)i, lat, (*points)[offset]);

  } else {

    double rotation_axis[3];
    crossproduct_d(a, b, rotation_axis);
    normalise_vector(rotation_axis);

    double d_angle = edge_length / (double)num_points_per_edge;

    for (size_t i = 0, offset = *num_points; i < num_points_per_edge;
         ++i, ++offset)
      rotate_vector(rotation_axis, d_angle * (double)i, a, (*points)[offset]);
  }

  *num_points += num_points_per_edge;
}

void yac_write_grid_cells_to_file(
  struct grid_cell * cells, size_t num_cells, char * name,
  size_t num_points_per_edge) {

  if (num_cells == 0) return;

  num_points_per_edge = MAX(num_points_per_edge, 3);

  unsigned * num_points_per_cell =
    xmalloc(num_cells * sizeof(*num_points_per_cell));

  //-------------------------------------
  // get the points for all cells
  //-------------------------------------

  double (*points)[3] = NULL;
  size_t points_array_size = 0;
  size_t num_points = 0;

  for (size_t i = 0; i < num_cells; ++i) {

    struct grid_cell * curr_cell = cells + i;
    size_t num_edges = curr_cell->num_corners;

    size_t prev_num_points = num_points;

    // for all edges of the current cell
    for (size_t j = 0; j < num_edges; ++j)
      get_edge_points(curr_cell->coordinates_xyz[j],
                      curr_cell->coordinates_xyz[(j + 1) % num_edges],
                      curr_cell->edge_type[j], &points, &points_array_size,
                      &num_points, num_points_per_edge);

    // get the number of points for the current cell
    num_points_per_cell[i] = num_points - prev_num_points;
  }

  //------------------------------------
  // generate cell data for the vtk file
  //------------------------------------

  unsigned * cell_data = xmalloc(num_points * sizeof(*cell_data));

  // write the cell data for the great circle
  for (size_t i = 0; i < num_points; ++i) cell_data[i] = (unsigned)i;

  //------------------------------------
  // generate the file name
  //------------------------------------

  char * filename = xcalloc((strlen(name) + 5), sizeof(*filename));
  strcpy(filename, name);
  strcat(filename, ".vtk");

  //------------------------------------
  // generate the actual vtk file
  //------------------------------------

  VTK_FILE * file = vtk_open(filename, "grid data");
  vtk_write_point_data(file, &(points[0][0]), (unsigned)num_points);
  vtk_write_cell_data(file, cell_data, num_points_per_cell, (unsigned)num_cells);
  vtk_close(file);

  //------------------------------------
  // some final cleanup
  //------------------------------------

  free(filename);
  free(cell_data);
  free(num_points_per_cell);
  free(points);
}

