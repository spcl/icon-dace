// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>

#include "read_cube_csv_grid.h"
#include "geometry.h"
#include "utils_common.h"

static FILE * open_cube_csv_file    ( const char *filename);

static void close_cube_csv_file   ( FILE * file );

#define _GET_NTH_ARG(_1, _2, _3, _4, _5, N, ...) N
#define EXPAND(x) x
#define COUNT_ARGS(...) EXPAND(_GET_NTH_ARG(__VA_ARGS__, 5, 4, 3, 2, 1))
#define READ_LINE(format, ...) \
  YAC_ASSERT_F( \
    fscanf(file, format, __VA_ARGS__) == COUNT_ARGS(__VA_ARGS__), \
    "ERROR(yac_read_cube_csv_grid_information): " \
    "failed while reading a line from file %s", filename)

void yac_read_cube_csv_grid_information(const char * filename, int * nbr_vertices,
                                        int * nbr_cells, int ** cell_to_vertex,
                                        double ** x_vertices, double ** y_vertices) {

  // open file

  FILE * file = open_cube_csv_file ( filename );

  // read number of vertices and cells
  READ_LINE("%d,%d\n", nbr_vertices, nbr_cells);

  // read coordinates of vertices

  *x_vertices = xmalloc (*nbr_vertices * sizeof(**x_vertices));
  *y_vertices = xmalloc (*nbr_vertices * sizeof(**y_vertices));

  for (int i = 0, dummy; i < *nbr_vertices; ++i) {
    READ_LINE("%d,%lf,%lf\n", &dummy, *x_vertices+i, *y_vertices+i);
    (*x_vertices)[i] *= YAC_RAD;
    (*y_vertices)[i] *= YAC_RAD;
  }

  // read indices of cell vertices

  *cell_to_vertex = xmalloc(4 * *nbr_cells * sizeof(**cell_to_vertex));

  for (int i = 0, dummy; i < *nbr_cells; ++i) {
    READ_LINE(
      "%d,%d,%d,%d,%d\n", &dummy, *cell_to_vertex+4*i+0,
      *cell_to_vertex+4*i+1, *cell_to_vertex+4*i+2,
      *cell_to_vertex+4*i+3);

    for (unsigned j = 0; j < 4; ++j)
      (*cell_to_vertex+4*i)[j]--;
  }

  close_cube_csv_file ( file );
}

struct yac_basic_grid_data yac_read_cube_csv_grid(char * filename) {

  int nbr_vertices;
  int nbr_cells;
  int * cell_to_vertex;

  double * x_vertices;
  double * y_vertices;

  yac_read_cube_csv_grid_information(filename, &nbr_vertices, &nbr_cells,
                                     &cell_to_vertex, &x_vertices, &y_vertices);

  int * num_vertices_per_cell =
    xmalloc(nbr_cells * sizeof(*num_vertices_per_cell));

  for (int i = 0; i < nbr_cells; ++i)
    num_vertices_per_cell[i] = 4;

  struct yac_basic_grid_data grid =
    yac_generate_basic_grid_data_unstruct(
      (size_t)nbr_vertices, (size_t)nbr_cells, num_vertices_per_cell,
      x_vertices, y_vertices, cell_to_vertex);
  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(cell_to_vertex);

  return grid;
}

static FILE * open_cube_csv_file (const char *filename) {

  FILE * file = xfopen(filename, "r");

  YAC_ASSERT_F(
    file, "ERROR(open_cube_csv_file): could not open file %s", filename);

  return file;
}

/* ---------------------------------------------------------------- */

static void close_cube_csv_file (FILE * file) {

   xfclose(file);
}

