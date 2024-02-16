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
#include <stdlib.h>
#include <stdio.h>

#include "grid.h"
#include "geometry.h"
#include "read_cube_csv_grid.h"
#include "utils.h"

static FILE * open_cube_csv_file    ( const char *filename);

static void close_cube_csv_file   ( FILE * file );

#define _GET_NTH_ARG(_1, _2, _3, _4, _5, N, ...) N
#define EXPAND(x) x
#define COUNT_ARGS(...) EXPAND(_GET_NTH_ARG(__VA_ARGS__, 5, 4, 3, 2, 1))
#define READ_LINE(format, ...) \
  YAC_ASSERT_F( \
    fscanf(file, format, __VA_ARGS__) == COUNT_ARGS(__VA_ARGS__), \
    "ERROR(read_cube_csv_grid_information): " \
    "failed while reading a line from file %s", filename)

void read_cube_csv_grid_information(const char * filename, int * nbr_vertices,
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

struct basic_grid_data read_cube_csv_grid(char * filename) {

  int nbr_vertices;
  int nbr_cells;
  int * cell_to_vertex;

  double * x_vertices;
  double * y_vertices;

  read_cube_csv_grid_information(filename, &nbr_vertices, &nbr_cells,
                                 &cell_to_vertex, &x_vertices, &y_vertices);

  int * num_vertices_per_cell =
    xmalloc(nbr_cells * sizeof(*num_vertices_per_cell));

  for (int i = 0; i < nbr_cells; ++i)
    num_vertices_per_cell[i] = 4;

  struct basic_grid_data grid =
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

