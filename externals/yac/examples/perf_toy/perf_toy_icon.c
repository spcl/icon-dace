// #define VERBOSE

/**
 * @file perf_toy_icon.c
 *
 * @copyright Copyright  (C)  2015 DKRZ, MPI-M
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

#include "yac_config.h"

#if defined YAC_NETCDF_ENABLED
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "utils.h"
#include "fields.h"
#include "yac_interface.h"
#include "geometry.h"
#include "yac_mpi.h"
#include "read_icon_grid.h"
#include "test_function.h"

// #define VTK_OUTPUT
#ifdef VTK_OUTPUT
#include "vtk_output.h"
#endif

/* ------------------------------------------------- */

/* For simplicity we define the same 8 fields that are in the
 * coupling configuration */

const char * fieldName[] = {"icon_out", "cube_out"};

// redefine YAC assert macros
#undef YAC_ASSERT
#define STR_USAGE "Usage: %s -c configFilename -g gridFilename\n"
#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, char const ** configFilename, char const ** gridFilename);

int main (int argc, char *argv[]) {

  // Initialisation of MPI

  MPI_Init (NULL, NULL);

  char const * configFilename = "perf_toy.yaml"; // default configuration file
  char const * gridFilename = "iconR2B06-grid.nc"; // default grid file
  parse_arguments(argc, argv, &configFilename, &gridFilename);

  xt_initialize(MPI_COMM_WORLD);

  /* The initialisation phase includes the reading of the
   * coupling configuration */
#ifdef VERBOSE
  printf (". main: calling yac_cinit\n");
#endif

  double tic, toc, time;

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  // move yaml file to node cache
  {
    FILE *f = fopen(configFilename, "rb");
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);  /* same as rewind(f); */

    char *string = xmalloc(fsize + 1);
    size_t dummy = fread(string, 1, fsize, f);
    UNUSED(dummy);
    fclose(f);
    free(string);
  }
  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  /* The usual component definition, here for two sequential components on the same process */

#ifdef VERBOSE
  printf (". main: calling yac_cdef_comp\n");
#endif

  int comp_id;
  char * comp_name = "ICON";

  yac_cdef_comp ( comp_name, &comp_id );
#ifdef VERBOSE
  printf ( ". main: defined %s with local comp ID %i \n", "ICON", comp_id );
#endif

  MPI_Comm local_comm;

  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;

  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);

  yac_mpi_call(MPI_Comm_free(&local_comm), MPI_COMM_WORLD);


  int cell_point_id;
  int corner_point_id;

  int field_ids[2];
  int grid_id;

  /* Grid definition for the first component (ICON) */

  int num_vertices;
  int num_cells;
  int * num_vertices_per_cell;
  int * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * cell_mask;
  int * global_cell_id;
  int * cell_core_mask;
  int * global_corner_id;
  int * corner_core_mask;

  read_part_icon_grid_information(gridFilename, &num_vertices, &num_cells,
                                  &num_vertices_per_cell, &cell_to_vertex,
                                  &x_vertices, &y_vertices, &x_cells,
                                  &y_cells, &global_cell_id,
                                  &cell_mask,
                                  &cell_core_mask, &global_corner_id,
                                  &corner_core_mask, rank, size);

  double * x_points, * y_points;

  x_points = x_vertices;
  y_points = y_vertices;

  double * x_center, * y_center;

  x_center = x_cells;
  y_center = y_cells;

  yac_cdef_grid_unstruct(
    "icon_grid", num_vertices, num_cells, num_vertices_per_cell,
    x_vertices, y_vertices, cell_to_vertex, &grid_id);

  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);

  yac_cdef_points_unstruct(
    grid_id, num_cells, YAC_LOCATION_CELL, x_center, y_center, &cell_point_id );
  yac_cdef_points_unstruct(
    grid_id, num_vertices, YAC_LOCATION_CORNER, x_points, y_points, &corner_point_id );

  /* Field definition for the first component (ICON-atmosphere) */

  yac_cdef_field(
    fieldName[0], comp_id, &cell_point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_ids[0]);
  yac_cdef_field(
    fieldName[1], comp_id, &cell_point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_ids[1]);

  toc=MPI_Wtime();
  time = toc-tic;
  printf ("ICON: Time for initialisation %f\n", time ); 

  /* Search. */

#ifdef VERBOSE
  printf (". main: calling yac_cenddef\n");
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  yac_cenddef ( );

  toc=MPI_Wtime();
  time = toc-tic;
  printf ("ICON: Time for search %f\n", time ); 

  double * conserv_in = xmalloc(num_cells * sizeof(*conserv_in));
  double * avg_in = xmalloc(num_vertices * sizeof(*avg_in));

  for (int i = 0; i < num_cells; ++i)
    conserv_in[i] = -10;
  for (int i = 0; i < num_vertices; ++i) avg_in[i] = -10;

  double * cell_out = xmalloc(num_cells * sizeof(*cell_out));
  double * corner_out = xmalloc(num_vertices * sizeof(*corner_out));

  int err, info;

  int cell_to_vertex_offset = 0;

  for (int i = 0; i < num_cells; ++i) {

    double middle_point[3] = {0, 0, 0};

    for (int j = 0; j < num_vertices_per_cell[i]; ++j) {

      double curr_point[3];

      LLtoXYZ(x_vertices[cell_to_vertex[cell_to_vertex_offset + j]],
              y_vertices[cell_to_vertex[cell_to_vertex_offset + j]],
              curr_point);

      middle_point[0] += curr_point[0];
      middle_point[1] += curr_point[1];
      middle_point[2] += curr_point[2];
    }

    double scale = 1.0 / sqrt(middle_point[0] * middle_point[0] + 
                              middle_point[1] * middle_point[1] + 
                              middle_point[2] * middle_point[2]);

    middle_point[0] *= scale;
    middle_point[1] *= scale;
    middle_point[2] *= scale;

    double lon, lat;

    XYZtoLL(middle_point, &lon, &lat);

    cell_to_vertex_offset += num_vertices_per_cell[i];

    cell_out[i] = test_func(lon, lat);
  }
  for (int i = 0; i < num_vertices; ++i)
    corner_out[i] = test_func(x_vertices[i], y_vertices[i]);

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  {
    double *point_set_data[1];
    double **collection_data[1] = {point_set_data};

    point_set_data[0] = cell_out;
    yac_cput(field_ids[0], 1, collection_data, &info, &err);
  }

  {
    double *collection_data[1] = {conserv_in};

    yac_cget(field_ids[1], 1, collection_data, &info, &err);
  }

  toc=MPI_Wtime();
  time = toc-tic;
  printf ("ICON: Time for ping-pong %f\n", time );

#ifdef VTK_OUTPUT
  //----------------------------------------------------------
  // write field to vtk output file
  //----------------------------------------------------------

  char vtk_filename[32];

  sprintf(vtk_filename, "perf_toy_icon_%d.vtk", rank);

  double point_data[num_vertices][3];
  for (int i = 0; i < num_vertices; ++i) {
   LLtoXYZ(x_vertices[i], y_vertices[i], point_data[i]);
  }

  VTK_FILE *vtk_file = vtk_open(vtk_filename, "unstruct_out");
  vtk_write_point_data(vtk_file, (double *)point_data, num_vertices);
  vtk_write_cell_data(vtk_file, (unsigned *)cell_to_vertex,
                      (unsigned*)num_vertices_per_cell, num_cells);
  vtk_write_point_scalars_int(
    vtk_file, corner_core_mask, num_vertices, "corner_core_mask");
  vtk_write_point_scalars_int(
    vtk_file, global_corner_id, num_vertices, "global_corner_id");
  vtk_write_cell_scalars_int(
    vtk_file, cell_core_mask, num_cells, "cell_core_mask");
  vtk_write_cell_scalars_int(
    vtk_file, global_cell_id, num_cells, "global_cell_id");

  vtk_write_cell_scalars_double(vtk_file, conserv_in, num_cells, "conserv_in");
  vtk_write_cell_scalars_double(vtk_file, cell_out, num_cells, "cell_out");
  vtk_write_point_scalars_double(vtk_file, avg_in, num_vertices, "avg_in");
  vtk_write_point_scalars_double(vtk_file, corner_out, num_vertices, "corner_out");

  vtk_close(vtk_file);

#endif // VTK_OUTPUT

  yac_cfinalize();

  xt_finalize();

  MPI_Finalize();

  free(num_vertices_per_cell);
  free(cell_to_vertex);
  free(x_cells);
  free(y_cells);
  free(x_vertices);
  free(y_vertices);
  free(cell_mask);
  free(global_cell_id);
  free(cell_core_mask);
  free(global_corner_id);
  free(corner_core_mask);
  free(conserv_in);
  free(cell_out);
  free(avg_in);
  free(corner_out);

  return EXIT_SUCCESS;
}



static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, char const ** gridFilename) {

  int opt;
  while ((opt = getopt(argc, argv, "c:g:")) != -1) {
    YAC_ASSERT((opt == 'c') || (opt == 'g'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
      case 'g':
        *gridFilename = optarg;
        break;
    }
  }
}

#else
#include <stdlib.h>
#include <stdio.h>
int main () {
  printf ("Examples requires compiling with NetCDF.\n");
  return EXIT_FAILURE;
}
#endif

