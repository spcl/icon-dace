//#define VERBOSE
/**
 * @file toy_icon_atm.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
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
#include "vtk_output.h"
#include "read_icon_grid.h"
#include "test_function.h"

/* ------------------------------------------------- */

const char * fieldName[] = {"runoff_field_atm2oce",
                            "runoff_field_oce2atm"};
const int num_fields = (int)(sizeof(fieldName)/sizeof(fieldName[0]));
const int num_steps = 1;

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

  MPI_Init (0, NULL);

  char const * configFilename = "toy_icon_runoff.yaml"; // default configuration file
  char const * gridFilename = "icon_grid_R02B03_G.nc";
  parse_arguments(argc, argv, &configFilename, &gridFilename);

  /* The initialisation phase includes the reading of the configuration file file */
#ifdef VERBOSE
  printf (". main: calling yac_cinit\n");
#endif

  double tic, toc, time;
  double time_min, time_max, time_ave;
  double time_min_acc, time_max_acc, time_ave_acc;

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  yac_cinit ();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_yaml(configFilename);

  /* The usual component definition, here for two sequential components on the same process */

#ifdef VERBOSE
  printf (". main: calling yac_cdef_comp\n");
#endif

  int comp_id;
  char * comp_name = "ATMOS";

  yac_cdef_comp ( comp_name, &comp_id );
#ifdef VERBOSE
  printf ( ". main: defined %s with local comp ID %i \n", "toy-icon-atmosphere", comp_id );
#endif

  MPI_Comm local_comm;

  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;

  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);

  int point_id;

  int field_ids[8];
  int grid_id;

  /* Grid definition for the first component (toy-icon-atmosphere) */

  int num_vertices;
  int num_cells;
  int * num_vertices_per_cell;
  int * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * cell_mask;
  int * cell_mask_own;
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

  double * x_center, * y_center;

  x_center = x_cells;
  y_center = y_cells;

  yac_cdef_grid_unstruct ( "atm_grid", num_vertices, num_cells, num_vertices_per_cell,
                           x_vertices, y_vertices, cell_to_vertex, &grid_id );

  yac_cset_global_index(global_cell_id, CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, CELL, grid_id);

  yac_cdef_points_unstruct ( grid_id, num_cells, CELL, x_center, y_center, &point_id );

  cell_mask_own = xmalloc (num_cells * sizeof(*cell_mask_own));

  // negative value represent sea, positive values represent land
  // -+2 inner, -+1 boundary
  //
  // Here we mask out land points:
  //
  // ICON mask values of 0 are left untouched!
  // In case the whole mask is 0 we set all values to true!

  int checksum_mask = 0;

  for ( int i = 0; i < num_cells; ++i )
     checksum_mask += abs(cell_mask[i]);

  if ( checksum_mask > 0 ) {
     for ( int i = 0; i < num_cells; ++i ) {
        cell_mask_own[i] = ( cell_mask[i] == -1 );
     }
  }
  else {
     for ( int i = 0; i < num_cells; ++i )
        cell_mask_own[i] = 1;
  }

#ifdef VERBOSE
  unsigned num_valid_cells = 0;
  for ( int i = 0; i < num_cells; ++i )
	  num_valid_cells += cell_mask_own[i] ;
  printf (". main: number of valid cells is %i\n", num_valid_cells );
#endif

  yac_cset_mask ( cell_mask_own, point_id );

  /* Field definition for the first component (toy-icon-atmosphere) */

  yac_cdef_field(
    fieldName[0], comp_id, &point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_ids[0]);
  yac_cdef_field(
    fieldName[1], comp_id, &point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_ids[1]);

  toc=MPI_Wtime();
  time = toc-tic;

  MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
  MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
  MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

  time_ave /= (double) size;

  if ( rank == 0 )
    printf ("toy-icon-atmosphere: Time for initialisation %f %f %f \n", time_min, time_max, time_ave );

  /* Search. */

#ifdef VERBOSE
  printf (". main: calling yac_cenddef\n");
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  yac_cenddef ( );

  toc=MPI_Wtime();
  time = toc-tic;

  MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
  MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
  MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

  time_ave /= (double) size;

  if ( rank == 0 )
    printf ("toy-icon-atmosphere: Time for search         %f %f %f \n", time_min, time_max, time_ave );

  double * cell_out_runoff_data = (double *) xmalloc(num_cells * sizeof(*cell_out_runoff_data));
  double * cell_in_runoff_data  = (double *) xmalloc(num_cells * sizeof(*cell_in_runoff_data));

  int err;
  int info;

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

    if ( cell_mask_own[i] ) {
      cell_in_runoff_data[i] = -10;
      cell_out_runoff_data[i] = test_func(lon, lat);
    } else {
      cell_in_runoff_data[i]  = -99999;
      cell_out_runoff_data[i] = -99999;
    }
  }

  double point_data[num_vertices][3];
  for (int i = 0; i < num_vertices; ++i) {
    LLtoXYZ(x_vertices[i], y_vertices[i], point_data[i]);
  }

  time_min_acc = 0.0;
  time_max_acc = 0.0;
  time_ave_acc = 0.0;

  for (int step = 0; step < num_steps; ++step) {

    MPI_Barrier(MPI_COMM_WORLD);

    tic=MPI_Wtime();

    {
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};

      point_set_data[0] = cell_out_runoff_data;
      yac_cput(field_ids[0], 1, collection_data, &info, &err);
      if ( info == 0 ) printf("WARNING: field is not coupled\n");
    }

    {
      double *collection_data[1];

      collection_data[0] = cell_in_runoff_data;
      yac_cget(field_ids[1], 1, collection_data, &info, &err);
      if ( info == 0 ) printf("WARNING: field is not coupled\n");
    }

    toc=MPI_Wtime();
    time = toc-tic;

    MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
    MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
    MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

    time_ave /= (double) size;

    if ( rank == 0 ) {
      time_min_acc += time_min;
      time_max_acc += time_max;
      time_ave_acc += time_ave;
    }

    //----------------------------------------------------------
    // write field to vtk output file
    //----------------------------------------------------------

    char vtk_filename[36];

    sprintf(vtk_filename, "toy-icon-atmosphere_out_%d_%d.vtk", rank, step);

    VTK_FILE *vtk_file = vtk_open(vtk_filename, "toy-icon-atmosphere_out");

    // mask out land points

    vtk_write_point_data(vtk_file, (double *)point_data, num_vertices);

    vtk_write_cell_data(vtk_file, (unsigned *)cell_to_vertex,
                        (unsigned*)num_vertices_per_cell, num_cells);
    vtk_write_point_scalars_int(vtk_file, corner_core_mask, num_vertices,
                                "corner_core_mask");
    vtk_write_point_scalars_int(vtk_file, global_corner_id, num_vertices,
                                "global_corner_id");
    vtk_write_cell_scalars_int(vtk_file, cell_core_mask, num_cells,
                               "cell_core_mask");
    vtk_write_cell_scalars_int(vtk_file, global_cell_id, num_cells,
                               "global_cell_id");

    vtk_write_cell_scalars_int(vtk_file, cell_mask, num_cells, "mask");
    vtk_write_cell_scalars_int(vtk_file, cell_mask_own, num_cells, "cell_mask_own");
    vtk_write_cell_scalars_double(vtk_file, cell_in_runoff_data, num_cells, "cell_in_runoff");
    vtk_write_cell_scalars_double(vtk_file, cell_out_runoff_data, num_cells, "cell_out_runoff");
    vtk_close(vtk_file);

    for ( int i = 0; i < num_cells; ++i ) {
      if ( cell_mask_own[i] ) {
        cell_out_runoff_data[i] = cell_in_runoff_data[i];
        cell_in_runoff_data[i] = -10;
      } else {
        cell_out_runoff_data[i] = -99999;
        cell_in_runoff_data[i]  = -99999;
      }
    }
  }

  if ( rank == 0 ) {
    time_min_acc /= (double) num_steps;
    time_max_acc /= (double) num_steps;
    time_ave_acc /= (double) num_steps;
    printf ("toy-icon-atmosphere: Time for ping-pong      %f %f %f \n", time_min_acc, time_max_acc, time_ave_acc );
  }

  free (cell_mask_own);

  delete_icon_grid_data ( &cell_mask,
                          &global_cell_id,
                          &cell_core_mask,
                          &num_vertices_per_cell,
                          &global_corner_id,
                          &corner_core_mask,
                          &cell_to_vertex,
                          &x_cells,
                          &y_cells,
                          &x_vertices,
                          &y_vertices );

  yac_cfinalize();

  MPI_Finalize();

  free(cell_in_runoff_data);
  free(cell_out_runoff_data);

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
