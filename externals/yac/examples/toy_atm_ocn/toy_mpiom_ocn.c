// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#undef VERBOSE

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "yac.h"
#include "yac_utils.h"

/* ------------------------------------------------- */

/* For simplicity we define the same 4 fields that are in the
 * coupling configuration */

#include "toy_common.h"

#define STR_USAGE "Usage: %s -c configFilename -g gridFilename\n"
#define YAC_ASSERT_ARGS(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, char const ** gridFilename);

int main (int argc, char *argv[]) {

  // Initialisation of MPI

  MPI_Init (0, NULL);

  /* The initialisation phase includes the reading of the
   * coupling configuration */
#ifdef VERBOSE
  printf (". main: calling yac_cinit\n");
#endif

  double tic, toc, time;
  double time_min, time_max, time_ave;
  double time_min_acc = 0.0, time_max_acc = 0.0, time_ave_acc = 0.0;
  
  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  char const * configFilename = "toy_atm_ocn.yaml"; // default configuration file
  //char const * gridFilename = "TP04_lsm.nc";
  //char const * gridFilename = "TP10_lsm.nc";
  //char const * gridFilename = "GR15_lsm.nc";
  char const * gridFilename = "GR30_lsm.nc";
  //char const * gridFilename = "TOY_lsm.nc";
  parse_arguments(argc, argv, &configFilename, &gridFilename);
  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  /* The usual component definition, here for two sequential components on the same process */

#ifdef VERBOSE
  printf (". main: calling yac_cdef_comp\n");
#endif

  int comp_id;
  char * comp_name = "OCEAN";

  yac_cdef_comp ( comp_name, &comp_id );
#ifdef VERBOSE
  printf ( ". main: defined %s with local comp ID %i \n", "toy-mpiom-ocn", comp_id );
#endif

  MPI_Comm local_comm;

  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;

  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);


  int point_id;

  int field_ids[4];
  int grid_id;

  /* Grid definition for the third component (toy-mpiom-ocean) */

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

  yac_read_part_mpiom_grid_information(gridFilename, &num_vertices, &num_cells,
                                       &num_vertices_per_cell, &cell_to_vertex,
                                       &x_vertices, &y_vertices, &x_cells,
                                       &y_cells, &global_cell_id,
                                       &cell_mask, &cell_core_mask,
                                       &global_corner_id, &corner_core_mask,
                                       rank, size);

  double * x_center, * y_center;

  x_center = x_cells;
  y_center = y_cells;

  yac_cdef_grid_unstruct ( "grid2", num_vertices, num_cells, num_vertices_per_cell,
                           x_vertices, y_vertices, cell_to_vertex, &grid_id );

  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);

  yac_cdef_points_unstruct(
    grid_id, num_cells, YAC_LOCATION_CELL, x_center, y_center, &point_id );

  cell_mask_own = malloc (num_cells * sizeof(*cell_mask_own));

  // water : 0 ; land  : 1

  for ( int i = 0; i < num_cells; ++i )
     cell_mask_own[i] = !cell_mask[i];

  yac_cset_mask ( cell_mask_own, point_id );

  /* Field definition for the ocean component (toy-mpiom-ocean) */

  for ( int i = 0; i < num_fields; i++ )
    yac_cdef_field(
      fieldName[i], comp_id, &point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
      &field_ids[i]);

  toc=MPI_Wtime();
  time = toc-tic;

  MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
  MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
  MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

  time_ave /= (double) size;

  if ( rank == 0 )
    printf ("toy-mpiom-ocn:       Time for initialisation %f %f %f \n", time_min, time_max, time_ave );

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

  time_ave /= size;

  if ( rank == 0 )
    printf ("toy-mpiom-ocn:       Time for search         %f %f %f \n", time_min, time_max, time_ave );

  double * cell_out_conserv_data = malloc(num_cells * sizeof(*cell_out_conserv_data));
  double * cell_out_hcsbb_data = malloc(num_cells * sizeof(*cell_out_hcsbb_data));
  double * cell_in_conserv_data = malloc(num_cells * sizeof(*cell_in_conserv_data));
  double * cell_in_hcsbb_data = malloc(num_cells * sizeof(*cell_in_hcsbb_data));

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

    cell_in_conserv_data[i]  = -10;
    cell_in_hcsbb_data[i]    = -10;
    cell_out_conserv_data[i] =
      (cell_out_hcsbb_data[i] = yac_test_func(lon, lat));
  }

  for ( int i = 0; i < num_cells; ++i ) {
    if ( !cell_mask_own[i] ) {
      cell_out_conserv_data[i] = -99999;
      cell_out_hcsbb_data[i]   = -99999;
      cell_in_conserv_data[i]  = -999;
      cell_in_hcsbb_data[i]    = -999;
    }
  }

  double point_data[num_vertices][3];
  for (int i = 0; i < num_vertices; ++i) {
   LLtoXYZ(x_vertices[i], y_vertices[i], point_data[i]);
  }

  for (int step = 0; step < num_steps; ++step) {

    MPI_Barrier(MPI_COMM_WORLD);

    tic=MPI_Wtime();

    {
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};

      point_set_data[0] = cell_out_conserv_data;
      yac_cput(field_ids[2], 1, collection_data, &info, &err);
      YAC_ASSERT(info, "check coupling period")
      point_set_data[0] = cell_out_hcsbb_data;
      yac_cput(field_ids[3], 1, collection_data, &info, &err);
      YAC_ASSERT(info, "check coupling period")
    }

    {
      double *collection_data[1];

      collection_data[0] = cell_in_conserv_data;
      yac_cget(field_ids[0], 1, collection_data, &info, &err);
      YAC_ASSERT(info, "check coupling period")
      collection_data[0] = cell_in_hcsbb_data;
      yac_cget(field_ids[1], 1, collection_data, &info, &err);
      YAC_ASSERT(info, "check coupling period")
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

    char vtk_filename[32];

    sprintf(vtk_filename, "toy-mpiom-ocean_out_%d_%d.vtk", rank, step);

    YAC_VTK_FILE *vtk_file = yac_vtk_open(vtk_filename, "toy-mpiom-ocean");
    yac_vtk_write_point_data(vtk_file, (double *)point_data, num_vertices);
    yac_vtk_write_cell_data(
      vtk_file, (unsigned *)cell_to_vertex, (unsigned*)num_vertices_per_cell,
      num_cells);
    yac_vtk_write_point_scalars_int(
      vtk_file, global_corner_id, num_vertices, "global_corner_id");
    yac_vtk_write_cell_scalars_int(
      vtk_file, global_cell_id, num_cells, "global_cell_id");
    yac_vtk_write_cell_scalars_int(
      vtk_file, cell_mask_own, num_cells, "cell_mask_own");
    yac_vtk_write_cell_scalars_int(
      vtk_file, cell_core_mask, num_cells, "cell_core_mask");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_conserv_data, num_cells, "cell_in_conserv");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_hcsbb_data, num_cells, "cell_in_hcsbb");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_out_conserv_data, num_cells, "cell_out_conserv");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_out_hcsbb_data, num_cells, "cell_out_hcsbb");
    yac_vtk_close(vtk_file);

    for ( int i = 0; i < num_cells; ++i ) {
      if ( !cell_mask_own[i] ) {
        cell_out_conserv_data[i] = -99999;
        cell_out_hcsbb_data[i]   = -99999;
        cell_in_conserv_data[i]  = -999;
        cell_in_hcsbb_data[i]    = -999;
      } else {
        cell_out_conserv_data[i] = cell_in_conserv_data[i];
        cell_out_hcsbb_data[i]   = cell_in_hcsbb_data[i];
        cell_in_conserv_data[i]  = -10;
        cell_in_hcsbb_data[i]    = -10;
      }
    }
  }

  if ( rank == 0 ) {
    time_min_acc /= (double) num_steps;
    time_max_acc /= (double) num_steps;
    time_ave_acc /= (double) num_steps;
    printf ("toy-mpiom-ocn:       Time for ping-pong      %f %f %f \n", time_min_acc, time_max_acc, time_ave_acc );
  }

  yac_cfinalize();

  MPI_Finalize();

  free(cell_mask_own);
  free(cell_core_mask);
  free(corner_core_mask);
  free(global_cell_id);
  free(global_corner_id);
  free(cell_in_conserv_data);
  free(cell_in_hcsbb_data);
  free(cell_out_conserv_data);
  free(cell_out_hcsbb_data);

  return EXIT_SUCCESS;
}

static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, char const ** gridFilename) {

  int opt;
  while ((opt = getopt(argc, argv, "c:g:")) != -1) {
    YAC_ASSERT_ARGS((opt == 'c') || (opt == 'g'), "invalid command argument")
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
