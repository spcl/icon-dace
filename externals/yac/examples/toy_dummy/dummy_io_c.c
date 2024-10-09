// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "yac.h"

#define NBR_CELLS 2
#define NBR_VERTICES 4
#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

#define STR_USAGE "Usage: %s -c configFilename\n"
#define YAC_ASSERT_ARGS(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, char const ** configFilename);

int main (int argc, char *argv[]) {

   int no_of_fields = 8;

   const char * fieldName[] = { "atmos_out1",
                                "atmos_out2",
                                "atmos_out3",
                                "atmos_out4",
                                "ocean_out1",
                                "ocean_out2",
                                "ocean_out3",
                                "ocean_out4"};

  char * comp_name = "dummy_io";

  char * grid_name[] = { "ocean_grid",
                         "atmos_grid"};

  const int nbr_cells = NBR_CELLS;
  const int nbr_vertices = NBR_VERTICES;
  int nbr_vertices_per_cell[NBR_CELLS];

  int i;
  int comp_id;

  int cell_point_ids[2];

  int grid_ids[2];

  int global_index[NBR_CELLS];
  int cell_core_mask[NBR_CELLS];

  double * buffer_lon;
  double * buffer_lat;
  int * cell_to_vertex;

  int * cell_mask;
  int * field_id;

  MPI_Comm local_comm;

  int size, rank;

  MPI_Init (0, NULL);

  // Initialise the coupler
  char const * configFilename = "toy_dummy.yaml"; // default configuration file
  parse_arguments(argc, argv, &configFilename);
  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  // Inform the coupler about what we are
  yac_cdef_comp ( comp_name, &comp_id );

  printf ( "YAC Version %s\n", yac_cget_version() );

  yac_cget_comp_comm ( comp_id, &local_comm );

  MPI_Comm_rank ( local_comm, &rank );
  MPI_Comm_size ( local_comm, &size );

  printf (" %s rank %d : local size %d \n", comp_name, rank, size );

  /* Here we define two grids (although numrically identical) to mimick
     two different Output server groups */

  for ( int igrid = 0; igrid < 2; igrid++ ) {

    buffer_lon = malloc ( nbr_vertices * sizeof(*buffer_lon) );
    buffer_lat = malloc ( nbr_vertices * sizeof(*buffer_lat) );
    cell_to_vertex = malloc ( 3 * nbr_cells * sizeof(*cell_to_vertex) );

    /* Define vertices

                0
               / \
              / o \
             /     \
            1-------2   Eq.
             \     /
              \ o /
               \ /
                3

    */

    buffer_lon[0] =  0.0 * YAC_RAD; buffer_lat[0] =  1.1 * YAC_RAD;
    buffer_lon[1] = -1.0 * YAC_RAD; buffer_lat[1] =  0.0 * YAC_RAD;
    buffer_lon[2] =  1.0 * YAC_RAD; buffer_lat[2] =  0.0 * YAC_RAD;
    buffer_lon[3] =  0.0 * YAC_RAD; buffer_lat[3] = -1.0 * YAC_RAD;

    // Connectivity

    cell_to_vertex[0] = 0; cell_to_vertex[1] = 1; cell_to_vertex[2] = 2; // cell 1
    cell_to_vertex[3] = 1; cell_to_vertex[4] = 3; cell_to_vertex[5] = 2; // cell 2

    for ( i = 0; i < nbr_cells; ++i ) nbr_vertices_per_cell[i] = 3;

    // Define unstructured grid

    yac_cdef_grid_unstruct ( grid_name[igrid],
                             nbr_vertices,
                             nbr_cells,
                             nbr_vertices_per_cell,
                             buffer_lon,
                             buffer_lat,
                             cell_to_vertex,
                             &grid_ids[igrid] );

  // Decomposition information

  for (int i = 0; i < nbr_cells; ++i ) {
    global_index[i] = i;
    cell_core_mask[i] = 1;
  }

  yac_cset_global_index(global_index, YAC_LOCATION_CELL, grid_ids[igrid]);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_ids[igrid]);

    // Center points in cells (needed e.g. for nearest neighbour)

    buffer_lon[0] = 0.0 * YAC_RAD; buffer_lat[0] =  0.5 * YAC_RAD;
    buffer_lon[1] = 0.0 * YAC_RAD; buffer_lat[1] = -0.5 * YAC_RAD;

    yac_cdef_points_unstruct ( grid_ids[igrid],
                               nbr_cells,
                               YAC_LOCATION_CELL,
                               buffer_lon,
                               buffer_lat,
                               &cell_point_ids[igrid] );

    free ( buffer_lon );
    free ( buffer_lat );

    // Mask generation

    cell_mask = malloc (nbr_cells * sizeof(*cell_mask));

    for ( i = 0; i < nbr_cells; ++i )
      cell_mask[i] = 1;

    yac_cset_mask ( cell_mask,
                    cell_point_ids[igrid] );

    free (cell_mask);
  }

  field_id = malloc ( no_of_fields * sizeof(*field_id));

  // atmosphere output

  for ( i = 0; i < no_of_fields-4; i++ )
    yac_cdef_field ( fieldName[i],
                     comp_id,
                     &cell_point_ids[1],
                     1,1,"1",YAC_TIME_UNIT_SECOND,
                     &field_id[i] );

  // ocean output

  for ( i = no_of_fields-4; i < no_of_fields; i++ )
    yac_cdef_field ( fieldName[i],
                     comp_id,
                     &cell_point_ids[0],
                     1,1,"1",YAC_TIME_UNIT_SECOND,
                     &field_id[i] );

  yac_cenddef ( );

  yac_cfinalize ();

  MPI_Finalize ();

  return 0;
}

static void parse_arguments(
  int argc, char ** argv, char const ** configFilename) {

  int opt;
  while ((opt = getopt(argc, argv, "c:")) != -1) {
    YAC_ASSERT_ARGS((opt == 'c'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
    }
  }
}
