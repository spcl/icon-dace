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

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "yac_interface.h"
#include "utils.h"

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

  const int no_of_fields = 14;

  const char * fieldName[] = { "surface_downward_eastward_stress",  // bundled field containing two components
                               "surface_downward_northward_stress", // bundled field containing two components
                               "surface_fresh_water_flux",          // bundled field containing three components
                               "surface_temperature",
                               "total_heat_flux",                   // bundled field containing four components
                               "atmosphere_sea_ice_bundle",         // bundled field containing four components
                               "sea_surface_temperature",
                               "eastward_sea_water_velocity",
                               "northward_sea_water_velocity",
                               "ocean_sea_ice_bundle",
                               "ocean_out1",                        // output field
                               "ocean_out2",                        // output field
                               "ocean_out3",                        // output field
                               "ocean_out4"};                       // output field
  int field_collection_size[] = {
    2,
    2,
    3,
    1,
    4,
    4,
    1,
    1,
    1,
    5,
    1,
    1,
    1,
    1};

  char * comp_name = "dummy_ocean";
  char * grid_name = "dummy_ocean_grid";

  const int nbr_cells = NBR_CELLS;
  const int nbr_vertices = NBR_VERTICES;
  int nbr_vertices_per_cell[NBR_CELLS];

  int i, info, ierror;
  int comp_id;

  int cell_point_ids[1];

  int grid_id;

  int global_index[NBR_CELLS];
  int cell_core_mask[NBR_CELLS];

  double * buffer;
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

  buffer_lon = xmalloc ( nbr_vertices * sizeof(*buffer_lon) );
  buffer_lat = xmalloc ( nbr_vertices * sizeof(*buffer_lat) );
  cell_to_vertex = xmalloc ( 3 * nbr_cells * sizeof(*cell_to_vertex) );

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

  buffer_lon[0] =  0.0 * YAC_RAD; buffer_lat[0] =  1.0 * YAC_RAD;
  buffer_lon[1] = -1.0 * YAC_RAD; buffer_lat[1] =  0.0 * YAC_RAD;
  buffer_lon[2] =  1.0 * YAC_RAD; buffer_lat[2] =  0.0 * YAC_RAD;
  buffer_lon[3] =  0.0 * YAC_RAD; buffer_lat[3] = -1.0 * YAC_RAD;

  // Connectivity

  cell_to_vertex[0] = 0; cell_to_vertex[1] = 1; cell_to_vertex[2] = 2; // cell 1
  cell_to_vertex[3] = 1; cell_to_vertex[4] = 3; cell_to_vertex[5] = 2; // cell 2

  for ( i = 0; i < nbr_cells; ++i ) nbr_vertices_per_cell[i] = 3;

  // Define unstructured grid

  yac_cdef_grid_unstruct(
    grid_name, nbr_vertices, nbr_cells, nbr_vertices_per_cell,
    buffer_lon, buffer_lat, cell_to_vertex, &grid_id);

  // Decomposition information

  for (int i = 0; i < nbr_cells; ++i ) {
    global_index[i] = i;
    cell_core_mask[i] = 1;
  }

  yac_cset_global_index(global_index, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);

  // Center points in cells (needed e.g. for nearest neighbour)

  buffer_lon[0] = 0.0 * YAC_RAD; buffer_lat[0] =  0.5 * YAC_RAD;
  buffer_lon[1] = 0.0 * YAC_RAD; buffer_lat[1] = -0.5 * YAC_RAD;

  yac_cdef_points_unstruct(
    grid_id, nbr_cells, YAC_LOCATION_CELL, buffer_lon, buffer_lat,
    &cell_point_ids[0]);

  free ( buffer_lon );
  free ( buffer_lat );

  // Mask generation

  cell_mask = xmalloc (nbr_cells * sizeof(*cell_mask));

  for (unsigned i = 0; i < nbr_cells; ++i) cell_mask[i] = 1;

  yac_cset_mask(cell_mask, cell_point_ids[0]);

  free (cell_mask);

  field_id = xmalloc ( no_of_fields * sizeof(*field_id));

  for ( i = 0; i < no_of_fields; ++i )
    yac_cdef_field ( fieldName[i],
                     comp_id,
                     cell_point_ids,
                     1,field_collection_size[i],
                     "1",YAC_TIME_UNIT_SECOND,
                     &field_id[i] );

  yac_cenddef ( );

  // these queries can only be fulfilled once the search has been completed.

  for ( i = 0; i < no_of_fields; ++i ) {
    const char* timestep_string = yac_cget_timestep_from_field_id ( field_id[i] );
    printf ( "Field ID %d: model step %s \n", field_id[i], timestep_string);
  }

  // Data exchange

  buffer = xmalloc (5 * nbr_cells * sizeof (*buffer));

  /* field_id[0] represents "TAUX"   wind stress component
     field_id[1] represents "TAUY"   wind stress component
     field_id[2] represents "SFWFLX" surface fresh water flux
     field_id[3] represents "SFTEMP" surface temperature
     field_id[4] represents "THFLX"  total heat flux
     field_id[5] represents "ICEATM" ice temperatures and melt potential

     field_id[6] represents "SST"    sea surface temperature
     field_id[7] represents "OCEANU" u component of ocean surface current
     field_id[8] represents "OCEANV" v component of ocean surface current
     field_id[9]represents "ICEOCE" ice thickness, concentration and temperatures
  */

  // Send fields from ocean to atmosphere
  // ------------------------------------

  {
    double *point_set_data[5];
    double **collection_data[5];

    for (int j = 0; j < 5; ++j) {
      point_set_data[j] = buffer + j * nbr_cells;
      collection_data[j] = &(point_set_data[j]);
    }

    // SST
    for ( i = 0; i < nbr_cells; ++i ) buffer[0*nbr_cells+i] = 110.1;

    yac_cput ( field_id[6], 1, collection_data, &info, &ierror );
    if ( info > 0 ) printf ( "ocean CPL SST %f \n", buffer[0*nbr_cells+0]);

    // zonal velocity

    for ( i = 0; i < nbr_cells; ++i ) buffer[0*nbr_cells+i] = 120.1;

    yac_cput ( field_id[7], 1, collection_data, &info, &ierror );
    if ( info > 0 ) printf ( "ocean CPL U %f \n", buffer[0*nbr_cells+0]);

    // meridional velocity

    for ( i = 0; i < nbr_cells; ++i ) buffer[0*nbr_cells+i] = 130.1;

    yac_cput ( field_id[8], 1, collection_data, &info, &ierror );
    if ( info > 0 ) printf ( "ocean CPL V %f \n", buffer[0*nbr_cells+0]);

    // Ice thickness, concentration, T1 and T2

    for ( i = 0; i < nbr_cells; ++i ) {
      buffer[0*nbr_cells+i] = 140.1;
      buffer[1*nbr_cells+i] = 140.2;
      buffer[2*nbr_cells+i] = 140.3;
      buffer[3*nbr_cells+i] = 140.4;
      buffer[4*nbr_cells+i] = 140.5;
    }
    point_set_data[0] = buffer;

    yac_cput ( field_id[9], 5, collection_data, &info, &ierror );
    if ( info > 0 )
      printf (
        "ocean CPL ice 1 %f \n"
        "ocean CPL ice 2 %f \n"
        "ocean CPL ice 3 %f \n"
        "ocean CPL ice 4 %f \n"
        "ocean CPL ice 5 %f \n",
        buffer[0*nbr_cells+0], buffer[1*nbr_cells+0], buffer[2*nbr_cells+0],
        buffer[3*nbr_cells+0], buffer[4*nbr_cells+0]);
  }

  // Receive fields from atmosphere
  // ------------------------------

  {
    double *collection_data[5];

    for (int j = 0; j < 5; ++j) {
      collection_data[j] = buffer + j * nbr_cells;
    }

    // zonal wind stress

    yac_cget ( field_id[0], 2, collection_data, &info, &ierror );

    // meridional wind stress

    yac_cget ( field_id[1], 2, collection_data, &info, &ierror );

    // freshwater flux

    yac_cget ( field_id[2], 3, collection_data, &info, &ierror );

    // surface air temperature
    yac_cget ( field_id[3], 1, collection_data, &info, &ierror );

    // total heat flux - 4 parts - record 5

    yac_cget ( field_id[4], 4, collection_data, &info, &ierror );

    // ice parameter

    yac_cget ( field_id[5], 4, collection_data, &info, &ierror );
  }

  free (buffer);
  free (field_id);

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
