/**
 * @file toy_coupling.c
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

#define EXACT

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "yac_interface.h"
#include "geometry.h"
#include "vtk_output.h"
#include "utils.h"

const char * fieldName[] = { "AtoB",
                             "BtoA"};

const int no_of_fields = 2;

int info, ierror;
int comp_id;
int comp_ids[1];

int cell_point_id;

int grid_id;

int nbr_fields;

double * recv_buffer;
double * send_buffer;

double * buffer = NULL;

int * global_index;
int * cell_core_mask;

int * field_id;

static void dummy_compA (char const * configFilename, char const * gridPath);
static void dummy_compB (char const * configFilename);

#define STR_USAGE "Usage: %s -c configFilename -g gridPath\n"
#define YAC_ASSERT_ARGS(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, char const ** configFilenamem, char const ** gridPath);

/* -------------------------------------------------------------------- */

int main (int argc, char *argv[]) {
  int size, rank;

  MPI_Init (0, NULL);

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  YAC_ASSERT(size == 2, "Too many processes have been launched")

  char const * configFilename = "toy_coupling.yaml"; // default configuration file
  char const * gridPath = ".";
  parse_arguments(argc, argv, &configFilename, &gridPath);

  if (rank == 0) dummy_compA(configFilename, gridPath);
  else           dummy_compB(configFilename);

  MPI_Finalize ();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------- */

struct point_with_index {

  struct {
    double lon, lat;
  } p;
  int i;
};

static int compare_point_with_index(const void * a,const void * b) {

   const struct point_with_index * a_ = (struct point_with_index*)a;
   const struct point_with_index * b_ = (struct point_with_index*)b;

  int lon_diff = fabs(a_->p.lon - b_->p.lon) > 1e-9;
  int lat_diff = fabs(a_->p.lat - b_->p.lat) > 1e-9;

  if (lon_diff) {

    if (a_->p.lon > b_->p.lon) return -1;
    else return 1;

  } else if (lat_diff) {

    if (a_->p.lat > b_->p.lat) return -1;
    else return 1;
  } else
    return 0;
}

static void remove_duplicated_vertices(double * temp_vertex_lon,
                                       double * temp_vertex_lat,
                                       int * temp_nbr_vertices,
                                       int * old_to_new_id) {

   struct point_with_index * sort_array =
      xmalloc(*temp_nbr_vertices * sizeof(*sort_array));

   for (int i = 0; i < *temp_nbr_vertices; ++i) {

      double curr_lon, curr_lat;

      curr_lon = temp_vertex_lon[i];
      curr_lat = temp_vertex_lat[i];

      while (curr_lon < 0.0) curr_lon += 2.0 * M_PI;
      while (curr_lon >= 2.0 * M_PI) curr_lon -= 2.0 * M_PI;

      sort_array[i].p.lon = curr_lon;
      sort_array[i].p.lat = curr_lat;

      sort_array[i].i = i;
   }

   yac_mergesort(sort_array, *temp_nbr_vertices, sizeof(*sort_array),
                 compare_point_with_index);

   old_to_new_id[sort_array[0].i] = 1;

   int last_unique_idx = sort_array[0].i;

   for (int i = 1; i < *temp_nbr_vertices; ++i) {

      if (compare_point_with_index(sort_array + i - 1, sort_array + i)) {

         old_to_new_id[sort_array[i].i] = 1;
         last_unique_idx = sort_array[i].i;

      } else {

         old_to_new_id[sort_array[i].i] = -last_unique_idx;
      }
   }

   free(sort_array);

   size_t new_nbr_vertices = 0;

   for (int i = 0; i < *temp_nbr_vertices; ++i) {

      if (old_to_new_id[i] == 1) {

         temp_vertex_lon[new_nbr_vertices] = temp_vertex_lon[i];
         temp_vertex_lat[new_nbr_vertices] = temp_vertex_lat[i];

         old_to_new_id[i] = new_nbr_vertices;

         new_nbr_vertices++;
      }
   }

   for (int i = 0; i < *temp_nbr_vertices; ++i)
      if (old_to_new_id[i] <= 0)
         old_to_new_id[i] = old_to_new_id[-old_to_new_id[i]];

   *temp_nbr_vertices = new_nbr_vertices;
}

static void read_compA_grid_data(char const * gridPath,
                                 int * nbr_vertices, int * nbr_cells,
                                 int ** num_vertices_per_cell, int ** cell_to_vertex,
                                 double ** buffer_lonv, double ** buffer_latv,
                                 double ** buffer_lon, double ** buffer_lat) {

  *nbr_cells = 20;
#ifdef EXACT
  int nbr_vertices_per_cell = 12;
#else
  int nbr_vertices_per_cell = 7;
#endif

  // Read in grid data

  char * gridFilename =
#ifdef EXACT
  "mesh_compA_exact.txt";
#else
  "mesh_compA_convex.txt";
#endif

  char * full_gridFilename = xmalloc(strlen(gridPath) + strlen(gridFilename) + 2);
  sprintf(full_gridFilename, "%s/%s", gridPath, gridFilename);

  char* inFormat0 = "%s %*s %*s";
  char* inFormat1 = "%*10c %lf";
#ifdef EXACT
  char* inFormat2 = "%*10c %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf";
  enum {num_entries_in_inFormat2 = 12};
#else
  char* inFormat2 = "%*10c %lf %lf %lf %lf %lf %lf %lf";
  enum {num_entries_in_inFormat2 = 7};
#endif

  FILE *inData = fopen(full_gridFilename, "r");
  YAC_ASSERT_F(
    inData, "Cannot open file \"%s\"", full_gridFilename);
  free(full_gridFilename);

  char dummy_string[32];

  //skip first three lines

  for ( int i = 0; i < 3; ++i )
    YAC_ASSERT(1 == fscanf (inData, inFormat0, dummy_string),
      "fscanf failed reading data");

  *nbr_vertices = nbr_vertices_per_cell * (*nbr_cells);

  *buffer_lon = xmalloc ( *nbr_cells * sizeof(*buffer_lon) );
  *buffer_lat = xmalloc ( *nbr_cells * sizeof(*buffer_lat) );
  *buffer_lonv = xmalloc ( (*nbr_vertices) * sizeof(*buffer_lonv) );
  *buffer_latv = xmalloc ( (*nbr_vertices) * sizeof(*buffer_latv) );

  for ( int i = 0; i < *nbr_cells; ++i ) {
    YAC_ASSERT(1 == fscanf ( inData, inFormat1, &((*buffer_lon)[i]) ),
               "fscanf failed reading data")
    printf ("point buffer_lat %i : %lf \n", i, (*buffer_lon)[i]);
  }

  for ( int i = 0; i < *nbr_cells; ++i ) {
    YAC_ASSERT(
      num_entries_in_inFormat2 ==
      fscanf ( inData, inFormat2,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  0,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  1,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  2,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  3,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  4,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  5,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  6,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  7,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  8,
               (*buffer_lonv) + i*nbr_vertices_per_cell +  9,
               (*buffer_lonv) + i*nbr_vertices_per_cell + 10,
               (*buffer_lonv) + i*nbr_vertices_per_cell + 11),
      "fscanf failed reading data")
  }

  for ( int i = 0; i < *nbr_cells; ++i ) {
    YAC_ASSERT (1 == fscanf ( inData, inFormat1, &((*buffer_lat)[i]) ),
                "fscanf failed reading data")
    printf ("point buffer_lat %i : %lf \n", i, (*buffer_lat)[i]);
  }

  for ( int i = 0; i < *nbr_cells; ++i ) {
    YAC_ASSERT(
      num_entries_in_inFormat2 ==
      fscanf ( inData, inFormat2,
               (*buffer_latv) + i*nbr_vertices_per_cell +  0,
               (*buffer_latv) + i*nbr_vertices_per_cell +  1,
               (*buffer_latv) + i*nbr_vertices_per_cell +  2,
               (*buffer_latv) + i*nbr_vertices_per_cell +  3,
               (*buffer_latv) + i*nbr_vertices_per_cell +  4,
               (*buffer_latv) + i*nbr_vertices_per_cell +  5,
               (*buffer_latv) + i*nbr_vertices_per_cell +  6,
               (*buffer_latv) + i*nbr_vertices_per_cell +  7,
               (*buffer_latv) + i*nbr_vertices_per_cell +  8,
               (*buffer_latv) + i*nbr_vertices_per_cell +  9,
               (*buffer_latv) + i*nbr_vertices_per_cell + 10,
               (*buffer_latv) + i*nbr_vertices_per_cell + 11),
      "fscanf failed reading data")
  }

  for ( int i = 0; i < *nbr_vertices; ++i )
    printf ("%i vertex buffer_lonv : %lf   buffer_latv : %lf \n", i,
            (*buffer_lonv)[i], (*buffer_latv)[i]);

  for ( int i = 0; i < *nbr_cells; ++i ) {
    (*buffer_lon)[i] *= YAC_RAD;
    (*buffer_lat)[i] *= YAC_RAD;
  }

  for ( int i = 0; i < *nbr_vertices; ++i ) {
    (*buffer_lonv)[i] *= YAC_RAD;
    (*buffer_latv)[i] *= YAC_RAD;
  }

  // Connectivity

  (*cell_to_vertex) = xmalloc ( (*nbr_cells) * nbr_vertices_per_cell *
                                sizeof(**cell_to_vertex) );
  (*num_vertices_per_cell) = xmalloc ( *nbr_cells * sizeof(**num_vertices_per_cell) );

  int * old_to_new_id = xmalloc((*nbr_vertices) * sizeof(*old_to_new_id));
  remove_duplicated_vertices(*buffer_lonv, *buffer_latv, nbr_vertices,
                             old_to_new_id);

  int count = 0;

  for (int i = 0; i < *nbr_cells; ++i) {

    (*cell_to_vertex)[count++] = old_to_new_id[i * nbr_vertices_per_cell];
    (*num_vertices_per_cell)[i] = 1;

    for (int j = 1; j < nbr_vertices_per_cell; ++j) {

      if (old_to_new_id[i * nbr_vertices_per_cell + j - 1] !=
          old_to_new_id[i * nbr_vertices_per_cell + j]) {

        (*cell_to_vertex)[count++] = old_to_new_id[i * nbr_vertices_per_cell + j];
        (*num_vertices_per_cell)[i]++;
      }
    }
  }

  fclose (inData);

  free(old_to_new_id);
}

static void dummy_compA (char const * configFilename, char const * gridPath) {

  int size, rank;

  MPI_Comm local_comm;

  char * comp_name = "dummy_compA";
  char * grid_name = "dummy_compA_grid";

  // Initialise the coupler
  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_yaml(configFilename);

  // Inform the coupler about what we are
  yac_cdef_comp ( comp_name, &comp_id );
  comp_ids[0] = comp_id;

  yac_cget_comp_comm ( comp_id, &local_comm );

  MPI_Comm_rank ( local_comm, &rank );
  MPI_Comm_size ( local_comm, &size );

  int nbr_vertices, nbr_cells, * num_vertices_per_cell, * cell_to_vertex;
  double * buffer_lonv, * buffer_latv, * buffer_lon, * buffer_lat;


  read_compA_grid_data(gridPath, &nbr_vertices, &nbr_cells, &num_vertices_per_cell,
                       &cell_to_vertex, &buffer_lonv, &buffer_latv,
                       &buffer_lon, &buffer_lat);

  // Define of unstructured grid

  yac_cdef_grid_unstruct(
    grid_name, nbr_vertices, nbr_cells, num_vertices_per_cell,
    buffer_lonv, buffer_latv, cell_to_vertex, &grid_id);

  // Decomposition information

  global_index = xmalloc ( nbr_cells * sizeof(*global_index) ) ;
  cell_core_mask  = xmalloc ( nbr_cells * sizeof(*cell_core_mask) );

  for ( int i = 0; i < nbr_cells; ++i ) {
    global_index[i] = i;
    cell_core_mask[i] = 1;
  }

  yac_cset_global_index(global_index, CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, CELL, grid_id);

  // Center points in cells

  yac_cdef_points_unstruct ( grid_id,
                             nbr_cells,
                             CELL,
                             buffer_lon,
                             buffer_lat,
                             &cell_point_id );

  // Define fields

  field_id = xmalloc ( no_of_fields * sizeof(*field_id));

  for ( int i = 0; i < no_of_fields; ++i )
    yac_cdef_field ( fieldName[i],
                     comp_id,
                     &cell_point_id,
                     1,1,"600",YAC_TIME_UNIT_SECOND,
                     &field_id[i] );

  // Invoke the search

  yac_cenddef ( );

  // Data exchange

  recv_buffer = xmalloc (nbr_cells * sizeof (*recv_buffer) );
  send_buffer = xmalloc (nbr_cells * sizeof (*send_buffer) );

  for ( int i = 0; i < nbr_cells; ++i ) recv_buffer[i] = 20.0;
  for ( int i = 0; i < nbr_cells; ++i ) send_buffer[i] = (double) i;

  // Send fields to compB
  // --------------------

  {
    double *point_set_data[1];
    double **collection_data[1];

    point_set_data[0] = send_buffer;
    collection_data[0] = &(point_set_data[0]);

    // AtoB

    yac_cput ( field_id[0], 1, collection_data, &info, &ierror );

  }

/*
  //
  // Receive fields from compB
  // -------------------------

  {
    double *collection_data[1];

    collection_data[0] = recv_buffer;

    // BtoA

    yac_cget ( field_id[1], 1, collection_data, &info, &ierror );

  }
*/

  //----------------------------------------------------------
  // write field to vtk output file
  //----------------------------------------------------------

  char vtk_filename[32];

  sprintf(vtk_filename, "compA_out_%d.vtk", rank);

  double point_data[nbr_vertices][3];
  for (int i = 0; i < nbr_vertices; ++i)
   LLtoXYZ(buffer_lonv[i], buffer_latv[i], point_data[i]);

  VTK_FILE *vtk_file = vtk_open(vtk_filename, "compA_out");

  vtk_write_point_data(vtk_file, (double *)point_data, nbr_vertices);

  vtk_write_cell_data(vtk_file, (unsigned *)cell_to_vertex,
                      (unsigned*)num_vertices_per_cell, nbr_cells);

  vtk_write_cell_scalars_int(vtk_file, global_index, nbr_cells,
                             "global_cell_id");

  vtk_write_cell_scalars_double(vtk_file, recv_buffer, nbr_cells, "cell_in");

  vtk_close(vtk_file);

  /* Free memory and finalise */

  free (cell_to_vertex);
  free (num_vertices_per_cell);

  free ( buffer );

  free ( buffer_lon );
  free ( buffer_lat );

  free ( buffer_lonv );
  free ( buffer_latv );

  free( global_index );
  free( cell_core_mask );

  free (recv_buffer);
  free (send_buffer);

  free (field_id);

  yac_cfinalize ();

}

/* -------------------------------------------------------------------- */

static void dummy_compB (char const * configFilename) {

  int size, rank;

  MPI_Comm local_comm;

  char * comp_name = "dummy_compB";
  char * grid_name = "dummy_compB_grid";

  // Initialise the coupler
  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_yaml(configFilename);

  // Inform the coupler about what we are
  yac_cdef_comp ( comp_name, &comp_id );
  comp_ids[0] = comp_id;

  yac_cget_comp_comm ( comp_id, &local_comm );

  MPI_Comm_rank ( local_comm, &rank );
  MPI_Comm_size ( local_comm, &size );

  /* Define compB grid */

  int nbr_vertices[2];
  int nbr_cells[2] = {200, 320};
  int cyclic[2] = {0,0};

  double dx = 0.01 * YAC_RAD;

  nbr_vertices[0] = nbr_cells[0] + 1;
  nbr_vertices[1] = nbr_cells[1] + 1;

  double * buffer_lonv = xmalloc ( nbr_vertices[0] * sizeof(*buffer_lonv) );
  double * buffer_latv = xmalloc ( nbr_vertices[1] * sizeof(*buffer_latv) );

  // Description of reg2d elements

  buffer_lonv[0] = 265.6 * YAC_RAD - 0.5 * dx;
  buffer_latv[0] =  58.5 * YAC_RAD - 0.5 * dx;

  for ( int i = 1; i < nbr_vertices[0]; ++i )
    buffer_lonv[i] = buffer_lonv[i-1] + dx;

  for ( int i = 1; i < nbr_vertices[1]; ++i )
    buffer_latv[i] = buffer_latv[i-1] + dx;

  yac_cdef_grid_reg2d ( grid_name,
                        nbr_vertices,
                        cyclic,
                        buffer_lonv,
                        buffer_latv,
                        &grid_id );

  // Decomposition information

  int total_nbr_cells =  nbr_cells[0] * nbr_cells[1];

  global_index = xmalloc ( total_nbr_cells * sizeof(*global_index) ) ;
  cell_core_mask  = xmalloc ( total_nbr_cells * sizeof(*cell_core_mask) );

  for ( int j = 0; j < nbr_cells[1]; ++j ) {
     for ( int i = 0; i < nbr_cells[0]; ++i ) {
       global_index[j*nbr_cells[0]+i] = j*nbr_cells[0]+i;
       cell_core_mask[j*nbr_cells[0]+i] = 1;
     }
  }

  yac_cset_global_index(global_index, CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, CELL, grid_id);

  // Center points in cells

  double * buffer_lon = xmalloc ( nbr_cells[0] * sizeof(*buffer_lon) );
  double * buffer_lat = xmalloc ( nbr_cells[1] * sizeof(*buffer_lat) );

  buffer_lon[0] = 265.6 * YAC_RAD;
  buffer_lat[0] =  58.5 * YAC_RAD;

  for ( int i = 1; i < nbr_cells[0]; ++i )
    buffer_lon[i] = buffer_lonv[i-1] + dx;

  for ( int i = 1; i < nbr_cells[1]; ++i )
    buffer_lat[i] = buffer_latv[i-1] + dx;

  yac_cdef_points_reg2d ( grid_id,
                          nbr_cells,
                          CELL,
                          buffer_lon,
                          buffer_lat,
                          &cell_point_id );

  // Define fields

  field_id = xmalloc ( no_of_fields * sizeof(*field_id));

  for ( int i = 0; i < 2; ++i )
    yac_cdef_field ( fieldName[i],
                     comp_id,
                     &cell_point_id,
                     1,1,"3600",YAC_TIME_UNIT_SECOND,
                     &field_id[i] );

  // Invoke the search

  yac_cenddef ( );

 // Data exchange

  recv_buffer = xmalloc ( total_nbr_cells * sizeof (*recv_buffer) );
  send_buffer = xmalloc ( total_nbr_cells * sizeof (*send_buffer) );

  for ( int j = 0; j < nbr_cells[1]; ++j ) {
     for ( int i = 0; i < nbr_cells[0]; ++i ) {
    	 send_buffer[j*nbr_cells[0]+i] = (float) j;
    }
  }

  for ( int i = 0; i < total_nbr_cells; ++i ) recv_buffer[i] = 20;

/*
  // Send fields from compB to compA
  // ------------------------------------

  {
    double *point_set_data[1];
    double **collection_data[1];

    point_set_data[0] = send_buffer;
    collection_data[0] = &(point_set_data[0]);

    // BtoA

    yac_cput ( field_id[1], 1, collection_data, &info, &ierror );

  }
*/

  // Receive fields from compA
  // ------------------------------

  {
    double *collection_data[1];

    collection_data[0] = recv_buffer;

    // AtoB

    yac_cget ( field_id[0], 1, collection_data, &info, &ierror );

  }

  //----------------------------------------------------------
  // write field to vtk output file
  //----------------------------------------------------------

  double * point_data = xmalloc(nbr_vertices[0] * nbr_vertices[1] * 3 * sizeof(*point_data));
  for (int i = 0; i < nbr_vertices[1]; ++i)
    for (int j = 0; j < nbr_vertices[0]; ++j)
      LLtoXYZ(buffer_lonv[j], buffer_latv[i],
              &point_data[3*(i * nbr_vertices[0] + j)]);

  unsigned * cell_corners = xmalloc(total_nbr_cells * 4 * sizeof(*cell_corners));
  unsigned * num_points_per_cell = xmalloc(total_nbr_cells * sizeof(*num_points_per_cell));

  for (int i = 0; i < total_nbr_cells; ++i) {

    unsigned x_index, y_index;

    y_index = i / nbr_cells[0];
    x_index = i - y_index * nbr_cells[0];

    cell_corners[i*4+0] =  y_index      * (nbr_cells[0] + 1) + x_index;
    cell_corners[i*4+1] =  y_index      * (nbr_cells[0] + 1) + x_index + 1;
    cell_corners[i*4+2] = (y_index + 1) * (nbr_cells[0] + 1) + x_index + 1;
    cell_corners[i*4+3] = (y_index + 1) * (nbr_cells[0] + 1) + x_index;

    num_points_per_cell[i] = 4;
  }

  char vtk_filename[32];

  sprintf(vtk_filename, "compB_out_%d.vtk", rank);

  VTK_FILE *vtk_file = vtk_open(vtk_filename, "compB_out");

  vtk_write_point_data(vtk_file, point_data, nbr_vertices[0]*nbr_vertices[1]);

  vtk_write_cell_data(vtk_file, cell_corners, num_points_per_cell, total_nbr_cells);

  vtk_write_cell_scalars_int(vtk_file, global_index, total_nbr_cells,
                             "global_cell_id");

  vtk_write_cell_scalars_double(vtk_file, recv_buffer, total_nbr_cells, "cell_in");

  vtk_close(vtk_file);

  /* Free memory */

  free( global_index );
  free( cell_core_mask );

  free(num_points_per_cell);
  free(cell_corners);
  free(point_data);

  free ( buffer );

  free ( buffer_lon );
  free ( buffer_lat );

  free ( buffer_lonv );
  free ( buffer_latv );

  free (send_buffer);
  free (recv_buffer);
  free (field_id);

  yac_cfinalize ();

}

/* -------------------------------------------------------------------- */



static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, char const ** gridPath) {

  int opt;
  while ((opt = getopt(argc, argv, "c:g:")) != -1) {
    YAC_ASSERT_ARGS((opt == 'c') || (opt == 'g'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
      case 'g':
        *gridPath = optarg;
        break;
    }
  }
}
