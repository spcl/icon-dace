// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "yac.h"
#include "yac_utils.h"

/* ------------------------------------------------- */

const char fieldName[] = "icon_to_cube";

#define STR_USAGE "Usage: %s -c configFilename -n cube edge length\n"
#define YAC_ASSERT_ARGS(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, char const ** configFilename, size_t * cube_n);

int main (int argc, char *argv[]) {

  // Initialisation of MPI
  MPI_Init (0, NULL);
  xt_initialize(MPI_COMM_WORLD);

  char const * configFilename = "toy_callback.yaml"; // default configuration file
  size_t cube_n = 202; // default cube edge length
  parse_arguments(argc, argv, &configFilename, &cube_n);
  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  int comp_id;
  char * comp_name = "CUBE";

  yac_cdef_comp ( comp_name, &comp_id );

  MPI_Comm local_comm;
  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;
  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);

  MPI_Comm_free(&local_comm);


  int corner_point_id;

  int field_id;
  int grid_id;

  unsigned nbr_vertices;
  unsigned nbr_cells;
  unsigned * num_vertices_per_cell;
  unsigned * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * global_cell_id;
  int * cell_core_mask;
  int * global_corner_id;
  int * corner_core_mask;

  yac_generate_part_cube_grid_information((unsigned)cube_n, &nbr_vertices, &nbr_cells,
                                          &num_vertices_per_cell, &cell_to_vertex,
                                          &x_vertices, &y_vertices, &x_cells,
                                          &y_cells, &global_cell_id,
                                          &cell_core_mask, &global_corner_id,
                                          &corner_core_mask, rank, size);

  double * x_points, * y_points;

  x_points = x_vertices;
  y_points = y_vertices;

  yac_cdef_grid_unstruct(
    "cube_grid", nbr_vertices, nbr_cells, (int*)num_vertices_per_cell,
    x_vertices, y_vertices, (int*)cell_to_vertex, &grid_id);


  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);

  yac_cdef_points_unstruct(
    grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_points, y_points, &corner_point_id);

  /* Field definition for the component */
  yac_cdef_field(
    fieldName, comp_id, &corner_point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_id);

  /* Search. */
  yac_cenddef( );

  double * callback_in = malloc(nbr_vertices * sizeof(*callback_in));
  for (unsigned i = 0; i < nbr_vertices; ++i) callback_in[i] = -10;

  int err;
  int info;
  {
    double *collection_data[1] = {callback_in};

    yac_cget(field_id, 1, collection_data, &info, &err);
  }

#ifdef VTK_OUTPUT
  //----------------------------------------------------------
  // write field to vtk output file
  //----------------------------------------------------------

  char vtk_filename[32];

  sprintf(vtk_filename, "toy_cube_callback_%d.vtk", rank);

  YAC_VTK_FILE *vtk_file = yac_vtk_open(vtk_filename, "cube_out");
  yac_vtk_write_point_data_ll(
    vtk_file, x_vertices, y_vertices, nbr_vertices);
  yac_vtk_write_cell_data(vtk_file, (unsigned *)cell_to_vertex,
                      (unsigned*)num_vertices_per_cell, nbr_cells);
  yac_vtk_write_point_scalars_int(
    vtk_file, corner_core_mask, nbr_vertices, "corner_core_mask");
  yac_vtk_write_point_scalars_int(
    vtk_file, global_corner_id, nbr_vertices, "global_corner_id");
  yac_vtk_write_cell_scalars_int(
    vtk_file, global_cell_id, nbr_cells, "global_cell_id");

  yac_vtk_write_point_scalars_double(
    vtk_file, callback_in, nbr_vertices, "callback_in");

  yac_vtk_close(vtk_file);
#endif // VTK_OUTPUT

  free(corner_core_mask);
  free(global_corner_id);
  free(cell_core_mask);
  free(global_cell_id);
  free(x_cells);
  free(y_cells);
  free(x_vertices);
  free(y_vertices);
  free(num_vertices_per_cell);
  free(cell_to_vertex);

  yac_cfinalize();

  xt_finalize();

  MPI_Finalize();

  free(callback_in);

  return EXIT_SUCCESS;
}


static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, size_t * cube_n) {

  int opt;
  while ((opt = getopt(argc, argv, "c:n:")) != -1) {
    YAC_ASSERT((opt == 'c') || (opt == 'n'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
      case 'n':
        *cube_n = atoi(optarg);
        YAC_ASSERT_ARGS(*cube_n > 0, "invalid cube edge length");
        break;
    }
  }
}
