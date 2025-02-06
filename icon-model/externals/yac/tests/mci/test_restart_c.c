// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#define EXACT

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "yac.h"
#include "utils_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "generate_cubed_sphere.h"
#include "test_function.h"
#include "tests.h"

#define DUMMY_VALUE (-1337.0)

static int point_id = INT_MAX;
double * field_out_data = NULL;
double * field_in_data = NULL;
size_t field_data_size = 0;
char const * config_filename = "test_restart_c.yaml";

static void generate_icon_grid(
  int comm_rank, int comm_size, char const * grid_dir);
static void generate_cube_grid(int comm_rank, int comm_size);

static void run_model_config(
  int icon_to_cube, int cube_to_icon, int config_from_file,
  char const * grid_dir) {

  yac_cinit ();
  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);

  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime("+1800-01-01T00:00:00.000", "+2100-01-01T00:00:00.000");

  if (global_size < 3) {
    fprintf(stderr, "Wrong number of processes (should be at least 3)\n");
    exit(EXIT_FAILURE);
  }

  int is_icon = global_rank < (global_size / 2);

  // register component(s)
  int comp_id;
  yac_cdef_comp((is_icon)?"icon":"cube", &comp_id);

  MPI_Comm local_comm;
  int local_rank, local_size;
  yac_cget_comp_comm(comp_id, &local_comm);
  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_size(local_comm, &local_size);
  MPI_Comm_free(&local_comm);

  // if no points have been registered yet
  if (point_id == INT_MAX) {
    if (is_icon) generate_icon_grid(local_rank, local_size, grid_dir);
    else         generate_cube_grid(local_rank, local_size);
  }

  // register grid
  int field_out_id, field_in_id;
  int receives_data;
  if (is_icon) {
    receives_data = cube_to_icon;
    yac_cdef_field(
      "icon_to_cube", comp_id, &point_id, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_out_id);
    yac_cdef_field(
      "cube_to_icon", comp_id, &point_id, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_in_id);
  } else {
    receives_data = icon_to_cube;
    yac_cdef_field(
      "cube_to_icon", comp_id, &point_id, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_out_id);
    yac_cdef_field(
      "icon_to_cube", comp_id, &point_id, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_in_id);
  }

  // end defintion phase
  yac_csync_def();

  if (config_from_file) {

    yac_cread_config_yaml(config_filename);
    yac_cenddef( );

  } else {
    // generate an interpolation interpolation stack
    int interp_stack_id;
    yac_cget_interp_stack_config(&interp_stack_id);
    yac_cadd_interp_stack_config_conservative(interp_stack_id, 1, 0, 0, 0);

    // generate couplings
    if (icon_to_cube && !is_icon)
      yac_cdef_couple(
        "icon", "icon", "icon_to_cube",
        "cube", "cube", "icon_to_cube",
        "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
        interp_stack_id, 0, 0);
    if (cube_to_icon && is_icon)
      yac_cdef_couple(
        "cube", "cube", "cube_to_icon",
        "icon", "icon", "cube_to_icon",
        "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
        interp_stack_id, 0, 0);

    yac_cfree_interp_stack_config(interp_stack_id);

    char * config;
    yac_cenddef_and_emit_config(YAC_YAML_EMITTER_DEFAULT, &config);

    if (global_rank == 0) {
      FILE * config_file = fopen(config_filename, "w");
      fputs(config, config_file);
      fclose(config_file);
    }

    free(config);
  }

  // do some ping-pongs
  for (int t = 0; t < 10; ++t) {

    {
      int info, err;
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};
      point_set_data[0] = field_out_data;
      yac_cput(field_out_id, 1, collection_data, &info, &err);
    }

    {
      for (size_t i = 0; i < field_data_size; ++i)
        field_in_data[i] = DUMMY_VALUE;

      int info, err;
      double *collection_data[1] = {field_in_data};
      yac_cget(field_in_id, 1, collection_data, &info, &err);

      for (size_t i = 0; i < field_data_size; ++i) {
        if (receives_data) {
          if (fabs(field_in_data[i] - field_out_data[i]) > 1e-3)
            PUT_ERR("data missmatch");
        } else {
          if (field_in_data[i] != DUMMY_VALUE)
            PUT_ERR("data missmatch");
        }
      }
    }
  }

  yac_ccleanup();
}

int main(int argc, char** argv) {

  if (argc != 2) {
    PUT_ERR("ERROR: missing grid file directory");
    return TEST_EXIT_CODE;
  }

  // run the models multiple time to check whether yac can be restarted
  for (int icon_to_cube = 0; icon_to_cube < 2; ++icon_to_cube)
    for (int cube_to_icon = 0; cube_to_icon < 2; ++cube_to_icon)
      for (int config_from_file = 0; config_from_file < 2; ++config_from_file)
        run_model_config(
          icon_to_cube, cube_to_icon, config_from_file, argv[1]);

  // clean-up
  free(field_out_data);
  free(field_in_data);

  int global_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  if (global_rank == 0) unlink(config_filename);

  yac_cfinalize();

  return TEST_EXIT_CODE;
}

static void generate_icon_grid(
  int comm_rank, int comm_size, char const * grid_dir) {

  int nbr_vertices;
  int nbr_cells;
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

  char const * grid_filename =
    strcat(
      strcpy(malloc(strlen(grid_dir) + 32), grid_dir),
      "icon_grid_0030_R02B03_G.nc");

  yac_read_part_icon_grid_information(
    grid_filename, &nbr_vertices, &nbr_cells, &num_vertices_per_cell,
    &cell_to_vertex, &x_vertices, &y_vertices, &x_cells, &y_cells,
    &global_cell_id, &cell_mask,
    &cell_core_mask, &global_corner_id, &corner_core_mask,
    comm_rank, comm_size);

  int grid_id;
  yac_cdef_grid_unstruct(
    "icon", nbr_vertices, nbr_cells, num_vertices_per_cell,
    x_vertices, y_vertices, cell_to_vertex, &grid_id);

  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);
  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);

  yac_cdef_points_unstruct(
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, &point_id);

  for (int i = 0; i < nbr_cells; ++i) cell_mask[i] = cell_mask[i] == 0;

  yac_cset_mask(cell_mask, point_id);

  field_out_data = xmalloc(nbr_cells * sizeof(*field_out_data));
  field_data_size = nbr_cells;
  for (int i = 0; i < nbr_cells; ++i)
    field_out_data[i] =
      (cell_core_mask[i])?
        (yac_test_harmonic(x_cells[i]*YAC_RAD, y_cells[i]*YAC_RAD)):
        DUMMY_VALUE;
  field_in_data = xmalloc(nbr_cells * sizeof(*field_in_data));

  yac_delete_icon_grid_data(&cell_mask,
                            &global_cell_id,
                            &cell_core_mask,
                            &num_vertices_per_cell,
                            &global_corner_id,
                            &corner_core_mask,
                            &cell_to_vertex,
                            &x_cells,
                            &y_cells,
                            &x_vertices,
                            &y_vertices);
}

static void generate_cube_grid(int comm_rank, int comm_size) {

  unsigned n = 50;

  unsigned nbr_vertices;
  unsigned nbr_cells;
  unsigned * num_vertices_per_cell;
  unsigned * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * cell_core_mask;
  int * corner_core_mask;
  int * global_cell_id;
  int * global_corner_id;

  yac_generate_part_cube_grid_information(n, &nbr_vertices, &nbr_cells,
                                          &num_vertices_per_cell, &cell_to_vertex,
                                          &x_vertices, &y_vertices, &x_cells,
                                          &y_cells, &global_cell_id,
                                          &cell_core_mask, &global_corner_id,
                                          &corner_core_mask, comm_rank, comm_size);

  int grid_id;
  yac_cdef_grid_unstruct(
    "cube", nbr_vertices, nbr_cells, (int*)num_vertices_per_cell,
    x_vertices, y_vertices, (int*)cell_to_vertex, &grid_id);

  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);
  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);

  yac_cdef_points_unstruct(
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, &point_id);

  field_out_data = xmalloc(nbr_cells * sizeof(*field_out_data));
  for (unsigned i = 0; i < nbr_cells; ++i)
    field_out_data[i] =
      (cell_core_mask[i])?
        (yac_test_harmonic(x_cells[i]*YAC_RAD, y_cells[i]*YAC_RAD)):
        DUMMY_VALUE;
  field_in_data = xmalloc(nbr_cells * sizeof(*field_in_data));
  field_data_size = nbr_cells;
  free(cell_core_mask);
  free(corner_core_mask);
  free(global_cell_id);
  free(global_corner_id);
  free(x_vertices);
  free(y_vertices);
  free(x_cells);
  free(y_cells);
  free(num_vertices_per_cell);
  free(cell_to_vertex);
}
