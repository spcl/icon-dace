// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#define EXACT

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yac.h"
#include "utils_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "generate_cubed_sphere.h"
#include "test_function.h"
#include "tests.h"

#define DUMMY_VALUE (-1337.0)

struct field_config {
  double * data;
  size_t data_size;
  int id;
  int dummy_id;
};

static void generate_icon_grid(
  int comm_rank, int comm_size, int comp_id, int * grid_id,
  struct field_config * in_config, struct field_config * out_config,
  char const * grid_dir);
static void generate_cube_grid(
  int comm_rank, int comm_size, int comp_id, int * grid_id,
  struct field_config * in_config, struct field_config * out_config);
static void init_in_field_data(struct field_config field_config);
static void check_results(
  struct field_config field_config, struct field_config ref_field_config);

int main(int argc, char** argv) {

  yac_cinit ();

  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);

  if (global_size != 3) {
    fprintf(stderr, "Wrong number of processes (should be 3)\n");
    exit(EXIT_FAILURE);
  }

  if (argc != 3) {
    PUT_ERR("ERROR: missing config/grid file directory");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  // register component(s)
  // (rank 0 => comp1; rank 1 => comp2; rank 2 => comp1 + comp2)
  int comp_ids[2];
  int num_comps;
  if ((global_rank == 0) || (global_rank == 1)) {
    yac_cdef_comp((global_rank == 0)?"comp1":"comp2", comp_ids);
    num_comps = 1;
  } else {
    char const * comp_names[2] = {"comp1", "comp2"};
    yac_cdef_comps(&(comp_names[0]), 2, comp_ids);
    num_comps = 2;
  }

  // register grids
  int grid_ids[2];
  struct field_config in_field_configs[2];
  struct field_config out_field_configs[2];
  if ((global_rank == 0) || (global_rank == 2))
    generate_icon_grid(
      (global_rank == 2), 2, comp_ids[0], &(grid_ids[0]),
      &(in_field_configs[0]),
      &(out_field_configs[0]), argv[2]);

  if ((global_rank == 1) || (global_rank == 2))
    generate_cube_grid(
      (global_rank == 2), 2, comp_ids[global_rank == 2],
      &(grid_ids[global_rank == 2]),
      &(in_field_configs[global_rank == 2]),
      &(out_field_configs[global_rank == 2]));

  // read configuration file
  char * yaml_filename =
    strcat(
      strcpy(
        malloc(strlen(argv[1]) + 32), argv[1]), "coupling_test2.yaml");
  yac_cread_config_yaml(yaml_filename);
  free(yaml_filename);

  yac_cenddef();

  // do some ping-pongs
  for (int t = 0; t < 100; ++t) {

    if (global_rank == 0) {

      {
        int send_info, recv_info, err;
        int id = out_field_configs[0].id;
        int dummy_id = out_field_configs[0].dummy_id;
        double *point_set_data[1];
        double **collection_data[1] = {point_set_data};
        point_set_data[0] = out_field_configs[0].data;

        yac_cexchange(
          id, dummy_id, 1, collection_data, NULL,
          &send_info, &recv_info, &err);

        if (send_info != YAC_ACTION_COUPLING) PUT_ERR("ERROR in yac_cexchange");
        if (recv_info != YAC_ACTION_NONE) PUT_ERR("ERROR in yac_cexchange");
        if (err)  PUT_ERR("ERROR in yac_cexchange");
      }

      {
        int send_info, recv_info, err;
        int id = in_field_configs[0].id;
        int dummy_id = in_field_configs[0].dummy_id;
        double *collection_data[1] = {in_field_configs[0].data};

        init_in_field_data(in_field_configs[0]);

        yac_cexchange(
          dummy_id, id, 1, NULL, collection_data,
          &send_info, &recv_info, &err);

        check_results(in_field_configs[0], out_field_configs[0]);
        if (send_info != YAC_ACTION_NONE) PUT_ERR("ERROR in yac_cexchange");
        if (recv_info != YAC_ACTION_COUPLING) PUT_ERR("ERROR in yac_cexchange");
        if (err)  PUT_ERR("ERROR in yac_cexchange");
      }
    } else if (global_rank == 1) {

      {
        int info, err;
        int id = in_field_configs[0].id;
        double *collection_data[1] = {in_field_configs[0].data};

        init_in_field_data(in_field_configs[0]);

        yac_cget(id, 1, collection_data, &info, &err);

        check_results(in_field_configs[0], out_field_configs[0]);
      }

      {
        int info, err;
        int id = out_field_configs[0].id;
        double *point_set_data[1];
        double **collection_data[1] = {point_set_data};
        point_set_data[0] = out_field_configs[0].data;

        yac_cput(id, 1, collection_data, &info, &err);
      }
    } else {

      for (int i = 0; i < 2; ++i) {
        int out_info, in_info, err;
        int out_id = out_field_configs[i].id;
        int in_id = in_field_configs[i^1].id;
        double *out_point_set_data[1];
        double **out_collection_data[1] = {out_point_set_data};
        out_point_set_data[0] = out_field_configs[i].data;
        double *in_collection_data[1] = {in_field_configs[i^1].data};

        init_in_field_data(in_field_configs[i^1]);

        yac_cexchange(
          out_id, in_id, 1, out_collection_data, in_collection_data,
          &out_info, &in_info, &err);

        check_results(in_field_configs[i^1], out_field_configs[i^1]);
      }
    }
  }

  // clean-up

  for (int i = 0; i < num_comps; ++i) {
    free(out_field_configs[i].data);
    free(in_field_configs[i].data);
  }

  yac_cfinalize();

  return TEST_EXIT_CODE;
}

static void generate_icon_grid(
  int comm_rank, int comm_size, int comp_id, int * grid_id,
  struct field_config * in_config, struct field_config * out_config,
  char const * grid_dir) {

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

  char * grid_filename =
    strcat(
      strcpy(
        malloc(strlen(grid_dir) + 32), grid_dir), "icon_grid_0030_R02B03_G.nc");
  yac_read_part_icon_grid_information(
    grid_filename, &nbr_vertices, &nbr_cells, &num_vertices_per_cell,
    &cell_to_vertex, &x_vertices, &y_vertices, &x_cells, &y_cells,
    &global_cell_id, &cell_mask,
    &cell_core_mask, &global_corner_id, &corner_core_mask,
    comm_rank, comm_size);
  free(grid_filename);

  yac_cdef_grid_unstruct(
    "icon", nbr_vertices, nbr_cells, num_vertices_per_cell,
    x_vertices, y_vertices, cell_to_vertex, grid_id);

  if (yac_cget_grid_size(YAC_LOCATION_CORNER, *grid_id) != (size_t)nbr_vertices)
    PUT_ERR("error in yac_cget_grid_size");
  if (yac_cget_grid_size(YAC_LOCATION_CELL, *grid_id) != (size_t)nbr_cells)
    PUT_ERR("error in yac_cget_grid_size");

  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, *grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, *grid_id);
  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, *grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, *grid_id);

  int corner_point_id, cell_point_id;
  yac_cdef_points_unstruct(
    *grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices,
    &corner_point_id);
  yac_cdef_points_unstruct(
    *grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells,
    &cell_point_id);

  if (yac_cget_points_size(corner_point_id) != (size_t)nbr_vertices)
    PUT_ERR("error in yac_cget_points_size");
  if (yac_cget_points_size(cell_point_id) != (size_t)nbr_cells)
    PUT_ERR("error in yac_cget_points_size");

  for (int i = 0; i < nbr_cells; ++i) cell_mask[i] = cell_mask[i] == 0;

  yac_cset_mask(cell_mask, cell_point_id);

  yac_cdef_field(
    "icon_to_cube", comp_id, &cell_point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND,
    &(out_config->id));
  yac_cdef_field(
    "icon_to_cube_dummy", comp_id, &cell_point_id, 1, 1, "1",
    YAC_TIME_UNIT_SECOND, &(out_config->dummy_id));
  yac_cdef_field(
    "cube_to_icon", comp_id, &cell_point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND,
    &(in_config->id));
  yac_cdef_field(
    "cube_to_icon_dummy", comp_id, &cell_point_id, 1, 1, "1",
    YAC_TIME_UNIT_SECOND, &(in_config->dummy_id));

  out_config->data = xmalloc(nbr_cells * sizeof(*(out_config->data)));
  for (int i = 0; i < nbr_cells; ++i)
    out_config->data[i] =
      (cell_core_mask[i])?
        (yac_test_harmonic(x_cells[i]*YAC_RAD, y_cells[i]*YAC_RAD)):
        DUMMY_VALUE;
  in_config->data = xmalloc(nbr_cells * sizeof(*(in_config->data)));

  in_config->data_size = nbr_cells;
  out_config->data_size = nbr_cells;

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

static void generate_cube_grid(
  int comm_rank, int comm_size, int comp_id, int * grid_id,
  struct field_config * in_config, struct field_config * out_config) {

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

  yac_cdef_grid_unstruct(
    "cube", nbr_vertices, nbr_cells, (int*)num_vertices_per_cell,
    x_vertices, y_vertices, (int*)cell_to_vertex, grid_id);

  if (yac_cget_grid_size(YAC_LOCATION_CORNER, *grid_id) != nbr_vertices)
    PUT_ERR("error in yac_cget_grid_size");
  if (yac_cget_grid_size(YAC_LOCATION_CELL, *grid_id) != nbr_cells)
    PUT_ERR("error in yac_cget_grid_size");

  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, *grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, *grid_id);
  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, *grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, *grid_id);

  int corner_point_id, cell_point_id;
  yac_cdef_points_unstruct(
    *grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices,
    &corner_point_id);
  yac_cdef_points_unstruct(
    *grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells,
    &cell_point_id);

  if (yac_cget_points_size(corner_point_id) != nbr_vertices)
    PUT_ERR("error in yac_cget_points_size");
  if (yac_cget_points_size(cell_point_id) != nbr_cells)
    PUT_ERR("error in yac_cget_points_size");

  yac_cdef_field("cube_to_icon", comp_id, &cell_point_id, 1,
    1, "1", YAC_TIME_UNIT_SECOND, &(out_config->id));
  yac_cdef_field("cube_to_icon_dummy", comp_id, &cell_point_id, 1,
    1, "1", YAC_TIME_UNIT_SECOND, &(out_config->dummy_id));
  yac_cdef_field("icon_to_cube", comp_id, &cell_point_id, 1,
    1, "1", YAC_TIME_UNIT_SECOND, &(in_config->id));
  yac_cdef_field("icon_to_cube_dummy", comp_id, &cell_point_id, 1,
    1, "1", YAC_TIME_UNIT_SECOND, &(in_config->dummy_id));

  out_config->data = xmalloc(nbr_cells * sizeof(*(out_config->data)));
  for (unsigned i = 0; i < nbr_cells; ++i)
    out_config->data[i] =
      (cell_core_mask[i])?
        (yac_test_harmonic(x_cells[i]*YAC_RAD, y_cells[i]*YAC_RAD)):
        DUMMY_VALUE;
  in_config->data = xmalloc(nbr_cells * sizeof(*(in_config->data)));
  out_config->data_size = nbr_cells;
  in_config->data_size = nbr_cells;
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

static void init_in_field_data(struct field_config field_config) {
  for (size_t i = 0; i < field_config.data_size; ++i)
    field_config.data[i] = DUMMY_VALUE;
}

static void check_results(
  struct field_config field_config, struct field_config ref_field_config) {

  if (field_config.data_size != ref_field_config.data_size) {
    fputs("ERROR(check_results): inconsistent data_size\n", stderr);
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < field_config.data_size; ++i) {
    if (fabs(field_config.data[i] - ref_field_config.data[i]) > 1e-3) {
      fputs("ERROR(check_results): data missmatch\n", stderr);
      exit(EXIT_FAILURE);
    }
  }
}
