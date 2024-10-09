// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "instance.h"
#include "yac.h"
#include "yac_mpi.h"
#include "dist_grid_utils.h"
#include "weight_file_common.h"
#include "grid_file_common.h"
#include "config_yaml.h"
#include "interp_method_callback.h"
#include "event.h"
#include "geometry.h"

static void compute_weights_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);

int main(int argc, char** argv) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size != 3) {

    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }

  if (argc != 2) {
    PUT_ERR("ERROR: missing config file directory");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  int collection_size = 1;
  char const * timestep_iso8601 = yac_time_to_ISO("1", C_SECOND);

  { // tests with instance_test_1_1.yaml

    { // no process defines a components
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      struct yac_couple_config * couple_config =
        yac_couple_config_new();
      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_1.yaml");
      yac_yaml_read_coupling(
        couple_config, yaml_filename, YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);
      yac_instance_set_couple_config(instance, couple_config);

      yac_instance_def_components(instance, NULL, 0);

      yac_instance_setup(instance, NULL, 0);

      yac_instance_delete(instance);
    }

    { // first two process define component comp_1 and last defines comp_2
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);
      yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

      char const * component_names[2] = {"comp_1", "comp_2"};
      yac_instance_def_components(
        instance, &(component_names[rank >> 1]), 1);

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[3][2] = {{0,1},{1,2},{0.5,1.5}};
      double coordinates_y[3][2] = {{0,1},{0,1},{0,1}};
      yac_int global_cell_ids [3][1] = {{0},{1},{0}};
      yac_int global_corner_ids[3][4] = {{0,1,3,4},{1,2,4,5},{0,1,2,3}};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x[rank], coordinates_y[rank]);
      grid_data.cell_ids = TO_POINTER(global_cell_ids[rank]);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids[rank]);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[2] = {"grid1", "grid2"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank>>1], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      yac_instance_add_field(
        instance, "field_1", component_names[rank >> 1], grid,
        interp_fields, 1, collection_size, timestep_iso8601);

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_1.yaml");
      char const * grid_filename = "instance_test_1_1_grid2.nc";
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      yac_couple_config_grid_set_output_filename(
        yac_instance_get_couple_config(instance), "grid2", grid_filename);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);

      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == 0) {

        unlink("comp_1.err");
        unlink("comp_2.err");
        unlink("comp_1.log");
        unlink("comp_2.log");

        check_grid_file(
          grid_filename, "grid2", 1, 4,
          (double[]){0,0,1,1}, (double[]){0.5,1.5,1.5,0.5}, NULL, NULL,
          (yac_int[]){0}, (int[]){1}, NULL, NULL, NULL, NULL);
        unlink(grid_filename);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  { // tests with instance_test_1_2.yaml

    { // each process has its own component and data exchange occurs in a
      // round robin fashion
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      char const * component_names[3] = {"comp_1", "comp_2", "comp_3"};
      yac_instance_def_components(
        instance, &(component_names[rank]), 1);

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[2] = {0,1};
      double coordinates_y[2] = {0,1};
      yac_int global_cell_ids [1] = {0};
      yac_int global_corner_ids[4] = {0,1,2,3};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x, coordinates_y);
      grid_data.cell_ids = TO_POINTER(global_cell_ids);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[3] = {"grid1", "grid2", "grid3"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      char * field_name[3] = {"field_1", "field_2", "field_3"};
      yac_instance_add_field(
        instance, field_name[rank], component_names[0], grid,
        interp_fields, 1, collection_size, timestep_iso8601);
      yac_instance_add_field(
        instance, field_name[(rank + 2) % 3], component_names[rank], grid,
        interp_fields, 1, collection_size, timestep_iso8601);

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_2.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_3.yaml

    { // each process has its own component and data exchange occurs in a
      // round robin fashion
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      char const * component_names[3] = {"comp_1", "comp_2", "comp_3"};
      yac_instance_def_components(
        instance, &(component_names[rank]), 1);

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[2] = {0,1};
      double coordinates_y[2] = {0,1};
      yac_int global_cell_ids [1] = {0};
      yac_int global_corner_ids[4] = {0,1,2,3};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x, coordinates_y);
      grid_data.cell_ids = TO_POINTER(global_cell_ids);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[3] = {"grid1", "grid2", "grid3"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      char * field_name[6] = {"field_1", "field_2", "field_3",
                              "field_4", "field_5", "field_6"};
      yac_instance_add_field(
        instance, field_name[rank], component_names[rank], grid,
        interp_fields, 1, collection_size, timestep_iso8601);
      yac_instance_add_field(
        instance, field_name[(rank + 2) % 3], component_names[rank], grid,
        interp_fields, 1, collection_size, timestep_iso8601);
      yac_instance_add_field(
        instance, field_name[rank + 3], component_names[rank], grid,
        interp_fields, 1, collection_size, timestep_iso8601);
      yac_instance_add_field(
        instance, field_name[((rank + 2) % 3) + 3], component_names[rank], grid,
        interp_fields, 1, collection_size, timestep_iso8601);

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_3.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_4.yaml

    { // each process has its own component
      // process one has one field that is to be sent to both other
      // processes in the same put
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      char const * component_names[3] = {"comp_1", "comp_2", "comp_3"};
      yac_instance_def_components(
        instance, &(component_names[rank]), 1);

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[2] = {0,1};
      double coordinates_y[2] = {0,1};
      yac_int global_cell_ids [1] = {0};
      yac_int global_corner_ids[4] = {0,1,2,3};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x, coordinates_y);
      grid_data.cell_ids = TO_POINTER(global_cell_ids);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[3] = {"grid1", "grid2", "grid3"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      yac_instance_add_field(
        instance, "field_1", component_names[rank], grid,
        interp_fields, 1, collection_size, timestep_iso8601);

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_4.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_5.yaml

    { // each process has its own component
      // process one has one field that is to be sent to both other
      // processes in the same put
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      char const * component_names[3] = {"comp_1", "comp_2", "comp_3"};
      yac_instance_def_components(
        instance, &(component_names[rank]), 1);

      size_t num_vertices[2] = {5,5};
      int cyclic[2] = {0,0};
      double coordinates_x[5] = {0,1,2,3,4};
      double coordinates_y[5] = {0,1,2,3,4};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x, coordinates_y);

      char * grid_name[3] = {"grid1", "grid2", "grid3"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      yac_instance_add_field(
        instance, "field_1", component_names[rank], grid,
        interp_fields, 1, collection_size, timestep_iso8601);

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_5.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_6.yaml

    { // first two process define component comp_1 and last defines comp_2
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      char const * component_names[2] = {"comp_1", "comp_2"};
      yac_instance_def_components(
        instance, &(component_names[rank >> 1]), 1);

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[3][2] = {{0,1},{1,2},{0.5,1.5}};
      double coordinates_y[3][2] = {{0,1},{0,1},{0,1}};
      yac_int global_cell_ids [3][1] = {{0},{1},{0}};
      yac_int global_corner_ids[3][4] = {{0,1,3,4},{1,2,4,5},{0,1,2,3}};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x[rank], coordinates_y[rank]);
      grid_data.cell_ids = TO_POINTER(global_cell_ids[rank]);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids[rank]);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[2] = {"grid1", "grid2"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank>>1], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      char * weight_file_name = "test_instance_parallel1_weight_file.nc";

      // delete weight file if it extists
      if (rank == 0) unlink(weight_file_name);

      yac_instance_add_field(
        instance, "field_1", component_names[rank >> 1], grid,
        interp_fields, 1, collection_size, timestep_iso8601);

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_6.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      // ensure that the weight file has been written
      MPI_Barrier(MPI_COMM_WORLD);

      // check whether weight file exists
      if (rank == 0) {
        if (access(weight_file_name, F_OK ) == -1)
          PUT_ERR("weight file is missing\n");
        // delete weight file if it extists
        unlink(weight_file_name);
      }

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_7.yaml

    { // first two process define component comp_1 and last defines comp_2
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      char const * component_names[2] = {"comp_1", "comp_2"};
      yac_instance_def_components(
        instance, &(component_names[rank >> 1]), 1);

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[3][2] = {{0,1},{1,2},{0.5,1.5}};
      double coordinates_y[3][2] = {{0,1},{0,1},{0,1}};
      yac_int global_cell_ids [3][1] = {{0},{1},{0}};
      yac_int global_corner_ids[3][4] = {{0,1,3,4},{1,2,4,5},{0,1,2,3}};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x[rank], coordinates_y[rank]);
      grid_data.cell_ids = TO_POINTER(global_cell_ids[rank]);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids[rank]);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[2] = {"grid1", "grid2"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank>>1], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      char * field_name[3] = {"field_1", "field_2", "field_3"};

      char * weight_file_name[3] =
        {NULL, "test_instance_parallel1_weight_file_1.nc",
         "test_instance_parallel1_weight_file_2.nc"};

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_7.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      // delete weight file if it exists
      for (unsigned i = 0; i < 3; ++i) {

        if ((rank == 0) && (weight_file_name[i] != NULL))
          unlink(weight_file_name[i]);

        yac_instance_add_field(
          instance, field_name[i], component_names[rank >> 1], grid,
          interp_fields, 1, collection_size, timestep_iso8601);
      }

      yac_instance_setup(instance, &grid, 1);

      // ensure that the weight file has been written
      MPI_Barrier(MPI_COMM_WORLD);

      // check whether weight file exists
      if (rank == 0) {
        for (unsigned i = 0; i < 3; ++i) {
          if (weight_file_name[i] == NULL) continue;
          if (access(weight_file_name[i], F_OK ) == -1)
            PUT_ERR("weight file is missing\n");
          // delete weight file if it extists
          unlink(weight_file_name[i]);
        }
      }

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_8.yaml

    { // all processes define comp_1, but the yaml lists comp_1 and comp_2
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      char const * component_names[1] = {"comp_1"};
      yac_instance_def_components(instance, component_names, 1);

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[3][2] = {{0,1},{1,2},{2,3}};
      double coordinates_y[3][2] = {{0,1},{0,1},{0,1}};
      yac_int global_cell_ids [3][1] = {{0},{1},{2}};
      yac_int global_corner_ids[3][4] = {{0,1,4,5},{1,2,5,6},{2,3,6,7}};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x[rank], coordinates_y[rank]);
      grid_data.cell_ids = TO_POINTER(global_cell_ids[rank]);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids[rank]);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[1] = {"grid1"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[0], grid_data);

      struct yac_interp_field interp_fields[1];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;

      char * field_name[3] = {"field_1", "field_2", "field_3"};

      for (unsigned i = 0; i < 3; ++i)
        yac_instance_add_field(
          instance, field_name[i], component_names[0], grid,
          interp_fields, 1, collection_size, timestep_iso8601);

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_8.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      yac_instance_delete(instance);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_9.yaml

    { // tests various interpolation methods
      struct yac_instance * instance =
        yac_instance_new(MPI_COMM_WORLD);

      int comp_idx = rank >> 1;
      char const * component_names[2] = {"comp_1", "comp_2"};
      yac_instance_def_components(
        instance, &(component_names[comp_idx]), 1);

      double coordinates_x[5] = {0.0,1.0,2.0,3.0,4.0};
      double coordinates_y[5] = {0.0,1.0,2.0,3.0,4.0};
      size_t const num_cells[2] = {4,4};
      size_t local_start[3][2] = {{0,0},{2,0}, {0,0}};
      size_t local_count[3][2] = {{2,4},{2,4}, {4,4}};
      int with_halo = 0;
      for (size_t i = 0; i < 5; ++i) coordinates_x[i] *= YAC_RAD;
      for (size_t i = 0; i < 5; ++i) coordinates_y[i] *= YAC_RAD;

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[rank], local_count[rank], with_halo);

      char * grid_name[2] = {"grid1", "grid2"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[comp_idx], grid_data);

      yac_coordinate_pointer cell_middle_points =
        xmalloc(grid_data.num_cells * sizeof(*cell_middle_points));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        for (size_t j = 0; j < 3; ++j) cell_middle_points[i][j] = 0.0;
        for (int j = 0; j < grid_data.num_vertices_per_cell[i]; ++j)
          for (size_t k = 0; k < 3; ++k)
            cell_middle_points[i][k] +=
              grid_data.vertex_coordinates[
                grid_data.cell_to_vertex[
                  grid_data.cell_to_vertex_offsets[i] + j]][k];
        normalise_vector(cell_middle_points[i]);
      }

      struct yac_interp_field interp_fields[2];
      interp_fields[0].location = YAC_LOC_CORNER;
      interp_fields[0].coordinates_idx = SIZE_MAX;
      interp_fields[0].masks_idx = SIZE_MAX;
      interp_fields[1].location = YAC_LOC_CELL;
      interp_fields[1].coordinates_idx =
        yac_basic_grid_add_coordinates(
          grid, YAC_LOC_CELL, cell_middle_points, grid_data.num_cells);
      interp_fields[1].masks_idx = SIZE_MAX;
      free(cell_middle_points);

      char const * weight_file_name = "instance_test_1_9.nc";
      if (rank == 0) {

        int src_indices[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        int tgt_indices[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        double weights[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        size_t num_links = 16;
        enum yac_location src_locations[1] = {YAC_LOC_CELL};
        enum yac_location tgt_location = YAC_LOC_CELL;
        unsigned num_src_fields = 1;
        int num_links_per_field[1] = {num_links};
        int * tgt_id_fixed = NULL;
        size_t num_fixed_tgt = 0;
        double * fixed_values = NULL;
        int * num_tgt_per_fixed_value = NULL;
        size_t num_fixed_values = 0;

        write_weight_file(
          weight_file_name, src_indices, tgt_indices, weights, num_links,
          src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
          num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
          num_fixed_values, tgt_location, grid_name[0], grid_name[1]);
      }

      // register weight computation callback routine
      if (comp_idx == 0)
        yac_cadd_compute_weights_callback(
          compute_weights_callback, NULL, "compute_weights_callback");

      char * field_name[] = {"AVG_ARITHMETIC",
                             "AVG_DIST",
                             "AVG_BARY",
                             "4NN_ARITHMETIC",
                             "4NN_DIST",
                             "4NN_GAUSS",
                             "HCSBB",
                             "RBF_4_GAUSS",
                             "FIXED",
                             "SPMAP",
                             "CONSERV_FRACAREA",
                             "CONSERV_DESTAREA",
                             "CONSERV2ND",
                             "USER_FILE",
                             "CREEP",
                             "USER_CALLBACK"};
      size_t const field_count = 16;

      for (size_t i = 0; i < 9; ++i)
        yac_instance_add_field(
          instance, field_name[i], component_names[comp_idx], grid,
          &(interp_fields[0]), 1, collection_size, timestep_iso8601);
      for (size_t i = 9; i < 15; ++i)
        yac_instance_add_field(
          instance, field_name[i], component_names[comp_idx], grid,
          &(interp_fields[1]), 1, collection_size, timestep_iso8601);

      // two source fields and one target field
      if (comp_idx == 0) {
        for (size_t i = 15; i < field_count; ++i)
          yac_instance_add_field(
            instance, field_name[i], component_names[comp_idx], grid,
            &(interp_fields[0]), 2, collection_size, timestep_iso8601);
      } else {
        for (size_t i = 15; i < field_count; ++i)
          yac_instance_add_field(
            instance, field_name[i], component_names[comp_idx], grid,
            &(interp_fields[0]), 1, collection_size, timestep_iso8601);
      }

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_9.yaml");
      yac_yaml_read_coupling(
        yac_instance_get_couple_config(instance), yaml_filename,
        YAC_YAML_PARSER_DEFAULT);
      free(yaml_filename);

      yac_instance_setup(instance, &grid, 1);

      yac_instance_delete(instance);
      if (rank == 0) unlink(weight_file_name);

      yac_basic_grid_delete(grid);
    }
  }

  { // tests with instance_test_1_10.yaml

    { // first two process define component comp_1 and last defines comp_2
      // configuration file contains a coupling, however only source field
      // is defined

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};
      double coordinates_x[3][2] = {{0,1},{1,2},{0.5,1.5}};
      double coordinates_y[3][2] = {{0,1},{0,1},{0,1}};
      yac_int global_cell_ids [3][1] = {{0},{1},{0}};
      yac_int global_corner_ids[3][4] = {{0,1,3,4},{1,2,4,5},{0,1,2,3}};
      int cell_core_mask[1] = {1};

      struct yac_basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x[rank], coordinates_y[rank]);
      grid_data.cell_ids = TO_POINTER(global_cell_ids[rank]);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids[rank]);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);

      char * grid_name[2] = {"grid1", "grid2"};
      struct yac_basic_grid * grid =
        yac_basic_grid_new(grid_name[rank>>1], grid_data);

      // define field only on source/target component
      for (int i = 0; i < 2; ++i) {

        struct yac_instance * instance =
          yac_instance_new(MPI_COMM_WORLD);

        char const * component_names[2] = {"comp_1", "comp_2"};
        yac_instance_def_components(
          instance, &(component_names[rank >> 1]), 1);

        if ((rank >> 1) == i) {

          struct yac_interp_field interp_fields[1];
          interp_fields[0].location = YAC_LOC_CORNER;
          interp_fields[0].coordinates_idx = SIZE_MAX;
          interp_fields[0].masks_idx = SIZE_MAX;

          yac_instance_add_field(
            instance, "field", component_names[rank >> 1], grid,
            interp_fields, 1, collection_size, timestep_iso8601);
        }

      char * yaml_filename =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + 32), argv[1]), "instance_test_1_10.yaml");
        yac_yaml_read_coupling(
          yac_instance_get_couple_config(instance), yaml_filename,
          YAC_YAML_PARSER_DEFAULT);
        free(yaml_filename);

        yac_instance_setup(instance, &grid, 1);

        yac_instance_delete(instance);
      }

      yac_basic_grid_delete(grid);
    }
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void compute_weights_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  UNUSED(tgt_coords);
  UNUSED(src_cell_id);
  UNUSED(src_cell_idx);
  UNUSED(user_data);

  for (size_t i = 0; i < 2; ++i) {
    global_results_points[i] = NULL;
    result_weights[i] = NULL;
    result_count[i] = 0;
  }
}
