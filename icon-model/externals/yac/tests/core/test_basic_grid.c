// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "basic_grid.h"
#include "geometry.h"

static enum yac_location locations[] =
  {YAC_LOC_CELL, YAC_LOC_CORNER, YAC_LOC_EDGE};
enum {
  NUM_LOCACTIONS = sizeof(locations) / sizeof(locations[0])
};

static void check_basic_grid(
  struct yac_basic_grid * grid, char const * ref_grid_name,
  struct yac_basic_grid_data ref_grid_data,
  int const * const * ref_core_masks,
  size_t const * ref_grid_data_sizes);

static void check_interp_field(
  struct yac_basic_grid * grid, struct yac_interp_field field,
  yac_const_coordinate_pointer ref_coords, int const * ref_mask);

// declaration of Fortran to C interfaces
size_t yac_basic_grid_get_data_size_f2c(
  struct yac_basic_grid * grid, int location);
size_t yac_basic_grid_add_coordinates_nocpy_f2c(
  struct yac_basic_grid * grid, int location, double * coordinates);
size_t yac_basic_grid_add_coordinates_f2c(
  struct yac_basic_grid * grid, int location,
  double * coordinates, size_t count);
size_t yac_basic_grid_add_mask_nocpy_f2c(
  struct yac_basic_grid * grid, int location,
  int const * mask, char const * mask_name);
size_t yac_basic_grid_add_mask_f2c(
  struct yac_basic_grid * grid, int location,
  int const * mask, size_t count, char const * mask_name);

int main (void) {

  char const * grid_name = "basic_grid_name";

  { // test yac_basic_grid_new

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0,1};
    double lat_vertices[] = {0,1,2};

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid * grid = yac_basic_grid_new(grid_name, grid_data);

    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_empty_new

    struct yac_basic_grid * grid = yac_basic_grid_empty_new(grid_name);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_unstruct(0, 0, NULL, NULL, NULL, NULL);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {0, 0, 0};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_reg_2d_new

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0.0,0.1};
    double lat_vertices[] = {0.0,0.1,0.2};

    struct yac_basic_grid * grid =
      yac_basic_grid_reg_2d_new(
        grid_name, nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_reg_2d(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_reg_2d_deg_new

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0.0,0.1};
    double lat_vertices[] = {0.0,0.1,0.2};

    struct yac_basic_grid * grid =
      yac_basic_grid_reg_2d_deg_new(
        grid_name, nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_curve_2d_new

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double lat_vertices[] = {0.0,0.0,0.1,0.1,0.2,0.2};

    struct yac_basic_grid * grid =
      yac_basic_grid_curve_2d_new(
        grid_name, nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_curve_2d(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_curve_2d_deg_new

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double lat_vertices[] = {0.0,0.0,0.1,0.1,0.2,0.2};

    struct yac_basic_grid * grid =
      yac_basic_grid_curve_2d_deg_new(
        grid_name, nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_curve_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_unstruct_new

    size_t nbr_vertices = 6;
    size_t nbr_cells = 2;
    int num_vertices_per_cell[] = {4, 4};
    double x_vertices[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double y_vertices[] = {0.0,0.0,0.1,0.1,0.2,0.2};
    int cell_to_vertex[] = {0,1,3,2, 2,3,5,4};

    struct yac_basic_grid * grid =
      yac_basic_grid_unstruct_new(
        grid_name, nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_unstruct(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_unstruct_deg_new

    size_t nbr_vertices = 6;
    size_t nbr_cells = 2;
    int num_vertices_per_cell[] = {4, 4};
    double x_vertices[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double y_vertices[] = {0.0,0.0,0.1,0.1,0.2,0.2};
    int cell_to_vertex[] = {0,1,3,2, 2,3,5,4};

    struct yac_basic_grid * grid =
      yac_basic_grid_unstruct_deg_new(
        grid_name, nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_unstruct_deg(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_unstruct_ll_new

    size_t nbr_vertices = 6;
    size_t nbr_cells = 2;
    int num_vertices_per_cell[] = {4, 4};
    double x_vertices[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double y_vertices[] = {0.0,0.0,0.1,0.1,0.2,0.2};
    int cell_to_vertex[] = {0,1,3,2, 2,3,5,4};

    struct yac_basic_grid * grid =
      yac_basic_grid_unstruct_ll_new(
        grid_name, nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_unstruct_ll(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_unstruct_ll_deg_new

    size_t nbr_vertices = 6;
    size_t nbr_cells = 2;
    int num_vertices_per_cell[] = {4, 4};
    double x_vertices[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double y_vertices[] = {0.0,0.0,0.1,0.1,0.2,0.2};
    int cell_to_vertex[] = {0,1,3,2, 2,3,5,4};

    struct yac_basic_grid * grid =
      yac_basic_grid_unstruct_ll_deg_new(
        grid_name, nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_unstruct_ll_deg(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {2, 6, 7};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_cloud_new

    size_t nbr_points = 6;
    double x_points[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double y_points[] = {0.0,0.0,0.1,0.1,0.2,0.2};

    struct yac_basic_grid * grid =
      yac_basic_grid_cloud_new(
        grid_name, nbr_points, x_points, y_points);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_cloud(nbr_points, x_points, y_points);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {0, 6, 0};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_basic_grid_cloud_deg_new

    size_t nbr_points = 6;
    double x_points[] = {0.0,0.1,0.0,0.1,0.0,0.1};
    double y_points[] = {0.0,0.0,0.1,0.1,0.2,0.2};

    struct yac_basic_grid * grid =
      yac_basic_grid_cloud_deg_new(
        grid_name, nbr_points, x_points, y_points);

    struct yac_basic_grid_data ref_grid_data =
      yac_generate_basic_grid_data_cloud_deg(nbr_points, x_points, y_points);
    int const * ref_core_masks[NUM_LOCACTIONS] = {NULL, NULL, NULL};
    size_t const ref_grid_data_sizes[NUM_LOCACTIONS] = {0, 6, 0};

    // some basic checks
    check_basic_grid(
      grid, grid_name, ref_grid_data, ref_core_masks, ref_grid_data_sizes);

    // cleanup
    yac_basic_grid_data_free(ref_grid_data);
    yac_basic_grid_delete(grid);
  }

  { // test yac_interp_fields

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0,1};
    double lat_vertices[] = {0,1,2};

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid * grid = yac_basic_grid_new(grid_name, grid_data);

    {
      double ref_cell_center_coordinates[2][3];
      int ref_cell_mask[] = {0,1};
      LLtoXYZ_deg(0.5, 0.5, ref_cell_center_coordinates[0]);
      LLtoXYZ_deg(0.5, 1.5, ref_cell_center_coordinates[1]);
      struct yac_interp_field cell_interp_field =
        {
          .location = YAC_LOC_CELL,
          .coordinates_idx =
            yac_basic_grid_add_coordinates(
              grid, YAC_LOC_CELL, ref_cell_center_coordinates, 2),
          .masks_idx =
            yac_basic_grid_add_mask(
              grid, YAC_LOC_CELL, ref_cell_mask, 2, NULL)
        };

      check_interp_field(
        grid, cell_interp_field,
        (yac_const_coordinate_pointer)ref_cell_center_coordinates, ref_cell_mask);
    }

    {
      struct yac_interp_field cell_interp_field =
        {
          .location = YAC_LOC_CELL,
          .coordinates_idx = SIZE_MAX,
          .masks_idx = SIZE_MAX
        };

      check_interp_field(
        grid, cell_interp_field, NULL, NULL);
    }

    // cleanup
    yac_basic_grid_delete(grid);
  }

  { // test named masks

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0,1};
    double lat_vertices[] = {0,1,2};

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid * grid = yac_basic_grid_new(grid_name, grid_data);

    int const cell_mask[] = {1,0,1,0,1,0};

    char const mask_name[] = "test_mask";

    size_t ref_mask_idx =
      yac_basic_grid_add_mask(
          grid, YAC_LOC_CORNER, cell_mask, 6, mask_name);

    if (ref_mask_idx !=
        yac_basic_grid_get_named_mask_idx(
          grid, YAC_LOC_CORNER, mask_name))
      PUT_ERR("ERROR in yac_basic_grid_get_named_mask_idx");

    if (SIZE_MAX !=
        yac_basic_grid_get_named_mask_idx(
          grid, YAC_LOC_CORNER, NULL))
      PUT_ERR("ERROR in yac_basic_grid_get_named_mask_idx");

    // cleanup
    yac_basic_grid_delete(grid);
  }

  { // test field data

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0,1};
    double lat_vertices[] = {0,1,2};

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid * grid = yac_basic_grid_new(grid_name, grid_data);

    {
      double cell_coordinates[2][2][3];
      double vertex_coordinates[1][6][3];
      int cell_mask[3][2] = {{1,0},{0,1},{0,0}};
      int edge_mask[2][7] = {{0,0,0,1,1,1,1},{1,1,1,0,0,0,0}};

      LLtoXYZ_deg(0.5, 0.5, cell_coordinates[0][0]);
      LLtoXYZ_deg(0.5, 1.5, cell_coordinates[0][1]);
      LLtoXYZ_deg(0.6, 0.6, cell_coordinates[1][0]);
      LLtoXYZ_deg(0.6, 1.6, cell_coordinates[1][1]);

      LLtoXYZ_deg(0, 0, vertex_coordinates[0][0]);
      LLtoXYZ_deg(1, 0, vertex_coordinates[0][1]);
      LLtoXYZ_deg(0, 1, vertex_coordinates[0][2]);
      LLtoXYZ_deg(1, 1, vertex_coordinates[0][3]);
      LLtoXYZ_deg(0, 2, vertex_coordinates[0][4]);
      LLtoXYZ_deg(1, 2, vertex_coordinates[0][5]);

      enum {
        NUM_CELL_COORDS =
          sizeof(cell_coordinates)/sizeof(cell_coordinates[0]),
        NUM_CELL_MASKS = sizeof(cell_mask)/sizeof(cell_mask[0]),
        NUM_CORNER_COORDS =
          sizeof(vertex_coordinates)/sizeof(vertex_coordinates[0]),
        NUM_CORNER_MASKS = 0,
        NUM_EDGE_COORDS = 0,
        NUM_EDGE_MASKS = sizeof(edge_mask)/sizeof(edge_mask[0]),
      };

      size_t ref_coordinates_idxs[3][2] = {
        {yac_basic_grid_add_coordinates(
           grid, YAC_LOC_CELL, cell_coordinates[0], 2),
         yac_basic_grid_add_coordinates(
           grid, YAC_LOC_CELL, cell_coordinates[1], 2)},
        {yac_basic_grid_add_coordinates(
           grid, YAC_LOC_CORNER, vertex_coordinates[0], 6)},
        {SIZE_MAX}
      };
      size_t ref_masks_idxs[3][3] = {
        {yac_basic_grid_add_mask(
           grid, YAC_LOC_CELL, cell_mask[0], 2, NULL),
         yac_basic_grid_add_mask(
           grid, YAC_LOC_CELL, cell_mask[1], 2, NULL),
         yac_basic_grid_add_mask(
           grid, YAC_LOC_CELL, cell_mask[2], 2, NULL)},
        {SIZE_MAX},
        {yac_basic_grid_add_mask(
           grid, YAC_LOC_EDGE, edge_mask[0], 7, NULL),
         yac_basic_grid_add_mask(
           grid, YAC_LOC_EDGE, edge_mask[1], 7, NULL)}
      };

      size_t ref_num_coords[3] =
        {NUM_CELL_COORDS, NUM_CORNER_COORDS, NUM_EDGE_COORDS};
      size_t ref_num_masks[3] =
        {NUM_CELL_MASKS, NUM_CORNER_MASKS, NUM_EDGE_MASKS};

      for (size_t loc_idx = 0; loc_idx < NUM_LOCACTIONS; ++loc_idx) {

        struct yac_field_data * field_data =
          yac_basic_grid_get_field_data(grid, locations[loc_idx]);

        if (ref_num_coords[loc_idx] !=
            yac_field_data_get_coordinates_count(field_data))
          PUT_ERR("ERROR in yac_field_data_get_coordinates_count");

        size_t data_size = yac_basic_grid_get_data_size(grid, locations[loc_idx]);

        for (size_t coord_idx = 0; coord_idx < ref_num_coords[loc_idx];
             ++coord_idx) {

          struct yac_interp_field interp_field =
            {.location = locations[loc_idx],
             .coordinates_idx = ref_coordinates_idxs[loc_idx][coord_idx],
             .masks_idx = SIZE_MAX};

          yac_const_coordinate_pointer ref_coords =
            yac_basic_grid_get_field_coordinates(grid, interp_field);
          yac_const_coordinate_pointer coords =
            yac_field_data_get_coordinates_data(field_data, coord_idx);

          for (size_t i = 0; i < data_size; ++i)
            if (get_vector_angle(ref_coords[i], coords[i]) > yac_angle_tol)
              PUT_ERR("ERROR in yac_basic_grid_get_field_data");
        }

        for (size_t masks_idx = 0; masks_idx < ref_num_masks[loc_idx];
             ++masks_idx) {

          struct yac_interp_field interp_field =
            {.location = locations[loc_idx],
             .coordinates_idx = SIZE_MAX,
             .masks_idx = ref_masks_idxs[loc_idx][masks_idx]};

          int const * ref_masks =
            yac_basic_grid_get_field_mask(grid, interp_field);
          int const * masks =
            yac_field_data_get_mask_data(field_data, masks_idx);

          for (size_t i = 0; i < data_size; ++i)
            if (ref_masks[i] != masks[i])
              PUT_ERR("ERROR in yac_basic_grid_get_field_data");
        }
      }
    }

    // cleanup
    yac_basic_grid_delete(grid);
  }

  { // test Fortran to C interfaces

    size_t nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double lon_vertices[] = {0,1};
    double lat_vertices[] = {0,1,2};

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);

    struct yac_basic_grid * grid = yac_basic_grid_new(grid_name, grid_data);

    {
      double ref_cell_center_coordinates[2][2][3];
      int ref_cell_mask[2][2] = {{0,1},{1,0}};
      LLtoXYZ_deg(0.5, 0.5, ref_cell_center_coordinates[0][0]);
      LLtoXYZ_deg(0.5, 1.5, ref_cell_center_coordinates[0][1]);
      LLtoXYZ_deg(0.6, 0.6, ref_cell_center_coordinates[1][0]);
      LLtoXYZ_deg(0.6, 1.6, ref_cell_center_coordinates[1][1]);

      struct yac_interp_field interp_fields[2] = {
        {
          .location = YAC_LOC_CELL,
          .coordinates_idx =
            yac_basic_grid_add_coordinates_nocpy_f2c(
              grid, YAC_LOC_CELL, TO_POINTER(ref_cell_center_coordinates[0])),
          .masks_idx =
            yac_basic_grid_add_mask_nocpy_f2c(
              grid, YAC_LOC_CELL, TO_POINTER(ref_cell_mask[0]), NULL)
        },
        {
          .location = YAC_LOC_CELL,
          .coordinates_idx =
            yac_basic_grid_add_coordinates_f2c(
              grid, YAC_LOC_CELL, &(ref_cell_center_coordinates[1][0][0]), 2),
          .masks_idx =
            yac_basic_grid_add_mask_f2c(
              grid, YAC_LOC_CELL, ref_cell_mask[1], 2, NULL)
        }
      };
      enum {NUM_FIELDS = sizeof(interp_fields) / sizeof(interp_fields[0])};

      for (size_t loc_idx = 0; loc_idx < NUM_LOCACTIONS; ++loc_idx) {
        if (yac_basic_grid_get_data_size(grid, locations[loc_idx]) !=
            yac_basic_grid_get_data_size_f2c(grid, (int)(locations[loc_idx])))
          PUT_ERR("ERROR in yac_basic_grid_get_data_size_f2c");
      }

      for (size_t i = 0; i < NUM_FIELDS; ++i)
        check_interp_field(
          grid, interp_fields[i],
          (yac_const_coordinate_pointer)(ref_cell_center_coordinates[i]),
          ref_cell_mask[i]);
    }

    // cleanup
    yac_basic_grid_delete(grid);
  }

   return TEST_EXIT_CODE;
}

static void check_basic_grid(
  struct yac_basic_grid * grid, char const * ref_grid_name,
  struct yac_basic_grid_data ref_grid_data,
  int const * const * ref_core_masks,
  size_t const * ref_grid_data_sizes) {

  // yac_basic_grid_get_name
  if (strcmp(ref_grid_name, yac_basic_grid_get_name(grid)))
    PUT_ERR("ERROR in yac_basic_grid_get_name");

  // yac_basic_grid_get_data
  check_basic_grid_data(
    ref_grid_data, *yac_basic_grid_get_data(grid), ref_grid_name);

  for (size_t loc_idx = 0; loc_idx < NUM_LOCACTIONS; ++loc_idx) {

    // yac_basic_grid_get_core_mask
    if (ref_core_masks[loc_idx] !=
        yac_basic_grid_get_core_mask(grid, locations[loc_idx]))
      PUT_ERR("ERROR in yac_basic_grid_get_core_mask");

    // yac_basic_grid_get_data_size
    if (ref_grid_data_sizes[loc_idx] !=
        yac_basic_grid_get_data_size(grid, locations[loc_idx]))
      PUT_ERR("ERROR in yac_basic_grid_get_data_size");
  }

  size_t num_cells = ref_grid_data.num_cells;
  double * cell_areas = malloc(num_cells * sizeof(*cell_areas));
  for (size_t i = 0; i < num_cells; ++i) cell_areas[i] = DBL_MAX;
  yac_basic_grid_compute_cell_areas(grid, cell_areas);
  for (size_t i = 0; i < num_cells; ++i)
    if (cell_areas[i] == DBL_MAX)
      PUT_ERR("ERROR in yac_basic_grid_compute_cell_areas");
  free(cell_areas);
}

static void check_interp_field(
  struct yac_basic_grid * grid, struct yac_interp_field field,
  yac_const_coordinate_pointer ref_coords, int const * ref_mask) {

  size_t data_size = yac_basic_grid_get_data_size(grid, field.location);

  yac_const_coordinate_pointer coords =
    yac_basic_grid_get_field_coordinates(grid, field);

  if ((ref_coords == NULL) != (coords == NULL))
    PUT_ERR("ERROR in yac_basic_grid_get_field_coordinates");

  if (ref_coords != NULL)
    for (size_t i = 0; i < data_size; ++i)
      if (get_vector_angle(ref_coords[i], coords[i]) > yac_angle_tol)
        PUT_ERR("ERROR in yac_basic_grid_get_field_coordinates");

  int const * mask =
    yac_basic_grid_get_field_mask(grid, field);

  if ((ref_mask == NULL) != (mask == NULL))
    PUT_ERR("ERROR in yac_basic_grid_get_field_mask");

  if (ref_mask != NULL)
    for (size_t i = 0; i < data_size; ++i)
      if (ref_mask[i] != mask[i])
        PUT_ERR("ERROR in yac_basic_grid_get_field_mask");
}
