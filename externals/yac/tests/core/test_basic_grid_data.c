// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "tests.h"
#include "basic_grid_data.h"
#include "geometry.h"
#include "ensure_array_size.h"
#include "area.h"
#include "test_common.h"

enum grid_type {

   REG2D,
   UNSTRUCT,
   CURVE2D,
};

static void test_grid_data_1x1(
  struct yac_basic_grid_data grid, char * grid_name, enum grid_type type) {

  double ref_coords_x[4] = {0.5,1.5,0.5,1.5};
  double ref_coords_y[4] = {0.5,0.5,1.5,1.5};
  double vertex_coordinates[4][3];
  for (size_t i = 0; i < 4; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  enum yac_edge_type reg2d_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE};
  enum yac_edge_type unstruct_edge_type[] =
    {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type * ref_edge_type =
    (type == REG2D)?reg2d_edge_type:unstruct_edge_type;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 1,
    .num_vertices = 4,
    .num_edges = 4,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4},
    .num_cells_per_vertex = (int[]){1,1,1,1},
    .cell_to_vertex = (size_t[]){0,1,3,2},
    .cell_to_vertex_offsets = (size_t[]){0},
    .cell_to_edge = (size_t[]){0,2,3,1},
    .cell_to_edge_offsets = (size_t[]){0},
    .vertex_to_cell = (size_t[]){0,0,0,0},
    .vertex_to_cell_offsets = (size_t[]){0,1,2,3},
    .edge_to_vertex = (size_t[][2]){{0,1},{0,2},{1,3},{2,3}},
    .edge_type = ref_edge_type,
    .num_total_cells = 1,
    .num_total_vertices = 4,
    .num_total_edges = 4};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_2x2(
  struct yac_basic_grid_data grid, char * grid_name, enum grid_type type) {

  double ref_coords_x[9] = {0,1,2,0,1,2,0,1,2};
  double ref_coords_y[9] = {0,0,0,1,1,1,2,2,2};
  double vertex_coordinates[9][3];
  for (size_t i = 0; i < 9; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  enum yac_edge_type reg2d_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE};
  enum yac_edge_type unstruct_edge_type[] =
    {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type * ref_edge_type =
    (type == REG2D)?reg2d_edge_type:unstruct_edge_type;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 4,
    .num_vertices = 9,
    .num_edges = 12,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4,4,4,4},
    .num_cells_per_vertex = (int[]){1,2,1, 2,4,2, 1,2,1},
    .cell_to_vertex = (size_t[]){0,1,4,3, 1,2,5,4, 3,4,7,6, 4,5,8,7},
    .cell_to_vertex_offsets = (size_t[]){0,4,8,12},
    .cell_to_edge = (size_t[]){0,3,5,1, 2,4,7,3, 5,8,10,6, 7,9,11,8},
    .cell_to_edge_offsets = (size_t[]){0,4,8,12},
    .vertex_to_cell = (size_t[]){0, 0,1, 1, 0,2, 0,1,2,3, 1,3, 2, 2,3, 3},
    .vertex_to_cell_offsets = (size_t[]){0, 1, 3, 4, 6, 10, 12, 13, 15, 16},
    .edge_to_vertex = (size_t[][2]){{0,1},{0,3},{1,2},{1,4},{2,5},
                                    {3,4},{3,6},{4,5},{4,7},{5,8},
                                    {6,7},{7,8}},
    .edge_type = ref_edge_type,
    .num_total_cells = 4,
    .num_total_vertices = 9,
    .num_total_edges = 12};

  double ref_cell_areas[4];
  {
    double cell_coords_x[2][4] = {{0,1,1,0},{1,2,2,1}};
    double cell_coords_y[2][4] = {{0,0,1,1},{1,1,2,2}};
    enum yac_edge_type edge_types[2][4] =
      {{YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
        YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE},
       {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
        YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE}};
    for (size_t i = 0, k = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j, ++k) {
        struct yac_grid_cell cell =
          generate_cell_deg(
            cell_coords_x[j], cell_coords_y[i], edge_types[type == REG2D], 4);
        ref_cell_areas[k] = yac_huiliers_area(cell);
        yac_free_grid_cell(&cell);
      }
    }
  }

  double cell_areas[4];
  yac_basic_grid_data_compute_cell_areas(grid, cell_areas);

  for (size_t i = 0; i < 4; ++i)
    if (fabs(ref_cell_areas[i] - cell_areas[i]) > YAC_AREA_TOL)
      PUT_ERR("ERROR in yac_basic_grid_data_compute_cell_areas");

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_2x3(
  struct yac_basic_grid_data grid, char * grid_name, enum grid_type type) {

  double ref_coords_x[12] = {0,1,2,0,1,2,0,1,2,0,1,2};
  double ref_coords_y[12] = {0,0,0,1,1,1,2,2,2,3,3,3};
  double vertex_coordinates[12][3];
  for (size_t i = 0; i < 12; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  enum yac_edge_type reg2d_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE};
  enum yac_edge_type unstruct_edge_type[] =
    {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type * ref_edge_type =
    (type == REG2D)?reg2d_edge_type:unstruct_edge_type;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 6,
    .num_vertices = 12,
    .num_edges = 17,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4,4,4,4,4,4},
    .num_cells_per_vertex = (int[]){1,2,1, 2,4,2, 2,4,2, 1,2,1},
    .cell_to_vertex = (size_t[]){0,1,4,3, 1,2,5,4, 3,4,7,6, 4,5,8,7,
                                 6,7,10,9, 7,8,11,10},
    .cell_to_vertex_offsets = (size_t[]){0,4,8,12,16,20},
    .cell_to_edge = (size_t[]){0,3,5,1, 2,4,7,3, 5,8,10,6,
                               7,9,12,8, 10,13,15,11, 12,14,16,13},
    .cell_to_edge_offsets = (size_t[]){0,4,8,12,16,20},
    .vertex_to_cell = (size_t[]){0, 0,1, 1, 0,2, 0,1,2,3, 1,3,
                                 2,4, 2,3,4,5, 3,5, 4, 4,5, 5},
    .vertex_to_cell_offsets = (size_t[]){0,1,3,4,6,10,12,14,18,20,21,23,24},
    .edge_to_vertex = (size_t[][2]){{0,1},{0,3},{1,2},{1,4},{2,5},
                                    {3,4},{3,6},{4,5},{4,7},{5,8},
                                    {6,7},{6,9},{7,8},{7,10},{8,11},
                                    {9,10},{10,11}},
    .edge_type = ref_edge_type,
    .num_total_cells = 6,
    .num_total_vertices = 12,
    .num_total_edges = 17};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_3x3(
  struct yac_basic_grid_data grid, char * grid_name, enum grid_type type) {

  double ref_coords_x[16] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
  double ref_coords_y[16] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
  double vertex_coordinates[16][3];
  for (size_t i = 0; i < 16; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  enum yac_edge_type reg2d_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE};
  enum yac_edge_type unstruct_edge_type[] =
    {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type * ref_edge_type =
    (type == REG2D)?reg2d_edge_type:unstruct_edge_type;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 9,
    .num_vertices = 16,
    .num_edges = 24,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4,4,4,4,4,4,4,4,4},
    .num_cells_per_vertex = (int[]){1,2,2,1, 2,4,4,2, 2,4,4,2, 1,2,2,1},
    .cell_to_vertex = (size_t[]){0,1,5,4, 1,2,6,5, 2,3,7,6,
                                 4,5,9,8, 5,6,10,9, 6,7,11,10,
                                 8,9,13,12, 9,10,14,13, 10,11,15,14},
    .cell_to_vertex_offsets = (size_t[]){0,4,8,12,16,20,24,28,32},
    .cell_to_edge = (size_t[]){0,3,7,1, 2,5,9,3, 4,6,11,5,
                               7,10,14,8, 9,12,16,10, 11,13,18,12,
                               14,17,21,15, 16,19,22,17, 18,20,23,19},
    .cell_to_edge_offsets = (size_t[]){0,4,8, 12,16,20, 24,28,32},
    .vertex_to_cell = (size_t[]){0, 0,1, 1,2, 2,
                                 0,3, 0,1,3,4, 1,2,4,5, 2,5,
                                 3,6, 3,4,6,7, 4,5,7,8, 5,8,
                                 6, 6,7, 7,8, 8},
    .vertex_to_cell_offsets = (size_t[]){0,1,3,5,
                                         6,8,12,16,
                                         18,20,24,28,
                                         30,31,33,35},
    .edge_to_vertex = (size_t[][2]){{0,1},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},
                                    {4,5},{4,8},{5,6},{5,9},{6,7},{6,10},{7,11},
                                    {8,9},{8,12},{9,10},{9,13},{10,11},{10,14},{11,15},
                                    {12,13},{13,14},{14,15}},
    .edge_type = ref_edge_type,
    .num_total_cells = 9,
    .num_total_vertices = 16,
    .num_total_edges = 24};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_3x3_cloud(
  struct yac_basic_grid_data grid, char * grid_name) {

  double ref_coords_x[16] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
  double ref_coords_y[16] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
  double vertex_coordinates[16][3];
  for (size_t i = 0; i < 16; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 0,
    .num_vertices = 16,
    .num_edges = 0,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = NULL,
    .num_cells_per_vertex = (int[]){0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
    .cell_to_vertex = NULL,
    .cell_to_vertex_offsets = NULL,
    .cell_to_edge = NULL,
    .cell_to_edge_offsets = NULL,
    .vertex_to_cell = NULL,
    .vertex_to_cell_offsets = (size_t[]){0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
    .edge_to_vertex = NULL,
    .edge_type = NULL,
    .num_total_cells = 0,
    .num_total_vertices = 16,
    .num_total_edges = 0};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_3x2(
  struct yac_basic_grid_data grid, char * grid_name, enum grid_type type) {

  double ref_coords_x[12] = {0,1,2,3,0,1,2,3,0,1,2,3};
  double ref_coords_y[12] = {0,0,0,0,1,1,1,1,2,2,2,2};
  double vertex_coordinates[12][3];
  for (size_t i = 0; i < 12; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  enum yac_edge_type reg2d_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE};
  enum yac_edge_type unstruct_edge_type[] =
    {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type * ref_edge_type =
    (type == REG2D)?reg2d_edge_type:unstruct_edge_type;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 6,
    .num_vertices = 12,
    .num_edges = 17,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4,4,4,4,4,4},
    .num_cells_per_vertex = (int[]){1,2,2,1, 2,4,4,2, 1,2,2,1},
    .cell_to_vertex = (size_t[]){0,1,5,4, 1,2,6,5, 2,3,7,6,
                                 4,5,9,8, 5,6,10,9, 6,7,11,10},
    .cell_to_vertex_offsets = (size_t[]){0,4,8,12,16,20},
    .cell_to_edge = (size_t[]){0,3,7,1, 2,5,9,3, 4,6,11,5,
                               7,10,14,8, 9,12,15,10, 11,13,16,12},
    .cell_to_edge_offsets = (size_t[]){0,4,8, 12,16,20},
    .vertex_to_cell = (size_t[]){0, 0,1, 1,2, 2,
                                 0,3, 0,1,3,4, 1,2,4,5, 2,5,
                                 3, 3,4, 4,5, 5},
    .vertex_to_cell_offsets = (size_t[]){0,1,3,5,
                                         6,8,12,16,
                                         18,19,21,23},
    .edge_to_vertex = (size_t[][2]){{0,1},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},
                                    {4,5},{4,8},{5,6},{5,9},{6,7},{6,10},{7,11},
                                    {8,9},{9,10},{10,11}},
    .edge_type = ref_edge_type,
    .num_total_cells = 6,
    .num_total_vertices = 12,
    .num_total_edges = 17};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_1x4(
  struct yac_basic_grid_data grid, char * grid_name, enum grid_type type) {

  double ref_coords_x[10] = {0,1,0,1,0,1,0,1,0,1};
  double ref_coords_y[10] = {0,0,1,1,2,2,3,3,4,4};
  double vertex_coordinates[10][3];
  for (size_t i = 0; i < 10; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  enum yac_edge_type reg2d_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE};
  enum yac_edge_type unstruct_edge_type[] =
    {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type * ref_edge_type =
    (type == REG2D)?reg2d_edge_type:unstruct_edge_type;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 4,
    .num_vertices = 10,
    .num_edges = 13,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4,4,4,4},
    .num_cells_per_vertex = (int[]){1,1, 2,2, 2,2, 2,2, 1,1},
    .cell_to_vertex = (size_t[]){0,1,3,2, 2,3,5,4, 4,5,7,6, 6,7,9,8},
    .cell_to_vertex_offsets = (size_t[]){0,4,8,12},
    .cell_to_edge = (size_t[]){0,2,3,1, 3,5,6,4, 6,8,9,7, 9,11,12,10},
    .cell_to_edge_offsets = (size_t[]){0,4,8,12},
    .vertex_to_cell = (size_t[]){0, 0, 0,1, 0,1, 1,2, 1,2, 2,3, 2,3, 3, 3},
    .vertex_to_cell_offsets = (size_t[]){0,1,2,4,6,8,10,12,14,15},
    .edge_to_vertex = (size_t[][2]){{0,1},{0,2},{1,3},
                                    {2,3},{2,4},{3,5},
                                    {4,5},{4,6},{5,7},
                                    {6,7},{6,8},{7,9},
                                    {8,9}},
    .edge_type = ref_edge_type,
    .num_total_cells = 4,
    .num_total_vertices = 10,
    .num_total_edges = 13};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_4x1(
  struct yac_basic_grid_data grid, char * grid_name, enum grid_type type) {

  double ref_coords_x[10] = {0,1,2,3,4,0,1,2,3,4};
  double ref_coords_y[10] = {0,0,0,0,0,1,1,1,1,1};
  double vertex_coordinates[10][3];
  for (size_t i = 0; i < 10; ++i)
    LLtoXYZ_deg(ref_coords_x[i], ref_coords_y[i], vertex_coordinates[i]);

  enum yac_edge_type reg2d_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE};
  enum yac_edge_type unstruct_edge_type[] =
    {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
     YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type * ref_edge_type =
    (type == REG2D)?reg2d_edge_type:unstruct_edge_type;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 4,
    .num_vertices = 10,
    .num_edges = 13,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4,4,4,4},
    .num_cells_per_vertex = (int[]){1,2,2,2,1, 1,2,2,2,1},
    .cell_to_vertex = (size_t[]){0,1,6,5, 1,2,7,6, 2,3,8,7, 3,4,9,8},
    .cell_to_vertex_offsets = (size_t[]){0,4,8,12},
    .cell_to_edge = (size_t[]){0,3,9,1, 2,5,10,3, 4,7,11,5, 6,8,12,7},
    .cell_to_edge_offsets = (size_t[]){0,4,8,12},
    .vertex_to_cell = (size_t[]){0, 0,1, 1,2, 2,3, 3,
                                 0, 0,1, 1,2, 2,3, 3},
    .vertex_to_cell_offsets = (size_t[]){0,1,3,5,7, 8,9,11,13,15},
    .edge_to_vertex = (size_t[][2]){{0,1},{0,5},{1,2},{1,6},{2,3},{2,7},{3,4},{3,8},{4,9},
                                    {5,6},{6,7},{7,8},{8,9}},
    .edge_type = ref_edge_type,
    .num_total_cells = 4,
    .num_total_vertices = 10,
    .num_total_edges = 13};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void test_grid_data_6x4_cyclic(
  struct yac_basic_grid_data grid, char * grid_name) {

  double ref_coords_x[6] = {0,60,120,180,240,300};
  double ref_coords_y[5] = {-89,-45,0,45,89};
  double vertex_coordinates[30][3];
  for (size_t i = 0, k = 0; i < 5; ++i)
    for (size_t j = 0; j < 6; ++j, ++k)
      LLtoXYZ_deg(ref_coords_x[j], ref_coords_y[i], vertex_coordinates[k]);

  enum yac_edge_type ref_edge_type[] =
    {YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LON_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,YAC_LON_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE,
     YAC_LAT_CIRCLE_EDGE,YAC_LAT_CIRCLE_EDGE};
  size_t cell_to_vertex_offsets[24];
  for (size_t i = 0; i < 24; ++i) cell_to_vertex_offsets[i] = 4 * i;
  size_t vertex_to_cell_offsets[30];
  for (size_t i = 0; i < 6; ++i) vertex_to_cell_offsets[i] = 2 * i;
  for (size_t i = 6; i < 24; ++i) vertex_to_cell_offsets[i] = 4 * (i - 6) + 6 * 2;
  for (size_t i = 24; i < 30; ++i) vertex_to_cell_offsets[i] = 2 * (i - 24) + 6 * 2 + 18 * 4 ;

  struct yac_basic_grid_data ref_grid_data = {
    .vertex_coordinates = vertex_coordinates,
    .cell_ids = NULL,
    .vertex_ids = NULL,
    .edge_ids = NULL,
    .num_cells = 24,
    .num_vertices = 30,
    .num_edges = 54,
    .core_cell_mask = NULL,
    .core_vertex_mask = NULL,
    .core_edge_mask = NULL,
    .num_vertices_per_cell = (int[]){4,4,4,4,4,4,
                                     4,4,4,4,4,4,
                                     4,4,4,4,4,4,
                                     4,4,4,4,4,4},
    .num_cells_per_vertex = (int[]){2,2,2,2,2,2,
                                    4,4,4,4,4,4,
                                    4,4,4,4,4,4,
                                    4,4,4,4,4,4,
                                    2,2,2,2,2,2},
    .cell_to_vertex = (size_t[]){0,1,7,6, 1,2,8,7, 2,3,9,8, 3,4,10,9, 4,5,11,10, 5,0,6,11,
                                 6,7,13,12, 7,8,14,13, 8,9,15,14, 9,10,16,15, 10,11,17,16, 11,6,12,17,
                                 12,13,19,18, 13,14,20,19, 14,15,21,20, 15,16,22,21, 16,17,23,22, 17,12,18,23,
                                 18,19,25,24, 19,20,26,25, 20,21,27,26, 21,22,28,27, 22,23,29,28, 23,18,24,29},
    .cell_to_vertex_offsets = cell_to_vertex_offsets,
    .cell_to_edge = (size_t[]){0,4,12,2, 3,6,15,4, 5,8,17,6, 7,10,19,8, 9,11,21,10, 1,2,13,11,
                               12,16,24,14, 15,18,27,16, 17,20,29,18, 19,22,31,20, 21,23,33,22, 13,14,25,23,
                               24,28,36,26, 27,30,39,28, 29,32,41,30, 31,34,43,32, 33,35,45,34, 25,26,37,35,
                               36,40,48,38, 39,42,50,40, 41,44,51,42, 43,46,52,44, 45,47,53,46, 37,38,49,47},
    .cell_to_edge_offsets = cell_to_vertex_offsets,
    .vertex_to_cell = (size_t[]){0,5, 0,1, 1,2, 2,3, 3,4, 4,5,
                                 0,6,11,5, 0,1,7,6, 1,2,8,7, 2,3,9,8, 3,4,10,9, 4,5,11,10,
                                 6,12,17,11, 6,7,13,12, 7,8,14,13, 8,9,15,14, 9,10,16,15, 10,11,17,16,
                                 12,18,23,17, 12,13,19,18, 13,14,20,19, 14,15,21,20, 15,16,22,21, 16,17,23,22,
                                 18,23, 18,19, 19,20, 20,21, 21,22, 22,23},
    .vertex_to_cell_offsets = vertex_to_cell_offsets,
    .edge_to_vertex = (size_t[][2]){{0,1},{0,5},{0,6},{1,2},{1,7},{2,3},{2,8},{3,4},{3,9},{4,5},{4,10},{5,11},
                                    {6,7},{6,11},{6,12},{7,8},{7,13},{8,9},{8,14},{9,10},{9,15},{10,11},{10,16},{11,17},
                                    {12,13},{12,17},{12,18},{13,14},{13,19},{14,15},{14,20},{15,16},{15,21},{16,17},{16,22},{17,23},
                                    {18,19},{18,23},{18,24},{19,20},{19,25},{20,21},{20,26},{21,22},{21,27},{22,23},{22,28},{23,29},
                                    {24,25},{24,29},{25,26},{26,27},{27,28},{28,29}},
    .edge_type = ref_edge_type,
    .num_total_cells = 24,
    .num_total_vertices = 30,
    .num_total_edges = 54};

  check_basic_grid_data(grid, ref_grid_data, grid_name);
}

static void check_cell_triangulation(
  double * coordinates_x, double * coordinates_y, size_t num_corners) {

  enum yac_edge_type tmp_edge_type[num_corners];
  double tmp_coordinates_xyz[num_corners][3];
  for (size_t i = 0; i < num_corners; ++i)
    tmp_edge_type[i] = YAC_GREAT_CIRCLE_EDGE;

  struct yac_grid_cell cell =
    generate_cell_deg(
      coordinates_x, coordinates_y, tmp_edge_type, num_corners);

  double cell_area = yac_huiliers_area(cell);

  struct yac_grid_cell triangles[num_corners-2];
  size_t triangle_indices[num_corners-2][3];
  for (size_t j = 0; j < num_corners-2; ++j)
    yac_init_grid_cell(&(triangles[j]));

  struct yac_grid_cell tmp_cell =
    {.coordinates_xyz = tmp_coordinates_xyz,
     .edge_type = tmp_edge_type,
     .num_corners = num_corners,
     .array_size = num_corners};
  size_t corner_indices[num_corners];

  // check both directions
  for (int order = -1; order < 2; order += 2) {

    // for all start corners
    for (int j = 0; j < (int)num_corners; ++j) {

      for (int k = 0; k < (int)num_corners; ++k) {
        int idx = (order*k+j+((int)num_corners))%((int)num_corners);
        tmp_cell.coordinates_xyz[k][0] = cell.coordinates_xyz[idx][0];
        tmp_cell.coordinates_xyz[k][1] = cell.coordinates_xyz[idx][1];
        tmp_cell.coordinates_xyz[k][2] = cell.coordinates_xyz[idx][2];
        corner_indices[k] = (size_t)idx;
      }

      // for all start corners of the triangulation
      for (size_t k = 0; k < num_corners; ++k) {

        // compute triangulation
        yac_triangulate_cell(tmp_cell, k, triangles);
        yac_triangulate_cell_indices(
          corner_indices, num_corners, k, triangle_indices);

        for (size_t l = 0; l < num_corners - 2; ++l)
          for (size_t m = 0; m < 3; ++m)
            if (memcmp(
                  cell.coordinates_xyz[triangle_indices[l][m]],
                  triangles[l].coordinates_xyz[m], 3 * sizeof(double)))
              PUT_ERR("error in yac_triangulate_cell_indices");

        double total_area = 0.0;

        // for all triangles
        for (size_t l = 0; l < num_corners - 2; ++l) {

          if (triangles[l].num_corners != 3)
            PUT_ERR("wrong number of corners in triangle\n");

          total_area += yac_huiliers_area(triangles[l]);
        }

        if (fabs(total_area - cell_area) > 1e-9)
          PUT_ERR("error in cell triangulation\n");
      }
    }
  }

  yac_free_grid_cell(&cell);
  for (size_t j = 0; j < num_corners-2; ++j)
    yac_free_grid_cell(&(triangles[j]));
}

int main (void) {

   { // test setting up of 2d regular 1x1 grids

      double coord_x[2] = {0.5,1.5};
      double coord_y[2] = {0.5,1.5};

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_1x1(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }

   { // test setting up of 2d unstructured 1x1 grids

      size_t num_vertices = 4;
      size_t num_cells = 1;
      double coord_x[4] = {0.5,1.5,0.5,1.5};
      double coord_y[4] = {0.5,0.5,1.5,1.5};
      int num_vertices_per_cell[] = {4};
      int cell_to_vertex[] = {0,1,3,2};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_1x1(unstruct_grid, "unstructured 2d grid", UNSTRUCT);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of 2d curvilinear 1x1 grids

      double coord_x[4] = {0.5,1.5,0.5,1.5};
      double coord_y[4] = {0.5,0.5,1.5,1.5};

      size_t num_vertices[2] = {2,2};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data curvi_grid =
        yac_generate_basic_grid_data_curve_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_1x1(curvi_grid, "curvilinear 2d grid", CURVE2D);

      yac_basic_grid_data_free(curvi_grid);
   }

   { // test setting up of 2d unstructured lonlat 1x1 grids

      size_t num_vertices = 4;
      size_t num_cells = 1;
      double coord_x[4] = {0.5,1.5,0.5,1.5};
      double coord_y[4] = {0.5,0.5,1.5,1.5};
      int num_vertices_per_cell[] = {4};
      int cell_to_vertex[] = {0,1,3,2};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_1x1(unstruct_grid, "unstructured lonlat 2d grid", REG2D);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of 2d regular 2x2 grids

      double coord_x[3] = {0,1,2};
      double coord_y[3] = {0,1,2};

      size_t num_vertices[2] = {3,3};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_2x2(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }

   { // test setting up of unstructured 2x2 grids

      size_t num_vertices = 9;
      size_t num_cells = 4;
      double coord_x[9] = {0,1,2,0,1,2,0,1,2};
      double coord_y[9] = {0,0,0,1,1,1,2,2,2};
      int num_vertices_per_cell[4] = {4,4,4,4};
      int cell_to_vertex[16] = {0,1,4,3, 1,2,5,4, 3,4,7,6, 4,5,8,7};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_2x2(unstruct_grid, "unstructured 2d grid", UNSTRUCT);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of curvilinear 2x2 grids

      double coord_x[9] = {0,1,2,0,1,2,0,1,2};
      double coord_y[9] = {0,0,0,1,1,1,2,2,2};

      size_t num_vertices[2] = {3,3};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data curvi_grid =
        yac_generate_basic_grid_data_curve_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_2x2(curvi_grid, "curvilinear 2d grid", CURVE2D);

      yac_basic_grid_data_free(curvi_grid);
   }

   { // test setting up of unstructured lonlat 2x2 grids

      size_t num_vertices = 9;
      size_t num_cells = 4;
      double coord_x[9] = {0,1,2,0,1,2,0,1,2};
      double coord_y[9] = {0,0,0,1,1,1,2,2,2};
      int num_vertices_per_cell[4] = {4,4,4,4};
      int cell_to_vertex[16] = {0,1,4,3, 1,2,5,4, 3,4,7,6, 4,5,8,7};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_2x2(unstruct_grid, "unstructured lonlat 2d grid", REG2D);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of 2d regular 2x3 grids

      double coord_x[3] = {0,1,2};
      double coord_y[4] = {0,1,2,3};

      size_t num_vertices[2] = {3,4};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_2x3(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }

   { // test setting up of unstructured 2x3 grids

      size_t num_vertices = 12;
      size_t num_cells = 6;
      double coord_x[12] = {0,1,2,0,1,2,0,1,2,0,1,2};
      double coord_y[12] = {0,0,0,1,1,1,2,2,2,3,3,3};
      int num_vertices_per_cell[6] = {4,4,4,4,4,4};
      int cell_to_vertex[24] = {0,1,4,3, 1,2,5,4,
                                3,4,7,6, 4,5,8,7,
                                6,7,10,9, 7,8,11,10};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_2x3(unstruct_grid, "unstructured 2d grid", UNSTRUCT);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of curvilinear 2x3 grids

      double coord_x[12] = {0,1,2,0,1,2,0,1,2,0,1,2};
      double coord_y[12] = {0,0,0,1,1,1,2,2,2,3,3,3};

      size_t num_vertices[2] = {3,4};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data curvi_grid =
        yac_generate_basic_grid_data_curve_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_2x3(curvi_grid, "curvilinear 2d grid", CURVE2D);

      yac_basic_grid_data_free(curvi_grid);
   }

   { // test setting up of unstructured lonlat 2x3 grids

      size_t num_vertices = 12;
      size_t num_cells = 6;
      double coord_x[12] = {0,1,2,0,1,2,0,1,2,0,1,2};
      double coord_y[12] = {0,0,0,1,1,1,2,2,2,3,3,3};
      int num_vertices_per_cell[6] = {4,4,4,4,4,4};
      int cell_to_vertex[24] = {0,1,4,3, 1,2,5,4,
                                3,4,7,6, 4,5,8,7,
                                6,7,10,9, 7,8,11,10};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_2x3(unstruct_grid, "unstructured lonlat 2d grid", REG2D);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of 2d regular 3x3 grids

      double coord_x[] = {0,1,2,3};
      double coord_y[] = {0,1,2,3};

      size_t num_vertices[2] = {4,4};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_3x3(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }

   { // test setting up of unstructured 3x3 grids

      size_t num_vertices = 16;
      size_t num_cells = 9;
      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
      int num_vertices_per_cell[] = {4,4,4, 4,4,4, 4,4,4};
      int cell_to_vertex[] = {0,1,5,4, 1,2,6,5, 2,3,7,6,
                              4,5,9,8, 5,6,10,9, 6,7,11,10,
                              8,9,13,12, 9,10,14,13, 10,11,15,14};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_3x3(unstruct_grid, "unstructured 2d grid", UNSTRUCT);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of curvilinear 3x3 grids

      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};

      size_t num_vertices[2] = {4,4};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data curvi_grid =
        yac_generate_basic_grid_data_curve_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_3x3(curvi_grid, "curvilinear 2d grid", CURVE2D);

      yac_basic_grid_data_free(curvi_grid);
   }

   { // test setting up of unstructured lonlat 3x3 grids

      size_t num_vertices = 16;
      size_t num_cells = 9;
      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
      int num_vertices_per_cell[] = {4,4,4, 4,4,4, 4,4,4};
      int cell_to_vertex[] = {0,1,5,4, 1,2,6,5, 2,3,7,6,
                              4,5,9,8, 5,6,10,9, 6,7,11,10,
                              8,9,13,12, 9,10,14,13, 10,11,15,14};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_3x3(unstruct_grid, "unstructured lonlat 2d grid", REG2D);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of 2d regular 3x2 grids

      double coord_x[] = {0,1,2,3};
      double coord_y[] = {0,1,2};

      size_t num_vertices[2] = {4,3};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_3x2(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }

   { // test setting up of unstructured 3x2 grids

      size_t num_vertices = 12;
      size_t num_cells = 6;
      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2};
      int num_vertices_per_cell[] = {4,4,4, 4,4,4};
      int cell_to_vertex[] = {0,1,5,4, 1,2,6,5, 2,3,7,6,
                              4,5,9,8, 5,6,10,9, 6,7,11,10};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_3x2(unstruct_grid, "unstructured 2d grid", UNSTRUCT);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of curvilinear 3x2 grids

      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2};

      size_t num_vertices[2] = {4,3};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data curvi_grid =
        yac_generate_basic_grid_data_curve_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_3x2(curvi_grid, "curvilinear 2d grid", CURVE2D);

      yac_basic_grid_data_free(curvi_grid);
   }

   { // test setting up of unstructured lonlat 3x2 grids

      size_t num_vertices = 12;
      size_t num_cells = 6;
      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2};
      int num_vertices_per_cell[] = {4,4,4, 4,4,4};
      int cell_to_vertex[] = {0,1,5,4, 1,2,6,5, 2,3,7,6,
                              4,5,9,8, 5,6,10,9, 6,7,11,10};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_3x2(unstruct_grid, "unstructured lonlat 2d grid", REG2D);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of 2d regular 1x4 grids

      double coord_x[] = {0,1};
      double coord_y[] = {0,1,2,3,4};

      size_t num_vertices[2] = {2,5};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_1x4(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }

   { // test setting up of unstructured 1x4 grids

      size_t num_vertices = 10;
      size_t num_cells = 4;
      double coord_x[] = {0,1,0,1,0,1,0,1,0,1};
      double coord_y[] = {0,0,1,1,2,2,3,3,4,4};
      int num_vertices_per_cell[] = {4,4,4,4};
      int cell_to_vertex[] = {0,1,3,2, 2,3,5,4, 4,5,7,6, 6,7,9,8};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_1x4(unstruct_grid, "unstructured 2d grid", UNSTRUCT);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of curvilinear 1x4 grids

      double coord_x[] = {0,1,0,1,0,1,0,1,0,1};
      double coord_y[] = {0,0,1,1,2,2,3,3,4,4};

      size_t num_vertices[2] = {2,5};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data curvi_grid =
        yac_generate_basic_grid_data_curve_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_1x4(curvi_grid, "curvilinear 2d grid", CURVE2D);

      yac_basic_grid_data_free(curvi_grid);
   }

   { // test setting up of unstructured lonlat 1x4 grids

      size_t num_vertices = 10;
      size_t num_cells = 4;
      double coord_x[] = {0,1,0,1,0,1,0,1,0,1};
      double coord_y[] = {0,0,1,1,2,2,3,3,4,4};
      int num_vertices_per_cell[] = {4,4,4,4};
      int cell_to_vertex[] = {0,1,3,2, 2,3,5,4, 4,5,7,6, 6,7,9,8};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_1x4(unstruct_grid, "unstructured lonlat 2d grid", REG2D);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of 2d regular 4x1 grids

      double coord_x[] = {0,1,2,3,4};
      double coord_y[] = {0,1};

      size_t num_vertices[2] = {5,2};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_4x1(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }

   { // test setting up of cyclic 2d regular 6x4 grids

      double coord_x[] = {0,60,120,180,240,300};
      double coord_y[] = {-89,-45,0,45,89};

      size_t num_vertices[2] = {6,5};
      int cyclic[2] = {1,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_6x4_cyclic(reg_grid, "regular 2d grid");

      yac_basic_grid_data_free(reg_grid);
   }
#ifdef NOT_YET_SUPPORTED
   { // test setting up of cyclic 2d regular 6x4 grids (that cover the pole)

      double coord_x[] = {0,60,120,180,240};
      double coord_y[] = {-90,45,0,45,90};

      size_t num_vertices[2] = {5,4};
      int cyclic[2] = {1,0};

      struct yac_basic_grid_data reg_grid =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      // test_grid_data_6x4_cyclic_pole(reg_grid, "regular 2d grid", REG2D);

      yac_basic_grid_data_free(reg_grid);
   }
#endif
   { // test setting up of unstructured 4x1 grids

      size_t num_vertices = 10;
      size_t num_cells = 4;
      double coord_x[] = {0,1,2,3,4, 0,1,2,3,4};
      double coord_y[] = {0,0,0,0,0, 1,1,1,1,1};
      int num_vertices_per_cell[] = {4,4,4,4};
      int cell_to_vertex[] = {0,1,6,5, 1,2,7,6, 2,3,8,7, 3,4,9,8};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_4x1(unstruct_grid, "unstructured 2d grid", UNSTRUCT);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of curvilinear 4x1 grids

      double coord_x[] = {0,1,2,3,4, 0,1,2,3,4};
      double coord_y[] = {0,0,0,0,0, 1,1,1,1,1};

      size_t num_vertices[2] = {5,2};
      int cyclic[2] = {0,0};

      struct yac_basic_grid_data curvi_grid =
        yac_generate_basic_grid_data_curve_2d_deg(
          num_vertices, cyclic, coord_x, coord_y);

      test_grid_data_4x1(curvi_grid, "curvilinear 2d grid", CURVE2D);

      yac_basic_grid_data_free(curvi_grid);
   }

   { // test setting up of unstructured lonlat 4x1 grids

      size_t num_vertices = 10;
      size_t num_cells = 4;
      double coord_x[] = {0,1,2,3,4, 0,1,2,3,4};
      double coord_y[] = {0,0,0,0,0, 1,1,1,1,1};
      int num_vertices_per_cell[] = {4,4,4,4};
      int cell_to_vertex[] = {0,1,6,5, 1,2,7,6, 2,3,8,7, 3,4,9,8};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_4x1(unstruct_grid, "unstructured lonlat 2d grid", REG2D);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test periodic coordinates

      size_t num_vertices = 9;
      size_t num_cells = 4;
      double coord_x[9] = {-1,0,1, 719,720,721, -1,0,1};
      double coord_y[9] = { 0,0,0, 1,1,1,        2,2,2};
      int num_vertices_per_cell[4] = {4,4,4,4};
      int cell_to_vertex[16] = {0,1,4,3, 1,2,5,4, 3,4,7,6, 4,5,8,7};

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_ll_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test triangulation of cells

      check_cell_triangulation((double[]){0,1,0}, (double[]){0,0,1}, 3);
      check_cell_triangulation((double[]){0,1,1,0}, (double[]){0,0,1,1}, 4);
      check_cell_triangulation(
        (double[]){0,1,1,0.5,0}, (double[]){0,0,1,1.5,1}, 5);
      check_cell_triangulation(
        (double[]){0,0.5,1,1,0.5,0}, (double[]){0,-0.5,0,1,1.5,1}, 6);
      check_cell_triangulation(
        (double[]){-0.25,0, 0.5,1,1.25,1,0},
        (double[]){ 0.5 ,0,-0.5,0,0.5 ,1,1}, 7);
      check_cell_triangulation(
        (double[]){-0.25,0, 0.5,1,1.25,1,0.5,0},
        (double[]){ 0.5 ,0,-0.5,0,0.5 ,1,1.5,1}, 8);
   }

   { // test setting up of unstructured 3x3 grids without cells

      size_t num_vertices = 16;
      size_t num_cells = 0;
      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
      int * num_vertices_per_cell = NULL;
      int * cell_to_vertex = NULL;

      struct yac_basic_grid_data unstruct_grid =
        yac_generate_basic_grid_data_unstruct_deg(
          num_vertices, num_cells, num_vertices_per_cell,
          coord_x, coord_y, cell_to_vertex);

      test_grid_data_3x3_cloud(unstruct_grid, "unstructured 2d grid cloud");

      yac_basic_grid_data_free(unstruct_grid);
   }

   { // test setting up of unstructured 3x3 grids without cells

      size_t num_coords = 16;
      double coord_x[] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
      double coord_y[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};

      struct yac_basic_grid_data cloud_grid =
        yac_generate_basic_grid_data_cloud_deg(num_coords, coord_x, coord_y);

      test_grid_data_3x3_cloud(cloud_grid, "cloud 2d grid");

      yac_basic_grid_data_free(cloud_grid);
   }

   return TEST_EXIT_CODE;
}

