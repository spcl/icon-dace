// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BASIC_GRID_H
#define BASIC_GRID_H

#include "basic_grid_data.h"
#include "location.h"
#include "field_data.h"

// YAC PUBLIC HEADER START

/** \example test_basic_grid.c
 */

struct yac_interp_field {
  enum yac_location location;
  size_t coordinates_idx;
  size_t masks_idx;
};

struct yac_basic_grid;

struct yac_basic_grid * yac_basic_grid_new(
  char const * name, struct yac_basic_grid_data grid_data);
struct yac_basic_grid * yac_basic_grid_empty_new(char const * name);
yac_const_coordinate_pointer yac_basic_grid_get_field_coordinates(
  struct yac_basic_grid * grid, struct yac_interp_field field);
int const * yac_basic_grid_get_field_mask(
  struct yac_basic_grid * grid, struct yac_interp_field field);
int const * yac_basic_grid_get_core_mask(
  struct yac_basic_grid * grid, enum yac_location location);
char const * yac_basic_grid_get_name(struct yac_basic_grid * grid);
struct yac_basic_grid_data * yac_basic_grid_get_data(
  struct yac_basic_grid * grid);
struct yac_field_data * yac_basic_grid_get_field_data(
  struct yac_basic_grid * grid, enum yac_location location);
size_t yac_basic_grid_get_data_size(
  struct yac_basic_grid * grid, enum yac_location location);
size_t yac_basic_grid_get_named_mask_idx(
  struct yac_basic_grid * grid, enum yac_location location,
  char const * mask_name);
size_t yac_basic_grid_add_coordinates(
  struct yac_basic_grid * grid, enum yac_location location,
  yac_coordinate_pointer coordinates, size_t count);
size_t yac_basic_grid_add_coordinates_nocpy(
  struct yac_basic_grid * grid, enum yac_location location,
  yac_coordinate_pointer coordinates);
size_t yac_basic_grid_add_mask(
  struct yac_basic_grid * grid, enum yac_location location,
  int const * mask, size_t count, char const * mask_name);
size_t yac_basic_grid_add_mask_nocpy(
  struct yac_basic_grid * grid, enum yac_location location,
  int const * mask, char const * mask_name);
void yac_basic_grid_delete(struct yac_basic_grid * grid);

struct yac_basic_grid * yac_basic_grid_reg_2d_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid * yac_basic_grid_reg_2d_deg_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid * yac_basic_grid_curve_2d_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid * yac_basic_grid_curve_2d_deg_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid * yac_basic_grid_unstruct_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex);

struct yac_basic_grid * yac_basic_grid_unstruct_deg_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex);

struct yac_basic_grid * yac_basic_grid_unstruct_ll_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex);

struct yac_basic_grid * yac_basic_grid_unstruct_ll_deg_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex);

struct yac_basic_grid * yac_basic_grid_cloud_new(
  char const * name, size_t nbr_points, double *x_points, double *y_points);

struct yac_basic_grid * yac_basic_grid_cloud_deg_new(
  char const * name, size_t nbr_points, double *x_points, double *y_points);

void yac_basic_grid_to_file_parallel(
  struct yac_basic_grid * grid, char const * filename, MPI_Comm comm);

void yac_basic_grid_compute_cell_areas(
  struct yac_basic_grid * grid, double * cell_areas);

// YAC PUBLIC HEADER STOP

#endif // BASIC_GRID_H
