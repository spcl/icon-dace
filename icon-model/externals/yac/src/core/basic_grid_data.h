// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BASIC_GRID_DATA_H
#define BASIC_GRID_DATA_H

#include "yac_types.h"
#include "grid_cell.h"

// YAC PUBLIC HEADER START

/** \example test_basic_grid_data.c
 */

struct yac_basic_grid_data {
  yac_coordinate_pointer vertex_coordinates;
  yac_int * cell_ids;
  yac_int * vertex_ids;
  yac_int * edge_ids;
  size_t num_cells; // number of local cells (owned by local process)
  size_t num_vertices; // number of local vertices (owned by local process)
  size_t num_edges; // number of local edges (owned by local process)
  int * core_cell_mask;
  int * core_vertex_mask;
  int * core_edge_mask;
  int * num_vertices_per_cell;
  int * num_cells_per_vertex;
  size_t * cell_to_vertex;
  size_t * cell_to_vertex_offsets;
  size_t * cell_to_edge;
  size_t * cell_to_edge_offsets;
  size_t * vertex_to_cell;
  size_t * vertex_to_cell_offsets;
  yac_size_t_2_pointer edge_to_vertex;
  enum yac_edge_type * edge_type;
  size_t num_total_cells; // number of locally stored cells
  size_t num_total_vertices; // number of locally stored vertices
  size_t num_total_edges; // number of locally stored edges
};

struct yac_basic_grid_data yac_generate_basic_grid_data_reg_2d(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid_data yac_generate_basic_grid_data_reg_2d_deg(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid_data yac_generate_basic_grid_data_curve_2d(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid_data yac_generate_basic_grid_data_curve_2d_deg(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct yac_basic_grid_data yac_generate_basic_grid_data_unstruct(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct yac_basic_grid_data yac_generate_basic_grid_data_unstruct_deg(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct yac_basic_grid_data yac_generate_basic_grid_data_unstruct_ll(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct yac_basic_grid_data yac_generate_basic_grid_data_unstruct_ll_deg(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct yac_basic_grid_data yac_generate_basic_grid_data_cloud(
  size_t nbr_points, double * x_points, double * y_points);

struct yac_basic_grid_data yac_generate_basic_grid_data_cloud_deg(
  size_t nbr_points, double * x_points, double * y_points);

void yac_basic_grid_data_compute_cell_areas(
  struct yac_basic_grid_data grid, double * cell_areas);

void yac_basic_grid_data_free(struct yac_basic_grid_data grid);

// YAC PUBLIC HEADER STOP

#endif // BASIC_GRID_DATA_H
