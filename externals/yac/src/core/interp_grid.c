// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <yaxt.h>

#include "dist_grid_internal.h"
#include "interp_grid_internal.h"
#include "yac_mpi_internal.h"
#include "geometry.h"
#include "utils_core.h"
#include "sphere_part.h"

struct yac_interp_grid {
  char * src_grid_name;
  char * tgt_grid_name;
  struct yac_dist_grid_pair * grid_pair;
  struct yac_interp_field tgt_field;
  size_t num_src_fields;
  struct yac_interp_field src_fields[];
};

struct yac_interp_grid * yac_interp_grid_new(
  struct yac_dist_grid_pair * grid_pair,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t num_src_fields, struct yac_interp_field const * src_fields,
  struct yac_interp_field const tgt_field) {

  struct yac_interp_grid * interp_grid =
    xmalloc(1 * sizeof(*interp_grid) + num_src_fields * sizeof(*src_fields));

  interp_grid->src_grid_name = strdup(src_grid_name);
  interp_grid->tgt_grid_name = strdup(tgt_grid_name);
  interp_grid->grid_pair = grid_pair;
  interp_grid->num_src_fields = num_src_fields;
  memcpy(&interp_grid->tgt_field, &tgt_field, 1 * sizeof(tgt_field));
  memcpy(
    interp_grid->src_fields, src_fields, num_src_fields * sizeof(*src_fields));

  return interp_grid;
}

struct yac_interp_grid * yac_interp_grid_new_f2c(
  struct yac_dist_grid_pair * grid_pair,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t num_src_fields, int const * src_field_locations,
  size_t const * src_field_coordinate_idxs,
  size_t const * src_field_masks_idxs,
  int tgt_field_location, size_t tgt_field_coordinate_idx,
  size_t tgt_field_masks_idx) {

  struct yac_interp_field src_fields[num_src_fields];
  struct yac_interp_field tgt_field;
  for (size_t i = 0; i < num_src_fields; ++i) {
    src_fields[i].location = yac_get_location(src_field_locations[i]);
    src_fields[i].coordinates_idx = src_field_coordinate_idxs[i];
    src_fields[i].masks_idx = src_field_masks_idxs[i];
  }
  tgt_field.location = yac_get_location(tgt_field_location);
  tgt_field.coordinates_idx = tgt_field_coordinate_idx;
  tgt_field.masks_idx = tgt_field_masks_idx;

  return
    yac_interp_grid_new(
      grid_pair, src_grid_name, tgt_grid_name, num_src_fields,
      src_fields, tgt_field);
}

char const * yac_interp_grid_get_src_grid_name(
  struct yac_interp_grid * interp_grid) {
  return interp_grid->src_grid_name;
}

char const * yac_interp_grid_get_tgt_grid_name(
  struct yac_interp_grid * interp_grid) {
  return interp_grid->tgt_grid_name;
}

void yac_interp_grid_get_src_points(
  struct yac_interp_grid * interp_grid, size_t src_field_idx,
  size_t ** src_indices, size_t * count) {

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);
  yac_dist_grid_get_local_unmasked_points(
    dist_grid, interp_grid->src_fields[src_field_idx],
    src_indices, count);
}

void yac_interp_grid_get_tgt_points(
  struct yac_interp_grid * interp_grid, size_t ** tgt_indices,
  size_t * count) {

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);
  yac_dist_grid_get_local_unmasked_points(
    dist_grid, interp_grid->tgt_field, tgt_indices, count);
}

struct remote_point * yac_interp_grid_get_src_remote_points2(
  struct yac_interp_grid * interp_grid, enum yac_location location,
  size_t * src_points, size_t count) {

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  return
    yac_dist_grid_get_remote_points(dist_grid, location, src_points, count);
}

struct remote_point * yac_interp_grid_get_src_remote_points(
  struct yac_interp_grid * interp_grid, size_t src_field_idx,
  size_t * src_points, size_t count) {

  return
    yac_interp_grid_get_src_remote_points2(
      interp_grid, interp_grid->src_fields[src_field_idx].location,
      src_points, count);
}

void yac_interp_grid_src_global_to_local(
  struct yac_interp_grid * interp_grid, size_t src_field_idx,
  yac_int * src_global_ids, size_t count, size_t * src_local_ids) {

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  yac_dist_grid_global_to_local(
    dist_grid, interp_grid->src_fields[src_field_idx].location,
    src_global_ids, count, src_local_ids);
}

void yac_interp_grid_tgt_global_to_local(
  struct yac_interp_grid * interp_grid, yac_int * tgt_global_ids,
  size_t count, size_t * tgt_local_ids) {

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);

  yac_dist_grid_global_to_local(
    dist_grid, interp_grid->tgt_field.location,
    tgt_global_ids, count, tgt_local_ids);
}

struct remote_point * yac_interp_grid_get_tgt_remote_points(
  struct yac_interp_grid * interp_grid, size_t * tgt_points, size_t count) {

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);

  return
    yac_dist_grid_get_remote_points(
      dist_grid, interp_grid->tgt_field.location, tgt_points, count);
}

static yac_int const * get_tgt_grid_global_ids(
  struct yac_interp_grid * interp_grid) {

  struct yac_const_basic_grid_data * yac_basic_grid_data =
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->tgt_grid_name));

  enum yac_location location = interp_grid->tgt_field.location;
  YAC_ASSERT(
    (location == YAC_LOC_CELL) ||
    (location == YAC_LOC_CORNER) ||
    (location == YAC_LOC_EDGE),
    "ERROR(get_tgt_grid_global_ids): invalide target location")

  return yac_basic_grid_data->ids[location];
}

enum yac_location yac_interp_grid_get_tgt_field_location(
  struct yac_interp_grid * interp_grid) {

  return interp_grid->tgt_field.location;
}

void yac_interp_grid_get_src_global_ids(
  struct yac_interp_grid * interp_grid, size_t * src_points, size_t count,
  size_t src_field_idx, yac_int * src_global_ids) {

  const_yac_int_pointer grid_global_ids =
    yac_interp_grid_get_src_field_global_ids(interp_grid, src_field_idx);

  for (size_t i = 0; i < count; ++i)
    src_global_ids[i] = grid_global_ids[src_points[i]];
}

void yac_interp_grid_get_tgt_global_ids(
  struct yac_interp_grid * interp_grid, size_t * tgt_points, size_t count,
  yac_int * tgt_global_ids) {

  yac_int const * grid_global_ids = get_tgt_grid_global_ids(interp_grid);

  for (size_t i = 0; i < count; ++i)
    tgt_global_ids[i] = grid_global_ids[tgt_points[i]];
}

static void yac_interp_grid_get_coordinates(
  size_t * points, size_t count, yac_coordinate_pointer coordinates,
  yac_const_coordinate_pointer grid_coordinates) {

  YAC_ASSERT(
    (grid_coordinates != NULL) || (count == 0),
    "ERROR(yac_interp_grid_get_coordinates): grid_coordinates == NULL")

  for (size_t i = 0; i < count; ++i)
    for (int j = 0; j < 3; ++j)
      coordinates[i][j] = grid_coordinates[points[i]][j];
}

void yac_interp_grid_get_src_coordinates(
  struct yac_interp_grid * interp_grid, size_t * src_points, size_t count,
  size_t src_field_idx, yac_coordinate_pointer src_coordinates) {

  yac_interp_grid_get_coordinates(
    src_points, count, src_coordinates,
    yac_interp_grid_get_src_field_coords(interp_grid, src_field_idx));
}

void yac_interp_grid_get_tgt_coordinates(
  struct yac_interp_grid * interp_grid, size_t * tgt_points, size_t count,
  yac_coordinate_pointer tgt_coordinates) {

  yac_interp_grid_get_coordinates(
    tgt_points, count, tgt_coordinates,
    yac_interp_grid_get_tgt_field_coords(interp_grid));
}

size_t yac_interp_grid_get_num_src_fields(
  struct yac_interp_grid * interp_grid) {

  return interp_grid->num_src_fields;
}

enum yac_location yac_interp_grid_get_src_field_location(
  struct yac_interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_location): invalid src_field_idx")

  return interp_grid->src_fields[src_field_idx].location;
}

const_yac_int_pointer yac_interp_grid_get_src_field_global_ids(
  struct yac_interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_global_ids): invalid src_field_idx")

  struct yac_const_basic_grid_data * yac_basic_grid_data =
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->src_grid_name));

  enum yac_location location =
    interp_grid->src_fields[src_field_idx].location;
  YAC_ASSERT(
    (location == YAC_LOC_CORNER) ||
    (location == YAC_LOC_CELL) ||
    (location == YAC_LOC_EDGE),
    "ERROR(yac_interp_grid_get_src_field_global_ids): "
    "invalid source field location")
  return yac_basic_grid_data->ids[location];
}

yac_const_coordinate_pointer yac_interp_grid_get_src_field_coords(
  struct yac_interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_coords): invalid src_field_idx")

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  return
    (yac_const_coordinate_pointer)
      yac_dist_grid_get_field_coords(
        dist_grid, interp_grid->src_fields[src_field_idx]);
}

yac_const_coordinate_pointer yac_interp_grid_get_tgt_field_coords(
  struct yac_interp_grid * interp_grid) {

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);

  return
    (yac_const_coordinate_pointer)
      yac_dist_grid_get_field_coords(dist_grid, interp_grid->tgt_field);
}

const_int_pointer yac_interp_grid_get_src_field_mask(
  struct yac_interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_mask): invalid src_field_idx")

  struct yac_dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  return
    yac_dist_grid_get_field_mask(
      dist_grid, interp_grid->src_fields[src_field_idx]);
}

void yac_interp_grid_do_points_search(
  struct yac_interp_grid * interp_grid, yac_coordinate_pointer search_coords,
  size_t count, size_t * src_cells) {

  yac_dist_grid_pair_do_point_search(
    interp_grid->grid_pair, interp_grid->src_grid_name, search_coords, count,
    src_cells);
}

void yac_interp_grid_do_points_search_gc(
  struct yac_interp_grid * interp_grid, yac_coordinate_pointer search_coords,
  size_t count, size_t * src_cells) {

  yac_dist_grid_pair_do_point_search_gc(
    interp_grid->grid_pair, interp_grid->src_grid_name, search_coords, count,
    src_cells);
}

void yac_interp_grid_do_nnn_search_src(
  struct yac_interp_grid * interp_grid, yac_coordinate_pointer search_coords,
  size_t count, size_t n, size_t * src_points, double max_search_distance) {

  YAC_ASSERT(
    interp_grid->num_src_fields == 1,
    "ERROR(yac_interp_grid_do_nnn_search_src): invalid number of source fields")

  yac_dist_grid_pair_do_nnn_search(
    interp_grid->grid_pair, interp_grid->src_grid_name, search_coords, count,
    src_points, n, interp_grid->src_fields[0], max_search_distance);
}

void yac_interp_grid_do_nnn_search_tgt(
  struct yac_interp_grid * interp_grid, yac_coordinate_pointer search_coords,
  size_t count, size_t n, size_t * tgt_points, double max_search_distance) {

  yac_dist_grid_pair_do_nnn_search(
    interp_grid->grid_pair, interp_grid->tgt_grid_name, search_coords, count,
    tgt_points, n, interp_grid->tgt_field, max_search_distance);
}

void yac_interp_grid_do_bnd_circle_search_src(
  struct yac_interp_grid * interp_grid,
  const_bounding_circle_pointer bnd_circles,
  size_t count, size_t src_field_idx, size_t ** src_cells,
  size_t * num_src_per_bnd_circle) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_do_bnd_circle_search_src): invalid src_field_idx")

  YAC_ASSERT(
    interp_grid->src_fields[src_field_idx].location == YAC_LOC_CELL,
    "ERROR(yac_interp_grid_do_bnd_circle_search_src): "
    "invalid source field location; has to be YAC_LOC_CELL")

  yac_dist_grid_pair_do_bnd_circle_search(
    interp_grid->grid_pair, interp_grid->src_grid_name, bnd_circles, count,
    src_cells, num_src_per_bnd_circle, interp_grid->src_fields[src_field_idx]);
}

void yac_interp_grid_do_bnd_circle_search_tgt(
  struct yac_interp_grid * interp_grid,
  const_bounding_circle_pointer bnd_circles,
  size_t count, size_t ** tgt_cells, size_t * num_tgt_per_bnd_circle) {

  YAC_ASSERT(
    interp_grid->tgt_field.location == YAC_LOC_CELL,
    "ERROR(yac_interp_grid_do_bnd_circle_search_tgt): "
    "invalid target field location; has to be YAC_LOC_CELL")

  yac_dist_grid_pair_do_bnd_circle_search(
    interp_grid->grid_pair, interp_grid->tgt_grid_name, bnd_circles, count,
    tgt_cells, num_tgt_per_bnd_circle, interp_grid->tgt_field);
}

void yac_interp_grid_do_cell_search_src(
  struct yac_interp_grid * interp_grid, size_t * tgt_cells, size_t count,
  size_t ** src_cells, size_t * num_src_per_tgt) {

  YAC_ASSERT(
    interp_grid->num_src_fields == 1,
    "ERROR(yac_interp_grid_do_cell_search_src): "
    "invalid number of source fields")

  YAC_ASSERT(
    interp_grid->src_fields[0].location == YAC_LOC_CELL,
    "ERROR(yac_interp_grid_do_cell_search_src): "
    "invalid source field location; has to be YAC_LOC_CELL")

  YAC_ASSERT(
    interp_grid->tgt_field.location == YAC_LOC_CELL,
    "ERROR(yac_interp_grid_do_cell_search_src): "
    "invalid target field location; has to be YAC_LOC_CELL")

  yac_dist_grid_pair_do_cell_search(
    interp_grid->grid_pair, interp_grid->tgt_grid_name,
    interp_grid->src_grid_name, tgt_cells, count, src_cells,
    num_src_per_tgt, interp_grid->src_fields[0]);
}

void yac_interp_grid_do_cell_search_tgt(
  struct yac_interp_grid * interp_grid, size_t * src_cells, size_t count,
  size_t ** tgt_cells, size_t * num_tgt_per_src) {

  YAC_ASSERT(
    interp_grid->num_src_fields == 1,
    "ERROR(yac_interp_grid_do_cell_search_tgt): "
    "invalid number of source fields")

  YAC_ASSERT(
    interp_grid->src_fields[0].location == YAC_LOC_CELL,
    "ERROR(yac_interp_grid_do_cell_search_tgt): "
    "invalid source field location; has to be YAC_LOC_CELL")

  YAC_ASSERT(
    interp_grid->tgt_field.location == YAC_LOC_CELL,
    "ERROR(yac_interp_grid_do_cell_search_tgt): "
    "invalid target field location; has to be YAC_LOC_CELL")

  yac_dist_grid_pair_do_cell_search(
    interp_grid->grid_pair, interp_grid->src_grid_name,
    interp_grid->tgt_grid_name, src_cells, count, tgt_cells,
    num_tgt_per_src, interp_grid->tgt_field);
}

MPI_Comm yac_interp_grid_get_MPI_Comm(struct yac_interp_grid * interp_grid) {

  return yac_dist_grid_pair_get_MPI_Comm(interp_grid->grid_pair);
}

struct yac_const_basic_grid_data * yac_interp_grid_get_basic_grid_data_src(
  struct yac_interp_grid * interp_grid) {

  return
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->src_grid_name));
}

struct yac_const_basic_grid_data * yac_interp_grid_get_basic_grid_data_tgt(
  struct yac_interp_grid * interp_grid) {

  return
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->tgt_grid_name));
}

void yac_interp_grid_get_src_cell_neighbours(
  struct yac_interp_grid * interp_grid, size_t * src_cells, size_t count,
  size_t * neighbours) {

  yac_dist_grid_pair_get_cell_neighbours(
    interp_grid->grid_pair, interp_grid->src_grid_name, src_cells, count,
    neighbours);
}

void yac_interp_grid_get_tgt_cell_neighbours(
  struct yac_interp_grid * interp_grid, size_t * tgt_cells, size_t count,
  size_t * neighbours) {

  yac_dist_grid_pair_get_cell_neighbours(
    interp_grid->grid_pair, interp_grid->tgt_grid_name, tgt_cells, count,
    neighbours);
}

void yac_interp_grid_get_src_corner_cells(
  struct yac_interp_grid * interp_grid,
  size_t * src_corners, size_t count, size_t ** src_cells,
  size_t * num_cells_per_corner) {

  yac_dist_grid_pair_get_corner_cells(
    interp_grid->grid_pair, interp_grid->src_grid_name,
    src_corners, count, src_cells, num_cells_per_corner);
}

void yac_interp_grid_get_tgt_corner_cells(
  struct yac_interp_grid * interp_grid,
  size_t * tgt_corners, size_t count, size_t ** tgt_cells,
  size_t * num_cells_per_corner) {

  yac_dist_grid_pair_get_corner_cells(
    interp_grid->grid_pair, interp_grid->tgt_grid_name,
    tgt_corners, count, tgt_cells, num_cells_per_corner);
}

void yac_interp_grid_get_aux_grid_src(
  struct yac_interp_grid * interp_grid, size_t * cells, size_t count,
  size_t ** vertex_to_cell, size_t ** vertex_to_cell_offsets,
  int ** num_cells_per_vertex) {

  YAC_ASSERT(
    interp_grid->num_src_fields == 1,
    "ERROR(yac_interp_grid_get_aux_grid_src): invalid number of source fields")

  yac_dist_grid_pair_get_aux_grid(
    interp_grid->grid_pair, interp_grid->src_grid_name, cells, count,
    vertex_to_cell, vertex_to_cell_offsets, num_cells_per_vertex,
    interp_grid->src_fields[0]);
}

void yac_interp_grid_relocate_src_tgt_pairs(
  struct yac_interp_grid * interp_grid, int to_tgt_owner,
  size_t src_field_idx, size_t ** src_points,
  size_t ** tgt_points, double ** weights, size_t * count) {

  yac_dist_grid_pair_relocate_point_pairs(
    interp_grid->grid_pair, !to_tgt_owner, 1,
    interp_grid->src_grid_name, src_points,
    interp_grid->src_fields[src_field_idx].location,
    interp_grid->tgt_grid_name, tgt_points,
    interp_grid->tgt_field.location, weights, count);
}

void yac_interp_grid_determine_dist_tgt_owners(
  struct yac_interp_grid * interp_grid, size_t * tgt_indices, size_t count,
  int * owners) {

  yac_dist_grid_pair_determine_dist_owner(
    interp_grid->grid_pair, interp_grid->tgt_grid_name,
    tgt_indices, count, interp_grid->tgt_field.location, owners);
}

void yac_interp_grid_get_tgt_vertex_neighbours(
  struct yac_interp_grid * interp_grid, size_t * vertices, size_t count,
  size_t ** neigh_vertices, int * num_neighs_per_vertex) {

  yac_dist_grid_pair_get_vertex_neighbours(
    interp_grid->grid_pair, interp_grid->tgt_grid_name,
    vertices, count, neigh_vertices, num_neighs_per_vertex,
    interp_grid->tgt_field);
}

void yac_interp_grid_relocate_src_tgt_pairs_orig(
  struct yac_interp_grid * interp_grid, int to_tgt_owner,
  enum yac_location src_location, size_t ** src_points,
  size_t ** tgt_points, double ** weights, size_t * count) {

  yac_dist_grid_pair_relocate_point_pairs(
    interp_grid->grid_pair, !to_tgt_owner, 0,
    interp_grid->src_grid_name, src_points, src_location,
    interp_grid->tgt_grid_name, tgt_points,
    interp_grid->tgt_field.location, weights, count);
}

void yac_interp_grid_delete(struct yac_interp_grid * interp_grid) {

  if (interp_grid == NULL) return;

  free(interp_grid->src_grid_name);
  free(interp_grid->tgt_grid_name);
  free(interp_grid);
}
