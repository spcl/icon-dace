/**
 * @file interp_grid.c
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <yaxt.h>

#include "dist_grid.h"
#include "dist_grid_utils.h"
#include "interp_grid.h"
#include "geometry.h"
#include "yac_mpi.h"
#include "utils.h"
#include "sphere_part.h"
#include "yac_interface.h"

struct interp_grid {
  char const src_grid_name[YAC_MAX_CHARLEN];
  char const tgt_grid_name[YAC_MAX_CHARLEN];
  struct dist_grid_pair * grid_pair;
  struct interp_field tgt_field;
  size_t num_src_fields;
  struct interp_field src_fields[];
};

struct interp_grid * yac_interp_grid_new(
  struct dist_grid_pair * grid_pair,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t num_src_fields, struct interp_field const * src_fields,
  struct interp_field const tgt_field) {

  struct interp_grid * interp_grid =
    xmalloc(1 * sizeof(*interp_grid) + num_src_fields * sizeof(*src_fields));

  YAC_ASSERT_F(
    strlen(src_grid_name) < YAC_MAX_CHARLEN,
    "ERROR(yac_interp_grid_new): source grid name \"%s\"is too long "
    "(maximum length is YAC_MAX_CHARLEN)", src_grid_name)
  YAC_ASSERT_F(
    strlen(tgt_grid_name) < YAC_MAX_CHARLEN,
    "ERROR(yac_interp_grid_new): target grid name \"%s\"is too long "
    "(maximum length is YAC_MAX_CHARLEN)", tgt_grid_name)

  strncpy((char *)interp_grid->src_grid_name, src_grid_name, YAC_MAX_CHARLEN-1);
  strncpy((char *)interp_grid->tgt_grid_name, tgt_grid_name, YAC_MAX_CHARLEN-1);
  interp_grid->grid_pair = grid_pair;
  interp_grid->num_src_fields = num_src_fields;
  memcpy(&interp_grid->tgt_field, &tgt_field, 1 * sizeof(tgt_field));
  memcpy(
    interp_grid->src_fields, src_fields, num_src_fields * sizeof(*src_fields));

  return interp_grid;
}

void yac_interp_grid_get_src_points(
  struct interp_grid * interp_grid, size_t src_field_idx,
  size_t ** src_indices, size_t * count) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);
  yac_dist_grid_get_local_unmasked_points(
    dist_grid, interp_grid->src_fields[src_field_idx],
    src_indices, count);
}

void yac_interp_grid_get_tgt_points(
  struct interp_grid * interp_grid, size_t ** tgt_indices, size_t * count) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);
  yac_dist_grid_get_local_unmasked_points(
    dist_grid, interp_grid->tgt_field, tgt_indices, count);
}

struct remote_point * yac_interp_grid_get_src_remote_points2(
  struct interp_grid * interp_grid, enum yac_location location,
  size_t * src_points, size_t count) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  return
    yac_dist_grid_get_remote_points(dist_grid, location, src_points, count);
}

struct remote_point * yac_interp_grid_get_src_remote_points(
  struct interp_grid * interp_grid, size_t src_field_idx,
  size_t * src_points, size_t count) {

  return
    yac_interp_grid_get_src_remote_points2(
      interp_grid, interp_grid->src_fields[src_field_idx].location,
      src_points, count);
}

void yac_interp_grid_src_global_to_local(
  struct interp_grid * interp_grid, size_t src_field_idx,
  yac_int * src_global_ids, size_t count, size_t * src_local_ids) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  return
    yac_dist_grid_global_to_local(
      dist_grid, interp_grid->src_fields[src_field_idx].location,
      src_global_ids, count, src_local_ids);
}

void yac_interp_grid_tgt_global_to_local(
  struct interp_grid * interp_grid, yac_int * tgt_global_ids,
  size_t count, size_t * tgt_local_ids) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);

  return
    yac_dist_grid_global_to_local(
      dist_grid, interp_grid->tgt_field.location,
      tgt_global_ids, count, tgt_local_ids);
}

struct remote_point * yac_interp_grid_get_tgt_remote_points(
  struct interp_grid * interp_grid, size_t * tgt_points, size_t count) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);

  return
    yac_dist_grid_get_remote_points(
      dist_grid, interp_grid->tgt_field.location, tgt_points, count);
}

static yac_int const * get_tgt_grid_global_ids(
  struct interp_grid * interp_grid) {

  struct const_basic_grid_data * basic_grid_data =
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->tgt_grid_name));

  enum yac_location location = interp_grid->tgt_field.location;
  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(get_tgt_grid_global_ids): invalide target location")

  return basic_grid_data->ids[location];
}

enum yac_location yac_interp_grid_get_tgt_field_location(
  struct interp_grid * interp_grid) {

  return interp_grid->tgt_field.location;
}

void yac_interp_grid_get_src_global_ids(
  struct interp_grid * interp_grid, size_t * src_points, size_t count,
  size_t src_field_idx, yac_int * src_global_ids) {

  const_yac_int_pointer grid_global_ids =
    yac_interp_grid_get_src_field_global_ids(interp_grid, src_field_idx);

  for (size_t i = 0; i < count; ++i)
    src_global_ids[i] = grid_global_ids[src_points[i]];
}

void yac_interp_grid_get_tgt_global_ids(
  struct interp_grid * interp_grid, size_t * tgt_points, size_t count,
  yac_int * tgt_global_ids) {

  yac_int const * grid_global_ids = get_tgt_grid_global_ids(interp_grid);

  for (size_t i = 0; i < count; ++i)
    tgt_global_ids[i] = grid_global_ids[tgt_points[i]];
}

static void yac_interp_grid_get_coordinates(
  size_t * points, size_t count, coordinate_pointer coordinates,
  const_coordinate_pointer grid_coordinates) {

  YAC_ASSERT(
    (grid_coordinates != NULL) || (count == 0),
    "ERROR(yac_interp_grid_get_coordinates): grid_coordinates == NULL")

  for (size_t i = 0; i < count; ++i)
    for (int j = 0; j < 3; ++j)
      coordinates[i][j] = grid_coordinates[points[i]][j];
}

void yac_interp_grid_get_src_coordinates(
  struct interp_grid * interp_grid, size_t * src_points, size_t count,
  size_t src_field_idx, coordinate_pointer src_coordinates) {

  yac_interp_grid_get_coordinates(
    src_points, count, src_coordinates,
    yac_interp_grid_get_src_field_coords(interp_grid, src_field_idx));
}

void yac_interp_grid_get_tgt_coordinates(
  struct interp_grid * interp_grid, size_t * tgt_points, size_t count,
  coordinate_pointer tgt_coordinates) {

  yac_interp_grid_get_coordinates(
    tgt_points, count, tgt_coordinates,
    yac_interp_grid_get_tgt_field_coords(interp_grid));
}

size_t yac_interp_grid_get_num_src_fields(struct interp_grid * interp_grid) {

  return interp_grid->num_src_fields;
}

enum yac_location yac_interp_grid_get_src_field_location(
  struct interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_location): invalid src_field_idx")

  return interp_grid->src_fields[src_field_idx].location;
}

const_yac_int_pointer yac_interp_grid_get_src_field_global_ids(
  struct interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_global_ids): invalid src_field_idx")

  struct const_basic_grid_data * basic_grid_data =
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->src_grid_name));

  enum yac_location location =
    interp_grid->src_fields[src_field_idx].location;
  YAC_ASSERT(
    (location == CORNER) ||
    (location == CELL) ||
    (location == EDGE),
    "ERROR(yac_interp_grid_get_src_field_global_ids): "
    "invalid source field location")
  return basic_grid_data->ids[location];
}

const_coordinate_pointer yac_interp_grid_get_src_field_coords(
  struct interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_coords): invalid src_field_idx")

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  return
    (const_coordinate_pointer)
      yac_dist_grid_get_field_coords(
        dist_grid, interp_grid->src_fields[src_field_idx]);
}

const_coordinate_pointer yac_interp_grid_get_tgt_field_coords(
  struct interp_grid * interp_grid) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->tgt_grid_name);

  return
    (const_coordinate_pointer)
      yac_dist_grid_get_field_coords(dist_grid, interp_grid->tgt_field);
}

const_int_pointer yac_interp_grid_get_src_field_mask(
  struct interp_grid * interp_grid, size_t src_field_idx) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_get_src_field_mask): invalid src_field_idx")

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(
      interp_grid->grid_pair, interp_grid->src_grid_name);

  return
    yac_dist_grid_get_field_mask(
      dist_grid, interp_grid->src_fields[src_field_idx]);
}

void yac_interp_grid_do_points_search(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t * src_cells) {

  yac_dist_grid_pair_do_point_search(
    interp_grid->grid_pair, interp_grid->src_grid_name, search_coords, count,
    src_cells);
}

void yac_interp_grid_do_points_search_gc(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t * src_cells) {

  yac_dist_grid_pair_do_point_search_gc(
    interp_grid->grid_pair, interp_grid->src_grid_name, search_coords, count,
    src_cells);
}

void yac_interp_grid_do_nnn_search_src(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t n, size_t * src_points) {

  YAC_ASSERT(
    interp_grid->num_src_fields == 1,
    "ERROR(yac_interp_grid_do_nnn_search_src): invalid number of source fields")

  yac_dist_grid_pair_do_nnn_search(
    interp_grid->grid_pair, interp_grid->src_grid_name, search_coords, count,
    src_points, n, interp_grid->src_fields[0]);
}

void yac_interp_grid_do_nnn_search_tgt(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t n, size_t * tgt_points) {

  yac_dist_grid_pair_do_nnn_search(
    interp_grid->grid_pair, interp_grid->tgt_grid_name, search_coords, count,
    tgt_points, n, interp_grid->tgt_field);
}

void yac_interp_grid_do_bnd_circle_search_src(
  struct interp_grid * interp_grid, const_bounding_circle_pointer bnd_circles,
  size_t count, size_t src_field_idx, size_t ** src_cells,
  size_t * num_src_per_bnd_circle) {

  YAC_ASSERT(
    src_field_idx < interp_grid->num_src_fields,
    "ERROR(yac_interp_grid_do_bnd_circle_search_src): invalid src_field_idx")

  YAC_ASSERT(
    interp_grid->src_fields[src_field_idx].location == CELL,
    "ERROR(yac_interp_grid_do_bnd_circle_search_src): "
    "invalid source field location; has to be CELL")

  yac_dist_grid_pair_do_bnd_circle_search(
    interp_grid->grid_pair, interp_grid->src_grid_name, bnd_circles, count,
    src_cells, num_src_per_bnd_circle, interp_grid->src_fields[src_field_idx]);
}

void yac_interp_grid_do_bnd_circle_search_tgt(
  struct interp_grid * interp_grid, const_bounding_circle_pointer bnd_circles,
  size_t count, size_t ** tgt_cells, size_t * num_tgt_per_bnd_circle) {

  YAC_ASSERT(
    interp_grid->tgt_field.location == CELL,
    "ERROR(yac_interp_grid_do_bnd_circle_search_tgt): "
    "invalid target field location; has to be CELL")

  yac_dist_grid_pair_do_bnd_circle_search(
    interp_grid->grid_pair, interp_grid->tgt_grid_name, bnd_circles, count,
    tgt_cells, num_tgt_per_bnd_circle, interp_grid->tgt_field);
}

void yac_interp_grid_do_cell_search_src(
  struct interp_grid * interp_grid, size_t * tgt_cells, size_t count,
  size_t ** src_cells, size_t * num_src_per_tgt) {

  YAC_ASSERT(
    interp_grid->num_src_fields == 1,
    "ERROR(yac_interp_grid_do_cell_search_src): "
    "invalid number of source fields")

  YAC_ASSERT(
    interp_grid->src_fields[0].location == CELL,
    "ERROR(yac_interp_grid_do_cell_search_src): "
    "invalid source field location; has to be CELL")

  YAC_ASSERT(
    interp_grid->tgt_field.location == CELL,
    "ERROR(yac_interp_grid_do_cell_search_src): "
    "invalid target field location; has to be CELL")

  yac_dist_grid_pair_do_cell_search(
    interp_grid->grid_pair, interp_grid->tgt_grid_name,
    interp_grid->src_grid_name, tgt_cells, count, src_cells,
    num_src_per_tgt, interp_grid->src_fields[0]);
}

void yac_interp_grid_do_cell_search_tgt(
  struct interp_grid * interp_grid, size_t * src_cells, size_t count,
  size_t ** tgt_cells, size_t * num_tgt_per_src) {

  YAC_ASSERT(
    interp_grid->num_src_fields == 1,
    "ERROR(yac_interp_grid_do_cell_search_tgt): "
    "invalid number of source fields")

  YAC_ASSERT(
    interp_grid->src_fields[0].location == CELL,
    "ERROR(yac_interp_grid_do_cell_search_tgt): "
    "invalid source field location; has to be CELL")

  YAC_ASSERT(
    interp_grid->tgt_field.location == CELL,
    "ERROR(yac_interp_grid_do_cell_search_tgt): "
    "invalid target field location; has to be CELL")

  yac_dist_grid_pair_do_cell_search(
    interp_grid->grid_pair, interp_grid->src_grid_name,
    interp_grid->tgt_grid_name, src_cells, count, tgt_cells,
    num_tgt_per_src, interp_grid->tgt_field);
}

MPI_Comm yac_interp_grid_get_MPI_Comm(struct interp_grid * interp_grid) {

  return yac_dist_grid_pair_get_MPI_Comm(interp_grid->grid_pair);
}

struct const_basic_grid_data * yac_interp_grid_get_basic_grid_data_src(
  struct interp_grid * interp_grid) {

  return
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->src_grid_name));
}

struct const_basic_grid_data * yac_interp_grid_get_basic_grid_data_tgt(
  struct interp_grid * interp_grid) {

  return
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(
        interp_grid->grid_pair, interp_grid->tgt_grid_name));
}

void yac_interp_grid_get_src_cell_neighbours(
  struct interp_grid * interp_grid, size_t * src_cells, size_t count,
  size_t * neighbours) {

  yac_dist_grid_pair_get_cell_neighbours(
    interp_grid->grid_pair, interp_grid->src_grid_name, src_cells, count,
    neighbours);
}

void yac_interp_grid_get_tgt_cell_neighbours(
  struct interp_grid * interp_grid, size_t * tgt_cells, size_t count,
  size_t * neighbours) {

  yac_dist_grid_pair_get_cell_neighbours(
    interp_grid->grid_pair, interp_grid->tgt_grid_name, tgt_cells, count,
    neighbours);
}

void yac_interp_grid_get_aux_grid_src(
  struct interp_grid * interp_grid, size_t * cells, size_t count,
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
  struct interp_grid * interp_grid, int to_tgt_owner,
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
  struct interp_grid * interp_grid, size_t * tgt_indices, size_t count,
  int * owners) {

  yac_dist_grid_pair_determine_dist_owner(
    interp_grid->grid_pair, interp_grid->tgt_grid_name,
    tgt_indices, count, interp_grid->tgt_field.location, owners);
}

void yac_interp_grid_get_tgt_vertex_neighbours(
  struct interp_grid * interp_grid, size_t * vertices, size_t count,
  size_t ** neigh_vertices, int * num_neighs_per_vertex) {

  yac_dist_grid_pair_get_vertex_neighbours(
    interp_grid->grid_pair, interp_grid->tgt_grid_name,
    vertices, count, neigh_vertices, num_neighs_per_vertex,
    interp_grid->tgt_field);
}

void yac_interp_grid_relocate_src_tgt_pairs_orig(
  struct interp_grid * interp_grid, int to_tgt_owner,
  enum yac_location src_location, size_t ** src_points,
  size_t ** tgt_points, double ** weights, size_t * count) {

  yac_dist_grid_pair_relocate_point_pairs(
    interp_grid->grid_pair, !to_tgt_owner, 0,
    interp_grid->src_grid_name, src_points, src_location,
    interp_grid->tgt_grid_name, tgt_points,
    interp_grid->tgt_field.location, weights, count);
}

void yac_interp_grid_delete(struct interp_grid * interp_grid) {

  if (interp_grid == NULL) return;

  free(interp_grid);
}
