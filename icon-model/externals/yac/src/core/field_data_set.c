// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "ppm/core.h"
#include "field_data_set.h"

struct yac_field_data_set {
  struct yac_field_data * cell;
  struct yac_field_data * vertex;
  struct yac_field_data * edge;
};

struct yac_field_data_set *
  yac_field_data_set_new(
    struct yac_field_data * cell_field_data,
    struct yac_field_data * vertex_field_data,
    struct yac_field_data * edge_field_data) {

  struct yac_field_data_set * field_data_set =
    xmalloc(1 * sizeof(*field_data_set));

  field_data_set->cell = cell_field_data;
  field_data_set->vertex = vertex_field_data;
  field_data_set->edge = edge_field_data;

  return field_data_set;
}

struct yac_field_data_set * yac_field_data_set_empty_new() {

  return
    yac_field_data_set_new(
      yac_field_data_empty_new(),
      yac_field_data_empty_new(),
      yac_field_data_empty_new());
}

struct yac_field_data * yac_field_data_set_get_field_data(
  struct yac_field_data_set * field_data_set, enum yac_location location) {

  YAC_ASSERT(
    (location == YAC_LOC_CELL) ||
    (location == YAC_LOC_CORNER) ||
    (location == YAC_LOC_EDGE),
    "ERROR(yac_field_data_set_get_field_data): invalid location")

  switch (location) {
    default:
    case (YAC_LOC_CELL): return field_data_set->cell;
    case (YAC_LOC_CORNER): return field_data_set->vertex;
    case (YAC_LOC_EDGE): return field_data_set->edge;
  };
}

size_t yac_field_data_set_add_mask_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask,
  char const * mask_name) {

  return
    yac_field_data_add_mask_nocpy(
      yac_field_data_set_get_field_data(
        field_data_set, location), mask, mask_name);
}

size_t yac_field_data_set_add_mask(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask,
  size_t count, char const * mask_name) {

  int * mask_cpy = xmalloc(count * sizeof(*mask_cpy));
  memcpy(mask_cpy, mask, count * sizeof(*mask));

  char * mask_name_cpy =
    (mask_name != NULL)?strdup(mask_name):NULL;

  return
    yac_field_data_set_add_mask_nocpy(
      field_data_set, location, mask_cpy,
      mask_name_cpy);
}

size_t yac_field_data_set_add_coordinates_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, yac_coordinate_pointer coordinates) {

  return
    yac_field_data_add_coordinates_nocpy(
      yac_field_data_set_get_field_data(
        field_data_set, location), coordinates);
}

size_t yac_field_data_set_add_coordinates(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, yac_coordinate_pointer coordinates,
  size_t count) {

  yac_coordinate_pointer coordinates_cpy =
    xmalloc(count * sizeof(*coordinates_cpy));
  memcpy(coordinates_cpy, coordinates, count * sizeof(*coordinates));

  return
    yac_field_data_set_add_coordinates_nocpy(
      field_data_set, location, coordinates_cpy);
}

void yac_field_data_set_delete(struct yac_field_data_set * field_data_set) {
  yac_field_data_delete(field_data_set->cell);
  yac_field_data_delete(field_data_set->vertex);
  yac_field_data_delete(field_data_set->edge);
  free(field_data_set);
}
