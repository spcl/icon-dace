// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>

#include "field_data.h"
#include "utils_common.h"

struct yac_field_data {
  struct {
    char * name;
    int * data;
  } * masks;
  size_t masks_count;
  yac_coordinate_pointer * coordinates;
  size_t coordinates_count;
};

struct yac_field_data * yac_field_data_empty_new() {

  struct yac_field_data * field_data = xmalloc(1 * sizeof(*field_data));
  field_data->masks_count = 0;
  field_data->masks = NULL;
  field_data->coordinates_count = 0;
  field_data->coordinates = NULL;
  return field_data;
}

size_t yac_field_data_add_mask_nocpy(
  struct yac_field_data * field_data, int const * mask,
  char const * mask_name) {

  size_t masks_idx = field_data->masks_count++;
  field_data->masks =
    xrealloc(field_data->masks,
             field_data->masks_count * sizeof(*(field_data->masks)));
  field_data->masks[masks_idx].name = (char *)mask_name;
  field_data->masks[masks_idx].data = (int*)mask;

  return masks_idx;
}

size_t yac_field_data_add_coordinates_nocpy(
  struct yac_field_data * field_data, yac_coordinate_pointer coordinates) {

  size_t coordinates_idx = field_data->coordinates_count++;
  field_data->coordinates =
    xrealloc(
      field_data->coordinates,
      field_data->coordinates_count * sizeof(*(field_data->coordinates)));
  field_data->coordinates[coordinates_idx] = coordinates;

  return coordinates_idx;
}

size_t yac_field_data_get_masks_count(struct yac_field_data * field_data) {

  return field_data->masks_count;
}

int const * yac_field_data_get_mask_data(
  struct yac_field_data * field_data, size_t mask_idx) {

  YAC_ASSERT(
    mask_idx < field_data->masks_count,
    "ERROR(yac_field_data_get_mask_data): invalid mask index");

  return field_data->masks[mask_idx].data;
}

void yac_field_data_set_mask_data(
  struct yac_field_data * field_data, size_t mask_idx, int * mask_data) {

  YAC_ASSERT(
    mask_idx < field_data->masks_count,
    "ERROR(yac_field_data_set_mask_data): invalid mask index");

  field_data->masks[mask_idx].data = mask_data;
}

char const * yac_field_data_get_mask_name(
  struct yac_field_data * field_data, size_t mask_idx) {

  YAC_ASSERT(
    mask_idx < field_data->masks_count,
    "ERROR(yac_field_data_get_mask_name): invalid mask index");

  return field_data->masks[mask_idx].name;
}

size_t yac_field_data_get_coordinates_count(struct yac_field_data * field_data) {

  return field_data->coordinates_count;
}

yac_const_coordinate_pointer yac_field_data_get_coordinates_data(
  struct yac_field_data * field_data, size_t coordinates_idx) {

  YAC_ASSERT(
    coordinates_idx < field_data->coordinates_count,
    "ERROR(yac_field_data_get_coordinates_data): invalid coordinates index");

  return
    (yac_const_coordinate_pointer)(field_data->coordinates[coordinates_idx]);
}

void yac_field_data_set_coordinates_data(
  struct yac_field_data * field_data, size_t coordinates_idx,
  yac_coordinate_pointer coordinates_data) {

  YAC_ASSERT(
    coordinates_idx < field_data->coordinates_count,
    "ERROR(yac_field_data_set_coordinates_data): invalid coordinates index");

  field_data->coordinates[coordinates_idx] = coordinates_data;
}

void yac_field_data_delete(struct yac_field_data * field_data) {

  for (size_t i = 0; i < field_data->masks_count; ++i) {
    free(field_data->masks[i].name);
    free(field_data->masks[i].data);
  }
  free(field_data->masks);
  for (size_t i = 0; i < field_data->coordinates_count; ++i)
    free(field_data->coordinates[i]);
  free(field_data->coordinates);
  free(field_data);
}
