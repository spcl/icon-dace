// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIELD_DATA_SET_H
#define FIELD_DATA_SET_H

#include "basic_grid.h"
#include "utils_common.h"
#include "location.h"
#include "field_data.h"

struct yac_field_data_set;

struct yac_field_data_set * yac_field_data_set_empty_new();
struct yac_field_data_set *
  yac_field_data_set_new(
    struct yac_field_data * cell_field_data,
    struct yac_field_data * vertex_field_data,
    struct yac_field_data * edge_field_data);
size_t yac_field_data_set_add_mask(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask, size_t count,
  char const * mask_name);
size_t yac_field_data_set_add_coordinates(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, yac_coordinate_pointer coordinates,
  size_t count);
size_t yac_field_data_set_add_mask_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask, char const * mask_name);
size_t yac_field_data_set_add_coordinates_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, yac_coordinate_pointer coordinates);
struct yac_field_data * yac_field_data_set_get_field_data(
  struct yac_field_data_set * field_data_set, enum yac_location location);
void yac_field_data_set_delete(struct yac_field_data_set * field_data_set);

#endif // FIELD_DATA_SET_H

