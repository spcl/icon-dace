// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include "yac_types.h"

// YAC PUBLIC HEADER START

struct yac_field_data;

struct yac_field_data * yac_field_data_empty_new();
size_t yac_field_data_add_mask_nocpy(
  struct yac_field_data * field_data, int const * mask,
  char const * mask_name);
size_t yac_field_data_add_coordinates_nocpy(
  struct yac_field_data * field_data, yac_coordinate_pointer coordinates);
size_t yac_field_data_get_masks_count(struct yac_field_data * field_data);
int const * yac_field_data_get_mask_data(
  struct yac_field_data * field_data, size_t mask_idx);
void yac_field_data_set_mask_data(
  struct yac_field_data * field_data, size_t mask_idx, int * mask_data);
char const * yac_field_data_get_mask_name(
  struct yac_field_data * field_data, size_t mask_idx);
size_t yac_field_data_get_coordinates_count(struct yac_field_data * field_data);
yac_const_coordinate_pointer yac_field_data_get_coordinates_data(
  struct yac_field_data * field_data, size_t coordinates_idx);
void yac_field_data_set_coordinates_data(
  struct yac_field_data * field_data, size_t coordinates_idx,
  yac_coordinate_pointer coordinates_data);
void yac_field_data_delete(struct yac_field_data * field_data);

// YAC PUBLIC HEADER STOP

#endif // FIELD_DATA_H
