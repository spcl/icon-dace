// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef GRID_REG2D_COMMON_H
#define GRID_REG2D_COMMON_H

#include "basic_grid_data.h"

struct yac_basic_grid_data yac_generate_basic_grid_data_reg2d_common(
  size_t nbr_vertices[2], int cyclic[2]);

#endif // GRID_REG2D_COMMON_H
