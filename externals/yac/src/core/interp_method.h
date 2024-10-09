// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_H
#define INTERP_METHOD_H

#include "interp_weights.h"
#include "interp_grid.h"
#include "yac_types.h"

// YAC PUBLIC HEADER START

struct interp_method;

struct yac_interp_weights * yac_interp_method_do_search(
  struct interp_method ** method, struct yac_interp_grid * grid);

void yac_interp_method_delete(struct interp_method ** method);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_H
