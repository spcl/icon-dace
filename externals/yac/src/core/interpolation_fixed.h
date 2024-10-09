// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERPOLATION_FIXED_H
#define INTERPOLATION_FIXED_H

#include "interpolation_internal.h"

struct yac_interpolation_type * yac_interpolation_fixed_new(
  size_t collection_size, double value, size_t count, size_t const * pos);

#endif // INTERPOLATION_FIXED_H
