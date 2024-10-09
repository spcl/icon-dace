// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_FIXED_H
#define INTERP_METHOD_FIXED_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

#define YAC_INTERP_FIXED_VALUE_DEFAULT (DBL_MAX)

struct interp_method * yac_interp_method_fixed_new(double value);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_FIXED_H
