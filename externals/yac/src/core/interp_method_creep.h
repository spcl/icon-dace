// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_CREEP_H
#define INTERP_METHOD_CREEP_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

/** \example test_interp_method_creep_parallel.c
 * A test for the parallel creep interpolation method.
 */

#define YAC_INTERP_CREEP_DISTANCE_DEFAULT (-1)

struct interp_method * yac_interp_method_creep_new(int creep_distance);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_CREEP_H
