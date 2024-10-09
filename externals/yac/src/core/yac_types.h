// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_TYPES_H
#define YAC_TYPES_H

// YAC PUBLIC HEADER START

#include <stddef.h>

#include <yaxt.h>

// idxtype from yaxt as global id type in yac
typedef Xt_int yac_int;
#define yac_int_dt Xt_int_dt

// types for 3D coordinate arrays
typedef double(*yac_coordinate_pointer)[3];
typedef double const (* yac_const_coordinate_pointer)[3] ;

// type for storing an array for two indices
typedef size_t (* yac_size_t_2_pointer)[2];

// YAC PUBLIC HEADER STOP

#endif // YAC_TYPES_H
