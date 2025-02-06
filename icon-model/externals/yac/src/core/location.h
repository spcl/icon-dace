// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LOCATION_H
#define LOCATION_H

// YAC PUBLIC HEADER START

#define YAC_MAX_LOC_STR_LEN 10

enum yac_location {

   YAC_LOC_CELL =   0,
   YAC_LOC_CORNER = 1,
   YAC_LOC_EDGE =   2,
   YAC_LOC_UNDEFINED = 3,
};

enum yac_location yac_str2loc(char const * location);
char const * yac_loc2str(enum yac_location location);
enum yac_location yac_get_location(int const location);

// YAC PUBLIC HEADER STOP

#endif // LOCATION_H

