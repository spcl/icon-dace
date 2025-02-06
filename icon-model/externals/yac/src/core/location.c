// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <string.h>

#include "location.h"
#include "utils_common.h"

struct yac_name_type_pair
  yac_location_name_type_pair[] =
  {{.name = "CELL",      .type = YAC_LOC_CELL},
   {.name = "CORNER",    .type = YAC_LOC_CORNER},
   {.name = "EDGE",      .type = YAC_LOC_EDGE},
   {.name = "UNDEFINED", .type = YAC_LOC_UNDEFINED}};
enum  {YAC_LOCATION_COUNT =
         sizeof(yac_location_name_type_pair) /
         sizeof(yac_location_name_type_pair[0])};

enum yac_location yac_str2loc(char const * location_str) {

  int location =
    yac_name_type_pair_get_type(
      yac_location_name_type_pair, YAC_LOCATION_COUNT, location_str);

  YAC_ASSERT(location != INT_MAX, "ERROR(yac_str2loc): invalid location")

  return (enum yac_location)location;
}

char const * yac_loc2str(enum yac_location location) {

  char const * location_str =
    yac_name_type_pair_get_name(
      yac_location_name_type_pair, YAC_LOCATION_COUNT, location);
  YAC_ASSERT(location_str,
    "ERROR(yac_loc2str): location must be one of "
    "YAC_LOC_CORNER/YAC_LOC_EDGE/YAC_LOC_CELL/YAC_LOC_UNDEFINED.")

  return location_str;
}

enum yac_location yac_get_location(int const location) {

  YAC_ASSERT(
    (location == YAC_LOC_CELL) || (location == YAC_LOC_CORNER) ||
    (location == YAC_LOC_EDGE) || (location == YAC_LOC_UNDEFINED),
    "ERROR(get_location): location must be one of "
    "YAC_LOC_CORNER/YAC_LOC_EDGE/YAC_LOC_CELL/YAC_LOC_UNDEFINED.")

  return (enum yac_location) location;
}
