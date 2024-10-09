// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "utils_core.h"

int yac_file_exists(const char * filename) {
  struct stat buffer;
  return !stat(filename,&buffer);
}

char const * yac_name_type_pair_get_name(
  struct yac_name_type_pair const * pairs, size_t count, int type) {

  char const * name = NULL;
  for (size_t i = 0; (i < count) && (name == NULL); ++i)
    if (pairs[i].type == type) name = pairs[i].name;
  return name;
}

int yac_name_type_pair_get_type(
  struct yac_name_type_pair const * pairs, size_t count, char const * name) {

  int type = INT_MAX;
  if (name != NULL)
    for (size_t i = 0; (i < count) && (type == INT_MAX); ++i)
      if (!strcmp(pairs[i].name, name)) type = pairs[i].type;
  return type;
}
