/**
 * @file grid.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <stdio.h>
#include <string.h>

#include "grid.h"
#include "grid_cell.h"
#include "yac_interface.h"

struct yac_name_type_pair
  yac_location_name_type_pair[] =
  {{.name = "CELL",      .type = CELL},
   {.name = "CORNER",    .type = CORNER},
   {.name = "EDGE",      .type = EDGE},
   {.name = "UNDEFINED", .type = LOC_UNDEFINED}};
enum  {YAC_LOCATION_COUNT =
         sizeof(yac_location_name_type_pair) /
         sizeof(yac_location_name_type_pair[0])};

struct yac_basic_grid {
  char name[YAC_MAX_CHARLEN];
  int is_empty;
  struct yac_field_data_set field_data_set;
  struct basic_grid_data data;
};

static struct basic_grid_data basic_grid_data_empty = {
  .vertex_coordinates = NULL,
  .cell_ids = NULL,
  .vertex_ids = NULL,
  .edge_ids = NULL,
  .num_cells = 0,
  .num_vertices = 0,
  .num_edges = 0,
  .core_cell_mask = NULL,
  .core_vertex_mask = NULL,
  .core_edge_mask = NULL,
  .num_vertices_per_cell = NULL,
  .num_cells_per_vertex = NULL,
  .cell_to_vertex = NULL,
  .cell_to_vertex_offsets = NULL,
  .cell_to_edge = NULL,
  .cell_to_edge_offsets = NULL,
  .vertex_to_cell = NULL,
  .vertex_to_cell_offsets = NULL,
  .edge_to_vertex = NULL,
  .edge_type = NULL,
  .num_total_cells = 0,
  .num_total_vertices = 0,
  .num_total_edges = 0
};

struct yac_basic_grid * yac_basic_grid_new(
  char const * name, struct basic_grid_data grid_data) {

  struct yac_basic_grid * grid = xmalloc(1 * sizeof(*grid));

  YAC_ASSERT_F(
    strlen(name) < YAC_MAX_CHARLEN,
    "ERROR(yac_basic_grid_new): grid name \"%s\"is too long "
    "(maximum length is YAC_MAX_CHARLEN)", name)

  strncpy(grid->name, name, YAC_MAX_CHARLEN-1);
  grid->is_empty = 0;
  grid->field_data_set = yac_field_data_set_init();
  grid->data = grid_data;

  return grid;
};

struct yac_basic_grid * yac_basic_grid_empty_new(char const * name) {
  struct yac_basic_grid * grid =
    yac_basic_grid_new(name, basic_grid_data_empty);
  grid->is_empty = 1;
  return grid;
}

void yac_basic_grid_delete(struct yac_basic_grid * grid) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_delete): "
    "NULL is not a valid value for argument grid")
  yac_field_data_set_free(grid->field_data_set);
  yac_basic_grid_data_free(grid->data);
  free(grid);
}

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
    "CORNER/EDGE/CELL/LOC_UNDEFINED.")

  return location_str;
}

enum yac_location yac_get_location(int const location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) ||
    (location == EDGE) || (location == LOC_UNDEFINED),
    "ERROR(get_location): location must be one of "
    "CORNER/EDGE/CELL/LOC_UNDEFINED.")

  return (enum yac_location) location;
}

coordinate_pointer yac_basic_grid_get_field_coordinates(
  struct yac_basic_grid * grid, struct interp_field field) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_field_coordinates): "
    "NULL is not a valid value for argument grid")

  coordinate_pointer coords = NULL;

  if (field.coordinates_idx != SIZE_MAX) {

    struct yac_field_data * data =
      yac_basic_grid_get_field_data(grid, field.location);

    YAC_ASSERT(
      field.coordinates_idx < data->coordinates_count,
      "ERROR(yac_basic_grid_get_field_coordinates): invalid coord_idx")

    coords = data->coordinates[field.coordinates_idx];
  }

  // if no field coordinates are defined, but the location is at the corners of
  // of the grid cells, return coordinates of them
  if ((coords == NULL) && (field.location == CORNER))
    coords = grid->data.vertex_coordinates;

  return coords;
}

int const * yac_basic_grid_get_core_mask(
  struct yac_basic_grid * grid, enum yac_location location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_basic_grid_get_core_mask): invalid location")

  switch (location) {
    default:
    case(CELL): return grid->data.core_cell_mask;
    case(CORNER): return grid->data.core_vertex_mask;
    case(EDGE): return grid->data.core_edge_mask;
  };
}

int const * yac_basic_grid_get_field_mask(
  struct yac_basic_grid * grid, struct interp_field field) {

  if (field.masks_idx == SIZE_MAX) return NULL;

  struct yac_field_data * data =
    yac_basic_grid_get_field_data(grid, field.location);

  YAC_ASSERT(
    field.masks_idx < data->masks_count,
    "ERROR(yac_basic_grid_get_field_mask): invalid mask_idx")

  return data->masks[field.masks_idx].data;
}

char const * yac_basic_grid_get_name(struct yac_basic_grid * grid) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_name): "
    "NULL is not a valid value for argument grid")

  return grid->name;
}

struct basic_grid_data * yac_basic_grid_get_data(struct yac_basic_grid * grid) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_data): "
    "NULL is not a valid value for argument grid")

  return &(grid->data);
}

static struct yac_field_data * get_field_data(
  struct yac_field_data_set * field_data_set, enum yac_location location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(get_field_data): invalid location")

  switch (location) {
    default:
    case (CELL): return &field_data_set->cell;
    case (CORNER): return &field_data_set->vertex;
    case (EDGE): return &field_data_set->edge;
  };
}

size_t yac_basic_grid_get_data_size(
  struct yac_basic_grid * grid, enum yac_location location) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_data_size): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_basic_grid_get_data_size): invalid location")

  switch (location) {
    default:
    case (CELL):
      return grid->data.num_cells;
    case (CORNER):
      return grid->data.num_vertices;
    case (EDGE):
      return grid->data.num_edges;
  };
}

size_t yac_basic_grid_get_named_mask_idx(
  struct yac_basic_grid * grid, enum yac_location location,
  char const * mask_name) {

  if (mask_name == NULL) return SIZE_MAX;

  struct yac_field_data * data =
    yac_basic_grid_get_field_data(grid, location);

  if (data == NULL) return SIZE_MAX;

  size_t mask_idx = SIZE_MAX;

  for (size_t i = 0; (i < data->masks_count) && (mask_idx == SIZE_MAX); ++i)
    if ((data->masks[i].name != NULL) &&
        (!strcmp(mask_name, data->masks[i].name)))
      mask_idx = i;

  YAC_ASSERT_F(
    mask_idx != SIZE_MAX,
    "ERROR(yac_basic_grid_get_named_mask_idx): grid \"%s\" does not contain "
    "%s-mask with the name \"%s\"",
    grid->name, yac_loc2str(location), mask_name)

  return mask_idx;
}

size_t yac_basic_grid_add_coordinates_nocpy(
  struct yac_basic_grid * grid,
  enum yac_location location, coordinate_pointer coordinates) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_coordinates_nocpy): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_coordinates_nocpy): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_coordinates_nocpy(
      &grid->field_data_set, location, coordinates);
}

size_t yac_basic_grid_add_coordinates(
  struct yac_basic_grid * grid, enum yac_location location,
  coordinate_pointer coordinates, size_t count) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_coordinates): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_coordinates): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_coordinates(
      &grid->field_data_set, location, coordinates, count);
}

size_t yac_basic_grid_add_mask_nocpy(
  struct yac_basic_grid * grid,
  enum yac_location location, int const * mask,
  char const * mask_name) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_mask_nocpy): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_mask_nocpy): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_mask_nocpy(
      &grid->field_data_set, location, mask, mask_name);
}

size_t yac_basic_grid_add_mask(
  struct yac_basic_grid * grid, enum yac_location location,
  int const * mask, size_t count, char const * mask_name) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_mask): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_mask): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_mask(
      &grid->field_data_set, location, mask, count, mask_name);
}

struct yac_field_data * yac_basic_grid_get_field_data(
  struct yac_basic_grid * grid, enum yac_location location) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_field_data): "
    "NULL is not a valid value for argument grid")

  if (grid->is_empty)
    return NULL;
  else
    return get_field_data(&(grid->field_data_set), location);
}

void yac_basic_grid_data_free(struct basic_grid_data grid) {

  free(grid.vertex_coordinates);
  free(grid.cell_ids);
  free(grid.vertex_ids);
  free(grid.edge_ids);
  free(grid.core_cell_mask);
  free(grid.core_vertex_mask);
  free(grid.core_edge_mask);
  free(grid.num_vertices_per_cell);
  free(grid.num_cells_per_vertex);
  free(grid.cell_to_vertex);
  free(grid.cell_to_vertex_offsets);
  free(grid.cell_to_edge);
  if (grid.cell_to_vertex_offsets != grid.cell_to_edge_offsets)
    free(grid.cell_to_edge_offsets);
  free(grid.vertex_to_cell);
  free(grid.vertex_to_cell_offsets);
  free(grid.edge_to_vertex);
  free(grid.edge_type);
}

struct yac_field_data yac_field_data_init() {
  struct yac_field_data field_data = {.masks_count = 0,
                                      .masks = NULL,
                                      .coordinates_count = 0,
                                      .coordinates = NULL};
  return field_data;
}

size_t yac_field_data_add_mask_nocpy(
  struct yac_field_data * field_data, int const * mask,
  char const * mask_name) {

  size_t masks_idx = field_data->masks_count++;
  field_data->masks =
    xrealloc(field_data->masks,
             field_data->masks_count * sizeof(*(field_data->masks)));
  field_data->masks[masks_idx].name = (char *)mask_name;
  field_data->masks[masks_idx].data = (int*)mask;

  return masks_idx;
}

size_t yac_field_data_add_coordinates_nocpy(
  struct yac_field_data * field_data, coordinate_pointer coordinates) {

  size_t coordinates_idx = field_data->coordinates_count++;
  field_data->coordinates =
    xrealloc(
      field_data->coordinates,
      field_data->coordinates_count * sizeof(*(field_data->coordinates)));
  field_data->coordinates[coordinates_idx] = coordinates;

  return coordinates_idx;
}

struct yac_field_data_set yac_field_data_set_init() {
  struct yac_field_data_set field_data_set =
    {.cell = yac_field_data_init(),
     .vertex = yac_field_data_init(),
     .edge = yac_field_data_init()};
  return field_data_set;
}

size_t yac_field_data_set_add_mask_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask,
  char const * mask_name) {

  return
    yac_field_data_add_mask_nocpy(
      get_field_data(
        field_data_set, location), mask, mask_name);
}

size_t yac_field_data_set_add_mask(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask,
  size_t count, char const * mask_name) {

  int * mask_cpy = xmalloc(count * sizeof(*mask_cpy));
  memcpy(mask_cpy, mask, count * sizeof(*mask));

  char * mask_name_cpy =
    (mask_name != NULL)?strdup(mask_name):NULL;

  return
    yac_field_data_set_add_mask_nocpy(
      field_data_set, location, mask_cpy,
      mask_name_cpy);
}

size_t yac_field_data_set_add_coordinates_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, coordinate_pointer coordinates) {

  return
    yac_field_data_add_coordinates_nocpy(
      get_field_data(field_data_set, location), coordinates);
}

size_t yac_field_data_set_add_coordinates(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, coordinate_pointer coordinates,
  size_t count) {

  coordinate_pointer coordinates_cpy =
    xmalloc(count * sizeof(*coordinates_cpy));
  memcpy(coordinates_cpy, coordinates, count * sizeof(*coordinates));

  return
    yac_field_data_set_add_coordinates_nocpy(
      field_data_set, location, coordinates_cpy);
}

static void yac_field_data_free(struct yac_field_data field_data) {

  for (size_t i = 0; i < field_data.masks_count; ++i) {
    free(field_data.masks[i].name);
    free(field_data.masks[i].data);
  }
  free(field_data.masks);
  for (size_t i = 0; i < field_data.coordinates_count; ++i)
    free(field_data.coordinates[i]);
  free(field_data.coordinates);
}

void yac_field_data_set_free(struct yac_field_data_set field_data_set) {
  yac_field_data_free(field_data_set.cell);
  yac_field_data_free(field_data_set.vertex);
  yac_field_data_free(field_data_set.edge);
}
