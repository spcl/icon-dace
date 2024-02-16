/**
 * @file grid.h
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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

#ifndef GRID_H
#define GRID_H

#include "utils.h"
#include "grid_cell.h"

#define YAC_MAX_LOC_STR_LEN 10

enum yac_location {

   CELL =   0,
   CORNER = 1,
   EDGE =   2,
   LOC_UNDEFINED = 3,
};

/** \example test_grid.c
 */

#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif

#ifndef  M_PI_2
#define  M_PI_2      1.57079632679489661923132169163975144  /* pi/2 */
#endif

#define YAC_EARTH_RADIUS (6371.2290)
#define YAC_EARTH_RADIUS2 ((6371.2290 * 6371.2290) * 0.5)
#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

typedef size_t (* size_t_2_pointer)[2];

struct yac_field_data {
  struct {
    char * name;
    int * data;
  } * masks;
  size_t masks_count;
  coordinate_pointer * coordinates;
  size_t coordinates_count;
};

struct yac_field_data_set {
  struct yac_field_data cell, vertex, edge;
};

struct interp_field {
  enum yac_location location;
  size_t coordinates_idx;
  size_t masks_idx;
};

struct basic_grid_data {
  coordinate_pointer vertex_coordinates;
  yac_int * cell_ids;
  yac_int * vertex_ids;
  yac_int * edge_ids;
  size_t num_cells; // number of local cells (owned by local process)
  size_t num_vertices; // number of local vertices (owned by local process)
  size_t num_edges; // number of local edges (owned by local process)
  int * core_cell_mask;
  int * core_vertex_mask;
  int * core_edge_mask;
  int * num_vertices_per_cell;
  int * num_cells_per_vertex;
  size_t * cell_to_vertex;
  size_t * cell_to_vertex_offsets;
  size_t * cell_to_edge;
  size_t * cell_to_edge_offsets;
  size_t * vertex_to_cell;
  size_t * vertex_to_cell_offsets;
  size_t_2_pointer edge_to_vertex;
  enum yac_edge_type * edge_type;
  size_t num_total_cells; // number of locally stored cells
  size_t num_total_vertices; // number of locally stored vertices
  size_t num_total_edges; // number of locally stored edges
};

struct yac_basic_grid;

struct basic_grid_data yac_generate_basic_grid_data_reg_2d(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct basic_grid_data yac_generate_basic_grid_data_reg_2d_deg(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct basic_grid_data yac_generate_basic_grid_data_curve_2d(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct basic_grid_data yac_generate_basic_grid_data_curve_2d_deg(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices);

struct basic_grid_data yac_generate_basic_grid_data_unstruct(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct basic_grid_data yac_generate_basic_grid_data_unstruct_deg(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct basic_grid_data yac_generate_basic_grid_data_unstruct_ll(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct basic_grid_data yac_generate_basic_grid_data_unstruct_ll_deg(
  size_t nbr_vertices, size_t nbr_cells, int *num_vertices_per_cell,
  double *x_vertices, double *y_vertices, int *cell_to_vertex);

struct yac_basic_grid * yac_basic_grid_new(
  char const * name, struct basic_grid_data grid_data);
struct yac_basic_grid * yac_basic_grid_empty_new(char const * name);
coordinate_pointer yac_basic_grid_get_field_coordinates(
  struct yac_basic_grid * grid, struct interp_field field);
int const * yac_basic_grid_get_field_mask(
  struct yac_basic_grid * grid, struct interp_field field);
int const * yac_basic_grid_get_core_mask(
  struct yac_basic_grid * grid, enum yac_location location);
char const * yac_basic_grid_get_name(struct yac_basic_grid * grid);
struct basic_grid_data * yac_basic_grid_get_data(struct yac_basic_grid * grid);
struct yac_field_data * yac_basic_grid_get_field_data(
  struct yac_basic_grid * grid, enum yac_location location);
size_t yac_basic_grid_get_data_size(
  struct yac_basic_grid * grid, enum yac_location location);
size_t yac_basic_grid_get_named_mask_idx(
  struct yac_basic_grid * grid, enum yac_location location,
  char const * mask_name);
size_t yac_basic_grid_add_coordinates(
  struct yac_basic_grid * grid, enum yac_location location,
  coordinate_pointer coordinates, size_t count);
size_t yac_basic_grid_add_coordinates_nocpy(
  struct yac_basic_grid * grid, enum yac_location location,
  coordinate_pointer coordinates);
size_t yac_basic_grid_add_mask(
  struct yac_basic_grid * grid, enum yac_location location,
  int const * mask, size_t count, char const * mask_name);
size_t yac_basic_grid_add_mask_nocpy(
  struct yac_basic_grid * grid, enum yac_location location,
  int const * mask, char const * mask_name);
void yac_basic_grid_delete(struct yac_basic_grid * grid);

void yac_basic_grid_data_free(struct basic_grid_data grid);

enum yac_location yac_str2loc(char const * location);
char const * yac_loc2str(enum yac_location location);
enum yac_location yac_get_location(int const location);

struct yac_field_data yac_field_data_init();
size_t yac_field_data_add_mask_nocpy(
  struct yac_field_data * field_data, int const * mask,
  char const * mask_name);
size_t yac_field_data_add_coordinates_nocpy(
  struct yac_field_data * field_data, coordinate_pointer coordinates);

struct yac_field_data_set yac_field_data_set_init();
size_t yac_field_data_set_add_mask(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask, size_t count,
  char const * mask_name);
size_t yac_field_data_set_add_coordinates(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, coordinate_pointer coordinates,
  size_t count);
size_t yac_field_data_set_add_mask_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, int const * mask, char const * mask_name);
size_t yac_field_data_set_add_coordinates_nocpy(
  struct yac_field_data_set * field_data_set,
  enum yac_location location, coordinate_pointer coordinates);
void yac_field_data_set_free(struct yac_field_data_set field_data_set);

#endif // GRID_H

