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
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <netcdf.h>

#include "grid.h"
#include "geometry.h"
#include "read_mpiom_grid.h"
#include "utils.h"
#include "io_utils.h"

static void remove_duplicated_vertices(
  double ** vertex_lon, double ** vertex_lat,
  size_t * nbr_vertices, int * old_to_new_id);

void read_scrip_basic_grid_information(
  const char * filename, const char * grid_name,
  size_t cell_mask_size, int * cell_mask,
  size_t * num_vertices_, size_t * num_cells_, int ** num_vertices_per_cell_,
  int ** cell_to_vertex_, double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells) {

  size_t grid_name_len = strlen(grid_name) + 1;
  char crn_dim_name[5 + grid_name_len];
  char x_dim_name[2 + grid_name_len];
  char y_dim_name[2 + grid_name_len];
  char cla_var_name[4 + grid_name_len];
  char clo_var_name[4 + grid_name_len];
  char lat_var_name[4 + grid_name_len];
  char lon_var_name[4 + grid_name_len];

  snprintf(crn_dim_name, 5 + grid_name_len, "crn_%s", grid_name);
  snprintf(x_dim_name, 2 + grid_name_len, "x_%s", grid_name);
  snprintf(y_dim_name, 2 + grid_name_len, "y_%s", grid_name);
  snprintf(cla_var_name, 4 + grid_name_len, "%s.cla", grid_name);
  snprintf(clo_var_name, 4 + grid_name_len, "%s.clo", grid_name);
  snprintf(lat_var_name, 4 + grid_name_len, "%s.lat", grid_name);
  snprintf(lon_var_name, 4 + grid_name_len, "%s.lon", grid_name);

  int ncid;
  yac_nc_open(filename, NC_NOWRITE, &ncid);

  // get dimension ids
  int crn_dim_id;
  int x_dim_id;
  int y_dim_id;
  HANDLE_ERROR(nc_inq_dimid(ncid, crn_dim_name, &crn_dim_id));
  HANDLE_ERROR(nc_inq_dimid(ncid, x_dim_name, &x_dim_id));
  HANDLE_ERROR(nc_inq_dimid(ncid, y_dim_name, &y_dim_id));

  // get dimension length
  size_t crn_dim_len;
  size_t x_dim_len;
  size_t y_dim_len;
  HANDLE_ERROR(nc_inq_dimlen(ncid, crn_dim_id, &crn_dim_len));
  HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, &x_dim_len));
  HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, &y_dim_len));

  size_t num_cells = x_dim_len * y_dim_len;
  size_t num_vertices = num_cells * crn_dim_len;

  YAC_ASSERT(
    num_cells == cell_mask_size,
    "ERROR(read_scrip_basic_grid_information): "
    "cell mask size is inconsistent with number of grid cells")

  // get variable ids
  int cla_var_id;
  int clo_var_id;
  int lat_var_id;
  int lon_var_id;
  yac_nc_inq_varid(ncid, cla_var_name, &cla_var_id);
  yac_nc_inq_varid(ncid, clo_var_name, &clo_var_id);
  yac_nc_inq_varid(ncid, lat_var_name, &lat_var_id);
  yac_nc_inq_varid(ncid, lon_var_name, &lon_var_id);

  // allocate variables
  double * cla = xmalloc(num_vertices * sizeof(*cla));
  double * clo = xmalloc(num_vertices * sizeof(*clo));
  double * lat = xmalloc(num_cells * sizeof(*lat));
  double * lon = xmalloc(num_cells * sizeof(*lon));

  //read variables
  HANDLE_ERROR(nc_get_var_double(ncid, cla_var_id, cla));
  HANDLE_ERROR(nc_get_var_double(ncid, clo_var_id, clo));
  HANDLE_ERROR(nc_get_var_double(ncid, lat_var_id, lat));
  HANDLE_ERROR(nc_get_var_double(ncid, lon_var_id, lon));

  HANDLE_ERROR(nc_close(ncid));

  size_t * reorder_idx = xmalloc(num_vertices * sizeof(*reorder_idx));
  for (size_t y = 0, l = 0; y < y_dim_len; ++y)
    for (size_t x = 0; x < x_dim_len; ++x)
      for (size_t n = 0; n < crn_dim_len; ++n, ++l)
        reorder_idx[x + y * x_dim_len + n * x_dim_len * y_dim_len] = l;

  // remove duplicated vertices
  int * cell_to_vertex = xmalloc(num_vertices * sizeof(*cell_to_vertex));
  remove_duplicated_vertices(&clo, &cla, &num_vertices, cell_to_vertex);

  *x_vertices = clo;
  *y_vertices = cla;
  *x_cells = lon;
  *y_cells = lat;

  // we have to reorder cell_to_vertex
  yac_quicksort_index_size_t_int(
    reorder_idx, num_cells * crn_dim_len, cell_to_vertex);
  free(reorder_idx);

  // determine number of vertices per cell and compact cell_to_vertex
  int * num_vertices_per_cell =
    xmalloc(num_cells * sizeof(*num_vertices_per_cell));
  size_t total_num_cell_vertices = 0;
  int * to_vertices = cell_to_vertex;
  int * from_vertices = cell_to_vertex;
  for (size_t i = 0; i < num_cells; ++i) {
    size_t curr_num_vertices = 0;
    if (cell_mask[i]) {
      int prev_vertex = from_vertices[crn_dim_len-1];
      for (size_t j = 0; j < crn_dim_len; ++j, ++from_vertices) {
        int curr_vertex = *from_vertices;
        if (prev_vertex != curr_vertex) {
          prev_vertex = curr_vertex;
          if (to_vertices != from_vertices) *to_vertices = curr_vertex;
          ++curr_num_vertices;
          ++to_vertices;
        }
      }
    } else {
      from_vertices += crn_dim_len;
    }
    num_vertices_per_cell[i] = (int)curr_num_vertices;
    total_num_cell_vertices += curr_num_vertices;
  }


  if (total_num_cell_vertices != num_cells * crn_dim_len)
    cell_to_vertex =
      xrealloc(
        cell_to_vertex, total_num_cell_vertices * sizeof(*cell_to_vertex));

  *num_vertices_ = num_vertices;
  *num_cells_ = num_cells;
  *num_vertices_per_cell_ = num_vertices_per_cell;
  *cell_to_vertex_ = cell_to_vertex;
}

void read_scrip_mask_information(
  const char * filename, const char * grid_name, size_t * num_cells_,
  int ** cell_mask) {

  size_t grid_name_len = strlen(grid_name) + 1;
  char x_dim_name[2 + grid_name_len];
  char y_dim_name[2 + grid_name_len];
  char msk_var_name[4 + grid_name_len];

  snprintf(x_dim_name, 2 + grid_name_len, "x_%s", grid_name);
  snprintf(y_dim_name, 2 + grid_name_len, "y_%s", grid_name);
  snprintf(msk_var_name, 4 + grid_name_len, "%s.msk", grid_name);

  int ncid;
  yac_nc_open(filename, NC_NOWRITE, &ncid);

  // get dimension ids
  int x_dim_id;
  int y_dim_id;
  HANDLE_ERROR(nc_inq_dimid(ncid, x_dim_name, &x_dim_id));
  HANDLE_ERROR(nc_inq_dimid(ncid, y_dim_name, &y_dim_id));

  // get dimension length
  size_t x_dim_len;
  size_t y_dim_len;
  HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, &x_dim_len));
  HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, &y_dim_len));

  size_t num_cells = x_dim_len * y_dim_len;

  // get variable id
  int msk_var_id;
  yac_nc_inq_varid(ncid, msk_var_name, &msk_var_id);

  // allocate variable
  *cell_mask = xmalloc(num_cells * sizeof(**cell_mask));

  //read variable
  HANDLE_ERROR(nc_get_var_int(ncid, msk_var_id, *cell_mask));

  HANDLE_ERROR(nc_close(ncid));

  *num_cells_ = num_cells;
}

void read_scrip_grid_information(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value,
  size_t * num_vertices, size_t * num_cells, int ** num_vertices_per_cell,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells,
  int ** cell_to_vertex, int ** cell_core_mask) {

  size_t cell_mask_size;
  int * cell_mask;
  read_scrip_mask_information(
    mask_filename, grid_name, &cell_mask_size, &cell_mask);

  for (size_t i = 0; i < cell_mask_size; ++i)
    cell_mask[i] = cell_mask[i] == valid_mask_value;

  read_scrip_basic_grid_information(
    grid_filename, grid_name, cell_mask_size, cell_mask,
    num_vertices, num_cells, num_vertices_per_cell, cell_to_vertex,
    x_vertices, y_vertices, x_cells, y_cells);

  for (size_t i = 0; i < *num_vertices; ++i) {
    (*x_vertices)[i] *= YAC_RAD;
    (*y_vertices)[i] *= YAC_RAD;
  }
  for (size_t i = 0; i < *num_cells; ++i) {
    (*x_cells)[i] *= YAC_RAD;
    (*y_cells)[i] *= YAC_RAD;
  }

  YAC_ASSERT(
    *num_vertices <= INT_MAX,
    "ERROR(read_scrip_grid_information): too man vertices in grid")
  YAC_ASSERT(
    *num_cells <= INT_MAX,
    "ERROR(read_scrip_grid_information): too man cells in grid")

  *cell_core_mask = cell_mask;
}

struct basic_grid_data read_scrip_grid(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value) {

  size_t num_vertices;
  size_t num_cells;

  int * num_vertices_per_cell;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;
  int * cell_to_vertex;
  int * cell_core_mask;

  read_scrip_grid_information(
    grid_filename, mask_filename, grid_name, valid_mask_value,
    &num_vertices, &num_cells, &num_vertices_per_cell,
    &x_vertices, &y_vertices, &x_cells, &y_cells,
    &cell_to_vertex, &cell_core_mask);

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_unstruct(
      num_vertices, num_cells, num_vertices_per_cell,
      x_vertices, y_vertices, cell_to_vertex);

  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(cell_to_vertex);
  free(x_cells);
  free(y_cells);

  free(grid_data.core_cell_mask);
  grid_data.core_cell_mask = cell_core_mask;

  return grid_data;
}

/* ---------------------------------------------------------------- */

struct point_with_index {

  int32_t lon, lat;
  int64_t i;
};

static inline int compare_point_with_index_coord(
  const void * a,const void * b) {

   return ( (*(int64_t*)a < *(int64_t*)b) - (*(int64_t*)a > *(int64_t*)b) );
}

static void remove_duplicated_vertices(
  double ** vertex_lon, double ** vertex_lat,
  size_t * nbr_vertices, int * old_to_new_id) {

  struct point_with_index * sort_array =
    xmalloc(*nbr_vertices * sizeof(*sort_array));

  double const scale = (double)(2 << 20);
  int32_t const periode = (int32_t)(360.0 * scale);

  for (size_t i = 0; i < *nbr_vertices; ++i) {

    int32_t curr_lon = (int32_t)round((*vertex_lon)[i] * scale);
    int32_t curr_lat = (int32_t)round((*vertex_lat)[i] * scale);

    if ((curr_lat == (int32_t)(90.0 * scale)) ||
        (curr_lat == (int32_t)(-90.0 * scale))) {

      curr_lon = 0;

    } else {
      if (curr_lon <= 0)
         curr_lon = curr_lon - (((curr_lon - periode) / periode) * periode);
      else if (curr_lon > periode)
         curr_lon = curr_lon - ((curr_lon / periode) * periode);
    }

    sort_array[i].lon = curr_lon;
    sort_array[i].lat = curr_lat;

    sort_array[i].i = (int64_t)i;
  }

  yac_mergesort(sort_array, *nbr_vertices, sizeof(*sort_array),
                compare_point_with_index_coord);

  struct point_with_index dummy =
    {.lon = INT32_MAX, .lat = INT32_MAX, .i = INT64_MAX};
  struct point_with_index * prev = &dummy, * curr;
  size_t new_nbr_vertices = 0;
  for (size_t i = 0; i < *nbr_vertices; ++i) {

    curr = sort_array + i;
    int64_t idx = curr->i;
    if (compare_point_with_index_coord(prev, curr)) {
      new_nbr_vertices++;
      prev = curr;
    } else {
      curr->i = INT64_MAX;
    }
    old_to_new_id[idx] = (int)(new_nbr_vertices - 1);
  }

  double * new_vertex_lat = xmalloc(new_nbr_vertices * sizeof(*new_vertex_lat));
  double * new_vertex_lon = xmalloc(new_nbr_vertices * sizeof(*new_vertex_lon));

  new_nbr_vertices = 0;
  for (size_t i = 0; i < *nbr_vertices; ++i) {

    int64_t idx = sort_array[i].i;
    if (idx != INT64_MAX) {
      new_vertex_lon[new_nbr_vertices] = (*vertex_lon)[(size_t)idx];
      new_vertex_lat[new_nbr_vertices] = (*vertex_lat)[(size_t)idx];
      new_nbr_vertices++;
    }
  }
  free(*vertex_lon);
  free(*vertex_lat);
  *vertex_lon = new_vertex_lon;
  *vertex_lat = new_vertex_lat;
  *nbr_vertices = new_nbr_vertices;

  free(sort_array);
}
