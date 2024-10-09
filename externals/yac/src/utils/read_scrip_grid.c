// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <mpi.h>

#include "read_scrip_grid.h"
#include "geometry.h"
#include "utils_common.h"
#include "io_utils.h"
#include "yac_mpi_internal.h"

#ifdef YAC_NETCDF_ENABLED

#include <netcdf.h>

struct point_with_index {

  int32_t lon, lat;
  int rank;
  yac_int id;
  size_t i;
  double dlon, dlat;
};

struct cell_vertices_with_index {
  int num_vertices;
  int * cell_to_vertex;
  size_t i;
};

struct cell_vertex_ids_with_index {
  yac_int cell_id;
  yac_int * vertex_ids;
  size_t num_vertices;
  size_t i;
};

static void remove_duplicated_vertices(
  double ** vertex_lon, double ** vertex_lat,
  size_t * nbr_vertices, int * old_to_new_id);
static void remove_duplicated_vertices_parallel(
  double ** vertex_lon, double ** vertex_lat,
  size_t * nbr_vertices, int * old_to_new_id,
  yac_int ** vertex_ids, MPI_Comm comm);
static void detect_duplicated_cells(
  int * cell_to_vertex, int * num_vertices_per_cell, int * cell_mask,
  size_t nbr_cells, size_t ** duplicated_cell_idx,
  size_t ** orig_cell_idx, size_t * nbr_duplicated_cells);
static void detect_duplicated_cells_parallel(
  int * cell_to_vertex, int * num_vertices_per_cell, int * cell_mask,
  yac_int * cell_ids, size_t num_cells, yac_int * vertex_ids,
  MPI_Comm comm, size_t ** duplicated_cell_idx,
  yac_int ** orig_cell_ids, size_t * num_duplicated_cells_);

static inline int compare_int(
  const void * a,const void * b) {
  return *(int const *)a - *(int const*)b;
}

static inline int compare_yac_int(
  const void * a, const void * b) {
  return (*(yac_int const *)a > *(yac_int const*)b) -
         (*(yac_int const *)a < *(yac_int const*)b);
}

static inline int compare_cell_vertex_ids_with_index_ids(
  const void * a_, const void * b_) {

  struct cell_vertex_ids_with_index const * a =
    (struct cell_vertex_ids_with_index const *)a_;
  struct cell_vertex_ids_with_index const * b =
    (struct cell_vertex_ids_with_index const *)b_;

  size_t num_vertices_a = a->num_vertices;
  size_t num_vertices_b = b->num_vertices;
  int ret = (num_vertices_a > num_vertices_b) -
            (num_vertices_a < num_vertices_b);

  for (size_t i = 0; (i < num_vertices_a) && !ret; ++i)
    ret = (a->vertex_ids[i] > b->vertex_ids[i]) -
          (a->vertex_ids[i] < b->vertex_ids[i]);

  return ret;
}

static inline int compare_cell_vertex_ids_with_index_ids_id(
  const void * a_, const void * b_) {

  int ret = compare_cell_vertex_ids_with_index_ids(a_, b_);
  if (ret) return ret;

  struct cell_vertex_ids_with_index const * a =
    (struct cell_vertex_ids_with_index const *)a_;
  struct cell_vertex_ids_with_index const * b =
    (struct cell_vertex_ids_with_index const *)b_;

  return (a->cell_id > b->cell_id) - (a->cell_id < b->cell_id);
}

static inline int compare_point_with_index_coord_id(
  const void * a_, const void * b_) {

  struct point_with_index * a = (struct point_with_index *)a_;
  struct point_with_index * b = (struct point_with_index *)b_;

  int ret = (a->lon > b->lon) - (a->lon < b->lon);
  if (ret) return ret;
  ret = (a->lat > b->lat) - (a->lat < b->lat);
  if (ret) return ret;
  return (a->id > b->id) - (a->id < b->id);
}

static inline int compare_point_with_index_coord(
  const void * a_, const void * b_) {

  struct point_with_index * a = (struct point_with_index *)a_;
  struct point_with_index * b = (struct point_with_index *)b_;

  int ret = (a->lon > b->lon) - (a->lon < b->lon);
  if (ret) return ret;
  return (a->lat > b->lat) - (a->lat < b->lat);
}

static inline int compare_point_with_index_rank(
  const void * a_, const void * b_) {

  struct point_with_index * a = (struct point_with_index *)a_;
  struct point_with_index * b = (struct point_with_index *)b_;

  return (a->rank > b->rank) - (a->rank < b->rank);
}

static inline int compare_point_with_index_i(
  const void * a_, const void * b_) {

  struct point_with_index * a = (struct point_with_index *)a_;
  struct point_with_index * b = (struct point_with_index *)b_;

  return (a->i > b->i) - (a->i < b->i);
}

static void yac_read_scrip_basic_grid_information(
  const char * filename, const char * grid_name,
  size_t cell_mask_size, int * cell_mask,
  size_t * num_vertices_, size_t * num_cells_, int ** num_vertices_per_cell_,
  int ** cell_to_vertex_, double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells,
  size_t ** duplicated_cell_idx_, size_t ** orig_cell_idx_,
  size_t * nbr_duplicated_cells_) {

  size_t grid_name_len = strlen(grid_name) + 1;
  char cla_var_name[4 + grid_name_len];
  char clo_var_name[4 + grid_name_len];
  char lat_var_name[4 + grid_name_len];
  char lon_var_name[4 + grid_name_len];

  snprintf(cla_var_name, 4 + grid_name_len, "%s.cla", grid_name);
  snprintf(clo_var_name, 4 + grid_name_len, "%s.clo", grid_name);
  snprintf(lat_var_name, 4 + grid_name_len, "%s.lat", grid_name);
  snprintf(lon_var_name, 4 + grid_name_len, "%s.lon", grid_name);

  int ncid;
  yac_nc_open(filename, NC_NOWRITE, &ncid);

  // get variable ids
  int cla_var_id;
  int clo_var_id;
  int lat_var_id;
  int lon_var_id;
  yac_nc_inq_varid(ncid, cla_var_name, &cla_var_id);
  yac_nc_inq_varid(ncid, clo_var_name, &clo_var_id);
  yac_nc_inq_varid(ncid, lat_var_name, &lat_var_id);
  yac_nc_inq_varid(ncid, lon_var_name, &lon_var_id);

  // get dimension ids
  int ndims;
  int dim_ids[NC_MAX_VAR_DIMS];
  YAC_HANDLE_ERROR(
    nc_inq_var(ncid, cla_var_id, NULL, NULL, &ndims, dim_ids, NULL));
  YAC_ASSERT_F(
    ndims == 3,
    "ERROR(yac_read_scrip_basic_grid_information): "
    "invalid number of dimensions for variable \"%s\" (has to be 3 not %d)",
    cla_var_name, ndims);
  int crn_dim_id = dim_ids[0];
  int x_dim_id = dim_ids[2];
  int y_dim_id = dim_ids[1];

  // get dimension length
  size_t crn_dim_len;
  size_t x_dim_len;
  size_t y_dim_len;
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, crn_dim_id, &crn_dim_len));
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, &x_dim_len));
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, &y_dim_len));

  size_t num_cells = x_dim_len * y_dim_len;
  size_t num_vertices = num_cells * crn_dim_len;

  YAC_ASSERT(
    num_cells == cell_mask_size,
    "ERROR(yac_read_scrip_basic_grid_information): "
    "cell mask size is inconsistent with number of grid cells")

  // allocate variables
  double * cla = xmalloc(num_vertices * sizeof(*cla));
  double * clo = xmalloc(num_vertices * sizeof(*clo));

  //read variables
  YAC_HANDLE_ERROR(nc_get_var_double(ncid, cla_var_id, cla));
  YAC_HANDLE_ERROR(nc_get_var_double(ncid, clo_var_id, clo));

  if ((x_cells != NULL) && (y_cells != NULL)) {
    double * lat = xmalloc(num_cells * sizeof(*lat));
    double * lon = xmalloc(num_cells * sizeof(*lon));
    YAC_HANDLE_ERROR(nc_get_var_double(ncid, lat_var_id, lat));
    YAC_HANDLE_ERROR(nc_get_var_double(ncid, lon_var_id, lon));
    *x_cells = lon;
    *y_cells = lat;
  }

  YAC_HANDLE_ERROR(nc_close(ncid));

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

  // make a copy of cell_to_vertex and sort the vertex indices of each cell,
  // this is required by the detection of duplicated cells
  int * cell_to_vertex_copy =
    xmalloc(total_num_cell_vertices * sizeof(*cell_to_vertex_copy));
  memcpy(
    cell_to_vertex_copy, cell_to_vertex,
    total_num_cell_vertices * sizeof(*cell_to_vertex_copy));
  for (size_t i = 0, offset = 0; i < num_cells; ++i) {
    qsort(
      cell_to_vertex_copy + offset, (size_t)(num_vertices_per_cell[i]),
      sizeof(*cell_to_vertex_copy), compare_int);
    offset += num_vertices_per_cell[i];
  }

  // detect duplicated cells
  size_t * duplicated_cell_idx;
  size_t * orig_cell_idx;
  size_t nbr_duplicated_cells;
  detect_duplicated_cells(
    cell_to_vertex_copy, num_vertices_per_cell, cell_mask, num_cells,
    &duplicated_cell_idx, &orig_cell_idx, &nbr_duplicated_cells);
  free(cell_to_vertex_copy);

  // mask out duplicated cells
  for (size_t i = 0; i < nbr_duplicated_cells; ++i)
    cell_mask[duplicated_cell_idx[i]] = 0;

  if ((duplicated_cell_idx_ != NULL) && (orig_cell_idx_ != NULL) &&
      (nbr_duplicated_cells_ != NULL)) {
    *duplicated_cell_idx_ = duplicated_cell_idx;
    *orig_cell_idx_ = orig_cell_idx;
    *nbr_duplicated_cells_ = nbr_duplicated_cells;
  } else {
    free(duplicated_cell_idx);
    free(orig_cell_idx);
  }

  *num_vertices_ = num_vertices;
  *num_cells_ = num_cells;
  *num_vertices_per_cell_ = num_vertices_per_cell;
  *cell_to_vertex_ = cell_to_vertex;
}

static void yac_read_part_scrip_basic_grid_information(
  const char * filename, const char * grid_name,
  size_t cell_mask_size, int * cell_mask,
  size_t * num_vertices_, size_t * num_cells_, int ** num_vertices_per_cell_,
  int ** cell_to_vertex_,
  double ** x_vertices_, double ** y_vertices_, yac_int ** vertex_ids_,
  double ** x_cells_, double ** y_cells_, yac_int ** cell_ids_,
  size_t ** duplicated_cell_idx_, yac_int ** orig_cell_ids_,
  size_t * nbr_duplicated_cells_,
  int io_rank_idx, int num_io_ranks, MPI_Comm comm) {

  size_t num_cells = 0;
  size_t num_vertices = 0;
  size_t * reorder_idx = NULL;
  size_t max_num_vertices_per_cell = 0;

  double * x_vertices = NULL;
  double * y_vertices = NULL;
  yac_int * vertex_ids = NULL;
  double * x_cells = NULL;
  double * y_cells = NULL;
  yac_int * cell_ids = NULL;

  if ((io_rank_idx >= 0) && (io_rank_idx < num_io_ranks)) {

    size_t grid_name_len = strlen(grid_name) + 1;
    char cla_var_name[4 + grid_name_len];
    char clo_var_name[4 + grid_name_len];
    char lat_var_name[4 + grid_name_len];
    char lon_var_name[4 + grid_name_len];

    snprintf(cla_var_name, 4 + grid_name_len, "%s.cla", grid_name);
    snprintf(clo_var_name, 4 + grid_name_len, "%s.clo", grid_name);
    snprintf(lat_var_name, 4 + grid_name_len, "%s.lat", grid_name);
    snprintf(lon_var_name, 4 + grid_name_len, "%s.lon", grid_name);

    int ncid;
    yac_nc_open(filename, NC_NOWRITE, &ncid);

    // get variable ids
    int cla_var_id, clo_var_id, lat_var_id, lon_var_id;
    yac_nc_inq_varid(ncid, cla_var_name, &cla_var_id);
    yac_nc_inq_varid(ncid, clo_var_name, &clo_var_id);
    yac_nc_inq_varid(ncid, lat_var_name, &lat_var_id);
    yac_nc_inq_varid(ncid, lon_var_name, &lon_var_id);

    // get dimension ids
    int ndims;
    int dim_ids[NC_MAX_VAR_DIMS];
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, cla_var_id, NULL, NULL, &ndims, dim_ids, NULL));
    YAC_ASSERT_F(
      ndims == 3,
      "ERROR(yac_read_part_scrip_basic_grid_information): "
      "invalid number of dimensions for variable \"%s\" (has to be 3 not %d)",
      cla_var_name, ndims);
    int crn_dim_id = dim_ids[0];
    int x_dim_id = dim_ids[2];
    int y_dim_id = dim_ids[1];

    // get dimension length
    size_t crn_dim_len, x_dim_len, y_dim_len;
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, crn_dim_id, &crn_dim_len));
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, &x_dim_len));
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, &y_dim_len));

    // decompose x dimension among io ranks
    size_t x_start_local =
      (size_t)
        (((unsigned long)x_dim_len * (unsigned long)io_rank_idx) /
        (unsigned long)num_io_ranks);
    size_t x_count_local =
      (size_t)
        (((unsigned long)x_dim_len * (io_rank_idx + 1)) /
        (unsigned long)num_io_ranks) - x_start_local;

    num_cells = x_count_local * y_dim_len;
    num_vertices = num_cells * crn_dim_len;
    max_num_vertices_per_cell = crn_dim_len;

    YAC_ASSERT(
      num_cells == cell_mask_size,
      "ERROR(yac_read_part_scrip_basic_grid_information): "
      "cell mask size is inconsistent with number of grid cells")

    // allocate variables
    x_vertices = xmalloc(num_vertices * sizeof(*x_vertices));
    y_vertices = xmalloc(num_vertices * sizeof(*y_vertices));

    //read variables
    size_t start[3] = {0, 0, x_start_local};
    size_t count[3] = {crn_dim_len, y_dim_len, x_count_local};
    YAC_HANDLE_ERROR(
      nc_get_vara_double(ncid, clo_var_id, start, count, x_vertices));
    YAC_HANDLE_ERROR(
      nc_get_vara_double(ncid, cla_var_id, start, count, y_vertices));

    if ((x_cells_ != NULL) && (y_cells_ != NULL)) {
      x_cells = xmalloc(num_cells * sizeof(*x_cells));
      y_cells = xmalloc(num_cells * sizeof(*y_cells));
      size_t start[2] = {0, x_start_local};
      size_t count[2] = {y_dim_len, x_count_local};
      YAC_HANDLE_ERROR(
        nc_get_vara_double(ncid, lon_var_id, start, count, x_cells));
      YAC_HANDLE_ERROR(
        nc_get_vara_double(ncid, lat_var_id, start, count, y_cells));
    }

    YAC_HANDLE_ERROR(nc_close(ncid));

    reorder_idx = xmalloc(num_vertices * sizeof(*reorder_idx));
    for (size_t y = 0, l = 0; y < y_dim_len; ++y)
      for (size_t x = 0; x < x_count_local; ++x)
        for (size_t n = 0; n < crn_dim_len; ++n, ++l)
          reorder_idx[x + y * x_count_local + n * x_count_local * y_dim_len] = l;

    // generate global cell ids
    cell_ids = xmalloc(num_cells * sizeof(*cell_ids));
    for (size_t y = 0, k = 0; y < y_dim_len; ++y)
      for (size_t x = 0; x < x_count_local; ++x, ++k)
        cell_ids[k] = y * x_dim_len + x + x_start_local;


    YAC_ASSERT_F(
      (y_dim_len * x_dim_len * crn_dim_len) <= XT_INT_MAX,
      "ERROR(yac_read_part_scrip_basic_grid_information): "
      "number of vertices exceed maximum global id (%zu)",
      (size_t)XT_INT_MAX);

    // generate global vertex ids
    vertex_ids = xmalloc(num_cells * crn_dim_len * sizeof(*vertex_ids));
    for (size_t n = 0, l = 0; n < crn_dim_len; ++n)
      for (size_t i = 0; i < num_cells; ++i, ++l)
        vertex_ids[l] = cell_ids[i] * crn_dim_len + n;
  }

  // remove duplicated vertices
  int * cell_to_vertex = xmalloc(num_vertices * sizeof(*cell_to_vertex));
  remove_duplicated_vertices_parallel(
    &x_vertices, &y_vertices, &num_vertices,
    cell_to_vertex, &vertex_ids, comm);

  // we have to reorder cell_to_vertex
  yac_quicksort_index_size_t_int(
    reorder_idx, num_cells * max_num_vertices_per_cell, cell_to_vertex);
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
      int prev_vertex = from_vertices[max_num_vertices_per_cell-1];
      for (size_t j = 0; j < max_num_vertices_per_cell; ++j, ++from_vertices) {
        int curr_vertex = *from_vertices;
        if (prev_vertex != curr_vertex) {
          prev_vertex = curr_vertex;
          if (to_vertices != from_vertices) *to_vertices = curr_vertex;
          ++curr_num_vertices;
          ++to_vertices;
        }
      }
    } else {
      from_vertices += max_num_vertices_per_cell;
    }
    num_vertices_per_cell[i] = (int)curr_num_vertices;
    total_num_cell_vertices += curr_num_vertices;
  }

  if (total_num_cell_vertices != num_cells * max_num_vertices_per_cell)
    cell_to_vertex =
      xrealloc(
        cell_to_vertex, total_num_cell_vertices * sizeof(*cell_to_vertex));

  // detect duplicated cells
  size_t * duplicated_cell_idx;
  yac_int * orig_cell_ids;
  size_t nbr_duplicated_cells;

  detect_duplicated_cells_parallel(
    cell_to_vertex, num_vertices_per_cell, cell_mask, cell_ids, num_cells,
    vertex_ids, comm, &duplicated_cell_idx, &orig_cell_ids,
    &nbr_duplicated_cells);

  // mask out duplicated cells
  for (size_t i = 0; i < nbr_duplicated_cells; ++i)
    cell_mask[duplicated_cell_idx[i]] = 0;

  if ((duplicated_cell_idx_ != NULL) && (orig_cell_ids_ != NULL) &&
      (nbr_duplicated_cells_ != NULL)) {
    *duplicated_cell_idx_ = duplicated_cell_idx;
    *orig_cell_ids_ = orig_cell_ids;
    *nbr_duplicated_cells_ = nbr_duplicated_cells;
  } else {
    free(duplicated_cell_idx);
    free(orig_cell_ids);
  }

  *num_vertices_ = num_vertices;
  *num_cells_ = num_cells;
  *num_vertices_per_cell_ = num_vertices_per_cell;
  *cell_to_vertex_ = cell_to_vertex;
  *x_vertices_ = x_vertices;
  *y_vertices_ = y_vertices;
  *vertex_ids_ = vertex_ids;
  if ((x_cells_ != NULL) && (y_cells_ != NULL)) {
    *x_cells_ = x_cells;
    *y_cells_ = y_cells;
  }
  *cell_ids_ = cell_ids;
}

static void yac_read_scrip_mask_information(
  const char * filename, const char * grid_name, size_t * num_cells_,
  int ** cell_mask) {

  size_t grid_name_len = strlen(grid_name) + 1;
  char msk_var_name[4 + grid_name_len];

  snprintf(msk_var_name, 4 + grid_name_len, "%s.msk", grid_name);

  int ncid;
  yac_nc_open(filename, NC_NOWRITE, &ncid);

  // get variable id
  int msk_var_id;
  yac_nc_inq_varid(ncid, msk_var_name, &msk_var_id);

  // get dimension ids
  int ndims;
  int dim_ids[NC_MAX_VAR_DIMS];
  YAC_HANDLE_ERROR(
    nc_inq_var(ncid, msk_var_id, NULL, NULL, &ndims, dim_ids, NULL));
  YAC_ASSERT_F(
    ndims == 2,
    "ERROR(yac_read_scrip_mask_information): "
    "invalid number of dimensions for variable \"%s\" (has to be 2 not %d)",
    msk_var_name, ndims);
  int x_dim_id = dim_ids[1];
  int y_dim_id = dim_ids[0];

  // get dimension length
  size_t x_dim_len;
  size_t y_dim_len;
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, &x_dim_len));
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, &y_dim_len));

  size_t num_cells = x_dim_len * y_dim_len;

  // allocate variable
  *cell_mask = xmalloc(num_cells * sizeof(**cell_mask));

  //read variable
  YAC_HANDLE_ERROR(nc_get_var_int(ncid, msk_var_id, *cell_mask));

  YAC_HANDLE_ERROR(nc_close(ncid));

  *num_cells_ = num_cells;
}

static void yac_read_part_scrip_mask_information(
  const char * filename, const char * grid_name, size_t * num_cells_,
  int ** cell_mask, int io_rank_idx, int num_io_ranks) {

  if ((io_rank_idx >= 0) && (io_rank_idx < num_io_ranks)) {

    size_t grid_name_len = strlen(grid_name) + 1;
    char msk_var_name[4 + grid_name_len];

    snprintf(msk_var_name, 4 + grid_name_len, "%s.msk", grid_name);

    int ncid;
    yac_nc_open(filename, NC_NOWRITE, &ncid);

    // get variable id
    int msk_var_id;
    yac_nc_inq_varid(ncid, msk_var_name, &msk_var_id);

    // get dimension ids
    int ndims;
    int dim_ids[NC_MAX_VAR_DIMS];
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, msk_var_id, NULL, NULL, &ndims, dim_ids, NULL));
    YAC_ASSERT_F(
      ndims == 2,
      "ERROR(yac_read_part_scrip_mask_information): "
      "invalid number of dimensions for variable \"%s\" (has to be 2 not %d)",
      msk_var_name, ndims);
    int x_dim_id = dim_ids[1];
    int y_dim_id = dim_ids[0];

    // get dimension length
    size_t x_dim_len;
    size_t y_dim_len;
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, &x_dim_len));
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, &y_dim_len));

    // decompose x dimension among io ranks
    size_t x_start_local =
      (size_t)
        (((unsigned long)x_dim_len * (unsigned long)io_rank_idx) /
        (unsigned long)num_io_ranks);
    size_t x_count_local =
      (size_t)
        (((unsigned long)x_dim_len * (io_rank_idx+1)) /
        (unsigned long)num_io_ranks) - x_start_local;

    size_t num_cells = x_count_local * y_dim_len;

    // allocate variable
    *cell_mask = xmalloc(num_cells * sizeof(**cell_mask));

    //read variable
    size_t start[2] = {0, x_start_local};
    size_t count[2] = {y_dim_len, x_count_local};
    YAC_HANDLE_ERROR(
      nc_get_vara_int(ncid, msk_var_id, start, count, *cell_mask));

    YAC_HANDLE_ERROR(nc_close(ncid));

    *num_cells_ = num_cells;

  } else {
    *cell_mask = NULL;
    *num_cells_ = 0;
  }
}

void yac_read_scrip_grid_information(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value,
  size_t * num_vertices, size_t * num_cells, int ** num_vertices_per_cell,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells,
  int ** cell_to_vertex, int ** cell_core_mask, size_t ** duplicated_cell_idx,
  size_t ** orig_cell_idx, size_t * nbr_duplicated_cells) {

  size_t cell_mask_size;
  int * cell_mask;
  yac_read_scrip_mask_information(
    mask_filename, grid_name, &cell_mask_size, &cell_mask);

  for (size_t i = 0; i < cell_mask_size; ++i)
    cell_mask[i] = cell_mask[i] == valid_mask_value;

  yac_read_scrip_basic_grid_information(
    grid_filename, grid_name, cell_mask_size, cell_mask,
    num_vertices, num_cells, num_vertices_per_cell, cell_to_vertex,
    x_vertices, y_vertices, x_cells, y_cells, duplicated_cell_idx,
    orig_cell_idx, nbr_duplicated_cells);

  for (size_t i = 0; i < *num_vertices; ++i) {
    (*x_vertices)[i] *= YAC_RAD;
    (*y_vertices)[i] *= YAC_RAD;
  }
  if ((x_cells != NULL) && (y_cells != NULL)) {
    for (size_t i = 0; i < *num_cells; ++i) {
      (*x_cells)[i] *= YAC_RAD;
      (*y_cells)[i] *= YAC_RAD;
    }
  }

  YAC_ASSERT(
    *num_vertices <= INT_MAX,
    "ERROR(yac_read_scrip_grid_information): "
    "too man vertices in grid")
  YAC_ASSERT(
    *num_cells <= INT_MAX,
    "ERROR(yac_read_scrip_grid_information): "
    "too man cells in grid")

  *cell_core_mask = cell_mask;
}

static void yac_read_part_scrip_grid_information(
  char const * grid_filename, char const * mask_filename,
  MPI_Comm comm, char const * grid_name, int valid_mask_value,
  size_t * num_vertices, size_t * num_cells, int ** num_vertices_per_cell,
  double ** x_vertices, double ** y_vertices, yac_int ** vertex_ids,
  double ** x_cells, double ** y_cells, yac_int ** cell_ids,
  int ** cell_to_vertex, int ** cell_core_mask, size_t ** duplicated_cell_idx,
  yac_int ** orig_cell_ids, size_t * nbr_duplicated_cells) {

  // if no communicator was provided
  if (comm == MPI_COMM_NULL) {

    size_t * orig_cell_idx = NULL;

    yac_read_scrip_grid_information(
      grid_filename, mask_filename, grid_name, valid_mask_value,
      num_vertices, num_cells, num_vertices_per_cell,
      x_vertices, y_vertices, x_cells, y_cells, cell_to_vertex,
      cell_core_mask, duplicated_cell_idx,
      (orig_cell_ids != NULL)?&orig_cell_idx:NULL,
      nbr_duplicated_cells);

    *vertex_ids = xmalloc(*num_vertices * sizeof(**vertex_ids));
    for (size_t i = 0; i < *num_vertices; ++i) (*vertex_ids)[i] = i;

    *cell_ids = xmalloc(*num_cells * sizeof(**cell_ids));
    for (size_t i = 0; i < *num_cells; ++i) (*cell_ids)[i] = i;

    if (orig_cell_ids != NULL) {
      *orig_cell_ids =
        xmalloc(*nbr_duplicated_cells * sizeof(**orig_cell_ids));
      for (size_t i = 0; i < *nbr_duplicated_cells; ++i)
        (*orig_cell_ids)[i] = (*cell_ids)[orig_cell_idx[i]];
      free(orig_cell_idx);
    }

  } else {

    // get io configuration
    int local_is_io, * io_ranks, num_io_ranks;
    yac_get_io_ranks(comm, &local_is_io, &io_ranks, &num_io_ranks);

    int comm_rank;
    yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

    int io_rank_idx = 0;
    while ((io_rank_idx < num_io_ranks) &&
            (comm_rank != io_ranks[io_rank_idx]))
      ++io_rank_idx;
    YAC_ASSERT(
      !local_is_io || (io_rank_idx < num_io_ranks),
      "ERROR(yac_read_part_scrip_grid_information): "
      "unable to determine io_rank_idx");

    int * cell_mask;
    size_t cell_mask_size;
    yac_read_part_scrip_mask_information(
      mask_filename, grid_name, &cell_mask_size, &cell_mask,
      io_rank_idx, num_io_ranks);

    for (size_t i = 0; i < cell_mask_size; ++i)
      cell_mask[i] = cell_mask[i] == valid_mask_value;

    yac_read_part_scrip_basic_grid_information(
      grid_filename, grid_name, cell_mask_size, cell_mask,
      num_vertices, num_cells, num_vertices_per_cell, cell_to_vertex,
      x_vertices, y_vertices, vertex_ids,
      x_cells, y_cells, cell_ids, duplicated_cell_idx,
      orig_cell_ids, nbr_duplicated_cells,
      io_rank_idx, num_io_ranks, comm);

    for (size_t i = 0; i < *num_vertices; ++i) {
      (*x_vertices)[i] *= YAC_RAD;
      (*y_vertices)[i] *= YAC_RAD;
    }
    if ((x_cells != NULL) && (y_cells != NULL)) {
      for (size_t i = 0; i < *num_cells; ++i) {
        (*x_cells)[i] *= YAC_RAD;
        (*y_cells)[i] *= YAC_RAD;
      }
    }

    YAC_ASSERT(
      *num_vertices <= INT_MAX,
      "ERROR(yac_read_part_scrip_grid_information): "
      "too man vertices in grid")
    YAC_ASSERT(
      *num_cells <= INT_MAX,
      "ERROR(yac_read_part_scrip_grid_information): "
      "too man cells in grid")
    YAC_ASSERT_F(
      cell_mask_size == *num_cells,
      "ERROR(yac_read_part_scrip_grid_information): "
      "mask and grid size do not match "
      "(mask_file: \"%s\" grid_file: \"%s\" grid_name: \"%s\"",
      mask_filename, grid_filename, grid_name)

    *cell_core_mask = cell_mask;

    free(io_ranks);
  }
}

static struct yac_basic_grid_data yac_read_scrip_basic_grid_data_(
  char const * grid_filename, char const * mask_filename,
  MPI_Comm comm, char const * grid_name, int valid_mask_value, int use_ll,
  double ** x_cells, double ** y_cells, size_t ** duplicated_cell_idx,
  yac_int ** orig_cell_global_ids, size_t * nbr_duplicated_cells) {

  size_t num_vertices;
  size_t num_cells;

  int * num_vertices_per_cell;
  double * x_vertices;
  double * y_vertices;
  yac_int * vertex_ids;
  yac_int * cell_ids;
  int * cell_to_vertex;
  int * cell_core_mask;

  yac_read_part_scrip_grid_information(
    grid_filename, mask_filename, comm, grid_name, valid_mask_value,
    &num_vertices, &num_cells, &num_vertices_per_cell,
    &x_vertices, &y_vertices, &vertex_ids, x_cells, y_cells, &cell_ids,
    &cell_to_vertex, &cell_core_mask, duplicated_cell_idx,
    orig_cell_global_ids, nbr_duplicated_cells);

  // check whether any cell has more then four vertices
  // --> write out a warning and switch to gc edges
  if (use_ll) {
    for (size_t i = 0; (i < num_cells) && use_ll; ++i)
      use_ll = !cell_core_mask[i] || (num_vertices_per_cell[i] <= 4);
    int comm_rank = 0;
    if (comm != MPI_COMM_NULL) {
      yac_mpi_call(
        MPI_Allreduce(
          MPI_IN_PLACE, &use_ll, 1, MPI_INT, MPI_MIN, comm), comm);
      yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
    }
    if (!use_ll && (comm_rank == 0))
      fprintf(
        stderr, "WARNING(yac_read_scrip_basic_grid_data_): "
        "grid \"%s\" from \"%s\": required LL but stored as GC "
        "(>4 vertices per cell)\n",
        grid_name, grid_filename);
  }

  struct yac_basic_grid_data
    (*generate_basic_grid_data_unstruct_ptr)(
      size_t, size_t, int *, double *, double *, int *) =
    (use_ll)?
      yac_generate_basic_grid_data_unstruct_ll:
      yac_generate_basic_grid_data_unstruct;

  struct yac_basic_grid_data grid_data =
    generate_basic_grid_data_unstruct_ptr(
      num_vertices, num_cells, num_vertices_per_cell,
      x_vertices, y_vertices, cell_to_vertex);

  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(cell_to_vertex);

  free(grid_data.core_cell_mask);
  grid_data.core_cell_mask = cell_core_mask;

  free(grid_data.vertex_ids);
  grid_data.vertex_ids = vertex_ids;
  free(grid_data.cell_ids);
  grid_data.cell_ids = cell_ids;

  return grid_data;
}

struct yac_basic_grid_data yac_read_scrip_basic_grid_data(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value, int use_ll_edges) {

  return
    yac_read_scrip_basic_grid_data_(
      grid_filename, mask_filename, MPI_COMM_NULL,
      grid_name, valid_mask_value, use_ll_edges,
      NULL, NULL, NULL, NULL, NULL);
}

struct yac_basic_grid * yac_read_scrip_basic_grid_parallel(
  char const * grid_filename, char const * mask_filename,
  MPI_Comm comm, char const * grid_name, int valid_mask_value,
  char const * name, int use_ll_edges, size_t * cell_coord_idx,
  size_t ** duplicated_cell_idx, yac_int ** orig_cell_global_ids,
  size_t * nbr_duplicated_cells) {

  // generate basic grid
  double * x_cells, * y_cells;
  struct yac_basic_grid * basic_grid =
    yac_basic_grid_new(
      name,
      yac_read_scrip_basic_grid_data_(
        grid_filename, mask_filename, comm, grid_name, valid_mask_value,
        use_ll_edges, (cell_coord_idx != NULL)?&x_cells:NULL,
        (cell_coord_idx != NULL)?&y_cells:NULL, duplicated_cell_idx,
        orig_cell_global_ids, nbr_duplicated_cells));

  if (cell_coord_idx != NULL) {

    // generate 3d cell center coordinates
    size_t num_cells =
      yac_basic_grid_get_data_size(basic_grid, YAC_LOC_CELL);
    yac_coordinate_pointer cell_coords =
      xmalloc(num_cells * sizeof(*cell_coords));
    for (size_t i = 0; i < num_cells; ++i)
      LLtoXYZ(x_cells[i], y_cells[i], cell_coords[i]);
    free(x_cells);
    free(y_cells);

    // register cell center coordinates in basic grid
    *cell_coord_idx =
      yac_basic_grid_add_coordinates_nocpy(
        basic_grid, YAC_LOC_CELL, cell_coords);
  }

  return basic_grid;
}

struct yac_basic_grid * yac_read_scrip_basic_grid_parallel_f2c(
  char const * grid_filename, char const * mask_filename,
  MPI_Fint comm, char const * grid_name, int valid_mask_value,
  char const * name, int use_ll_edges, size_t * cell_coord_idx,
  size_t ** duplicated_cell_idx, yac_int ** orig_cell_global_ids,
  size_t * nbr_duplicated_cells) {

  return
    yac_read_scrip_basic_grid_parallel(
      grid_filename, mask_filename, MPI_Comm_f2c(comm), grid_name,
      valid_mask_value, name, use_ll_edges, cell_coord_idx,
      duplicated_cell_idx, orig_cell_global_ids, nbr_duplicated_cells);
}

struct yac_basic_grid * yac_read_scrip_basic_grid(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value, char const * name,
  int use_ll_edges, size_t * cell_coord_idx,
  size_t ** duplicated_cell_idx, yac_int ** orig_cell_global_ids,
  size_t * nbr_duplicated_cells) {

  return
    yac_read_scrip_basic_grid_parallel(
      grid_filename, mask_filename, MPI_COMM_NULL, grid_name, valid_mask_value,
      name, use_ll_edges, cell_coord_idx, duplicated_cell_idx,
      orig_cell_global_ids, nbr_duplicated_cells);
}

/* ---------------------------------------------------------------- */

static inline int compare_cell_vertices_with_index(
  const void * a,const void * b) {

  int ret = ((struct cell_vertices_with_index const *)a)->num_vertices -
            ((struct cell_vertices_with_index const *)b)->num_vertices;
  if (ret) return ret;
  int count = ((struct cell_vertices_with_index const *)a)->num_vertices;
  for (int i = 0; (i < count) && !ret; ++i)
    ret = ((struct cell_vertices_with_index const *)a)->cell_to_vertex[i] -
          ((struct cell_vertices_with_index const *)b)->cell_to_vertex[i];
  return ret;
}

static struct point_with_index * remove_duplicated_vertices_(
  double ** vertex_lon, double ** vertex_lat, yac_int ** vertex_ids,
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
    sort_array[i].dlon = (*vertex_lon)[i];
    sort_array[i].dlat = (*vertex_lat)[i];
    sort_array[i].id = (vertex_ids != NULL)?(*vertex_ids)[i]:XT_INT_MAX;

    sort_array[i].i = i;
  }

  yac_mergesort(sort_array, *nbr_vertices, sizeof(*sort_array),
                compare_point_with_index_coord_id);

  struct point_with_index dummy =
    {.lon = INT32_MAX, .lat = INT32_MAX, .i = SIZE_MAX};
  struct point_with_index * prev = &dummy, * curr = sort_array;
  size_t new_nbr_vertices = 0;
  for (size_t i = 0; i < *nbr_vertices; ++i, ++curr) {

    if (compare_point_with_index_coord(prev, curr)) {

      (*vertex_lon)[new_nbr_vertices] = curr->dlon;
      (*vertex_lat)[new_nbr_vertices] = curr->dlat;
      if (vertex_ids != NULL)
        (*vertex_ids)[new_nbr_vertices] = curr->id;
      sort_array[new_nbr_vertices] = *curr;
      new_nbr_vertices++;
      prev = curr;
    }
    old_to_new_id[curr->i] = (int)(new_nbr_vertices - 1);
  }

  (*vertex_lon) = xrealloc(*vertex_lon, new_nbr_vertices * sizeof(**vertex_lon));
  (*vertex_lat) = xrealloc(*vertex_lat, new_nbr_vertices * sizeof(**vertex_lat));
  if (vertex_ids != NULL)
    (*vertex_ids) = xrealloc(*vertex_ids, new_nbr_vertices * sizeof(**vertex_ids));
  *nbr_vertices = new_nbr_vertices;

  return sort_array;
}

static void remove_duplicated_vertices(
  double ** vertex_lon, double ** vertex_lat,
  size_t * nbr_vertices, int * old_to_new_id) {

  free(
    remove_duplicated_vertices_(
      vertex_lon, vertex_lat, NULL, nbr_vertices, old_to_new_id));
}

// djb2 hash function
unsigned long djb2_hash(unsigned char * values, size_t count) {

    unsigned long hash = 5381;

    for (size_t i = 0; i < count; ++i) {
      unsigned int value = values[i];
      hash = ((hash << 5) + hash) + value; /* hash * 33 + value */
    }

    return hash;
}

static int compute_bucket(
  void * data, size_t data_size, int comm_size) {

  return
    djb2_hash((unsigned char*)data, data_size) % (unsigned long)comm_size;
}

static MPI_Datatype get_point_with_index_mpi_datatype(MPI_Comm comm) {

  struct point_with_index dummy;
  MPI_Datatype point_with_index_dt;
  int array_of_blocklengths[] = {1, 1, 1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.lon) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.lat) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.dlon) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.dlat) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_INT32_T, MPI_INT32_T, yac_int_dt, MPI_DOUBLE, MPI_DOUBLE};
  yac_mpi_call(
    MPI_Type_create_struct(5, array_of_blocklengths, array_of_displacements,
                           array_of_types, &point_with_index_dt), comm);
  return yac_create_resized(point_with_index_dt, sizeof(dummy), comm);
}

static void remove_duplicated_vertices_parallel(
  double ** vertex_lon, double ** vertex_lat,
  size_t * nbr_vertices_, int * old_to_new_id, yac_int ** vertex_ids,
  MPI_Comm comm) {

  // remove the locally duplicated vertices
  struct point_with_index * sort_array =
    remove_duplicated_vertices_(
      vertex_lon, vertex_lat, vertex_ids, nbr_vertices_, old_to_new_id);

  sort_array = xrealloc(sort_array, *nbr_vertices_ * sizeof(*sort_array));

  size_t nbr_vertices = *nbr_vertices_;

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // set up counts and displs for redistribution of points
  for (size_t i = 0; i < nbr_vertices; ++i) {
    int rank =
      compute_bucket(
        &(sort_array[i].lon), 2 * sizeof(sort_array[i].lon), comm_size);
    sendcounts[rank]++;
    sort_array[i].rank = rank;
    sort_array[i].i = i;
  }

  // sort points by bucket ranks
  yac_mergesort(sort_array, nbr_vertices, sizeof(*sort_array),
                compare_point_with_index_rank);

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t recv_count = recvcounts[comm_size - 1] + rdispls[comm_size - 1];
  struct point_with_index * dist_sort_array =
    xmalloc(recv_count * sizeof(*dist_sort_array));

  // redistribute point information
  MPI_Datatype point_with_index_dt =
    get_point_with_index_mpi_datatype(comm);
  yac_mpi_call(MPI_Type_commit(&point_with_index_dt), comm);
  yac_alltoallv_p2p(
    sort_array, sendcounts, sdispls+1,
    dist_sort_array, recvcounts, rdispls,
    sizeof(*sort_array), point_with_index_dt, comm);

  for (size_t i = 0; i < recv_count; ++i) dist_sort_array[i].i = i;

  // sort received point data based on coordinates
  yac_mergesort(dist_sort_array, recv_count, sizeof(*dist_sort_array),
                compare_point_with_index_coord_id);

  struct point_with_index dummy =
    {.lon = INT32_MAX, .lat = INT32_MAX};
  struct point_with_index * prev = &dummy, * curr = dist_sort_array;
  for (size_t i = 0; i < recv_count; ++i, ++curr) {
    if (compare_point_with_index_coord(prev, curr)) {

      prev = curr;

    } else {

      curr->id = prev->id;
      curr->dlon = prev->dlon;
      curr->dlat = prev->dlat;
    }
  }

  // bring array into original order
  yac_mergesort(dist_sort_array, recv_count, sizeof(*dist_sort_array),
                compare_point_with_index_i);

  // redistribute point information
  yac_alltoallv_p2p(
    dist_sort_array, recvcounts, rdispls,
    sort_array, sendcounts, sdispls+1,
    sizeof(*sort_array), point_with_index_dt, comm);
  free(dist_sort_array);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  yac_mpi_call(MPI_Type_free(&point_with_index_dt), comm);

  // bring vertex indices into original order
  for (size_t i = 0; i < nbr_vertices; ++i) {

    (*vertex_lon)[sort_array[i].i] = sort_array[i].dlon;
    (*vertex_lat)[sort_array[i].i] = sort_array[i].dlat;
    (*vertex_ids)[sort_array[i].i] = sort_array[i].id;
  }

  free(sort_array);
}

static void detect_duplicated_cells(
  int * cell_to_vertex, int * num_vertices_per_cell, int * cell_mask,
  size_t nbr_cells, size_t ** duplicated_cell_idx,
  size_t ** orig_cell_idx, size_t * nbr_duplicated_cells) {

  // count number of unmasked cells
  size_t nbr_unmasked_cells = 0;
  for (size_t i = 0; i < nbr_cells; ++i)
    if (cell_mask[i]) ++nbr_unmasked_cells;

  struct cell_vertices_with_index * sort_array =
    xmalloc(nbr_unmasked_cells * sizeof(*sort_array));

  // get all unmasked cells
  for (size_t i = 0, offset = 0, j = 0; i < nbr_cells; ++i) {

    if (cell_mask[i]) {

      sort_array[j].num_vertices = num_vertices_per_cell[i];
      sort_array[j].cell_to_vertex = cell_to_vertex + offset;
      sort_array[j].i = i;
      ++j;
    }

    offset += (size_t)(num_vertices_per_cell[i]);
  }

  // sort cells
  yac_mergesort(sort_array, nbr_unmasked_cells, sizeof(*sort_array),
                compare_cell_vertices_with_index);

  // count number of duplicated cells
  *nbr_duplicated_cells = 0;
  for (size_t i = 1; i < nbr_unmasked_cells; ++i)
    if (!compare_cell_vertices_with_index(sort_array + (i-1), sort_array + i))
      ++*nbr_duplicated_cells;

  // get duplicated cells
  *duplicated_cell_idx =
    xmalloc(*nbr_duplicated_cells * sizeof(**duplicated_cell_idx));
  *orig_cell_idx =
    xmalloc(*nbr_duplicated_cells * sizeof(**orig_cell_idx));
  struct cell_vertices_with_index *prev = sort_array, *curr = sort_array + 1;
  for (size_t i = 1, j = 0; i < nbr_unmasked_cells; ++i, ++curr) {
    if (compare_cell_vertices_with_index(prev, curr)) {
      prev = curr;
    } else {
      (*duplicated_cell_idx)[j] = (size_t)(curr->i);
      (*orig_cell_idx)[j] = (size_t)(prev->i);
      ++j;
    }
  }

  free(sort_array);
}

static MPI_Datatype get_cell_to_vertex_ids_mpi_datatype(
  MPI_Comm comm, size_t count) {

  MPI_Datatype cell_to_vertex_ids_dt;
  yac_mpi_call(
    MPI_Type_contiguous(
      (int)count, yac_int_dt, &cell_to_vertex_ids_dt), comm);
  return cell_to_vertex_ids_dt;
}

static void detect_duplicated_cells_parallel(
  int * cell_to_vertex, int * num_vertices_per_cell, int * cell_mask,
  yac_int * cell_ids, size_t num_cells, yac_int * vertex_ids,
  MPI_Comm comm, size_t ** duplicated_cell_idx_,
  yac_int ** orig_cell_ids_, size_t * num_duplicated_cells_) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // determine the local maximum number of vertices per cell and count
  // number of unmasked cells
  size_t max_num_vertices_per_cell = 0;
  size_t local_num_cells = 0;
  for (size_t i = 0; i < num_cells; ++i) {
    if (cell_mask[i]) {
      if ((size_t)(num_vertices_per_cell[i]) > max_num_vertices_per_cell)
        max_num_vertices_per_cell = (size_t)(num_vertices_per_cell[i]);
      local_num_cells++;
    }
  }

  // determine the global maximum number of vertices per cell
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &max_num_vertices_per_cell, 1, YAC_MPI_SIZE_T, MPI_MAX,
      comm), comm);

  // generate local cell to vertex ids array
  yac_int * local_cell_ids =
    xmalloc(local_num_cells * (1 + max_num_vertices_per_cell) *
            sizeof(*local_cell_ids));
  yac_int * local_cell_to_vertex_ids = local_cell_ids + local_num_cells;

  // fill local cell to vertex ids array
  yac_int * curr_local_cell_ids = local_cell_ids;
  yac_int * curr_local_cell_to_vertex_ids = local_cell_to_vertex_ids;
  int * curr_cell_to_vertex = cell_to_vertex;
  for (size_t i = 0; i < num_cells; ++i) {

    size_t curr_num_vertices_per_cell = (size_t)(num_vertices_per_cell[i]);
    if (cell_mask[i]) {

      *curr_local_cell_ids = cell_ids[i];

      // get global vertex ids for current cell
      size_t j = 0;
      for (; j < curr_num_vertices_per_cell; ++j)
        curr_local_cell_to_vertex_ids[j] =
          vertex_ids[curr_cell_to_vertex[j]];

      // sort global vertex ids
      qsort(curr_local_cell_to_vertex_ids, curr_num_vertices_per_cell,
            sizeof(*curr_local_cell_to_vertex_ids), compare_yac_int);

      // fill remaining entries
      for (; j < max_num_vertices_per_cell; ++j)
        curr_local_cell_to_vertex_ids[j] = XT_INT_MAX;

      curr_local_cell_ids++;
      curr_local_cell_to_vertex_ids += max_num_vertices_per_cell;
    }
    curr_cell_to_vertex += curr_num_vertices_per_cell;
  }

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  // determine distributed owners for all cells
  int * dist_cell_ranks = xmalloc(local_num_cells * sizeof(*dist_cell_ranks));
  for (size_t i = 0; i < local_num_cells; ++i) {
    int rank =
      compute_bucket(
        local_cell_to_vertex_ids + i * max_num_vertices_per_cell,
        max_num_vertices_per_cell * sizeof(*local_cell_to_vertex_ids),
        comm_size);
    sendcounts[rank]++;
    dist_cell_ranks[i] = rank;
  }

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t dist_num_cells = recvcounts[comm_size-1] + rdispls[comm_size-1];
  yac_int * yac_int_buffer =
    xmalloc((local_num_cells + dist_num_cells) *
            (1 + max_num_vertices_per_cell) * sizeof(*yac_int_buffer));
  yac_int * dist_cell_ids = yac_int_buffer;
  yac_int * dist_cell_vertex_ids = yac_int_buffer + dist_num_cells;
  yac_int * send_local_cell_ids =
    yac_int_buffer + dist_num_cells * (1 + max_num_vertices_per_cell);
  yac_int * send_local_cell_to_vertex_ids =
    yac_int_buffer + local_num_cells +
    dist_num_cells * (1 + max_num_vertices_per_cell);

  // pack send data
  for (size_t i = 0; i < local_num_cells; ++i) {
    size_t pos = sdispls[dist_cell_ranks[i] + 1]++;
    send_local_cell_ids[pos] = local_cell_ids[i];
    memcpy(
      send_local_cell_to_vertex_ids + pos * max_num_vertices_per_cell,
      local_cell_to_vertex_ids + i * max_num_vertices_per_cell,
      max_num_vertices_per_cell * sizeof(*send_local_cell_to_vertex_ids));
  }
  local_cell_ids =
    xrealloc(local_cell_ids, local_num_cells * sizeof(*local_cell_ids));

  // redistribute cell ids and cell vertex ids
  MPI_Datatype cell_to_vertex_ids_dt =
    get_cell_to_vertex_ids_mpi_datatype(comm, max_num_vertices_per_cell);
  yac_mpi_call(MPI_Type_commit(&cell_to_vertex_ids_dt), comm);
  yac_alltoallv_yac_int_p2p(
    send_local_cell_ids, sendcounts, sdispls,
    dist_cell_ids, recvcounts, rdispls, comm);
  yac_alltoallv_p2p(
    send_local_cell_to_vertex_ids, sendcounts, sdispls,
    dist_cell_vertex_ids, recvcounts, rdispls,
    max_num_vertices_per_cell * sizeof(*dist_cell_vertex_ids),
    cell_to_vertex_ids_dt, comm);
  yac_mpi_call(MPI_Type_free(&cell_to_vertex_ids_dt), comm);

  // free some memory
  yac_int_buffer =
    xrealloc(
      yac_int_buffer,
      dist_num_cells * (1 + max_num_vertices_per_cell) *
      sizeof(*yac_int_buffer));
  dist_cell_ids = yac_int_buffer;
  dist_cell_vertex_ids = yac_int_buffer + dist_num_cells;

  // generate sortable data structure
  struct cell_vertex_ids_with_index * sort_array =
    xmalloc(dist_num_cells * sizeof(*sort_array));
  for (size_t i = 0; i < dist_num_cells; ++i) {
    sort_array[i].cell_id = dist_cell_ids[i];
    sort_array[i].vertex_ids =
      dist_cell_vertex_ids + i * max_num_vertices_per_cell;
    sort_array[i].num_vertices = max_num_vertices_per_cell;
    sort_array[i].i = i;
  }

  // sort cells by vertex ids
  yac_mergesort(
    sort_array, dist_num_cells, sizeof(*sort_array),
    compare_cell_vertex_ids_with_index_ids_id);

  // determine duplicated cells
  struct cell_vertex_ids_with_index * prev = sort_array;
  struct cell_vertex_ids_with_index * curr = sort_array + 1;
  for (size_t i = 1; i < dist_num_cells; ++i, ++curr) {

    if (compare_cell_vertex_ids_with_index_ids(prev, curr)) {
      prev = curr;
    } else {
      curr->cell_id = prev->cell_id;
    }
  }

  // pack cell ids
  yac_int_buffer =
    xrealloc(
      yac_int_buffer,
      (local_num_cells + dist_num_cells) * sizeof(*yac_int_buffer));
  yac_int * recv_local_cell_ids = yac_int_buffer;
  dist_cell_ids = yac_int_buffer + local_num_cells;
  for (size_t i = 0; i < dist_num_cells; ++i)
    dist_cell_ids[sort_array[i].i] = sort_array[i].cell_id;
  free(sort_array);

  // redistribute cell ids
  yac_alltoallv_yac_int_p2p(
    dist_cell_ids, recvcounts, rdispls,
    recv_local_cell_ids, sendcounts, sdispls, comm);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  // bring received local cell ids into original order and
  // count number of duplicated cells
  size_t num_duplicated_cells = 0;
  for (size_t i = 0; i < local_num_cells; ++i) {

    size_t pos = sdispls[dist_cell_ranks[i]]++;

    if (local_cell_ids[i] != recv_local_cell_ids[pos]) {
      num_duplicated_cells++;
      local_cell_ids[i] = recv_local_cell_ids[pos];
    } else {
      local_cell_ids[i] = XT_INT_MAX;
    }
  }
  free(yac_int_buffer);
  free(dist_cell_ranks);

  size_t * duplicated_cell_idx =
    xmalloc(num_duplicated_cells * sizeof(*duplicated_cell_idx));
  yac_int * orig_cell_ids =
    xmalloc(num_duplicated_cells * sizeof(*orig_cell_ids));

  // get duplicated cells
  for (size_t i = 0, j = 0, k = 0; i < num_cells; ++i) {

    if (cell_mask[i]) {

      if (local_cell_ids[j] != XT_INT_MAX) {
        duplicated_cell_idx[k] = i;
        orig_cell_ids[k] = local_cell_ids[j];
        ++k;
      }
      ++j;
    }
  }
  free(local_cell_ids);

  *duplicated_cell_idx_ = duplicated_cell_idx;
  *orig_cell_ids_ = orig_cell_ids;
  *num_duplicated_cells_ = num_duplicated_cells;
}

#else

struct yac_basic_grid_data yac_read_scrip_basic_grid_data(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value, int use_ll_edges) {

   UNUSED(grid_filename);
   UNUSED(mask_filename);
   UNUSED(grid_name);
   UNUSED(valid_mask_value);
   UNUSED(use_ll_edges);
   die(
     "ERROR(yac_read_scrip_basic_grid_data): "
     "YAC is built without the NetCDF support");

   return
    yac_generate_basic_grid_data_reg_2d(
      (size_t[]){0,0}, (int[]){0,0}, NULL, NULL);
}

struct yac_basic_grid * yac_read_scrip_basic_grid(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value, char const * name,
  int use_ll_edges, size_t * cell_coord_idx,
  size_t ** duplicated_cell_idx, yac_int ** orig_cell_global_ids,
  size_t * nbr_duplicated_cells) {

   UNUSED(grid_filename);
   UNUSED(mask_filename);
   UNUSED(grid_name);
   UNUSED(valid_mask_value);
   UNUSED(name);
   UNUSED(use_ll_edges);
   UNUSED(cell_coord_idx);
   UNUSED(duplicated_cell_idx);
   UNUSED(orig_cell_global_ids);
   UNUSED(nbr_duplicated_cells);
   die(
     "ERROR(yac_read_scrip_basic_grid): "
     "YAC is built without the NetCDF support");

  return NULL;
}

struct yac_basic_grid * yac_read_scrip_basic_grid_parallel(
  char const * grid_filename, char const * mask_filename,
  MPI_Comm comm, char const * grid_name, int valid_mask_value,
  char const * name, int use_ll_edges, size_t * cell_coord_idx,
  size_t ** duplicated_cell_idx, yac_int ** orig_cell_global_ids,
  size_t * nbr_duplicated_cells) {

   UNUSED(grid_filename);
   UNUSED(mask_filename);
   UNUSED(comm);
   UNUSED(grid_name);
   UNUSED(valid_mask_value);
   UNUSED(name);
   UNUSED(use_ll_edges);
   UNUSED(cell_coord_idx);
   UNUSED(duplicated_cell_idx);
   UNUSED(orig_cell_global_ids);
   UNUSED(nbr_duplicated_cells);
   die(
     "ERROR(yac_read_scrip_basic_grid_parallel): "
     "YAC is built without the NetCDF support");

  return NULL;
}

void yac_read_scrip_grid_information(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value,
  size_t * num_vertices, size_t * num_cells, int ** num_vertices_per_cell,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells,
  int ** cell_to_vertex, int ** cell_core_mask, size_t ** duplicated_cell_idx,
  size_t ** orig_cell_idx, size_t * nbr_duplicated_cells) {

   UNUSED(grid_filename);
   UNUSED(mask_filename);
   UNUSED(grid_name);
   UNUSED(valid_mask_value);
   UNUSED(num_vertices);
   UNUSED(num_cells);
   UNUSED(num_vertices_per_cell);
   UNUSED(x_vertices);
   UNUSED(y_vertices);
   UNUSED(x_cells);
   UNUSED(y_cells);
   UNUSED(cell_to_vertex);
   UNUSED(cell_core_mask);
   UNUSED(duplicated_cell_idx);
   UNUSED(orig_cell_idx);
   UNUSED(nbr_duplicated_cells);
   die(
     "ERROR(yac_read_scrip_grid_information): "
     "YAC is built without the NetCDF support");
}

#endif // YAC_NETCDF_ENABLED
