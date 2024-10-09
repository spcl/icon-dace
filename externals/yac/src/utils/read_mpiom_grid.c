// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "read_mpiom_grid.h"
#include "geometry.h"
#include "utils_common.h"
#include "io_utils.h"

#ifdef YAC_NETCDF_ENABLED

#include <netcdf.h>


static void get_mpiom_vertices(int ncid,
                               double **vertex_lon, // longitude of vertex coordinates
                               double **vertex_lat, // latitude of vertex coordinates
                               size_t *nbr_vertices);

static void get_mpiom_cell_center(int ncid,
                                  double **cell_lon, // longitude of cell coordinates
                                  double **cell_lat, // latitude of cell coordinates
                                  size_t *nbr_cells);

static void get_mpiom_cell_mask(int ncid,
                                int **cell_mask, // data for mask
                                size_t *nbr_cells );

static void remove_duplicated_coords(double * temp_coords_lon,
                                     double * temp_coords_lat,
                                     int * temp_coords_mask,
                                     size_t * temp_nbr_coords,
                                     int * old_to_new_id);

void yac_read_mpiom_grid_information(const char * filename, int * nbr_vertices,
                                     int * nbr_cells, int ** num_vertices_per_cell,
                                     int ** cell_to_vertex, double ** x_vertices,
                                     double ** y_vertices, double ** x_cells,
                                     double ** y_cells, int ** cell_mask) {

  int ncid;

  /* Open file */

  yac_nc_open(filename, NC_NOWRITE, &ncid);

  /* Get vertex longitudes and latitudes of cells and
  relations between vertices and cells*/

  {
    double * temp_vertex_lon;
    double * temp_vertex_lat;
    size_t temp_nbr_vertices;

    get_mpiom_vertices(
      ncid, &temp_vertex_lon, &temp_vertex_lat, &temp_nbr_vertices);

    *nbr_cells = (int)(temp_nbr_vertices / 4);

    int * old_to_new_id = xmalloc(temp_nbr_vertices * sizeof(*old_to_new_id));

    remove_duplicated_coords(
      temp_vertex_lon, temp_vertex_lat, NULL,
      &temp_nbr_vertices, old_to_new_id);

    *nbr_vertices = temp_nbr_vertices;

    *x_vertices =
      xrealloc(temp_vertex_lon, temp_nbr_vertices * sizeof(*temp_vertex_lon));
    *y_vertices =
      xrealloc(temp_vertex_lat, temp_nbr_vertices * sizeof(*temp_vertex_lat));

    for (size_t i = 0; i < temp_nbr_vertices; ++i) {
      (*x_vertices)[i] *= YAC_RAD;
      (*y_vertices)[i] *= YAC_RAD;
    }

    *num_vertices_per_cell = xmalloc(*nbr_cells * sizeof(**num_vertices_per_cell));
    for (int i = 0; i < *nbr_cells; ++i) (*num_vertices_per_cell)[i] = 4;

    *cell_to_vertex = xmalloc(*nbr_cells * 4 * sizeof(**cell_to_vertex));

    // Unfortunately the data is only available in Fortran order
    for (int i = 0; i < *nbr_cells; ++i)
      for (int j = 0; j < 4; ++j)
        (*cell_to_vertex)[4*i+j] = old_to_new_id[4*i+j];
    free(old_to_new_id);
  }

  {
    /* Get longitudes and latitudes of cell center points */

    double * temp_cell_lon;
    double * temp_cell_lat;
    int *  temp_cell_mask;

    size_t temp_nbr_cells;
    size_t temp_nbr_masks;

    get_mpiom_cell_center(
      ncid, &temp_cell_lon, &temp_cell_lat, &temp_nbr_cells);

    get_mpiom_cell_mask(ncid, &temp_cell_mask, &temp_nbr_masks);

    YAC_ASSERT(
      temp_nbr_cells == temp_nbr_masks,
      "ERROR(yac_read_mpiom_grid_information): "
      "missmatch between cell count and mask value count")

    int * old_to_new_id = xmalloc(temp_nbr_cells * sizeof(*old_to_new_id));

    remove_duplicated_coords(
      temp_cell_lon, temp_cell_lat, temp_cell_mask,
      &temp_nbr_cells, old_to_new_id);

    *nbr_cells = temp_nbr_cells;

    *x_cells = xrealloc(temp_cell_lon, temp_nbr_cells * sizeof(**x_cells));
    *y_cells = xrealloc(temp_cell_lat, temp_nbr_cells * sizeof(**y_cells));

    for (size_t i = 0; i < temp_nbr_cells; ++i) {
      (*x_cells)[i] *= YAC_RAD;
      (*y_cells)[i] *= YAC_RAD;
    }

    *cell_mask = xrealloc(temp_cell_mask,temp_nbr_cells * sizeof(int));
    for (size_t i = 0; i < temp_nbr_cells; ++i)
      (*cell_mask)[i] = !(*cell_mask)[i];

    free(old_to_new_id);
  }
  YAC_HANDLE_ERROR(nc_close(ncid));
}

void yac_read_part_mpiom_grid_information(
  const char * filename, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, double ** x_vertices,
  double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size) {

  // read the global grid
  yac_read_mpiom_grid_information(filename, num_vertices, num_cells,
                                  num_vertices_per_cell, cell_to_vertex,
                                  x_vertices, y_vertices, x_cells, y_cells, cell_mask);

  // determine local rank
  int num_cells_per_process = (*num_cells + size - 1)/size;
  int local_start = rank * num_cells_per_process;
  int num_local_cells = MIN(num_cells_per_process,
                            *num_cells - local_start);

  // mask for required vertices and cells
  int * required_vertices = xcalloc(*num_vertices, sizeof(*required_vertices));
  int * required_cells = xcalloc(*num_cells, sizeof(*required_cells));

  int offset = 0;
  for (int i = 0; i < local_start; ++i)
    offset += (*num_vertices_per_cell)[i];

  // mark all local cells and their vertices as required
  for (int i = local_start; i < local_start + num_local_cells; ++i) {

    required_cells[i] = 2;
    for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j)
      required_vertices[(*cell_to_vertex)[offset+j]] = 2;

    offset += (*num_vertices_per_cell)[i];
  }

  // mark all halo cells as required
  offset = 0;
  for (int i = 0; i < *num_cells; ++i) {

    if (!required_cells[i]) {

      for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j) {

        if (required_vertices[(*cell_to_vertex)[offset+j]]) {

          required_cells[i] = 1;
          break;
        }
      }
    }
    offset += (*num_vertices_per_cell)[i];
  }

  // mark all halo vertices as required
  offset = 0;
  for (int i = 0; i < *num_cells; ++i) {

    if (required_cells[i] == 1) {

      for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j) {

        if (!required_vertices[(*cell_to_vertex)[offset+j]]) {

          required_vertices[(*cell_to_vertex)[offset+j]] = 1;
        }
      }
    }

    offset += (*num_vertices_per_cell)[i];
  }

  // count the number of cells and vertices
  int part_num_vertices = 0;
  int part_num_cells = 0;
  for (int i = 0; i < *num_vertices; ++i)
    if (required_vertices[i])
      part_num_vertices++;
  for (int i = 0; i < *num_cells; ++i)
    if(required_cells[i])
      part_num_cells++;

  *global_cell_id = xmalloc(part_num_cells * sizeof(**global_cell_id));
  *cell_core_mask = xmalloc(part_num_cells * sizeof(**cell_core_mask));
  *global_corner_id = xmalloc(part_num_vertices * sizeof(**global_corner_id));
  *corner_core_mask = xmalloc(part_num_vertices * sizeof(**corner_core_mask));

  // generate final vertex data
  part_num_vertices = 0;
  int * global_to_local_vertex = xmalloc(*num_vertices * sizeof(*global_to_local_vertex));
  for (int i = 0; i < *num_vertices; ++i) {

    if (required_vertices[i]) {

      (*global_corner_id)[part_num_vertices] = i;
      (*corner_core_mask)[part_num_vertices] = required_vertices[i] == 2;
      (*x_vertices)[part_num_vertices] = (*x_vertices)[i];
      (*y_vertices)[part_num_vertices] = (*y_vertices)[i];
      global_to_local_vertex[i] = part_num_vertices;
      part_num_vertices++;
    }
  }

  *x_vertices = xrealloc(*x_vertices, part_num_vertices * sizeof(**x_vertices));
  *y_vertices = xrealloc(*y_vertices, part_num_vertices * sizeof(**y_vertices));
  *num_vertices = part_num_vertices;
  free(required_vertices);

  // generate final cell data
  int num_cell_vertex_dependencies = 0;
  part_num_cells = 0;
  offset = 0;
  for (int i = 0; i < *num_cells; ++i) {

    if (required_cells[i]) {

      (*global_cell_id)[part_num_cells] = i;
      (*cell_core_mask)[part_num_cells] = required_cells[i] == 2;
      (*x_cells)[part_num_cells] = (*x_cells)[i];
      (*y_cells)[part_num_cells] = (*y_cells)[i];
      (*cell_mask)[part_num_cells] = (*cell_mask)[i];

      for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j)
        (*cell_to_vertex)[num_cell_vertex_dependencies++] =
          global_to_local_vertex[(*cell_to_vertex)[offset+j]];

      (*num_vertices_per_cell)[part_num_cells] = (*num_vertices_per_cell)[i];

      part_num_cells++;
    }

    offset += (*num_vertices_per_cell)[i];
  }

  *x_cells = xrealloc(*x_cells, part_num_cells * sizeof(**x_cells));
  *y_cells = xrealloc(*y_cells, part_num_cells * sizeof(**y_cells));
  *cell_mask = xrealloc(*cell_mask, part_num_cells * sizeof(**cell_mask));

  *num_vertices_per_cell = xrealloc(*num_vertices_per_cell, part_num_cells *
                                    sizeof(**num_vertices_per_cell));
  *cell_to_vertex = xrealloc(*cell_to_vertex, num_cell_vertex_dependencies *
                            sizeof(**cell_to_vertex));
  *num_cells = part_num_cells;
  free(required_cells);
  free(global_to_local_vertex);
}


struct yac_basic_grid_data yac_read_mpiom_basic_grid_data(
  char const * filename) {

  int nbr_vertices;
  int nbr_cells;
  int * num_vertices_per_cell = NULL;
  int * cell_to_vertex = NULL;

  int * cell_mask = NULL;

  double * x_vertices = NULL;
  double * y_vertices = NULL;
  double * x_cells = NULL;
  double * y_cells = NULL;

  yac_read_mpiom_grid_information(filename, &nbr_vertices, &nbr_cells,
                                  &num_vertices_per_cell, &cell_to_vertex,
                                  &x_vertices, &y_vertices,
                                  &x_cells, &y_cells, &cell_mask);

  struct yac_basic_grid_data grid =
    yac_generate_basic_grid_data_unstruct(
      (size_t)nbr_vertices, (size_t)nbr_cells, num_vertices_per_cell,
      x_vertices, y_vertices, cell_to_vertex);
  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(cell_to_vertex);
  free(x_cells);
  free(y_cells);
  free(cell_mask);

  return grid;
}

struct yac_basic_grid * yac_read_mpiom_basic_grid(
  char const * filename, char const * gridname) {

  return
    yac_basic_grid_new(
      gridname, yac_read_mpiom_basic_grid_data(filename));
}

/* ---------------------------------------------------------------- */

static void get_mpiom_vertices(
  int ncid, double ** vertex_lon, double ** vertex_lat,
  size_t * nbr_vertices) {

  int glon_id;
  int glat_id;

  // get variable ids
  if (nc_inq_varid(ncid, "lon_bounds", &glon_id) != NC_NOERR)
    yac_nc_inq_varid(ncid, "lon_bnds", &glon_id);

  if (nc_inq_varid(ncid, "lat_bounds", &glat_id) != NC_NOERR)
    yac_nc_inq_varid(ncid, "lat_bnds", &glat_id);

  int glon_ndims;
  int glat_ndims;

  // get number of dimensions
  YAC_HANDLE_ERROR(nc_inq_varndims(ncid, glon_id, &glon_ndims));
  YAC_HANDLE_ERROR(nc_inq_varndims(ncid, glat_id, &glat_ndims));
  YAC_ASSERT(
    (glon_ndims == 3) && (glat_ndims == 3),
    "ERROR(get_mpiom_vertices): "
    "number of dimensions for lon and lat does not match")

  int glon_dimids[3];
  int glat_dimids[3];

  // get dimension ids
  YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, glon_id, glon_dimids));
  YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, glat_id, glat_dimids));

  size_t dimlen_lat[3];
  size_t dimlen_lon[3];

  // get size of arrays
  for (int i = 0; i < 3; ++i) {
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, glon_dimids[i], &dimlen_lon[i]));
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, glat_dimids[i], &dimlen_lat[i]));
    YAC_ASSERT(
      dimlen_lon[i] == dimlen_lat[i],
      "ERROR(get_mpiom_vertices): dimlen_lon[i] != dimlen_lat[i]")
  }
  YAC_ASSERT(
    dimlen_lon[2] == 4, "ERROR(get_mpiom_vertices): dimlen_lon[2] has to be 4")

  *nbr_vertices = dimlen_lon[0] * dimlen_lon[1] * dimlen_lon[2];

  *vertex_lon = xmalloc(*nbr_vertices * sizeof(**vertex_lon));
  *vertex_lat = xmalloc(*nbr_vertices * sizeof(**vertex_lat));

  YAC_HANDLE_ERROR(nc_get_var_double(ncid, glon_id, *vertex_lon));
  YAC_HANDLE_ERROR(nc_get_var_double(ncid, glat_id, *vertex_lat));
}

/* ---------------------------------------------------------------- */

static void get_mpiom_cell_center ( int ncid,
                                    double **cell_lon,
                                    double **cell_lat,
                                    size_t *nbr_cells ) {

  int glon_id;
  int glat_id;

  // get variable ids
  yac_nc_inq_varid (ncid, "lon", &glon_id);
  yac_nc_inq_varid (ncid, "lat", &glat_id);

  int glon_ndims;
  int glat_ndims;

  // get number of dimensions
  YAC_HANDLE_ERROR(nc_inq_varndims(ncid, glon_id, &glon_ndims));
  YAC_HANDLE_ERROR(nc_inq_varndims(ncid, glat_id, &glat_ndims));
  YAC_ASSERT(
    (glon_ndims == 2) && (glat_ndims == 2),
    "ERROR(get_mpiom_cell_center): (glon_ndims != 2) || (glat_ndims != 2)")

  int glon_dimids[2];
  int glat_dimids[2];

  // get dimension ids
  YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, glon_id, glon_dimids));
  YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, glat_id, glat_dimids));

  size_t dimlen_lon[2];
  size_t dimlen_lat[2];

  // get size of arrays
  for (int i = 0; i < 2; ++i) {
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, glon_dimids[i], &dimlen_lon[i]));
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, glat_dimids[i], &dimlen_lat[i]));
    YAC_ASSERT(
      dimlen_lon[i] == dimlen_lat[i],
      "ERROR(get_mpiom_cell_center): dimlen_lon[i] != dimlen_lat[i]")
  }

  *nbr_cells = dimlen_lon[0] * dimlen_lon[1];

  *cell_lon = xmalloc(*nbr_cells * sizeof(**cell_lon));
  *cell_lat = xmalloc(*nbr_cells * sizeof(**cell_lat));

  YAC_HANDLE_ERROR(nc_get_var_double(ncid, glon_id, *cell_lon));
  YAC_HANDLE_ERROR(nc_get_var_double(ncid, glat_id, *cell_lat));
}

/* ---------------------------------------------------------------- */

static void get_mpiom_cell_mask(
  int ncid, int **cell_mask, size_t *nbr_cells ) {

  int status;

  // get variable id
  int mask_id;
  status = nc_inq_varid (ncid, "var1", &mask_id);
  if (status != NC_NOERR) yac_nc_inq_varid(ncid, "zo", &mask_id);

  // get number of mask dimensions
  int mask_ndims;
  YAC_HANDLE_ERROR(nc_inq_varndims(ncid, mask_id, &mask_ndims));
  YAC_ASSERT(
    (mask_ndims == 3) || (mask_ndims == 4),
    "ERROR(get_mpiom_cell_mask): invalid number of dimensions")

  // get dimension ids
  int mask_dimids[4];
  YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, mask_id, mask_dimids));

  // get size of mask dimensions
  size_t dimlen_mask[4] = {0, 0, 0, 0};
  for (int i = 0; i < mask_ndims; ++i)
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, mask_dimids[i], &dimlen_mask[i]));

  *nbr_cells = dimlen_mask[mask_ndims - 2] * dimlen_mask[mask_ndims - 1];

  int * temp_cell_mask = xmalloc(*nbr_cells * sizeof(*temp_cell_mask));

  YAC_HANDLE_ERROR(nc_get_var_int(ncid, mask_id, temp_cell_mask));

  if (mask_ndims == 4)
    for (size_t i = 0; i < *nbr_cells; ++i)
      temp_cell_mask[i] = !(temp_cell_mask[i] == 1);

  *cell_mask = temp_cell_mask;
}

/* ---------------------------------------------------------------- */

struct point_with_index {

  int32_t lon, lat;
  int64_t i;
};

static inline int compare_point_with_index(const void * a,const void * b) {

   return ( (*(int64_t*)a < *(int64_t*)b) - (*(int64_t*)a > *(int64_t*)b) );
}

static void remove_duplicated_coords(
  double * temp_coords_lon, double * temp_coords_lat, int * temp_coords_mask,
  size_t * temp_nbr_coords, int * old_to_new_id) {

  struct point_with_index * sort_array =
    xmalloc(*temp_nbr_coords * sizeof(*sort_array));

  double const scale = (double)(2 << 20);
  int32_t const periode = (int32_t)(360.0 * scale);

  for (size_t i = 0; i < *temp_nbr_coords; ++i) {

    int32_t curr_lon, curr_lat;

    curr_lon = (int32_t)(temp_coords_lon[i] * scale);
    curr_lat = (int32_t)(temp_coords_lat[i] * scale);

    if (curr_lon <= 0)
      curr_lon = curr_lon - (((curr_lon - periode) / periode) * periode);
    else if (curr_lon > periode)
      curr_lon = curr_lon - ((curr_lon / periode) * periode);

    sort_array[i].lon = curr_lon;
    sort_array[i].lat = curr_lat;

    sort_array[i].i = i;
  }

  qsort(sort_array, *temp_nbr_coords, sizeof(*sort_array),
        compare_point_with_index);

  old_to_new_id[sort_array[0].i] = 1;

  int last_unique_idx = (int)(sort_array[0].i);

  for (size_t i = 1; i < *temp_nbr_coords; ++i) {

    if (compare_point_with_index(sort_array + i - 1, sort_array + i)) {

      old_to_new_id[sort_array[i].i] = 1;
      last_unique_idx = (int)(sort_array[i].i);

    } else {

      old_to_new_id[sort_array[i].i] = -last_unique_idx;
    }
  }

  free(sort_array);

  size_t new_nbr_coords = 0;

  for (size_t i = 0; i < *temp_nbr_coords; ++i) {

    if (old_to_new_id[i] == 1) {

      temp_coords_lon[new_nbr_coords] = temp_coords_lon[i];
      temp_coords_lat[new_nbr_coords] = temp_coords_lat[i];

      if (temp_coords_mask != NULL)
        temp_coords_mask[new_nbr_coords] = temp_coords_mask[i];

      old_to_new_id[i] = new_nbr_coords;

      new_nbr_coords++;
    }
  }

  for (size_t i = 0; i < *temp_nbr_coords; ++i)
    if (old_to_new_id[i] <= 0)
      old_to_new_id[i] = old_to_new_id[-old_to_new_id[i]];

  *temp_nbr_coords = new_nbr_coords;
}

#else

struct yac_basic_grid_data yac_read_mpiom_basic_grid_data(
  char const * filename) {

   UNUSED(filename);
   die(
     "ERROR(yac_read_mpiom_basic_grid_data): "
     "YAC is built without the NetCDF support");

   return
    yac_generate_basic_grid_data_reg_2d(
      (size_t[]){0,0}, (int[]){0,0}, NULL, NULL);
}

struct yac_basic_grid * yac_read_mpiom_basic_grid(
  char const * filename, char const * gridname) {

   UNUSED(filename);
   UNUSED(gridname);
   die(
     "ERROR(yac_read_mpiom_basic_grid): "
     "YAC is built without the NetCDF support");

   return NULL;
}

void yac_read_mpiom_grid_information(const char * filename, int * nbr_vertices,
                                     int * nbr_cells, int ** num_vertices_per_cell,
                                     int ** cell_to_vertex, double ** x_vertices,
                                     double ** y_vertices, double ** x_cells,
                                     double ** y_cells, int ** cell_mask) {

   UNUSED(filename);
   UNUSED(nbr_vertices);
   UNUSED(nbr_cells);
   UNUSED(num_vertices_per_cell);
   UNUSED(cell_to_vertex);
   UNUSED(x_vertices);
   UNUSED(y_vertices);
   UNUSED(x_cells);
   UNUSED(y_cells);
   UNUSED(cell_mask);
   die(
     "ERROR(yac_read_mpiom_grid_information): "
     "YAC is built without the NetCDF support");
}

void yac_read_part_mpiom_grid_information(
  const char * filename, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, double ** x_vertices,
  double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size) {

   UNUSED(filename);
   UNUSED(num_vertices);
   UNUSED(num_cells);
   UNUSED(num_vertices_per_cell);
   UNUSED(cell_to_vertex);
   UNUSED(x_vertices);
   UNUSED(y_vertices);
   UNUSED(x_cells);
   UNUSED(y_cells);
   UNUSED(global_cell_id);
   UNUSED(cell_mask);
   UNUSED(cell_core_mask);
   UNUSED(global_corner_id);
   UNUSED(corner_core_mask);
   UNUSED(rank);
   UNUSED(size);
   die(
     "ERROR(yac_read_part_mpiom_grid_information): "
     "YAC is built without the NetCDF support");
}

#endif // YAC_NETCDF_ENABLED
