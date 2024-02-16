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

#include <netcdf.h>

#include "grid.h"
#include "geometry.h"
#include "read_fesom_grid.h"
#include "utils.h"
#include "io_utils.h"

static void get_fesom_vertices ( int ncid,
                                 double **vertex_lon,
                                 double **vertex_lat,
                                 int * nbr_cells,
                                 int * nbr_vertices_per_cell);

static void get_fesom_cell_center ( int ncid,
                                    double **cell_lon,
                                    double **cell_lat,
                                    size_t *nbr_cells );

static void remove_duplicated_vertices(double * temp_vertex_lon,
                                       double * temp_vertex_lat,
                                       int * temp_nbr_vertices,
                                       int * old_to_new_id);

void read_fesom_grid_information(const char * filename, int * nbr_vertices,
                                 int * nbr_cells, int ** num_vertices_per_cell,
                                 int ** cell_to_vertex, double ** x_vertices,
                                 double ** y_vertices, double ** x_cells,
                                 double ** y_cells) {

   int ncid;

   /* Open file */

   yac_nc_open(filename, NC_NOWRITE, &ncid);

   /* Get vertex longitudes and latitudes of cells and
   relations between vertices and cells*/

   {
      double * temp_vertex_lon;
      double * temp_vertex_lat;
      int nbr_vertices_per_cell;
      int temp_nbr_vertices;

      get_fesom_vertices (ncid, &temp_vertex_lon, &temp_vertex_lat,
                          nbr_cells, &nbr_vertices_per_cell);

      temp_nbr_vertices = *nbr_cells * nbr_vertices_per_cell;

      int * old_to_new_id = xmalloc(temp_nbr_vertices * sizeof(*old_to_new_id));

      remove_duplicated_vertices(temp_vertex_lon, temp_vertex_lat,
                                 &temp_nbr_vertices, old_to_new_id);

      *nbr_vertices = temp_nbr_vertices;

      *x_vertices = xrealloc(temp_vertex_lon,
                             temp_nbr_vertices * sizeof(**x_vertices));
      *y_vertices = xrealloc(temp_vertex_lat,
                             temp_nbr_vertices * sizeof(**y_vertices));

      *num_vertices_per_cell = xmalloc(*nbr_cells * sizeof(**num_vertices_per_cell));
      for (int i = 0; i < *nbr_cells; ++i)
         (*num_vertices_per_cell)[i] = nbr_vertices_per_cell;

      *cell_to_vertex = xmalloc(*nbr_cells * nbr_vertices_per_cell *
                               sizeof(**cell_to_vertex));

      // Unfortunately the data is only available in Fortran order
      for (int i = 0, k = 0; i < *nbr_cells; ++i)
         for (int j = 0; j < nbr_vertices_per_cell; ++j, ++k)
            (*cell_to_vertex)[k] = old_to_new_id[k];

      int new_size_cell_to_vertex = 0;

      // remove duplicated corners in each cell
      for (int i = 0; i < *nbr_cells; ++i) {
        (*cell_to_vertex)[new_size_cell_to_vertex++] =
              (*cell_to_vertex)[i * nbr_vertices_per_cell];
        for (int j = 1; j < nbr_vertices_per_cell; ++j)
          if ((*cell_to_vertex)[i * nbr_vertices_per_cell + j] !=
              (*cell_to_vertex)[i * nbr_vertices_per_cell + j - 1])
            (*cell_to_vertex)[new_size_cell_to_vertex++] =
              (*cell_to_vertex)[i * nbr_vertices_per_cell + j];
          else
            (*num_vertices_per_cell)[i]--;
      }

      *cell_to_vertex = xrealloc(*cell_to_vertex, new_size_cell_to_vertex *
                                 sizeof(**cell_to_vertex));

      free(old_to_new_id);

   }

   {
      /* Get longitudes and latitudes of cell center points */

      size_t temp_nbr_cells;

      get_fesom_cell_center (ncid, x_cells, y_cells, &temp_nbr_cells);
   }
   HANDLE_ERROR(nc_close(ncid));
}

struct basic_grid_data read_fesom_grid(char * filename) {

  int nbr_vertices;
  int nbr_cells;
  int * num_vertices_per_cell = NULL;
  int * cell_to_vertex = NULL;

  double * x_vertices = NULL;
  double * y_vertices = NULL;
  double * x_cells = NULL;
  double * y_cells = NULL;

  read_fesom_grid_information(filename, &nbr_vertices, &nbr_cells,
                             &num_vertices_per_cell, &cell_to_vertex, 
                             &x_vertices, &y_vertices,
                             &x_cells, &y_cells);

  free(x_cells);
  free(y_cells);

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_unstruct(
      (size_t)nbr_vertices, (size_t)nbr_cells, num_vertices_per_cell,
      x_vertices, y_vertices, cell_to_vertex);
  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(cell_to_vertex);

  return grid_data;
}

/* ---------------------------------------------------------------- */

static void get_fesom_vertices ( int ncid,
                                 double **vertex_lon,
                                 double **vertex_lat,
                                 int * nbr_cells,
                                 int * nbr_vertices_per_cell) {

  int glat_id;                        // Various NetCDF IDs
  int glon_id;

  int glon_ndims;                     // number of dimensions in file
  int glat_ndims;

  int glon_dimids[NC_MAX_VAR_DIMS];   // dimension NetCDF IDs
  int glat_dimids[NC_MAX_VAR_DIMS];

  nc_type glon_type;                  // variable type
  nc_type glat_type;

  int nbr_atts;                       // number of attributes

  size_t nbr_vertices;

  yac_nc_inq_varid (ncid, "lon_vertices", &glon_id);
  yac_nc_inq_varid (ncid, "lat_vertices", &glat_id);

  HANDLE_ERROR(
    nc_inq_var (ncid, glon_id, 0, &glon_type, &glon_ndims, glon_dimids, &nbr_atts));

  HANDLE_ERROR(
    nc_inq_var (ncid, glat_id, 0, &glat_type, &glat_ndims, glat_dimids, &nbr_atts));

  /* Allocate memory to read in coordinates */
  YAC_ASSERT(
    glon_type == glat_type,
    "ERROR(get_fesom_vertices): lon and lat datatypes do not match");
  YAC_ASSERT(
    glon_ndims == 2 && glat_ndims == 2,
    "ERROR(get_fesom_vertices): unsupported number of dimensions for lon or lat");

  /* get size of arrays */

  size_t dimlen_lat[2]; // dimension size for latitude
  size_t dimlen_lon[2]; // dimension size for longitude

  for (int i = 0; i < 2; ++i) {
    HANDLE_ERROR(nc_inq_dimlen(ncid, glon_dimids[i], dimlen_lon + i));
    HANDLE_ERROR(nc_inq_dimlen(ncid, glat_dimids[i], dimlen_lat + i));
    YAC_ASSERT(
      dimlen_lon[i] == dimlen_lat[i],
      "ERROR(get_fesom_vertices): mismatching dimension size for lon and lat");
  }

  nbr_vertices = dimlen_lon[0] * dimlen_lon[1];
  *nbr_cells = dimlen_lon[0];
  *nbr_vertices_per_cell = dimlen_lon[1];

  /* read coordinates and convert radians into degrees */


  *vertex_lon = (double * ) xmalloc ( nbr_vertices * sizeof ( double ) );
  *vertex_lat = (double * ) xmalloc ( nbr_vertices * sizeof ( double ) );

  HANDLE_ERROR(nc_get_var_double (ncid, glon_id, (*vertex_lon)));
  HANDLE_ERROR(nc_get_var_double (ncid, glat_id, (*vertex_lat)));

  for (size_t i = 0; i < nbr_vertices; ++i) {
    (*vertex_lon)[i] *= YAC_RAD;
    (*vertex_lat)[i] *= YAC_RAD;
  }
}

/* ---------------------------------------------------------------- */

static void get_fesom_cell_center ( int ncid,
                                    double **cell_lon,
                                    double **cell_lat,
                                    size_t *nbr_cells ) {

  int glat_id;                        // Various NetCDF IDs
  int glon_id;

  int glon_ndims;                     // number of dimensions in file
  int glat_ndims;

  int glon_dimids[NC_MAX_VAR_DIMS];   // dimension NetCDF IDs
  int glat_dimids[NC_MAX_VAR_DIMS];

  nc_type glon_type;                  // variable type
  nc_type glat_type;

  int nbr_atts;                       // number of attributes

  yac_nc_inq_varid (ncid, "lon", &glon_id);
  yac_nc_inq_varid (ncid, "lat", &glat_id);

  HANDLE_ERROR(
    nc_inq_var (ncid, glon_id, 0, &glon_type, &glon_ndims, glon_dimids, &nbr_atts));
  HANDLE_ERROR(
    nc_inq_var (ncid, glat_id, 0, &glat_type, &glat_ndims, glat_dimids, &nbr_atts));

  /* Allocate memory to read in coordinates */
  YAC_ASSERT(
    glon_type == glat_type,
    "ERROR(get_fesom_cell_center): lon and lat datatypes do not match");
  YAC_ASSERT(
    glon_ndims == 1 && glat_ndims == 1,
    "ERROR(get_fesom_cell_center): unsupported number of dimensions for lon or lat");

  /* get size of arrays */

  size_t dimlen_lat; // dimension size for latitude
  size_t dimlen_lon; // dimension size for longitude

  HANDLE_ERROR(nc_inq_dimlen(ncid, glon_dimids[0], &dimlen_lon));
  HANDLE_ERROR(nc_inq_dimlen(ncid, glat_dimids[0], &dimlen_lat));
  YAC_ASSERT(
    dimlen_lon == dimlen_lat,
    "ERROR(get_fesom_cell_center): mismatching dimension size for lon and lat");

  *nbr_cells = dimlen_lon;

  /* read coordinates and convert radians into degrees */

  *cell_lon = (double * ) xmalloc ( *nbr_cells * sizeof ( **cell_lon ) );
  *cell_lat = (double * ) xmalloc ( *nbr_cells * sizeof ( **cell_lat ) );

  HANDLE_ERROR(nc_get_var_double (ncid, glon_id, *cell_lon));
  HANDLE_ERROR(nc_get_var_double (ncid, glat_id, *cell_lat));

  for (size_t i = 0; i < *nbr_cells; ++i) {
    (*cell_lon)[i] *= YAC_RAD;
    (*cell_lat)[i] *= YAC_RAD;
  }
}

/* ---------------------------------------------------------------- */

struct point_with_index {

  struct {
    double lon, lat;
  } p;
  int i;
};

static int compare_point_with_index(const void * a,const void * b) {

   const struct point_with_index * a_ = (struct point_with_index*)a;
   const struct point_with_index * b_ = (struct point_with_index*)b;

  int lon_diff = fabs(a_->p.lon - b_->p.lon) > 1e-9;
  int lat_diff = fabs(a_->p.lat - b_->p.lat) > 1e-9;

  if (lon_diff) {

    if (a_->p.lon > b_->p.lon) return -1;
    else return 1;

  } else if (lat_diff) {

    if (a_->p.lat > b_->p.lat) return -1;
    else return 1;
  } else
    return 0;
}

static void remove_duplicated_vertices(double * temp_vertex_lon,
                                       double * temp_vertex_lat,
                                       int * temp_nbr_vertices,
                                       int * old_to_new_id) {

   struct point_with_index * sort_array =
      xmalloc(*temp_nbr_vertices * sizeof(*sort_array));

   for (int i = 0; i < *temp_nbr_vertices; ++i) {

      double curr_lon, curr_lat;

      curr_lon = ((double*)temp_vertex_lon)[i];
      curr_lat = ((double*)temp_vertex_lat)[i];

      while (curr_lon < 0.0) curr_lon += 360.0;
      while (curr_lon >= 360) curr_lon -= 360.0;

      sort_array[i].p.lon = curr_lon;
      sort_array[i].p.lat = curr_lat;

      sort_array[i].i = i;
   }

   yac_mergesort(sort_array, *temp_nbr_vertices, sizeof(*sort_array),
                 compare_point_with_index);

   old_to_new_id[sort_array[0].i] = 1;

   int last_unique_idx = sort_array[0].i;

   for (int i = 1; i < *temp_nbr_vertices; ++i) {

      if (compare_point_with_index(sort_array + i - 1, sort_array + i)) {

         old_to_new_id[sort_array[i].i] = 1;
         last_unique_idx = sort_array[i].i;

      } else {

         old_to_new_id[sort_array[i].i] = -last_unique_idx;
      }
   }

   free(sort_array);

   size_t new_nbr_vertices = 0;

   for (int i = 0; i < *temp_nbr_vertices; ++i) {

      if (old_to_new_id[i] == 1) {

         temp_vertex_lon[new_nbr_vertices] = temp_vertex_lon[i];
         temp_vertex_lat[new_nbr_vertices] = temp_vertex_lat[i];

         old_to_new_id[i] = new_nbr_vertices;

         new_nbr_vertices++;
      }
   }

   for (int i = 0; i < *temp_nbr_vertices; ++i)
      if (old_to_new_id[i] <= 0)
         old_to_new_id[i] = old_to_new_id[-old_to_new_id[i]];

   *temp_nbr_vertices = new_nbr_vertices;
}
