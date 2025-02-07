/**
 * @file toy_ocean.c - Toy model for ocean to test coupling with YAC.
 *
 * Copyright (C) 2022 DKRZ, MPI-M
 * SPDX-License-Identifier: BSD-3-Clause
 * See ./LICENSES/ for license information
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "yac_interface.h"

static void read_icon_grid(const char * filename, int * nbr_vertices,
                           int * nbr_cells, int ** num_vertices_per_cell,
                           int ** cell_to_vertex, double ** x_vertices,
                           double ** y_vertices, double ** x_cells,
                           double ** y_cells, int ** cell_mask);

/* ------------------------------------------------- */

int rdc2ocn_field_id;

int main () {

  // Initialisation of MPI
  MPI_Init (0, NULL);

  // initialise YAC
  yac_cinit ();

  // configuration file is read in by HD model

  int comp_id;
  char * comp_name = "OCEAN";

  // define ocean component
  yac_cdef_comp ( comp_name, &comp_id );

  MPI_Comm local_comm;

  // get communicator for ocean component
  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;
  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);

  if (size != 1) {
    fputs("toy_ocean: only a single toy process is supported", stderr);
    exit(EXIT_FAILURE);
  }

  char const * filename =
    "/pool/data/ICON/grids/public/mpim/0020/icon_grid_0020_R02B05_O.nc";

  int num_vertices;
  int num_cells;
  int * num_vertices_per_cell;
  int * cell_to_vertex;
  double * x_vertices, * y_vertices;
  double * x_cells, * y_cells;
  int * cell_mask;

  read_icon_grid(
    filename, &num_vertices, &num_cells, &num_vertices_per_cell,
    &cell_to_vertex, &x_vertices, &y_vertices, &x_cells, &y_cells, &cell_mask);

  // define icon grid
  int grid_id;
  yac_cdef_grid_unstruct(
    "ocean_grid", num_vertices, num_cells, num_vertices_per_cell,
    x_vertices, y_vertices, cell_to_vertex, &grid_id);

  // set global cell ids and core mask
  int * global_cell_id = malloc((size_t)num_cells * sizeof(*global_cell_id));
  for (int i = 0; i < num_cells; ++i) global_cell_id[i] = i;
  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);

  int cell_point_id;

  // define field locations (at cell centers)
  yac_cdef_points_unstruct(
    grid_id, num_cells, YAC_LOCATION_CELL, x_cells, y_cells, &cell_point_id);

  free(x_cells);
  free(y_cells);

  int coast_cell_mask_id;

  // define a costal mask
  int * coast_cell_mask = malloc(num_cells * sizeof(coast_cell_mask));
  for (int i = 0; i < num_cells; ++i) coast_cell_mask[i] = (cell_mask[i] == 1);
  yac_cdef_mask(
    grid_id, num_cells, YAC_LOCATION_CELL, coast_cell_mask,
    &coast_cell_mask_id);

  int const nlev = 1;

  // define fields
  char const * timestep = "24";
  yac_cdef_field_mask(
    "RDC2NEMO", comp_id, &cell_point_id, &coast_cell_mask_id, 1,
    nlev, timestep, YAC_TIME_UNIT_HOUR, &rdc2ocn_field_id);

  free(coast_cell_mask);
  free(cell_mask);

  free(x_vertices);
  free(y_vertices);

  // setup coupling
  yac_cenddef ( );

  // allocate field arrays
  double * river_disc  = malloc(num_cells * sizeof(*river_disc));
  for (int i = 0; i < num_cells; ++i) river_disc[i] = 0.0;

  int t = 0;
  int err;
  int info = 0;

  while ((info == 0) || (info == 1)) {

    { // receive river discharge
      double *collection_data[1];

      collection_data[0] = river_disc;
      yac_cget(rdc2ocn_field_id, 1, collection_data, &info, &err);
    }

    if (info == 3) fputs("toy_ocean received last get\n", stdout);

    t++;
  }

  free(river_disc);

  yac_cfinalize();

  MPI_Finalize();

  return EXIT_SUCCESS;
}

static void inline handle_error(int status) {
  if (status != NC_NOERR) {
    const char *err_string;
    err_string = nc_strerror(status);
    printf ("%s \n", err_string);
    exit(EXIT_FAILURE);
  }
}

static size_t check_dimension(int ncid, int varids[2]) {

  for (int i = 0; i < 2; ++i) {
    int ndims;
    handle_error(nc_inq_varndims(ncid, varids[i], &ndims));
    if(ndims != 1) {
      fputs(
        "ERROR(check_dimension): coordinate array has more than one dimension",
        stderr);
      exit(EXIT_FAILURE);
    }
  }

  int dimids[2];
  for (int i = 0; i < 2; ++i)
    handle_error(nc_inq_vardimid(ncid, varids[i], &(dimids[i])));

  if(dimids[0] != dimids[1]) {
    fputs("ERROR(check_dimension): "
          "lon lat coordinate arrays have differing dimensions", stderr);
    exit(EXIT_FAILURE);
  }

  size_t dimlen;
  handle_error(nc_inq_dimlen(ncid, dimids[0], &dimlen));
  return dimlen;
}

static int check_coord_units(int ncid, int varid) {

  int is_degree = 0;
  nc_type att_type;
  size_t att_len;
  int status = nc_inq_att(ncid, varid, "units", &att_type, &att_len);
  // if the status is not "attribute not found"
  if (status != NC_ENOTATT) {
    handle_error(status);
    // if the attribute is not a string or too long
    if ((att_type != NC_CHAR) || (att_len > 8)) {
      fputs("ERROR(read_coord): invalid units type or len", stderr);
      exit(EXIT_FAILURE);
    }
    char units[8];
    memset(units, 0, 8 * sizeof(units[0]));
    handle_error(nc_get_att_text(ncid, varid, "units", units));
    is_degree = !strcmp(units, "degree");
    if (!is_degree && strcmp(units, "radian")) {
      fputs("ERROR(read_coord): unsupported units type", stderr);
      exit(EXIT_FAILURE);
    }
  }
  return is_degree;
}

static double * read_coord(int ncid, int varid, size_t varlen) {

  int is_degree = check_coord_units(ncid, varid);

  double * coord = malloc(varlen * sizeof(*coord));
  handle_error(nc_get_var_double (ncid, varid, coord));

  // convert to radiant if necessary
  if (is_degree)
    for (size_t i = 0; i < varlen; ++i)
      coord[i] *= 0.017453292519943295769;

  return coord;
}

static void read_coords(
  int ncid, char const * lon_name, char const * lat_name,
  double ** lon, double ** lat, size_t * len) {

  int vlonid, vlatid;
  handle_error(nc_inq_varid(ncid, lon_name, &vlonid));
  handle_error(nc_inq_varid(ncid, lat_name, &vlatid));

  size_t varlen = (*len = check_dimension(ncid, (int[]){vlonid, vlatid}));
  *lon = read_coord(ncid, vlonid, varlen);
  *lat = read_coord(ncid, vlatid, varlen);
}

static int * get_icon_connect(int ncid, size_t nbr_cells) {

  // get variable id
  int conn_id;
  handle_error(nc_inq_varid(ncid, "vertex_of_cell", &conn_id));

  // check number of dimension (has to be 2)
  int ndims;
  handle_error(nc_inq_varndims(ncid, conn_id, &ndims));
  if (ndims != 2) {
    fputs("ERROR(get_icon_connect): "
          "connectivity array has invalid number of dimensions", stderr);
    exit(EXIT_FAILURE);
  }

  // get ids of dimensions
  int dimids[2];
  handle_error(nc_inq_vardimid(ncid, conn_id, dimids));

  // check size of dimensions
  size_t dimlen;
  handle_error(nc_inq_dimlen(ncid, dimids[0], &dimlen));
  if (dimlen != 3) {
    fputs("ERROR(get_icon_connect): invalid size of first dimension of "
          "connectivity array (has to be 3)", stderr);
    exit(EXIT_FAILURE);
  }
  handle_error(nc_inq_dimlen(ncid, dimids[1], &dimlen));
  if (dimlen != nbr_cells) {
    fputs("ERROR(get_icon_connect): invalid size of second dimension of "
          "connectivity array (has to be nbr_cells)", stderr);
    exit(EXIT_FAILURE);
  }

  // get and return connectivity array
  int * vertex_of_cell = malloc(3 * nbr_cells * sizeof(*vertex_of_cell));
  handle_error(nc_get_var_int (ncid, conn_id, vertex_of_cell));
  return vertex_of_cell;
}

static int * get_icon_cell_mask ( int ncid, size_t nbr_cells ) {

  // get variable id
  int mask_id;
  handle_error(nc_inq_varid (ncid, "cell_sea_land_mask", &mask_id));

  // check number of dimension (has to be 1)
  int ndims;
  handle_error(nc_inq_varndims(ncid, mask_id, &ndims));
  if(ndims != 1) {
    fputs(
      "ERROR(get_icon_cell_mask): mask array has more than one dimension",
      stderr);
    exit(EXIT_FAILURE);
  }

  // get id of dimension
  int dimid;
  handle_error(nc_inq_vardimid(ncid, mask_id, &dimid));

  // check size of mask (has to be equal to nbr_cells)
  size_t dimlen;
  handle_error(nc_inq_dimlen(ncid, dimid, &dimlen));
  if(dimlen != nbr_cells) {
    fputs("ERROR(get_icon_cell_mask): invalid size of mask array", stderr);
    exit(EXIT_FAILURE);
  }

  // check mask type (has to be NC_INT)
  nc_type mask_type;
  handle_error(nc_inq_vartype(ncid, mask_id, &mask_type));
  if(mask_type != NC_INT) {
    fputs("ERROR(get_icon_cell_mask): invalid mask type", stderr);
    exit(EXIT_FAILURE);
  }

  // get and return mask
  int * cell_mask = malloc(nbr_cells * sizeof(*cell_mask));
  handle_error(nc_get_var_int (ncid, mask_id, cell_mask));
  return cell_mask;
}

static void read_icon_grid(const char * filename, int * nbr_vertices,
                           int * nbr_cells, int ** num_vertices_per_cell,
                           int ** cell_to_vertex, double ** x_vertices,
                           double ** y_vertices, double ** x_cells,
                           double ** y_cells, int ** cell_mask) {

   /* Open file */
   int ncid;
   handle_error(nc_open(filename, NC_NOWRITE, &ncid));

   /* Get vertex longitudes and latitudes of cells */
   size_t nbr_vertices_;
   read_coords(ncid, "vlon", "vlat", x_vertices, y_vertices, &nbr_vertices_);
   *nbr_vertices = (int)nbr_vertices_;

   /* Get cell center longitudes and latitudes of cells */
   size_t nbr_cells_;
   read_coords(ncid, "clon", "clat", x_cells, y_cells, &nbr_cells_);
   *nbr_cells = (int)nbr_cells_;

   /* Get mask of cells */
   *cell_mask = get_icon_cell_mask ( ncid, *nbr_cells );

   /* Get relations between vertices and cells */
   int * vertex_of_cell = get_icon_connect(ncid, *nbr_cells);

   /* Close file */
   handle_error(nc_close(ncid));

   //-------------------------------------------------------------------------//

   *num_vertices_per_cell = malloc(*nbr_cells * sizeof(**num_vertices_per_cell));
   for (int i = 0; i < *nbr_cells; (*num_vertices_per_cell)[i++] = 3);

   *cell_to_vertex = malloc(*nbr_cells * 3 * sizeof(**cell_to_vertex));

   // Unfortunately the data is only available in Fortran order
   for (int i = 0; i < *nbr_cells; ++i) {

      (*cell_to_vertex)[3*i+0] = vertex_of_cell[i+0*(*nbr_cells)] - 1;
      (*cell_to_vertex)[3*i+1] = vertex_of_cell[i+1*(*nbr_cells)] - 1;
      (*cell_to_vertex)[3*i+2] = vertex_of_cell[i+2*(*nbr_cells)] - 1;
   }

   free(vertex_of_cell);
}
