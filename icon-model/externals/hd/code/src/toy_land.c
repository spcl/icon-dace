/**
 * @file toy_land.c - Toy model for Land to test coupling with YAC
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

/* ------------------------------------------------- */

int runoff_s_field_id;
int runoff_g_field_id;

struct forcing_file {
  int ncid;
  size_t curr_timestep;
  size_t ntimesteps;
  size_t nlon, nlat;
  int qs_varid, qsb_varid;
};

static struct forcing_file open_input_file(
  char const * filename, double ** lon, double ** lat, int ** cell_mask);
static void read_timestep(
  struct forcing_file * file, double * qs, double * qsb);
static void close_input_file(struct forcing_file file);
static void receive_HD_domain(double * HD_domain_corners);

int main () {

  // Initialisation of MPI
  MPI_Init (0, NULL);

  // initialise YAC
  yac_cinit ();

  // configuration file is read in by HD model

  int comp_id;
  char * comp_name = "LAND";

  // define land component
  yac_cdef_comp ( comp_name, &comp_id );

  // double HD_domain_corners[2][2];
  // receive_HD_domain(&HD_domain_corners[0][0]);

  MPI_Comm local_comm;

  // get communicator for land component
  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;
  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);

  if (size != 1) {
    fputs("toy_land: only a single toy process is supported", stderr);
    exit(EXIT_FAILURE);
  }

  char * filename =
    "/work/gg0302/g260122/HD/forcing/55056/55056_daily_1979.nc";

  double * x_center, * y_center;
  int * cell_mask;
  struct forcing_file file =
    open_input_file(filename, &x_center, &y_center, &cell_mask);

  int grid_id;
  int num_cells[2] = {file.nlon, file.nlat};
  int num_vertices[2] = {num_cells[0] + 1, num_cells[1] + 1};
  int cyclic[2] = {0, 0};

  double * x_vertices = malloc(num_vertices[0] * sizeof(*x_vertices));
  double * y_vertices = malloc(num_vertices[1] * sizeof(*y_vertices));
  double dx = x_center[1] - x_center[0];
  double dy = y_center[1] - y_center[0];
  x_vertices[0] = x_center[0] - 0.5 * dx;
  y_vertices[0] = y_center[0] - 0.5 * dy;
  for (size_t i = 1; i < num_vertices[0]; ++i)
    x_vertices[i] = x_vertices[0] + dx * (double)i;
  for (size_t i = 1; i < num_vertices[1]; ++i)
    y_vertices[i] = y_vertices[0] + dy * (double)i;

  // for (size_t i = 0, k = 0; i < num_cells[1]; ++i) {
    // int within_lat_bounts = (y_center[i] >= HD_domain_corners[0][1]) &&
                            // (y_center[i] <= HD_domain_corners[1][1]);
    // for (size_t j = 0; j < num_cells[0]; ++j, ++k) {
      // int within_lon_bounts = (x_center[j] >= HD_domain_corners[0][0]) &&
                              // (x_center[j] <= HD_domain_corners[1][0]);
      // cell_mask[k] = cell_mask[k] && within_lon_bounts && within_lat_bounts;
    // }
  // }

  // define icon grid
  yac_cdef_grid_reg2d(
    "land_grid", num_vertices, cyclic, x_vertices, y_vertices, &grid_id);

  // set global cell ids and core mask
  int * global_cell_id =
    malloc(num_cells[0] * num_cells[1] * sizeof(*global_cell_id));
  for (int i = 0; i < num_cells[0] * num_cells[1]; ++i)
    global_cell_id[i] = i;
  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  // yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);

  int cell_point_id;

  // define field locations (at cell centers)
  yac_cdef_points_reg2d(
    grid_id, num_cells, YAC_LOCATION_CELL, x_center, y_center, &cell_point_id);

  // define field mask
  yac_cset_mask(cell_mask, cell_point_id);

  int const nlev = 1;

  // define fields
  char const * timestep = "24";
  yac_cdef_field(
    "RUNOFF_S", comp_id, &cell_point_id, 1, nlev,
    timestep, YAC_TIME_UNIT_HOUR, &runoff_s_field_id);
  yac_cdef_field(
    "RUNOFF_G", comp_id, &cell_point_id, 1, nlev,
    timestep, YAC_TIME_UNIT_HOUR, &runoff_g_field_id);

  free(cell_mask);
  free(x_vertices);
  free(y_vertices);

  // setup coupling
  yac_cenddef ( );

  // allocate field arrays
  double * runoff_s = malloc(num_cells[0] * num_cells[1] * sizeof(*runoff_s));
  double * runoff_g = malloc(num_cells[0] * num_cells[1] * sizeof(*runoff_g));

  int err;
  int info = 0;

  while ((info != 7) && (file.curr_timestep < file.ntimesteps)) {

    read_timestep(&file, runoff_s, runoff_g);

    { // send surface water runoff
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};

      point_set_data[0] = runoff_s;
      yac_cput(runoff_s_field_id, nlev, collection_data, &info, &err);
    }

    { // send soil water runoff
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};

      point_set_data[0] = runoff_g;
      yac_cput(runoff_g_field_id, nlev, collection_data, &info, &err);
    }
  }

  free(runoff_s);
  free(runoff_g);

  close_input_file(file);

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

static struct forcing_file open_input_file(
  char const * filename, double ** lon, double ** lat, int ** cell_mask) {

  struct forcing_file file;

  file.curr_timestep = 0;

  // open file
  handle_error(nc_open(filename, NC_NOWRITE, &file.ncid));

  // get dimensions
  int lon_dimid, lat_dimid, time_dimid;
  handle_error(nc_inq_dimid(file.ncid, "lon", &lon_dimid));
  handle_error(nc_inq_dimid(file.ncid, "lat", &lat_dimid));
  handle_error(nc_inq_dimid(file.ncid, "time", &time_dimid));

  // get length of dimensions
  handle_error(nc_inq_dimlen(file.ncid, lon_dimid, &(file.nlon)));
  handle_error(nc_inq_dimlen(file.ncid, lat_dimid, &(file.nlat)));
  handle_error(nc_inq_dimlen(file.ncid, time_dimid, &(file.ntimesteps)));

  // get variables
  int lon_varid, lat_varid;
  handle_error(nc_inq_varid(file.ncid, "lon", &lon_varid));
  handle_error(nc_inq_varid(file.ncid, "lat", &lat_varid));
  handle_error(nc_inq_varid(file.ncid, "qs", &file.qs_varid));
  handle_error(nc_inq_varid(file.ncid, "qsb", &file.qsb_varid));

  // check dimensions of variable
  int lon_ndims, lat_ndims, qs_ndims, qsb_ndims;
  handle_error(nc_inq_varndims(file.ncid, lon_varid, &lon_ndims));
  handle_error(nc_inq_varndims(file.ncid, lat_varid, &lat_ndims));
  handle_error(nc_inq_varndims(file.ncid, file.qs_varid, &qs_ndims));
  handle_error(nc_inq_varndims(file.ncid, file.qsb_varid, &qsb_ndims));

  if ((lon_ndims != 1) || (lat_ndims != 1)) {
    fputs(
      "ERROR(open_input_file): invalid number of dimensions (lon lat)", stderr);
    exit(EXIT_FAILURE);
  }
  if ((qs_ndims != 3) || (qsb_ndims != 3)) {
    fputs(
      "ERROR(open_input_file): invalid number of dimensions (qs qsb)", stderr);
    exit(EXIT_FAILURE);
  }

  int lon_var_dimid, lat_var_dimid, qs_var_dimids[3], qsb_var_dimids[3];
  handle_error(nc_inq_vardimid(file.ncid, lon_varid, &lon_var_dimid));
  handle_error(nc_inq_vardimid(file.ncid, lat_varid, &lat_var_dimid));
  handle_error(nc_inq_vardimid(file.ncid, file.qs_varid, qs_var_dimids));
  handle_error(nc_inq_vardimid(file.ncid, file.qsb_varid, qsb_var_dimids));

  if ((lon_var_dimid != lon_dimid) || (lat_var_dimid != lat_dimid)) {
    fputs(
      "ERROR(open_input_file): dimension mismatch(lon lat)", stderr);
    exit(EXIT_FAILURE);
  }
  if ((qs_var_dimids[0] != time_dimid) ||
      (qs_var_dimids[1] != lat_dimid) ||
      (qs_var_dimids[2] != lon_dimid) ||
      (qsb_var_dimids[0] != time_dimid) ||
      (qsb_var_dimids[1] != lat_dimid) ||
      (qsb_var_dimids[2] != lon_dimid)) {
    fputs(
      "ERROR(open_input_file): dimension mismatch(qs qsb)", stderr);
    exit(EXIT_FAILURE);
  }

  // get grid dimensions
  *lon = malloc(file.nlon * sizeof(**lon));
  *lat = malloc(file.nlat * sizeof(**lat));

  handle_error(nc_get_var_double(file.ncid, lon_varid, *lon));
  handle_error(nc_get_var_double(file.ncid, lat_varid, *lat));

#define RAD (0.01745329251994329576923690768489)
  for (size_t i = 0; i < file.nlon; ++i) (*lon)[i] *= RAD;
  for (size_t i = 0; i < file.nlat; ++i) (*lat)[i] *= RAD;

  size_t start[3] = {0, 0, 0};
  size_t count[3] = {1, file.nlat, file.nlon};
  double * qs = malloc(file.nlat * file.nlon * sizeof(*qs));
  handle_error(nc_get_vara_double(file.ncid, file.qs_varid, start, count, qs));

  *cell_mask = malloc(file.nlat * file.nlon * sizeof(**cell_mask));
  for (size_t i = 0; i < file.nlat * file.nlon; ++i)
    (*cell_mask)[i] = qs[i] < 10.0;

  free(qs);

  return file;
}

static void read_timestep(
  struct forcing_file * file, double * qs, double * qsb) {

  if (file->curr_timestep >= file->ntimesteps) {
    fputs(
      "ERROR(read_timestep): no more data available", stderr);
    exit(EXIT_FAILURE);
  }

  size_t start[3] = {file->curr_timestep, 0, 0};
  size_t count[3] = {1, file->nlat, file->nlon};

  handle_error(
    nc_get_vara_double(file->ncid, file->qs_varid, start, count, qs));
  handle_error(
    nc_get_vara_double(file->ncid, file->qsb_varid, start, count, qsb));

  file->curr_timestep++;
}

static void close_input_file(struct forcing_file file) {
  handle_error(nc_close(file.ncid));
}

static void receive_HD_domain(double * HD_domain_corners) {

  char const * comp_names[2] = {"LAND", "HD"};
  MPI_Comm comm;
  yac_cget_comps_comm(comp_names, 2, &comm);

  int root = INT_MAX;
  MPI_Allreduce(MPI_IN_PLACE, &root, 1, MPI_INT, MPI_MIN, comm);
  MPI_Bcast(HD_domain_corners, 4, MPI_DOUBLE, root, comm);

  if (HD_domain_corners[0] > HD_domain_corners[2]) {
    double temp = HD_domain_corners[0];
    HD_domain_corners[0] = HD_domain_corners[2];
    HD_domain_corners[2] = temp;
  }
  if (HD_domain_corners[1] > HD_domain_corners[3]) {
    double temp = HD_domain_corners[1];
    HD_domain_corners[1] = HD_domain_corners[3];
    HD_domain_corners[3] = temp;
  }
}
