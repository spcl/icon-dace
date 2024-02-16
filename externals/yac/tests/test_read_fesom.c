/**
 * @file test_read_fesom.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
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

#include <netcdf.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "tests.h"
#include "read_fesom_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

static void write_dummy_grid_file(char * name);
static void check_grid(struct basic_grid_data grid);

int main(void) {

   write_dummy_grid_file("fesom_grid.nc");

   struct basic_grid_data fesom_grid = read_fesom_grid("fesom_grid.nc");

   unlink("fesom_grid.nc");

   check_grid(fesom_grid);

#define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
   yac_write_basic_grid_to_file(&fesom_grid, "fesom");
#endif // WRITE_VTK_GRID_FILE

   yac_basic_grid_data_free(fesom_grid);

   return TEST_EXIT_CODE;
}

static void write_dummy_grid_file(char * file_name) {

  int ncid;

  // create file
  yac_nc_create(file_name, NC_CLOBBER, &ncid);

  int dim_ncells_id;
  int dim_nv_id;
  int dim_ids[2];

  // define dimensions
  HANDLE_ERROR(nc_def_dim(ncid, "ncells", 8, &dim_ncells_id));
  HANDLE_ERROR(nc_def_dim(ncid, "nv", 9, &dim_nv_id));

  dim_ids[0] = dim_ncells_id;
  dim_ids[1] = dim_nv_id;

  int var_lon_id, var_lon_v_id, var_lat_id, var_lat_v_id;

  // define grid arrays
  HANDLE_ERROR(
    nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dim_ncells_id, &var_lon_id));
  HANDLE_ERROR(
    nc_def_var(ncid, "lon_vertices", NC_DOUBLE, 2, dim_ids, &var_lon_v_id));
  HANDLE_ERROR(
    nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dim_ncells_id, &var_lat_id));
  HANDLE_ERROR(
    nc_def_var(ncid, "lat_vertices", NC_DOUBLE, 2, dim_ids, &var_lat_v_id));

  // end definition
  HANDLE_ERROR(nc_enddef(ncid));

  // write grid data

  double lon[8] = {1,2,3,4,3,3,5.5,5.5};
  double lon_vertices[8*9] = {0,1,2,2,1,0,0,0,0,
                              1,3,2,2,2,2,2,2,2,
                              3,4,2,2,2,2,2,2,2,
                              3,5,4,4,4,4,4,4,4,
                              2,4,4,2,2,2,2,2,2,
                              1,2,4,5,4,2,2,2,2,
                              4,6,7,6,5,4,4,4,4,
                              5,6,4,4,4,4,4,4,4};
  double lat[8] = {1.5,0.5,0.5,0.5,1.5,3,1.5,0.5};
  double lat_vertices[8*9] = {1,0,1,2,3,2,2,2,2,
                              0,0,1,1,1,1,1,1,1,
                              0,1,1,1,1,1,1,1,1,
                              0,0,1,1,1,1,1,1,1,
                              1,1,2,2,2,2,2,2,2,
                              3,2,2,3,4,4,4,4,4,
                              1,1,2,3,3,2,2,2,2,
                              0,1,1,1,1,1,1,1,1};

  HANDLE_ERROR(nc_put_var_double(ncid, var_lon_id, lon));
  HANDLE_ERROR(nc_put_var_double(ncid, var_lon_v_id, lon_vertices));
  HANDLE_ERROR(nc_put_var_double(ncid, var_lat_id, lat));
  HANDLE_ERROR(nc_put_var_double(ncid, var_lat_v_id, lat_vertices));

  HANDLE_ERROR(nc_close(ncid));
}

static void check_grid(struct basic_grid_data grid) {

  size_t ref_num_cells = 8;

  if (grid.num_cells != ref_num_cells)
    PUT_ERR("wrong number of grid cells");

  int ref_num_corners_per_cell[8] = {6,3,3,3,4,6,6,3};

  for (size_t i = 0; i < ref_num_cells; ++i)
    if (grid.num_vertices_per_cell[i] != ref_num_corners_per_cell[i])
      PUT_ERR("wrong number of corners per CELL");
}

