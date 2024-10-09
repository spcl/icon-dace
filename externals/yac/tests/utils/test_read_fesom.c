// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <netcdf.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "tests.h"
#include "read_fesom_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

static void write_dummy_grid_file(char * name);
static void check_grid(struct yac_basic_grid_data grid);

int main(void) {

   write_dummy_grid_file("fesom_grid.nc");

   struct yac_basic_grid_data fesom_grid =
     yac_read_fesom_basic_grid_data("fesom_grid.nc");

   unlink("fesom_grid.nc");

   check_grid(fesom_grid);

// #define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
   yac_write_basic_grid_data_to_file(&fesom_grid, "fesom");
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
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "ncells", 8, &dim_ncells_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "nv", 9, &dim_nv_id));

  dim_ids[0] = dim_ncells_id;
  dim_ids[1] = dim_nv_id;

  int var_lon_id, var_lon_v_id, var_lat_id, var_lat_v_id;

  // define grid arrays
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dim_ncells_id, &var_lon_id));
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "lon_vertices", NC_DOUBLE, 2, dim_ids, &var_lon_v_id));
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dim_ncells_id, &var_lat_id));
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "lat_vertices", NC_DOUBLE, 2, dim_ids, &var_lat_v_id));

  // end definition
  YAC_HANDLE_ERROR(nc_enddef(ncid));

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

  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lon_id, lon));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lon_v_id, lon_vertices));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lat_id, lat));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lat_v_id, lat_vertices));

  YAC_HANDLE_ERROR(nc_close(ncid));
}

static void check_grid(struct yac_basic_grid_data grid) {

  size_t ref_num_cells = 8;

  if (grid.num_cells != ref_num_cells)
    PUT_ERR("wrong number of grid cells");

  int ref_num_corners_per_cell[8] = {6,3,3,3,4,6,6,3};

  for (size_t i = 0; i < ref_num_cells; ++i)
    if (grid.num_vertices_per_cell[i] != ref_num_corners_per_cell[i])
      PUT_ERR("wrong number of corners per cell");
}

