// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <netcdf.h>

#include "tests.h"
#include "test_common.h"
#include "read_woa_data.h"
#include "io_utils.h"

#define LON (16)
#define LAT (8)
#define DEPTH (4)
#define TIME (2)
#define NV (2)

static void write_test_data_file(char const * file_name);

int main(void) {

  char const * test_data_file_name = "test_read_woa_data.nc";

  write_test_data_file(test_data_file_name);

  {
    // open woa data file
    char * woa_field_name = "s_an";
    int woa_file = yac_open_woa_output(test_data_file_name);

    struct yac_fieldMetadata field_info;
    yac_read_woa_dimensions(woa_file, woa_field_name, &field_info);

    // check meta data
    if (field_info.nbrTimeSteps != TIME)
      PUT_ERR("ERROR in meta data (nbrTimeSteps)");
    if (field_info.nbrLevels != DEPTH)
      PUT_ERR("ERROR in meta data (nbrLevels)");
    if (field_info.nbrLatPoints != LAT)
      PUT_ERR("ERROR in meta data (nbrLatPoints)");
    if (field_info.nbrLonPoints != LON)
      PUT_ERR("ERROR in meta data (nbrLonPoints)");

    double * global_salinity = yac_get_woa_memory(field_info);

    for (int time = 0; time < TIME; ++time) {
      for (int depth = 0; depth < DEPTH; ++depth) {

        yac_read_woa_timestep_level(
          woa_file, global_salinity, field_info, time+1, depth+1);

        // check retrieved field data
        for (int i = 0; i < LON * LAT; ++i)
          if (fabs(global_salinity[i] - (double)(time * 100 + depth)) > 1e-6)
            PUT_ERR("ERROR in field data");
      }
    }

    yac_free_woa_memory(global_salinity);

    yac_close_woa_output(woa_file);
  }

  unlink(test_data_file_name);

  return TEST_EXIT_CODE;
}

static void write_test_data_file(char const * file_name) {

  int ncid;

  // create file
  yac_nc_create(file_name, NC_CLOBBER, &ncid);

  int dim_lon_id, dim_lat_id, dim_depth_id, dim_time_id, dim_nv_id;

  // define dimensions
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "lon", LON, &dim_lon_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "lat", LAT, &dim_lat_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "depth", DEPTH, &dim_depth_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "time", TIME, &dim_time_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "nv", NV, &dim_nv_id));

  int var_s_an_id;

  // define variables
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "s_an", NC_FLOAT, 4,
      (int[]){dim_time_id,dim_depth_id,dim_lat_id,dim_lon_id},
      &var_s_an_id));

  // end definition
  YAC_HANDLE_ERROR(nc_enddef(ncid));

  // write dummy data
  float dummy_data[TIME][DEPTH][LAT][LON];
  for (int time = 0; time < TIME; ++time)
    for (int depth = 0; depth < DEPTH; ++depth)
      for (int lat = 0; lat < LAT; ++lat)
        for (int lon = 0; lon < LON; ++lon)
          dummy_data[time][depth][lat][lon] = (float)(time * 100 + depth);
  YAC_HANDLE_ERROR(nc_put_var_float(ncid, var_s_an_id, &dummy_data[0][0][0][0]));

  // close file
  YAC_HANDLE_ERROR(nc_close(ncid));
}
