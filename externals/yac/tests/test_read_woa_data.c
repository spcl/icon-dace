/**
 * @file test_read_woa_data.c
 *
 * @copyright Copyright  (C)  2022 DKRZ, MPI-M
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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <netcdf.h>

#include "tests.h"
#include "test_common.h"
#include "yac_interface.h"
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
    int woa_file = open_woa_output(test_data_file_name);

    struct fieldMetadata field_info;
    read_woa_dimensions(woa_file, woa_field_name, &field_info);

    // check meta data
    if (field_info.nbrTimeSteps != TIME)
      PUT_ERR("ERROR in meta data (nbrTimeSteps)");
    if (field_info.nbrLevels != DEPTH)
      PUT_ERR("ERROR in meta data (nbrLevels)");
    if (field_info.nbrLatPoints != LAT)
      PUT_ERR("ERROR in meta data (nbrLatPoints)");
    if (field_info.nbrLonPoints != LON)
      PUT_ERR("ERROR in meta data (nbrLonPoints)");

    double * global_salinity = get_woa_memory(field_info);

    for (int time = 0; time < TIME; ++time) {
      for (int depth = 0; depth < DEPTH; ++depth) {

        read_woa_timestep_level(
          woa_file, global_salinity, field_info, time+1, depth+1);

        // check retrieved field data
        for (int i = 0; i < LON * LAT; ++i)
          if (fabs(global_salinity[i] - (double)(time * 100 + depth)) > 1e-6)
            PUT_ERR("ERROR in field data");
      }
    }

    free_woa_memory(global_salinity);

    close_woa_output(woa_file);
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
  HANDLE_ERROR(nc_def_dim(ncid, "lon", LON, &dim_lon_id));
  HANDLE_ERROR(nc_def_dim(ncid, "lat", LAT, &dim_lat_id));
  HANDLE_ERROR(nc_def_dim(ncid, "depth", DEPTH, &dim_depth_id));
  HANDLE_ERROR(nc_def_dim(ncid, "time", TIME, &dim_time_id));
  HANDLE_ERROR(nc_def_dim(ncid, "nv", NV, &dim_nv_id));

  int var_s_an_id;

  // define variables
  HANDLE_ERROR(
    nc_def_var(
      ncid, "s_an", NC_FLOAT, 4,
      (int[]){dim_time_id,dim_depth_id,dim_lat_id,dim_lon_id},
      &var_s_an_id));

  // end definition
  HANDLE_ERROR(nc_enddef(ncid));

  // write dummy data
  float dummy_data[TIME][DEPTH][LAT][LON];
  for (int time = 0; time < TIME; ++time)
    for (int depth = 0; depth < DEPTH; ++depth)
      for (int lat = 0; lat < LAT; ++lat)
        for (int lon = 0; lon < LON; ++lon)
          dummy_data[time][depth][lat][lon] = (float)(time * 100 + depth);
  HANDLE_ERROR(nc_put_var_float(ncid, var_s_an_id, &dummy_data[0][0][0][0]));

  // close file
  HANDLE_ERROR(nc_close(ncid));
}
