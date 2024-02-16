/**
 * @file test_read_scrip.c
 *
 * @copyright Copyright  (C)  2020 DKRZ, MPI-M
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
#include <string.h>
#include <unistd.h>

#include "tests.h"
#include "test_common.h"
#include "read_scrip_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

static void write_dummy_grid_file(
  char * grid_name, char * grid_filename, char * mask_filename);

int main(void) {

// #define HAVE_OASIS_FILES
#ifdef HAVE_OASIS_FILES
  char * grid_filename = "../examples/OASIS-grid/grids.nc";
  char * mask_filename = "../examples/OASIS-grid/masks_no_atm.nc";
  // char * gridname = "bggd";
  // char * gridname = "icoh";
  // char * gridname = "icos";
  // char * gridname = "nogt";
  char * gridname = "sse7";
  // char * gridname = "torc";
  // char * gridname = "ssea";
#else
  char * grid_filename = "test_read_scrip_grids.nc";
  char * mask_filename = "test_read_scrip_masks.nc";
  char * gridname = "dummy_grid";
  write_dummy_grid_file(gridname, grid_filename, mask_filename);
#endif

  int valid_mask_value = 0;

  if (!yac_file_exists(grid_filename)) return EXIT_FAILURE;
  if (!yac_file_exists(mask_filename)) return EXIT_FAILURE;

  struct basic_grid_data scrip_grid =
    read_scrip_grid(grid_filename, mask_filename, gridname, valid_mask_value);

// #define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
   yac_write_basic_grid_to_file(&scrip_grid, gridname);
#endif // WRITE_VTK_GRID_FILE

  yac_basic_grid_data_free(scrip_grid);

#ifndef HAVE_OASIS_FILES
  unlink(grid_filename);
  unlink(mask_filename);
#endif

  return TEST_EXIT_CODE;
}

static void write_dummy_grid_file(
  char * grid_name, char * grid_filename, char * mask_filename) {

  { // grid file
    int ncid;

    // create file
    yac_nc_create(grid_filename, NC_CLOBBER, &ncid);

    char crn_dim_name[128];
    char x_dim_name[128];
    char y_dim_name[128];

    sprintf(crn_dim_name, "crn_%s", grid_name);
    sprintf(x_dim_name, "x_%s", grid_name);
    sprintf(y_dim_name, "y_%s", grid_name);

    int dim_crn_id;
    int dim_x_id;
    int dim_y_id;

    // define dimensions
    HANDLE_ERROR(nc_def_dim(ncid, crn_dim_name, 4, &dim_crn_id));
    HANDLE_ERROR(nc_def_dim(ncid, x_dim_name, 20, &dim_x_id));
    HANDLE_ERROR(nc_def_dim(ncid, y_dim_name, 10, &dim_y_id));

    char cla_var_name[128];
    char clo_var_name[128];
    char lat_var_name[128];
    char lon_var_name[128];

    sprintf(cla_var_name, "%s.cla", grid_name);
    sprintf(clo_var_name, "%s.clo", grid_name);
    sprintf(lat_var_name, "%s.lat", grid_name);
    sprintf(lon_var_name, "%s.lon", grid_name);

    int corner_dim_ids[3] = {dim_crn_id, dim_y_id, dim_x_id};
    int cell_dim_ids[2] = {dim_y_id, dim_x_id};

    int var_cla_id;
    int var_clo_id;
    int var_lat_id;
    int var_lon_id;

    char degree[] = "degree";
    char title[] = "This is a reg lon-lat dummy grid";

    // define variable
    HANDLE_ERROR(
      nc_def_var(
        ncid, cla_var_name, NC_DOUBLE, 3, corner_dim_ids, &var_cla_id));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_cla_id, "units", strlen(degree), degree));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_cla_id, "title", strlen(title), title));

    HANDLE_ERROR(
      nc_def_var(
        ncid, clo_var_name, NC_DOUBLE, 3, corner_dim_ids, &var_clo_id));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_clo_id, "units", strlen(degree), degree));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_clo_id, "title", strlen(title), title));

    HANDLE_ERROR(
      nc_def_var(
        ncid, lat_var_name, NC_DOUBLE, 2, cell_dim_ids, &var_lat_id));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_lat_id, "units", strlen(degree), degree));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_lat_id, "title", strlen(title), title));

    HANDLE_ERROR(
      nc_def_var(
        ncid, lon_var_name, NC_DOUBLE, 2, cell_dim_ids, &var_lon_id));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_lon_id, "units", strlen(degree), degree));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_lon_id, "title", strlen(title), title));


    // end definition
    HANDLE_ERROR(nc_enddef(ncid));

    // write grid data

    double cla[4][10][20];
    double clo[4][10][20];
    double lat[10][20];
    double lon[10][20];

    for (int i = 0; i < 20; ++i) {
      for (int j = 0; j < 10; ++j) {
        cla[0][j][i] = (double)(j+0);
        cla[1][j][i] = (double)(j+0);
        cla[2][j][i] = (double)(j+1);
        cla[3][j][i] = (double)(j+1);
        clo[0][j][i] = (double)(i+0);
        clo[1][j][i] = (double)(i+1);
        clo[2][j][i] = (double)(i+1);
        clo[3][j][i] = (double)(i+0);
        lat[j][i] = 0.5 + (double)j;
        lon[j][i] = 0.5 + (double)i;
      }
    }

    HANDLE_ERROR(nc_put_var_double(ncid, var_cla_id, &cla[0][0][0]));
    HANDLE_ERROR(nc_put_var_double(ncid, var_clo_id, &clo[0][0][0]));
    HANDLE_ERROR(nc_put_var_double(ncid, var_lat_id, &lat[0][0]));
    HANDLE_ERROR(nc_put_var_double(ncid, var_lon_id, &lon[0][0]));

    HANDLE_ERROR(nc_close(ncid));
  }

  { // mask file
    int ncid;

    // create file
    yac_nc_create(mask_filename, NC_CLOBBER, &ncid);

    char x_dim_name[128];
    char y_dim_name[128];

    sprintf(x_dim_name, "x_%s", grid_name);
    sprintf(y_dim_name, "y_%s", grid_name);

    int dim_x_id;
    int dim_y_id;

    // define dimensions
    HANDLE_ERROR(nc_def_dim(ncid, x_dim_name, 20, &dim_x_id));
    HANDLE_ERROR(nc_def_dim(ncid, y_dim_name, 10, &dim_y_id));

    char frc_var_name[128];
    char msk_var_name[128];

    sprintf(frc_var_name, "%s.frc", grid_name);
    sprintf(msk_var_name, "%s.msk", grid_name);

    int dim_ids[2] = {dim_y_id, dim_x_id};

    int var_frc_id;
    int var_msk_id;

    char adim[] = "adim";

    // define variable
    HANDLE_ERROR(
      nc_def_var(
        ncid, frc_var_name, NC_DOUBLE, 2, dim_ids, &var_frc_id));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_frc_id, "units", strlen(adim), adim));

    HANDLE_ERROR(
      nc_def_var(
        ncid, msk_var_name, NC_INT, 2, dim_ids, &var_msk_id));
    HANDLE_ERROR(
      nc_put_att_text(ncid, var_msk_id, "units", strlen(adim), adim));


    // end definition
    HANDLE_ERROR(nc_enddef(ncid));

    // write grid data

    double frc[10][20];
    int msk[10][20];

    for (int i = 0; i < 20; ++i) {
      for (int j = 0; j < 10; ++j) {
        frc[j][i] = 1;
        msk[j][i] = 0;
      }
    }

    HANDLE_ERROR(nc_put_var_double(ncid, var_frc_id, &frc[0][0]));
    HANDLE_ERROR(nc_put_var_int(ncid, var_msk_id, &msk[0][0]));

    HANDLE_ERROR(nc_close(ncid));
  }
}
