// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <netcdf.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "read_scrip_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

// #define HAVE_OASIS_FILES
// #define WRITE_VTK_GRID_FILE
#define DUMMY_GRID_NAME ("dummy_grid")

#ifndef HAVE_OASIS_FILES
static void write_dummy_grid_file(
  char * grid_name, char * grid_filename, char * mask_filename);
#endif

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_size, comm_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if ((comm_size != 4) && (comm_rank == 0)) {
    fputs("wrong number of processes (has to be four)\n", stderr);
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_OASIS_FILES
  char * grid_filename = "../examples/OASIS-grid/grids.nc";
  char * mask_filename = "../examples/OASIS-grid/masks_no_atm.nc";
#else
  char * grid_filename = "test_read_scrip_parallel_grids.nc";
  char * mask_filename = "test_read_scrip_parallel_masks.nc";
  if (comm_rank == 0)
    write_dummy_grid_file(DUMMY_GRID_NAME, grid_filename, mask_filename);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int valid_mask_value = 0;

  if (!yac_file_exists(grid_filename)) return EXIT_FAILURE;
  if (!yac_file_exists(mask_filename)) return EXIT_FAILURE;

  { // read unstructured grid
#ifdef HAVE_OASIS_FILES
    // char * gridname = "bggd";
    // char * gridname = "icoh";
    // char * gridname = "icos";
    // char * gridname = "nogt";
    char * gridname = "sse7";
    // char * gridname = "torc";
    // char * gridname = "ssea";
#else
    char * gridname = DUMMY_GRID_NAME;
#endif

    size_t * duplicated_cell_idx = NULL;
    yac_int * orig_cell_global_ids = NULL;
    size_t nbr_duplicated_cells = 0;

    size_t cell_coord_idx;

    set_even_io_rank_list(MPI_COMM_WORLD);
    struct yac_basic_grid * scrip_grid =
      yac_read_scrip_basic_grid_parallel(
        grid_filename, mask_filename, MPI_COMM_WORLD,
        gridname, valid_mask_value, gridname,
        0, &cell_coord_idx, &duplicated_cell_idx,
        &orig_cell_global_ids, &nbr_duplicated_cells);

#ifdef WRITE_VTK_GRID_FILE
    char vtk_gridname[64];
    snprintf(vtk_gridname, sizeof(vtk_gridname), "%s_%d", gridname, comm_rank);
    yac_write_basic_grid_to_file(scrip_grid, vtk_gridname);
#endif // WRITE_VTK_GRID_FILE

    yac_basic_grid_delete(scrip_grid);
    clear_yac_io_env();
  }

  { // read lon/lat grid
#ifdef HAVE_OASIS_FILES
    // char * gridname = "bggd";
    // char * gridname = "icoh";
    // char * gridname = "icos";
    // char * gridname = "nogt";
    // char * gridname = "sse7";
    char * gridname = "torc";
    // char * gridname = "ssea";
#else
    char * gridname = DUMMY_GRID_NAME;
#endif

    set_even_io_rank_list(MPI_COMM_WORLD);
    struct yac_basic_grid * scrip_grid =
      yac_read_scrip_basic_grid_parallel(
        grid_filename, mask_filename, MPI_COMM_WORLD,
        gridname, valid_mask_value, gridname,
        1, NULL, NULL, NULL, NULL);

#ifdef WRITE_VTK_GRID_FILE
    char vtk_gridname[64];
    snprintf(
      vtk_gridname, sizeof(vtk_gridname), "%s_ll_%d", gridname, comm_rank);
    yac_write_basic_grid_to_file(scrip_grid, vtk_gridname);
#endif // WRITE_VTK_GRID_FILE

    yac_basic_grid_delete(scrip_grid);
    clear_yac_io_env();
  }

#ifndef HAVE_OASIS_FILES
  MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank == 0) {
    unlink(grid_filename);
    unlink(mask_filename);
  }
#endif

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

#ifndef HAVE_OASIS_FILES
static void write_dummy_grid_file(
  char * grid_name, char * grid_filename, char * mask_filename) {

  enum {NUM_LAT = 10, NUM_LON = 361};

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
    YAC_HANDLE_ERROR(nc_def_dim(ncid, crn_dim_name, 4, &dim_crn_id));
    YAC_HANDLE_ERROR(nc_def_dim(ncid, x_dim_name, NUM_LON, &dim_x_id));
    YAC_HANDLE_ERROR(nc_def_dim(ncid, y_dim_name, NUM_LAT, &dim_y_id));

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
    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, cla_var_name, NC_DOUBLE, 3, corner_dim_ids, &var_cla_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_cla_id, "units", strlen(degree), degree));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_cla_id, "title", strlen(title), title));

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, clo_var_name, NC_DOUBLE, 3, corner_dim_ids, &var_clo_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_clo_id, "units", strlen(degree), degree));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_clo_id, "title", strlen(title), title));

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, lat_var_name, NC_DOUBLE, 2, cell_dim_ids, &var_lat_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lat_id, "units", strlen(degree), degree));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lat_id, "title", strlen(title), title));

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, lon_var_name, NC_DOUBLE, 2, cell_dim_ids, &var_lon_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lon_id, "units", strlen(degree), degree));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lon_id, "title", strlen(title), title));


    // end definition
    YAC_HANDLE_ERROR(nc_enddef(ncid));

    // write grid data

    double cla[4][NUM_LAT][NUM_LON];
    double clo[4][NUM_LAT][NUM_LON];
    double lat[NUM_LAT][NUM_LON];
    double lon[NUM_LAT][NUM_LON];

    for (int i = 0; i < NUM_LON; ++i) {
      for (int j = 0; j < NUM_LAT; ++j) {
        cla[0][j][i] = (double)(j+0-90);
        cla[1][j][i] = (double)(j+0-90);
        cla[2][j][i] = (double)(j+1-90);
        cla[3][j][i] = (double)(j+1-90);
        clo[0][j][i] = (double)(i+0);
        clo[1][j][i] = (double)(i+1);
        clo[2][j][i] = (double)(i+1);
        clo[3][j][i] = (double)(i+0);
        lat[j][i] = 0.5 + (double)(j-90);
        lon[j][i] = 0.5 + (double)i;
      }
    }

    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_cla_id, &cla[0][0][0]));
    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_clo_id, &clo[0][0][0]));
    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lat_id, &lat[0][0]));
    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lon_id, &lon[0][0]));

    YAC_HANDLE_ERROR(nc_close(ncid));
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
    YAC_HANDLE_ERROR(nc_def_dim(ncid, x_dim_name, NUM_LON, &dim_x_id));
    YAC_HANDLE_ERROR(nc_def_dim(ncid, y_dim_name, NUM_LAT, &dim_y_id));

    char frc_var_name[128];
    char msk_var_name[128];

    sprintf(frc_var_name, "%s.frc", grid_name);
    sprintf(msk_var_name, "%s.msk", grid_name);

    int dim_ids[2] = {dim_y_id, dim_x_id};

    int var_frc_id;
    int var_msk_id;

    char adim[] = "adim";

    // define variable
    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, frc_var_name, NC_DOUBLE, 2, dim_ids, &var_frc_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_frc_id, "units", strlen(adim), adim));

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, msk_var_name, NC_INT, 2, dim_ids, &var_msk_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_msk_id, "units", strlen(adim), adim));


    // end definition
    YAC_HANDLE_ERROR(nc_enddef(ncid));

    // write grid data

    double frc[NUM_LAT][NUM_LON];
    int msk[NUM_LAT][NUM_LON];

    for (int i = 0; i < NUM_LON; ++i) {
      for (int j = 0; j < NUM_LAT; ++j) {
        frc[j][i] = 1;
        msk[j][i] = 0;
      }
    }

    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_frc_id, &frc[0][0]));
    YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_msk_id, &msk[0][0]));

    YAC_HANDLE_ERROR(nc_close(ncid));
  }
}
#endif // HAVE_OASIS_FILES
