// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#include "yac_core.h"

static void get_grid_dimensions(
  char const * grids_file_name,
  char const * grid_names[], size_t count, size_t (*grid_sizes)[2]);
static int create_mask_file(
  char const * ocean_grid_name, char const **atmos_grid_names,
  size_t atmo_count, size_t (*atmo_grid_sizes)[2], double threshold,
  int (*var_ids)[2]);
static double * read_interp_field(
  char const * path, char const * ocean_grid_name,
  char const * atmo_grid_name, size_t atmo_size, double const threshold);
static void write_mask(
  int ncid, int var_ids[2], double * frac, int * mask);

int main(void) {

  char const vtk_path[] = "./output";
  char const grids_file_name[] = "./input/grids.nc";

  char const * ocean_names[] = {
    "nogt", "torc"
  };
  char const * atmo_names[] = {
    "bggd", "icos", "sse7", "icoh"
  };
  enum {
    OCEAN_COUNT = sizeof(ocean_names) / sizeof(ocean_names[0]),
    ATMO_COUNT = sizeof(atmo_names) / sizeof(atmo_names[0])
  };

  double const threshold = 0.001;

  // get dimensions from atmosphere grid file
  size_t atmo_sizes_2d[ATMO_COUNT][2];
  get_grid_dimensions(
    grids_file_name, atmo_names, ATMO_COUNT, atmo_sizes_2d);

  for (size_t i = 0; i < OCEAN_COUNT; ++i) {

    // create mask file
    int var_ids[ATMO_COUNT][2];
    int mask_file_ncid =
      create_mask_file(
        ocean_names[i], atmo_names, ATMO_COUNT, atmo_sizes_2d, threshold,
        var_ids);

    for (size_t j = 0; j < ATMO_COUNT; ++j) {

      size_t atmo_size = atmo_sizes_2d[j][0] * atmo_sizes_2d[j][1];

      // read interpolated field
      double * frac =
        read_interp_field(vtk_path, ocean_names[i], atmo_names[j], atmo_size, threshold);

      // generate integer mask
      int * mask = malloc(atmo_size * sizeof(*mask));
      for (size_t k = 0; k < atmo_size; ++k)
        mask[k] = frac[k] > 0.0 ? 0 : 1;

      // write masks
      write_mask(mask_file_ncid, var_ids[j], frac, mask);

      free(mask);
      free(frac);
    }

    // close mask file
    nc_close(mask_file_ncid);
  }

  return EXIT_SUCCESS;
}

static void get_grid_dimensions(
  char const * grids_file_name,
  char const * grid_names[], size_t count, size_t (*grid_sizes)[2]) {

  int ncid;
  yac_nc_open(grids_file_name, NC_NOWRITE, &ncid);

  for (size_t i = 0; i < count; ++i) {

    char const * grid_name = grid_names[i];
    size_t * grid_size = grid_sizes[i];
    size_t grid_name_len = strlen(grid_name) + 1;
    char x_dim_name[2 + grid_name_len];
    char y_dim_name[2 + grid_name_len];

    snprintf(x_dim_name, 2 + grid_name_len, "x_%s", grid_name);
    snprintf(y_dim_name, 2 + grid_name_len, "y_%s", grid_name);

    // get dimension ids
    int x_dim_id;
    int y_dim_id;
    yac_nc_inq_dimid(ncid, x_dim_name, &x_dim_id);
    yac_nc_inq_dimid(ncid, y_dim_name, &y_dim_id);

    // get dimension length
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, grid_size + 0));
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, grid_size + 1));
  }

  YAC_HANDLE_ERROR(nc_close(ncid));
}

static int create_mask_file(
  char const * ocean_grid_name, char const **atmos_grid_names,
  size_t atmo_count, size_t (*atmo_grid_sizes)[2], double threshold,
  int (*var_ids)[2]) {

  char mask_file_name[128];
  snprintf(
    mask_file_name, sizeof(mask_file_name), "./output/masks_%s_yac.nc",
    ocean_grid_name);

  int ncid;
  int dim_ids[atmo_count][2];

  YAC_HANDLE_ERROR(nc_create(mask_file_name, NC_CLOBBER, &ncid));

  // write dimensions
  for (size_t i = 0; i < atmo_count; ++i) {
    char const * grid_name = atmos_grid_names[i];
    size_t grid_name_len = strlen(grid_name) + 1;
    char y_dim_name[2 + grid_name_len];
    char x_dim_name[2 + grid_name_len];

    snprintf(y_dim_name, 2 + grid_name_len, "y_%s", grid_name);
    snprintf(x_dim_name, 2 + grid_name_len, "x_%s", grid_name);

    YAC_HANDLE_ERROR(
      nc_def_dim(ncid, y_dim_name, atmo_grid_sizes[i][1], &(dim_ids[i][0])));
    YAC_HANDLE_ERROR(
      nc_def_dim(ncid, x_dim_name, atmo_grid_sizes[i][0], &(dim_ids[i][1])));
  }

  // define variables
  for (size_t i = 0; i < atmo_count; ++i) {
    char const * grid_name = atmos_grid_names[i];
    size_t grid_name_len = strlen(grid_name) + 1;
    char frac_var_name[5 + grid_name_len];
    char mask_var_name[5 + grid_name_len];

    snprintf(frac_var_name, 5 + grid_name_len, "%s.frc", grid_name);
    snprintf(mask_var_name, 5 + grid_name_len, "%s.msk", grid_name);

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, frac_var_name, NC_DOUBLE, 2, dim_ids[i], &(var_ids[i][0])));
    YAC_HANDLE_ERROR(
      nc_put_att_text(
        ncid, var_ids[i][0], "coherent_with_grid",
        strlen(ocean_grid_name), ocean_grid_name));
    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, mask_var_name, NC_INT, 2, dim_ids[i], &(var_ids[i][1])));
    YAC_HANDLE_ERROR(
      nc_put_att_text(
        ncid, var_ids[i][1], "coherent_with_grid",
        strlen(ocean_grid_name), ocean_grid_name));
    YAC_HANDLE_ERROR(
      nc_put_att_double(
        ncid, var_ids[i][1], "threshold", NC_DOUBLE, 1, &threshold));
  }

  // write global attributes
  char str_description[] = "Created by YAC";
  char str_remapper[] = "YAC";
  YAC_HANDLE_ERROR(
    nc_put_att_text(
      ncid, NC_GLOBAL, "description",
      strlen(str_description), str_description));
  YAC_HANDLE_ERROR(
    nc_put_att_text(
      ncid, NC_GLOBAL, "remapper", strlen(str_remapper), str_remapper));

  YAC_HANDLE_ERROR(nc_enddef(ncid));

  return ncid;
}

static double * read_interp_field(
  char const * path, char const * ocean_grid_name,
  char const * atmo_grid_name, size_t atmo_size, double const threshold) {

  char vtk_file_name[1024];
  snprintf(vtk_file_name, 1024, "%s/%s.vtk", path, atmo_grid_name);

  char field_name[1024];
  snprintf(field_name, 1024, "test_one_%s_CONSERV_DESTAREA", ocean_grid_name);

  double * frac = malloc(atmo_size * sizeof(*frac));

  FILE * vtk_file = fopen(vtk_file_name, "r");

  if (vtk_file == NULL) {
    for (size_t i = 0; i < atmo_size; ++i) frac[i] = 1.0;
    return frac;
  }
  char line_buffer[1024];
  char * line;

  // search for cell declaration
  while(((line = fgets(line_buffer, 1024, vtk_file))) != NULL)
    if (strstr(line, "CELLS") != NULL) break;

  if (line == NULL) {
    fprintf(stderr, "ERROR could not find CELLS declaration in file %s\n",
            vtk_file_name);
    exit(EXIT_FAILURE);
  }

  unsigned num_cells, dummy;
  if (sscanf(line, "CELLS %u %u\n", &num_cells, &dummy) != 2) {
    fprintf(stderr, "ERROR could not read grid size from file %s\n",
            vtk_file_name);
    exit(EXIT_FAILURE);
  }

  if ((size_t)num_cells != atmo_size) {
    fprintf(stderr, "ERROR grid dimension missmatch in file %s\n"
            "value in file:  %u\n"
            "expected value: %zu\n", vtk_file_name, num_cells, atmo_size);
    exit(EXIT_FAILURE);
  }

  // search for field declaration
  while(((line = fgets(line_buffer, 1024, vtk_file))) != NULL)
    if (strstr(line, field_name) != NULL) break;
  if (line == NULL) {
    fprintf(stderr, "ERROR could not find field \"%s\" in file %s\n",
            field_name, vtk_file_name);
    exit(EXIT_FAILURE);
  }
  // skip one line
  line = fgets(line_buffer, 1024, vtk_file);

  for (size_t i = 0; i < atmo_size; ++i) {
    if (fscanf(vtk_file, "%lf\n", frac + i) != 1) {
      fprintf(stderr, "ERROR while reading field %s from file %s\n",
              field_name, vtk_file_name);
      exit(EXIT_FAILURE);
    }
    if (frac[i] < threshold) frac[i] = 0.0;
    if (frac[i] > 1.0) frac[i] = 1.0;
  }

  fclose(vtk_file);

  return frac;
}


static void write_mask(
  int ncid, int var_ids[2], double * frac, int * mask) {

  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_ids[0], frac));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_ids[1], mask));
}
