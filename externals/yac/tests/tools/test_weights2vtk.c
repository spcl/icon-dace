// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <unistd.h> 
#include <stdlib.h>
#include <libgen.h>
#include <string.h>
#include <netcdf.h>

#include "tests.h"
#include "yac_core.h"
#include "weight_file_common.h"

static void write_scrip_grid(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name,
  double * clo, double * cla, double * lon, double * lat,
  size_t crn, size_t x, size_t y, int create_file);

int main (int argc, char* argv[]) {

  if (argc != 4) {
    PUT_ERR("ERROR: wrong number of arguments");
    return TEST_EXIT_CODE;
  }

  char const * weights2vtk_exec = argv[1];

  char const * src_grid_name = "test_weights2vtk_src_grid_name";
  char const * tgt_grid_name = "test_weights2vtk_tgt_grid_name";
  char const * weight_file_name = "test_weights2vtk_weight_file_name.nc";
  char const * vtk_file_name = "test_weights2vtk.vtk";

  // test with cubed grid
  {
    // write dummy weight file
    {
      // these links form a message when viewed in 3D
      int src_indices[] = {316,344,346,384,422,413,448,448,469,507,418,83,214};
      int tgt_indices[] = {93,115,70,85,115,108,108,139,169,167,117,50,160};
      double weights[] = {0,1,1,1,1,2,2,2,2,2,3,0,3};
      size_t num_links = 13;
      enum yac_location src_locations[2] = {YAC_LOC_CORNER, YAC_LOC_CELL};
      enum yac_location tgt_location = YAC_LOC_CELL;
      unsigned num_src_fields = 2;
      int num_links_per_field[2] = {11,2};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        weight_file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    {
      char cmd[1024];
      sprintf(
        cmd, "%s -S C -T C -s 20 -t 15 -w %s -o %s",
        weights2vtk_exec, weight_file_name, vtk_file_name);
      if (system(cmd)) PUT_ERR("failed to execute");
      unlink(weight_file_name);
      unlink(vtk_file_name);
    }
  }

  // test with curved and unstructed grid, but without links in the weight file
  {
    // write dummy weight file
    {
      int * src_indices = NULL;
      int * tgt_indices = NULL;
      double * weights = NULL;
      size_t num_links = 0;
      enum yac_location src_locations[1] = {YAC_LOC_CELL};
      enum yac_location tgt_location = YAC_LOC_CELL;
      unsigned num_src_fields = 1;
      int num_links_per_field[1] = {0};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        weight_file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    {
      char const * format =
        "%s -S m -T i -s %sGR30_lsm.nc -t %sicon_grid_0030_R02B03_G.nc -w %s -o %s";
      char * cmd =
        malloc(strlen(format) + strlen(weights2vtk_exec) + 2 * strlen(argv[2]) +
               strlen(weight_file_name) + strlen(vtk_file_name));
      sprintf(
        cmd, format,
        weights2vtk_exec, argv[2], argv[3], weight_file_name, vtk_file_name);
      if (system(cmd)) PUT_ERR("failed to execute");
      unlink(weight_file_name);
      unlink(vtk_file_name);
      free(cmd);
    }
  }

  // test with cubed grid and field data on the edges
  {
    // write dummy weight file (just some dummy data)
    {
      int src_indices[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
      int tgt_indices[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      double weights[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
      size_t num_links = 20;
      enum yac_location src_locations[1] = {YAC_LOC_EDGE};
      enum yac_location tgt_location = YAC_LOC_EDGE;
      unsigned num_src_fields = 1;
      int num_links_per_field[1] = {20};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        weight_file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    {
      char cmd[1024];
      sprintf(
        cmd, "%s -S c -T c -s 5 -t 5 -w %s -o %s",
        weights2vtk_exec, weight_file_name, vtk_file_name);
      if (system(cmd)) PUT_ERR("failed to execute");
      unlink(weight_file_name);
      unlink(vtk_file_name);
    }
  }

  // test with gaussian grids
  {
    // write dummy weight file (just some dummy data)
    {
      int src_indices[] = {1,2,2,3,3, 1,2,2,3,3, 5,6,6,7,7, 5,6,6,7,7, 9,10,10,11,11,};
      int tgt_indices[] = {4,4,5,5,6, 8,8,9,9,10, 8,8,9,9,10, 12,12,13,13,14, 12,12,13,13,14};
      double weights[] = {0.5,0.25,0.25,0.25,0.25, 0.5,0.25,0.25,0.25,0.25,
                          0.5,0.25,0.25,0.25,0.25, 0.5,0.25,0.25,0.25,0.25,
                          1.0,0.5,0.5,0.5,0.5};
      size_t num_links = 25;
      enum yac_location src_locations[1] = {YAC_LOC_CELL};
      enum yac_location tgt_location = YAC_LOC_CELL;
      unsigned num_src_fields = 1;
      int num_links_per_field[1] = {25};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        weight_file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    {
      char cmd[1024];
      sprintf(
        cmd, "%s -S g -T g -s -2.75,-1.25,1.25,2.75,4,4 -t -1.25,-2.75,2.75,1.25,4,4 -w %s -o %s",
        weights2vtk_exec, weight_file_name, vtk_file_name);
      if (system(cmd)) PUT_ERR("failed to execute");
      unlink(weight_file_name);
      unlink(vtk_file_name);
    }
  }

  // test with scrip formated grid files
  {
    char const * grid_filename = "test_weights2vtk_grids.nc";
    char const * mask_filename = "test_weights2vtk_masks.nc";
    // write dummy source grid
    {
      enum {
        CRN = 4,
        Y = 4,
        X = 5,
      };
      double clo[CRN][Y][X] =
        {{{0,1,2,3,4},{0,1,2,3,4},{0,1,2,3,4},{0,1,2,3,4}},
         {{1,2,3,4,5},{1,2,3,4,5},{1,2,3,4,5},{1,2,3,4,5}},
         {{1,2,3,4,5},{1,2,3,4,5},{1,2,3,4,5},{1,2,3,4,5}},
         {{0,1,2,3,4},{0,1,2,3,4},{0,1,2,3,4},{0,1,2,3,4}}};
      double cla[CRN][Y][X] =
        {{{0,0,0,0,0},{1,1,1,1,1},{2,2,2,2,2},{3,3,3,3,3}},
         {{0,0,0,0,0},{1,1,1,1,1},{2,2,2,2,2},{3,3,3,3,3}},
         {{1,1,1,1,1},{2,2,2,2,2},{3,3,3,3,3},{4,4,4,4,4}},
         {{1,1,1,1,1},{2,2,2,2,2},{3,3,3,3,3},{4,4,4,4,4}}};
      double lon[Y][X] =
        {{0.5,1.5,2.5,3.5,4.5},
         {0.5,1.5,2.5,3.5,4.5},
         {0.5,1.5,2.5,3.5,4.5},
         {0.5,1.5,2.5,3.5,4.5}};
      double lat[Y][X] =
        {{0.5,0.5,0.5,0.5,0.5},
         {1.5,1.5,1.5,1.5,1.5},
         {2.5,2.5,2.5,2.5,2.5},
         {3.5,3.5,3.5,3.5,3.5}};
      write_scrip_grid(
        grid_filename, mask_filename, src_grid_name,
        &clo[0][0][0], &cla[0][0][0], &lon[0][0], &lat[0][0], CRN, X, Y, 1);
    }
    // write dummy target grid
    {
      enum {
        CRN = 3,
        Y = 1,
        X = 8,
      };
      double clo[CRN][Y][X] =
        {{{0  ,0  ,2.5,5  ,5  ,5  ,2.5,0  }},
         {{2.5,2.5,5  ,5  ,5  ,2.5,0  ,0  }},
         {{0  ,2.5,2.5,2.5,2.5,2.5,2.5,2.5}}};
      double cla[CRN][Y][X] =
        {{{0,0,0,0,2,4,4,4}},
         {{2,0,0,2,4,4,4,2}},
         {{2,2,2,2,2,2,2,2}}};
      double lon[Y][X] =
        {{0.5,2,3,4.5,4.5,3,2,0.5}};
      double lat[Y][X] =
        {{1.5,0.5,0.5,1.5,2.5,3.5,3.5,2.5}};
      write_scrip_grid(
        grid_filename, mask_filename, tgt_grid_name,
        &clo[0][0][0], &cla[0][0][0], &lon[0][0], &lat[0][0], CRN, X, Y, 0);
    }
    // write dummy weight file (just some dummy data)
    {
      int src_indices[] =
        {0,5,6,7, 0,1,2,6,7, 2,3,4,7,8, 4,7,8,9,
         12,13,14,19, 12,13,17,18,19, 11,12,15,16,17, 10,11,12,15};
      int tgt_indices[] =
        {0,0,0,0, 1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,
         4,4,4,4, 5,5,5,5,5, 6,6,6,6,6, 7,7,7,7};
      double weights[] =
        {0.5,1.0,0.5,0.125, 0.5,1.0,0.5,0.5,0.375,
         0.5,1.0,0.5,0.375,0.5, 0.5,0.125,0.5,1.0,
         0.125,0.5,1.0,0.5, 0.375,0.5,0.5,1.0,0.5,
         0.5,0.375,0.5,1.0,0.5, 1.0,0.5,0.125,0.5};
      size_t num_links = sizeof(src_indices)/sizeof(src_indices[0]);
      enum yac_location src_locations[1] = {YAC_LOC_CELL};
      enum yac_location tgt_location = YAC_LOC_CELL;
      unsigned num_src_fields = 1;
      int num_links_per_field[1] = {sizeof(src_indices)/sizeof(src_indices[0])};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        weight_file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    {
      char cmd[1024];
      sprintf(
        cmd, "%s -S s -T s -s %s,%s,%s -t %s,%s,%s -w %s -o %s",
        weights2vtk_exec,
        grid_filename, mask_filename, src_grid_name,
        grid_filename, mask_filename, tgt_grid_name,
        weight_file_name, vtk_file_name);
      if (system(cmd)) PUT_ERR("failed to execute");
      unlink(grid_filename);
      unlink(mask_filename);
      unlink(weight_file_name);
      unlink(vtk_file_name);
    }
  }

  return TEST_EXIT_CODE;
}

static void write_scrip_grid(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name,
  double * clo, double * cla, double * lon, double * lat,
  size_t crn, size_t x, size_t y, int create_file) {

  char buffer[128];

  int ncid;

  if (create_file) {
    yac_nc_create(grid_filename, NC_CLOBBER, &ncid);
  } else {
    yac_nc_open(grid_filename, NC_WRITE, &ncid);
    YAC_HANDLE_ERROR(nc_redef(ncid));
  }

  int dim_ids[3], cla_var_id, clo_var_id, lat_var_id, lon_var_id, msk_var_id;
  YAC_HANDLE_ERROR(
    nc_def_dim(
      ncid, strcat(strcpy(buffer, "crn_"), grid_name), crn, &dim_ids[0]));
  YAC_HANDLE_ERROR(
    nc_def_dim(
      ncid, strcat(strcpy(buffer, "y_"), grid_name), y, &dim_ids[1]));
  YAC_HANDLE_ERROR(
    nc_def_dim(
      ncid, strcat(strcpy(buffer, "x_"), grid_name), x, &dim_ids[2]));

  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, strcat(strcpy(buffer, grid_name), ".cla"), NC_DOUBLE, 3, dim_ids,
      &cla_var_id));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, cla_var_id, "units", strlen("degrees"), "degrees"));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, strcat(strcpy(buffer, grid_name), ".clo"), NC_DOUBLE, 3, dim_ids,
      &clo_var_id));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, clo_var_id, "units", strlen("degrees"), "degrees"));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, strcat(strcpy(buffer, grid_name), ".lat"), NC_DOUBLE, 2,
      &dim_ids[1], &lat_var_id));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, lat_var_id, "units", strlen("degrees"), "degrees"));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, strcat(strcpy(buffer, grid_name), ".lon"), NC_DOUBLE, 2,
      &dim_ids[1], &lon_var_id));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, lon_var_id, "units", strlen("degrees"), "degrees"));

  YAC_HANDLE_ERROR(nc_enddef(ncid));

  YAC_HANDLE_ERROR(nc_put_var_double(ncid, cla_var_id, cla));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, clo_var_id, clo));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, lat_var_id, lat));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, lon_var_id, lon));

  YAC_HANDLE_ERROR(nc_close(ncid));

  if (create_file) {
    yac_nc_create(mask_filename, NC_CLOBBER, &ncid);
  } else {
    yac_nc_open(mask_filename, NC_WRITE, &ncid);
    YAC_HANDLE_ERROR(nc_redef(ncid));
  }

  YAC_HANDLE_ERROR(
    nc_def_dim(
      ncid, strcat(strcpy(buffer, "y_"), grid_name), y, &dim_ids[0]));
  YAC_HANDLE_ERROR(
    nc_def_dim(
      ncid, strcat(strcpy(buffer, "x_"), grid_name), x, &dim_ids[1]));

  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, strcat(strcpy(buffer, grid_name), ".msk"), NC_INT, 2, dim_ids,
      &msk_var_id));

  YAC_HANDLE_ERROR(nc_enddef(ncid));

  int * msk = calloc(x * y, sizeof(*msk));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, msk_var_id, msk));
  free(msk);

  YAC_HANDLE_ERROR(nc_close(ncid));
}
