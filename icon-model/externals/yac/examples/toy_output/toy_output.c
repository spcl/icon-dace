// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include "yac.h"

#define YAC_HANDLE_ERROR(exp) \
  do { \
    int handle_error_status = (exp); \
    if (handle_error_status != NC_NOERR) { \
      fprintf(stderr, "Error: %s\n", nc_strerror(handle_error_status)); \
      exit(handle_error_status); \
    } \
  } while(0)


int main(int argc, char** argv){
  assert(argc == 4);
  const char* source_comp = argv[1];
  const char* source_grid = argv[2];
  const char* field_name = argv[3];

  yac_cinit();

  char filename[32];
  int ncid;
  sprintf(filename, "%s.nc", field_name);
  YAC_HANDLE_ERROR(nc_create(filename, NC_CLOBBER, &ncid));

  printf("Writing file %s", filename);

  int comp_id;
  char comp_name[256];
  sprintf(comp_name, "toy_output_%s_%s_%s", source_comp, source_grid, field_name);
  yac_cdef_comp(comp_name, &comp_id);

  int grid_id;
  int nbr_vertices[] = {360, 181};
  int cyclic[] = {1, 0};
  double* x_vertices = malloc(nbr_vertices[0]*sizeof(*x_vertices));
  for(int i = 0; i<nbr_vertices[0]; ++i){
    x_vertices[i] = -M_PI + i*2*M_PI/nbr_vertices[0];
  }
  double* y_vertices = malloc(nbr_vertices[1]*sizeof(*y_vertices));
  for(int i = 0; i<nbr_vertices[1]; ++i){
    y_vertices[i] = -0.5*M_PI + i*M_PI/nbr_vertices[1];
  }

  const char * grid_name = "toy_output_grid";
  yac_cdef_grid_reg2d ( grid_name,
                        nbr_vertices,
                        cyclic,
                        x_vertices,
                        y_vertices,
                        &grid_id);

  free(x_vertices);
  free(y_vertices);

  int nbr_cells[] = {nbr_vertices[0], nbr_vertices[1]-1};
  double* x_cells = malloc(nbr_cells[0]*sizeof(*x_cells));
  for(int i = 0; i<nbr_cells[0]; ++i){
    x_cells[i] = -M_PI + (((double)i) + 0.5)*2*M_PI/nbr_cells[0];
  }
  double* y_cells = malloc(nbr_cells[1]*sizeof(*y_cells));
  for(int i = 0; i<nbr_cells[1]; ++i){
    y_cells[i] = -0.5*M_PI + (((double)i) + 0.5)*M_PI/nbr_vertices[1];
  }
  int point_id;
  yac_cdef_points_reg2d( grid_id,
                         nbr_cells,
                         YAC_LOCATION_CELL,
                         x_cells,
                         y_cells,
                         &point_id );

  yac_csync_def();

  const char* dt = yac_cget_field_timestep(source_comp, source_grid, field_name);
  printf("toy_output: timestep for %s is %s\n", field_name, dt);

  int collection_size = yac_cget_field_collection_size(source_comp, source_grid, field_name);
  printf("toy_output: collection_size for %s is %d\n", field_name, collection_size);

  int field_id;
  yac_cdef_field ( field_name,
                   comp_id,
                   &point_id,
                   1,
                   /*collection_size*/ collection_size,
                   dt,
                   YAC_TIME_UNIT_ISO_FORMAT,
                   &field_id );

  int lat_dimid, lon_dimid, t_dimid, z_dimid;
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "lat", nbr_cells[1], &lat_dimid));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "lon", nbr_cells[0], &lon_dimid));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dimid));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "z", collection_size, &z_dimid));

  free(x_cells);
  free(y_cells);

  int varid;
  int dimids[] = {t_dimid, z_dimid, lat_dimid, lon_dimid};
  YAC_HANDLE_ERROR(nc_def_var(ncid, field_name, NC_DOUBLE, 4, dimids, &varid));

  int interp_stack_config_id;
  yac_cget_interp_stack_config(&interp_stack_config_id);
  yac_cadd_interp_stack_config_nnn(
    interp_stack_config_id, YAC_NNN_AVG, 1, 0.0, 1.0);

  yac_cdef_couple( source_comp, source_grid, field_name,
                   comp_name, grid_name, field_name,
                   dt, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                   interp_stack_config_id, 0, 1);

  yac_cfree_interp_stack_config(interp_stack_config_id);

  YAC_HANDLE_ERROR(nc_enddef(ncid));

  yac_cenddef();

  double* data = malloc(nbr_cells[0]*nbr_cells[1]*sizeof(data)*collection_size);
  int info, ierror;
  int time_counter = 0;
  while(1){
    const char* t = yac_cget_field_datetime(field_id);
    printf("receiving %s at %s\n", field_name, t);
    yac_cget_ ( field_id,
               collection_size,
               data,
               &info,
               &ierror );
    size_t start[] = {time_counter, 0, 0, 0};
    size_t count[] = {1, collection_size, nbr_cells[1], nbr_cells[0]};
    YAC_HANDLE_ERROR(nc_put_vara_double (ncid, varid, start, count, data ));
    yac_cget_action(field_id, &info);
    if (info == YAC_ACTION_OUT_OF_BOUND)
      break;
    time_counter++;
  }
  printf("done, wrote %d timesteps for %s\n", time_counter, field_name);

  free(data);

  YAC_HANDLE_ERROR(nc_close(ncid));

  yac_cfinalize();

  return 0;
}
