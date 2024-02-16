/**
 * @file toy_output.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *         Nils-Arne Dreier <dreier@dkrz.de>
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

#include "yac_config.h"

#if defined YAC_NETCDF_ENABLED

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include "yac_interface.h"

#define HANDLE_ERROR(exp) \
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
  HANDLE_ERROR(nc_create(filename, NC_CLOBBER, &ncid));

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
  HANDLE_ERROR(nc_def_dim(ncid, "lat", nbr_cells[1], &lat_dimid));
  HANDLE_ERROR(nc_def_dim(ncid, "lon", nbr_cells[0], &lon_dimid));
  HANDLE_ERROR(nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dimid));
  HANDLE_ERROR(nc_def_dim(ncid, "z", collection_size, &z_dimid));

  free(x_cells);
  free(y_cells);

  int varid;
  int dimids[] = {t_dimid, z_dimid, lat_dimid, lon_dimid};
  HANDLE_ERROR(nc_def_var(ncid, field_name, NC_DOUBLE, 4, dimids, &varid));

  int interp_stack_config_id;
  yac_cget_interp_stack_config(&interp_stack_config_id);
  yac_cadd_interp_stack_config_nnn(interp_stack_config_id, YAC_NNN_AVG, 1, 1.0);

  yac_cdef_couple( source_comp, source_grid, field_name,
                   comp_name, grid_name, field_name,
                   dt, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                   interp_stack_config_id, 0, 0);

  yac_cfree_interp_stack_config(interp_stack_config_id);

  HANDLE_ERROR(nc_enddef(ncid));

  yac_cenddef();

  double* data = malloc(nbr_cells[0]*nbr_cells[1]*sizeof(data)*collection_size);
  int info, ierror;
  int time_counter = 0;
  const char* end = yac_cget_end_datetime ( );
  const char* t = yac_cget_field_datetime(field_id);
  while(1){
    printf("receiving %s at %s\n", field_name, t);
    yac_cget_ ( field_id,
               collection_size,
               data,
               &info,
               &ierror );
    size_t start[] = {time_counter, 0, 0, 0};
    size_t count[] = {1, collection_size, nbr_cells[1], nbr_cells[0]};
    HANDLE_ERROR(nc_put_vara_double (ncid, varid, start, count, data ));
    t = yac_cget_field_datetime(field_id);
    if (strcmp(end, t) == 0)
      break;
    time_counter++;
  }
  printf("done, wrote %d timesteps for %s\n", time_counter, field_name);

  free(data);

  HANDLE_ERROR(nc_close(ncid));

  yac_cfinalize();

  return 0;
}

#else
#include <stdlib.h>
#include <stdio.h>
int main () {
  printf ("Example requires compiling with NetCDF.\n");
  return EXIT_FAILURE;
}
#endif
