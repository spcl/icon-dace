/**
 * @file toy_woa_ocn.c
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

#include "yac_config.h"

#undef VERBOSE

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "utils.h"
#include "fields.h"
#include "yac_interface.h"
#include "geometry.h"
#include "read_woa_data.h"
#include "vtk_output.h"
#include "generate_reg2d.h"

/* ------------------------------------------------- */

/* For simplicity we define the same 4 fields that are in the
 * coupling configuration */

#include "toy_common.h"

#define STR_USAGE "Usage: %s -c configFilename -x num_cells_x -y num_cells_y\n"
#define YAC_ASSERT_ARGS(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, size_t * num_cells_x, size_t * num_cells_y);

int main (int argc, char *argv[]) {

#ifdef YAC_NETCDF_ENABLED

  // Initialisation of MPI

  MPI_Init (0, NULL);

  /* The initialisation phase includes the reading of the
   * coupling configuration */
#ifdef VERBOSE
  printf (". main: calling yac_cinit\n");
#endif

  double tic, toc, time;
  double time_min, time_max, time_ave;
  double time_min_acc, time_max_acc, time_ave_acc;

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  char const * configFilename = "toy_atm_ocn.yaml"; // default configuration file
  size_t num_cells_x = 360; // default horizontal resolution
  size_t num_cells_y = 180; // default vertical resolution
  parse_arguments(argc, argv, &configFilename, &num_cells_x, &num_cells_y);
  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  /* The usual component definition, here for two sequential components on the same process */

#ifdef VERBOSE
  printf (". main: calling yac_cdef_comp\n");
#endif

  int comp_id;
  char * comp_name = "WOA";

  yac_cdef_comp ( comp_name, &comp_id );
#ifdef VERBOSE
  printf ( ". main: defined %s with local comp ID %i \n", "toy-woa-ocean", comp_id );
#endif

  MPI_Comm local_comm;

  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;

  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);


  int point_id;

  int field_ids[4];
  int grid_id;

  /* Grid definition for WOA 2009 data set  */

  int total_nbr_cells[2];
  int num_procs[2];
  total_nbr_cells[0] = (int)num_cells_x;
  total_nbr_cells[1] = (int)num_cells_y;

  generate_reg2d_decomp(total_nbr_cells, size, num_procs);

  int block_pos[2];
  int block_size[2];
  int nbr_cells[2];
  int cyclic[2];
  int nbr_points[2];
  block_pos[0] = rank%num_procs[0];
  block_pos[1] = rank/num_procs[0];
  block_size[0] = (total_nbr_cells[0] + num_procs[0] - 1)/num_procs[0];
  block_size[1] = (total_nbr_cells[1] + num_procs[1] - 1)/num_procs[1];
  nbr_cells[0] = MIN(block_size[0], total_nbr_cells[0] - block_size[0] * block_pos[0]);
  nbr_cells[1] = MIN(block_size[1], total_nbr_cells[1] - block_size[1] * block_pos[1]);
  cyclic[0] = num_procs[0] == 1;
  cyclic[1] = 0;

  // halo
  if (!cyclic[0])
   nbr_cells[0] += 2;
  if (num_procs[1] > 1) {

    nbr_cells[1] += 1;

    if ((block_pos[1] != 0) && (block_pos[1] != num_procs[1]-1))
      nbr_cells[1] += 1;
  }

  if (cyclic[0])
    nbr_points[0] = nbr_cells[0];
  else
    nbr_points[0] = nbr_cells[0]+1;

  nbr_points[1] = nbr_cells[1] + 1;

  // the grid

  double * global_x_vertices = xmalloc((total_nbr_cells[0]+3) * sizeof(*global_x_vertices));
  double * global_y_vertices = xmalloc((total_nbr_cells[1]+1) * sizeof(*global_y_vertices));;
  double * x_vertices;
  double * y_vertices;
  double * x_points;
  double * y_points;

  for (int i = 0; i < total_nbr_cells[0]+3; ++i)
    global_x_vertices[i] = 0.0 + (2.0 * M_PI / ((double)total_nbr_cells[0])) * (double)(i-1);
  for (int i = 0; i < total_nbr_cells[1]; ++i)
    global_y_vertices[i] = -M_PI_2 + (M_PI / ((double)total_nbr_cells[1])) * (double)i;
  global_y_vertices[total_nbr_cells[1]] = M_PI_2;

  if (cyclic[0])
    x_vertices = x_points = global_x_vertices + 1;
  else
    x_vertices = x_points = global_x_vertices + block_size[0] * block_pos[0];

  y_vertices = y_points = global_y_vertices + block_size[1] * block_pos[1];

  if (block_pos[1] != 0) {
    y_vertices -= 1;
    y_points -= 1;
  }

  // the data and mask

  unsigned num_cells = (unsigned)(nbr_cells[0] * nbr_cells[1]);

  double * global_temperature;
  double * global_salinity;

  // salinity data


  {
    char * fileName = "salinity_monthly_1deg.nc";
    char * woaFieldName = "s_an";
    struct fieldMetadata fieldInfo;

    int woaFile = open_woa_output ( fileName );
    read_woa_dimensions ( woaFile, woaFieldName, &fieldInfo );
    global_salinity = get_woa_memory ( fieldInfo);
    read_woa_timestep_level ( woaFile, global_salinity, fieldInfo, 1, 1 );
    close_woa_output ( woaFile );
  }

  // temperature data

  {
    char * fileName = "temperature_monthly_1deg.nc";
    char * woaFieldName = "t_an";
    struct fieldMetadata fieldInfo;

    int woaFile = open_woa_output ( fileName );
    read_woa_dimensions ( woaFile, woaFieldName, &fieldInfo );
    global_temperature = get_woa_memory ( fieldInfo);
    read_woa_timestep_level ( woaFile, global_temperature, fieldInfo, 1, 1 );
    close_woa_output ( woaFile );
  }

#ifdef VERBOSE
  printf (". main: got woa data\n");
#endif

  double * salinity = xmalloc(num_cells * sizeof(*salinity));
  double * temperature = xmalloc(num_cells * sizeof(*temperature));

  unsigned iii = 0;
  unsigned ioff = block_size[0] * block_pos[0];
  unsigned joff = block_size[1] * block_pos[1];

  int * cell_mask = xmalloc (num_cells * sizeof(*cell_mask));

  for ( int j = 0; j < nbr_cells[1]; ++j ) {
    for ( int i = 0; i < nbr_cells[0]; ++i ) {
      unsigned ii = ( joff + j ) * total_nbr_cells[0]  + ioff + i;
      salinity[iii] = global_salinity[ii];
      temperature[iii] = global_temperature[ii];

      if ( temperature[iii] > 100.0 ) {
        cell_mask[iii] = 0;
        temperature[iii] = 999.0;
        salinity[iii] = 999.0;
      }
      else cell_mask[iii] = 1;
      ++iii;
    }
  }

  free_woa_memory (global_salinity);
  free_woa_memory (global_temperature);

  // now get everything into yac

  yac_cdef_grid_reg2d ( "grid2", nbr_points, cyclic, x_vertices, y_vertices, &grid_id);

  int * global_global_cell_id = xmalloc(total_nbr_cells[1] * (total_nbr_cells[0] + 2) * sizeof(*global_global_cell_id));
  int * global_global_cell_id_rank = xmalloc(total_nbr_cells[1] * (total_nbr_cells[0] + 2) * sizeof(*global_global_cell_id_rank));
  int * global_global_corner_id = xmalloc((total_nbr_cells[1] + 1) * (total_nbr_cells[0] + 3) * sizeof(*global_global_corner_id));

  for (int j = 0, id = 0; j < total_nbr_cells[1]; ++j) {

    for (int i = 1; i <= total_nbr_cells[0]; ++i) {

      global_global_cell_id[j * (total_nbr_cells[0] + 2) + i] = id;
      global_global_cell_id_rank[j * (total_nbr_cells[0] + 2) + i] =
         (id%total_nbr_cells[0]) / block_size[0] + ((id/total_nbr_cells[0]) / block_size[1]) * num_procs[0];
      id++;
    }
  }
  for (int i = 0; i < total_nbr_cells[1]; ++i) {
    global_global_cell_id[i * (total_nbr_cells[0] + 2) + 0] = global_global_cell_id[i * (total_nbr_cells[0] + 2) + total_nbr_cells[0]];
    global_global_cell_id_rank[i * (total_nbr_cells[0] + 2) + 0] = global_global_cell_id_rank[i * (total_nbr_cells[0] + 2) + total_nbr_cells[0]];
    global_global_cell_id[i * (total_nbr_cells[0] + 2) + total_nbr_cells[0]+1] = global_global_cell_id[i * (total_nbr_cells[0] + 2) + 1];
    global_global_cell_id_rank[i * (total_nbr_cells[0] + 2) + total_nbr_cells[0]+1] = global_global_cell_id_rank[i * (total_nbr_cells[0] + 2) + 1];
  }

  if (num_procs[0] == 1) {

    for (int j = 0, id = 0; j < total_nbr_cells[1] + 1; ++j) {

      for (int i = 1; i <= total_nbr_cells[0]; ++i)
        global_global_corner_id[j * (total_nbr_cells[0] + 3) + i] = id++;
    }

  } else {

    for (int j = 0, id = 0; j < total_nbr_cells[1] + 1; ++j) {

      for (int i = 1; i <= total_nbr_cells[0] + 1; ++i)
        global_global_corner_id[j * (total_nbr_cells[0] + 3) + i] = id++;
    }

    for (int i = 0; i < total_nbr_cells[1] + 1; ++i) {

      global_global_corner_id[i * (total_nbr_cells[0] + 3) + 0] = global_global_corner_id[i * (total_nbr_cells[0] + 3) + total_nbr_cells[0] + 1];
      global_global_corner_id[i * (total_nbr_cells[0] + 3) + total_nbr_cells[0] + 2] = global_global_corner_id[i * (total_nbr_cells[0] + 3) + 1];
    }
  }

  for (int i = 0; i < (total_nbr_cells[1])*(total_nbr_cells[0] + 2); ++i)
    if(global_global_cell_id_rank[i] == rank) global_global_cell_id_rank[i] = -1;

  int * local_global_cell_id = xmalloc(nbr_cells[1] * nbr_cells[0] * sizeof(*local_global_cell_id));
  int * local_global_cell_id_rank = xmalloc(nbr_cells[1] * nbr_cells[0] * sizeof(*local_global_cell_id_rank));
  int * local_global_corner_id = xmalloc(nbr_points[1] * nbr_points[0] * sizeof(*local_global_corner_id));
  int * local_global_corner_id_rank = xmalloc(nbr_points[1] * nbr_points[0] * sizeof(*local_global_corner_id_rank));

  int offset_x, offset_y;

  if (cyclic[0])
    offset_x = 1;
  else
    offset_x = block_size[0] * block_pos[0];
  offset_y = block_size[1] * block_pos[1] - (block_pos[1] != 0);

  for (int j = 0; j < nbr_cells[1]; ++j) {
    for (int i = 0; i < nbr_cells[0]; ++i) {
      local_global_cell_id[j * nbr_cells[0] + i] = global_global_cell_id[(j+offset_y) * (total_nbr_cells[0] + 2) + i+offset_x];
      local_global_cell_id_rank[j * nbr_cells[0] + i] = global_global_cell_id_rank[(j+offset_y) * (total_nbr_cells[0] + 2) +i+offset_x];
    }
  }

  free(global_global_cell_id_rank);
  free(global_global_cell_id);

  for (int j = 0; j < nbr_points[1]; ++j)
    for (int i = 0; i < nbr_points[0]; ++i)
      local_global_corner_id[j * nbr_points[0] + i] = global_global_corner_id[(j+offset_y) * (total_nbr_cells[0] + 3) +i+offset_x];

   free(global_global_corner_id);

  for (int j = 0; j < nbr_points[1]; ++j)
    for (int i = 0; i < nbr_points[0]; ++i)
      local_global_corner_id_rank[j * nbr_points[0] + i] = -1;

  for (int i = 0; i < nbr_points[0]-2; ++i) {
    local_global_corner_id_rank[0 * nbr_points[0] + i] = local_global_cell_id_rank[0 * nbr_cells[0] + i];
    local_global_corner_id_rank[(nbr_points[1]-1) * nbr_points[0] + i] = local_global_cell_id_rank[(nbr_cells[1]-1) * nbr_cells[0] + i];
  }
  for (int j = 0; j < nbr_points[1]-2; ++j) {
    local_global_corner_id_rank[j * nbr_points[0] + 0] = local_global_cell_id_rank[j * nbr_cells[0] + 0];
    local_global_corner_id_rank[j * nbr_points[0] + nbr_points[0]-1] = local_global_cell_id_rank[j * nbr_cells[0] + nbr_cells[0]-1];
  }
  local_global_corner_id_rank[(nbr_points[1]-1)* nbr_points[0] + nbr_points[0]-1] = local_global_cell_id_rank[(nbr_cells[1]-1) * nbr_cells[0] + nbr_cells[0]-1];
  local_global_corner_id_rank[(nbr_points[1]-2)* nbr_points[0] + 0] = local_global_cell_id_rank[(nbr_cells[1]-2) * nbr_cells[0] + 0];
  local_global_corner_id_rank[0 * nbr_points[0] + nbr_points[0]-2] = local_global_cell_id_rank[0 * nbr_cells[0] + nbr_cells[0]-2];
  local_global_corner_id_rank[(nbr_points[1]-2) * nbr_points[0] + nbr_points[0]-1] = local_global_cell_id_rank[(nbr_cells[1]-2) * nbr_cells[0] + nbr_cells[0]-1];
  local_global_corner_id_rank[(nbr_points[1]-1) * nbr_points[0] + nbr_points[0]-2] = local_global_cell_id_rank[(nbr_cells[1]-1) * nbr_cells[0] + nbr_cells[0]-2];

  int * cell_core_mask = xmalloc(nbr_cells[1] * nbr_cells[0] * sizeof(*local_global_corner_id_rank));
  int * corner_core_mask = xmalloc(nbr_points[1] * nbr_points[0] * sizeof(*local_global_corner_id_rank));

  for (int i = 0; i < nbr_cells[1]*nbr_cells[0]; ++i)
    cell_core_mask[i] = local_global_cell_id_rank[i] == -1;
  for (int i = 0; i < nbr_points[1]*nbr_points[0]; ++i)
    corner_core_mask[i] = local_global_corner_id_rank[i] == -1;

  yac_cset_global_index(local_global_cell_id, CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, CELL, grid_id);
  yac_cset_global_index(local_global_corner_id, CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, CORNER, grid_id);

  yac_cdef_points_reg2d ( grid_id, nbr_cells, CELL, x_points, y_points, &point_id );

  yac_cset_mask ( cell_mask, point_id );

  /* Field definition for the second component (ocean) */

  for ( int i = 0; i < num_fields; i++ )
    yac_cdef_field(
      fieldName[i], comp_id, &point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
      &field_ids[i] );

  toc=MPI_Wtime();
  time = toc-tic;

  MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
  MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
  MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

  time_ave /= (double) size;

  if ( rank == 0 )
    printf ("toy-woa-ocean: Time for initialisation %f %f %f \n", time_min, time_max, time_ave );

  /* Search. */

#ifdef VERBOSE
  printf (". main: calling yac_cenddef\n");
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  yac_cenddef ( );

  toc=MPI_Wtime();
  time = toc-tic;

  MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
  MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
  MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

  time_ave /= (double) size;

  if ( rank == 0 )
    printf ("toy-icon-ocean: Time for search         %f %f %f \n", time_min, time_max, time_ave );

  double * cell_data_field = xmalloc(num_cells * sizeof(*cell_data_field));
  double * cell_out_conserv_data = xmalloc(num_cells * sizeof(*cell_out_conserv_data));
  double * cell_out_hcsbb_data = xmalloc(num_cells * sizeof(*cell_out_hcsbb_data));
  double * cell_in_conserv_data = xmalloc(num_cells * sizeof(*cell_in_conserv_data));
  double * cell_in_hcsbb_data = xmalloc(num_cells * sizeof(*cell_in_hcsbb_data));
  double * cell_in_conserv_err_abs = xmalloc(num_cells * sizeof(*cell_in_conserv_err_abs));
  double * cell_in_hcsbb_err_abs = xmalloc(num_cells * sizeof(*cell_in_hcsbb_err_abs));
  double * cell_in_conserv_err_rel = xmalloc(num_cells * sizeof(*cell_in_conserv_err_rel));
  double * cell_in_hcsbb_err_rel = xmalloc(num_cells * sizeof(*cell_in_hcsbb_err_rel));

  int err;
  int info;

  unsigned * cell_corners = xmalloc(num_cells * 4 * sizeof(*cell_corners));
  unsigned * num_points_per_cell = xmalloc(num_cells * sizeof(*num_points_per_cell));

  for (unsigned i = 0; i < num_cells; ++i) {

    {
      unsigned x_index, y_index;

      y_index = i / nbr_cells[0];
      x_index = i - y_index * nbr_cells[0];

      if (!cyclic[0]) {

        cell_corners[i*4+0] =  y_index      * (nbr_cells[0] + 1) + x_index;
        cell_corners[i*4+1] =  y_index      * (nbr_cells[0] + 1) + x_index + 1;
        cell_corners[i*4+2] = (y_index + 1) * (nbr_cells[0] + 1) + x_index + 1;
        cell_corners[i*4+3] = (y_index + 1) * (nbr_cells[0] + 1) + x_index;

      } else {

        cell_corners[i*4+0] = y_index * nbr_cells[0] + x_index;
        if (x_index + 1 != nbr_cells[0]) {
          cell_corners[i*4+1] =  y_index      * nbr_cells[0] + x_index + 1;
          cell_corners[i*4+2] = (y_index + 1) * nbr_cells[0] + x_index + 1;
        } else {
          cell_corners[i*4+1] =  y_index      * nbr_cells[0];
          cell_corners[i*4+2] = (y_index + 1) * nbr_cells[0];
        }
        cell_corners[i*4+3] = (y_index + 1) * nbr_cells[0] + x_index;
      }
    }
    num_points_per_cell[i] = 4;
  }

  for (unsigned i = 0; i < num_cells; ++i) {
    cell_in_conserv_data[i] = -999;
    cell_in_hcsbb_data[i]   = -999;
  }

  double * point_data = xmalloc(nbr_points[0] * nbr_points[1] * 3 * sizeof(*point_data));
  for (int i = 0; i < nbr_points[1]; ++i) {
    for (int j = 0; j < nbr_points[0]; ++j) {

      LLtoXYZ(x_vertices[j], y_vertices[i], &point_data[3*(i * nbr_points[0] + j)]);
    }
  }

  time_min_acc = 0.0;
  time_max_acc = 0.0;
  time_ave_acc = 0.0;

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  {
    double *point_set_data[1];
    double **collection_data[1] = {point_set_data};

    point_set_data[0] = salinity;
    yac_cput(field_ids[2], 1, collection_data, &info, &err);
    YAC_ASSERT(info, "check coupling period")
    point_set_data[0] = temperature;
    yac_cput(field_ids[3], 1, collection_data, &info, &err);
    YAC_ASSERT(info, "check coupling period")
  }

  {
    double *collection_data[1];

    collection_data[0] = cell_in_conserv_data;
    yac_cget(field_ids[0], 1, collection_data, &info, &err);
    YAC_ASSERT(info, "check coupling period")
    collection_data[0] = cell_in_hcsbb_data;
    yac_cget(field_ids[1], 1, collection_data, &info, &err);
    YAC_ASSERT(info, "check coupling period")
  }

  toc=MPI_Wtime();
  time = toc-tic;

  MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
  MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
  MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

  time_ave /= (double) size;

  if ( rank == 0 ) {
    time_min_acc += time_min;
    time_max_acc += time_max;
    time_ave_acc += time_ave;
  }

  //----------------------------------------------------------
  // compute error
  //----------------------------------------------------------

  for ( unsigned i = 0; i < num_cells; ++i ) {

    cell_in_conserv_err_abs[i] = fabs(cell_in_conserv_data[i] - cell_data_field[i]);
    cell_in_hcsbb_err_abs[i] = fabs(cell_in_hcsbb_data[i] - cell_data_field[i]);
    cell_in_conserv_err_rel[i] = fabs(cell_in_conserv_err_abs[i] / cell_data_field[i]);
    cell_in_hcsbb_err_rel[i] = fabs(cell_in_hcsbb_err_abs[i] / cell_data_field[i]);
  }

  //----------------------------------------------------------
  // write field to vtk output file
  //----------------------------------------------------------

  char vtk_filename[32];

  sprintf(vtk_filename, "toy-woa-ocean_out_%d_%d.vtk", rank, 1);

  VTK_FILE *vtk_file = vtk_open(vtk_filename, "toy-woa-ocean_out");

  vtk_write_point_data(vtk_file, point_data, nbr_points[0]*nbr_points[1]);

  vtk_write_cell_data(vtk_file, cell_corners, num_points_per_cell, num_cells);

  vtk_write_point_scalars_int(vtk_file, corner_core_mask,
                              nbr_points[0]*nbr_points[1], "corner_core_mask");
  vtk_write_point_scalars_int(vtk_file, local_global_corner_id,
                              nbr_points[0]*nbr_points[1], "global_corner_id");
  vtk_write_cell_scalars_int(vtk_file, cell_core_mask, num_cells,
                             "cell_core_mask");
  vtk_write_cell_scalars_int(vtk_file, local_global_cell_id, num_cells,
                             "global_cell_id");
  vtk_write_cell_scalars_int(vtk_file, cell_mask, num_cells, "cell_mask");

  vtk_write_cell_scalars_double(vtk_file, cell_in_conserv_data, num_cells, "cell_in_conserv");
  vtk_write_cell_scalars_double(vtk_file, cell_in_hcsbb_data, num_cells, "cell_in_hcsbb");
  vtk_write_cell_scalars_double(vtk_file, salinity, num_cells, "salinity");
  vtk_write_cell_scalars_double(vtk_file, temperature, num_cells, "temperature");
  vtk_close(vtk_file);

  for ( unsigned i = 0; i < num_cells; ++i ) {
      cell_out_conserv_data[i] = cell_in_conserv_data[i];
      cell_out_hcsbb_data[i]   = cell_in_hcsbb_data[i];
      cell_in_conserv_data[i] = -999;
      cell_in_hcsbb_data[i]   = -999;
  }

  if ( rank == 0 )
    printf ("toy-icon-atmosphere: Time for ping-pong      %f %f %f \n", time_min_acc, time_max_acc, time_ave_acc );

  yac_cfinalize();

  MPI_Finalize();

  free(cell_in_hcsbb_err_rel);
  free(cell_in_conserv_err_rel);
  free(cell_in_hcsbb_err_abs);
  free(cell_in_conserv_err_abs);
  free(cell_in_conserv_data);
  free(cell_in_hcsbb_data);
  free(cell_out_conserv_data);
  free(cell_out_hcsbb_data);

  free(num_points_per_cell);
  free(cell_corners);
  free(point_data);
  free(corner_core_mask);
  free(local_global_corner_id_rank);
  free(local_global_corner_id);
  free(cell_core_mask);
  free(local_global_cell_id_rank);
  free(local_global_cell_id);
  free(global_x_vertices);
  free(global_y_vertices);
  free(cell_mask);

  return EXIT_SUCCESS;

#else

  fprintf(stderr, "ERROR: toy-woa-ocean requires NETCDF\n");

  return EXIT_FAILURE;

#endif // YAC_NETCDF_ENABLED
}

static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, size_t * num_cells_x, size_t * num_cells_y) {

  int opt;
  while ((opt = getopt(argc, argv, "c:x:y:")) != -1) {
    YAC_ASSERT((opt == 'c') || (opt == 'x') || (opt == 'y'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
      case 'x':
        *num_cells_x = atoi(optarg);
        YAC_ASSERT_ARGS(*num_cells_x > 0, "invalid horizontal resolution");
        break;
      case 'y':
        *num_cells_y = atoi(optarg);
        YAC_ASSERT_ARGS(*num_cells_y > 0, "invalid vertical resolution");
        break;
    }
  }
}
