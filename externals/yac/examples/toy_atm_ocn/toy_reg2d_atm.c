// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#undef VERBOSE

#define T106

/* supported is
 *
 *     T21, T31, T42, T63, T85, T106, T127, T159, T255
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "yac.h"
#include "yac_core.h"
#include "yac_utils.h"

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

#if defined T21

#define NUM_CELLS_X (64)
#define NUM_CELLS_Y (32)

#elif defined T31

#define NUM_CELLS_X (96)
#define NUM_CELLS_Y (48)

#elif defined T42

#define NUM_CELLS_X (128)
#define NUM_CELLS_Y (64)

#elif defined T63

#define NUM_CELLS_X (192)
#define NUM_CELLS_Y (96)

#elif defined T85

#define NUM_CELLS_X (256)
#define NUM_CELLS_Y (128)

#elif defined T106

#define NUM_CELLS_X (320)
#define NUM_CELLS_Y (160)

#elif defined T127

#define NUM_CELLS_X (384)
#define NUM_CELLS_Y (192)

#elif defined T159

#define NUM_CELLS_X (480)
#define NUM_CELLS_Y (240)

#elif defined T255

#define NUM_CELLS_X (768)
#define NUM_CELLS_Y (384)

#else

#define NUM_CELLS_X (256)
#define NUM_CELLS_Y (128)

#endif

  char const * configFilename = "toy_atm_ocn.yaml"; // default configuration file
  size_t num_cells_x = NUM_CELLS_X; // default horizontal resolution
  size_t num_cells_y = NUM_CELLS_Y; // default vertical resolution
  parse_arguments(argc, argv, &configFilename, &num_cells_x, &num_cells_y);
  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  /* The usual component definition, here for two sequential components on the same process */

#ifdef VERBOSE
  printf (". main: calling yac_cdef_comp\n");
#endif

  int comp_id;
  char * comp_name = "ATMOS";

  yac_cdef_comp ( comp_name, &comp_id );
#ifdef VERBOSE
  printf ( ". main: defined %s with local comp ID %i \n", "toy-reg2d-atmosphere", comp_id );
#endif

  MPI_Comm local_comm;

  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;

  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);


  int point_id;

  int field_ids[4];
  int grid_id;

  /* Grid definition for toy-reg2d-atmosphere */

  int total_nbr_cells[2];
  int num_procs[2];
  total_nbr_cells[0] = (int)num_cells_x;
  total_nbr_cells[1] = (int)num_cells_y;

  yac_generate_reg2d_decomp(total_nbr_cells, size, num_procs);

  int block_pos[2];
  int block_size[2];
  int nbr_cells[2];
  int cyclic[2];
  int nbr_points[2];
  block_pos[0] = rank%num_procs[0];
  block_pos[1] = rank/num_procs[0];
  block_size[0] = (NUM_CELLS_X + num_procs[0] - 1)/num_procs[0];
  block_size[1] = (NUM_CELLS_Y + num_procs[1] - 1)/num_procs[1];
  nbr_cells[0] = MIN(block_size[0], NUM_CELLS_X - block_size[0] * block_pos[0]);
  nbr_cells[1] = MIN(block_size[1], NUM_CELLS_Y - block_size[1] * block_pos[1]);
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

  double * global_x_vertices = malloc((NUM_CELLS_X+3) * sizeof(*global_x_vertices));
  double * global_y_vertices = malloc((NUM_CELLS_Y+1) * sizeof(*global_y_vertices));;
  double * x_vertices;
  double * y_vertices;
  double * x_points;
  double * y_points;

  for (int i = 0; i < NUM_CELLS_X+3; ++i)
    global_x_vertices[i] =
      0.0 + (2.0 * M_PI / ((double)NUM_CELLS_X)) * (double)(i-1);
  for (int i = 0; i < NUM_CELLS_Y; ++i)
    global_y_vertices[i] =
      -M_PI_2 + (M_PI / ((double)NUM_CELLS_Y)) * (double)i;
  global_y_vertices[NUM_CELLS_Y] = M_PI_2;

  if (cyclic[0])
    x_vertices = x_points = global_x_vertices + 1;
  else
    x_vertices = x_points = global_x_vertices + block_size[0] * block_pos[0];

  y_vertices = y_points = global_y_vertices + block_size[1] * block_pos[1];

  if (block_pos[1] != 0) {
    y_vertices -= 1;
    y_points -= 1;
  }

  yac_cdef_grid_reg2d ( "grid1", nbr_points, cyclic, x_vertices, y_vertices, &grid_id);

  int * global_global_cell_id = malloc(NUM_CELLS_Y * (NUM_CELLS_X + 2) * sizeof(*global_global_cell_id));
  int * global_global_cell_id_rank = malloc(NUM_CELLS_Y * (NUM_CELLS_X + 2) * sizeof(*global_global_cell_id_rank));
  int * global_global_corner_id = malloc((NUM_CELLS_Y + 1) * (NUM_CELLS_X + 3) * sizeof(*global_global_corner_id));

  for (int j = 0, id = 0; j < NUM_CELLS_Y; ++j) {

    for (int i = 1; i <= NUM_CELLS_X; ++i) {

      global_global_cell_id[j * (NUM_CELLS_X + 2) + i] = id;
      global_global_cell_id_rank[j * (NUM_CELLS_X + 2) + i] =
         (id%NUM_CELLS_X) / block_size[0] + ((id/NUM_CELLS_X) / block_size[1]) * num_procs[0];
      id++;
    }
  }
  for (int i = 0; i < NUM_CELLS_Y; ++i) {
    global_global_cell_id[i * (NUM_CELLS_X + 2) + 0] = global_global_cell_id[i * (NUM_CELLS_X + 2) + NUM_CELLS_X];
    global_global_cell_id_rank[i * (NUM_CELLS_X + 2) + 0] = global_global_cell_id_rank[i * (NUM_CELLS_X + 2) + NUM_CELLS_X];
    global_global_cell_id[i * (NUM_CELLS_X + 2) + NUM_CELLS_X+1] = global_global_cell_id[i * (NUM_CELLS_X + 2) + 1];
    global_global_cell_id_rank[i * (NUM_CELLS_X + 2) + NUM_CELLS_X+1] = global_global_cell_id_rank[i * (NUM_CELLS_X + 2) + 1];
  }

  if (num_procs[0] == 1) {

    for (int j = 0, id = 0; j < NUM_CELLS_Y + 1; ++j) {

      for (int i = 1; i <= NUM_CELLS_X; ++i)
        global_global_corner_id[j * (NUM_CELLS_X + 3) + i] = id++;
    }

  } else {

    for (int j = 0, id = 0; j < NUM_CELLS_Y + 1; ++j) {

      for (int i = 1; i <= NUM_CELLS_X + 1; ++i)
        global_global_corner_id[j * (NUM_CELLS_X + 3) + i] = id++;
    }

    for (int i = 0; i < NUM_CELLS_Y + 1; ++i) {

      global_global_corner_id[i * (NUM_CELLS_X + 3) + 0] = global_global_corner_id[i * (NUM_CELLS_X + 3) + NUM_CELLS_X + 1];
      global_global_corner_id[i * (NUM_CELLS_X + 3) + NUM_CELLS_X + 2] = global_global_corner_id[i * (NUM_CELLS_X + 3) + 1];
    }
  }

  for (int i = 0; i < (NUM_CELLS_Y)*(NUM_CELLS_X + 2); ++i)
    if(global_global_cell_id_rank[i] == rank) global_global_cell_id_rank[i] = -1;

  int * local_global_cell_id = malloc(nbr_cells[1] * nbr_cells[0] * sizeof(*local_global_cell_id));
  int * local_global_cell_id_rank = malloc(nbr_cells[1] * nbr_cells[0] * sizeof(*local_global_cell_id_rank));
  int * local_global_corner_id = malloc(nbr_points[1] * nbr_points[0] * sizeof(*local_global_corner_id));
  int * local_global_corner_id_rank = malloc(nbr_points[1] * nbr_points[0] * sizeof(*local_global_corner_id_rank));

  int offset_x, offset_y;

  if (cyclic[0])
    offset_x = 1;
  else
    offset_x = block_size[0] * block_pos[0];
  offset_y = block_size[1] * block_pos[1] - (block_pos[1] != 0);

  for (int j = 0; j < nbr_cells[1]; ++j) {
    for (int i = 0; i < nbr_cells[0]; ++i) {
      local_global_cell_id[j * nbr_cells[0] + i] = global_global_cell_id[(j+offset_y) * (NUM_CELLS_X + 2) + i+offset_x];
      local_global_cell_id_rank[j * nbr_cells[0] + i] = global_global_cell_id_rank[(j+offset_y) * (NUM_CELLS_X + 2) +i+offset_x];
    }
  }

  free(global_global_cell_id_rank);
  free(global_global_cell_id);

  for (int j = 0; j < nbr_points[1]; ++j)
    for (int i = 0; i < nbr_points[0]; ++i)
      local_global_corner_id[j * nbr_points[0] + i] = global_global_corner_id[(j+offset_y) * (NUM_CELLS_X + 3) +i+offset_x];

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

  int * cell_core_mask = malloc(nbr_cells[1] * nbr_cells[0] * sizeof(*local_global_corner_id_rank));
  int * corner_core_mask = malloc(nbr_points[1] * nbr_points[0] * sizeof(*local_global_corner_id_rank));

  for (int i = 0; i < nbr_cells[1]*nbr_cells[0]; ++i)
    cell_core_mask[i] = local_global_cell_id_rank[i] == -1;
  for (int i = 0; i < nbr_points[1]*nbr_points[0]; ++i)
    corner_core_mask[i] = local_global_corner_id_rank[i] == -1;

  yac_cset_global_index(local_global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(local_global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);

  yac_cdef_points_reg2d(
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_points, y_points, &point_id);

  /* Field definition for toy-reg2d-atmosphere */

  for ( int i = 0; i < num_fields; i++ )
    yac_cdef_field(
      fieldName[i], comp_id, &point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
      &field_ids[i]);

  toc=MPI_Wtime();
  time = toc-tic;

  MPI_Reduce ( &time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, local_comm);
  MPI_Reduce ( &time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, local_comm);
  MPI_Reduce ( &time, &time_ave, 1, MPI_DOUBLE, MPI_SUM, 0, local_comm);

  time_ave /= (double) size;

  if ( rank == 0 )
    printf ("toy-reg2d-atmosphere: Time for initialisation %f %f %f \n", time_min, time_max, time_ave );

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
    printf ("toy-icon-atmosphere: Time for search         %f %f %f \n", time_min, time_max, time_ave );

  unsigned num_cells = nbr_cells[0] * nbr_cells[1];

  double * cell_data_field = malloc(num_cells * sizeof(*cell_data_field));
  double * cell_out_conserv_data = malloc(num_cells * sizeof(*cell_out_conserv_data));
  double * cell_out_hcsbb_data = malloc(num_cells * sizeof(*cell_out_hcsbb_data));
  double * cell_in_conserv_data = malloc(num_cells * sizeof(*cell_in_conserv_data));
  double * cell_in_hcsbb_data = malloc(num_cells * sizeof(*cell_in_hcsbb_data));
  double * cell_in_conserv_err_abs = malloc(num_cells * sizeof(*cell_in_conserv_err_abs));
  double * cell_in_hcsbb_err_abs = malloc(num_cells * sizeof(*cell_in_hcsbb_err_abs));
  double * cell_in_conserv_err_rel = malloc(num_cells * sizeof(*cell_in_conserv_err_rel));
  double * cell_in_hcsbb_err_rel = malloc(num_cells * sizeof(*cell_in_hcsbb_err_rel));

  int err;
  int info;

  unsigned * cell_corners = malloc(num_cells * 4 * sizeof(*cell_corners));
  unsigned * num_points_per_cell = malloc(num_cells * sizeof(*num_points_per_cell));

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
        if (x_index + 1 != (unsigned)(nbr_cells[0])) {
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

    double curr_x, curr_y;

    curr_x = (x_vertices[cell_corners[i*4+0]%nbr_points[0]] +
              x_vertices[cell_corners[i*4+1]%nbr_points[0]]) * 0.5;
    curr_y = (y_vertices[cell_corners[i*4+1]/nbr_points[0]] +
              y_vertices[cell_corners[i*4+2]/nbr_points[0]]) * 0.5;

    cell_data_field[i] = yac_test_func(curr_x, curr_y);
  }

  for (unsigned i = 0; i < num_cells; ++i) {

    cell_out_conserv_data[i] = cell_data_field[i];
    cell_out_hcsbb_data[i] = cell_data_field[i];

    cell_in_conserv_data[i] = -999;
    cell_in_hcsbb_data[i]   = -999;
  }

  double * point_data = malloc(nbr_points[0] * nbr_points[1] * 3 * sizeof(*point_data));
  for (int i = 0; i < nbr_points[1]; ++i) {
    for (int j = 0; j < nbr_points[0]; ++j) {

      LLtoXYZ(x_vertices[j], y_vertices[i], &point_data[3*(i * nbr_points[0] + j)]);
    }
  }

  time_min_acc = 0.0;
  time_max_acc = 0.0;
  time_ave_acc = 0.0;

  for (int step = 0; step < num_steps; ++step) {

    MPI_Barrier(MPI_COMM_WORLD);

    tic=MPI_Wtime();

    {
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};

      point_set_data[0] = cell_out_conserv_data;
      yac_cput(field_ids[0], 1, collection_data, &info, &err);
      YAC_ASSERT(info, "check coupling period")

      point_set_data[0] = cell_out_hcsbb_data;
      yac_cput(field_ids[1], 1, collection_data, &info, &err);
      YAC_ASSERT(info, "check coupling period")
    }

    {
      double *collection_data[1];

      collection_data[0] = cell_in_conserv_data;
      yac_cget(field_ids[2], 1, collection_data, &info, &err);
      YAC_ASSERT(info, "check coupling period")

      collection_data[0] = cell_in_hcsbb_data;
      yac_cget(field_ids[3], 1, collection_data, &info, &err);
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

#if 0

    double f1err_1 = 0.0;
    double f1err_100 = 0.0;
    int if1err_1 = 0;
    int if1err_100 = 0;

    double f2err_1 = 0.0;
    double f2err_100 = 0.0;
    int if2err_1 = 0;
    int if2err_100 = 0;

    for ( int i = 0; i < num_cells; ++i ) {

      if ( cell_in_conserv_err_rel[i] < 100.0 ) {
    	  f1err_1 += cell_in_conserv_err_rel[i];
    	  if1err_1++;
      } else {
    	  f1err_100 += cell_in_conserv_err_rel[i];
    	  if1err_100++;
      }

      if ( cell_in_hcsbb_err_rel[i] < 100.0 ) {
    	  f2err_1 += cell_in_hcsbb_err_rel[i];
    	  if2err_1++;
      } else {
    	  f2err_100 += cell_in_hcsbb_err_rel[i];
    	  if2err_100++;
      }
    }

    if (if1err_1 > 0.0) f1err_1 = f1err_1 / (double)if1err_1;
    if (if1err_1 > 0.0) f2err_1 = f2err_1 / (double)if2err_1;

    if (if1err_100 > 0.0) f1err_100 = f1err_100 / (double)if1err_100;
    if (if2err_100 > 0.0) f2err_100 = f2err_100 / (double)if2err_100;

    printf ("avg. rel err. for 1st field %f\n", f1err_1);
    printf ("avg. rel err. for 2nd field %f\n", f2err_1);

    printf ("extreme avg. rel err. for 1st field %f\n", f1err_100);
    printf ("extreme avg. rel err. for 2nd field %f\n", f2err_100);
#endif

    //----------------------------------------------------------
    // write field to vtk output file
    //----------------------------------------------------------

    char vtk_filename[36];

    sprintf(vtk_filename, "toy-reg2d-atmosphere_out_%d_%d.vtk", rank, step);

    YAC_VTK_FILE *vtk_file =
      yac_vtk_open(vtk_filename, "toy-reg2d-atmosphere_out");

    yac_vtk_write_point_data(
      vtk_file, point_data, nbr_points[0]*nbr_points[1]);

    yac_vtk_write_cell_data(
      vtk_file, cell_corners, num_points_per_cell, num_cells);

    yac_vtk_write_point_scalars_int(
      vtk_file, corner_core_mask, nbr_points[0]*nbr_points[1],
      "corner_core_mask");
    yac_vtk_write_point_scalars_int(
      vtk_file, local_global_corner_id, nbr_points[0]*nbr_points[1],
      "global_corner_id");
    yac_vtk_write_cell_scalars_int(
      vtk_file, cell_core_mask, num_cells, "cell_core_mask");
    yac_vtk_write_cell_scalars_int(
      vtk_file, local_global_cell_id, num_cells, "global_cell_id");

    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_conserv_data, num_cells, "cell_in_conserv");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_hcsbb_data, num_cells, "cell_in_hcsbb");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_out_conserv_data, num_cells, "cell_out_conserv");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_out_hcsbb_data, num_cells, "cell_out_hcsbb");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_conserv_err_abs, num_cells, "cell_in_conserv_err_abs");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_hcsbb_err_abs, num_cells, "cell_in_hcsbb_err_abs");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_conserv_err_rel, num_cells, "cell_in_conserv_err_rel");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_in_hcsbb_err_rel, num_cells, "cell_in_hcsbb_err_rel");

    yac_vtk_close(vtk_file);

    for ( unsigned i = 0; i < num_cells; ++i ) {
      cell_out_conserv_data[i] = cell_in_conserv_data[i];
      cell_out_hcsbb_data[i]   = cell_in_hcsbb_data[i];
      cell_in_conserv_data[i] = -10;
      cell_in_hcsbb_data[i]   = -10;
    }
  }

  if ( rank == 0 ) {
    time_min_acc /= (double) num_steps;
    time_max_acc /= (double) num_steps;
    time_ave_acc /= (double) num_steps;
    printf ("toy-icon-atmosphere: Time for ping-pong      %f %f %f \n", time_min_acc, time_max_acc, time_ave_acc );
  }

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
  free(cell_data_field);

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

  return EXIT_SUCCESS;
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
