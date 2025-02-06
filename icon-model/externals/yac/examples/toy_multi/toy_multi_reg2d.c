// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

// #define VERBOSE

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "yac.h"
#include "yac_utils.h"

#include "toy_multi_common.h"

// redefine YAC assert macros
#undef YAC_ASSERT
#define STR_USAGE "Usage: %s -x num_cells_x -y num_cells_y\n"
#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, size_t * num_cells_x, size_t * num_cells_y);

int main (int argc, char *argv[]) {

  double tic;

  // Initialisation of MPI

  MPI_Init (0, NULL);

  MPI_Barrier(MPI_COMM_WORLD);
  tic=MPI_Wtime();

  size_t num_cells_x = 256; // default horizontal resolution
  size_t num_cells_y = 128; // default vertical resolution
  parse_arguments(argc, argv, &num_cells_x, &num_cells_y);
  yac_cinit ();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime("2008-03-09T16:05:07", "2008-03-10T16:05:07");

  int comp_id;
  char * comp_name = "ECHAM";

  yac_cdef_comp(comp_name, &comp_id);

  MPI_Comm comp_comm;

  yac_cget_comp_comm(comp_id, &comp_comm);

  int rank, comp_rank, comp_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_rank(comp_comm,&comp_rank);
  MPI_Comm_size(comp_comm,&comp_size);
  MPI_Comm_free(&comp_comm);

  int total_nbr_cells[2] = {num_cells_x, num_cells_y};
  int num_procs[2];

  yac_generate_reg2d_decomp(total_nbr_cells, comp_size, num_procs);

  int block_pos[2];
  int block_size[2];
  int nbr_cells[2];
  int cyclic[2];
  int nbr_points[2];
  block_pos[0] = comp_rank%num_procs[0];
  block_pos[1] = comp_rank/num_procs[0];
  block_size[0] = (num_cells_x + num_procs[0] - 1)/num_procs[0];
  block_size[1] = (num_cells_y + num_procs[1] - 1)/num_procs[1];
  nbr_cells[0] = MIN(block_size[0], (int)num_cells_x - block_size[0] * block_pos[0]);
  nbr_cells[1] = MIN(block_size[1], (int)num_cells_y - block_size[1] * block_pos[1]);
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

  double * global_x_vertices = malloc((num_cells_x+3) * sizeof(*global_x_vertices));
  double * global_y_vertices = malloc((num_cells_y+1) * sizeof(*global_y_vertices));
  double * x_vertices;
  double * y_vertices;
  double * x_corner_points;
  double * y_corner_points;

  for (size_t i = 0; i < num_cells_x+3; ++i)
    global_x_vertices[i] =
      0.0 + (2.0 * M_PI / ((double)num_cells_x)) * (double)(i-1);
  for (size_t i = 0; i < num_cells_y; ++i)
    global_y_vertices[i] =
      -M_PI_2 + (M_PI / ((double)num_cells_y)) * (double)i;
  global_y_vertices[num_cells_y] = M_PI_2;

  if (cyclic[0])
    x_vertices = x_corner_points = global_x_vertices + 1;
  else
    x_vertices = x_corner_points = global_x_vertices + block_size[0] * block_pos[0];

  y_vertices = y_corner_points = global_y_vertices + block_size[1] * block_pos[1];

  if (block_pos[1] != 0) {
    y_vertices -= 1;
    y_corner_points -= 1;
  }

  double * x_cell_points = malloc(nbr_cells[0] * sizeof(*x_cell_points));
  double * y_cell_points = malloc(nbr_cells[1] * sizeof(*y_cell_points));

  for (int i = 0; i < nbr_cells[0]; ++i)
    x_cell_points[i] = (x_corner_points[i] + x_corner_points[i + 1]) * 0.5;
  for (int i = 0; i < nbr_cells[1]; ++i)
    y_cell_points[i] = (y_corner_points[i] + y_corner_points[i + 1]) * 0.5;

  int grid_id;
  char grid_name[256];
  sprintf(grid_name, "%s_grid", comp_name);

  yac_cdef_grid_reg2d(
    grid_name, nbr_points, cyclic, x_vertices, y_vertices, &grid_id);

  int * global_global_cell_id =
    malloc(num_cells_y * (num_cells_x + 2) * sizeof(*global_global_cell_id));
  int * global_global_corner_id =
    malloc((num_cells_y + 1) * (num_cells_x + 3) * sizeof(*global_global_corner_id));

  for (size_t j = 0, id = 0; j < num_cells_y; ++j)
    for (size_t i = 1; i <= num_cells_x; ++i)
      global_global_cell_id[j * (num_cells_x + 2) + i] = id++;
  for (size_t i = 0; i < num_cells_y; ++i) {
    global_global_cell_id[i * (num_cells_x + 2) + 0] =
      global_global_cell_id[i * (num_cells_x + 2) + num_cells_x];
    global_global_cell_id[i * (num_cells_x + 2) + num_cells_x+1] =
      global_global_cell_id[i * (num_cells_x + 2) + 1];
  }

  if (num_procs[0] == 1) {

    for (size_t j = 0, id = 0; j < num_cells_y + 1; ++j)
      for (size_t i = 1; i <= num_cells_x; ++i)
        global_global_corner_id[j * (num_cells_x + 3) + i] = id++;

  } else {

    for (size_t j = 0, id = 0; j < num_cells_y + 1; ++j)
      for (size_t i = 0; i < num_cells_x; ++i)
        global_global_corner_id[j * (num_cells_x + 3) + (i + 1)] = id++;

    for (size_t i = 0; i < num_cells_y + 1; ++i) {

      global_global_corner_id[i * (num_cells_x + 3) + 0] =
        global_global_corner_id[i * (num_cells_x + 3) + num_cells_x];
      global_global_corner_id[i * (num_cells_x + 3) + num_cells_x + 1] =
        global_global_corner_id[i * (num_cells_x + 3) + 1];
      global_global_corner_id[i * (num_cells_x + 3) + num_cells_x + 2] =
        global_global_corner_id[i * (num_cells_x + 3) + 2];
    }
  }

  int * local_global_cell_id = malloc(nbr_cells[1] * nbr_cells[0] * sizeof(*local_global_cell_id));
  int * cell_core_mask = malloc(nbr_cells[1] * nbr_cells[0] * sizeof(*cell_core_mask));
  int * local_global_corner_id = malloc(nbr_points[1] * nbr_points[0] * sizeof(*local_global_corner_id));
  int * corner_core_mask = malloc(nbr_points[1] * nbr_points[0] * sizeof(*corner_core_mask));

  int offset_x, offset_y;

  if (cyclic[0])
    offset_x = 1;
  else
    offset_x = block_size[0] * block_pos[0];
  offset_y = block_size[1] * block_pos[1] - (block_pos[1] != 0);

  for (int j = 0; j < nbr_cells[1]; ++j)
    for (int i = 0; i < nbr_cells[0]; ++i)
      local_global_cell_id[j * nbr_cells[0] + i] =
        global_global_cell_id[(j+offset_y) * (num_cells_x + 2) + i+offset_x];

  free(global_global_cell_id);

  for (int j = 0; j < nbr_points[1]; ++j)
    for (int i = 0; i < nbr_points[0]; ++i)
      local_global_corner_id[j * nbr_points[0] + i] =
        global_global_corner_id[(j+offset_y) * (num_cells_x + 3) +i+offset_x];

  free(global_global_corner_id);

  for (int j = 0; j < nbr_cells[1]; ++j)
    for (int i = 0; i < nbr_cells[0]; ++i)
      cell_core_mask[j * nbr_cells[0] + i] = 1;
  if (num_procs[1] > 1) {
    if (block_pos[1] > 0)
      for (int i = 0; i < nbr_cells[0]; ++i)
        cell_core_mask[i] = 0;
    if (block_pos[1] < num_procs[1] - 1)
      for (int i = 0; i < nbr_cells[0]; ++i)
        cell_core_mask[(nbr_cells[1] - 1) * nbr_cells[0] + i] = 0;
  }
  if (num_procs[0] > 1) {
    for (int j = 0; j < nbr_cells[1]; ++j) {
      cell_core_mask[j * nbr_cells[0]] = 0;
      cell_core_mask[j * nbr_cells[0] + nbr_cells[0] - 1] = 0;
    }
  }

  for (int j = 0; j < nbr_points[1]; ++j)
    for (int i = 0; i < nbr_points[0]; ++i)
      corner_core_mask[j * nbr_points[0] + i] = 1;
  if (num_procs[1] > 1) {
    if (block_pos[1] > 0)
      for (int i = 0; i < nbr_points[0]; ++i)
        corner_core_mask[i] = 0;
    if (block_pos[1] < num_procs[1] - 1)
      for (int i = 0; i < nbr_points[0]; ++i)
        corner_core_mask[(nbr_points[1] - 1) * nbr_points[0] + i] = 0;
  }
  if (num_procs[0] > 1) {
    for (int j = 0; j < nbr_points[1]; ++j) {
      corner_core_mask[j * nbr_points[0]] = 0;
      corner_core_mask[j * nbr_points[0] + nbr_points[0] - 1] = 0;
    }
  }

  yac_cset_global_index(local_global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(local_global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);

  int cell_point_id, corner_point_id;

  yac_cdef_points_reg2d(
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_cell_points, y_cell_points,
    &cell_point_id);
  yac_cdef_points_reg2d(
    grid_id, nbr_points, YAC_LOCATION_CORNER, x_corner_points, y_corner_points,
    &corner_point_id);

  unsigned num_cells = nbr_cells[0] * nbr_cells[1];
  unsigned num_corners = nbr_points[0] * nbr_points[1];

  double * cell_point_data = malloc(num_cells * sizeof(*cell_point_data));
  double * corner_point_data = malloc(num_corners * sizeof(*corner_point_data));

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

    cell_point_data[i] = yac_test_func(curr_x, curr_y);
  }

  for (unsigned i = 0; i < num_corners; ++i)
    corner_point_data[i] = yac_test_func(x_corner_points[i%nbr_points[0]],
                                     y_corner_points[i/nbr_points[0]]);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    printf(
      "toy_multi_common(%s): Time for initialisation %f\n",
      comp_name, MPI_Wtime()-tic);

  // initialize VTK file
  double * point_data =
    malloc(nbr_points[0] * nbr_points[1] * 3 * sizeof(*point_data));
  for (int i = 0; i < nbr_points[1]; ++i)
    for (int j = 0; j < nbr_points[0]; ++j)
      LLtoXYZ(
        x_vertices[j], y_vertices[i], &point_data[3*(i * nbr_points[0] + j)]);

  char vtk_filename[32];

  sprintf(vtk_filename, "%s_out_%d.vtk", comp_name, comp_rank);

  YAC_VTK_FILE *vtk_file = yac_vtk_open(vtk_filename, comp_name);
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

  int ret =
    run_toy_multi_common(
      comp_name, comp_id, grid_id, cell_point_id, corner_point_id,
      cell_point_data, corner_point_data, vtk_file);

  free(num_points_per_cell);
  free(cell_corners);
  free(point_data);
  free(corner_point_data);
  free(cell_point_data);
  free(corner_core_mask);
  free(local_global_corner_id);
  free(cell_core_mask);
  free(local_global_cell_id);
  free(x_cell_points);
  free(y_cell_points);
  free(global_x_vertices);
  free(global_y_vertices);

  return ret;
}

static void parse_arguments(
  int argc, char ** argv, size_t * num_cells_x, size_t * num_cells_y) {

  int opt;
  while ((opt = getopt(argc, argv, "x:y:")) != -1) {
    YAC_ASSERT((opt == 'c') || (opt == 'x') || (opt == 'y'), "invalid command argument")
    switch (opt) {
      default:
      case 'x':
        *num_cells_x = atoi(optarg);
        YAC_ASSERT(*num_cells_x > 0, "invalid horizontal resolution");
        break;
      case 'y':
        *num_cells_y = atoi(optarg);
        YAC_ASSERT(*num_cells_y > 0, "invalid vertical resolution");
        break;
    }
  }
}
