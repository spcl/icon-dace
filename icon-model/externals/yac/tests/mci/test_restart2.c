// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "yac.h"
#include "utils_common.h"
#include "geometry.h"
#include "generate_cubed_sphere.h"
#include "test_function.h"

#define DUMMY_VALUE (-1337.0)

static void run_setup(char local_comp_id, char comp_groups[2][2]);

int main (void) {

  // This test checks the restarting of YAC using varying MPI communicators.
  // This setup could also be executed with a single YAC initilisation.

  MPI_Init(NULL, NULL);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_size != 15) {
    fputs("wrong number of processes (has to be 15)\n", stderr);
    exit(EXIT_FAILURE);
  }

  char local_comp_id;
  if (comm_rank < 6)       local_comp_id = 'A';
  else if (comm_rank < 10) local_comp_id = 'B';
  else if (comm_rank < 13) local_comp_id = 'C';
  else                     local_comp_id = 'D';

  // couple A-B and C-D (each in their own world comm)
  run_setup(local_comp_id, (char[2][2]){{'A','B'},{'C','D'}});

  // couple A-C and B-D (each in their own world comm)
  run_setup(local_comp_id, (char[2][2]){{'A','C'},{'B','D'}});

  // couple A-D and B-C (each in their own world comm)
  run_setup(local_comp_id, (char[2][2]){{'A','D'},{'B','C'}});

  MPI_Finalize();

  exit(EXIT_SUCCESS);
}

static void def_grid(
  char const * grid_name, int comm_rank, int comm_size,
  int * grid_id, int * cell_point_id, int local_test_func_idx,
  int remote_test_func_idx, double ** field_data_out,
  double ** ref_field_data_in, size_t * field_data_size) {

  unsigned n = 50;

  unsigned nbr_vertices;
  unsigned nbr_cells;
  unsigned * num_vertices_per_cell;
  unsigned * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * cell_core_mask;
  int * corner_core_mask;
  int * global_cell_id;
  int * global_corner_id;

  yac_generate_part_cube_grid_information(
    n, &nbr_vertices, &nbr_cells, &num_vertices_per_cell, &cell_to_vertex,
    &x_vertices, &y_vertices, &x_cells, &y_cells, &global_cell_id,
    &cell_core_mask, &global_corner_id, &corner_core_mask,
    comm_rank, comm_size);

  yac_cdef_grid_unstruct(
    grid_name, nbr_vertices, nbr_cells, (int*)num_vertices_per_cell,
    x_vertices, y_vertices, (int*)cell_to_vertex, grid_id);

  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, *grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, *grid_id);

  yac_cdef_points_unstruct(
    *grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells,
    cell_point_id);

  double (*test_func[4])(double, double) =
    {&yac_test_ana_fcos,
     &yac_test_vortex,
     &yac_test_gulfstream,
     &yac_test_harmonic};

  *field_data_out = xmalloc(nbr_cells * sizeof(**field_data_out));
  for (unsigned i = 0; i < nbr_cells; ++i)
    (*field_data_out)[i] =
      (cell_core_mask[i])?
        (*(test_func[local_test_func_idx]))(x_cells[i], y_cells[i]):
        DUMMY_VALUE;
  *ref_field_data_in = xmalloc(nbr_cells * sizeof(**ref_field_data_in));
  for (unsigned i = 0; i < nbr_cells; ++i)
    (*ref_field_data_in)[i] =
      (cell_core_mask[i])?
        (*(test_func[remote_test_func_idx]))(x_cells[i], y_cells[i]):
        DUMMY_VALUE;
  *field_data_size = nbr_cells;
  free(cell_core_mask);
  free(corner_core_mask);
  free(global_cell_id);
  free(global_corner_id);
  free(x_vertices);
  free(y_vertices);
  free(x_cells);
  free(y_cells);
  free(num_vertices_per_cell);
  free(cell_to_vertex);
}

static void run_component(
  char local_comp_id, char remote_comp_id, MPI_Comm comm) {

  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  // init YAC
  yac_cinit_comm(comm);
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

  // add configuration (one out-going field)
  yac_cdef_datetime("2008-03-09T16:05:07", "2008-03-10T16:05:07");
  char local_comp_name[16],  local_grid_name[16];
  char remote_comp_name[16], remote_grid_name[16];
  char field_out_name[16], field_in_name[16];
  sprintf(local_comp_name, "comp%c", local_comp_id);
  sprintf(local_grid_name, "grid%c", local_comp_id);
  sprintf(remote_comp_name, "comp%c", remote_comp_id);
  sprintf(remote_grid_name, "grid%c", remote_comp_id);
  sprintf(field_out_name, "field_%c_to_%c", local_comp_id, remote_comp_id);
  sprintf(field_in_name, "field_%c_to_%c", remote_comp_id, local_comp_id);
  int interp_stack_config_id;
  yac_cget_interp_stack_config(&interp_stack_config_id);
  yac_cadd_interp_stack_config_nnn(
    interp_stack_config_id, YAC_NNN_AVG, 1, 0.0, -1.0);
  yac_cdef_couple(
    local_comp_name, local_grid_name, field_out_name,
    remote_comp_name, remote_grid_name, field_out_name,
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
    interp_stack_config_id, 0, 0);
  yac_cfree_interp_stack_config(interp_stack_config_id);

  // define component
  int comp_id;
  yac_cdef_comp(local_comp_name, &comp_id);

  MPI_Comm comp_comm;
  int comp_rank, comp_size;
  yac_cget_comp_comm(comp_id, &comp_comm);
  MPI_Comm_rank(comp_comm, &comp_rank);
  MPI_Comm_size(comp_comm, &comp_size);
  MPI_Comm_free(&comp_comm);

  // define grid
  int grid_id, cell_point_id;
  double * field_data_out, * ref_field_data_in;
  size_t field_data_size;
  def_grid(local_grid_name, comp_rank, comp_size, &grid_id, &cell_point_id,
           (int)(local_comp_id - 'A'), (int)(remote_comp_id - 'A'),
           &field_data_out, &ref_field_data_in, &field_data_size);

  // define fields
  int field_out_id, field_in_id;
  yac_cdef_field(
    field_out_name, comp_id, &cell_point_id, 1, 1, "1",
    YAC_TIME_UNIT_SECOND, &field_out_id);
  yac_cdef_field(
    field_in_name, comp_id, &cell_point_id, 1, 1, "1",
    YAC_TIME_UNIT_SECOND, &field_in_id);

  yac_cenddef();

  double * field_data_in = xmalloc(field_data_size * sizeof(*field_data_in));

  // do some ping-pongs
  for (int t = 0; t < 10; ++t) {

    {
      int info, err;
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};
      point_set_data[0] = field_data_out;
      yac_cput(field_out_id, 1, collection_data, &info, &err);
    }

    {
      for (size_t i = 0; i < field_data_size; ++i)
        field_data_in[i] = DUMMY_VALUE;

      int info, err;
      double *collection_data[1] = {field_data_in};
      yac_cget(field_in_id, 1, collection_data, &info, &err);

      for (size_t i = 0; i < field_data_size; ++i) {
        if(fabs(field_data_in[i] - ref_field_data_in[i]) > 1e-3) {
          fputs("data data_mismatch\n", stderr);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  free(field_data_in);
  free(field_data_out);
  free(ref_field_data_in);

  yac_cfinalize();

  MPI_Barrier(comm);
}


static void run_setup(char local_comp_id, char comp_groups[2][2]) {

  int group_idx = (local_comp_id == comp_groups[1][0]) ||
                  (local_comp_id == comp_groups[1][1]);

  MPI_Comm group_comm;
  MPI_Comm_split(MPI_COMM_WORLD, group_idx, 0, &group_comm);

  char remote_comp_id =
    comp_groups[group_idx][comp_groups[group_idx][0] == local_comp_id];

  run_component(local_comp_id, remote_comp_id, group_comm);

  MPI_Comm_free(&group_comm);
}
