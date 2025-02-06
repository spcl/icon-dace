// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

// #define VERBOSE

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "yac.h"
#include "yac_core.h"
#include "yac_utils.h"

#include "toy_multi_common.h"

/* ------------------------------------------------- */

struct {
  char const * name;
  int const * location;
} interpolations[] =
  {
    {.name = "conserv",
     .location = &YAC_LOCATION_CELL},
    {.name = "2nd_conserv",
     .location = &YAC_LOCATION_CELL},
    {.name = "avg",
     .location = &YAC_LOCATION_CORNER},
    {.name = "hcsbb",
     .location = &YAC_LOCATION_CELL},
    {.name = "rbf",
     .location = &YAC_LOCATION_CELL},
  };
size_t nbr_interpolations = sizeof(interpolations) / sizeof(interpolations[0]);

int run_toy_multi_common (
  char const * comp_name, int comp_id, int grid_id,
  int cell_point_id, int corner_point_id,
  double * cell_point_data, double * corner_point_data,
  YAC_VTK_FILE * vtk_file) {

  int info, err;
  double tic;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Barrier(MPI_COMM_WORLD);
  tic=MPI_Wtime();

  // get list of all registered components
  int nbr_comps = yac_cget_nbr_comps();
  char const ** comp_names = malloc((size_t)nbr_comps * sizeof(*comp_names));
  yac_cget_comp_names(nbr_comps, comp_names);
  size_t local_comp_idx = SIZE_MAX;
  for (int i = 0; i < nbr_comps; ++i) comp_names[i] = strdup(comp_names[i]);

  // register fields
  int * field_ids =
    malloc((size_t)nbr_comps * nbr_interpolations * sizeof(*field_ids));
  for (int comp_idx = 0, i = 0; comp_idx < nbr_comps; ++comp_idx) {
    if (!strcmp(comp_name, comp_names[comp_idx]))
      local_comp_idx = (size_t)comp_idx;
    for (size_t interp_idx = 0; interp_idx < nbr_interpolations;
         ++interp_idx, ++i) {
      char field_name[256];
      sprintf(
        field_name, "%s_%s_field_out", interpolations[interp_idx].name,
        comp_names[comp_idx]);
      int * point_id_ptr;
      YAC_ASSERT(
        (*(interpolations[interp_idx].location) == YAC_LOCATION_CELL) ||
        (*(interpolations[interp_idx].location) == YAC_LOCATION_CORNER),
        "invalid location");
      if (*(interpolations[interp_idx].location) == YAC_LOCATION_CELL)
        point_id_ptr = &cell_point_id;
      else
        point_id_ptr = &corner_point_id;
      yac_cdef_field(
        field_name, comp_id, point_id_ptr, 1, 1, "2",
        YAC_TIME_UNIT_SECOND, &field_ids[i]);
    }
  }

  // define couples
  if (rank == 0) {
    for (size_t interp_idx = 0; interp_idx < nbr_interpolations;
         ++interp_idx) {

      int interp_config;
      yac_cget_interp_stack_config(&interp_config);
      char const * interp_name = interpolations[interp_idx].name;
      if (!strcmp(interp_name, "conserv")) {
         yac_cadd_interp_stack_config_conservative(
           interp_config, 1, 0, 0, YAC_CONSERV_DESTAREA);
      } else if (!strcmp(interp_name, "2nd_conserv")) {
         yac_cadd_interp_stack_config_conservative(
           interp_config, 2, 0, 0, YAC_CONSERV_DESTAREA);
      } else if (!strcmp(interp_name, "avg")) {
        yac_cadd_interp_stack_config_average(
          interp_config, YAC_AVG_BARY, 1);
      } else if (!strcmp(interp_name, "hcsbb")) {
        yac_cadd_interp_stack_config_hcsbb(interp_config);
      } else if (!strcmp(interp_name, "rbf")) {
        yac_cadd_interp_stack_config_nnn(
          interp_config, YAC_NNN_RBF, 9,
          YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT, 1.487973e+01);
      } else {
        fprintf(
          stderr, "unsupported interpolation type: \"%s\"\n", interp_name);
        exit(EXIT_FAILURE);
      }
      yac_cadd_interp_stack_config_fixed(interp_config, -1.0);

      for (int src_comp_idx = 0; src_comp_idx < nbr_comps; ++src_comp_idx) {

        char src_grid_name[256];
        sprintf(src_grid_name, "%s_grid", comp_names[src_comp_idx]);

        for (int tgt_comp_idx = 0; tgt_comp_idx < nbr_comps; ++tgt_comp_idx) {

          if (src_comp_idx == tgt_comp_idx) continue;

          char tgt_grid_name[256];
          sprintf(tgt_grid_name, "%s_grid", comp_names[tgt_comp_idx]);

          char field_name[256];
          sprintf(
            field_name, "%s_%s_field_out", interp_name,
            comp_names[src_comp_idx]);

          yac_cdef_couple(
            comp_names[src_comp_idx], src_grid_name, field_name,
            comp_names[tgt_comp_idx], tgt_grid_name, field_name,
            "2", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
            interp_config, 0, 0);
        }
      }

      yac_cfree_interp_stack_config(interp_config);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    printf(
      "toy_multi_common(%s): Time for field/interpolation definition %f\n",
      comp_name, MPI_Wtime()-tic);

  // finalize initialisation phase

  MPI_Barrier(MPI_COMM_WORLD);
  tic=MPI_Wtime();

  yac_cenddef ( );

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    printf(
      "toy_multi_common(%s): Time for search %f\n",
      comp_name, MPI_Wtime()-tic);

  // initialise data

  double ** field_buffers =
    malloc((size_t)nbr_comps * nbr_interpolations * sizeof(*field_buffers));
  for (int comp_idx = 0, i = 0; comp_idx < nbr_comps; ++comp_idx) {
    for (size_t interp_idx = 0; interp_idx < nbr_interpolations;
         ++interp_idx, ++i) {

      if ((size_t)comp_idx == local_comp_idx) {
        field_buffers[i] =
          (*(interpolations[interp_idx].location) == YAC_LOCATION_CELL)?
            cell_point_data:corner_point_data;
      } else {
        int point_id =
          (*(interpolations[interp_idx].location) == YAC_LOCATION_CELL)?
            cell_point_id:corner_point_id;
        size_t data_size = yac_cget_points_size(point_id);
        field_buffers[i] = malloc(data_size * sizeof(**field_buffers));
        for (size_t j = 0; j < data_size; ++j)
          field_buffers[i][j] = -10.0;
      }
    }
  }

  // exchange/interpolate fields

  MPI_Barrier(MPI_COMM_WORLD);
  tic=MPI_Wtime();

  for (size_t interp_idx = 0, offset = local_comp_idx * nbr_interpolations;
       interp_idx < nbr_interpolations; ++interp_idx) {

    double *point_set_data[1] = {field_buffers[offset + interp_idx]};
    double **collection_data[1] = {point_set_data};

    yac_cput(field_ids[offset + interp_idx], 1, collection_data, &info, &err);
  }

  for (int comp_idx = 0, i = 0; comp_idx < nbr_comps; ++comp_idx) {
    for (size_t interp_idx = 0; interp_idx < nbr_interpolations;
         ++interp_idx, ++i) {

      if ((size_t)comp_idx == local_comp_idx) continue;

      double *collection_data[1];

      collection_data[0] = field_buffers[i];
      yac_cget(field_ids[i], 1, collection_data, &info, &err);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    printf(
      "toy_multi_common(%s): Time for exchange %f\n",
      comp_name, MPI_Wtime()-tic);

  // write fields to vtk output file

  yac_vtk_write_cell_scalars_double(
    vtk_file, cell_point_data,
    yac_cget_grid_size(YAC_LOCATION_CELL, grid_id), "cell_out");
  yac_vtk_write_point_scalars_double(
    vtk_file, corner_point_data,
    yac_cget_grid_size(YAC_LOCATION_CORNER, grid_id), "corner_out");

  for (int comp_idx = 0, i = 0; comp_idx < nbr_comps; ++comp_idx) {
    for (size_t interp_idx = 0; interp_idx < nbr_interpolations;
         ++interp_idx, ++i) {

      if ((size_t)comp_idx != local_comp_idx) {

        char field_name[256];
        sprintf(
          field_name, "%s_in_%s", interpolations[interp_idx].name,
          comp_names[comp_idx]);

        if (*(interpolations[interp_idx].location) == YAC_LOCATION_CELL)
          yac_vtk_write_cell_scalars_double(
            vtk_file, field_buffers[i],
            yac_cget_grid_size(
              *(interpolations[interp_idx].location), grid_id), field_name);
        else
          yac_vtk_write_point_scalars_double(
            vtk_file, field_buffers[i],
            yac_cget_grid_size(
              *(interpolations[interp_idx].location), grid_id), field_name);
      }
    }
  }

  // finalize

  yac_vtk_close(vtk_file);
  yac_cfinalize();
  MPI_Finalize();

  // clean up

  for (int comp_idx = 0, i = 0; comp_idx < nbr_comps; ++comp_idx)
    for (size_t interp_idx = 0; interp_idx < nbr_interpolations;
         ++interp_idx, ++i)
      if ((size_t)comp_idx != local_comp_idx) free(field_buffers[i]);
  free(field_buffers);

  return EXIT_SUCCESS;
}
