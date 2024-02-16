/**
 * @file test_dummy_coupling3.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
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

#define EXACT

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yac_interface.h"
#include "read_icon_grid.h"
#include "read_mpiom_grid.h"
#include "generate_cubed_sphere.h"
#include "test_function.h"
#include "utils.h"

#include "tests.h"

#define DUMMY_VALUE (-1337.0)

struct grid_data {

  char const * name;

  int nbr_vertices;
  int nbr_cells;
  int * num_vertices_per_cell;
  int * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * cell_mask;
  int * global_cell_id;
  int * cell_core_mask;
  int * global_corner_id;
  int * corner_core_mask;
};

struct intra_comp_yac_info {
  int yac_id;
  int is_comp_root;
  int out_field_id;
  int in_field_id;
};

struct comp_pair_yac_info {
  int yac_id;
  int nbr_cells[2];
  int out_field_ids[2];
  int in_field_ids[2];
};

static struct intra_comp_yac_info setup_intra_comp_yac(
  MPI_Comm comp_comm, char const * comp_name,
  int base_point_id, int regional_point_id);

static struct comp_pair_yac_info setup_comp_pair_yac(
  char const * comp_pair_names[2], int point_ids[2], int nbr_cells[2],
  int local_comp_indices[2]);

static struct grid_data generate_icon_grid_data(MPI_Comm comm);
static struct grid_data generate_cube_grid_data(MPI_Comm comm);
static struct grid_data generate_mpiom_grid_data(MPI_Comm comm);
static void grid_data_free(struct grid_data grid_data);

static void def_default_fields(
  int comp_id, char const * comp_name, int grid_id, int cell_point_id,
  int * in_field_ids, int * out_field_id);

/**
 * This tests checks the support for multiple parallel YAC instances.
 * There are three main components and one dummy component. The dummy component
 * does not take part in any coupling.
 *
 * All four components are registered in the default YAC instance.
 *
 * The three active components couple to each other.
 *
 * Each active component shares one process with each of the two other active
 * components.
 *
 * Each active component pair generates their own YAC instance.
 *
 * Each active component has its own YAC instance in which one of its processes
 * generates its own component.
 *
 * The whole setup is run twice in order to simulate the restarting of YAC
 * instances.
 */
int main (void) {
#ifdef YAC_NETCDF_ENABLED

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);

  if (global_size != 10) {
    fprintf(stderr, "Wrong number of processes (should be 10)\n");
    exit(EXIT_FAILURE);
  }

  // ranks 8 and 9 do not initialise YAC, but still need to call the
  // respective dummy initilialisation and finalise routines
  if (global_rank > 7) {

    yac_cinit_dummy();
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  // determine the local components
  char const * global_comp_names[] = {"comp_a", "comp_b", "comp_c", "dummy"};
  char const * local_comp_names[2];
  int num_local_comps;
  int is_dummy = 0;
  switch (global_rank) {
    default:
    case(0): {
      local_comp_names[0] = global_comp_names[0];
      num_local_comps = 1;
      break;
    }
    case(1): {
      local_comp_names[0] = global_comp_names[1];
      num_local_comps = 1;
      break;
    }
    case(2): {
      local_comp_names[0] = global_comp_names[2];
      num_local_comps = 1;
      break;
    }
    case(4): {
      local_comp_names[0] = global_comp_names[0];
      local_comp_names[1] = global_comp_names[1];
      num_local_comps = 2;
      break;
    }
    case(5): {
      local_comp_names[0] = global_comp_names[0];
      local_comp_names[1] = global_comp_names[2];
      num_local_comps = 2;
      break;
    }
    case(6): {
      local_comp_names[0] = global_comp_names[1];
      local_comp_names[1] = global_comp_names[2];
      num_local_comps = 2;
      break;
    }
    case(3):
    case(7): {
      local_comp_names[0] = global_comp_names[3];
      num_local_comps = 1;
      is_dummy = 1;
      break;
    }
  }

  // initialise default YAC instance
  yac_cinit ();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_yaml("coupling_test3_default.yaml");

  // register local components in default YAC instance
  int default_comp_ids[2];
  int default_num_comps = num_local_comps;
  yac_cdef_comps(local_comp_names, default_num_comps, default_comp_ids);

  MPI_Comm default_comp_comms[2];
  for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx)
    yac_cget_comp_comm(
      default_comp_ids[comp_idx], &(default_comp_comms[comp_idx]));

  // generate and register grid, global ids, core mask, and points
  int grid_ids[2], cell_point_ids[2];
  int regional_grid_id = -1, regional_point_id = -1;
  int nbr_vertices_regional = 0;
  struct grid_data comp_grid_data[2];
  if (!is_dummy) {
    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx) {
      switch(local_comp_names[comp_idx][5]) {
        case('a'): {
          comp_grid_data[comp_idx] =
          generate_icon_grid_data(default_comp_comms[comp_idx]);
          break;
        }
        case('b'): {
          comp_grid_data[comp_idx] =
            generate_cube_grid_data(default_comp_comms[comp_idx]);
          break;
        }
        case('c'): {
          comp_grid_data[comp_idx] =
            generate_mpiom_grid_data(default_comp_comms[comp_idx]);
          break;
        }
        default:
          exit(EXIT_FAILURE);
      };
      yac_cdef_grid_unstruct(
        comp_grid_data[comp_idx].name,
        comp_grid_data[comp_idx].nbr_vertices,
        comp_grid_data[comp_idx].nbr_cells,
        comp_grid_data[comp_idx].num_vertices_per_cell,
        comp_grid_data[comp_idx].x_vertices,
        comp_grid_data[comp_idx].y_vertices,
        comp_grid_data[comp_idx].cell_to_vertex,
        &(grid_ids[comp_idx]));
      yac_cset_global_index(
        comp_grid_data[comp_idx].global_corner_id, YAC_LOCATION_CORNER,
        grid_ids[comp_idx]);
      yac_cset_core_mask(
        comp_grid_data[comp_idx].corner_core_mask, YAC_LOCATION_CORNER,
        grid_ids[comp_idx]);
      yac_cset_global_index(
        comp_grid_data[comp_idx].global_cell_id, YAC_LOCATION_CELL,
        grid_ids[comp_idx]);
      yac_cset_core_mask(
        comp_grid_data[comp_idx].cell_core_mask, YAC_LOCATION_CELL,
        grid_ids[comp_idx]);

      yac_cdef_points_unstruct(
        grid_ids[comp_idx], comp_grid_data[comp_idx].nbr_cells,
        YAC_LOCATION_CELL, comp_grid_data[comp_idx].x_cells,
        comp_grid_data[comp_idx].y_cells, &(cell_point_ids[comp_idx]));
      if (comp_grid_data[comp_idx].cell_mask != NULL)
        yac_cset_mask(
          comp_grid_data[comp_idx].cell_mask, cell_point_ids[comp_idx]);
    }

    // regional grid on root processes of each component
    if ((global_rank == 0) || (global_rank == 1) || (global_rank == 2)) {

      // register regional grid on rank 0
      double regional_coord_x[] = {
        -20,-18,-16,-14,-12,-10, -8, -6, -4, -2,
          0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20};
      double regional_coord_y[] = {
        -20,-18,-16,-14,-12,-10, -8, -6, -4, -2,
          0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20};
      int num_vertices_regional[2] =
        {sizeof(regional_coord_x) / sizeof(regional_coord_x[0]),
        sizeof(regional_coord_y) / sizeof(regional_coord_y[0])};
      nbr_vertices_regional =
        num_vertices_regional[0] * num_vertices_regional[1];
      int cyclic[2] = {0, 0};
      for (int i = 0; i < num_vertices_regional[0]; ++i)
        regional_coord_x[i] *= YAC_RAD;
      for (int i = 0; i < num_vertices_regional[1]; ++i)
        regional_coord_y[i] *= YAC_RAD;
      yac_cdef_grid_reg2d(
        "regional_grid", num_vertices_regional, cyclic,
        regional_coord_x, regional_coord_y, &regional_grid_id);

      yac_cdef_points_reg2d(
        regional_grid_id, num_vertices_regional, YAC_LOCATION_CORNER,
        regional_coord_x, regional_coord_y,  &regional_point_id);
    }
  }

  // setup yac instance which couples a component to one of its own ranks
  struct intra_comp_yac_info yac_intra_infos[2];
  if (!is_dummy)
    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx)
      yac_intra_infos[comp_idx] =
        setup_intra_comp_yac(
          default_comp_comms[comp_idx], local_comp_names[comp_idx],
          cell_point_ids[comp_idx], regional_point_id);

  // register in and out fields
  int default_in_field_ids[2][2];// = {{INT_MAX, INT_MAX}, {INT_MAX, INT_MAX}};
  int default_out_field_ids[2];// = {INT_MAX, INT_MAX};
  if (!is_dummy)
    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx)
      def_default_fields(
        default_comp_ids[comp_idx], local_comp_names[comp_idx],
        grid_ids[comp_idx], cell_point_ids[comp_idx],
        default_in_field_ids[comp_idx], &(default_out_field_ids[comp_idx]));

  yac_csync_def();

  if (strcmp(yac_cget_start_datetime(), "+1800-01-01T00:00:00.000"))
    PUT_ERR("error in yac_cget_start_datetime");

  if (strcmp(yac_cget_end_datetime(), "+2100-01-01T00:00:00.000"))
    PUT_ERR("error in yac_cget_end_datetime");

  // generate interpolations for the default YAC instance
  yac_cenddef();

  // setup yac instances that only couple all component combinations
  struct comp_pair_yac_info yac_comp_pair_infos[3];
  for (int i = 0, k = 0; i < 3; ++i) {

    int local_comp_idx[2];

    local_comp_idx[0] = INT_MAX;
    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx)
      if (!strcmp(global_comp_names[i], local_comp_names[comp_idx]))
        local_comp_idx[0] = comp_idx;

    for (int j = i + 1; j < 3; ++j, ++k) {

      local_comp_idx[1] = INT_MAX;
      for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx)
        if (!strcmp(global_comp_names[j], local_comp_names[comp_idx]))
          local_comp_idx[1] = comp_idx;

      char const * comp_pair_names[2] =
        {global_comp_names[i], global_comp_names[j]};
      int nbr_cells[2] =
        {comp_grid_data[0].nbr_cells, comp_grid_data[1].nbr_cells};

      yac_comp_pair_infos[k] =
        setup_comp_pair_yac(
          comp_pair_names, cell_point_ids, nbr_cells, local_comp_idx);
    }
  }

  // if this is not the dummy process, which is not involved in any coupling
  if (!is_dummy) {

    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx) {

      int collection_size;
      char field_name[64];
      sprintf(field_name, "%s_out", local_comp_names[comp_idx]);
      const char* timestep_string = yac_cget_timestep_from_field_id(
        default_out_field_ids[comp_idx]);
      if (strcmp(timestep_string, "PT01.000S"))
        PUT_ERR("error in yac_cget_model_timestep_id");
      collection_size = yac_cget_collection_size_from_field_id(
        default_out_field_ids[comp_idx]);
      if (collection_size != 1)
        PUT_ERR("error in yac_cget_collection_size_from_field_id");
    }

    double * out_field_data[2];
    double * in_field_data[2];
    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx) {
      out_field_data[comp_idx] =
        xmalloc(comp_grid_data[comp_idx].nbr_cells * sizeof(double));
      for (int i = 0; i < comp_grid_data[comp_idx].nbr_cells; ++i)
        out_field_data[comp_idx][i] =
          (comp_grid_data[comp_idx].cell_core_mask[i])?
            (test_harmonic(
               comp_grid_data[comp_idx].x_cells[i],
               comp_grid_data[comp_idx].y_cells[i])):DUMMY_VALUE;
      in_field_data[comp_idx] =
        xmalloc(comp_grid_data[comp_idx].nbr_cells * sizeof(double));
    }

    // do some ping-pongs in the default yac instance
    for (int t = 0; t < 10; ++t) {

      for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx) {
        int info, err;
        int id = default_out_field_ids[comp_idx];
        double *point_set_data[1];
        double **collection_data[1] = {point_set_data};
        point_set_data[0] = out_field_data[comp_idx];
        yac_cput(id, 1, collection_data, &info, &err);
      }

      for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx) {
        for (int field_idx = 0; field_idx < 2; ++field_idx) {
          int info, err;
          int id = default_in_field_ids[comp_idx][field_idx];
          double *collection_data[1] = {in_field_data[comp_idx]};

          for (int i = 0; i < comp_grid_data[comp_idx].nbr_cells; ++i)
            in_field_data[comp_idx][i] = DUMMY_VALUE;
          yac_cget(id, 1, collection_data, &info, &err);
        }
      }
    }

    // interpolate received results to the regional grid defined on
    // rank 0 of the current component
    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx) {

      if (yac_intra_infos[comp_idx].is_comp_root) {
        int out_info, in_info, err;
        int out_id = yac_intra_infos[comp_idx].out_field_id;
        double *out_point_set_data[1];
        double **out_collection_data[1] = {out_point_set_data};
        out_point_set_data[0] = in_field_data[comp_idx];

        int in_id = yac_intra_infos[comp_idx].in_field_id;
        double *in_collection_data[1] = {
          xmalloc(nbr_vertices_regional * sizeof(double))};
        for (int i = 0; i < nbr_vertices_regional; ++i)
          in_collection_data[0][i] = DUMMY_VALUE;

        yac_cexchange(
          out_id, in_id, 1, out_collection_data, in_collection_data,
          &out_info, &in_info, &err);

        free(in_collection_data[0]);
      } else {
        int info, err;
        int id = yac_intra_infos[comp_idx].out_field_id;
        double *point_set_data[1];
        double **collection_data[1] = {point_set_data};
        point_set_data[0] = in_field_data[comp_idx];
        yac_cput(id, 1, collection_data, &info, &err);
      }
    }

    // clean-up

    for (size_t comp_idx = 0; comp_idx < num_local_comps; ++comp_idx) {
      free(out_field_data[comp_idx]);
      free(in_field_data[comp_idx]);
      grid_data_free(comp_grid_data[comp_idx]);
    }

    // do some exchanges between within the component pair YAC instances
    // for all available component pairs
    for (int i = 0; i < 3; ++i) {

      // if the local process is not a part of the current component pair
      if (yac_comp_pair_infos[i].yac_id == INT_MAX) continue;

      // if the local process has both components of the current pair
      if ((yac_comp_pair_infos[i].nbr_cells[0] > 0) &&
          (yac_comp_pair_infos[i].nbr_cells[1] > 0)) {

        double * out_field_data[2];
        double * in_field_data[2];

        for (int j = 0; j < 2; ++j) {
          int nbr_cells = yac_comp_pair_infos[i].nbr_cells[j];
          out_field_data[j] = xmalloc(nbr_cells * sizeof(double));
          for (int k = 0; k < nbr_cells; ++k) out_field_data[j][k] = 0.0;
          in_field_data[j] = xmalloc(nbr_cells * sizeof(double));
          for (int k = 0; k < nbr_cells; ++k) in_field_data[j][k] = DUMMY_VALUE;
        }

        for (int j = 0; j < 2; ++j) {

          int out_info, in_info, err;
          int out_id = yac_comp_pair_infos[i].out_field_ids[j];
          double *out_point_set_data[1];
          double **out_collection_data[1] = {out_point_set_data};
          out_point_set_data[0] = out_field_data[j];

          int in_id = yac_comp_pair_infos[i].in_field_ids[j^1];
          double *in_collection_data[1] = {in_field_data[j^1]};

          yac_cexchange(
            out_id, in_id, 1, out_collection_data, in_collection_data,
            &out_info, &in_info, &err);
        }
        for (int j = 0; j < 2; ++j) {
          free(out_field_data[j]);
          free(in_field_data[j]);
        }

      } else {

        for (int j = 0; j < 2; ++j) {

          // if the source component is available locally
          if (yac_comp_pair_infos[i].nbr_cells[j] > 0) {

            int nbr_cells = yac_comp_pair_infos[i].nbr_cells[j];

            double * out_field_data = xmalloc(nbr_cells * sizeof(double));
            // just some dummy data
            for (int j = 0; j < nbr_cells; ++j) out_field_data[j] = 0.0;

            int info, err;
            int id = yac_comp_pair_infos[i].out_field_ids[j];
            double *point_set_data[1];
            double **collection_data[1] = {point_set_data};
            point_set_data[0] = out_field_data;
            yac_cput(id, 1, collection_data, &info, &err);

            free(out_field_data);

          } else {

            int nbr_cells = yac_comp_pair_infos[i].nbr_cells[j^1];

            double * in_field_data = xmalloc(nbr_cells * sizeof(double));

            int info, err;
            int id = yac_comp_pair_infos[i].in_field_ids[j^1];
            double *collection_data[1] = {in_field_data};

            for (int j = 0; j < nbr_cells; ++j) in_field_data[j] = DUMMY_VALUE;
            yac_cget(id, 1, collection_data, &info, &err);

            free(in_field_data);
          }
        }
      }
    }
  }

  for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx)
    MPI_Comm_free(&(default_comp_comms[comp_idx]));

  if (!is_dummy) {
    for (int i = 0; i < 3; ++i)
      if (yac_comp_pair_infos[i].yac_id != INT_MAX)
        yac_ccleanup_instance(yac_comp_pair_infos[i].yac_id);
    for (int comp_idx = 0; comp_idx < num_local_comps; ++comp_idx)
      yac_ccleanup_instance(yac_intra_infos[comp_idx].yac_id);
  }
  yac_cfinalize();

  xt_finalize();
  MPI_Finalize();

#else

  PUT_ERR("YAC configured without YAML or NetCDF\n");

#endif

  return TEST_EXIT_CODE;
}

static struct intra_comp_yac_info setup_intra_comp_yac(
  MPI_Comm comp_comm, char const * comp_name,
  int base_point_id, int regional_point_id) {

  // initialise intra component YAC instance
  int yac_id;
  char yaml_file_name[128];
  sprintf(yaml_file_name, "coupling_test3_%s_intra.yaml", comp_name);
  yac_cinit_comm_instance(comp_comm, &yac_id);
  yac_cread_config_yaml_instance(yac_id, yaml_file_name);

  // register component(s) (all processes are part of the base component
  //                        process 0 has an additional component)
  char const * comp_names[2] = {comp_name, "comp_regional"};
  int num_comps = (regional_point_id != -1)?2:1;
  int comp_ids[2];
  if (regional_point_id != -1)
    yac_cdef_comps_instance(yac_id, comp_names, num_comps, comp_ids);
  else
    yac_cdef_comp_instance(yac_id, comp_names[0], &comp_ids[0]);

  // register field (from base component to regional one)
  int field_ids[2];
  yac_cdef_field(
    "base_to_regional", comp_ids[0], &base_point_id, 1, 1, "1",
    YAC_TIME_UNIT_SECOND, &(field_ids[0]));
  if (regional_point_id != -1)
    yac_cdef_field(
      "base_to_regional", comp_ids[1], &regional_point_id, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &(field_ids[1]));

  yac_csync_def_instance(yac_id);

  if (strcmp(yac_cget_start_datetime_instance(yac_id),
            "+1800-01-01T00:00:00.000"))
    PUT_ERR("error in yac_cget_start_datetime_instance");

  if (strcmp(yac_cget_end_datetime_instance(yac_id),
            "+2100-01-01T00:00:00.000"))
    PUT_ERR("error in yac_cget_end_datetime_instance");

  // generate interpolations for the default YAC instance
  yac_cenddef_instance(yac_id);

  struct intra_comp_yac_info yac_intra_info;
  yac_intra_info.yac_id = yac_id;
  yac_intra_info.is_comp_root = regional_point_id != -1;
  yac_intra_info.out_field_id = field_ids[0];
  yac_intra_info.in_field_id = field_ids[1];

  int collection_size;
  const char* timestep_string = yac_cget_timestep_from_field_id(field_ids[0]);
  if (strcmp(timestep_string, "PT01.000S"))
    PUT_ERR("error in yac_cget_model_timestep_id");
  collection_size = yac_cget_collection_size_from_field_id(field_ids[0]);
  if (collection_size != 1)
    PUT_ERR("error in yac_cget_collection_size_from_field_id");
  return yac_intra_info;
}

static struct comp_pair_yac_info setup_comp_pair_yac(
  char const * comp_pair_names[2], int point_ids[2], int nbr_cells[2],
  int local_comp_indices[2]) {

  if ((local_comp_indices[0] == INT_MAX) &&
      (local_comp_indices[1] == INT_MAX))
    return (struct comp_pair_yac_info){.yac_id = INT_MAX};

  // get the communicator for the component pair
  MPI_Comm comp_pair_comm;
  yac_cget_comps_comm(
    (char const *[]){comp_pair_names[0], comp_pair_names[1]}, 2,
    &comp_pair_comm);

  char yaml_file_name[256];
  sprintf(yaml_file_name, "coupling_test3_%s_%s.yaml",
         comp_pair_names[0], comp_pair_names[1]);

  // initialise component pair YAC instance
  int yac_id;
  yac_cinit_comm_instance(comp_pair_comm, &yac_id);
  yac_cread_config_yaml_instance(yac_id, yaml_file_name);

  // register component(s) the local process can have one or two components
  // associated to the YAC instance
  char const * comp_names[2];
  int num_comps = 0;
  int comp_ids[2];
  for (int i = 0; i < 2; ++i)
    if (local_comp_indices[i] != INT_MAX)
      comp_names[num_comps++] = comp_pair_names[i];
  yac_cdef_comps_instance(yac_id, comp_names, num_comps, comp_ids);

  int field_ids[2][2];
  for (int i = 0, k = 0; i < 2; ++i) {

    int comp_idx = local_comp_indices[i];
    if (comp_idx == INT_MAX) continue;

    // register field
    for (int j = 0; j < 2; ++j) {
      char field_name[16];
      sprintf(field_name, "comp_%c_out", comp_pair_names[j][5]);
      yac_cdef_field(
        field_name, comp_ids[k], &(point_ids[comp_idx]), 1, 1, "1",
        YAC_TIME_UNIT_SECOND, &(field_ids[i][j]));
    }
    ++k;
  }

  // generate interpolations for the component pair YAC instance
  yac_cenddef_instance(yac_id);

  MPI_Comm_free(&comp_pair_comm);

  struct comp_pair_yac_info yac_info;

  yac_info.yac_id = yac_id;
  for (int i = 0; i < 2; ++i) {
    int comp_idx = local_comp_indices[i];
    if (comp_idx == INT_MAX) {
      yac_info.nbr_cells[i] = 0;
      yac_info.out_field_ids[i] = INT_MAX;
      yac_info.in_field_ids[i] = INT_MAX;
    } else {
      yac_info.nbr_cells[i] = nbr_cells[comp_idx];
      int roles[2];
      for (int j = 0; j < 2; ++j) {
        roles[j] = yac_cget_role_from_field_id(field_ids[i][j]);
        switch (roles[j]) {
          case(0):
          default:
            PUT_ERR("invalid role");
            break;
          case(1):
            yac_info.out_field_ids[i] = field_ids[i][j];
            break;
          case(2):
            yac_info.in_field_ids[i] = field_ids[i][j];
            break;
        }
      }
      if (roles[0] == roles[1]) PUT_ERR("roles should not match");
    }
  }

  return yac_info;
}

static void def_default_fields(
  int comp_id, char const * comp_name, int grid_id, int cell_point_id,
  int * in_field_ids, int * out_field_id) {

  int num_in_fields = 0;

  for (int i = 0; i < 3; ++i) {
    char const * field_names[3] = {"comp_a_out", "comp_b_out", "comp_c_out"};
    int field_id;
    char const * curr_field_name = field_names[i];
    yac_cdef_field(
      curr_field_name, comp_id, &cell_point_id, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_id);

    if (comp_name[5] == curr_field_name[5])
      *out_field_id = field_id;
    else
      in_field_ids[num_in_fields++] = field_id;
  }
}

static struct grid_data generate_icon_grid_data(MPI_Comm comm) {

  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  struct grid_data grid_data;

  static char const * icon_grid_name = "icon_grid";
  grid_data.name = icon_grid_name;

  read_part_icon_grid_information(
    "icon_grid_R02B02.nc", &grid_data.nbr_vertices, &grid_data.nbr_cells,
    &grid_data.num_vertices_per_cell, &grid_data.cell_to_vertex,
    &grid_data.x_vertices, &grid_data.y_vertices,
    &grid_data.x_cells, &grid_data.y_cells,
    &grid_data.global_cell_id, &grid_data.cell_mask,
    &grid_data.cell_core_mask, &grid_data.global_corner_id,
    &grid_data.corner_core_mask, comm_rank, comm_size);

  for (size_t i = 0; i < grid_data.nbr_cells; ++i)
    grid_data.cell_mask[i] = grid_data.cell_mask[i] == 0;

  return grid_data;
}

static struct grid_data generate_cube_grid_data(MPI_Comm comm) {

  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  struct grid_data grid_data;

  static char const * cube_grid_name = "cube_grid";
  grid_data.name = cube_grid_name;

  unsigned n = 50;

  generate_part_cube_grid_information(
    n, (unsigned*)&grid_data.nbr_vertices, (unsigned*)&grid_data.nbr_cells,
    (unsigned**)&grid_data.num_vertices_per_cell,
    (unsigned**)&grid_data.cell_to_vertex,
    &grid_data.x_vertices, &grid_data.y_vertices,
    &grid_data.x_cells, &grid_data.y_cells, &grid_data.global_cell_id,
    &grid_data.cell_core_mask, &grid_data.global_corner_id,
    &grid_data.corner_core_mask, comm_rank, comm_size);

  grid_data.cell_mask = NULL;

  return grid_data;
}

static struct grid_data generate_mpiom_grid_data(MPI_Comm comm) {

  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  struct grid_data grid_data;

  static char const * mpiom_grid_name = "mpiom_grid";
  grid_data.name = mpiom_grid_name;

  read_part_mpiom_grid_information(
    "GR30_lsm.nc", &grid_data.nbr_vertices, &grid_data.nbr_cells,
    &grid_data.num_vertices_per_cell, &grid_data.cell_to_vertex,
    &grid_data.x_vertices, &grid_data.y_vertices,
    &grid_data.x_cells, &grid_data.y_cells,
    &grid_data.global_cell_id, &grid_data.cell_mask,
    &grid_data.cell_core_mask, &grid_data.global_corner_id,
    &grid_data.corner_core_mask, comm_rank, comm_size);

  for (size_t i = 0; i < grid_data.nbr_cells; ++i)
    grid_data.cell_mask[i] = grid_data.cell_mask[i] == 0;

  return grid_data;
}

static void grid_data_free(struct grid_data grid_data) {

  free(grid_data.num_vertices_per_cell);
  free(grid_data.cell_to_vertex);
  free(grid_data.x_vertices);
  free(grid_data.y_vertices);
  free(grid_data.x_cells);
  free(grid_data.y_cells);

  free(grid_data.cell_mask);
  free(grid_data.global_cell_id);
  free(grid_data.cell_core_mask);
  free(grid_data.global_corner_id);
  free(grid_data.corner_core_mask);
}
