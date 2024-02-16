/**
 * @file OASIS3_MCT_example.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
 * This example is based on a example presented in the following paper:
 * "Development and performance of a new version of the OASIS coupler, OASIS3-MCT_3.0"
 * (http://www.geosci-model-dev-discuss.net/gmd-2017-64/)
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
#include <unistd.h>
#include <string.h>

#include "yac_interface.h"
#include "read_icon_grid.h"
#include "event.h"
#include "test_function.h"
#include "utils.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

struct {

  unsigned component_idx[3];
  unsigned num_components;

  unsigned grid_idx[3];

} config_per_global_rank[38];

char const * component_names[] = {"comp1-grid1", "comp2-grid2", "comp3",
                                  "comp3-grid3", "comp3-grid4", "comp3-grid5",
                                  "comp4"};
char const * grid_names[] = {"grid1", "grid2", "grid3", "grid4", "grid5"};

struct {

  char const * field_names[8];
  unsigned num_fields;
  char const * file_name;

} grid_config[] = {{.field_names = {"A_1", "D_1", "E_1", "J_1",
                                    "A_2", "D_2", "E_2", "J_2"},
                    .num_fields = 8,
                    .file_name = "iconR2B05-grid.nc"},
                   {.field_names = {"A_1", "B_1", "F_1", "I_1",
                                    "A_2", "B_2", "F_2", "I_2"},
                    .num_fields = 8,
                    .file_name = "iconR2B05-grid.nc"},
                   {.field_names = {"B_1", "C_1", "D_1", "G_1",
                                    "B_2", "C_2", "D_2", "G_2"},
                    .num_fields = 8,
                    .file_name = "iconR2B06-grid.nc"},
                   {.field_names = {"C_1", "E_1", "F_1", "H_1",
                                    "C_2", "E_2", "F_2", "H_2"},
                    .num_fields = 8,
                    .file_name = "iconR2B05-grid.nc"},
                   {.field_names = {"G_1", "H_1", "I_1", "J_1",
                                    "G_2", "H_2", "I_2", "J_2"},
                    .num_fields = 8,
                    .file_name = "iconR2B07-grid.nc"}};

/* -------------------------------------------------------------------- */

static void setup_process_config();

/* -------------------------------------------------------------------- */

// redefine YAC assert macros
#undef YAC_ASSERT
#define STR_USAGE "Usage: %s -c configFilename -x num_cells_x -y num_cells_y\n"
#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, char const ** configFilename);

int main (int argc, char *argv[]) {
#ifdef YAC_NETCDF_ENABLED

  setup_process_config();

  yac_cinit ();

  char const * configFilename = "toy_OASIS3_MCT.yaml"; // default configuration file
  parse_arguments(argc, argv, &configFilename);
  yac_cread_config_yaml(configFilename);
  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);

  if (global_size != 38) {
    fprintf(stderr, "Wrong number of processes (should be 38)\n");
    exit(EXIT_FAILURE);
  }

  // register component

  char const * local_component_names[3];
  int component_ids[3];
  unsigned num_components = config_per_global_rank[global_rank].num_components;
  for (unsigned i = 0; i < num_components; ++i)
    local_component_names[i] =
      component_names[config_per_global_rank[global_rank].component_idx[i]];
  yac_cdef_comps(local_component_names, num_components, component_ids);

  int component_ranks[3], component_sizes[3];
  for (unsigned i = 0; i < num_components; ++i) {
    MPI_Comm component_comm;
    yac_cget_comp_comm(component_ids[i], &component_comm);
    MPI_Comm_rank(component_comm, &component_ranks[i]);
    MPI_Comm_size(component_comm, &component_sizes[i]);
  }

  // register grid and fields

  int field_ids[3*8];
  unsigned num_fields_per_grid[3];
  unsigned num_fields = 0;

  double * field_data[3][8];

  for (unsigned i = 0; i < num_components; ++i) {

    if (config_per_global_rank[global_rank].grid_idx[i] == -1) {
      num_fields_per_grid[i] = 0;
      continue;
    }

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

    read_part_icon_grid_information(
      grid_config[config_per_global_rank[global_rank].grid_idx[i]].file_name,
      &nbr_vertices, &nbr_cells, &num_vertices_per_cell, &cell_to_vertex,
      &x_vertices, &y_vertices, &x_cells, &y_cells, &global_cell_id, &cell_mask,
      &cell_core_mask, &global_corner_id, &corner_core_mask,
      component_ranks[i], component_sizes[i]);

    int grid_id;
    yac_cdef_grid_unstruct(
      grid_names[config_per_global_rank[global_rank].grid_idx[i]],
      nbr_vertices, nbr_cells, num_vertices_per_cell,
      x_vertices, y_vertices, cell_to_vertex, &grid_id);

    yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
    yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);
    yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
    yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);

    int corner_point_id, cell_point_id;
    yac_cdef_points_unstruct(
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices,
      &corner_point_id);
    yac_cdef_points_unstruct(
      grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells,
      &cell_point_id);

    yac_cset_mask(cell_mask, cell_point_id);

    num_fields_per_grid[i] =
      grid_config[config_per_global_rank[global_rank].grid_idx[i]].num_fields;

    for (unsigned j = 0; j < num_fields_per_grid[i]; ++j) {

      yac_cdef_field(
        grid_config[
          config_per_global_rank[global_rank].grid_idx[i]].field_names[j],
          component_ids[i], &cell_point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND,
          &field_ids[num_fields]);

      field_data[i][j] = xmalloc(nbr_cells * sizeof(*(field_data[i][j])));
      for (int k = 0; k < nbr_cells; ++k)
        field_data[i][j][k] = test_func_deg(x_cells[k], y_cells[k]);

      num_fields++;
    }

    delete_icon_grid_data ( &cell_mask,
                            &global_cell_id,
                            &cell_core_mask,
                            &num_vertices_per_cell,
                            &global_corner_id,
                            &corner_core_mask,
                            &cell_to_vertex,
                            &x_cells,
                            &y_cells,
                            &x_vertices,
                            &y_vertices );
  }

  yac_cenddef( );

  // do some ping-pongs

  for (unsigned t = 0; t < 100; ++t) {

    // call put for all local fields
    {
      double *point_set_data[1];
      double **collection_data[1] = {point_set_data};

      for (unsigned i = 0, k = 0; i < num_components; ++i) {
        for (unsigned j = 0; j < num_fields_per_grid[i]; ++j, ++k) {

          int err, info;
          point_set_data[0] = field_data[i][j];
          yac_cput(field_ids[k], 1, collection_data, &info, &err);
        }
      }
    }

    // call get for all local fields
    {
      for (unsigned i = 0, k = 0; i < num_components; ++i) {
        for (unsigned j = 0; j < num_fields_per_grid[i]; ++j, ++k) {

          int err, info;
          yac_cget(field_ids[k], 1, &(field_data[i][j]), &info, &err);
        }
      }
    }
  }

  // clean-up

  for (unsigned i = 0; i < num_components; ++i)
    for (unsigned j = 0; j < num_fields_per_grid[i]; ++j)
      free(field_data[i][j]);

  yac_cfinalize();

#else

  fprintf(stderr, "YAC configured without NetCDF\n");
  exit(EXIT_FAILURE);

#endif
}

/* -------------------------------------------------------------------- */

static void setup_process_config() {

  unsigned rank = 0;

  for (; rank < 6; ++rank) {
    config_per_global_rank[rank].component_idx[0] = 0;
    config_per_global_rank[rank].num_components = 1;
    config_per_global_rank[rank].grid_idx[0] = 0;
  }
  for (; rank < 12; ++rank) {
    config_per_global_rank[rank].component_idx[0] = 1;
    config_per_global_rank[rank].num_components = 1;
    config_per_global_rank[rank].grid_idx[0] = 1;
  }
  for (; rank < 22; ++rank) {
    config_per_global_rank[rank].component_idx[0] = 2;
    config_per_global_rank[rank].component_idx[1] = 3;
    config_per_global_rank[rank].component_idx[2] = 5;
    config_per_global_rank[rank].num_components = 3;
    config_per_global_rank[rank].grid_idx[0] = -1;
    config_per_global_rank[rank].grid_idx[1] = 2;
    config_per_global_rank[rank].grid_idx[2] = 4;
  }
  for (; rank < 27; ++rank) {
    config_per_global_rank[rank].component_idx[0] = 2;
    config_per_global_rank[rank].component_idx[1] = 4;
    config_per_global_rank[rank].component_idx[2] = 5;
    config_per_global_rank[rank].num_components = 3;
    config_per_global_rank[rank].grid_idx[0] = -1;
    config_per_global_rank[rank].grid_idx[1] = 3;
    config_per_global_rank[rank].grid_idx[2] = 4;
  }
  for (; rank < 31; ++rank) {
    config_per_global_rank[rank].component_idx[0] = 2;
    config_per_global_rank[rank].component_idx[1] = 4;
    config_per_global_rank[rank].num_components = 2;
    config_per_global_rank[rank].grid_idx[0] = -1;
    config_per_global_rank[rank].grid_idx[1] = 3;
  }
  for (; rank < 34; ++rank) {
    config_per_global_rank[rank].component_idx[0] = 2;
    config_per_global_rank[rank].num_components = 1;
    config_per_global_rank[rank].grid_idx[0] = -1;
  }
  for (; rank < 38; ++rank) {
    config_per_global_rank[rank].component_idx[0] = 6;
    config_per_global_rank[rank].num_components = 1;
    config_per_global_rank[rank].grid_idx[0] = -1;
  }
}

static void parse_arguments(
  int argc, char ** argv, char const ** configFilename) {

  int opt;
  while ((opt = getopt(argc, argv, "c:")) != -1) {
    YAC_ASSERT((opt == 'c'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
    }
  }
}
