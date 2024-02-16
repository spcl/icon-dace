/**
 * @file test_interpolation_parallel4.c
 *
 * @copyright Copyright  (C)  2020 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <mpi.h>

#include "tests.h"

#include "tests.h"
#include "interp_method.h"
#include "interp_method_file.h"
#include "interp_method_avg.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"
#include "interpolation.h"
#include "test_common.h"
#include "weight_file_common.h"

// needs to be a multiple of 2
#define NUM_CELLS_X 128
#define NUM_CELLS_Y 128
#define NUM_CELLS_X_2 64
#define NUM_CELLS_Y_2 64

static void generate_input_weights();

static void generate_ref_weights();

static void submain_src(MPI_Comm global_comm, MPI_Comm comp_comm,
                        enum interp_weights_reorder_type reorder_type);

static void submain_tgt(MPI_Comm global_comm, MPI_Comm comp_comm,
                        enum interp_weights_reorder_type reorder_type);

// link data for input weight file
unsigned num_links_file = 4*NUM_CELLS_X_2*NUM_CELLS_Y_2;
int src_address_file[NUM_CELLS_Y_2][NUM_CELLS_X_2][4];
int tgt_address_file[NUM_CELLS_Y_2][NUM_CELLS_X_2][4];
double weights_file[NUM_CELLS_Y_2][NUM_CELLS_X_2][4];

char const weight_file_in[] = "weight_file_in.nc";
char const weight_file_out[] = "weight_file_out.nc";
char const src_grid_name[] = "src_grid";
char const tgt_grid_name[] = "tgt_grid";

// reference link data for output weight file
unsigned ref_num_links = 4*NUM_CELLS_X*NUM_CELLS_Y;

int ref_src_address[NUM_CELLS_Y][NUM_CELLS_X][4];
int ref_tgt_address[NUM_CELLS_Y][NUM_CELLS_X][4];
double ref_weights[NUM_CELLS_Y][NUM_CELLS_X][4];

int const * ref_tgt_address_fixed = NULL;
double const * ref_fixed_values = NULL;
int const * ref_num_tgt_per_fixed_value = NULL;
unsigned ref_num_fixed_values = 0;

// grid information (is the same for source and target)
unsigned cyclic[2] = {0,0};
double global_coordinates_x[NUM_CELLS_X + 1];
double global_coordinates_y[NUM_CELLS_X + 1];
int with_halo = 1;

int main(int argc, char *argv[]) {

  if (argc != 2) {
    PUT_ERR("wrong number of arguments\n");
    return TEST_EXIT_CODE;
  }

  enum interp_weights_reorder_type reorder_type =
   (strcmp(argv[1], "src") == 0)?MAPPING_ON_SRC:MAPPING_ON_TGT;

  if ((reorder_type != MAPPING_ON_SRC) && strcmp(argv[1], "tgt")) {
    PUT_ERR("invalid argument (has to be either \"src\" or \"tgt\")\n");
    return TEST_EXIT_CODE;
  }

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  generate_input_weights();
  generate_ref_weights();

  for (unsigned i = 0; i < NUM_CELLS_X + 1; ++i)
    global_coordinates_x[i] = (double)i * YAC_RAD;

  for (unsigned i = 0; i < NUM_CELLS_Y + 1; ++i)
    global_coordinates_y[i] = (double)i * YAC_RAD;

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 5) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  // split processes into source an target

  int comp_flag = comm_rank < 4;

  MPI_Comm split_comm;
  yac_mpi_call(
    MPI_Comm_split(
      MPI_COMM_WORLD, comp_flag, 0, &split_comm), MPI_COMM_WORLD);

  if (comp_flag) submain_src(MPI_COMM_WORLD, split_comm, reorder_type);
  else           submain_tgt(MPI_COMM_WORLD, split_comm, reorder_type);

  yac_mpi_call(MPI_Comm_free(&split_comm), MPI_COMM_WORLD);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void submain_tgt(MPI_Comm global_comm, MPI_Comm comp_comm,
                        enum interp_weights_reorder_type reorder_type) {

  size_t local_start[2] = {0,0};
  size_t local_count[2] = {NUM_CELLS_X, NUM_CELLS_Y};

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      global_coordinates_x, global_coordinates_y, local_count,
      local_start, local_count, with_halo);

  double cell_coordinates_x[NUM_CELLS_X];
  double cell_coordinates_y[NUM_CELLS_Y];

  for (unsigned i = 0; i < NUM_CELLS_X; ++i)
    cell_coordinates_x[i] = ((double)i + 0.5)*YAC_RAD;
  for (unsigned i = 0; i < NUM_CELLS_Y; ++i)
    cell_coordinates_y[i] = ((double)i + 0.5)*YAC_RAD;

  double cell_field_coords[NUM_CELLS_Y][NUM_CELLS_X][3];
  for (size_t i = 0; i < NUM_CELLS_Y; ++i)
    for (size_t j = 0; j < NUM_CELLS_X; ++j)
      LLtoXYZ(
        cell_coordinates_x[j], cell_coordinates_y[i], cell_field_coords[i][j]);

  struct yac_basic_grid * tgt_grid =
    yac_basic_grid_new(tgt_grid_name, grid_data);
  yac_basic_grid_add_coordinates(
    tgt_grid, CELL, &(cell_field_coords[0][0]), grid_data.num_cells);
  struct yac_basic_grid * src_grid =
    yac_basic_grid_empty_new(src_grid_name);

  struct dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

  { // test interpolation

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_file_new(
         weight_file_in, src_grid_name, tgt_grid_name),
       yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       NULL};

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    struct interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

    //---------------------
    // do the interpolation
    //---------------------

    // target_data dimensions [collection_idx]
    //                        [local_idx]

    double target_field[NUM_CELLS_X][NUM_CELLS_Y];
    double * target_field_ = &target_field[0][0];
    double ** target_data = &target_field_;

    for (unsigned i = 0; i < NUM_CELLS_Y; ++i)
      for (unsigned j = 0; j < NUM_CELLS_X; ++j)
        target_field[i][j] = -1;

    yac_interpolation_execute_get(interpolation, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    yac_interp_weights_write_to_file(
      weights, weight_file_out, src_grid_name, tgt_grid_name);

    double ref_target_field[NUM_CELLS_X*NUM_CELLS_Y];

    for (unsigned i = 0; i < NUM_CELLS_X*NUM_CELLS_Y; ++i)
      ref_target_field[i] = 0;

    for (unsigned i = 0; i < NUM_CELLS_Y; ++i)
      for (unsigned j = 0; j < NUM_CELLS_X; ++j)
        for (unsigned k = 0; k < 4; ++k)
          ref_target_field[ref_tgt_address[i][j][k]] +=
            ((double)ref_src_address[i][j][k]) * ref_weights[i][j][k];

    for (unsigned i = 0; i < NUM_CELLS_Y; ++i)
      for (unsigned j = 0; j < NUM_CELLS_X; ++j)
        if (fabs(target_field[i][j] - ref_target_field[i * NUM_CELLS_X + j]) >
            1e-10)
          PUT_ERR("wrong interpolation result\n")

    //--------
    // cleanup
    //--------

    yac_interpolation_delete(interpolation);
    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
  }

  //--------
  // cleanup
  //--------

  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
}

static void submain_src(MPI_Comm global_comm, MPI_Comm comp_comm,
                        enum interp_weights_reorder_type reorder_type) {

  int my_rank;
  yac_mpi_call(MPI_Comm_rank(comp_comm, &my_rank), comp_comm);

   // create a weight file
  if (my_rank == 0) {
    enum yac_location src_locations[1] = {CORNER};
    enum yac_location tgt_location = CELL;
    int * tgt_id_fixed = NULL;
    unsigned num_fixed_tgt = 0;
    double * fixed_values = NULL;
    int * num_tgt_per_fixed_value = NULL;
    unsigned num_fixed_values = 0;
    write_weight_file(weight_file_in, &(src_address_file[0][0][0]),
                      &(tgt_address_file[0][0][0]), &(weights_file[0][0][0]),
                      num_links_file, src_locations, 1, (int*)&num_links_file,
                      tgt_id_fixed, num_fixed_tgt, fixed_values,
                      num_tgt_per_fixed_value, num_fixed_values,
                      tgt_location, src_grid_name, tgt_grid_name);
  }

  size_t local_start[4][2] =
    {{0,0},{NUM_CELLS_X_2,0},{NUM_CELLS_X_2, NUM_CELLS_Y_2},{0,NUM_CELLS_Y_2}};
  size_t local_count[2] = {NUM_CELLS_X_2,NUM_CELLS_Y_2};
  size_t global_num_cells[2] = {NUM_CELLS_X,NUM_CELLS_Y};

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      global_coordinates_x, global_coordinates_y, global_num_cells,
      local_start[my_rank], local_count, with_halo);
  size_t num_vertices = grid_data.num_vertices;

  struct yac_basic_grid * src_grid =
    yac_basic_grid_new(src_grid_name, grid_data);
  yac_basic_grid_add_coordinates_nocpy(src_grid, CELL, NULL);
  struct yac_basic_grid * tgt_grid =
    yac_basic_grid_empty_new(tgt_grid_name);

  struct dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

  { // test interpolation

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_file_new(
         weight_file_in, src_grid_name, tgt_grid_name),
       yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       NULL};

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    struct interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

    //---------------------
    // do the interpolation
    //---------------------

    // source_data dimensions [collection_idx]
    //                        [pointset_idx]
    //                        [local_idx]
    double * source_data_field =
      xmalloc(num_vertices * sizeof(*source_data_field));

    double * source_data_pointset[1] = {source_data_field};  // num_pointset == 1
    double ** source_data[1] = {source_data_pointset}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i)
      source_data_field[i] =
        (grid_data.core_vertex_mask[i])?
          ((double)(grid_data.vertex_ids[i])):(-1.0);

    yac_interpolation_execute_put(interpolation, source_data);

    //------------------------------
    // check the written weight file
    //------------------------------

    yac_interp_weights_write_to_file(
      weights, weight_file_out, src_grid_name, tgt_grid_name);

    if (my_rank == 0) {
      enum yac_location ref_src_locations[1] = {CORNER};
      enum yac_location ref_tgt_location = CELL;
      check_weight_file(weight_file_out, &ref_src_address[0][0][0],
                        &ref_tgt_address[0][0][0], &ref_weights[0][0][0],
                        ref_num_links, ref_src_locations, 1, (int*)&ref_num_links,
                        ref_tgt_address_fixed, ref_fixed_values,
                        ref_num_tgt_per_fixed_value, ref_num_fixed_values,
                        ref_tgt_location, src_grid_name, tgt_grid_name);
    }

    //--------
    // cleanup
    //--------

    free(source_data_field);
    yac_interpolation_delete(interpolation);
    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
  }

  //--------
  // cleanup
  //--------

  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(tgt_grid);
  yac_basic_grid_delete(src_grid);

  // delete weight file
  if (my_rank == 0) {
    unlink(weight_file_in);
    unlink(weight_file_out);
  }
}

static void
generate_input_weights() {

  for (unsigned i = 0; i < NUM_CELLS_Y_2; ++i) {
    for (unsigned j = 0; j < NUM_CELLS_X_2; ++j) {
      src_address_file[i][j][0] =
        0 + (i + NUM_CELLS_Y_2 / 2) * (NUM_CELLS_X + 1) + (j + NUM_CELLS_X_2 / 2);
      src_address_file[i][j][1] =
        1 + (i + NUM_CELLS_Y_2 / 2) * (NUM_CELLS_X + 1) + (j + NUM_CELLS_X_2 / 2);
      src_address_file[i][j][2] =
        0 + (i + 1 + NUM_CELLS_Y_2 / 2) * (NUM_CELLS_X + 1) + (j + NUM_CELLS_X_2 / 2);
      src_address_file[i][j][3] =
        1 + (i + 1 + NUM_CELLS_Y_2 / 2) * (NUM_CELLS_X + 1) + (j + NUM_CELLS_X_2 / 2);
      for (unsigned k = 0; k < 4; ++k)
        tgt_address_file[i][j][k] =
          (i + NUM_CELLS_Y_2 / 2) * NUM_CELLS_X + (j + NUM_CELLS_X_2 / 2);
      for (unsigned k = 0; k < 4; ++k)
        weights_file[i][j][k] = (double)(k + 1) * 0.1;
    }
  }
}

static void
generate_ref_weights() {

  for (unsigned i = 0, idx = 0; i < NUM_CELLS_Y; ++i) {
    for (unsigned j = 0; j < NUM_CELLS_X; ++j, ++idx) {
      ref_src_address[i][j][0] = 0 + i * (NUM_CELLS_X + 1) + j;
      ref_src_address[i][j][1] = 1 + i * (NUM_CELLS_X + 1) + j;
      ref_src_address[i][j][2] = 0 + (i + 1) * (NUM_CELLS_X + 1) + j;
      ref_src_address[i][j][3] = 1 + (i + 1) * (NUM_CELLS_X + 1) + j;
      for (unsigned k = 0; k < 4; ++k) {
        ref_tgt_address[i][j][k] = idx;
        ref_weights[i][j][k] = 0.25;
      }
    }
  }
  for (unsigned i = 0; i < NUM_CELLS_Y_2; ++i)
    for (unsigned j = 0; j < NUM_CELLS_X_2; ++j)
      for (unsigned k = 0; k < 4; ++k)
        ref_weights
          [tgt_address_file[i][j][k]/NUM_CELLS_X]
          [tgt_address_file[i][j][k]%NUM_CELLS_X][k] = weights_file[i][j][k];
}

