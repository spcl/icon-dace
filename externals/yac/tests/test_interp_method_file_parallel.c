/**
 * @file test_interp_method_file_parallel.c
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
#include <unistd.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "interp_method.h"
#include "interp_method_file.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"
#include "weight_file_common.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

double const tol = 1e-7;
double const err_tol = 1e-14;

char const src_grid_name[] = "src_grid";
char const tgt_grid_name[] = "tgt_grid";
char const file_name[] = "weights.nc";
char const file_name_2[] = "weights_2.nc";

static void target_main(MPI_Comm global_comm, MPI_Comm target_comm);
static void source_main(MPI_Comm global_comm, MPI_Comm source_comm);

int main (void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 3) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  int tgt_flag = comm_rank < 1;

  MPI_Comm split_comm;
  yac_mpi_call(
    MPI_Comm_split(
      MPI_COMM_WORLD, tgt_flag, 0, &split_comm), MPI_COMM_WORLD);

  if (tgt_flag) target_main(MPI_COMM_WORLD, split_comm);
  else          source_main(MPI_COMM_WORLD, split_comm);

  yac_mpi_call(MPI_Comm_free(&split_comm), MPI_COMM_WORLD);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void source_main(MPI_Comm global_comm,
                        MPI_Comm source_comm) {

  int my_source_rank;
  yac_mpi_call(MPI_Comm_rank(source_comm, &my_source_rank), source_comm);

 { // simple weight file interpolation (one link per target point) with two
   // source processes and a single target process the source fields have two
   // struct points per grid (the test check how the interpolation
   // methods handle NULL target coordinates )

    // the global source grid has 3x2 cells:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    //---------------
    // setup
    //---------------

    // create weight file
    if (my_source_rank == 0) {

      int src_indices[] = {0,1,2,3,4,5,6,7,8,9,10,11};
      int tgt_indices[] = {0,1,2,3,4,5,6,7,8,9,10,11};
      double weights[] = {0,1,2,3,4,5,6,7,8,9,10,11};
      size_t num_links = 12;
      enum yac_location src_locations[2] = {CORNER, CELL};
      enum yac_location tgt_location = CORNER;
      unsigned num_src_fields = 2;
      int num_links_per_field[2] = {num_links, 0};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0};
    size_t const num_global_cells[2] = {3,2};
    size_t local_start[2][2] = {{0,0},{2,0}};
    size_t local_count[2][2] = {{2,2},{1,2}};
    int with_halo = 1;
    int global_corner_mask[3][4] = {
      {1,1,1,1}, {1,1,0,0}, {0,0,0,0}};
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    int * src_corner_mask =
      xmalloc(grid_data.num_vertices * sizeof(*src_corner_mask));
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      src_corner_mask[i] =
        ((int*)(&(global_corner_mask[0][0])))[grid_data.vertex_ids[i]];

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    yac_basic_grid_add_mask_nocpy(src_grid, CORNER, src_corner_mask, NULL);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_empty_new(tgt_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {
      for (size_t collection_size = 1; collection_size < 4;
           collection_size += 2) {

        struct interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        //------------------------------
        // check generated interpolation
        //------------------------------
        {
          double *** src_data = xmalloc(collection_size * sizeof(*src_data));

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            src_data[collection_idx] = xmalloc(1 * sizeof(**src_data));
            src_data[collection_idx][0] =
              xmalloc(grid_data.num_vertices * sizeof(***src_data));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              src_data[collection_idx][0][i] =
                (double)(grid_data.vertex_ids[i]) +
                (double)(collection_idx * 12);
          }

          yac_interpolation_execute_put(interpolation, src_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            free(src_data[collection_idx][0]);
            free(src_data[collection_idx]);
          }
          free(src_data);
        }

        //---------------
        // cleanup
        //---------------

        yac_interpolation_delete(interpolation);
      }
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (my_source_rank == 0) unlink(file_name);
  }

  { // simple weight file interpolation (multiple links per target point)
    // with two source processes and a single target process
    // there are two source fields

    // the global grid has 4x2 cells:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    //---------------
    // setup
    //---------------

    // create weight file
    if (my_source_rank == 0) {

      int tgt_indices[] = { 0, 0, 0, 0,
                            1, 1, 1, 1,
                            2, 2, 2, 2,
                            3, 3, 3, 3,
                            4, 4, 4, 4,
                            5, 5, 5, 5};
      int src_indices[] = { 0, 1, 4, 5,
                            1, 2, 5, 6,
                            2, 3, 6, 7,
                            4, 5, 8, 9,
                            5, 6, 9,10,
                            6, 7,10,11};
      double weights[] = {0.1,0.2,0.3,0.4,
                          0.5,0.6,0.7,0.8,
                          0.9,1.0,1.1,1.2,
                          1.3,1.4,1.5,1.6,
                          1.7,1.8,1.9,2.0,
                          2.1,2.2,2.3,2.4};
      size_t num_links = 24;
      enum yac_location src_locations[2] = {CORNER, CELL};
      enum yac_location tgt_location = CORNER;
      unsigned num_src_fields = 2;
      int num_links_per_field[2] = {num_links, 0};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0};
    size_t const num_global_cells[2] = {3,2};
    size_t local_start[2][2] = {{0,0},{2,0}};
    size_t local_count[2][2] = {{2,2},{1,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_empty_new(tgt_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------
      {
        double ** src_fields = xmalloc(1 * sizeof(*src_fields));
        double *** src_data = &src_fields;

        double * src_field =
          xmalloc(grid_data.num_vertices * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_vertices; ++i)
          src_field[i] = (double)(grid_data.vertex_ids[i]);
        src_fields[0] = src_field;

        yac_interpolation_execute_put(interpolation, src_data);

        free(src_fields[0]);
        free(src_fields);
      }

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (my_source_rank == 0) unlink(file_name);
  }

  { // simple weight file interpolation (multiple links per target point)
    // with two source processes and a single target process
    // there are two source fields and a source mask

    // the global grid has 4x2 cells:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    // source mask:
    //
    //       0-------1-------1-------1
    //       |       |       |       |
    //       |   0   |   1   |   1   |
    //       |       |       |       |
    //       1-------1-------1-------1
    //       |       |       |       |
    //       |   1   |   1   |   0   |
    //       |       |       |       |
    //       1-------1-------1-------0
    //
    //---------------
    // setup
    //---------------

    // create weight file
    if (my_source_rank == 0) {

      int tgt_indices[] = { 0, 0, 0, 0,
                            1, 1, 1, 1,
                            2, 2, 2, 2,
                            3, 3, 3, 3,
                            4, 4, 4, 4,
                            5, 5, 5, 5};
      int src_indices[] = { 0, 1, 4, 5,
                            1, 2, 5, 6,
                            2, 3, 6, 7,
                            4, 5, 8, 9,
                            5, 6, 9,10,
                            6, 7,10,11};
      double weights[] = {0.1,0.2,0.3,0.4,
                          0.5,0.6,0.7,0.8,
                          0.9,1.0,1.1,1.2,
                          1.3,1.4,1.5,1.6,
                          1.7,1.8,1.9,2.0,
                          2.1,2.2,2.3,2.4};
      size_t num_links = 24;
      enum yac_location src_locations[2] = {CORNER, CELL};
      enum yac_location tgt_location = CORNER;
      unsigned num_src_fields = 2;
      int num_links_per_field[2] = {num_links, 0};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0};
    size_t const num_global_cells[2] = {3,2};
    size_t local_start[2][2] = {{0,0},{2,0}};
    size_t local_count[2][2] = {{2,2},{1,2}};
    int global_src_vertex_mask[] = {1, 1, 1, 0,
                                    1, 1, 1, 1,
                                    0, 1, 1, 1};
    int global_src_cell_mask[] = {1, 1, 0,
                                  0, 1, 1};
    int with_halo = 1;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    int * src_vertex_mask = xmalloc(12 * sizeof(*src_vertex_mask));
    int * src_cell_mask = xmalloc(6 * sizeof(*src_cell_mask));;
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      src_vertex_mask[i] =
        global_src_vertex_mask[grid_data.vertex_ids[i]];
    for (size_t i = 0; i < grid_data.num_cells; ++i)
      src_cell_mask[i] =
        global_src_cell_mask[grid_data.cell_ids[i]];

    struct yac_basic_grid * src_grid = yac_basic_grid_new(src_grid_name, grid_data);
    yac_basic_grid_add_mask_nocpy(src_grid, CORNER, src_vertex_mask, NULL);
    yac_basic_grid_add_mask_nocpy(src_grid, CELL, src_cell_mask, NULL);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_empty_new(tgt_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);
    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------
      {
        double ** src_fields = xmalloc(1 * sizeof(*src_fields));
        double *** src_data = &src_fields;

        double * src_field =
          xmalloc(grid_data.num_vertices * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_vertices; ++i)
          src_field[i] = (double)(grid_data.vertex_ids[i]);
        src_fields[0] = src_field;

        yac_interpolation_execute_put(interpolation, src_data);

        free(src_fields[0]);
        free(src_fields);
      }

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (my_source_rank == 0) unlink(file_name);
  }

  { // simple weight file interpolation (weight file contains no weights and
    // not fixed target points) with two source processes and
    // a single target process the source fields have two struct points per
    // grid

    // the global grid has 4x2 cells:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    //---------------
    // setup
    //---------------

    // create weight file
    if (my_source_rank == 0) {

      int * tgt_indices = NULL;
      int * src_indices = NULL;
      double * weights = NULL;
      size_t num_links = 0;
      enum yac_location src_locations[2] = {CORNER, CELL};
      enum yac_location tgt_location = CORNER;
      unsigned num_src_fields = 2;
      int num_links_per_field[2] = {num_links, 0};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0};
    size_t const num_global_cells[2] = {3,2};
    size_t local_start[2][2] = {{0,0},{2,0}};
    size_t local_count[2][2] = {{2,2},{1,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_empty_new(tgt_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------
      {
        double ** src_fields = xmalloc(1 * sizeof(*src_fields));
        double *** src_data = &src_fields;

        double * src_field =
          xmalloc(grid_data.num_vertices * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_vertices; ++i)
          src_field[i] = (double)(grid_data.vertex_ids[i]);
        src_fields[0] = src_field;

        yac_interpolation_execute_put(interpolation, src_data);

        free(src_fields[0]);
        free(src_fields);
      }

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (my_source_rank == 0) unlink(file_name);
  }

  { // simple weight file interpolation (weight file contains only fixed target
    // points) with two source processes and a single target process

    // the global grid has 4x2 cells:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    //---------------
    // setup
    //---------------

    // create weight file
    if (my_source_rank == 0) {

      int * tgt_indices = NULL;
      int * src_indices = NULL;
      double * weights = NULL;
      size_t num_links = 0;
      enum yac_location src_locations[2] = {CORNER, CELL};
      enum yac_location tgt_location = CORNER;
      unsigned num_src_fields = 2;
      int num_links_per_field[2] = {num_links, 0};
      int tgt_id_fixed[] = {0, 2, 4, 6};
      size_t num_fixed_tgt = 4;
      double fixed_values[] = {-1.0, -2.0};
      int num_tgt_per_fixed_value[] = {1, 3};
      size_t num_fixed_values = 2;

      write_weight_file(
        file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0};
    size_t const num_global_cells[2] = {3,2};
    size_t local_start[2][2] = {{0,0},{2,0}};
    size_t local_count[2][2] = {{2,2},{1,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_empty_new(tgt_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------
      {
        double ** src_fields = xmalloc(1 * sizeof(*src_fields));
        double *** src_data = &src_fields;

        double * src_field =
          xmalloc(grid_data.num_vertices * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_vertices; ++i)
          src_field[i] = (double)(grid_data.vertex_ids[i]);
        src_fields[0] = src_field;

        yac_interpolation_execute_put(interpolation, src_data);

        free(src_fields[0]);
        free(src_fields);
      }

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (my_source_rank == 0) unlink(file_name);
  }

  { // simple weight file interpolation (weight file contains only fixed target
    // points) with two source processes and a single target process

    // the global grid has 4x2 cells:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    //---------------
    // setup
    //---------------

    // create weight file
    if (my_source_rank == 0) {

      int * tgt_indices = NULL;
      int * src_indices = NULL;
      double * weights = NULL;
      size_t num_links = 0;
      enum yac_location src_locations[] = {CORNER};
      enum yac_location tgt_location = EDGE;
      unsigned num_src_fields = 1;
      int num_links_per_field[1] = {num_links};
      int tgt_id_fixed[] = {1,3,5,7,9,11,0,2,4,6,8,10};
      size_t num_fixed_tgt = 12;
      double fixed_values[] = {-1.0, -2.0};
      int num_tgt_per_fixed_value[] = {6, 6};
      size_t num_fixed_values = 2;

      write_weight_file(
        file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0};
    size_t const num_global_cells[2] = {3,2};
    size_t local_start[2][2] = {{0,0},{2,0}};
    size_t local_count[2][2] = {{2,2},{1,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_empty_new(tgt_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = EDGE, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------
      {
        double ** src_fields = xmalloc(1 * sizeof(*src_fields));
        double *** src_data = &src_fields;

        double * src_field =
          xmalloc(grid_data.num_vertices * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_vertices; ++i)
          src_field[i] = (double)(grid_data.vertex_ids[i]);
        src_fields[0] = src_field;

        yac_interpolation_execute_put(interpolation, src_data);

        free(src_fields[0]);
        free(src_fields);
      }

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (my_source_rank == 0) unlink(file_name);
  }

  { // multi source field weight file interpolation
    // with two source processes and a single target process

    // the global grid has 2x2 cells:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    // weight_type == 0 => weights have different values
    // weight_type == 1 => all weights have the value 1.0
    // weight_type == 2 => each target gets assigned a value
    //                     from a single source point (from
    //                     varying source field)
    for (int weight_type = 0; weight_type <= 2; ++weight_type) {

      // create weight file
      if (my_source_rank == 0) {

        int src_indices[3][36] = {{0,1,2,3,
                                   0,1,3,4, 1,2,4,5, 3,4,6,7, 4,5,7,8,
                                   0,1,3,5, 2,3,4,7, 5,6,8,10, 7,8,9,11},
                                  {0,1,2,3,
                                   0,1,3,4, 1,2,4,5, 3,4,6,7, 4,5,7,8,
                                   0,1,3,5, 2,3,4,7, 5,6,8,10, 7,8,9,11},
                                  {1,2,
                                   4,
                                   11}};
        int tgt_indices[3][36] = {{0,1,2,3,
                                   0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3,
                                   0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3},
                                  {0,1,2,3,
                                   0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3,
                                   0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3},
                                  {1,2,
                                   0,
                                   3}};
        double weights[3][36] = {{1,1,1,1,
                                  0.25,0.25,0.25,0.25, 0.25,0.25,0.25,0.25,
                                  0.25,0.25,0.25,0.25, 0.25,0.25,0.25,0.25,
                                  0.25,0.25,0.25,0.25, 0.25,0.25,0.25,0.25,
                                  0.25,0.25,0.25,0.25, 0.25,0.25,0.25,0.25},
                                 {1,1,1,1,
                                  1,1,1,1, 1,1,1,1,
                                  1,1,1,1, 1,1,1,1,
                                  1,1,1,1, 1,1,1,1,
                                  1,1,1,1, 1,1,1,1},
                                 {1,1,1,1}};
        size_t num_links[3] = {1*4 + 4*4 + 4*4,
                               1*4 + 4*4 + 4*4,
                               1*2 + 1*1 + 1*1};
        enum yac_location src_locations[3] = {CELL, CORNER, EDGE};
        enum yac_location tgt_location = CORNER;
        unsigned num_src_fields = 3;
        int num_links_per_field[3][3] = {{1*4, 4*4, 4*4},
                                         {1*4, 4*4, 4*4},
                                         {2, 1, 1}};
        int * tgt_id_fixed = NULL;
        size_t num_fixed_tgt = 0;
        double * fixed_values = NULL;
        int * num_tgt_per_fixed_value = NULL;
        size_t num_fixed_values = 0;

        write_weight_file(
          file_name, src_indices[weight_type], tgt_indices[weight_type],
          weights[weight_type], num_links[weight_type], src_locations,
          num_src_fields, num_links_per_field[weight_type], tgt_id_fixed,
          num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
          num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
      }

      double coordinates_x[] = {0.0, 1.0, 2.0};
      double coordinates_y[] = {0.0, 1.0, 2.0};
      size_t const num_global_cells[2] = {2,2};
      size_t local_start[2][2] = {{0,0},{1,0}};
      size_t local_count[2][2] = {{1,2},{1,2}};
      int with_halo = 1;
      for (size_t i = 0; i <= num_global_cells[0]; ++i)
        coordinates_x[i] *= YAC_RAD;
      for (size_t i = 0; i <= num_global_cells[1]; ++i)
        coordinates_y[i] *= YAC_RAD;

      struct basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_global_cells,
          local_start[my_source_rank], local_count[my_source_rank], with_halo);

      struct yac_basic_grid * src_grid =
        yac_basic_grid_new(src_grid_name, grid_data);
      struct yac_basic_grid * tgt_grid =
        yac_basic_grid_empty_new(tgt_grid_name);

      struct dist_grid_pair * grid_pair =
        yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

      struct interp_field src_fields[] =
        {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
         {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
         {.location = EDGE, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

      struct interp_grid * interp_grid =
        yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                            num_src_fields, src_fields, tgt_field);

      //----------------------------------------
      // test generation of interpolation method
      //----------------------------------------

      //-----------------
      // generate weights
      //-----------------

      struct interp_method * method_stack[2] = {
        yac_interp_method_file_new(
          file_name, src_grid_name, tgt_grid_name), NULL};
      struct interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);
      yac_interp_method_delete(method_stack);

      // write weights to file

      yac_interp_weights_write_to_file(
        weights, file_name_2, src_grid_name, tgt_grid_name);

      // read weights
      struct interp_method * method_stack_file[2] = {
        yac_interp_method_file_new(
          file_name_2, src_grid_name, tgt_grid_name), NULL};
      struct interp_weights * weights_from_file =
        yac_interp_method_do_search(method_stack_file, interp_grid);
      yac_interp_method_delete(method_stack_file);

      enum interp_weights_reorder_type reorder_type[2] =
        {MAPPING_ON_SRC, MAPPING_ON_TGT};

      for (int from_file = 0; from_file <= 1; ++from_file) {

        for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
             ++i) {

          for (size_t collection_size = 1; collection_size <= 16;
               collection_size *= 2) {

            struct interpolation * interpolation =
              yac_interp_weights_get_interpolation(
                (from_file)?weights_from_file:weights,
                reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

            //------------------------------
            // check generated interpolation
            //------------------------------
            {
              double *** src_data = xmalloc(collection_size * sizeof(*src_data));

              for (size_t collection_idx = 0; collection_idx < collection_size;
                   ++collection_idx) {

                double ** src_fields = xmalloc(3 * sizeof(*src_fields));
                src_data[collection_idx] = src_fields;

                {
                  double * src_field =
                    xmalloc(grid_data.num_cells * sizeof(*src_field));
                  for (size_t i = 0; i < grid_data.num_cells; ++i)
                    src_field[i] =
                      (double)(grid_data.cell_ids[i] + 10 * collection_idx) + 0.0;
                  src_fields[0] = src_field;
                }
                {
                  double * src_field =
                    xmalloc(grid_data.num_vertices * sizeof(*src_field));
                  for (size_t i = 0; i < grid_data.num_vertices; ++i)
                    src_field[i] =
                      (double)(grid_data.vertex_ids[i] + 10 * collection_idx) + 0.3;
                  src_fields[1] = src_field;
                }
                {
                  double * src_field =
                    xmalloc(grid_data.num_edges * sizeof(*src_field));
                  for (size_t i = 0; i < grid_data.num_edges; ++i)
                    src_field[i] =
                      (double)(grid_data.edge_ids[i] + 10 * collection_idx) + 0.7;
                  src_fields[2] = src_field;
                }
              }

              yac_interpolation_execute_put(interpolation, src_data);

              for (size_t collection_idx = 0; collection_idx < collection_size;
                   ++collection_idx) {
                for (size_t i = 0; i < 3; ++i) free(src_data[collection_idx][i]);
                free(src_data[collection_idx]);
              }
              free(src_data);
            }

            //---------------
            // cleanup
            //---------------

            yac_interpolation_delete(interpolation);
          }
        }
      }

      //---------------
      // cleanup
      //---------------

      yac_interp_weights_delete(weights_from_file);
      yac_interp_weights_delete(weights);
      yac_interp_grid_delete(interp_grid);
      yac_dist_grid_pair_delete(grid_pair);
      yac_basic_grid_delete(tgt_grid);
      yac_basic_grid_delete(src_grid);

      if (my_source_rank == 0) {
        unlink(file_name);
        unlink(file_name_2);
      }
    }
  }

  { // multi source field weight file interpolation
    // with two source processes and a single target process

    // the global source grid is a 4x4 grid:
    //
    //       20--36--21--37--22--38--23--39--24
    //       |       |       |       |       |
    //       28  12  30  13  32  14  34  15  35
    //       |       |       |       |       |
    //       15--27--16--29--17--31--18--33--19
    //       |       |       |       |       |
    //       19  08  21  09  23  10  25  11  26
    //       |       |       |       |       |
    //       10--18--11--20--12--22--13--24--14
    //       |       |       |       |       |
    //       10  04  12  05  14  06  16  07  17
    //       |       |       |       |       |
    //       05--09--06--11--07--13--08--15--09
    //       |       |       |       |       |
    //       01  00  03  01  05  02  07  03  08
    //       |       |       |       |       |
    //       00--00--01--02--02--04--03--06--04
    //
    //---------------
    // setup
    //---------------

    // create weight file
    if (my_source_rank == 0) {

      int src_indices[16] = {0,1,4,8,13,14,
                             8,12,13,17,24,
                             8,17,25,26,30};
      int tgt_indices[16] = {0,1,4,8,13,14,
                             2,5,6,9,15,
                             3,7,10,11,12};
      double weights[16] = {1,1,1,1,1,1,
                            1,1,1,1,1,
                            1,1,1,1,1};
      size_t num_links = 16;
      enum yac_location src_locations[3] = {CELL, CORNER, EDGE};
      enum yac_location tgt_location = CELL;
      unsigned num_src_fields = 3;
      int num_links_per_field[3] = {6, 5, 5};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
    }

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double coordinates_y[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    size_t const num_global_cells[2] = {4,4};
    size_t local_start[2][2] = {{0,0},{2,0}};
    size_t local_count[2][2] = {{2,4},{2,4}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_empty_new(tgt_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = EDGE, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    //-----------------
    // generate weights
    //-----------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};
    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);
    yac_interp_method_delete(method_stack);

    // write weights to file

    yac_interp_weights_write_to_file(
      weights, file_name_2, src_grid_name, tgt_grid_name);

    // read weights
    struct interp_method * method_stack_file[2] = {
      yac_interp_method_file_new(
        file_name_2, src_grid_name, tgt_grid_name), NULL};
    struct interp_weights * weights_from_file =
      yac_interp_method_do_search(method_stack_file, interp_grid);
    yac_interp_method_delete(method_stack_file);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (int from_file = 0; from_file <= 1; ++from_file) {

      for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
           ++i) {

        for (size_t collection_size = 1; collection_size <= 16;
             collection_size *= 2) {

          struct interpolation * interpolation =
            yac_interp_weights_get_interpolation(
              (from_file)?weights_from_file:weights,
              reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

          //------------------------------
          // check generated interpolation
          //------------------------------
          {
            double *** src_data = xmalloc(collection_size * sizeof(*src_data));

            for (size_t collection_idx = 0; collection_idx < collection_size;
                 ++collection_idx) {

              double ** src_fields = xmalloc(3 * sizeof(*src_fields));
              src_data[collection_idx] = src_fields;

              {
                double * src_field =
                  xmalloc(grid_data.num_cells * sizeof(*src_field));
                for (size_t i = 0; i < grid_data.num_cells; ++i)
                  src_field[i] =
                    (double)(grid_data.cell_ids[i] + 10 * collection_idx) + 0.0;
                src_fields[0] = src_field;
              }
              {
                double * src_field =
                  xmalloc(grid_data.num_vertices * sizeof(*src_field));
                for (size_t i = 0; i < grid_data.num_vertices; ++i)
                  src_field[i] =
                    (double)(grid_data.vertex_ids[i] + 10 * collection_idx) + 0.3;
                src_fields[1] = src_field;
              }
              {
                double * src_field =
                  xmalloc(grid_data.num_edges * sizeof(*src_field));
                for (size_t i = 0; i < grid_data.num_edges; ++i)
                  src_field[i] =
                    (double)(grid_data.edge_ids[i] + 10 * collection_idx) + 0.7;
                src_fields[2] = src_field;
              }
            }

            yac_interpolation_execute_put(interpolation, src_data);

            for (size_t collection_idx = 0; collection_idx < collection_size;
                 ++collection_idx) {
              for (size_t i = 0; i < 3; ++i) free(src_data[collection_idx][i]);
              free(src_data[collection_idx]);
            }
            free(src_data);
          }

          //---------------
          // cleanup
          //---------------

          yac_interpolation_delete(interpolation);
        }
      }
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights_from_file);
    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (my_source_rank == 0) {
      unlink(file_name);
      unlink(file_name_2);
    }
  }
}

static void target_main(MPI_Comm global_comm,
                        MPI_Comm target_comm) {

  int my_target_rank;
  yac_mpi_call(MPI_Comm_rank(target_comm, &my_target_rank), target_comm);

  { // simple weight file interpolation (one link per target point) with two
    // source processes and a single target process the source fields have two
    // struct points per grid (the test check how the interpolation
    // methods handle NULL target coordinates )

    // the global target grid has 2x2 cells:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    double coordinates_x[] = {0.5, 1.5, 2.5};
    double coordinates_y[] = {0.5, 1.5, 2.5};
    size_t const num_global_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    struct yac_basic_grid * src_grid =
      yac_basic_grid_empty_new(src_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {
      for (size_t collection_size = 1; collection_size < 4;
           collection_size += 2) {

        struct interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        //------------------------------
        // check generated interpolation
        //------------------------------

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            tgt_data[collection_idx][i] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            if (i >= 6) {
              if (tgt_data[collection_idx][i] != -1.0)
                PUT_ERR("wrong interpolation result");
            } else if (fabs((double)(i * (i + 12 * collection_idx)) -
                         tgt_data[collection_idx][i]) > 1e-9) {
                PUT_ERR("wrong interpolation result");
            }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);

        //---------------
        // cleanup
        //---------------

        yac_interpolation_delete(interpolation);
      }
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  { // simple weight file interpolation (multiple links per target point)
    // with two source processes and a single target process
    // there are two source fields

    // the global target grid has 2x2 cells:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    double coordinates_x[] = {0.5, 1.5, 2.5};
    double coordinates_y[] = {0.5, 1.5, 2.5};
    size_t const num_global_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    struct yac_basic_grid * src_grid =
      yac_basic_grid_empty_new(src_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------

      double ** tgt_data = xmalloc(1 * sizeof(*tgt_data));
      tgt_data[0] =
        xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        tgt_data[0][i] = -1;

      double ref_tgt_field[9] =
        {0.1*0+0.2*1+0.3*4+0.4*5,
         0.5*1+0.6*2+0.7*5+0.8*6,
         0.9*2+1.0*3+1.1*6+1.2*7,
         1.3*4+1.4*5+1.5*8+1.6*9,
         1.7*5+1.8*6+1.9*9+2.0*10,
         2.1*6+2.2*7+2.3*10+2.4*11,
         -1,-1,-1};

      yac_interpolation_execute_get(interpolation, tgt_data);

      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        if (i >= 6) {
          if (tgt_data[0][i] != -1.0)
            PUT_ERR("wrong interpolation result");
        } else if (fabs(ref_tgt_field[i] - tgt_data[0][i]) > 1e-9) {
            PUT_ERR("wrong interpolation result");
        }

      free(tgt_data[0]);
      free(tgt_data);

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  { // simple weight file interpolation (multiple links per target point)
    // with two source processes and a single target process
    // there are two source fields and source masks

    // the global target grid has 2x2 cells:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    double coordinates_x[] = {0.5, 1.5, 2.5};
    double coordinates_y[] = {0.5, 1.5, 2.5};
    size_t const num_global_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    struct yac_basic_grid * src_grid =
      yac_basic_grid_empty_new(src_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------

      double ** tgt_data = xmalloc(1 * sizeof(*tgt_data));
      tgt_data[0] =
        xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        tgt_data[0][i] = -1;

      double ref_tgt_field[9] =
        {0.1*0+0.2*1+0.3*4+0.4*5,
         0.5*1+0.6*2+0.7*5+0.8*6,
         -1,
         -1,
         1.7*5+1.8*6+1.9*9+2.0*10,
         2.1*6+2.2*7+2.3*10+2.4*11,
         -1,-1,-1};

      yac_interpolation_execute_get(interpolation, tgt_data);

      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        if (i >= 6) {
          if (tgt_data[0][i] != -1.0)
            PUT_ERR("wrong interpolation result");
        } else if (fabs(ref_tgt_field[i] - tgt_data[0][i]) > 1e-9) {
            PUT_ERR("wrong interpolation result");
        }

      free(tgt_data[0]);
      free(tgt_data);

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  { // simple weight file interpolation (weight file contains only fixed target
    // points) with two source processes and a single target process

    // the global target grid has 2x2 cells:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    double coordinates_x[] = {0.5, 1.5, 2.5};
    double coordinates_y[] = {0.5, 1.5, 2.5};
    size_t const num_global_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    struct yac_basic_grid * src_grid =
      yac_basic_grid_empty_new(src_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------

      double ** tgt_data = xmalloc(1 * sizeof(*tgt_data));
      tgt_data[0] =
        xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        tgt_data[0][i] = -1;

      yac_interpolation_execute_get(interpolation, tgt_data);

      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        if (tgt_data[0][i] != -1.0)
          PUT_ERR("wrong interpolation result");

      free(tgt_data[0]);
      free(tgt_data);

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  { // simple weight file interpolation (weight file contains only fixed target
    // points) with two source processes and a single target process

    // the global target grid has 2x2 cells:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    double coordinates_x[] = {0.5, 1.5, 2.5};
    double coordinates_y[] = {0.5, 1.5, 2.5};
    size_t const num_global_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    struct yac_basic_grid * src_grid =
      yac_basic_grid_empty_new(src_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------

      double ** tgt_data = xmalloc(1 * sizeof(*tgt_data));
      tgt_data[0] =
        xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
      for (size_t i = 0; i < grid_data.num_vertices; ++i) tgt_data[0][i] = -3;

      double ref_tgt_field[9] =
        {-1,-3,-2, -3,-2,-3, -2,-3,-3};

      yac_interpolation_execute_get(interpolation, tgt_data);

      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        if (ref_tgt_field[i] != tgt_data[0][i])
          PUT_ERR("wrong interpolation result");

      free(tgt_data[0]);
      free(tgt_data);

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  { // simple weight file interpolation (weight file contains only fixed target
    // points) with two source processes and a single target process

    // the global target grid has 2x2 cells:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    double coordinates_x[] = {0.5, 1.5, 2.5};
    double coordinates_y[] = {0.5, 1.5, 2.5};
    size_t const num_global_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    struct yac_basic_grid * src_grid =
      yac_basic_grid_empty_new(src_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = EDGE, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};

    //-----------------
    // generate weights
    //-----------------

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
         ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      //------------------------------
      // check generated interpolation
      //------------------------------

      double ** tgt_data = xmalloc(1 * sizeof(*tgt_data));
      tgt_data[0] = xmalloc(grid_data.num_edges * sizeof(**tgt_data));
      for (size_t i = 0; i < grid_data.num_edges; ++i) tgt_data[0][i] = -3;

      double ref_tgt_field[12] = {-2,-1,-2,-1,-2,-1,-2,-1,-2,-1,-2,-1};

      yac_interpolation_execute_get(interpolation, tgt_data);

      for (size_t i = 0; i < grid_data.num_edges; ++i)
        if (ref_tgt_field[i] != tgt_data[0][i])
          PUT_ERR("wrong interpolation result");

      free(tgt_data[0]);
      free(tgt_data);

      //---------------
      // cleanup
      //---------------

      yac_interpolation_delete(interpolation);
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  { // multi source field weight file interpolation
    // with two source processes and a single target process

    // the global target grid has 1x1 cells:
    //
    //       02--03--03
    //       |       |
    //       01  00  02
    //       |       |
    //       00--00--01
    //
    //---------------
    // setup
    //---------------

    // weight_type == 0 => weights have different values
    // weight_type == 1 => all weights have the value 1.0
    // weight_type == 2 => each target gets assigned a value
    //                     from a single source point (from
    //                     varying source field)
    for (int weight_type = 0; weight_type <= 2; ++weight_type) {

      double coordinates_x[] = {0.5, 1.5};
      double coordinates_y[] = {0.5, 1.5};
      size_t const num_global_cells[2] = {1,1};
      size_t local_start[2] = {0,0};
      size_t local_count[2] = {1,1};
      int with_halo = 0;
      for (size_t i = 0; i <= num_global_cells[0]; ++i)
        coordinates_x[i] *= YAC_RAD;
      for (size_t i = 0; i <= num_global_cells[1]; ++i)
        coordinates_y[i] *= YAC_RAD;

      struct basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_global_cells,
          local_start, local_count, with_halo);

      struct yac_basic_grid * tgt_grid =
        yac_basic_grid_new(tgt_grid_name, grid_data);
      struct yac_basic_grid * src_grid =
        yac_basic_grid_empty_new(src_grid_name);

      struct dist_grid_pair * grid_pair =
        yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

      struct interp_field src_fields[] =
        {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
         {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
         {.location = EDGE, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

      struct interp_grid * interp_grid =
        yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                            num_src_fields, src_fields, tgt_field);

      //----------------------------------------
      // test generation of interpolation method
      //----------------------------------------

      //-----------------
      // generate weights
      //-----------------

      struct interp_method * method_stack[2] = {
        yac_interp_method_file_new(
          file_name, src_grid_name, tgt_grid_name), NULL};
      struct interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);
      yac_interp_method_delete(method_stack);

      // write weights to file
      yac_interp_weights_write_to_file(
        weights, file_name_2, src_grid_name, tgt_grid_name);

      // read weights
      struct interp_method * method_stack_file[2] = {
        yac_interp_method_file_new(
          file_name_2, src_grid_name, tgt_grid_name), NULL};
      struct interp_weights * weights_from_file =
        yac_interp_method_do_search(method_stack_file, interp_grid);
      yac_interp_method_delete(method_stack_file);

      enum interp_weights_reorder_type reorder_type[2] =
        {MAPPING_ON_SRC, MAPPING_ON_TGT};

      for (int from_file = 0; from_file <= 1; ++from_file) {

        for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
             ++i) {

          for (size_t collection_size = 1; collection_size <= 16;
               collection_size *= 2) {

            struct interpolation * interpolation =
              yac_interp_weights_get_interpolation(
                (from_file)?weights_from_file:weights,
                reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

            //------------------------------
            // check generated interpolation
            //------------------------------

            double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
            for (size_t collection_idx = 0; collection_idx < collection_size;
                 ++collection_idx) {
              tgt_data[collection_idx] =
                xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
              for (size_t i = 0; i < grid_data.num_vertices; ++i)
                tgt_data[collection_idx][i] = -3;
            }

            double ref_tgt_field[3][4] =
              {{0.0 + 0.25 * (0.3+1.3+3.3+4.3) + 0.25 * (0.7+1.7+3.7+ 5.7),
                1.0 + 0.25 * (1.3+2.3+4.3+5.3) + 0.25 * (2.7+3.7+4.7+ 7.7),
                2.0 + 0.25 * (3.3+4.3+6.3+7.3) + 0.25 * (5.7+6.7+8.7+10.7),
                3.0 + 0.25 * (4.3+5.3+7.3+8.3) + 0.25 * (7.7+8.7+9.7+11.7)},
               {0.0 + (0.3+1.3+3.3+4.3) + (0.7+1.7+3.7+ 5.7),
                1.0 + (1.3+2.3+4.3+5.3) + (2.7+3.7+4.7+ 7.7),
                2.0 + (3.3+4.3+6.3+7.3) + (5.7+6.7+8.7+10.7),
                3.0 + (4.3+5.3+7.3+8.3) + (7.7+8.7+9.7+11.7)},
               {4.3, 1.0, 2.0, 11.7}};
            double collection_factor[3] = {30.0, 90.0, 10.0};

            yac_interpolation_execute_get(interpolation, tgt_data);

            for (size_t collection_idx = 0; collection_idx < collection_size;
                 collection_idx++)
              for (size_t i = 0; i < grid_data.num_vertices; ++i)
                if (fabs((ref_tgt_field[weight_type][i] +
                          collection_factor[weight_type] *
                          (double)collection_idx) -
                         tgt_data[collection_idx][i]) > 1e-9)
                  PUT_ERR("wrong interpolation result");

            for (size_t collection_idx = 0; collection_idx < collection_size;
                 ++collection_idx) free(tgt_data[collection_idx]);
            free(tgt_data);

            //---------------
            // cleanup
            //---------------

            yac_interpolation_delete(interpolation);
          }
        }
      }

      //---------------
      // cleanup
      //---------------

      yac_interp_weights_delete(weights_from_file);
      yac_interp_weights_delete(weights);
      yac_interp_grid_delete(interp_grid);
      yac_dist_grid_pair_delete(grid_pair);
      yac_basic_grid_delete(src_grid);
      yac_basic_grid_delete(tgt_grid);
    }
  }

  { // multi source field weight file interpolation
    // with two source processes and a single target process

    // the global target grid is a 4x4 grid:
    //
    //       20--36--21--37--22--38--23--39--24
    //       |       |       |       |       |
    //       28  12  30  13  32  14  34  15  35
    //       |       |       |       |       |
    //       15--27--16--29--17--31--18--33--19
    //       |       |       |       |       |
    //       19  08  21  09  23  10  25  11  26
    //       |       |       |       |       |
    //       10--18--11--20--12--22--13--24--14
    //       |       |       |       |       |
    //       10  04  12  05  14  06  16  07  17
    //       |       |       |       |       |
    //       05--09--06--11--07--13--08--15--09
    //       |       |       |       |       |
    //       01  00  03  01  05  02  07  03  08
    //       |       |       |       |       |
    //       00--00--01--02--02--04--03--06--04
    //
    //---------------
    // setup
    //---------------

    // weight_type == 0 => weights have different values
    // weight_type == 1 => all weights have the value 1.0
    // weight_type == 2 => each target gets assigned a value
    //                     from a single source point (from
    //                     varying source field)

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double coordinates_y[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    size_t const num_global_cells[2] = {4,4};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {4,4};
    int with_halo = 0;
    for (size_t i = 0; i <= num_global_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_global_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_global_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    struct yac_basic_grid * src_grid =
      yac_basic_grid_empty_new(src_grid_name);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = EDGE, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    //-----------------
    // generate weights
    //-----------------

    struct interp_method * method_stack[2] = {
      yac_interp_method_file_new(
        file_name, src_grid_name, tgt_grid_name), NULL};
    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);
    yac_interp_method_delete(method_stack);

    // write weights to file
    yac_interp_weights_write_to_file(
      weights, file_name_2, src_grid_name, tgt_grid_name);

    // read weights
    struct interp_method * method_stack_file[2] = {
      yac_interp_method_file_new(
        file_name_2, src_grid_name, tgt_grid_name), NULL};
    struct interp_weights * weights_from_file =
      yac_interp_method_do_search(method_stack_file, interp_grid);
    yac_interp_method_delete(method_stack_file);

    enum interp_weights_reorder_type reorder_type[2] =
      {MAPPING_ON_SRC, MAPPING_ON_TGT};

    for (int from_file = 0; from_file <= 1; ++from_file) {

      for (size_t i = 0; i < sizeof(reorder_type) / sizeof(reorder_type[0]);
           ++i) {

        for (size_t collection_size = 1; collection_size <= 16;
             collection_size *= 2) {

          struct interpolation * interpolation =
            yac_interp_weights_get_interpolation(
              (from_file)?weights_from_file:weights,
              reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

          //------------------------------
          // check generated interpolation
          //------------------------------

          double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            tgt_data[collection_idx] =
              xmalloc(grid_data.num_cells * sizeof(**tgt_data));
            for (size_t i = 0; i < grid_data.num_cells; ++i)
              tgt_data[collection_idx][i] = -3;
          }

          double ref_tgt_field[16] =
            { 0.0,  1.0,  8.3,  8.7,
              4.0, 12.3, 13.3, 17.7,
              8.0, 17.3, 25.7, 26.7,
             30.7, 13.0, 14.0, 24.3};

          yac_interpolation_execute_get(interpolation, tgt_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               collection_idx++)
            for (size_t i = 0; i < grid_data.num_cells; ++i)
              if (fabs((ref_tgt_field[i] + 10.0 * (double)collection_idx) -
                       tgt_data[collection_idx][i]) > 1e-9)
                PUT_ERR("wrong interpolation result");

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) free(tgt_data[collection_idx]);
          free(tgt_data);

          //---------------
          // cleanup
          //---------------

          yac_interpolation_delete(interpolation);
        }
      }
    }

    //---------------
    // cleanup
    //---------------

    yac_interp_weights_delete(weights_from_file);
    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }
}
