/**
 * @file test_interpolation_parallel3.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
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

static void
generate_ref_weights();

static void submain_src(MPI_Comm global_comm, MPI_Comm comp_comm,
                        enum interp_weights_reorder_type reorder_type);

static void submain_tgt(MPI_Comm global_comm, MPI_Comm comp_comm,
                        enum interp_weights_reorder_type reorder_type);

// link data for input weight file
unsigned const num_links_file = 48;
int src_address_file[48] = {
  9,10,16,17, 10,11,17,18, 15,16,22,23, 16,17,23,24, 17,18,24,25, 18,19,25,26,
  22,23,29,30, 23,24,30,31, 24,25,31,32, 25,26,32,33, 30,31,37,38, 31,32,38,39};
int tgt_address_file[48] = {
  8,8,8,8, 9,9,9,9, 13,13,13,13, 14,14,14,14, 15,15,15,15, 16,16,16,16,
  19,19,19,19, 20,20,20,20, 21,21,21,21, 22,22,22,22, 26,26,26,26, 27,27,27,27};
double weights_file[48] = {
  0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4,
  0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4,
  0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4, 0.1,0.2,0.3,0.4};
char const weight_file_in[] = "weight_file_in.nc";
char const weight_file_out[] = "weight_file_out.nc";
char const src_grid_name[] = "src_grid";
char const tgt_grid_name[] = "tgt_grid";

// reference link data for output weight file
unsigned ref_num_links = 4 * 36;

int ref_src_address[4 * 36];
int ref_tgt_address[4 * 36];
double ref_weights[4 * 36];

int const * ref_tgt_address_fixed = NULL;
double const * ref_fixed_values = NULL;
int const * ref_num_tgt_per_fixed_value = NULL;
unsigned ref_num_fixed_values = 0;

/*
 * The grid is distributed among 4 processes
 *
 * The global grid has 6x6 cells:
 *
 *    42--79--43--81--44--83--45--85--46--87--47--89--48
 *    |       |       |       |       |       |       |
 *    66  30  68  31  70  32  72  33  74  34  76  35  78
 *    |       |       |       |       |       |       |
 *    35--65--36--67--37--69--38--71--39--73--40--75--41
 *    |       |       |       |       |       |       |
 *    53  24  55  25  57  26  59  27  61  28  63  29  64
 *    |       |       |       |       |       |       |
 *    28--52--29--54--30--56--31--58--32--60--33--62--34
 *    |       |       |       |       |       |       |
 *    40  18  42  19  44  20  46  21  48  22  50  23  51
 *    |       |       |       |       |       |       |
 *    21--39--22--41--23--43--24--45--25--47--26--49--27
 *    |       |       |       |       |       |       |
 *    27  12  29  13  31  14  33  15  35  16  37  17  38
 *    |       |       |       |       |       |       |
 *    14--26--15--28--16--30--17--32--18--34--19--36--20
 *    |       |       |       |       |       |       |
 *    14  06  16  07  18  08  20  09  22  10  24  11  25
 *    |       |       |       |       |       |       |
 *    07--13--08--15--09--17--10--19--11--21--12--23--13
 *    |       |       |       |       |       |       |
 *    01  00  03  01  05  02  07  03  09  04  11  05  12
 *    |       |       |       |       |       |       |
 *    00--00--01--02--02--04--03--06--04--08--05--10--06
 */

// grid information (is the same for source and target)
size_t num_cells[2] = {6,6};
double coordinates_x[] = {0,1,2,3,4,5,6};
double coordinates_y[] = {0,1,2,3,4,5,6};
double cell_coordinates_x[] = {0.5,1.5,2.5,3.5,4.5,5.5};
double cell_coordinates_y[] = {0.5,1.5,2.5,3.5,4.5,5.5};
double cell_coords[36][3];
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

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 8) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;
  for (size_t i = 0; i < 6; ++i)
    for (size_t j = 0; j < 6; ++j)
      LLtoXYZ_deg(
        cell_coordinates_x[j], cell_coordinates_y[i], cell_coords[6 * i + j]);

  generate_ref_weights();

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

  int my_rank;
  yac_mpi_call(MPI_Comm_rank(comp_comm, &my_rank), comp_comm);

  size_t local_start[4][2] = {{0,0},{3,0},{0,3},{3,3}};
  size_t local_count[4][2] = {{3,3},{3,3},{3,3},{3,3}};

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells,
      local_start[my_rank], local_count[my_rank], with_halo);
  size_t num_cells = grid_data.num_cells;

  coordinate_pointer cell_field_coords =
    xmalloc(num_cells * sizeof(*cell_field_coords));
  for (size_t i = 0; i < num_cells; ++i)
    memcpy(cell_field_coords[i], cell_coords[grid_data.cell_ids[i]],
           3 * sizeof(double));

  struct yac_basic_grid * tgt_grid =
    yac_basic_grid_new(tgt_grid_name, grid_data);
  yac_basic_grid_add_coordinates_nocpy(tgt_grid, CELL, cell_field_coords);
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
    double * target_data_field =
      xmalloc(num_cells * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_cells; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_get(interpolation, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    yac_interp_weights_write_to_file(
      weights, weight_file_out, src_grid_name, tgt_grid_name);

    double ref_global_target_data[6*6];

    for (size_t i = 0; i < 6*6; ++i) {
      ref_global_target_data[i] = 0;
      for (size_t j = 0; j < 4; ++j)
        ref_global_target_data[i] +=
          ref_weights[4 * i + j] * (double)ref_src_address[4 * i + j];
    }

    for (size_t i = 0; i < num_cells; ++i) {
      if (grid_data.core_cell_mask[i]) {
        if (fabs(
              target_data[0][i] -
              ref_global_target_data[grid_data.cell_ids[i]]) > 1e-10)
          PUT_ERR("error in interpolated data on target side\n");
      } else {
        if (target_data[0][i] != -1.0)
          PUT_ERR("error in interpolated data on target side\n");
      }
    }

    //--------
    // cleanup
    //--------

    free(target_data_field);
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
    write_weight_file(weight_file_in, src_address_file, tgt_address_file,
                      weights_file, num_links_file, src_locations,
                      1, (int*)&num_links_file, tgt_id_fixed, num_fixed_tgt,
                      fixed_values, num_tgt_per_fixed_value, num_fixed_values,
                      tgt_location, src_grid_name, tgt_grid_name);
  }

  size_t local_start[4][2] = {{0,0},{1,0},{3,0},{5,0}};
  size_t local_count[4][2] = {{1,6},{2,6},{2,6},{1,6}};

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells,
      local_start[my_rank], local_count[my_rank], with_halo);
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
      check_weight_file(weight_file_out, ref_src_address, ref_tgt_address,
                        ref_weights, ref_num_links, ref_src_locations,
                        1, (int*)&ref_num_links, ref_tgt_address_fixed,
                        ref_fixed_values, ref_num_tgt_per_fixed_value,
                        ref_num_fixed_values, ref_tgt_location,
                        src_grid_name, tgt_grid_name);
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
generate_ref_weights() {

  for (unsigned i = 0; i < 6; ++i) {
    for (unsigned j = 0; j < 6; ++j) {
      ref_src_address[4 * (i * 6 + j) + 0] = 0 + i * 7 + j;
      ref_src_address[4 * (i * 6 + j) + 1] = 1 + i * 7 + j;
      ref_src_address[4 * (i * 6 + j) + 2] = 7 + i * 7 + j;
      ref_src_address[4 * (i * 6 + j) + 3] = 8 + i * 7 + j;
    }
  }
  for (unsigned i = 0; i < ref_num_links; ++i) ref_tgt_address[i] = i / 4;
  for (unsigned i = 0; i < ref_num_links; ++i) ref_weights[i] = 0.25;
  for (unsigned i = 0; i < num_links_file; ++i)
    ref_weights[tgt_address_file[i] * 4 + (i & 3)] = weights_file[i];
}

