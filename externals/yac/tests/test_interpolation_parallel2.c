/**
 * @file test_interpolation_parallel2.c
 *
 * @copyright Copyright  (C)  2020 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "tests.h"
#include "interp_method.h"
#include "interp_method_avg.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"
#include "interpolation.h"
#include "test_common.h"


static void submain_1(MPI_Comm global_comm, MPI_Comm comp_comm,
                      enum interp_weights_reorder_type reorder_type);

static void submain_2(MPI_Comm global_comm, MPI_Comm comp_comm,
                      enum interp_weights_reorder_type reorder_type);

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

  if (comm_size != 4) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  // split processes into source an target

  int comp_flag = comm_rank < 2;

  MPI_Comm split_comm;
  yac_mpi_call(
    MPI_Comm_split(
      MPI_COMM_WORLD, comp_flag, 0, &split_comm), MPI_COMM_WORLD);

  if (comp_flag) submain_1(MPI_COMM_WORLD, split_comm, reorder_type);
  else           submain_2(MPI_COMM_WORLD, split_comm, reorder_type);

  yac_mpi_call(MPI_Comm_free(&split_comm), MPI_COMM_WORLD);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

/*
 * The grid is distributed among 2 processes.
 *
 * The global grid has 4x4 cells:
 *
 *    20--36--21--37--22--38--23--39--24
 *    |       |       |       |       |
 *    28  12  30  13  32  14  34  15  35
 *    |       |       |       |       |
 *    15--27--16--29--17--31--18--33--19
 *    |       |       |       |       |
 *    19  08  21  09  23  10  25  11  26
 *    |       |       |       |       |
 *    10--18--11--20--12--22--13--24--14
 *    |       |       |       |       |
 *    10  04  12  05  14  06  16  07  17
 *    |       |       |       |       |
 *    05--09--06--11--07--13--08--15--09
 *    |       |       |       |       |
 *    01  00  03  01  05  02  07  03  08
 *    |       |       |       |       |
 *    00--00--01--02--02--04--03--06--04
 *
 * The mask looks as follows (# = masked points)
 *
 *    +-------+-------+-------+-------+
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    #---#---#---#---#-------+-------+
 *    |       |       |       |       |
 *    #   #   #   #   #       |       |
 *    |       |       |       |       |
 *    #---#---#---#---#-------+-------+
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    +-------+-------+-------+-------+
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    +-------+-------+-------+-------+
 */
static void submain_2(MPI_Comm global_comm, MPI_Comm comp_comm,
                      enum interp_weights_reorder_type reorder_type) {

  char const local_grid_name[] = "grid2";
  char const remote_grid_name[] = "grid1";

  int my_rank;
  yac_mpi_call(MPI_Comm_rank(comp_comm, &my_rank), comp_comm);

  double coordinates_x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  double coordinates_y[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  double edge_coordinates_x[] = {
    0.25,0.0,1.25,1.0,2.25,2.0,3.25,3.0,4.0,
    0.25,0.0,1.25,1.0,2.25,2.0,3.25,3.0,4.0,
    0.25,0.0,1.25,1.0,2.25,2.0,3.25,3.0,4.0,
    0.25,0.0,1.25,1.0,2.25,2.0,3.25,3.0,4.0,
    0.25,1.25,2.25,3.25};
  double edge_coordinates_y[] = {
    0.0,0.25,0.0,0.25,0.0,0.25,0.0,0.25,0.25,
    1.0,1.25,1.0,1.25,1.0,1.25,1.0,1.25,1.25,
    2.0,2.25,2.0,2.25,2.0,2.25,2.0,2.25,2.25,
    3.0,3.25,3.0,3.25,3.0,3.25,3.0,3.25,3.25,
    4.0,4.0,4.0,4.0};
  double edge_coords[40][3];
  size_t const num_cells[2] = {4,4};
  size_t local_start[2][2] = {{0,0},{2,0}};
  size_t local_count[2][2] = {{2,4},{2,4}};
  int with_halo = 1;
  for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;
  for (size_t i = 0; i < 40; ++i)
    LLtoXYZ_deg(edge_coordinates_x[i], edge_coordinates_y[i], edge_coords[i]);

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells,
      local_start[my_rank], local_count[my_rank], with_halo);
  size_t num_vertices = grid_data.num_vertices;
  size_t num_edges = grid_data.num_edges;

  struct yac_basic_grid * local_grid =
    yac_basic_grid_new(local_grid_name, grid_data);
  struct yac_basic_grid * remote_grid =
    yac_basic_grid_empty_new(remote_grid_name);

  yac_int masked_corner_ids[] = {10,11,12,15,16,17};
  size_t num_masked_corners =
    sizeof(masked_corner_ids)/sizeof(masked_corner_ids[0]);

  int * corner_mask = xmalloc(num_vertices * sizeof(*corner_mask));
  for (size_t i = 0; i < num_vertices; ++i) {
    corner_mask[i] = 1;
    for (size_t j = 0; j < num_masked_corners; ++j)
      if (grid_data.vertex_ids[i] == masked_corner_ids[j]) corner_mask[i] = 0;
  }

  yac_int masked_edge_ids[] = {18,19,20,21,23,27,29};
  size_t num_masked_edges =
    sizeof(masked_edge_ids)/sizeof(masked_edge_ids[0]);

  int * edge_mask = xmalloc(num_edges * sizeof(*edge_mask));
  for (size_t i = 0; i < num_edges; ++i) {
    edge_mask[i] = 1;
    for (size_t j = 0; j < num_masked_edges; ++j)
      if (grid_data.edge_ids[i] == masked_edge_ids[j]) edge_mask[i] = 0;
  }

  coordinate_pointer edge_field_coords =
    xmalloc(num_edges * sizeof(*edge_field_coords));
  for (size_t i = 0; i < num_edges; ++i)
    memcpy(edge_field_coords[i], edge_coords[grid_data.edge_ids[i]],
           3 * sizeof(double));

  yac_basic_grid_add_mask_nocpy(local_grid, CORNER, corner_mask, NULL);
  yac_basic_grid_add_mask_nocpy(local_grid, EDGE, edge_mask, NULL);
  yac_basic_grid_add_coordinates_nocpy(local_grid, EDGE, edge_field_coords);

  struct dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(local_grid, remote_grid, MPI_COMM_WORLD);

  { // test interpolation without using a mask

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid_out =
      yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                          num_src_fields, src_fields, tgt_field);
    struct interp_grid * interp_grid_in =
      yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);
    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);

    yac_interp_method_delete(method_stack_in);
    yac_interp_method_delete(method_stack_out);

    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {1338,1338,1338,1338,1338,
                                       1338,   3,   4,   5,   6,
                                       1338,   8,   9,  10,  11,
                                       1338,  13,  14,  15,  16,
                                       1338,  18,  19,  20,  21};

    for (size_t i = 0; i < num_vertices; ++i) {
      if (grid_data.core_vertex_mask[i]) {
        if (double_are_unequal(
              target_data[0][i],
              ref_global_target_data[grid_data.vertex_ids[i]]))
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid_out =
      yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                          num_src_fields, src_fields, tgt_field);
    struct interp_grid * interp_grid_in =
      yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);
    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);

    yac_interp_method_delete(method_stack_in);
    yac_interp_method_delete(method_stack_out);

    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {1338,1338,1338,1338,1338,
                                       1338,   3,1338,1338,1338,
                                       1338,   8,1338,1338,1338,
                                       1338,  13,1338,1338,1338,
                                       1338,  18,  19,  20,  21};

    for (size_t i = 0; i < num_vertices; ++i) {
      if (grid_data.core_vertex_mask[i]) {
        if (double_are_unequal(
              target_data[0][i],
              ref_global_target_data[grid_data.vertex_ids[i]]))
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask and allowing partial coverage

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid_out =
      yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                          num_src_fields, src_fields, tgt_field);
    struct interp_grid * interp_grid_in =
      yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);
    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);

    yac_interp_method_delete(method_stack_in);
    yac_interp_method_delete(method_stack_out);

    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {1338,1338,    1338,1338,1338,
                                       1338,   3,       3, 2.5, 3.5,
                                       1338,   8,     8.5,1338,1338,
                                       1338,  13,44.0/3.0,17.5,18.5,
                                       1338,  18,      19,  20,  21};

    for (size_t i = 0; i < num_vertices; ++i) {
      if (grid_data.core_vertex_mask[i]) {
        if (fabs(target_data[0][i] -
                 ref_global_target_data[grid_data.vertex_ids[i]]) > 1e-9)
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask and allowing partial coverage and
    // use a target mask

    struct interp_grid * interp_grid_in, * interp_grid_out;

    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

      interp_grid_out =
        yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }

    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

      interp_grid_in =
        yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }

    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);
    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);

    yac_interp_method_delete(method_stack_in);
    yac_interp_method_delete(method_stack_out);

    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {1338,1338,1338,1338,1338,
                                       1338,   3,   3, 2.5, 3.5,
                                         -1,  -1,  -1,1338,1338,
                                         -1,  -1,  -1,17.5,18.5,
                                       1338,  18,  19,  20,  21};

    for (size_t i = 0; i < num_vertices; ++i) {
      if ((grid_data.core_vertex_mask[i]) && corner_mask[i]) {
        if (fabs(target_data[0][i] -
                 ref_global_target_data[grid_data.vertex_ids[i]]) > 1e-9)
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask and allowing partial coverage and
    // use a target mask (target points on edges)

    struct interp_grid * interp_grid_out, * interp_grid_in;
    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CELL, .coordinates_idx = 0, .masks_idx = 0};
      interp_grid_out =
        yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }
    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = EDGE, .coordinates_idx = 0, .masks_idx = 0};
      interp_grid_in =
        yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }

    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);
    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);

    yac_interp_method_delete(method_stack_in);
    yac_interp_method_delete(method_stack_out);

    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_edges * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_edges; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {
        1338,1338,1338,1338,    1338,    1338,1338,1338,1338,
        1338,1338,   3,   3,       3,       3, 2.5, 2.5, 3.5,
          -1,  -1,  -1,  -1,     8.5,      -1,1338,1338,1338,
          -1,1338,  -1,  13,44.0/3.0,44.0/3.0,17.5,17.5,18.5,
        1338,  18,  19,  20};

    for (size_t i = 0; i < num_edges; ++i) {
      if ((grid_data.core_edge_mask[i]) && edge_mask[i]) {
        if (fabs(target_data[0][i] -
                 ref_global_target_data[grid_data.edge_ids[i]]) > 1e-9)
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  //--------
  // cleanup
  //--------

  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(remote_grid);
  yac_basic_grid_delete(local_grid);
}

/*
 * The grid is distributed among 2 processes.
 *
 * The global grid has 4x4 cells:
 *
 *    20--36--21--37--22--38--23--39--24
 *    |       |       |       |       |
 *    28  12  30  13  32  14  34  15  35
 *    |       |       |       |       |
 *    15--27--16--29--17--31--18--33--19
 *    |       |       |       |       |
 *    19  08  21  09  23  10  25  11  26
 *    |       |       |       |       |
 *    10--18--11--20--12--22--13--24--14
 *    |       |       |       |       |
 *    10  04  12  05  14  06  16  07  17
 *    |       |       |       |       |
 *    05--09--06--11--07--13--08--15--09
 *    |       |       |       |       |
 *    01  00  03  01  05  02  07  03  08
 *    |       |       |       |       |
 *    00--00--01--02--02--04--03--06--04
 *
 * The mask looks as follows (# = masked points)
 *
 *    +-------+-------+-------+-------+
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    +-------+-------+-------+-------+
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    +-------+-------#---#---#---#---#
 *    |       |       |       |       |
 *    |       |       #   #   #   #   #
 *    |       |       |       |       |
 *    +-------+-------#---#---#---#---#
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    |       |       |       |       |
 *    +-------+-------+-------+-------+
 */
static void submain_1(MPI_Comm global_comm, MPI_Comm comp_comm,
                      enum interp_weights_reorder_type reorder_type) {

  char const local_grid_name[] = "grid1";
  char const remote_grid_name[] = "grid2";

  int my_rank;
  yac_mpi_call(MPI_Comm_rank(comp_comm, &my_rank), comp_comm);

  double vertex_coordinates_x[] = {0.5, 1.5, 2.5, 3.5, 4.5};
  double vertex_coordinates_y[] = {0.5, 1.5, 2.5, 3.5, 4.5};
  double cell_coordinates_x[] = {0.75, 1.75, 2.75, 3.75};
  double cell_coordinates_y[] = {0.75, 1.75, 2.75, 3.75};
  double cell_coords[16][3];
  size_t const num_global_cells[2] = {4,4};
  size_t local_start[2][2] = {{0,0},{0,2}};
  size_t local_count[2][2] = {{4,2},{4,2}};
  int with_halo = 1;
  for (size_t i = 0; i <= num_global_cells[0]; ++i)
    vertex_coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= num_global_cells[1]; ++i)
    vertex_coordinates_y[i] *= YAC_RAD;
  for (size_t i = 0, k = 0; i < num_global_cells[1]; ++i)
    for (size_t j = 0; j < num_global_cells[0]; ++j, ++k)
      LLtoXYZ_deg(cell_coordinates_x[j], cell_coordinates_y[i], cell_coords[k]);

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      vertex_coordinates_x, vertex_coordinates_y, num_global_cells,
      local_start[my_rank], local_count[my_rank], with_halo);
  size_t num_vertices = grid_data.num_vertices;
  size_t num_cells = grid_data.num_cells;

  struct yac_basic_grid * local_grid =
    yac_basic_grid_new(local_grid_name, grid_data);
  struct yac_basic_grid * remote_grid =
    yac_basic_grid_empty_new(remote_grid_name);

  yac_int masked_corner_ids[] = {7,8,9,12,13,14};
  size_t num_masked_corners =
    sizeof(masked_corner_ids)/sizeof(masked_corner_ids[0]);

  int * corner_mask = xmalloc(num_vertices * sizeof(*corner_mask));
  for (size_t i = 0; i < num_vertices; ++i) {
    corner_mask[i] = 1;
    for (size_t j = 0; j < num_masked_corners; ++j)
      if (grid_data.vertex_ids[i] == masked_corner_ids[j]) corner_mask[i] = 0;
  }

  yac_int masked_cell_ids[] = {6,7};
  size_t num_masked_cells =
    sizeof(masked_cell_ids)/sizeof(masked_cell_ids[0]);

  int * cell_mask = xmalloc(num_cells * sizeof(*cell_mask));
  for (size_t i = 0; i < num_cells; ++i) {
    cell_mask[i] = 1;
    for (size_t j = 0; j < num_masked_cells; ++j)
      if (grid_data.cell_ids[i] == masked_cell_ids[j]) cell_mask[i] = 0;
  }

  coordinate_pointer cell_field_coords =
    xmalloc(num_cells * sizeof(*cell_field_coords));
  for (size_t i = 0; i < num_cells; ++i)
    memcpy(cell_field_coords[i], cell_coords[grid_data.cell_ids[i]],
           3 * sizeof(double));

  yac_basic_grid_add_mask_nocpy(local_grid, CORNER, corner_mask, NULL);
  yac_basic_grid_add_mask_nocpy(local_grid, CELL, cell_mask, NULL);
  yac_basic_grid_add_coordinates_nocpy(local_grid, CELL, cell_field_coords);

  struct dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(local_grid, remote_grid, MPI_COMM_WORLD);

  { // test interpolation without using a mask

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid_in =
      yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                          num_src_fields, src_fields, tgt_field);
    struct interp_grid * interp_grid_out =
      yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);
    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);

    yac_interp_method_delete(method_stack_out);
    yac_interp_method_delete(method_stack_in);

    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {   3,   4,   5,   6,1337,
                                          8,   9,  10,  11,1337,
                                         13,  14,  15,  16,1337,
                                         18,  19,  20,  21,1337,
                                       1337,1337,1337,1337,1337};

    for (size_t i = 0; i < num_vertices; ++i) {
      if (grid_data.core_vertex_mask[i]) {
        if (double_are_unequal(
              target_data[0][i],
              ref_global_target_data[grid_data.vertex_ids[i]]))
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid_in =
      yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                          num_src_fields, src_fields, tgt_field);
    struct interp_grid * interp_grid_out =
      yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);
    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);

    yac_interp_method_delete(method_stack_out);
    yac_interp_method_delete(method_stack_in);

    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {   3,   4,   5,   6,1337,
                                       1337,1337,1337,  11,1337,
                                       1337,1337,1337,  16,1337,
                                       1337,1337,1337,  21,1337,
                                       1337,1337,1337,1337,1337};

    for (size_t i = 0; i < num_vertices; ++i) {
      if (grid_data.core_vertex_mask[i]) {
        if (double_are_unequal(
              target_data[0][i],
              ref_global_target_data[grid_data.vertex_ids[i]]))
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask and allowing partial coverage

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid_in =
      yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                          num_src_fields, src_fields, tgt_field);
    struct interp_grid * interp_grid_out =
      yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);
    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);

    yac_interp_method_delete(method_stack_out);
    yac_interp_method_delete(method_stack_in);

    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {   3,   4,       5,   6,1337,
                                        5.5, 6.5,28.0/3.0,  11,1337,
                                       1337,1337,    15.5,  16,1337,
                                       20.5,21.5,      21,  21,1337,
                                       1337,1337,    1337,1337,1337};

    for (size_t i = 0; i < num_vertices; ++i) {
      if (grid_data.core_vertex_mask[i]) {
        if (fabs(target_data[0][i] -
                 ref_global_target_data[grid_data.vertex_ids[i]]) > 1e-9)
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask and allowing partial coverage and
    // use a target mask

    struct interp_grid * interp_grid_in, * interp_grid_out;

    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

      interp_grid_in =
        yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }

    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

      interp_grid_out =
        yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }

    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);
    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);

    yac_interp_method_delete(method_stack_out);
    yac_interp_method_delete(method_stack_in);

    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_vertices * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_vertices; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {   3,   4,   5,   6,1337,
                                        5.5, 6.5,  -1,  -1,  -1,
                                       1337,1337,  -1,  -1,  -1,
                                       20.5,21.5,  21,  21,1337,
                                       1337,1337,1337,1337,1337};

    for (size_t i = 0; i < num_vertices; ++i) {
      if ((grid_data.core_vertex_mask[i]) && corner_mask[i]) {
        if (fabs(target_data[0][i] -
                 ref_global_target_data[grid_data.vertex_ids[i]]) > 1e-9)
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  { // test interpolation using a source mask and allowing partial coverage and
    // use a target mask (target points in cell centers)

    struct interp_grid * interp_grid_in, * interp_grid_out;
    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CELL, .coordinates_idx = 0, .masks_idx = 0};
      interp_grid_in =
        yac_interp_grid_new(grid_pair, remote_grid_name, local_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }
    {
      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = EDGE, .coordinates_idx = 0, .masks_idx = 0};
      interp_grid_out =
        yac_interp_grid_new(grid_pair, local_grid_name, remote_grid_name,
                            num_src_fields, src_fields, tgt_field);
    }

    struct interp_method * method_stack_in[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337),
       NULL};
    struct interp_method * method_stack_out[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 1),
       yac_interp_method_fixed_new(1337 + 1),
       NULL};

    struct interp_weights * weights_in =
      yac_interp_method_do_search(method_stack_in, interp_grid_in);
    struct interp_weights * weights_out =
      yac_interp_method_do_search(method_stack_out, interp_grid_out);

    yac_interp_method_delete(method_stack_out);
    yac_interp_method_delete(method_stack_in);

    struct interpolation * interpolation_in =
      yac_interp_weights_get_interpolation(
        weights_in, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct interpolation * interpolation_out =
      yac_interp_weights_get_interpolation(
        weights_out, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

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
    // target_data dimensions [collection_idx]
    //                        [local_idx]
    double * target_data_field =
      xmalloc(num_cells * sizeof(*target_data_field));
    double * target_data[1] = {target_data_field}; // collection_size == 1
    for (size_t i = 0; i < num_cells; ++i) target_data_field[i] = -1.0;

    yac_interpolation_execute_put(interpolation_out, source_data);
    yac_interpolation_execute_get(interpolation_in, target_data);

    //----------------------------
    // check interpolation results
    //----------------------------

    double ref_global_target_data[] = {   3,   4,   5,   6,
                                        5.5, 6.5,  -1,  -1,
                                       1337,1337,15.5,  16,
                                       20.5,21.5,  21,  21};

    for (size_t i = 0; i < num_cells; ++i) {
      if ((grid_data.core_cell_mask[i]) && cell_mask[i]) {
        if (fabs(target_data[0][i] -
                 ref_global_target_data[grid_data.cell_ids[i]]) > 1e-9)
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
    free(source_data_field);
    yac_interpolation_delete(interpolation_out);
    yac_interpolation_delete(interpolation_in);
    yac_interp_weights_delete(weights_out);
    yac_interp_weights_delete(weights_in);
    yac_interp_grid_delete(interp_grid_out);
    yac_interp_grid_delete(interp_grid_in);
  }

  //---------
  // clean up
  //---------

  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(remote_grid);
  yac_basic_grid_delete(local_grid);
}

