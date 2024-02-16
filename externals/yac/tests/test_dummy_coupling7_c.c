/**
 * @file test_dummy_coupling7_c.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
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

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <yaxt.h>
#include "tests.h"
#include "yac_interface.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

static void run_comp_a(int with_field_mask, int exchange_type);
static void run_comp_b(int with_field_mask, int exchange_type);

int main (void) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (int exchange_type = 0; exchange_type < 2; ++exchange_type) {
    for (int with_field_mask = 0; with_field_mask < 2; ++with_field_mask) {
      switch(rank) {
        case(0): {
          MPI_Comm yac_comm;
          const char* group_name = "yac";
          yac_cmpi_handshake(MPI_COMM_WORLD, 1, &group_name, &yac_comm);
          yac_cinit_comm_dummy(yac_comm);
          MPI_Comm_free(&yac_comm);
          break;
        }
        case(1):
          run_comp_a(with_field_mask, exchange_type);
          break;
        case(2):
          run_comp_b(with_field_mask, exchange_type);
          break;
        default:
          PUT_ERR("wrong number of processes (has to be 3)");
          break;
      }
    }
  }

  xt_finalize();
  MPI_Finalize();
  return TEST_EXIT_CODE;
}

static void run_comp_a(int with_field_mask, int exchange_type) {

  // initialise YAC default instance
  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime("2000-01-01T00:00:00", "2000-01-01T00:00:12");

  // define local component
  int comp_id;
  yac_cdef_comp("comp_A", &comp_id);

  // get communicator that contains both components
  MPI_Comm pair_comm;
  yac_cget_comps_comm((char const*[]){"comp_A", "comp_B"}, 2, &pair_comm);

  // check the pair_comm
  {
    char * sendbuf = "A";
    char recvbuf[2];
    MPI_Allgather(sendbuf, 1, MPI_CHAR, recvbuf, 1, MPI_CHAR, pair_comm);

    if ((recvbuf[0] != 'A') || (recvbuf[1] != 'B'))
      PUT_ERR("ERROR in yac_cget_comps_instance");
  }
  MPI_Comm_free(&pair_comm);

  int nbr_vertices = 4;
  int nbr_cells = 2;
  int nbr_vertices_per_cell[2] = {3,3};
  double x_vertices[4] =
    {0.0 * YAC_RAD, -1.0 * YAC_RAD, 1.0 * YAC_RAD, 0.0 * YAC_RAD};
  double y_vertices[4] =
    {1.0 * YAC_RAD, 0.0 * YAC_RAD, 0.0 * YAC_RAD, -1.0 * YAC_RAD};
  int cell_to_vertex[6] = {0,1,2, 1,3,2};

  /* define local grid

      0
     / \
    / o \
   /     \
  1-------2   Eq.
   \     /
    \ o /
     \ /
      3

  */
  int grid_id;
  yac_cdef_grid_unstruct(
    "grid_A", nbr_vertices, nbr_cells, nbr_vertices_per_cell,
    x_vertices, y_vertices, cell_to_vertex, &grid_id);

  double x_cells[2] = {0.0 * YAC_RAD, 0.0 * YAC_RAD};
  double y_cells[2] = {0.5 * YAC_RAD, -0.5 * YAC_RAD};

  // center points in cells (needed e.g. for nearest neighbour)
  int cell_point_id, vertex_point_id;
  yac_cdef_points_unstruct(
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, &cell_point_id);
  // vertex points
  yac_cdef_points_unstruct(
    grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices,
    &vertex_point_id);

  // masks
  int cell_mask_id, vertex_mask_id;
  int cell_mask[2] = {1, 0};
  int vertex_mask[4] = {1, 1, 1, 1};
  yac_cdef_mask(
    grid_id, nbr_cells, YAC_LOCATION_CELL, cell_mask, &cell_mask_id);
  yac_cdef_mask(
    grid_id, nbr_vertices, YAC_LOCATION_CORNER, vertex_mask, &vertex_mask_id);

  // define field
  int field_id, multi_field_id, dummy_field_id;
  double frac_mask_fallback_value = -1.0;
  double frac_mask_fallback_value_multi = strtod("inf", NULL);
  int collection_size = 1;
  if (with_field_mask) {
    yac_cdef_field_mask(
      "A_to_B_src", comp_id, &cell_point_id, &cell_mask_id, 1,
      collection_size, "1", YAC_TIME_UNIT_SECOND, &field_id);
  } else {
    yac_cdef_field(
      "A_to_B_src", comp_id, &cell_point_id, 1,
      collection_size, "1", YAC_TIME_UNIT_SECOND, &field_id);
  }
  yac_cenable_field_frac_mask(
    "comp_A", "grid_A", "A_to_B_src", frac_mask_fallback_value);
  int multi_field_point_ids[2], multi_field_mask_ids[2];
  multi_field_point_ids[0] = cell_point_id;
  multi_field_point_ids[1] = vertex_point_id;
  multi_field_mask_ids[0] = cell_mask_id;
  multi_field_mask_ids[1] = vertex_mask_id;
  if (with_field_mask) {
    yac_cdef_field_mask(
      "A_to_B_multi", comp_id, multi_field_point_ids, multi_field_mask_ids,
      2, collection_size, "1", YAC_TIME_UNIT_SECOND, &multi_field_id);
  } else {
    yac_cdef_field(
      "A_to_B_multi", comp_id, multi_field_point_ids,
      2, collection_size, "1", YAC_TIME_UNIT_SECOND, &multi_field_id);
  }
  yac_cenable_field_frac_mask(
    "comp_A", "grid_A", "A_to_B_multi", frac_mask_fallback_value_multi);
  yac_cdef_field(
    "dummy_field", comp_id, &cell_point_id, 1,
    collection_size, "1", YAC_TIME_UNIT_SECOND, &dummy_field_id);

  // generate coupling
  int ierror;
  yac_cenddef();

  if (yac_cget_field_frac_mask_fallback_value(
        "comp_A", "grid_A", "A_to_B_multi") != frac_mask_fallback_value_multi)
    PUT_ERR("ERROR in yac_cget_field_frac_mask_fallback_value");

  // move data from comp_A to comp_B
  {
    double send_field_data[2] = {3.0, 4.0};
    double * send_field_[1] = {&send_field_data[0]};
    double ** send_field[1] = {&send_field_[0]};
    double send_frac_mask_data[2] = {1.0, 1.0};
    double * send_frac_mask_[1] = {&send_frac_mask_data[0]};
    double ** send_frac_mask[1] = {&send_frac_mask_[0]};
    int info;
    if (exchange_type) {
      yac_cput_frac(
        field_id, collection_size, send_field, send_frac_mask, &info, &ierror);
      yac_cget(dummy_field_id, collection_size, NULL, &info, &ierror);
    } else {
      yac_cexchange_frac(
        field_id, dummy_field_id, collection_size,
        send_field, send_frac_mask, NULL, &info, &info, &ierror);
    }
  }
  {
    double send_cell_field_data[2] = {3.0, 4.0};
    double send_vertex_field_data[4] = {3.0, 4.0, 5.0, 6.0};
    double * send_field_[2] =
      {&send_cell_field_data[0], &send_vertex_field_data[0]};
    double ** send_field[1] = {&send_field_[0]};
    double send_cell_frac_mask_data[2] = {1.0, 1.0};
    double send_vertex_frac_mask_data[4] = {1.0, 1.0, 1.0, 1.0};
    double * send_frac_mask_[2] =
      {&send_cell_frac_mask_data[0], &send_vertex_frac_mask_data[0]};
    double ** send_frac_mask[1] = {&send_frac_mask_[0]};
    int info;
    yac_cput_frac(
      multi_field_id, collection_size, send_field, send_frac_mask,
      &info, &ierror);
  }

  // finalise YAC default instance
  yac_cfinalize();
}

static void run_comp_b(int with_field_mask, int exchange_type) {

  double frac_mask_fallback_value = strtod("inf", NULL);

  // initialise YAC default instance
  int default_instance_id;
  yac_cinit_instance(&default_instance_id);
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime_instance(
    default_instance_id, "2000-01-01T00:00:00", "2000-01-01T00:00:12");

  // define local component
  int comp_id;
  yac_cdef_comp_instance(default_instance_id, "comp_B", &comp_id);

  // get communicator that contains both components
  MPI_Comm pair_comm;
  yac_cget_comps_comm_instance(
    default_instance_id, (char const *[]){"comp_B", "comp_A"}, 2, &pair_comm);

  // check the pair_comm
  {
    char * sendbuf = "B";
    char recvbuf[2];
    MPI_Allgather(sendbuf, 1, MPI_CHAR, recvbuf, 1, MPI_CHAR, pair_comm);
    if ((recvbuf[0] != 'A') || (recvbuf[1] != 'B'))
      PUT_ERR("ERROR in yac_cget_comps_comm");
  }
  MPI_Comm_free(&pair_comm);

  int nbr_vertices[2] = {2,3};
  int cyclic[2] = {0,0};
  double x_vertices[2] = {-0.5 * YAC_RAD, 0.5 * YAC_RAD};
  double y_vertices[3] = {-1.0 * YAC_RAD, 0.0 * YAC_RAD, 1.0 * YAC_RAD};

  /* define local grid

  4-------5
  |       |
  |   o   |
  |       |
  2-------3
  |       |
  |   o   |
  |       |
  0-------1

  */
  int grid_id_B;
  yac_cdef_grid_reg2d(
    "grid_B", nbr_vertices, cyclic, x_vertices, y_vertices, &grid_id_B);

  int nbr_cells[2] = {1,2};
  double x_cells[1] = {0.0 * YAC_RAD};
  double y_cells[2] = {-0.5 * YAC_RAD, 0.5 * YAC_RAD};

  // center points in cells (needed e.g. for nearest neighbour)
  int point_id_B;
  yac_cdef_points_reg2d(
    grid_id_B, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, &point_id_B);

  // define field
  int field_id, multi_field_id, dummy_field_id;
  int collection_size = 1;
  yac_cdef_field(
    "A_to_B_tgt", comp_id, &point_id_B, 1, collection_size, "1",
    YAC_TIME_UNIT_SECOND, &field_id);
  yac_cdef_field(
    "A_to_B_multi", comp_id, &point_id_B, 1, collection_size, "1",
    YAC_TIME_UNIT_SECOND, &multi_field_id);
  yac_cdef_field(
    "dummy_field", comp_id, &point_id_B, 1, collection_size, "1",
    YAC_TIME_UNIT_SECOND, &dummy_field_id);

  // define interpolation stacks
  int interp_stack_nnn, interp_stack_fixed;
  yac_cget_interp_stack_config(&interp_stack_nnn);
  yac_cadd_interp_stack_config_nnn(
    interp_stack_nnn, YAC_NNN_AVG, 1, 0.0);
  yac_cget_interp_stack_config(&interp_stack_fixed);
  yac_cadd_interp_stack_config_fixed(interp_stack_fixed, -2.0);

  // define couplings
  yac_cdef_couple_instance(
    default_instance_id,
    "comp_A", "grid_A", "A_to_B_src",
    "comp_B", "grid_B", "A_to_B_tgt",
    "PT01.000S", YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
    interp_stack_nnn, 0, 0);
  yac_cdef_couple_instance(
    default_instance_id,
    "comp_A", "grid_A", "A_to_B_multi",
    "comp_B", "grid_B", "A_to_B_multi",
    "PT01.000S", YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
    interp_stack_fixed, 0, 0);
  yac_cfree_interp_stack_config(interp_stack_nnn);
  yac_cfree_interp_stack_config(interp_stack_fixed);

  // generate coupling
  int ierror;
  yac_cenddef_instance(default_instance_id);

  if (yac_cget_field_frac_mask_fallback_value_instance(
        default_instance_id, "comp_A", "grid_A", "A_to_B_multi") !=
      frac_mask_fallback_value)
    PUT_ERR("ERROR in yac_cget_field_frac_mask_fallback_value_instance");

  // get data from comp_A
  double recv_field_data[2] = {-1.0, -1.0};
  double * recv_field[1] = {&recv_field_data[0]};
  int info;
  yac_cget(multi_field_id, collection_size, recv_field, &info, &ierror);
  if ((recv_field_data[0] != -2.0) || (recv_field_data[1] != -2.0))
    PUT_ERR("ERROR: wrong results from multi field exchange");
  if (exchange_type) {
    yac_cput(
      dummy_field_id, collection_size, NULL, &info, &ierror);
    yac_cexchange(
      dummy_field_id, field_id, collection_size, NULL, recv_field,
      &info, &info, &ierror);
  } else {
    yac_cget(field_id, collection_size, recv_field, &info, &ierror);
  }

  { // run internal instance

    MPI_Comm internal_comm;
    yac_cget_comp_comm(comp_id, &internal_comm);

    // initialise internal YAC instance
    int internal_instance_id;
    yac_cinit_comm_instance(internal_comm, &internal_instance_id);
    yac_cdef_datetime_instance(
      internal_instance_id, "2000-01-01T00:00:00", "2000-01-01T00:00:12");

    // define internal component
    int internal_comp_id;
    yac_cdef_comp_instance(internal_instance_id, "comp_B", &internal_comp_id);

    int nbr_vertices[2] = {2,3};
    int cyclic[2] = {0,0};
    double x_vertices[2*3] = {-0.4 * YAC_RAD, 0.4 * YAC_RAD,
                              -0.4 * YAC_RAD, 0.4 * YAC_RAD,
                              -0.4 * YAC_RAD, 0.4 * YAC_RAD};
    double y_vertices[2*3] = {-0.7 * YAC_RAD, -0.7 * YAC_RAD,
                               0.0 * YAC_RAD,  0.0 * YAC_RAD,
                               0.8 * YAC_RAD,  0.8 * YAC_RAD};

    /* define local grid

    4-------5
    |       |
    |       |
    |       |
    2-------3
    |       |
    |       |
    |       |
    0-------1

    */
    int grid_id_C;
    yac_cdef_grid_curve2d(
      "grid_C", nbr_vertices, cyclic, x_vertices, y_vertices, &grid_id_C);

    // vertices points (needed e.g. for nearest neighbour)
    int point_id_C;
    yac_cdef_points_curve2d(
      grid_id_C, nbr_vertices, YAC_LOCATION_CORNER,
      x_vertices, y_vertices, &point_id_C);

    // define field
    int field_id_src, field_id_tgt;
    yac_cdef_field(
      "B_to_C", internal_comp_id, &point_id_B, 1,
      1, "1", YAC_TIME_UNIT_SECOND, &field_id_src);
    yac_cenable_field_frac_mask_instance(
      internal_instance_id, "comp_B", "grid_B", "B_to_C",
      frac_mask_fallback_value);
    yac_cdef_field(
      "B_to_C", internal_comp_id, &point_id_C, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_id_tgt);

    // define interpolation stack
    int interp_stack_nnn;
    yac_cget_interp_stack_config(&interp_stack_nnn);
    yac_cadd_interp_stack_config_nnn(
      interp_stack_nnn, YAC_NNN_AVG, 1, 0.0);

    // define coupling
    int ext_couple_config;
    double const scale_factor = 2.0;
    double const scale_summand = 0.5;
    yac_cget_ext_couple_config(&ext_couple_config);
    yac_cset_ext_couple_config_scale_factor(
      ext_couple_config, scale_factor);
    yac_cset_ext_couple_config_scale_summand(
      ext_couple_config, scale_summand);
    yac_cdef_couple_custom_instance(
      internal_instance_id,
      "comp_B", "grid_B", "B_to_C",
      "comp_B", "grid_C", "B_to_C",
      "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
      interp_stack_nnn, 0, 0, ext_couple_config);
    yac_cfree_ext_couple_config(ext_couple_config);
    yac_cfree_interp_stack_config(interp_stack_nnn);

    // generate internal coupling
    int ierror;
    yac_cenddef_instance(internal_instance_id);

    // move data from grid B to grid C
    int send_info, recv_info;
    int collection_size = 1;
    double * send_field_[1] = {&recv_field_data[0]};
    double ** send_field[1] = {&send_field_[0]};
    double send_frac_mask_data[2] = {1.0, 0.0};
    double * send_frac_mask_[1] = {send_frac_mask_data};
    double ** send_frac_mask[1] = {&send_frac_mask_[0]};
    double recv_field_data_C[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    double * recv_field[1] = {&recv_field_data_C[0]};

    if (exchange_type)
      yac_cexchange_frac(
        field_id_src, field_id_tgt, collection_size,
        send_field, send_frac_mask, recv_field,
        &send_info, &recv_info, &ierror);
    else
      yac_cexchange_frac_ptr_(
        field_id_src, field_id_tgt, collection_size,
        send_field_, send_frac_mask_, recv_field,
        &send_info, &recv_info, &ierror);

    double ref_recv_field_data_C[2][6] =
      {{4.0, 4.0, 4.0, 4.0,
        frac_mask_fallback_value, frac_mask_fallback_value},
      {3.0, 3.0, 3.0, 3.0,
        frac_mask_fallback_value, frac_mask_fallback_value}};
    for (int i = 0; i < 6; ++i)
      if (recv_field_data_C[i] !=
          scale_factor * ref_recv_field_data_C[with_field_mask][i] +
          scale_summand)
        PUT_ERR("data missmatch");

    yac_ccleanup_instance(internal_instance_id);

    MPI_Comm_free(&internal_comm);
  }

  // finalise YAC default instance
  yac_cfinalize_instance(default_instance_id);
}
