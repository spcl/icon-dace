/**
 * @file test_dummy_coupling6_c.c
 *
 * @copyright Copyright  (C)  2022 DKRZ, MPI-M
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

static void run_comp_a();
static void run_comp_b();

int main (void) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
      run_comp_a();
      break;
    case(2):
      run_comp_b();
      break;
    default:
      PUT_ERR("wrong number of processes (has to be 3)");
      break;
  }

  xt_finalize();
  MPI_Finalize();
  return TEST_EXIT_CODE;
}

static void run_comp_a() {

  // initialise YAC default instance
  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_json("coupling_test6.json");

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
  int cell_mask[2] = {1, 1};
  int vertex_mask[4] = {1, 1, 1, 1};
  yac_cdef_mask(
    grid_id, nbr_cells, YAC_LOCATION_CELL, cell_mask, &cell_mask_id);
  yac_cdef_mask(
    grid_id, nbr_vertices, YAC_LOCATION_CORNER, vertex_mask, &vertex_mask_id);

  // define field
  int field_id, multi_field_id;
  int collection_size = 1;
  yac_cdef_field(
    "A_to_B_src", comp_id, &cell_point_id, 1, collection_size, "1",
    YAC_TIME_UNIT_SECOND, &field_id);
  int multi_field_point_ids[2], multi_field_mask_ids[2];
  multi_field_point_ids[0] = cell_point_id;
  multi_field_point_ids[1] = vertex_point_id;
  multi_field_mask_ids[0] = cell_mask_id;
  multi_field_mask_ids[1] = vertex_mask_id;
  yac_cdef_field_mask(
    "A_to_B_multi", comp_id, multi_field_point_ids, multi_field_mask_ids,
    2, collection_size, "1", YAC_TIME_UNIT_SECOND, &multi_field_id);

  // generate coupling
  int ierror;
  yac_cenddef();

  // move data from comp_A to comp_B
  {
    double send_field_data[2] = {3.0, 4.0};
    double * send_field_[1] = {&send_field_data[0]};
    double ** send_field[1] = {&send_field_[0]};
    int info;
    yac_cput(field_id, collection_size, send_field, &info, &ierror);
  }
  {
    double send_cell_field_data[2] = {3.0, 4.0};
    double send_vertex_field_data[4] = {3.0, 4.0, 5.0, 6.0};
    double * send_field_[2] =
      {&send_cell_field_data[0], &send_vertex_field_data[0]};
    double ** send_field[1] = {&send_field_[0]};
    int info;
    yac_cput(multi_field_id, collection_size, send_field, &info, &ierror);
  }

  // finalise YAC default instance
  yac_cfinalize();
}

static void run_comp_b() {

  // initialise YAC default instance
  int default_instance_id;
  yac_cinit_instance(&default_instance_id);
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_json_instance(
    default_instance_id, "coupling_test6.json");

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
  int field_id, multi_field_id;
  int collection_size = 1;
  yac_cdef_field(
    "A_to_B_tgt", comp_id, &point_id_B, 1, collection_size, "1",
    YAC_TIME_UNIT_SECOND, &field_id);
  yac_cdef_field(
    "A_to_B_multi", comp_id, &point_id_B, 1, collection_size, "1",
    YAC_TIME_UNIT_SECOND, &multi_field_id);

  // generate coupling
  int ierror;
  yac_cenddef_instance(default_instance_id);

  // get data from comp_A
  double recv_field_data[2] = {-1.0, -1.0};
  double * recv_field[1] = {&recv_field_data[0]};
  int info;
  yac_cget(multi_field_id, collection_size, recv_field, &info, &ierror);
  if ((recv_field_data[0] != -2.0) || (recv_field_data[1] != -2.0))
    PUT_ERR("ERROR: wrong results from multi field exchange");
  yac_cget(field_id, collection_size, recv_field, &info, &ierror);

  { // run internal instance

    MPI_Comm internal_comm;
    yac_cget_comp_comm(comp_id, &internal_comm);

    // initialise internal YAC instance
    int internal_instance_id;
    yac_cinit_comm_instance(internal_comm, &internal_instance_id);
    yac_cread_config_yaml_instance(internal_instance_id,
      "coupling_test6_local.yaml");

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
      "B_to_C", internal_comp_id, &point_id_B, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_id_src);
    yac_cdef_field(
      "B_to_C", internal_comp_id, &point_id_C, 1, 1, "1",
      YAC_TIME_UNIT_SECOND, &field_id_tgt);

    // generate internal coupling
    int ierror;
    yac_cenddef_instance(internal_instance_id);

    // move data from grid B to grid C
    int send_info, recv_info;
    int collection_size = 1;
    double * send_field_[1] = {&recv_field_data[0]};
    double ** send_field[1] = {&send_field_[0]};
    double recv_field_data_C[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    double * recv_field[1] = {&recv_field_data_C[0]};
    yac_cexchange(
      field_id_src, field_id_tgt, collection_size, send_field, recv_field,
      &send_info, &recv_info, &ierror);

    double ref_recv_field_data_C[6] = {4.0, 4.0, 4.0, 4.0, 3.0, 3.0};
    for (int i = 0; i < 6; ++i)
      if (recv_field_data_C[i] != ref_recv_field_data_C[i])
        PUT_ERR("data missmatch");

    yac_ccleanup_instance(internal_instance_id);

    MPI_Comm_free(&internal_comm);
  }

  // finalise YAC default instance
  yac_cfinalize_instance(default_instance_id);
}
