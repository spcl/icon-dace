/**
 * @file test_dummy_coupling4_c.c
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

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "tests.h"
#include "test_common.h"
#include "yac_interface.h"

enum {
  NUM_FIELDS = 5,
  COLLECTION_SIZE = 3,
  NUM_POINTSETS = 1,
  NUM_POINTS = 9
};

static void init_ref_recv_field(
  double ref_recv_field[][COLLECTION_SIZE][NUM_POINTS]) {

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < COLLECTION_SIZE; ++j)
      for (int k = 0; k < NUM_POINTS; ++k)
        ref_recv_field[i][j][k] = 0.0;
  for (int j = 0; j < COLLECTION_SIZE; ++j) {
    for (int k = 0; k < NUM_POINTS; ++k) {
      ref_recv_field[3][j][k] = DBL_MAX;
      ref_recv_field[4][j][k] = -DBL_MAX;
    }
  }
}

int main (void) {

  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_yaml("coupling_test4.yaml");

  int size, rank;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  if (size != 2) {
    fputs("wrong number of processes (has to be 2)\n", stderr);
    exit(EXIT_FAILURE);
  }

  int is_target = rank == 1;

  // define local component
  int comp_id;
  yac_cdef_comp((is_target)?"target_comp":"source_comp", &comp_id);

  // define grid (both components use an identical grid
  int grid_id;
  yac_cdef_grid_reg2d(
    (is_target)?"target_grid":"source_grid", (int[2]){3,3}, (int[2]){0,0},
    (double[]){0,1,2}, (double[]){0,1,2}, &grid_id);

  // define points at the vertices of the grid
  int point_id;
  yac_cdef_points_reg2d(
    grid_id, (int[2]){3,3}, YAC_LOCATION_CORNER,
    (double[]){0,1,2}, (double[]){0,1,2}, &point_id);

  // define fields
  int field_ids[NUM_FIELDS];
  const char * fieldName[NUM_FIELDS] =
    {"time_op_accu_field",
     "time_op_avg_field",
     "time_op_none_field",
     "time_op_min_field",
     "time_op_max_field"};
  for (int i = 0; i < NUM_FIELDS; ++i )
    yac_cdef_field(
      fieldName[i], comp_id, &point_id, 1, 3,
        (is_target)?"2":"1", YAC_TIME_UNIT_SECOND,
        &field_ids[i]);

  yac_cenddef ( );

  double send_field[COLLECTION_SIZE][NUM_POINTSETS][NUM_POINTS] =
    {{{ 1, 2, 3, 4, 5, 6, 7, 8, 9}},
     {{10,11,12,13,14,15,16,17,18}},
     {{19,20,21,22,23,24,25,26,27}}};
  double ref_recv_field[NUM_FIELDS][COLLECTION_SIZE][NUM_POINTS];
  init_ref_recv_field(ref_recv_field);

  // do time steps
  for (int t = 0; t < 32; ++t) {

    if (!is_target) {

      for (int i = 0; i < NUM_FIELDS; ++i) {

        int info, ierror;
        yac_cput_(
          field_ids[i], COLLECTION_SIZE, &send_field[0][0][0], &info, &ierror);

        int ref_send_info =
          (t%4)?
            ((i == 2)?YAC_ACTION_NONE:YAC_ACTION_REDUCTION):
            YAC_ACTION_COUPLING;
        if (info != ref_send_info) PUT_ERR("error in yac_cput_: wrong info");
        if (ierror) PUT_ERR("error in yac_cput_: wrong ierror");
      }
    }

    double scale = ((t == 0)?1.0:0.25);
    for (int i = 0; i < COLLECTION_SIZE; ++i)
      for (int j = 0; j < NUM_POINTS; ++j) {

        // update ref_recv_field (accu)
        ref_recv_field[0][i][j] += send_field[i][0][j];

        // update ref_recv_field (avg)
        ref_recv_field[1][i][j] += send_field[i][0][j] * scale;

        // update ref_recv_field (none)
        ref_recv_field[2][i][j] = send_field[i][0][j];

        // update ref_recv_field (minimum)
        ref_recv_field[3][i][j] =
          MIN(ref_recv_field[3][i][j], send_field[i][0][j]);

        // update ref_recv_field (maximum)
        ref_recv_field[4][i][j] =
          MAX(ref_recv_field[4][i][j], send_field[i][0][j]);
    }

    if (is_target) {

      // target calls get every second timestep
      if ((t % 2) == 0) {
        for (int i = 0; i < NUM_FIELDS; ++i) {

          // initialise recv_field
          double recv_field[COLLECTION_SIZE][NUM_POINTS];
          for (int j = 0; j < COLLECTION_SIZE; ++j)
            for (int k = 0; k < NUM_POINTS; ++k)
              recv_field[j][k] = -1;

          int info, ierror;
          yac_cget_(
            field_ids[i], COLLECTION_SIZE, &recv_field[0][0], &info, &ierror);

          int ref_recv_info = (t%4)?YAC_ACTION_NONE:YAC_ACTION_COUPLING;
          if (info != ref_recv_info) PUT_ERR("error in yac_cget_: wrong info");
          if (ierror) PUT_ERR("error in yac_cget_: wrong ierror");

          if (info == YAC_ACTION_COUPLING) {
            for (int j = 0; j < COLLECTION_SIZE; ++j)
              for (int k = 0; k < NUM_POINTS; ++k)
                if (fabs(recv_field[j][k] - ref_recv_field[i][j][k]) > 1e-6)
                  PUT_ERR("error in yac_cget_: wrong recv_field");
          } else {
            for (int j = 0; j < COLLECTION_SIZE; ++j)
              for (int k = 0; k < NUM_POINTS; ++k)
                if (recv_field[j][k] != -1)
                  PUT_ERR("error in yac_cget_: wrong recv_field");
          }
        }
      }
    }

    // update send_field
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < NUM_POINTS; ++j)
        send_field[i][0][j] += i - 1;

    // clean ref_recv_field at every coupling timestep
    if ((t % 4) == 0) init_ref_recv_field(ref_recv_field);
  }

  yac_cfinalize ();

  return TEST_EXIT_CODE;
}
