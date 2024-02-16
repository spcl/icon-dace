/**
 * @file test_dummy_coupling8_c.c
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
  NUM_POINTS = 9,
  COUPLING_DT = 4,
  TGT_DT = 2,
  SRC_DT = 1,
};

#define FRAC_MASK_VALUE (1337.0)

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

static void init_temp_frac_mask(
  double temp_frac_mask[][COLLECTION_SIZE][NUM_POINTSETS][NUM_POINTS]) {

  for (int i = 0; i < NUM_FIELDS; ++i)
    for (int j = 0; j < COLLECTION_SIZE; ++j)
      for (int k = 0; k < NUM_POINTS; ++k)
        temp_frac_mask[i][j][0][k] = 0.0;
}

int main (void) {

  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime("1850-01-01T00:00:00", "1850-01-03T00:00:00");

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
  for (int field_idx = 0; field_idx < NUM_FIELDS; ++field_idx) {
    yac_cdef_field(
      fieldName[field_idx], comp_id, &point_id, NUM_POINTSETS,
      COLLECTION_SIZE, (is_target)?"2":"1", YAC_TIME_UNIT_SECOND,
      &field_ids[field_idx]);
    yac_cenable_field_frac_mask(
      (is_target)?"target_comp":"source_comp",
      (is_target)?"target_grid":"source_grid",
      fieldName[field_idx], FRAC_MASK_VALUE);
  }

  // define interpolation stacks
  int interp_stack_nnn;
  yac_cget_interp_stack_config(&interp_stack_nnn);
  yac_cadd_interp_stack_config_nnn(
    interp_stack_nnn, YAC_NNN_AVG, 1, 0.0);

  // define couplings
  int reduction_type[NUM_FIELDS] =
    {YAC_REDUCTION_TIME_ACCUMULATE,
     YAC_REDUCTION_TIME_AVERAGE,
     YAC_REDUCTION_TIME_NONE,
     YAC_REDUCTION_TIME_MINIMUM,
     YAC_REDUCTION_TIME_MAXIMUM};
  for (int field_idx = 0; field_idx < NUM_FIELDS; ++field_idx)
    yac_cdef_couple(
      "source_comp", "source_grid", fieldName[field_idx],
      "target_comp", "target_grid", fieldName[field_idx],
      "4", YAC_TIME_UNIT_SECOND, reduction_type[field_idx],
      interp_stack_nnn, 0, 0);
  yac_cfree_interp_stack_config(interp_stack_nnn);

  yac_cenddef ( );

  double send_field[COLLECTION_SIZE][NUM_POINTSETS][NUM_POINTS] =
    {{{ 1, 2, 3, 4, 5, 6, 7, 8, 9}},
     {{10,11,12,13,14,15,16,17,18}},
     {{19,20,21,22,23,24,25,26,27}}};
  double frac_mask[COUPLING_DT][COLLECTION_SIZE][NUM_POINTSETS][NUM_POINTS] =
    {{{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}},
      {{1.0, 1.0, 0.7, 0.7, 0.5, 0.3, 0.3, 0.0, 0.0}},
      {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0}}},
     {{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0}},
      {{0.7, 0.7, 0.7, 0.7, 0.5, 0.3, 0.3, 0.3, 0.3}},
      {{1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0}}},
     {{{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0}},
      {{0.3, 0.3, 0.3, 0.3, 0.5, 0.7, 0.7, 0.7, 0.7}},
      {{1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}},
     {{{0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}},
      {{0.0, 0.0, 0.3, 0.3, 0.5, 0.7, 0.7, 1.0, 1.0}},
      {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}}};
  double temp_frac_mask[NUM_FIELDS][COLLECTION_SIZE][NUM_POINTSETS][NUM_POINTS];
  double ref_recv_field[NUM_FIELDS][COLLECTION_SIZE][NUM_POINTS];
  init_ref_recv_field(ref_recv_field);
  init_temp_frac_mask(temp_frac_mask);

  // do time steps
  for (int t = 0; t < 8 * COUPLING_DT; ++t) {

    if (!is_target) {

      for (int field_idx = 0; field_idx < NUM_FIELDS; ++field_idx) {

        int info, ierror;
        yac_cput_frac_(
          field_ids[field_idx], COLLECTION_SIZE, &send_field[0][0][0],
          &frac_mask[t%COUPLING_DT][0][0][0], &info, &ierror);

        int ref_send_info =
          (t%COUPLING_DT)?
            ((field_idx == 2)?YAC_ACTION_NONE:YAC_ACTION_REDUCTION):
            YAC_ACTION_COUPLING;
        if (info != ref_send_info) PUT_ERR("error in yac_cput_: wrong info");
        if (ierror) PUT_ERR("error in yac_cput_: wrong ierror");
      }
    }

    // the first timestep is coupled directly
    double scale = ((t == 0)?1.0:0.25);
    for (int i = 0; i < COLLECTION_SIZE; ++i) {
      for (int j = 0; j < NUM_POINTS; ++j) {

        double frac_mask_value = frac_mask[t%COUPLING_DT][i][0][j];
        double frac_send_field_value =
          frac_mask_value * send_field[i][0][j];

        if (frac_mask_value != 0.0) {

          // update ref_recv_field (accu)
          ref_recv_field[0][i][j] += frac_send_field_value;
          temp_frac_mask[0][i][0][j] += frac_mask_value * scale;

          // update ref_recv_field (avg)
          ref_recv_field[1][i][j] += frac_send_field_value * scale;
          temp_frac_mask[1][i][0][j] += frac_mask_value * scale;

          // update ref_recv_field (minimum)
          if (ref_recv_field[3][i][j] > frac_send_field_value) {
            ref_recv_field[3][i][j] = frac_send_field_value;
            temp_frac_mask[3][i][0][j] = frac_mask_value;
          }

          // update ref_recv_field (maximum)
          if (ref_recv_field[4][i][j] < frac_send_field_value) {
            ref_recv_field[4][i][j] = frac_send_field_value;
            temp_frac_mask[4][i][0][j] = frac_mask_value;
          }

        }

        // update ref_recv_field (none)
        ref_recv_field[2][i][j] = frac_send_field_value;
        temp_frac_mask[2][i][0][j] = frac_mask_value;
      }
    }

    if (is_target) {

      // target calls get every second timestep
      if ((t % TGT_DT) == 0) {
        for (int field_idx = 0; field_idx < NUM_FIELDS; ++field_idx) {

          // initialise recv_field
          double recv_field[COLLECTION_SIZE][NUM_POINTS];
          for (int j = 0; j < COLLECTION_SIZE; ++j)
            for (int k = 0; k < NUM_POINTS; ++k)
              recv_field[j][k] = -1;

          int info, ierror;
          yac_cget_(
            field_ids[field_idx], COLLECTION_SIZE, &recv_field[0][0],
            &info, &ierror);

          int ref_recv_info =
            (t%COUPLING_DT)?YAC_ACTION_NONE:YAC_ACTION_COUPLING;
          if (info != ref_recv_info) PUT_ERR("error in yac_cget_: wrong info");
          if (ierror) PUT_ERR("error in yac_cget_: wrong ierror");

          if (info == YAC_ACTION_COUPLING) {
            for (int j = 0; j < COLLECTION_SIZE; ++j) {
              for (int k = 0; k < NUM_POINTS; ++k) {
                if (temp_frac_mask[field_idx][j][0][k] != 0.0) {
                  if (fabs(recv_field[j][k] -
                           (ref_recv_field[field_idx][j][k] /
                            temp_frac_mask[field_idx][j][0][k])) > 1e-6)
                    PUT_ERR("error in yac_cget_: wrong recv_field (unmasked)");
                } else {
                  if (recv_field[j][k] != FRAC_MASK_VALUE)
                    PUT_ERR("error in yac_cget_: wrong recv_field (masked)");
                }
              }
            }
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
    for (int i = 0; i < COLLECTION_SIZE; ++i)
      for (int j = 0; j < NUM_POINTS; ++j)
        send_field[i][0][j] += i - 1;

    // clean ref_recv_field at every coupling timestep
    if ((t % COUPLING_DT) == 0) {
      init_ref_recv_field(ref_recv_field);
      init_temp_frac_mask(temp_frac_mask);
    }
  }

  yac_cfinalize ();

  return TEST_EXIT_CODE;
}
