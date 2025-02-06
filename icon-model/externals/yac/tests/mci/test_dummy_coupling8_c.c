// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "tests.h"
#include "test_common.h"
#include "yac.h"

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

static void test_aggregation();
static void test_changing_frac_mask();

int main (void) {


  test_aggregation();
  test_changing_frac_mask();

  yac_cfinalize ();

  return TEST_EXIT_CODE;
}

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

static void test_aggregation() {

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

  // define grid (both components use an identical grid)
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
    interp_stack_nnn, YAC_NNN_AVG, 1, 0.0, 0.0);

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

  yac_ccleanup();
}

static void test_changing_frac_mask() {

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

  // define grid (source has a 2x2 grid, while the target has a 1x1 grid)
  int grid_id;
  if (is_target)
    yac_cdef_grid_reg2d(
      "target_grid_2", (int[2]){2,2}, (int[2]){0,0},
      (double[]){0,2}, (double[]){0,2}, &grid_id);
  else
    yac_cdef_grid_reg2d(
      "source_grid_2", (int[2]){3,3}, (int[2]){0,0},
      (double[]){0,1,2}, (double[]){0,1,2}, &grid_id);

  // define points at the vertices of the grid
  int point_id;
  if (is_target)
    yac_cdef_points_reg2d(
      grid_id, (int[2]){1,1}, YAC_LOCATION_CELL,
      (double[]){1}, (double[]){1}, &point_id);
  else
    yac_cdef_points_reg2d(
      grid_id, (int[2]){2,2}, YAC_LOCATION_CELL,
      (double[]){0.5,1.5}, (double[]){0.5,1.5}, &point_id);

  // define field
  int field_id;
  yac_cdef_field(
    "field", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_id);
    yac_cenable_field_frac_mask(
      (is_target)?"target_comp":"source_comp",
      (is_target)?"target_grid_2":"source_grid_2",
      "field", FRAC_MASK_VALUE);

  // define interpolation stacks
  int interp_stack_nnn;
  yac_cget_interp_stack_config(&interp_stack_nnn);
  yac_cadd_interp_stack_config_nnn(
    interp_stack_nnn, YAC_NNN_AVG, 4, 0.0, 0.0);

  // define couplings
  yac_cdef_couple(
    "source_comp", "source_grid_2", "field",
    "target_comp", "target_grid_2", "field",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
    interp_stack_nnn, 0, 0);
  yac_cfree_interp_stack_config(interp_stack_nnn);

  yac_cenddef();

  enum {NUM_TIMESTEP = 5, SRC_FIELD_SIZE = 4};

  double send_field[NUM_TIMESTEP][SRC_FIELD_SIZE] =
    {{ 1, 2, 3, 4}, {10,11,12,13}, {20,21,22,23},
     {30,31,32,33}, {40,41,42,43}};
  double frac_mask[NUM_TIMESTEP][SRC_FIELD_SIZE] =
    {{1.0,1.0,0.0,0.0}, {1.0,0.0,1.0,0.0}, {1.0,1.0,0.0,0.0},
     {1.0,1.0,1.0,1.0}, {0.0,0.0,0.0,0.0}};

  for (int t = 0; t < NUM_TIMESTEP; ++t) {

    if (!is_target) {

      int info, ierror;
      yac_cput_frac_(
        field_id, 1, &send_field[t][0], &frac_mask[t][0], &info, &ierror);
      if (info != YAC_ACTION_COUPLING) PUT_ERR("error in yac_cput_: wrong info");
      if (ierror) PUT_ERR("error in yac_cput_: wrong ierror");
    }

    if (is_target) {

        double const weight = 0.25;
        double ref_recv_field = 0.0;
        double frac_weight_sum = 0.0;
        for (int i = 0; i < SRC_FIELD_SIZE; ++i) {
          ref_recv_field += send_field[t][i] * frac_mask[t][i] * weight;
          frac_weight_sum += weight * frac_mask[t][i];
        }
        ref_recv_field =
          (frac_weight_sum < 1.e-9)?
            FRAC_MASK_VALUE:(ref_recv_field/frac_weight_sum);

          double recv_field;
          int info, ierror;
          yac_cget_(field_id, 1, &recv_field, &info, &ierror);

          if (info != YAC_ACTION_COUPLING)
            PUT_ERR("error in yac_cget_: wrong info");
          if (ierror) PUT_ERR("error in yac_cget_: wrong ierror");
          if (fabs(recv_field - ref_recv_field) > 1.e-9)
            PUT_ERR("error in yac_cget_: wrong recv_field");
    }
  }

  yac_ccleanup();
}
