/**
 * @file test_dummy_coupling9_c.c
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
  COLLECTION_SIZE = 3,
  NUM_POINTSETS = 1,
  NUM_POINTS = 9
};

static char const * mask_types[3] =
  {"even", "odd", "none"};

typedef int (*func_check_src_ptr)(
  double value, int idx, int collection_idx, int with_field_mask);
typedef int (*func_check_tgt_ptr)(
  double value, int idx, int collection_idx, int with_field_mask,
  func_check_src_ptr check_src);

static int check_src_no_mask(
  double value, int idx, int collection_idx, int with_field_mask);
static int check_src_even_mask(
  double value, int idx, int collection_idx, int with_field_mask);
static int check_src_odd_mask(
  double value, int idx, int collection_idx, int with_field_mask);
static int check_tgt_no_mask(
  double value, int idx, int collection_idx, int with_field_mask,
  func_check_src_ptr check_src);
static int check_tgt_even_mask(
  double value, int idx, int collection_idx, int with_field_mask,
  func_check_src_ptr check_src);
static int check_tgt_odd_mask(
  double value, int idx, int collection_idx, int with_field_mask,
  func_check_src_ptr check_src);

func_check_src_ptr check_src[3] = {
  check_src_even_mask, check_src_odd_mask, check_src_no_mask};
func_check_tgt_ptr check_tgt[3] = {
  check_tgt_even_mask, check_tgt_odd_mask, check_tgt_no_mask};

int main (void) {

  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cread_config_yaml("coupling_test9.yaml");

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
  char const * src_comp_name = "source_comp";
  char const * tgt_comp_name = "target_comp";
  yac_cdef_comp((is_target)?tgt_comp_name:src_comp_name, &comp_id);

  // define grid (both components use an identical grid)
  int grid_id;
  char const * src_grid_name = "source_grid";
  char const * tgt_grid_name = "target_grid";
  yac_cdef_grid_reg2d(
    (is_target)?tgt_grid_name:src_grid_name, (int[2]){3,3}, (int[2]){0,0},
    (double[]){0,1,2}, (double[]){0,1,2}, &grid_id);

  // define points at the vertices of the grid
  int point_id;
  yac_cdef_points_reg2d(
    grid_id, (int[2]){3,3}, YAC_LOCATION_CORNER,
    (double[]){0,1,2}, (double[]){0,1,2}, &point_id);

  // define masks for vertices
  int default_mask_id, dummy_mask_id;
  yac_cdef_mask(
    grid_id, 3*3, YAC_LOCATION_CORNER,
    (int[3*3]){1,0,0, 0,0,0, 0,0,0}, &default_mask_id);
  char const * even_mask_name = "even_mask";
  yac_cdef_mask_named(
    grid_id, 3*3, YAC_LOCATION_CORNER,
    (int[3*3]){1,0,1, 0,1,0, 1,0,1}, even_mask_name, &dummy_mask_id);
  char const * odd_mask_name = "odd_mask";
  yac_cdef_mask_named(
    grid_id, 3*3, YAC_LOCATION_CORNER,
    (int[3*3]){0,1,0, 1,0,1, 0,1,0}, odd_mask_name, &dummy_mask_id);

  // define field
  int field_ids[2][2][3][3];
  for (int config_from_file = 0; config_from_file < 2; ++config_from_file) {
    for (int with_field_mask = 0; with_field_mask < 2; ++with_field_mask) {
      for (int source_mask_type = 0; source_mask_type < 3; ++source_mask_type) {
        for (int target_mask_type = 0; target_mask_type < 3; ++target_mask_type) {
          char field_name[128];
          sprintf(
            field_name, "src2tgt_%s_%s_field_mask_%s_src_mask_%s_tgt_mask",
            config_from_file?"yaml":"manual",
            with_field_mask?"with":"without",
            mask_types[source_mask_type],
            mask_types[target_mask_type]);
          if (with_field_mask) {
            yac_cdef_field_mask(
              field_name, comp_id, &point_id, &default_mask_id, NUM_POINTSETS,
              COLLECTION_SIZE, "1", YAC_TIME_UNIT_SECOND,
              &field_ids[config_from_file]
                        [with_field_mask]
                        [source_mask_type]
                        [target_mask_type]);
          } else {
            yac_cdef_field(
              field_name, comp_id, &point_id, NUM_POINTSETS,
              COLLECTION_SIZE, "1", YAC_TIME_UNIT_SECOND,
              &field_ids[config_from_file]
                        [with_field_mask]
                        [source_mask_type]
                        [target_mask_type]);
          }
        }
      }
    }
  }

  // define manual couples
  int interp_stack_config;
  yac_cget_interp_stack_config(&interp_stack_config);
  yac_cadd_interp_stack_config_nnn(interp_stack_config, YAC_NNN_AVG, 1, 0.0);
  int ext_couple_config;
  yac_cget_ext_couple_config(&ext_couple_config);
  for (int with_field_mask = 0; with_field_mask < 2; ++with_field_mask) {
    for (int source_mask_type = 0; source_mask_type < 3; ++source_mask_type) {
      for (int target_mask_type = 0; target_mask_type < 3; ++target_mask_type) {
        if (source_mask_type != 2)
          yac_cset_ext_couple_config_src_mask_names(
            ext_couple_config, 1, source_mask_type?&odd_mask_name:&even_mask_name);
        else
          yac_cset_ext_couple_config_src_mask_names(ext_couple_config, 0, NULL);
        if (target_mask_type != 2)
          yac_cset_ext_couple_config_tgt_mask_name(
            ext_couple_config, target_mask_type?odd_mask_name:even_mask_name);
        else
          yac_cset_ext_couple_config_tgt_mask_name(ext_couple_config, NULL);

        size_t num_src_maks_names;
        char const * const * src_mask_names;
        yac_cget_ext_couple_config_src_mask_names(
          ext_couple_config, &num_src_maks_names, &src_mask_names);
        if ((num_src_maks_names != (source_mask_type != 2)) ||
            ((source_mask_type == 2) && (src_mask_names != NULL)) ||
            ((source_mask_type != 2) &&
             strcmp(src_mask_names[0],
                    source_mask_type?odd_mask_name:even_mask_name)))
          PUT_ERR("ERROR in yac_cget_ext_couple_config_src_mask_names");
        char const * tgt_mask_name;
        yac_cget_ext_couple_config_mask_name(ext_couple_config, &tgt_mask_name);
        if (((target_mask_type == 2) && (tgt_mask_name != NULL)) ||
            ((target_mask_type != 2) &&
             strcmp(tgt_mask_name, target_mask_type?odd_mask_name:even_mask_name)))
          PUT_ERR("ERROR in yac_cget_ext_couple_config_tgt_mask_name");

        char field_name[128];
          sprintf(
            field_name, "src2tgt_manual_%s_field_mask_%s_src_mask_%s_tgt_mask",
            with_field_mask?"with":"without",
            mask_types[source_mask_type],
            mask_types[target_mask_type]);
        yac_cdef_couple_custom(
          src_comp_name, src_grid_name, field_name,
          tgt_comp_name, tgt_grid_name, field_name,
          "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
        interp_stack_config, 0, 0, ext_couple_config);
      }
    }
  }
  yac_cfree_ext_couple_config(ext_couple_config);
  yac_cfree_interp_stack_config(interp_stack_config);

  yac_cenddef ( );

  double send_field[COLLECTION_SIZE][NUM_POINTSETS][NUM_POINTS] =
    {{{ 0, 1, 2, 3, 4, 5, 6, 7, 8}},
     {{10,11,12,13,14,15,16,17,18}},
     {{20,21,22,23,24,25,26,27,28}}};

  // do time steps
  for (int t = 0; t < 4; ++t) {

    if (!is_target) {

      for (int config_from_file = 0; config_from_file < 2; ++config_from_file) {
        for (int with_field_mask = 0; with_field_mask < 2; ++with_field_mask) {
          for (int source_mask_type = 0; source_mask_type < 3; ++source_mask_type) {
            for (int target_mask_type = 0; target_mask_type < 3; ++target_mask_type) {

              int info, ierror;
              yac_cput_(
                field_ids[config_from_file]
                         [with_field_mask]
                         [source_mask_type]
                         [target_mask_type],
                COLLECTION_SIZE, &send_field[0][0][0], &info, &ierror);
              if (info != YAC_ACTION_COUPLING)
                PUT_ERR("error in yac_cput_: wrong info");
              if (ierror) PUT_ERR("error in yac_cput_: wrong ierror");
            }
          }
        }
      }
    } else {

      for (int config_from_file = 0; config_from_file < 2; ++config_from_file) {
        for (int with_field_mask = 0; with_field_mask < 2; ++with_field_mask) {
          for (int source_mask_type = 0; source_mask_type < 3; ++source_mask_type) {
            for (int target_mask_type = 0; target_mask_type < 3; ++target_mask_type) {

              // initialise recv_field
              double recv_field[COLLECTION_SIZE][NUM_POINTS];
              for (int j = 0; j < COLLECTION_SIZE; ++j)
                for (int k = 0; k < NUM_POINTS; ++k)
                  recv_field[j][k] = -1;

              int info, ierror;
              yac_cget_(
                field_ids[config_from_file]
                         [with_field_mask]
                         [source_mask_type]
                         [target_mask_type],
                COLLECTION_SIZE, &recv_field[0][0], &info, &ierror);

              if (info != YAC_ACTION_COUPLING)
                PUT_ERR("error in yac_cget_: wrong info");
              if (ierror) PUT_ERR("error in yac_cget_: wrong ierror");

              for (int j = 0; j < COLLECTION_SIZE; ++j)
                for (int k = 0; k < NUM_POINTS; ++k)
                  if (check_tgt[target_mask_type](
                        recv_field[j][k], k, j, with_field_mask,
                        check_src[source_mask_type]))
                    PUT_ERR("error in yac_cget_: wrong recv_field");
            }
          }
        }
      }
    }
  }

  yac_cfinalize ();

  return TEST_EXIT_CODE;
}

static int check_src_even_mask(
  double value, int idx, int collection_idx, int with_field_mask) {

  UNUSED(idx);
  UNUSED(collection_idx);
  UNUSED(with_field_mask);
  return (int)value & 1;
}

static int check_src_odd_mask(
  double value, int idx, int collection_idx, int with_field_mask) {

  UNUSED(idx);
  UNUSED(collection_idx);
  UNUSED(with_field_mask);
  return !((int)value & 1);
}

static int check_src_no_mask(
  double value, int idx, int collection_idx, int with_field_mask) {

  return
    with_field_mask?
      ((int)value != (10 * collection_idx)):
      ((int)value != idx + 10 * collection_idx);
}

static int check_tgt_even_mask(
  double value, int idx, int collection_idx, int with_field_mask,
  func_check_src_ptr check_src) {

  return
    (idx&1)?
      (value != -1.0):check_src(value, idx, collection_idx, with_field_mask);
}

static int check_tgt_odd_mask(
  double value, int idx, int collection_idx, int with_field_mask,
  func_check_src_ptr check_src) {

  return
    (idx&1)?
      check_src(value, idx, collection_idx, with_field_mask):(value != -1.0);
}

static int check_tgt_no_mask(
  double value, int idx, int collection_idx, int with_field_mask,
  func_check_src_ptr check_src) {

  return
    (with_field_mask && idx)?
      (value != -1.0):check_src(value, idx, collection_idx, with_field_mask);
}
