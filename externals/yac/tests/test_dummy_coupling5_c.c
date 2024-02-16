/**
 * @file test_dummy_coupling5_c.c
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

// if this is included in the wrong order SNAN is not defined
#include "feenableexcept.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tests.h"
#include "test_common.h"
#include "yac_interface.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

int main (void) {

  MPI_Init(NULL, NULL);

  feenableexcept(FE_INVALID);

  for (int with_core_mask = 0; with_core_mask < 2; ++with_core_mask) {
    for (int with_field_mask = 0; with_field_mask < 2; ++with_field_mask) {

      int instance_id;
      yac_cinit_instance(&instance_id);
      yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
      yac_cread_config_yaml_instance(instance_id, "coupling_test5.yaml");

      int size, rank;
      MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
      MPI_Comm_size ( MPI_COMM_WORLD, &size );

      if (size != 3) {
        fputs("wrong number of processes (has to be 3)\n", stderr);
        exit(EXIT_FAILURE);
      }

      int is_target = rank < 2;
      int is_source = rank > 0;
      YAC_ASSERT((is_target || is_source), "internal error")

      // define local component
      int comp_id;
      char const * comp_name = "main_comp";
      yac_cdef_comps_instance(instance_id, &comp_name, 1, &comp_id);

      // define grid (both components use an identical grids)
      int grid_ids[2];
      double x_vertices[2][3] = {{0.0*YAC_RAD,1.0*YAC_RAD,2.0*YAC_RAD},
                                {1.0*YAC_RAD,2.0*YAC_RAD,3.0*YAC_RAD}};
      double y_vertices[3] = {0.0*YAC_RAD,1.0*YAC_RAD,2.0*YAC_RAD};
      if (is_target)
        yac_cdef_grid_reg2d(
          "target_grid", (int[2]){3,3}, (int[2]){0,0},
          x_vertices[rank], y_vertices, &grid_ids[0]);
      if (is_source)
        yac_cdef_grid_reg2d(
          "source_grid", (int[2]){3,3}, (int[2]){0,0},
          x_vertices[rank-1], y_vertices, &grid_ids[1]);

      yac_int cell_global_index[2][4] = {{0,1,3,4}, {1,2,4,5}};
      yac_int corner_global_index[2][9] =
        {{0,1,2, 4,5,6, 8,9,10}, {1,2,3, 5,6,7, 9,10,11}};
      yac_int edge_global_index[2][12] =
        {{0,1,2,3,5, 7,8,9,10,12, 14,15}, {2,3,4,5,6, 9,10,11,12,13, 15,16}};
      int cell_core_mask[2][4] = {{1,1,1,0}, {0,1,1,1}};
      int corner_core_mask[2][9] = {{1,1,0, 0,1,0, 1,1,0}, {0,1,1, 0,1,1, 0,1,1}};
      int edge_core_mask[2][12] = {{1,1,1,1,0, 1,1,1,1,0, 1,0},
                                    {0,0,1,1,1, 0,0,1,1,1, 1,1}};
      int corner_field_mask[2][9] = {{1,1,1, 0,1,1, 1,1,1}, {1,1,1, 1,1,1, 1,1,1}};
      if (is_target) {
        yac_cset_global_index(
          cell_global_index[rank], YAC_LOCATION_CELL, grid_ids[0]);
        yac_cset_global_index(
          corner_global_index[rank], YAC_LOCATION_CORNER, grid_ids[0]);
        yac_cset_global_index(
          edge_global_index[rank], YAC_LOCATION_EDGE, grid_ids[0]);
        if (with_core_mask) {
          yac_cset_core_mask(
            cell_core_mask[rank], YAC_LOCATION_CELL, grid_ids[0]);
          yac_cset_core_mask(
            corner_core_mask[rank], YAC_LOCATION_CORNER, grid_ids[0]);
          yac_cset_core_mask(
            edge_core_mask[rank], YAC_LOCATION_EDGE, grid_ids[0]);
        }
      }
      if (is_source) {
        yac_cset_global_index(
          cell_global_index[rank-1], YAC_LOCATION_CELL, grid_ids[1]);
        yac_cset_global_index(
          corner_global_index[rank-1], YAC_LOCATION_CORNER, grid_ids[1]);
        yac_cset_global_index(
          edge_global_index[rank-1], YAC_LOCATION_EDGE, grid_ids[1]);
        if (with_core_mask) {
          yac_cset_core_mask(
            cell_core_mask[rank-1], YAC_LOCATION_CELL, grid_ids[1]);
          yac_cset_core_mask(
            corner_core_mask[rank-1], YAC_LOCATION_CORNER, grid_ids[1]);
          yac_cset_core_mask(
            edge_core_mask[rank-1], YAC_LOCATION_EDGE, grid_ids[1]);
        }
      }

      // define points at the vertices of the grid
      int point_ids[2];
      if (is_target) {
        yac_cdef_points_reg2d(
          grid_ids[0], (int[2]){3,3}, YAC_LOCATION_CORNER,
          x_vertices[rank], y_vertices, &point_ids[0]);
        if (with_field_mask)
          yac_cset_mask(corner_field_mask[rank], point_ids[0]);
      }
      if (is_source) {
        yac_cdef_points_reg2d(
          grid_ids[1], (int[2]){3,3}, YAC_LOCATION_CORNER,
          x_vertices[rank-1], y_vertices, &point_ids[1]);
        if (with_field_mask)
          yac_cset_mask(corner_field_mask[rank-1], point_ids[1]);
      }

      // define fields
      int field_ids[2][2];
      int dummy_field_ids[2];
      if (is_target) {
        yac_cdef_field(
          "field_a", comp_id, &point_ids[0], 1, 1, "3", YAC_TIME_UNIT_SECOND,
          &field_ids[0][0]);
        yac_cdef_field(
          "field_b", comp_id, &point_ids[0], 1, 1, "3", YAC_TIME_UNIT_SECOND,
          &field_ids[1][0]);
        yac_cdef_field(
          "dummy_field", comp_id, &point_ids[0], 1, 1, "3", YAC_TIME_UNIT_SECOND,
          &dummy_field_ids[0]);
      }
      if (is_source) {
        yac_cdef_field(
          "field_a", comp_id, &point_ids[1], 1, 1, "2", YAC_TIME_UNIT_SECOND,
          &field_ids[0][1]);
        yac_cdef_field(
          "field_b", comp_id, &point_ids[1], 1, 1, "2", YAC_TIME_UNIT_SECOND,
          &field_ids[1][1]);
        yac_cdef_field(
          "dummy_field", comp_id, &point_ids[1], 1, 1, "2", YAC_TIME_UNIT_SECOND,
          &dummy_field_ids[1]);
      }

      yac_cenddef_instance ( instance_id );

      double source_data[10][2][9] = // [put_idx][rank][field_idx]
        {{{ 0, 1, 2, 4, 5, 6, 8, 9,10}, { 1, 2, 3, 5, 6, 7, 9,10,11}}, // t =  0
         {{ 1, 2, 3, 5, 6, 7, 9,10,11}, { 2, 3, 4, 6, 7, 8,10,11,12}}, // t =  2
         {{ 2, 3, 4, 6, 7, 8,10,11,12}, { 3, 4, 5, 7, 8, 9,11,12,13}}, // t =  4
         {{ 3, 4, 5, 7, 8, 9,11,12,13}, { 4, 5, 6, 8, 9,10,12,13,14}}, // t =  6
         {{ 4, 5, 6, 8, 9,10,12,13,14}, { 5, 6, 7, 9,10,11,13,14,15}}, // t =  8
         {{ 5, 6, 7, 9,10,11,13,14,15}, { 6, 7, 8,10,11,12,14,15,16}}, // t = 10
         {{ 6, 7, 8,10,11,12,14,15,16}, { 7, 8, 9,11,12,13,15,16,17}}, // t = 12
         {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t = 14
         {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t = 16
         {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}};// t = 18
      double ref_recv_field[2][7][2][9] = // [reduction_idx][get_idx][rank][field_idx]
        {{{{ 0, 1, 2, 4, 5, 6, 8, 9,10}, { 1, 2, 3, 5, 6, 7, 9,10,11}}, // t =  0
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t =  3
          {{ 6, 9,12,18,21,24,30,33,36}, { 9,12,15,21,24,27,33,36,39}}, // t =  6
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t =  9
          {{15,18,21,27,30,33,39,42,45}, {18,21,24,30,33,36,42,45,48}}, // t = 12
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t = 15
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}},// t = 18
         {{{ 0, 1, 2, 4, 5, 6, 8, 9,10}, { 1, 2, 3, 5, 6, 7, 9,10,11}}, // t =  0
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t =  3
          {{ 3, 4, 5, 7, 8, 9,11,12,13}, { 4, 5, 6, 8, 9,10,12,13,14}}, // t =  6
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t =  9
          {{ 6, 7, 8,10,11,12,14,15,16}, { 7, 8, 9,11,12,13,15,16,17}}, // t = 12
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}, // t = 15
          {{-1,-1,-1,-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1,-1,-1,-1}}}};// t = 18

#ifdef SNAN
      if (is_source)
        for (size_t i = 0; i < 10; ++i)
          for (size_t j = 0; j < 2; ++j)
            for (size_t k = 0; k < 9; ++k)
              if ((with_core_mask && !corner_core_mask[rank-1][k]) ||
                  (with_field_mask && !corner_field_mask[rank-1][k]))
                source_data[i][j][k] = SNAN;
#endif
      if (is_target)
        for (size_t i = 0; i < 2; ++i)
          for (size_t j = 0; j < 7; ++j)
            for (size_t k = 0; k < 2; ++k)
              for (size_t l = 0; l < 9; ++l)
                if ((with_core_mask && !corner_core_mask[rank][l]) ||
                    (with_field_mask && !corner_field_mask[rank][l]))
                  ref_recv_field[i][j][k][l] = -1;

      int ref_send_info[2][10] =
        {{YAC_ACTION_COUPLING,        // t =  0
          YAC_ACTION_REDUCTION,       // t =  2
          YAC_ACTION_REDUCTION,       // t =  4
          YAC_ACTION_COUPLING,        // t =  6
          YAC_ACTION_REDUCTION,       // t =  8
          YAC_ACTION_REDUCTION,       // t = 10
          YAC_ACTION_PUT_FOR_RESTART, // t = 12
          YAC_ACTION_OUT_OF_BOUND,    // t = 14
          YAC_ACTION_OUT_OF_BOUND,    // t = 16
          YAC_ACTION_OUT_OF_BOUND},   // t = 18
        {YAC_ACTION_COUPLING,        // t =  0
          YAC_ACTION_NONE,            // t =  2
          YAC_ACTION_NONE,            // t =  4
          YAC_ACTION_COUPLING,        // t =  6
          YAC_ACTION_NONE,            // t =  8
          YAC_ACTION_NONE,            // t = 10
          YAC_ACTION_PUT_FOR_RESTART, // t = 12
          YAC_ACTION_OUT_OF_BOUND,    // t = 14
          YAC_ACTION_OUT_OF_BOUND,    // t = 16
          YAC_ACTION_OUT_OF_BOUND}};  // t = 18
      int ref_recv_info[7] =
        {YAC_ACTION_COUPLING,        // t =  0
         YAC_ACTION_NONE,            // t =  3
         YAC_ACTION_COUPLING,        // t =  6
         YAC_ACTION_NONE,            // t =  9
         YAC_ACTION_GET_FOR_RESTART, // t = 12
         YAC_ACTION_OUT_OF_BOUND,    // t = 15
         YAC_ACTION_OUT_OF_BOUND};   // t = 18

      int const source_period = 2;
      int const target_period = 3;
      int const test_time = 18;

      // do time steps
      for (int t = 0, put_idx = 0, get_idx = 0; t <= test_time; ++t) {

        int is_even_timestep = !(t%2);
        int is_source_timestep = is_source && (!(t%source_period));
        int is_target_timestep = is_target && (!(t%target_period));

        int const collection_size = 1;
        double * send_field_vertex =
          (is_source_timestep)?(source_data[put_idx][rank-1]):NULL;
        double * send_fields[1] = {send_field_vertex};
        double ** send_field_collection_data[1] = {send_fields};
        double recv_field_data[9] = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
        double * recv_field_vertex =
          (is_target_timestep)?recv_field_data:NULL;
        double * recv_field_collection_data[1] = {recv_field_vertex};

        int send_info, recv_info, send_action, recv_action;
        int dummy_send_ierr, dummy_recv_ierr;

        if (is_source) {
          yac_cget_action(dummy_field_ids[1], &send_action);
          if (send_action != YAC_ACTION_NONE) PUT_ERR("ERROR in dummy_send_action");
          if (is_even_timestep) {
            yac_cput(
              dummy_field_ids[1], collection_size, send_field_collection_data,
              &send_info, &dummy_send_ierr);
            if (dummy_send_ierr != 0) PUT_ERR("ERROR in dummy_send_ierr");
            if (send_info != YAC_ACTION_NONE) PUT_ERR("ERROR in dummy_send_info");
          } else {
            yac_cupdate(dummy_field_ids[1]);
          }
        }
        if (is_target) {
          yac_cget_action(dummy_field_ids[0], &recv_action);
          if (recv_action != YAC_ACTION_NONE) PUT_ERR("ERROR in dummy_recv_action");
          if (is_even_timestep) {
            yac_cget(
              dummy_field_ids[0], collection_size, recv_field_collection_data,
              &recv_info, &dummy_recv_ierr);
            if (dummy_recv_ierr != 0) PUT_ERR("ERROR in dummy_recv_ierr");
            if (recv_info != YAC_ACTION_NONE) PUT_ERR("ERROR in dummy_recv_info");
          } else {
            yac_cupdate(dummy_field_ids[0]);
          }
        }

        for (int i = 0; i < 2; ++i) {
          int ierr = 0;
          if (is_source_timestep && !is_target_timestep) {
            yac_cget_action(field_ids[i][1], &send_action);
            if (is_even_timestep && (send_action == YAC_ACTION_NONE)) {
              yac_cupdate(field_ids[i][1]);
              send_info = YAC_ACTION_NONE;
            } else {
              yac_cput(
                field_ids[i][1], collection_size, send_field_collection_data,
                &send_info, &ierr);
            }
          } else if (!is_source_timestep && is_target_timestep) {
            yac_cget_action(field_ids[i][0], &recv_action);
            if (is_even_timestep && (recv_action == YAC_ACTION_NONE)) {
              yac_cupdate(field_ids[i][0]);
              recv_info = YAC_ACTION_NONE;
            } else {
              yac_cget(
                field_ids[i][0], collection_size, recv_field_collection_data,
                &recv_info, &ierr);
            }
          } else if (is_source_timestep && is_target_timestep) {
            yac_cget_action(field_ids[i][1], &send_action);
            yac_cget_action(field_ids[i][0], &recv_action);
            yac_cexchange(
              field_ids[i][1], field_ids[i][0], collection_size,
              send_field_collection_data, recv_field_collection_data,
              &send_info, &recv_info, &ierr);
          }

          if (ierr != 0) PUT_ERR("ERROR in ierr");
          if (is_source_timestep) {
            if (send_info != ref_send_info[i][put_idx])
              PUT_ERR("ERROR in send_info");
            if (send_action != ref_send_info[i][put_idx])
              PUT_ERR("ERROR in send_action");
          }
          if (is_target_timestep) {
            if (recv_info != ref_recv_info[get_idx])
              PUT_ERR("ERROR in recv_info");
            if (recv_action != ref_recv_info[get_idx])
              PUT_ERR("ERROR in recv_action");
            for (int j = 0; j < 9; ++j)
              if (fabs(recv_field_data[j] -
                      ref_recv_field[i][get_idx][rank][j]) > 1e-9)
                PUT_ERR("ERROR in recv_field_data");
          }
        }

        if (is_source_timestep) ++put_idx;
        if (is_target_timestep) ++get_idx;
      }

      yac_cfinalize_instance(instance_id);
    } // with field mask
  } // with core mask

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
