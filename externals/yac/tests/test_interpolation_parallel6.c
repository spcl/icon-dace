/**
 * @file test_interpolation_parallel6.c
 *
 * @copyright Copyright  (C)  2022 DKRZ, MPI-M
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
#include <limits.h>
#include <math.h>

#include <mpi.h>
#include "tests.h"
#include "test_common.h"
#include "interpolation.h"

#define UNSET (1337.0)
#define FRAC_MASK_VALUE (-1337.0)

#define SCALE(X) \
  ((X) * scaling_factors[scaling_factor_idx] + \
   scaling_summand[scaling_summand_idx])

static int check_interpolation(
  struct interpolation * interp,
  double * src_data_raw, size_t num_src_fields, size_t src_field_size,
  double * ref_tgt_data_raw, size_t tgt_field_size, size_t collection_size);
static int check_interpolation_frac(
  struct interpolation * interp,
  double * src_data_raw, double * src_frac_mask_raw,
  size_t num_src_fields, size_t src_field_size,
  double * ref_tgt_data_raw, size_t tgt_field_size, size_t collection_size);

int main(int argc, char *argv[]) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  double scaling_factors[] = {1.0, 0.1, -1.0, 0.5, 10.0};
  double scaling_summand[] = {0.0, -1.0, 1.0, 10.0};
  enum {
    num_procs = 4,
    NUM_SCALING_FACTORS = sizeof(scaling_factors) / sizeof(scaling_factors[0]),
    NUM_SCALING_SUMMAND = sizeof(scaling_summand) / sizeof(scaling_summand[0]),
  };

  if (comm_size != num_procs) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  { // testing yac_interpolation_add_fixed
    // 16 target points, all odd positions receive -1 and the others 1

    enum {num_fixed_values = 2};
    double value[num_fixed_values] = {-1.0, 1.0};
    size_t count[num_fixed_values] = {8, 8};
    size_t pos[num_fixed_values][8] = {{0, 2, 4, 6, 8, 10, 12, 14},
                                       {1, 3, 5, 7, 9, 11, 13, 15}};

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {

        enum {collection_size = 1};
        struct interpolation * interp =
          yac_interpolation_new(
            collection_size, YAC_FRAC_MASK_NO_VALUE,
            scaling_factors[scaling_factor_idx],
            scaling_summand[scaling_summand_idx]);
        for (int i = 0; i < num_fixed_values; ++i)
          yac_interpolation_add_fixed(
            interp, value[i], count[i], pos[i]);

        enum {num_src_fields = 1};
        enum {src_field_size = 0};
        double * src_data_raw = NULL;
        enum {tgt_field_size = 16};
        double ref_tgt_data_raw[collection_size][tgt_field_size] =
          {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1}};

        if (check_interpolation(
              interp, src_data_raw, num_src_fields, src_field_size,
              &ref_tgt_data_raw[0][0], tgt_field_size, collection_size))
          PUT_ERR("ERROR in yac_interpolation_add_fixed");

        yac_interpolation_delete(interp);
      }
    }
  }

  { // testing yac_interpolation_add_direct
    // rank 0 pos {0,2} -> rank 1 pos {0,1}
    // rank 0 pos {0,2} -> rank 2 pos {0,2}
    // rank 1 pos {0,2} -> rank 2 pos {1,3}
    // (rank 0 is only source
    //  rank 1 is source and target
    //  rank 2 is only target
    //  rank 3 is neither source nor target)
    Xt_int src_indices[num_procs][4] =
      {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    size_t num_src_indices[num_procs] = {4,4,0,0};
    Xt_int dst_indices[num_procs][4] = {{-1},{0,2},{0,4,2,6},{-1}};
    size_t num_dst_indices[num_procs] = {0,2,4,0};
    Xt_idxlist src_idxlist =
      xt_idxvec_new(src_indices[comm_rank], num_src_indices[comm_rank]);
    Xt_idxlist dst_idxlist =
      xt_idxvec_new(dst_indices[comm_rank], num_dst_indices[comm_rank]);
    Xt_xmap xmap
      = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {

        enum {collection_size = 1};
        struct interpolation * interp =
          yac_interpolation_new(
            collection_size, YAC_FRAC_MASK_NO_VALUE,
            scaling_factors[scaling_factor_idx],
            scaling_summand[scaling_summand_idx]);
        yac_interpolation_add_direct(interp, redist);

        enum {num_src_fields = 1};
        size_t const src_field_size[num_procs] = {4,4,0,0};
        double src_data_raw[num_procs][num_src_fields][4] =
          {{{0,1,2,3}},{{4,5,6,7}},{{-1}},{{-1}}};
        size_t const tgt_field_size[num_procs] = {0,4,4,0};
        double ref_tgt_data_raw[num_procs][collection_size][4] =
          {{{UNSET}},
           {{SCALE(0),SCALE(2),UNSET,UNSET}},
           {{SCALE(0),SCALE(4),SCALE(2),SCALE(6)}},
           {{UNSET}}};

        if (check_interpolation(
              interp,
              (src_field_size[comm_rank] > 0)?
                (&src_data_raw[comm_rank][0][0]):NULL,
              num_src_fields, src_field_size[comm_rank],
              (tgt_field_size[comm_rank] > 0)?
                (&ref_tgt_data_raw[comm_rank][0][0]):NULL,
              tgt_field_size[comm_rank], collection_size))
          PUT_ERR("ERROR in yac_interpolation_add_direct");

        yac_interpolation_delete(interp);
      }
    }
    xt_redist_delete(redist);
  }

  { // testing yac_interpolation_add_direct with dynamic fractional masking
    // rank 0 pos {0,2} -> rank 1 pos {0,1} (pos 2 is masked)
    // rank 0 pos {0,2} -> rank 2 pos {0,2} (pos 2 is masked)
    // rank 1 pos {0,2} -> rank 2 pos {1,3}
    // (rank 0 is only source
    //  rank 1 is source and target
    //  rank 2 is only target
    //  rank 3 is neither source nor target)
    Xt_int src_indices[num_procs][4] =
      {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    size_t num_src_indices[num_procs] = {4,4,0,0};
    Xt_int dst_indices[num_procs][4] = {{-1},{0,2},{0,4,2,6},{-1}};
    size_t num_dst_indices[num_procs] = {0,2,4,0};
    Xt_idxlist src_idxlist =
      xt_idxvec_new(src_indices[comm_rank], num_src_indices[comm_rank]);
    Xt_idxlist dst_idxlist =
      xt_idxvec_new(dst_indices[comm_rank], num_dst_indices[comm_rank]);
    Xt_xmap xmap
      = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {

        enum {collection_size = 1};
        struct interpolation * interp =
          yac_interpolation_new(
            collection_size, FRAC_MASK_VALUE,
            scaling_factors[scaling_factor_idx],
            scaling_summand[scaling_summand_idx]);
        yac_interpolation_add_direct(interp, redist);

        enum {num_src_fields = 1};
        size_t const src_field_size[num_procs] = {4,4,0,0};
        double src_data_raw[num_procs][num_src_fields][4] =
          {{{0,1,2,3}},{{4,5,6,7}},{{-1}},{{-1}}};
        double src_frac_mask_raw[num_procs][num_src_fields][4] =
          {{{1.0,1.0,0.0,1.0}},{{0.5,0.5,0.5,0.5}},{{-1}},{{-1}}};
        size_t const tgt_field_size[num_procs] = {0,4,4,0};
        double ref_tgt_data_raw[num_procs][collection_size][4] =
          {{{UNSET}},
           {{SCALE(0),FRAC_MASK_VALUE,UNSET,UNSET}},
           {{SCALE(0),SCALE(4),FRAC_MASK_VALUE,SCALE(6)}},
           {{UNSET}}};

        if (check_interpolation_frac(
              interp,
              (src_field_size[comm_rank] > 0)?
                (&src_data_raw[comm_rank][0][0]):NULL,
              (src_field_size[comm_rank] > 0)?
                (&src_frac_mask_raw[comm_rank][0][0]):NULL,
              num_src_fields, src_field_size[comm_rank],
              (tgt_field_size[comm_rank] > 0)?
                (&ref_tgt_data_raw[comm_rank][0][0]):NULL,
              tgt_field_size[comm_rank], collection_size))
          PUT_ERR("ERROR in yac_interpolation_add_direct");

        yac_interpolation_delete(interp);
      }
    }
    xt_redist_delete(redist);
  }

  { // testing yac_interpolation_add_direct_mf
    // rank 0 field 0 id {0,2} field 1 id {1,3} -> rank 1 pos {0,2,3,1}
    // rank 0 field 0 id {1}   field 1 id {2}   -> rank 2 pos {0,2}
    // rank 1 field 0 id {4}   field 1 id {5}   -> rank 2 pos {1,3}
    // (rank 0 is only source
    //  rank 1 is source and target
    //  rank 2 is only target
    //  rank 3 is neither source nor target)
    enum {num_src_fields = 2};
    enum {num_src_indices = 4};
    Xt_int src_indices[num_procs][num_src_indices] =
      {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    Xt_int dst_indices[num_procs][2][2] = {{{-1},{-1}},
                                           {{0,2},{1,3}},
                                           {{1,4},{2,5}},
                                           {{-1},{-1}}};
    int src_offsets[num_procs] = {0,1,2,3};
    int dst_offsets[num_procs][num_src_fields][2] = {{{-1},{-1}},
                                                     {{0,2},{3,1}},
                                                     {{0,1},{2,3}},
                                                     {{-1},{-1}}};
    size_t num_dst_indices[num_procs][num_src_fields] =
      {{0,0},{2,2},{2,2},{0,0}};
    Xt_redist redists[num_src_fields];
    for (int i = 0; i < num_src_fields; ++i) {
      Xt_idxlist src_idxlist =
        xt_idxvec_new(src_indices[comm_rank], num_src_indices);
      Xt_idxlist dst_idxlist =
        xt_idxvec_new(dst_indices[comm_rank][i], num_dst_indices[comm_rank][i]);
      Xt_xmap xmap =
        xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
      redists[i] =
        xt_redist_p2p_off_new(
          xmap, src_offsets, dst_offsets[comm_rank][i], MPI_DOUBLE);
      xt_xmap_delete(xmap);
      xt_idxlist_delete(dst_idxlist);
      xt_idxlist_delete(src_idxlist);
    }

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {

        enum {collection_size = 1};
        struct interpolation * interp =
          yac_interpolation_new(
            collection_size, YAC_FRAC_MASK_NO_VALUE,
            scaling_factors[scaling_factor_idx],
            scaling_summand[scaling_summand_idx]);
        yac_interpolation_add_direct_mf(interp, redists, num_src_fields);

        size_t const src_field_size[num_procs] = {4,4,4,4};
        double src_data_raw[num_procs][num_src_fields][num_src_indices] =
          {{{0,1,2,3},{00,10,20,30}},
          {{4,5,6,7},{40,50,60,70}},
          {{8,9,10,11},{80,90,100,110}},
          {{12,13,14,15},{120,130,140,150}}};
        size_t const tgt_field_size[num_procs] = {0,4,4,0};
        double ref_tgt_data_raw[num_procs][collection_size][4] =
          {{{UNSET}},
           {{SCALE(0),SCALE(30),SCALE(2),SCALE(10)}},
           {{SCALE(1),SCALE(4),SCALE(20),SCALE(50)}},
           {{UNSET}}};

        if (check_interpolation(
              interp,
              (src_field_size[comm_rank] > 0)?
                (&src_data_raw[comm_rank][0][0]):NULL,
              num_src_fields, src_field_size[comm_rank],
              (tgt_field_size[comm_rank] > 0)?
                (&ref_tgt_data_raw[comm_rank][0][0]):NULL,
              tgt_field_size[comm_rank], collection_size))
          PUT_ERR("ERROR in yac_interpolation_add_direct_mf");

        yac_interpolation_delete(interp);
      }
    }
    xt_redist_delete(redists[1]);
    xt_redist_delete(redists[0]);
  }

  { // testing yac_interpolation_add_direct_mf with dynamic fractional mask
    // rank 0 field 0 id {0,2} field 1 id {1,3} -> rank 1 pos {0,2,3,1}
    // rank 0 field 0 id {1}   field 1 id {2}   -> rank 2 pos {0,2}
    // rank 1 field 0 id {4}   field 1 id {5}   -> rank 2 pos {1,3}
    // (rank 0 is only source
    //  rank 1 is source and target
    //  rank 2 is only target
    //  rank 3 is neither source nor target)
    enum {num_src_fields = 2};
    enum {num_src_indices = 4};
    Xt_int src_indices[num_procs][num_src_indices] =
      {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    Xt_int dst_indices[num_procs][2][2] = {{{-1},{-1}},
                                           {{0,2},{1,3}},
                                           {{1,4},{2,5}},
                                           {{-1},{-1}}};
    int src_offsets[num_procs] = {0,1,2,3};
    int dst_offsets[num_procs][num_src_fields][2] = {{{-1},{-1}},
                                                     {{0,2},{3,1}},
                                                     {{0,1},{2,3}},
                                                     {{-1},{-1}}};
    size_t num_dst_indices[num_procs][num_src_fields] =
      {{0,0},{2,2},{2,2},{0,0}};
    Xt_redist redists[num_src_fields];
    for (int i = 0; i < num_src_fields; ++i) {
      Xt_idxlist src_idxlist =
        xt_idxvec_new(src_indices[comm_rank], num_src_indices);
      Xt_idxlist dst_idxlist =
        xt_idxvec_new(dst_indices[comm_rank][i], num_dst_indices[comm_rank][i]);
      Xt_xmap xmap =
        xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
      redists[i] =
        xt_redist_p2p_off_new(
          xmap, src_offsets, dst_offsets[comm_rank][i], MPI_DOUBLE);
      xt_xmap_delete(xmap);
      xt_idxlist_delete(dst_idxlist);
      xt_idxlist_delete(src_idxlist);
    }

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {

        enum {collection_size = 1};
        struct interpolation * interp =
          yac_interpolation_new(
            collection_size, FRAC_MASK_VALUE,
            scaling_factors[scaling_factor_idx],
            scaling_summand[scaling_summand_idx]);
        yac_interpolation_add_direct_mf(interp, redists, num_src_fields);

        size_t const src_field_size[num_procs] = {4,4,4,4};
        double src_data_raw[num_procs][num_src_fields][num_src_indices] =
          {{{0,1,2,3},{00,10,20,30}},
          {{4,5,6,7},{40,50,60,70}},
          {{8,9,10,11},{80,90,100,110}},
          {{12,13,14,15},{120,130,140,150}}};
        double src_frac_mask_raw[num_procs][num_src_fields][num_src_indices] =
          {{{1.0,0.0,1.0,1.0},{1.0,1.0,1.0,1.0}},
          {{0.0,1.0,1.0,1.0},{1.0,0.0,1.0,1.0}},
          {{1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0}},
          {{1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0}}};
        size_t const tgt_field_size[num_procs] = {0,4,4,0};
        double ref_tgt_data_raw[num_procs][collection_size][4] =
          {{{UNSET}},
          {{SCALE(0),SCALE(30),SCALE(2),SCALE(10)}},
          {{FRAC_MASK_VALUE,FRAC_MASK_VALUE,SCALE(20),FRAC_MASK_VALUE}},
          {{UNSET}}};

        if (check_interpolation_frac(
              interp,
              (src_field_size[comm_rank] > 0)?
                (&src_data_raw[comm_rank][0][0]):NULL,
              (src_field_size[comm_rank] > 0)?
                (&src_frac_mask_raw[comm_rank][0][0]):NULL,
              num_src_fields, src_field_size[comm_rank],
              (tgt_field_size[comm_rank] > 0)?
                (&ref_tgt_data_raw[comm_rank][0][0]):NULL,
              tgt_field_size[comm_rank], collection_size))
          PUT_ERR("ERROR in yac_interpolation_add_direct_mf");

        yac_interpolation_delete(interp);
      }
    }
    xt_redist_delete(redists[1]);
    xt_redist_delete(redists[0]);
  }

  { // test yac_interpolation_add_sum_at_src and
    //      yac_interpolation_add_weight_sum_mvp_at_src
    // distribution of source and destination (only field_id 0) ids
    //        rank |       0        |       1        |       2        |       3
    // ------------|----------------|----------------|----------------|---------------
    //  field_id 0 |  0,  1,  2,  3 |  4,  5,  6,  7 |  8,  9, 10, 11 | 12, 13, 14, 15
    //  field_id 1 | 16, 17, 18, 19 | 20, 21, 22, 23 | 24, 25, 26, 27 | 28, 29, 30, 31
    //
    // tgt id | src ids      | generate on rank | optional weights
    // -------|--------------|------------------|--------------------
    //    4   | 0, 17, 2, 19 | 0                | 0.1, 0.2, 0.3, 0.4
    //    5   | 2, 19, 20, 5 | 1                | 0.5, 0.6, 0.7, 0.8
    //    6   | 6, 22, 7, 23 | 1                | 0.9, 1.0, 1.1, 1.2
    //    7   | 23           | 1                | 1.3
    //    8   | 2, 19, 22, 7 | 0                | 1.4, 1.5, 1.6, 1.7
    //    9   | 0, 17, 2, 19 | 0                | 0.1, 0.2, 0.3, 0.4
    //   10   | 5, 6         | 1                | 2.2, 2.3
    //   11   | 21, 22       | 1                | 2.4, 2.5

    // at first the data for all stencils is gathered on a single process
    // this is done using the halo_redists
    enum {num_src_fields = 2};
    Xt_redist halo_redists[num_src_fields];
    {
      enum {num_src_indices = 4};
      Xt_int src_indices[num_procs][num_src_fields][num_src_indices] =
        {{{0,1,2,3},{16,17,18,19}},
         {{4,5,6,7},{20,21,22,23}},
         {{8,9,10,11},{24,25,26,27}},
         {{12,13,14,15},{28,29,30,31}}};
      Xt_int dst_indices[num_procs][num_src_fields][1] =
        {{{7},{22}},
         {{2},{19}},
         {{-1},{-1}},
         {{-1},{-1}}};
      size_t num_dst_indices[num_procs][num_src_fields] =
        {{1,1},{1,1},{0,0},{0,0}};
      for (int field_id = 0; field_id < num_src_fields; ++field_id) {
        Xt_idxlist src_idxlist =
          xt_idxvec_new(src_indices[comm_rank][field_id], num_src_indices);
        Xt_idxlist dst_idxlist =
          xt_idxvec_new(
            dst_indices[comm_rank][field_id],
            num_dst_indices[comm_rank][field_id]);
        Xt_xmap xmap =
          xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
        halo_redists[field_id] = xt_redist_p2p_new(xmap, MPI_DOUBLE);
        xt_xmap_delete(xmap);
        xt_idxlist_delete(dst_idxlist);
        xt_idxlist_delete(src_idxlist);
      }
    }

    // number of stencils to be computed locally (tgt id 4 and 9 have the same
    // stencil -> we compute it only once)
    size_t tgt_count[num_procs] = {2,5,0,0};
    // number of source points per stencil
    size_t num_src_per_tgt[num_procs][5] = {{4,4},{4,4,1,2,2},{-1},{-1}};
    double weights[num_procs][13] =
      {{0.1,0.2,0.3,0.4,
        1.4,1.5,1.6,1.7},
       {0.5,0.6,0.7,0.8,
        0.9,1.0,1.1,1.2,
        1.3,
        2.2,2.3,
        2.4,2.5},{-1},{-1}};
    // source field index for each source point used in the stencils (SIZE_MAX
    // indicates that the respective source point is in the halo buffer)
    size_t src_field_idx[num_procs][13] =
      {{0,1,0,1,
        0,1,SIZE_MAX,SIZE_MAX},
       {SIZE_MAX,SIZE_MAX,1,0,
        0,1,0,1,
        1,
        0,0,
        1,1},{-1},{-1}};
    // source index; position of the respective source point within the source
    // field or halo buffer
    size_t src_idx[num_procs][13] =
      {{0,1,2,3,
        2,3,1,0},
       {0,1,0,1,
        2,2,3,3,
        3,
        1,2,
        1,2},{-1},{-1}};

    // the results from the stencil computation get send directly to their
    // final location in the destination buffer on the respective process using
    // result_redist
    // (tgt id 4 and 9 have the same stencil)
    Xt_redist result_redist;
    {
      Xt_int src_indices[num_procs][5] = {{4,8},{5,6,7,10,11},{-1},{-1}};
      size_t num_src_indices[num_procs] = {2,5,0,0};
      Xt_int dst_indices[num_procs][4] = {{-1},{4,5,6,7},{8,4,10,11},{-1}};
      size_t num_dst_indices[num_procs] = {0,4,4,0};
      Xt_idxlist src_idxlist =
        xt_idxvec_new(src_indices[comm_rank], num_src_indices[comm_rank]);
      Xt_idxlist dst_idxlist =
        xt_idxvec_new(dst_indices[comm_rank], num_dst_indices[comm_rank]);
      Xt_xmap xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
      result_redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
      xt_xmap_delete(xmap);
      xt_idxlist_delete(dst_idxlist);
      xt_idxlist_delete(src_idxlist);
    }

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {
        for (int with_weights = 0; with_weights < 2; ++with_weights) {

          enum {collection_size = 1};
          struct interpolation * interp =
            yac_interpolation_new(
              collection_size, YAC_FRAC_MASK_NO_VALUE,
              scaling_factors[scaling_factor_idx],
              scaling_summand[scaling_summand_idx]);
          if (with_weights) {
            yac_interpolation_add_weight_sum_mvp_at_src(
              interp, halo_redists, tgt_count[comm_rank],
                num_src_per_tgt[comm_rank], weights[comm_rank],
                src_field_idx[comm_rank], src_idx[comm_rank], num_src_fields,
                result_redist);
          } else {
            yac_interpolation_add_sum_at_src(
              interp, halo_redists, tgt_count[comm_rank],
              num_src_per_tgt[comm_rank], src_field_idx[comm_rank],
              src_idx[comm_rank], num_src_fields, result_redist);
          }

          enum {src_field_size = 4};
          double src_data_raw[num_procs][2][src_field_size] =
            {{{0,1,2,3},    {16,17,18,19}},
            {{4,5,6,7},    {20,21,22,23}},
            {{8,9,10,11},  {24,25,26,27}},
            {{12,13,14,15},{28,29,30,31}}};
          size_t const tgt_field_size[num_procs] = {0,4,4,0};
          double ref_tgt_data_raw[2][num_procs][collection_size][4] =
            {{{{UNSET}},
              {{SCALE(0+17+2+19),SCALE(2+19+20+5),SCALE(6+22+7+23),SCALE(23)}},
              {{SCALE(2+19+22+7),SCALE(0+17+2+19),SCALE(5+6),SCALE(21+22)}},
              {{UNSET}}},
            {{{UNSET}},
              {{SCALE(0.1*0+0.2*17+0.3*2+0.4*19),
                SCALE(0.5*2+0.6*19+0.7*20+0.8*5),
                SCALE(0.9*6+1.0*22+1.1*7+1.2*23),
                SCALE(1.3*23)}},
              {{SCALE(1.4*2+1.5*19+1.6*22+1.7*7),
                SCALE(0.1*0+0.2*17+0.3*2+0.4*19),
                SCALE(2.2*5+2.3*6),
                SCALE(2.4*21+2.5*22)}},
              {{UNSET}}}};

          if (check_interpolation(
                interp, &src_data_raw[comm_rank][0][0],
                num_src_fields, src_field_size,
                &ref_tgt_data_raw[with_weights][comm_rank][0][0],
                tgt_field_size[comm_rank], collection_size))
            PUT_ERR("ERROR in yac_interpolation_add_sum_at_src/"
                    "yac_interpolation_add_weight_sum_mvp_at_src");

          yac_interpolation_delete(interp);
        }
      }
    }
    xt_redist_delete(result_redist);
    xt_redist_delete(halo_redists[1]);
    xt_redist_delete(halo_redists[0]);
  }

  { // test yac_interpolation_add_sum_at_src and
    //      yac_interpolation_add_weight_sum_mvp_at_src with dynamic fractional masking
    // distribution of source and destination (only field_id 0) ids
    //        rank |         0          |         1          |         2          |         3
    // ------------|--------------------|--------------------|--------------------|-------------------
    //  field_id 0 |   0,   1,   2,   3 |   4,   5,   6,   7 |   8,   9,  10,  11 |  12,  13,  14,  15
    //  frac mask  | 1.0, 1.0, 1.0, 1.0 | 1.0, 0.5, 0.0, 0.0 | 1.0, 1.0, 1.0, 1.0 | 1.0, 1.0, 1.0, 1.0
    //  field_id 1 |  16,  17,  18,  19 |  20,  21,  22,  23 |  24,  25,  26,  27 |  28,  29,  30,  31
    //  frac_mask  | 1.0, 1.0, 1.0, 1.0 | 0.5, 0.0, 0.0, 0.0 | 1.0, 1.0, 1.0, 1.0 | 1.0, 1.0, 1.0, 1.0
    //
    // tgt id | src ids      | generate on rank | optional weights
    // -------|--------------|------------------|--------------------
    //    4   | 0, 17, 2, 19 | 0                | 0.1, 0.2, 0.3, 0.4
    //    5   | 2, 19, 20, 5 | 1                | 0.5, 0.6, 0.7, 0.8
    //    6   | 6, 22, 7, 23 | 1                | 0.9, 1.0, 1.1, 1.2
    //    7   | 23           | 1                | 1.3
    //    8   | 2, 19, 22, 7 | 0                | 1.4, 1.5, 1.6, 1.7
    //    9   | 0, 17, 2, 19 | 0                | 0.1, 0.2, 0.3, 0.4
    //   10   | 5, 6         | 1                | 2.2, 2.3
    //   11   | 21, 22       | 1                | 2.4, 2.5

    // at first the data for all stencils is gathered on a single process
    // this is done using the halo_redists
    enum {num_src_fields = 2};
    Xt_redist halo_redists[num_src_fields];
    {
      enum {num_src_indices = 4};
      Xt_int src_indices[num_procs][num_src_fields][num_src_indices] =
        {{{0,1,2,3},{16,17,18,19}},
         {{4,5,6,7},{20,21,22,23}},
         {{8,9,10,11},{24,25,26,27}},
         {{12,13,14,15},{28,29,30,31}}};
      Xt_int dst_indices[num_procs][num_src_fields][1] =
        {{{7},{22}},
         {{2},{19}},
         {{-1},{-1}},
         {{-1},{-1}}};
      size_t num_dst_indices[num_procs][num_src_fields] =
        {{1,1},{1,1},{0,0},{0,0}};
      for (int field_id = 0; field_id < num_src_fields; ++field_id) {
        Xt_idxlist src_idxlist =
          xt_idxvec_new(src_indices[comm_rank][field_id], num_src_indices);
        Xt_idxlist dst_idxlist =
          xt_idxvec_new(
            dst_indices[comm_rank][field_id],
            num_dst_indices[comm_rank][field_id]);
        Xt_xmap xmap =
          xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
        halo_redists[field_id] = xt_redist_p2p_new(xmap, MPI_DOUBLE);
        xt_xmap_delete(xmap);
        xt_idxlist_delete(dst_idxlist);
        xt_idxlist_delete(src_idxlist);
      }
    }

    // number of stencils to be computed locally (tgt id 4 and 9 have the same
    // stencil -> we compute it only once)
    size_t tgt_count[num_procs] = {2,5,0,0};
    // number of source points per stencil
    size_t num_src_per_tgt[num_procs][5] = {{4,4},{4,4,1,2,2},{-1},{-1}};
    double weights[num_procs][13] =
      {{0.1,0.2,0.3,0.4,
        1.4,1.5,1.6,1.7},
       {0.5,0.6,0.7,0.8,
        0.9,1.0,1.1,1.2,
        1.3,
        2.2,2.3,
        2.4,2.5},{-1},{-1}};
    // source field index for each source point used in the stencils (SIZE_MAX
    // indicates that the respective source point is in the halo buffer)
    size_t src_field_idx[num_procs][13] =
      {{0,1,0,1,
        0,1,SIZE_MAX,SIZE_MAX},
       {SIZE_MAX,SIZE_MAX,1,0,
        0,1,0,1,
        1,
        0,0,
        1,1},{-1},{-1}};
    // source index; position of the respective source point within the source
    // field or halo buffer
    size_t src_idx[num_procs][13] =
      {{0,1,2,3,
        2,3,1,0},
       {0,1,0,1,
        2,2,3,3,
        3,
        1,2,
        1,2},{-1},{-1}};

    // the results from the stencil computation get send directly to their
    // final location in the destination buffer on the respective process using
    // result_redist
    // (tgt id 4 and 9 have the same stencil)
    Xt_redist result_redist;
    {
      Xt_int src_indices[num_procs][5] = {{4,8},{5,6,7,10,11},{-1},{-1}};
      size_t num_src_indices[num_procs] = {2,5,0,0};
      Xt_int dst_indices[num_procs][4] = {{-1},{4,5,6,7},{8,4,10,11},{-1}};
      size_t num_dst_indices[num_procs] = {0,4,4,0};
      Xt_idxlist src_idxlist =
        xt_idxvec_new(src_indices[comm_rank], num_src_indices[comm_rank]);
      Xt_idxlist dst_idxlist =
        xt_idxvec_new(dst_indices[comm_rank], num_dst_indices[comm_rank]);
      Xt_xmap xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
      result_redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
      xt_xmap_delete(xmap);
      xt_idxlist_delete(dst_idxlist);
      xt_idxlist_delete(src_idxlist);
    }

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {
        for (int with_weights = 0; with_weights < 2; ++with_weights) {

          enum {collection_size = 1};
          struct interpolation * interp =
            yac_interpolation_new(
              collection_size, FRAC_MASK_VALUE,
              scaling_factors[scaling_factor_idx],
              scaling_summand[scaling_summand_idx]);
          if (with_weights) {
            yac_interpolation_add_weight_sum_mvp_at_src(
              interp, halo_redists, tgt_count[comm_rank],
                num_src_per_tgt[comm_rank], weights[comm_rank],
                src_field_idx[comm_rank], src_idx[comm_rank], num_src_fields,
                result_redist);
          } else {
            yac_interpolation_add_sum_at_src(
              interp, halo_redists, tgt_count[comm_rank],
              num_src_per_tgt[comm_rank], src_field_idx[comm_rank],
              src_idx[comm_rank], num_src_fields, result_redist);
          }

          enum {src_field_size = 4};
          double src_data_raw[num_procs][2][src_field_size] =
            {{{0,1,2,3},    {16,17,18,19}},
            {{4,5,6,7},    {20,21,22,23}},
            {{8,9,10,11},  {24,25,26,27}},
            {{12,13,14,15},{28,29,30,31}}};
          double src_frac_mask_raw[num_procs][2][src_field_size] =
            {{{1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}},
            {{1.0, 0.5, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0}},
            {{1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}},
            {{1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}}};
          size_t const tgt_field_size[num_procs] = {0,4,4,0};
          double ref_tgt_data_raw[2][num_procs][collection_size][4] =
            {{{{UNSET}},
              {{SCALE((1.0*0+1.0*17+1.0*2+1.0*19)/(1.0+1.0+1.0+1.0)),
                SCALE((1.0*2+1.0*19+0.5*20+0.5*5)/(1.0+1.0+0.5+0.5)),
                FRAC_MASK_VALUE,FRAC_MASK_VALUE}},
              {{SCALE((1.0*2+1.0*19+0.0*22+0.0*7)/(1.0+1.0+0.0+0.0)),
                SCALE((1.0*0+1.0*17+1.0*2+1.0*19)/(1.0+1.0+1.0+1.0)),
                SCALE((0.5*5+0.0*6)/(0.5+0.0)),
                FRAC_MASK_VALUE}},
              {{UNSET}}},
            {{{UNSET}},
              {{SCALE((1.0*0.1*0+1.0*0.2*17+1.0*0.3*2+1.0*0.4*19)/
                      (1.0*0.1+1.0*0.2+1.0*0.3+1.0*0.4)),
                SCALE((1.0*0.5*2+1.0*0.6*19+0.5*0.7*20+0.5*0.8*5)/
                      (1.0*0.5+1.0*0.6+0.5*0.7+0.5*0.8)),
                FRAC_MASK_VALUE, FRAC_MASK_VALUE}},
              {{SCALE((1.0*1.4*2+1.0*1.5*19+0.0*1.6*22+0.0*1.7*7)/
                      (1.0*1.4+1.0*1.5+0.0*1.6+0.0*1.7)),
                SCALE((1.0*0.1*0+1.0*0.2*17+1.0*0.3*2+1.0*0.4*19)/
                      (1.0*0.1+1.0*0.2+1.0*0.3+1.0*0.4)),
                SCALE((0.5*2.2*5+0.0*2.3*6)/(0.5*2.2+0.0*2.3)),
                FRAC_MASK_VALUE}},
              {{UNSET}}}};

          if (check_interpolation_frac(
                interp, &src_data_raw[comm_rank][0][0],
                &src_frac_mask_raw[comm_rank][0][0],
                num_src_fields, src_field_size,
                &ref_tgt_data_raw[with_weights][comm_rank][0][0],
                tgt_field_size[comm_rank], collection_size))
            PUT_ERR("ERROR in yac_interpolation_add_sum_at_src/"
                    "yac_interpolation_add_weight_sum_mvp_at_src frac");

          yac_interpolation_delete(interp);
        }
      }
    }
    xt_redist_delete(result_redist);
    xt_redist_delete(halo_redists[1]);
    xt_redist_delete(halo_redists[0]);
  }

  { // test yac_interpolation_add_sum_at_tgt and
    //      yac_interpolation_add_weight_sum_mvp_at_tgt
    // distribution of source and destination (only field_id 0) ids
    //        rank |       0        |       1        |       2        |       3
    // ------------|----------------|----------------|----------------|---------------
    //  field_id 0 |  0,  1,  2,  3 |  4,  5,  6,  7 |  8,  9, 10, 11 | 12, 13, 14, 15
    //  field_id 1 | 16, 17, 18, 19 | 20, 21, 22, 23 | 24, 25, 26, 27 | 28, 29, 30, 31
    //
    // tgt id | src ids      | optional weights
    // -------|--------------|--------------------
    //    4   | 0, 17, 2, 19 | 0.1, 0.2, 0.3, 0.4
    //    5   | 2, 19, 20, 5 | 0.5, 0.6, 0.7, 0.8
    //    6   | 6, 22, 7, 23 | 0.9, 1.0, 1.1, 1.2
    //    7   | 23           | 1.3
    //    8   | 2, 19, 22, 7 | 1.4, 1.5, 1.6, 1.7
    //    9   | 0, 17, 2, 19 | 0.1, 0.2, 0.3, 0.4
    //   10   | 5, 6         | 2.2, 2.3
    //   11   | 21, 22       | 2.4, 2.5

    // at first the source points required for each target point have to be
    // transfered to the respective target processes using src_redists
    enum {num_src_fields = 2};
    Xt_redist src_redists[num_src_fields];
    {
      enum {num_src_indices = 4};
      Xt_int src_indices[num_procs][num_src_fields][num_src_indices] =
        {{{0,1,2,3},{16,17,18,19}},
         {{4,5,6,7},{20,21,22,23}},
         {{8,9,10,11},{24,25,26,27}},
         {{12,13,14,15},{28,29,30,31}}};

      Xt_int dst_indices[num_procs][num_src_fields][5] =
        {{{-1},{-1}},
         {{0,2},{17,19}},
         {{0,2,5,6,7},{17,19,21,22}},
         {{-1},{-1}}};
      size_t num_dst_indices[num_procs][num_src_fields] =
        {{0,0},{2,2},{5,4},{0,0}};
      for (int field_id = 0; field_id < num_src_fields; ++field_id) {
        Xt_idxlist src_idxlist =
          xt_idxvec_new(src_indices[comm_rank][field_id], num_src_indices);
        Xt_idxlist dst_idxlist =
          xt_idxvec_new(
            dst_indices[comm_rank][field_id],
            num_dst_indices[comm_rank][field_id]);
        Xt_xmap xmap =
          xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
        src_redists[field_id] = xt_redist_p2p_new(xmap, MPI_DOUBLE);
        xt_xmap_delete(xmap);
        xt_idxlist_delete(dst_idxlist);
        xt_idxlist_delete(src_idxlist);
      }
    }

    // positions of target points to be written
    size_t tgt_pos[num_procs][4] = {{-1},{0,1,2,3},{0,1,2,3},{-1}};
    // number of entries in tgt_pos
    size_t tgt_count[num_procs] = {0,4,4,0};
    // number of source points per stencil
    size_t num_src_per_tgt[num_procs][4] = {{-1},{4,4,4,1},{4,4,2,2},{-1}};
    double weights[num_procs][13] =
      {{-1},
       {0.1,0.2,0.3,0.4,
        0.5,0.6,0.7,0.8,
        0.9,1.0,1.1,1.2,
        1.3},
       {1.4,1.5,1.6,1.7,
        0.1,0.2,0.3,0.4,
        2.2,2.3,
        2.4,2.5},
       {-1}};
    // source field index for each source point used in the stencils (SIZE_MAX
    // indicates that the respective source point is in the halo buffer)
    size_t src_field_idx[num_procs][13] =
      {{-1},
       {SIZE_MAX,SIZE_MAX,SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX,1,0,
        0,1,0,1,
        1},
       {SIZE_MAX,SIZE_MAX,SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX,SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX},
       {-1}};
    // source index; position of the respective source point within the source
    // field or halo buffer
    size_t src_idx[num_procs][13] =
      {{-1},
       {0,2,1,3, 1,3,0,1, 2,2,3,3, 3},
       {1,6,8,4, 0,5,1,6, 2,3, 7,8},
       {-1}};

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {
        for (int with_weights = 0; with_weights < 2; ++with_weights) {

          enum {collection_size = 1};
          struct interpolation * interp =
            yac_interpolation_new(
              collection_size, YAC_FRAC_MASK_NO_VALUE,
              scaling_factors[scaling_factor_idx],
              scaling_summand[scaling_summand_idx]);
          if (with_weights) {
            yac_interpolation_add_weight_sum_mvp_at_tgt(
              interp, src_redists, tgt_pos[comm_rank], tgt_count[comm_rank],
              num_src_per_tgt[comm_rank], weights[comm_rank],
              src_field_idx[comm_rank], src_idx[comm_rank], num_src_fields);
          } else {
            yac_interpolation_add_sum_at_tgt(
              interp, src_redists, tgt_pos[comm_rank], tgt_count[comm_rank],
              num_src_per_tgt[comm_rank], src_field_idx[comm_rank],
              src_idx[comm_rank], num_src_fields);
          }

          enum {src_field_size = 4};
          double src_data_raw[num_procs][num_src_fields][src_field_size] =
            {{{0,1,2,3},    {16,17,18,19}},
            {{4,5,6,7},    {20,21,22,23}},
            {{8,9,10,11},  {24,25,26,27}},
            {{12,13,14,15},{28,29,30,31}}};
          size_t const tgt_field_size[num_procs] = {0,4,4,0};
          double ref_tgt_data_raw[2][num_procs][collection_size][4] =
          {{{{UNSET}},
            {{SCALE(0+17+2+19),SCALE(2+19+20+5),SCALE(6+22+7+23),SCALE(23)}},
            {{SCALE(2+19+22+7),SCALE(0+17+2+19),SCALE(5+6),SCALE(21+22)}},
            {{UNSET}}},
          {{{UNSET}},
            {{SCALE(0.1*0+0.2*17+0.3*2+0.4*19),
              SCALE(0.5*2+0.6*19+0.7*20+0.8*5),
              SCALE(0.9*6+1.0*22+1.1*7+1.2*23),
              SCALE(1.3*23)}},
            {{SCALE(1.4*2+1.5*19+1.6*22+1.7*7),
              SCALE(0.1*0+0.2*17+0.3*2+0.4*19),
              SCALE(2.2*5+2.3*6),
              SCALE(2.4*21+2.5*22)}},
            {{UNSET}}}};

          if (check_interpolation(
                interp, &src_data_raw[comm_rank][0][0],
                num_src_fields, src_field_size,
                &ref_tgt_data_raw[with_weights][comm_rank][0][0],
                tgt_field_size[comm_rank], collection_size))
            PUT_ERR("ERROR in yac_interpolation_add_sum_at_tgt/"
                    "yac_interpolation_add_weight_sum_mvp_at_tgt");

          yac_interpolation_delete(interp);
        }
      }
    }
    xt_redist_delete(src_redists[1]);
    xt_redist_delete(src_redists[0]);
  }

  { // test yac_interpolation_add_sum_at_tgt and
    //      yac_interpolation_add_weight_sum_mvp_at_tgt with dynamic fractional masking
    // distribution of source and destination (only field_id 0) ids
    //        rank |         0          |         1          |         2          |         3
    // ------------|--------------------|--------------------|--------------------|-------------------
    //  field_id 0 |   0,   1,   2,   3 |   4,   5,   6,   7 |   8,   9,  10,  11 |  12,  13,  14,  15
    //  frac mask  | 1.0, 1.0, 1.0, 1.0 | 1.0, 0.5, 0.0, 0.0 | 1.0, 1.0, 1.0, 1.0 | 1.0, 1.0, 1.0, 1.0
    //  field_id 1 |  16,  17,  18,  19 |  20,  21,  22,  23 |  24,  25,  26,  27 |  28,  29,  30,  31
    //  frac_mask  | 1.0, 1.0, 1.0, 1.0 | 0.5, 0.0, 0.0, 0.0 | 1.0, 1.0, 1.0, 1.0 | 1.0, 1.0, 1.0, 1.0
    //
    // tgt id | src ids      | optional weights
    // -------|--------------|--------------------
    //    4   | 0, 17, 2, 19 | 0.1, 0.2, 0.3, 0.4
    //    5   | 2, 19, 20, 5 | 0.5, 0.6, 0.7, 0.8
    //    6   | 6, 22, 7, 23 | 0.9, 1.0, 1.1, 1.2
    //    7   | 23           | 1.3
    //    8   | 2, 19, 22, 7 | 1.4, 1.5, 1.6, 1.7
    //    9   | 0, 17, 2, 19 | 0.1, 0.2, 0.3, 0.4
    //   10   | 5, 6         | 2.2, 2.3
    //   11   | 21, 22       | 2.4, 2.5

    // at first the source points required for each target point have to be
    // transfered to the respective target processes using src_redists
    enum {num_src_fields = 2};
    Xt_redist src_redists[num_src_fields];
    {
      enum {num_src_indices = 4};
      Xt_int src_indices[num_procs][num_src_fields][num_src_indices] =
        {{{0,1,2,3},{16,17,18,19}},
         {{4,5,6,7},{20,21,22,23}},
         {{8,9,10,11},{24,25,26,27}},
         {{12,13,14,15},{28,29,30,31}}};

      Xt_int dst_indices[num_procs][num_src_fields][5] =
        {{{-1},{-1}},
         {{0,2},{17,19}},
         {{0,2,5,6,7},{17,19,21,22}},
         {{-1},{-1}}};
      size_t num_dst_indices[num_procs][num_src_fields] =
        {{0,0},{2,2},{5,4},{0,0}};
      for (int field_id = 0; field_id < num_src_fields; ++field_id) {
        Xt_idxlist src_idxlist =
          xt_idxvec_new(src_indices[comm_rank][field_id], num_src_indices);
        Xt_idxlist dst_idxlist =
          xt_idxvec_new(
            dst_indices[comm_rank][field_id],
            num_dst_indices[comm_rank][field_id]);
        Xt_xmap xmap =
          xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
        src_redists[field_id] = xt_redist_p2p_new(xmap, MPI_DOUBLE);
        xt_xmap_delete(xmap);
        xt_idxlist_delete(dst_idxlist);
        xt_idxlist_delete(src_idxlist);
      }
    }

    // positions of target points to be written
    size_t tgt_pos[num_procs][4] = {{-1},{0,1,2,3},{0,1,2,3},{-1}};
    // number of entries in tgt_pos
    size_t tgt_count[num_procs] = {0,4,4,0};
    // number of source points per stencil
    size_t num_src_per_tgt[num_procs][4] = {{-1},{4,4,4,1},{4,4,2,2},{-1}};
    double weights[num_procs][13] =
      {{-1},
       {0.1,0.2,0.3,0.4,
        0.5,0.6,0.7,0.8,
        0.9,1.0,1.1,1.2,
        1.3},
       {1.4,1.5,1.6,1.7,
        0.1,0.2,0.3,0.4,
        2.2,2.3,
        2.4,2.5},
       {-1}};
    // source field index for each source point used in the stencils (SIZE_MAX
    // indicates that the respective source point is in the halo buffer)
    size_t src_field_idx[num_procs][13] =
      {{-1},
       {SIZE_MAX,SIZE_MAX,SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX,1,0,
        0,1,0,1,
        1},
       {SIZE_MAX,SIZE_MAX,SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX,SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX,
        SIZE_MAX,SIZE_MAX},
       {-1}};
    // source index; position of the respective source point within the source
    // field or halo buffer
    size_t src_idx[num_procs][13] =
      {{-1},
       {0,2,1,3, 1,3,0,1, 2,2,3,3, 3},
       {1,6,8,4, 0,5,1,6, 2,3, 7,8},
       {-1}};

    for (int scaling_factor_idx = 0;
         scaling_factor_idx < NUM_SCALING_FACTORS; ++scaling_factor_idx) {
      for (int scaling_summand_idx = 0;
           scaling_summand_idx < NUM_SCALING_SUMMAND; ++scaling_summand_idx) {
        for (int with_weights = 0; with_weights < 2; ++with_weights) {

          enum {collection_size = 1};
          struct interpolation * interp =
            yac_interpolation_new(
              collection_size, FRAC_MASK_VALUE,
              scaling_factors[scaling_factor_idx],
              scaling_summand[scaling_summand_idx]);
          if (with_weights) {
            yac_interpolation_add_weight_sum_mvp_at_tgt(
              interp, src_redists, tgt_pos[comm_rank], tgt_count[comm_rank],
              num_src_per_tgt[comm_rank], weights[comm_rank],
              src_field_idx[comm_rank], src_idx[comm_rank], num_src_fields);
          } else {
            yac_interpolation_add_sum_at_tgt(
              interp, src_redists, tgt_pos[comm_rank], tgt_count[comm_rank],
              num_src_per_tgt[comm_rank], src_field_idx[comm_rank],
              src_idx[comm_rank], num_src_fields);
          }

          enum {src_field_size = 4};
          double src_data_raw[num_procs][num_src_fields][src_field_size] =
            {{{0,1,2,3},    {16,17,18,19}},
            {{4,5,6,7},    {20,21,22,23}},
            {{8,9,10,11},  {24,25,26,27}},
            {{12,13,14,15},{28,29,30,31}}};
          double src_frac_mask_raw[num_procs][2][src_field_size] =
            {{{1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}},
            {{1.0, 0.5, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0}},
            {{1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}},
            {{1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}}};
          size_t const tgt_field_size[num_procs] = {0,4,4,0};
          double ref_tgt_data_raw[2][num_procs][collection_size][4] =
            {{{{UNSET}},
              {{SCALE((1.0*0+1.0*17+1.0*2+1.0*19)/(1.0+1.0+1.0+1.0)),
                SCALE((1.0*2+1.0*19+0.5*20+0.5*5)/(1.0+1.0+0.5+0.5)),
                FRAC_MASK_VALUE,FRAC_MASK_VALUE}},
              {{SCALE((1.0*2+1.0*19+0.0*22+0.0*7)/(1.0+1.0+0.0+0.0)),
                SCALE((1.0*0+1.0*17+1.0*2+1.0*19)/(1.0+1.0+1.0+1.0)),
                SCALE((0.5*5+0.0*6)/(0.5+0.0)),
                FRAC_MASK_VALUE}},
              {{UNSET}}},
            {{{UNSET}},
              {{SCALE((1.0*0.1*0+1.0*0.2*17+1.0*0.3*2+1.0*0.4*19)/
                      (1.0*0.1+1.0*0.2+1.0*0.3+1.0*0.4)),
                SCALE((1.0*0.5*2+1.0*0.6*19+0.5*0.7*20+0.5*0.8*5)/
                      (1.0*0.5+1.0*0.6+0.5*0.7+0.5*0.8)),
                FRAC_MASK_VALUE, FRAC_MASK_VALUE}},
              {{SCALE((1.0*1.4*2+1.0*1.5*19+0.0*1.6*22+0.0*1.7*7)/
                      (1.0*1.4+1.0*1.5+0.0*1.6+0.0*1.7)),
                SCALE((1.0*0.1*0+1.0*0.2*17+1.0*0.3*2+1.0*0.4*19)/
                      (1.0*0.1+1.0*0.2+1.0*0.3+1.0*0.4)),
                SCALE((0.5*2.2*5+0.0*2.3*6)/(0.5*2.2+0.0*2.3)),
                FRAC_MASK_VALUE}},
              {{UNSET}}}};

          if (check_interpolation_frac(
                interp,
                &src_data_raw[comm_rank][0][0],
                &src_frac_mask_raw[comm_rank][0][0],
                num_src_fields, src_field_size,
                &ref_tgt_data_raw[with_weights][comm_rank][0][0],
                tgt_field_size[comm_rank], collection_size))
            PUT_ERR("ERROR in yac_interpolation_add_sum_at_tgt/"
                    "yac_interpolation_add_weight_sum_mvp_at_tgt frac");

          yac_interpolation_delete(interp);
        }
      }
    }
    xt_redist_delete(src_redists[1]);
    xt_redist_delete(src_redists[0]);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void init_tgt_data(
  double ** tgt_data, size_t tgt_field_size, size_t collection_size) {

  if (tgt_data == NULL) return;

  for (size_t i = 0; i < collection_size; ++i)
    for (size_t j = 0; j < tgt_field_size; ++j)
      tgt_data[i][j] = UNSET;
}

static int check_tgt_data(
  double ** tgt_data, double * ref_tgt_data_raw,
  size_t tgt_field_size, size_t collection_size) {

  int diff_count = 0;
  if (tgt_data != NULL)
    for (size_t i = 0; i < collection_size;
        ++i, ref_tgt_data_raw += tgt_field_size)
      for (size_t j = 0; j < tgt_field_size; ++j)
        if (fabs(tgt_data[i][j] - ref_tgt_data_raw[j]) > 1e-9)
          ++diff_count;
  return diff_count;
}

static int check_interpolation_execute(
  struct interpolation * interp, double *** src_data, double ** tgt_data,
  double * ref_tgt_data_raw, size_t tgt_field_size, size_t collection_size) {

  int diff_count = 0;

  // there should be not active interpolation exchange
  if (!yac_interpolation_execute_test(interp))
    PUT_ERR("ERROR in yac_interpolation_execute_test");

  // synchronous exchange
  init_tgt_data(tgt_data, tgt_field_size, collection_size);
  yac_interpolation_execute(interp, src_data, tgt_data);
  diff_count +=
    check_tgt_data(tgt_data, ref_tgt_data_raw, tgt_field_size, collection_size);

  MPI_Barrier(MPI_COMM_WORLD);

  // there should be not active interpolation exchange
  if (!yac_interpolation_execute_test(interp))
    PUT_ERR("ERROR in yac_interpolation_execute_test");

  // first put then get
  if (src_data != NULL) {
    yac_interpolation_execute_put(interp, src_data);

    // this should return at some point
    if (tgt_data == NULL)
      while (!yac_interpolation_execute_test(interp));
  }
  if (tgt_data != NULL) yac_interpolation_execute_get(interp, tgt_data);
  diff_count +=
    check_tgt_data(tgt_data, ref_tgt_data_raw, tgt_field_size, collection_size);

  return diff_count;
}

static int check_interpolation_execute_frac(
  struct interpolation * interp, double *** src_data, double *** src_frac_mask,
  double ** tgt_data, double * ref_tgt_data_raw, size_t tgt_field_size,
  size_t collection_size) {

  int diff_count = 0;

  // there should be not active interpolation exchange
  if (!yac_interpolation_execute_test(interp))
    PUT_ERR("ERROR in yac_interpolation_execute_test");

  // synchronous exchange
  init_tgt_data(tgt_data, tgt_field_size, collection_size);
  yac_interpolation_execute_frac(interp, src_data, src_frac_mask, tgt_data);
  diff_count +=
    check_tgt_data(tgt_data, ref_tgt_data_raw, tgt_field_size, collection_size);

  MPI_Barrier(MPI_COMM_WORLD);

  // there should be not active interpolation exchange
  if (!yac_interpolation_execute_test(interp))
    PUT_ERR("ERROR in yac_interpolation_execute_test");

  // first put then get
  if (src_data != NULL) {
    yac_interpolation_execute_put_frac(interp, src_data, src_frac_mask);
    // this should return at some point
    if (tgt_data == NULL)
      while (!yac_interpolation_execute_test(interp));
  }
  if (tgt_data != NULL) yac_interpolation_execute_get(interp, tgt_data);
  diff_count +=
    check_tgt_data(tgt_data, ref_tgt_data_raw, tgt_field_size, collection_size);

  return diff_count;
}

static int check_interpolation(
  struct interpolation * interp,
  double * src_data_raw, size_t num_src_fields, size_t src_field_size,
  double * ref_tgt_data_raw, size_t tgt_field_size, size_t collection_size) {

  // generate source input data pointers
  double *** src_data = NULL;
  if (src_data_raw != NULL) {
    src_data = xmalloc(collection_size * sizeof(src_data));
    double ** src_field_data =
      xmalloc(collection_size * num_src_fields * sizeof(*src_field_data));
    for (size_t i = 0; i < collection_size;
         ++i, src_field_data += num_src_fields) {
      src_data[i] = src_field_data;
      for (size_t j = 0; j < num_src_fields;
           ++j, src_data_raw += src_field_size)
        src_data[i][j] = src_data_raw;
    }
  }

  // generate target input data pointers
  double ** tgt_data = NULL;
  if (ref_tgt_data_raw != NULL) {
    double * tgt_data_raw =
      xmalloc(collection_size * tgt_field_size * sizeof(*tgt_data_raw));
    for (size_t i = 0; i < collection_size * tgt_field_size; ++i)
      tgt_data_raw[i] = UNSET;
    tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
    for (size_t i = 0; i < collection_size; ++i, tgt_data_raw += tgt_field_size)
      tgt_data[i] = tgt_data_raw;
  }

  // check results
  int diff_count =
    check_interpolation_execute(
      interp, src_data, tgt_data, ref_tgt_data_raw,
      tgt_field_size, collection_size);

  struct interpolation * interp_copy = yac_interpolation_copy(interp);
  diff_count +=
    check_interpolation_execute(
      interp, src_data, tgt_data, ref_tgt_data_raw,
      tgt_field_size, collection_size);
  yac_interpolation_delete(interp_copy);

  // clean up
  if (ref_tgt_data_raw != NULL) {
    free(tgt_data[0]);
    free(tgt_data);
  }
  if (src_data_raw != NULL) {
    free(src_data[0]);
    free(src_data);
  }

  return diff_count;
}

static int check_interpolation_frac(
  struct interpolation * interp,
  double * src_data_raw, double * src_frac_mask_raw,
  size_t num_src_fields, size_t src_field_size,
  double * ref_tgt_data_raw, size_t tgt_field_size, size_t collection_size) {

  // generate source input data pointers
  double *** src_data = NULL;
  if (src_data_raw != NULL) {
    src_data = xmalloc(collection_size * sizeof(src_data));
    double ** src_field_data =
      xmalloc(collection_size * num_src_fields * sizeof(*src_field_data));
    for (size_t i = 0; i < collection_size;
         ++i, src_field_data += num_src_fields) {
      src_data[i] = src_field_data;
      for (size_t j = 0; j < num_src_fields;
           ++j, src_data_raw += src_field_size)
        src_data[i][j] = src_data_raw;
    }
  }

  // generate source fractional mask data pointers
  double *** src_frac_data = NULL;
  if (src_frac_mask_raw != NULL) {
    src_frac_data = xmalloc(collection_size * sizeof(src_frac_data));
    double ** src_frac_mask_data =
      xmalloc(collection_size * num_src_fields * sizeof(*src_frac_mask_data));
    for (size_t i = 0; i < collection_size;
         ++i, src_frac_mask_data += num_src_fields) {
      src_frac_data[i] = src_frac_mask_data;
      for (size_t j = 0; j < num_src_fields;
           ++j, src_frac_mask_raw += src_field_size)
        src_frac_data[i][j] = src_frac_mask_raw;
    }
  }

  // apply fractional mask to source data
  if ((src_data != NULL) && (src_frac_data != NULL))
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_src_fields; ++j)
        for (size_t k = 0; k < src_field_size; ++k)
          src_data[i][j][k] *= src_frac_data[i][j][k];

  // generate target input data pointers
  double ** tgt_data = NULL;
  if (ref_tgt_data_raw != NULL) {
    double * tgt_data_raw =
      xmalloc(collection_size * tgt_field_size * sizeof(*tgt_data_raw));
    for (size_t i = 0; i < collection_size * tgt_field_size; ++i)
      tgt_data_raw[i] = UNSET;
    tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
    for (size_t i = 0; i < collection_size; ++i, tgt_data_raw += tgt_field_size)
      tgt_data[i] = tgt_data_raw;
  }

  // check results
  int diff_count =
    check_interpolation_execute_frac(
      interp, src_data, src_frac_data, tgt_data, ref_tgt_data_raw,
      tgt_field_size, collection_size);

  struct interpolation * interp_copy = yac_interpolation_copy(interp);
  diff_count +=
    check_interpolation_execute_frac(
      interp, src_data, src_frac_data, tgt_data, ref_tgt_data_raw,
      tgt_field_size, collection_size);
  yac_interpolation_delete(interp_copy);

  // clean up
  if (ref_tgt_data_raw != NULL) {
    free(tgt_data[0]);
    free(tgt_data);
  }
  if (src_frac_mask_raw != NULL) {
    free(src_frac_data[0]);
    free(src_frac_data);
  }
  if (src_data_raw != NULL) {
    free(src_data[0]);
    free(src_data);
  }

  return diff_count;
}
