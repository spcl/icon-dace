//#define VERBOSE

/**
 * @file test_clipping.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "grid.h"
#include "clipping.h"
#include "geometry.h"
#include "tests.h"
#include "area.h"
#include "test_common.h"

static void test_clipping(struct grid_cell cells[2],
                          struct grid_cell * ref_cells,
                          unsigned num_ref_cells);

static int compare_cells(struct grid_cell a, struct grid_cell b);

static void check_weight_correction(
  double * weights, size_t count, double * ref_weights);

int main (void) {

  enum yac_edge_type gc_edges[] = {
    GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
    GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE};
  enum yac_edge_type latlon_edges[] = {
    LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE};

  struct grid_cell overlap_cell[2];

  yac_init_grid_cell(overlap_cell);
  yac_init_grid_cell(overlap_cell+1);

#ifdef VERBOSE
  printf ("\n 1.) Simple test case with two quadrangle\n\n");
#endif
  /* 1.) Simple test case with two quadrangle

         source cell test data

              0.0 (lon) [2]
              0.7 (lat)
              / \
             /   \
   -0.5 [1] /     \  0.5  [3]
    0.0     \     /  0.0
             \   /
              \ /
              0.0  [0]
             -0.5
  */

  /* target cell test data

            0.0 (lon)  [0]
            0.6 (lat)
             / \
            /   \
  -0.6 [3] /     \  0.6  [1]
   0.0     \     /  0.0
            \   /
             \ /
             0.0  [2]
            -0.6
  */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, 0.0, 0.6, 0.6, 0.0,
                   GREAT_CIRCLE_EDGE, 0.0, 0.7, 0.5, 0.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0, 0.5, 0.0, -0.5},
                          (double[]){-0.5, 0.0, 0.7, 0.0}, gc_edges, 4),
        generate_cell_deg((double[]){0.0, 0.6, 0.0, -0.6},
                          (double[]){0.6, 0.0, -0.6, 0.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0, intersection_lon, 0.5, 0.0, -0.5, -intersection_lon},
                          (double[]){0.6, intersection_lat, 0.0, -0.5, 0.0, intersection_lat}, gc_edges, 6)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 2.a) Simple test case with one quadrangle and one triangle \n\n");
#endif
  /* 2.a) Simple test case with one quadrangle and one triangle

         source cell test data

              0.0 (lon) [2]
              0.7 (lat)
              / \
             /   \
   -0.5 [1] /     \  0.5  [3]
    0.0     \     /  0.0
             \   /
              \ /
              0.0  [0]
             -0.5
  */

  /* target cell test data

            0.0 (lon)
            0.6 (lat)
             / \
            /   \
  -0.6     /     \   0.6
   0.0    /_______\  0.0

  */

  {
    double intersection[3], intersection_lon, intersection_lat;
    if (!intersect(GREAT_CIRCLE_EDGE, 0.0, 0.6, 0.6, 0.0,
                   GREAT_CIRCLE_EDGE, 0.0, 0.7, 0.5, 0.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 0.0,-0.5,0.0,0.5},
                          (double[]){-0.5,0.0,0.7,0.0}, gc_edges, 4),
        generate_cell_deg((double[]){0.6,0.0,-0.6},
                          (double[]){0.0,0.6,0.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,intersection_lon,0.5,-0.5,-intersection_lon},
                          (double[]){0.6,intersection_lat,0.0,0.0,intersection_lat}, gc_edges, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 2.b) Simple test case with one quadrangle and one triangle \n\n");
#endif
  /* 2.b) Simple test case with one quadrangle and one triangle

         source cell test data

              0.0 (lon) [2]
              0.7 (lat)
              / \
             /   \
   -0.5 [1] /     \  0.5  [3]
    0.0     \     /  0.0
             \   /
              \ /
              0.0  [0]
             -0.5
  */

  /* target cell test data

            0.6 (lon)
            0.6 (lat)
             / \
            /   \
   0.0     /     \   1.0
   0.1    /_______\  0.1

  */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 0.0, 0.1, 1.0, 0.1,
                   GREAT_CIRCLE_EDGE, 0.0, 0.7, 0.5, 0.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 0.0, 0.1, 0.6, 0.6,
                   GREAT_CIRCLE_EDGE, 0.0, 0.7, 0.5, 0.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.5,0.0,0.5},
                          (double[]){-0.5,0.0,0.7,0.0}, gc_edges, 4),
        generate_cell_deg((double[]){1.0,0.6,0.0},
                          (double[]){0.1,0.6,0.1}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,intersection_lon[0],intersection_lon[1]},
                          (double[]){0.1,intersection_lat[0],intersection_lat[1]}, gc_edges, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 2.c) Simple test case with one quadrangle and one triangle \n\n");
#endif
  /* 2.c) Simple test case with one quadrangle and one triangle

         source cell test data

              0.0 (lon) [2]
              0.7 (lat)
              / \
             /   \
   -0.5 [1] /     \  0.5  [3]
    0.0     \     /  0.0
             \   /
              \ /
              0.0  [0]
             -0.5
  */

  /* target cell test data

            0.0 (lon)
            0.6 (lat)
             / \
            /   \
  -0.5     /     \   0.5
   0.1    /_______\  0.1

  */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,0.6,0.5,0.1,
                   GREAT_CIRCLE_EDGE, 0.0,0.7,0.5,0.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -0.5,0.1,0.5,0.1,
                   GREAT_CIRCLE_EDGE,  0.0,0.7,0.5,0.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.5,0.0,0.5},
                          (double[]){-0.5,0.0,0.7,0.0}, gc_edges, 4),
        generate_cell_deg((double[]){0.5,0.0,-0.5},
                          (double[]){0.1,0.6,0.1}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,intersection_lon[0],intersection_lon[1],-intersection_lon[1],-intersection_lon[0]},
                          (double[]){0.6,intersection_lat[0],intersection_lat[1],intersection_lat[1],intersection_lat[0]}, gc_edges, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 3.) Two source cells overlapping one target cell\n\n");
#endif

  /* 3.) Two source cells overlapping one target cell

     target cell test data

            0.0 (lon)
            0.7 (lat)
             / \
            1   0
  -0.5     /     \  0.5
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.7

  */

  /* source cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \   0.6
   0.0    /___2___\  0.0
          _________
          \   0   /
           \     /
            1   2
             \ /
             0.0
            -0.5
  */

  // generate reference intersection cells
  {
    double intersection[3], intersection_lon, intersection_lat;
    if (!intersect(GREAT_CIRCLE_EDGE, 0.0, 0.5, 0.6, 0.0,
                   GREAT_CIRCLE_EDGE, 0.0, 0.7, 0.5, 0.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    struct grid_cell ref_cells[] = {
      generate_cell_deg((double[]){0.0,intersection_lon,0.5,-0.5,-intersection_lon},
                        (double[]){0.5,intersection_lat,0.0,0.0,intersection_lat},
                        gc_edges, 5),
      generate_cell_deg((double[]){0.0,-intersection_lon,-0.5,0.5,intersection_lon},
                        (double[]){-0.5,-intersection_lat,0.0,0.0,-intersection_lat},
                        gc_edges, 5)};

    for (int src_order_a = -1; src_order_a <= 1; src_order_a += 2) {
      for (int src_order_b = -1; src_order_b <= 1; src_order_b += 2) {
        for (int tgt_order = -1; tgt_order <= 1; tgt_order += 2) {
          for (int src_start_a = 0; src_start_a < 3; ++src_start_a) {
            for (int src_start_b = 0; src_start_b < 3; ++src_start_b) {
              for (int tgt_start = 0; tgt_start < 4; ++tgt_start) {

                struct grid_cell SourceCells[] = {
                  generate_cell_deg((double[]){0.6,0.0,-0.6},
                                  (double[]){0.0,0.5,0.0}, gc_edges, 3),
                  generate_cell_deg((double[]){0.6,-0.6,0.0},
                                  (double[]){0.0,0.0,-0.5}, gc_edges, 3)};
                struct grid_cell TargetCell =
                  generate_cell_deg((double[]){0.5,0.0,-0.5,0.0},
                                  (double[]){0.0,0.7,0.0,-0.7}, gc_edges, 4);

                yac_cell_clipping ( 2, SourceCells, TargetCell, overlap_cell );
                if (compare_cells(overlap_cell[0], ref_cells[0]) ||
                    compare_cells(overlap_cell[1], ref_cells[1]))
                  PUT_ERR("ERROR: wrong clipping cell\n");
                yac_free_grid_cell(&SourceCells[0]);
                yac_free_grid_cell(&SourceCells[1]);
                yac_free_grid_cell(&TargetCell);
              }
            }
          }
        }
      }
    }
    yac_free_grid_cell(&ref_cells[0]);
    yac_free_grid_cell(&ref_cells[1]);
  }

#ifdef VERBOSE
  printf ("\n 4.a.1) Simple cases: identity \n\n");
#endif

  /* 4.a.1) Simple cases: identity

     target cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \  0.6
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.5

  */

  /* source cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \   0.6
   0.0    /___2___\  0.0
          _________
          \   0   /
           \     /
            1   2
             \ /
             0.0
            -0.5
  */

  // generate reference intersection cells
  {
    struct grid_cell ref_cells[] = {
      generate_cell_deg((double[]){0.0,0.6,-0.6},
                        (double[]){0.5,0.0,0.0}, gc_edges, 3),
      generate_cell_deg((double[]){0.0,0.6,-0.6},
                        (double[]){-0.5,0.0,0.0}, gc_edges, 3)};

    for (int src_order_a = -1; src_order_a <= 1; src_order_a += 2) {
      for (int src_order_b = -1; src_order_b <= 1; src_order_b += 2) {
        for (int tgt_order = -1; tgt_order <= 1; tgt_order += 2) {
          for (int src_start_a = 0; src_start_a < 3; ++src_start_a) {
            for (int src_start_b = 0; src_start_b < 3; ++src_start_b) {
              for (int tgt_start = 0; tgt_start < 4; ++tgt_start) {

                struct grid_cell SourceCells[] = {
                  generate_cell_deg((double[]){0.6,0.0,-0.6},
                                  (double[]){0.0,0.5,0.0}, gc_edges, 3),
                  generate_cell_deg((double[]){0.6,-0.6,0.0},
                                  (double[]){0.0,0.0,-0.5}, gc_edges, 3)};
                struct grid_cell TargetCell =
                  generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                                  (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4);

                yac_cell_clipping ( 2, SourceCells, TargetCell, overlap_cell );
                if (compare_cells(overlap_cell[0], ref_cells[0]) ||
                    compare_cells(overlap_cell[1], ref_cells[1]))
                  PUT_ERR("ERROR: wrong clipping cell\n");
                yac_free_grid_cell(overlap_cell);
                yac_free_grid_cell(overlap_cell+1);
                yac_free_grid_cell(&SourceCells[0]);
                yac_free_grid_cell(&SourceCells[1]);
                yac_free_grid_cell(&TargetCell);
              }
            }
          }
        }
      }
    }
    yac_free_grid_cell(&ref_cells[0]);
    yac_free_grid_cell(&ref_cells[1]);
  }

#ifdef VERBOSE
  printf ("\n 4.a.2) Simple cases: identity \n\n");
#endif

  /* 4.a.2) Simple cases: identity

     target cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \  0.6
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.5

  */

  /* source cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \  0.6
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.5
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                          (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4),
        generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                          (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                          (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 4.b) Simple cases: inside cell\n\n");
#endif

  /* 4.b) Simple cases: source inside target cell

     target cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \  0.6
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.5

  */

  /* source cell test data

            0.0 (lon)
            0.4 (lat)
             / \
            1   0
  -0.5     /     \  0.5
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.4
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.5,0.0,0.5},
                          (double[]){0.4,0.0,-0.4,0.0}, gc_edges, 4),
        generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                          (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.5,0.0,0.5},
                          (double[]){0.4,0.0,-0.4,0.0}, gc_edges, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 4.c) Simple cases: target inside source cell\n\n");
#endif

  /* 4.c) Simple cases: target inside source cell

     target cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \  0.6
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.5

  */

  /* source cell test data

            0.0 (lon)
            0.6 (lat)
             / \
            1   0
  -0.7     /     \  0.7
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.6
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.7,0.0,0.7},
                          (double[]){0.6,0.0,-0.6,0.0}, gc_edges, 4),
        generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                          (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                          (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 4.d) Simple cases: degenerated cell \n\n");
#endif

  /* 4.d) Simple cases: degenerated cell

     target cell test data

            0.0 (lon)
            0.5 (lat)
             / \
            1   0
  -0.6     /     \  0.6
   0.0     \     /  0.0
            2   3
             \ /
             0.0
            -0.5

  */

  /* source cell test data

       0.0     0.0 (lon)
       0.5     0.6 (lat)
        --------
         \      \
          \      \
           \      \
            --------
           0.6     1.0
           0.0     0.0
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,0.6,1.0,0.5},
                          (double[]){0.5,0.0,0.0,0.5}, gc_edges, 4),
        generate_cell_deg((double[]){0.0,-0.6,0.0,0.6},
                          (double[]){0.5,0.0,-0.5,0.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 5.) Special case from ICON_toy \n\n");
#endif
  /* 5.) source cell test data

     +------+
     |      |
     |      |
     |      |
     +------+
  */

  /* target cell test data

     /|
    / |
   /  |
  /   |
  \   |
   \  |
    \ |
     \|
  */

  {
    double intersection[4][3], intersection_lon[4], intersection_lat[4];

    if (!intersect(GREAT_CIRCLE_EDGE, 108.0000, -7.2167, 88.4979, -19.2921,
                   GREAT_CIRCLE_EDGE, 90.0, 0.0, 90.0, -22.5, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 88.4979, -19.2921, 108.0000, -26.5650,
                   GREAT_CIRCLE_EDGE, 90.0, 0.0, 90.0, -22.5, intersection[1]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 88.4979, -19.2921, 108.0000, -26.5650,
                   GREAT_CIRCLE_EDGE, 90.0, -22.5, 112.5, -22.5, intersection[2]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 108.0000, -26.5650, 108.0000, -7.2167,
                   GREAT_CIRCLE_EDGE, 90.0, -22.5, 112.5, -22.5, intersection[3]))
      return EXIT_FAILURE;
    for (int i = 0; i < 4; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){90.0,112.5,112.5,90.0},
                          (double[]){-22.5,-22.5,0.0,0.0}, gc_edges, 4),
        generate_cell_deg((double[]){108.0000,88.4979,108.0000},
                          (double[]){-7.2167,-19.2921,-26.5650}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){108.0000,
                                   intersection_lon[0],
                                   intersection_lon[1],
                                   intersection_lon[2],
                                   intersection_lon[3]},
                          (double[]){-7.2167,
                                   intersection_lat[0],
                                   intersection_lat[1],
                                   intersection_lat[2],
                                   intersection_lat[3]}, gc_edges, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 6.) Special case from ICON_toy \n\n");
#endif
  /* 6.) source cell test data

     +------+  0.5
     |      |
     |      |
     |      |
     +------+ -0.5
   -0.5    0.5
  */

  /* target cell test data

    -0.75   0.0
        ______  0.0
        \    |
         \   |
          \  |
           \ |
            \| -0.75
  */

  {
    double intersection[4][3], intersection_lon[4], intersection_lat[4];

    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,0.0,-0.75,0.0,
                   GREAT_CIRCLE_EDGE, -0.5,0.5,-0.5,-0.5, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -0.75,0.0,0.0,-0.75,
                   GREAT_CIRCLE_EDGE, -0.5,0.5,-0.5,-0.5, intersection[1]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -0.75,0.0,0.0,-0.75,
                   GREAT_CIRCLE_EDGE, -0.5,-0.5,0.5,-0.5, intersection[2]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,-0.75,0.0,0.0,
                   GREAT_CIRCLE_EDGE, -0.5,-0.5,0.5,-0.5, intersection[3]))
      return EXIT_FAILURE;
    for (int i = 0; i < 4; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-0.5,0.5,0.5,-0.5},
                          (double[]){-0.5,-0.5,0.5,0.5}, gc_edges, 4),
        generate_cell_deg((double[]){0.0,-0.75,0.0},
                          (double[]){0.0,0.0,-0.75}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,
                                     intersection_lon[0],
                                     intersection_lon[1],
                                     intersection_lon[2],
                                     intersection_lon[3]},
                          (double[]){0.0,
                                     intersection_lat[0],
                                     intersection_lat[1],
                                     intersection_lat[2],
                                     intersection_lat[3]}, gc_edges, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 7.) Special case in which no corner on any square is within the other \n\n");
#endif
  /* 7.) source cell test data

     +------+  0.5
     |      |
     |      |
     |      |
     +------+ -0.5
   -0.5    0.5
  */

  /* target cell test data

   -0.7  0.0  0.7
          /\     0.7
         /  \
        /    \
       /      \  0.0
       \      /
        \    /
         \  /
          \/    -0.7
  */

  {
    double intersection[3], intersection_lon[8], intersection_lat[8];
    double src_data[2][4] = {{-0.5,-0.5, 0.5, 0.5},
                           { 0.5,-0.5,-0.5, 0.5}};
    double tgt_data[2][4] = {{ 0.0,-0.7, 0.0, 0.7},
                           { 0.7, 0.0,-0.7, 0.0}};

    for (int i = 0; i < 8; ++i) {
      if (!intersect(GREAT_CIRCLE_EDGE, src_data[0][i/2],
                                 src_data[1][i/2],
                                 src_data[0][(i/2+1)%4],
                                 src_data[1][(i/2+1)%4],
                     GREAT_CIRCLE_EDGE, tgt_data[0][((i+1)/2)%4],
                                 tgt_data[1][((i+1)/2)%4],
                                 tgt_data[0][((i+1)/2+1)%4],
                                 tgt_data[1][((i+1)/2+1)%4], intersection))
        return EXIT_FAILURE;
      XYZtoLL(intersection, &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-0.5,0.5,0.5,-0.5},
                          (double[]){-0.5,-0.5,0.5,0.5}, gc_edges, 4),
        generate_cell_deg((double[]){0.0,0.7,0.0,-0.7},
                          (double[]){-0.7,0.0,0.7,0.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(intersection_lon, intersection_lat, gc_edges, 8)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 8.1) Test non intersecting cells \n\n");
#endif
  /* 8.) source cell test data

     +------+ 1
     |      |
     |      |
     |      |
     +------+ 0
     0      1
  */

  /* target cell test data

     +------+ 1
     |      |
     |      |
     |      |
     +------+ 0
     2      3
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,1.0,1.0,0.0},
                          (double[]){0.0,0.0,1.0,1.0}, gc_edges, 4),
        generate_cell_deg((double[]){2.0,3.0,3.0,2.0},
                          (double[]){0.0,0.0,1.0,1.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 8.2) Test non intersecting cells \n\n");
#endif
  /* 8.) source cell test data

     +------+ 1
     |      |
     |      |
     |      |
     +------+ 0
     0      1
  */

  /* target cell test data

     +------+ 3
     |      |
     |      |
     |      |
     +------+ 2
     1      2
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,1.0,1.0,0.0},
                          (double[]){0.0,0.0,1.0,1.0}, gc_edges, 4),
        generate_cell_deg((double[]){1.0,2.0,2.0,1.0},
                          (double[]){2.0,2.0,3.0,3.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 9.) triangle-square intersection \n\n");
#endif
  /* 9.) source cell test data

    -1.5     1.5
      --------  0.5
      \      /
       \    /
        \  /
         \/    -0.5
  */

  /* target cell test data

     +------+  1
     |      |
     |      |
     |      |
     +------+ -1
    -1      1
  */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 1.0,-1.0,1.0,1.0,
                   GREAT_CIRCLE_EDGE, 0.0,-0.5,1.5,0.5, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 1.0,-1.0,1.0,1.0,
                   GREAT_CIRCLE_EDGE, 1.5,0.5,-1.5,0.5, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){1.5,-1.5,0.0},
                          (double[]){0.5,0.5,-0.5}, gc_edges, 3),
        generate_cell_deg((double[]){-1.0,1.0,1.0,-1.0},
                          (double[]){-1.0,-1.0,1.0,1.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,
                                   intersection_lon[0],
                                   intersection_lon[1],
                                   -intersection_lon[1],
                                   -intersection_lon[0]},
                          (double[]){-0.5,
                                   intersection_lat[0],
                                   intersection_lat[1],
                                   intersection_lat[1],
                                   intersection_lat[0]}, gc_edges, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 10.1) touching edges of two triangles \n\n");
#endif
  /* 10.1) source cell test data

       /\     1
      /  \
     /    \
    /      \
    --------  0
    0      2
  */
  /* target cell test data

      --------  0.0
      \      /
       \    /
        \  /
         \/    -0.5
     0  0.5  1
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,2.0,1.0},
                          (double[]){0.0,0.0,1.0}, gc_edges, 3),
        generate_cell_deg((double[]){0.0,1.0, 0.5},
                          (double[]){0.0,0.0,-0.5}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 10.2) touching edges \n\n");
#endif

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){160.0,180.0,180.0},
                          (double[]){ 65.0, 65.0, 90.0}, gc_edges, 3),
        generate_cell_deg((double[]){180.0,200.0,  0.0},
                          (double[]){ 70.0, 70.0, 90.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 11.1) test from test_interpolation_method_conserv.x \n\n");
#endif
  /* 11.) source cell test data

     +------+  0
     |      |
     |      |
     |      |
     +------+ -1
     0      1
  */
  /* target cell test data

     +------+  0.5
     |      |
     |      |
     |      |
     +------+ -0.5
   -0.5    0.5
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 0.0, 1.0,1.0,0.0},
                          (double[]){-1.0,-1.0,0.0,0.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-0.5, 0.5,0.5,-0.5},
                          (double[]){-0.5,-0.5,0.5, 0.5}, latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,0.5, 0.5, 0.0},
                          (double[]){0.0,0.0,-0.5,-0.5}, latlon_edges, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 11.2) test from test_interpolation_method_conserv.x \n\n");
#endif
  /* 11.) source cell test data

    -1      0
     +------+  0
     |      |
     |      |
     |      |
     +------+ -1
  */
  /* target cell test data

     +------+  0.5
     |      |
     |      |
     |      |
     +------+ -0.5
   -0.5    0.5
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 0.0,-1.0,-1.0,0.0},
                          (double[]){-1.0,-1.0, 0.0,0.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-0.5, 0.5,0.5,-0.5},
                          (double[]){-0.5,-0.5,0.5, 0.5}, latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,-0.5,-0.5, 0.0},
                          (double[]){0.0, 0.0,-0.5,-0.5}, latlon_edges, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.1) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

        --------  39
        \      /
         \    /
          \  /
           \/     22
     -20   0   20
    */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(LAT_CIRCLE_EDGE,   -20.0,40.0, 0.0,40.0,
                   GREAT_CIRCLE_EDGE, -20.0,39.0,20.0,39.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0, 60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){ 39.0,22.0,39.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon,-intersection_lon},
                          (double[]){intersection_lat,intersection_lat},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 2),
        generate_cell_deg((double[]){intersection_lon,-intersection_lon},
                          (double[]){intersection_lat,intersection_lat},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 2)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 12.2) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

  40  --....__  39
      \      /
       \    /
        \  /
         \/     22
   -20   0   20
  */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(LAT_CIRCLE_EDGE,   20.0,40.0,-19.0,40.0,
                   GREAT_CIRCLE_EDGE, 20.0,39.0,-20.0,40.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0, 60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){ 40.0,22.0,39.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,intersection_lon},
                          (double[]){ 40.0,intersection_lat},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 2),
        generate_cell_deg((double[]){-20.0,intersection_lon},
                          (double[]){ 40.0,intersection_lat},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 2)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 12.3) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

         /\     80
        /  \
       /    \
      /      \
      --------  59
   -19   0   19
  */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE,  19.0,59.0, 0.0,80.0,
                   LAT_CIRCLE_EDGE, -20.0,60.0,20.0,60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 19.0, 59.0,-19.0, 59.0,
                   LAT_CIRCLE_EDGE, -20.0, 60.0, 0.0, 60.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0, 60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-19.0, 0.0,19.0},
                          (double[]){ 59.0,80.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){19.0,intersection_lon[0],
                                   fabs(intersection_lon[1]),
                                   -fabs(intersection_lon[1]),
                                   -intersection_lon[0],-19.0},
                          (double[]){59.0,60.0,60.0,60.0,60.0,59.0},
                          (enum yac_edge_type[]){
                          GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                          LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 6),
        generate_cell_deg((double[]){19.0,intersection_lon[0],-intersection_lon[0],-19.0,-fabs(intersection_lon[1]),fabs(intersection_lon[1])},
                          (double[]){59.0,60.0,60.0,59.0,60.0,60.0},
                          (enum yac_edge_type[]){
                          GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                          GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 6)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 12.4) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

            /| 80
          /  |
        /    |
      /      |
  60  --...__|
               59
   -20       21
  */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, -20.0,60.0,21.0,59.0,
                   LAT_CIRCLE_EDGE,   -19.0,60.0,20.0,60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -20.0,60.0,21.0,59.0,
                   LON_CIRCLE_EDGE,    20.0,40.0,20.0,60.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0, 60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0,21.0,21.0},
                          (double[]){ 60.0,80.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){20.0,intersection_lon[0],intersection_lon[1]},
                          (double[]){60.0,intersection_lat[0],intersection_lat[1]},
                          (enum yac_edge_type[]){
                          LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 3),
        generate_cell_deg((double[]){20.0,intersection_lon[0],-20.0,intersection_lon[1]},
                          (double[]){60.0,intersection_lat[0], 60.0,intersection_lat[1]},
                          (enum yac_edge_type[]){
                          LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 4)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 12.5) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

 60 --...__  59
    \      |
      \    |
        \  |
          \| 40

  -20     21

  */

  {
    double intersection[3][3], intersection_lon[3], intersection_lat[3];

    if (!intersect(GREAT_CIRCLE_EDGE, -20.0,60.0,21.0,59.0,
                   LAT_CIRCLE_EDGE,   -19.0,60.0,20.0,60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -20.0,60.0,21.0,59.0,
                   LON_CIRCLE_EDGE,    20.0,40.0,20.0,60.0, intersection[1]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -20.0,60.0,21.0,40.0,
                   LON_CIRCLE_EDGE,    20.0,40.0,20.0,60.0, intersection[2]))
      return EXIT_FAILURE;
    for (int i = 0; i < 3; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0, 60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0,21.0,21.0},
                          (double[]){ 60.0,40.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon[0],intersection_lon[1],intersection_lon[2],-20.0},
                          (double[]){intersection_lat[0],intersection_lat[1],intersection_lat[2], 60.0},
                          (enum yac_edge_type[]){
                          GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.6) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
       x      x+40
    */

    /* cell b

        --------  39
        \      /
         \    /
          \  /
           \/     20
      -20  0   20
    */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, -20.0,39.0,20.0,39.0,
                   LAT_CIRCLE_EDGE,   -20.0,40.0, 0.0,40.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon,intersection_lon+40.0,
                                     intersection_lon+40.0,intersection_lon},
                          (double[]){40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0,20.0, 0.0},
                          (double[]){ 39.0,39.0,20.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon,-intersection_lon},
                          (double[]){40.0,40.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 2),
        generate_cell_deg((double[]){intersection_lon,-intersection_lon},
                          (double[]){40.0,40.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 2)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 12.7) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
       x      x+40
    */

    /* cell b

      |\        50
      |  \
      |    \
      |      \
      --------  39
   -20   0   20
    */

  {
    double intersection[3][3], intersection_lon[3], intersection_lat[3];

    if (!intersect(GREAT_CIRCLE_EDGE, -20.0,39.0,20.0,39.0,
                   LAT_CIRCLE_EDGE,   -20.0,40.0, 0.0,40.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE,  20.0,39.0,-20.0,50.0,
                   LAT_CIRCLE_EDGE,   -20.0,40.0, 20.0,40.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }
    if (!intersect(GREAT_CIRCLE_EDGE, 20.0,39.0,-20.0,50.0,
                   LON_CIRCLE_EDGE,   intersection_lon[0],40.0,
                                 intersection_lon[0],60.0, intersection[2]))
      return EXIT_FAILURE;
    XYZtoLL(intersection[2], &intersection_lon[2], &intersection_lat[2]);
    intersection_lon[2] /= YAC_RAD;
    intersection_lat[2] /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon[0],
                                     intersection_lon[0]+40.0,
                                     intersection_lon[0]+40.0,
                                     intersection_lon[0]},
                          (double[]){40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0,20.0,-20.0},
                          (double[]){ 39.0,39.0, 50.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon[0],-intersection_lon[0],
                                   intersection_lon[1],intersection_lon[2]},
                          (double[]){40.0,40.0,40.0,intersection_lat[2]},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.8) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

        --------  40
        \      /
         \    /
          \  /
           \/     20
     -20   0   20
    */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0, 60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){ 40.0,20.0,40.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0},
                          (double[]){ 40.0,40.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 2),
        generate_cell_deg((double[]){-20.0,20.0},
                          (double[]){ 40.0,40.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 2)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
    printf ("\n 12.9) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

         /\     50
        /  \
       /    \
      /      \
      --------  40
   -20   0   20
    */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){ 40.0,50.0,40.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,0.0},
                          (double[]){40.0,40.0,50.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.10) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

        --------  60
        \      /
         \    /
          \  /
           \/     50
     -20   0   20
    */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){ 60.0,50.0,60.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,0.0},
                          (double[]){60.0,60.0,50.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.11) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

         /\     70
        /  \
       /    \
      /      \
      --------  60
   -20   0   20
    */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){ 60.0,70.0,60.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0},
                          (double[]){60.0,60.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 2),
        generate_cell_deg((double[]){-20.0,20.0},
                          (double[]){60.0,60.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 2)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
    printf ("\n 12.12) some mean tests \n\n");
#endif
    /* cell a

       +------+ -40
       |      |
       |      |
       |      |
       +------+ -60
      -20     20
    */

    /* cell b

        -------- -40
        \      /
         \    /
          \  /
           \/    -50
     -20   0   20
    */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){-40.0,-40.0,-60.0,-60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){-40.0,-50.0,-40.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,0.0},
                          (double[]){-40.0,-40.0,-50.0}, gc_edges, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.13) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+  40
      -20     20
    */

    /* cell b

     59 ___..--- 60
        \      /
         \    /
          \  /
           \/    45
     -15   0   15
    */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(LAT_CIRCLE_EDGE,   -20.0,60.0,14.9,60.0,
                   GREAT_CIRCLE_EDGE, -15.0,59.0,15.0,60.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){60.0,60.0,40.0,40.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-15.0, 0.0,15.0},
                          (double[]){59.0,45.0,60.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon,15.0,0.0,-15.0},
                          (double[]){intersection_lat,60.0,45.0,59.0},
                          (enum yac_edge_type[]){
                              LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                              GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.14) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       +------+  59
      -20     20
    */

    /* cell b

        --------  58.5
        \      /
         \    /
          \  /
           \/     40
     -20   0   20
    */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(LAT_CIRCLE_EDGE,   -20.0, 59.0, 0.0, 59.0,
                   GREAT_CIRCLE_EDGE, -20.0, 58.5,20.0, 58.5, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(LAT_CIRCLE_EDGE,   -20.0, 60.0, 0.0, 60.0,
                   GREAT_CIRCLE_EDGE, -20.0, 58.5,20.0, 58.5, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){60.0,60.0,59.0,59.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){58.5,40.0,58.5}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ intersection_lon[0],
                                      intersection_lon[1],
                                     -intersection_lon[1],
                                     -intersection_lon[0]},
                          (double[]){59.0,60.0,60.0,59.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.15) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

         /\     80
        /  \
       /    \
      /      \
      --------  59
   -20   0   20
  */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 20.0,59.0,0.0,80.0,
                   LAT_CIRCLE_EDGE,  -20.0, 60.0,20.0,60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 20.0, 59.0,-20.0, 59.0,
                   LAT_CIRCLE_EDGE,  -20.0, 60.0, 0.0, 60.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){59.0,80.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){20.0,intersection_lon[0],
                                     fabs(intersection_lon[1]),
                                     -fabs(intersection_lon[1]),
                                     -intersection_lon[0],-20.0},
                          (double[]){59.0,60.0,60.0,60.0,60.0,59.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 6),
        generate_cell_deg((double[]){20.0,intersection_lon[0],
                                     -intersection_lon[0],-20.0,
                                     -fabs(intersection_lon[1]),
                                     fabs(intersection_lon[1])},
                          (double[]){59.0,60.0,60.0,59.0,60.0,60.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 6)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 12.16) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

        --------  59
        \      /
         \    /
          \  /
           \/     50
     -20   0   20
  */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, 20.0,59.0,-20.0,59.0,
                   LAT_CIRCLE_EDGE,  -20.0,60.0,0.0,60.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0, 0.0,20.0},
                          (double[]){59.0,50.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon,-intersection_lon,
                                     20.0,0.0,-20.0},
                          (double[]){60.0,60.0,59.0,50.0,59.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.17) some mean tests \n\n");
#endif
  /* cell a

     +------+  20
     |      |
     |      |
     |      |
     +------+ -20
    -20     20
  */

  /* cell b

     +------+  19
     |      |
     |      |
     |      |
     +------+ -19
    -19     19
  */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, -19.0,19.0,19.0,19.0,
                   LAT_CIRCLE_EDGE,   -20.0,20.0, 0.0,20.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 20.0,20.0,-20.0,-20.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-19.0,-19.0,19.0,19.0},
                          (double[]){ 19.0,-19.0,-19.0,19.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-19.0,intersection_lon,
                                     -intersection_lon,19.0,19.0,
                                     -intersection_lon,intersection_lon,-19.0},
                          (double[]){19.0,20.0,20.0,19.0,
                                     -19.0,-20.0,-20.0,-19.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 8)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.18) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

         /\     55
        /  \
       /    \
      /      \
      --------  45
   -40  -20   0
    */

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, -40.0, 45.0, 0.0, 45.0,
                   LON_CIRCLE_EDGE,   -20.0, 40.0,-20.0, 60.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-40.0, 0.0,-20.0},
                          (double[]){ 45.0,45.0,55.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,-20.0,0.0},
                          (double[]){intersection_lat,55.0,45.0},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.20) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

            /| 50
          /  |
        /    | 40
        \    |
          \  |
            \| 30
      -10   10
    */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-10.0,10.0,10.0},
                          (double[]){ 40.0,30.0,50.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-10.0,10.0,10.0},
                          (double[]){40.0,50.0,40.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.21) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

            /| 80
          /  |
        /    |
      /      |
  61  --...__|
               59
   -20       21
  */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, -20.0, 61.0,21.0, 59.0,
                   LAT_CIRCLE_EDGE,   -20.0, 60.0,20.0, 60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -20.0, 61.0,21.0, 59.0,
                   LON_CIRCLE_EDGE,    20.0,40.0,20.0,60.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0,21.0,21.0},
                          (double[]){61.0,80.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){20.0,intersection_lon[0],
                                     intersection_lon[1]},
                          (double[]){60.0,intersection_lat[0],
                                     intersection_lat[1]},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.22) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

         /\     x (slightly above 40)
        /  \
       /    \
      /      \
      --------  30
   -20   0   20
    */

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(LON_CIRCLE_EDGE,     0.0,60.0,0.0,30.0,
                   GREAT_CIRCLE_EDGE, -20.0, 40.0,20.0, 40.0, intersection[0]))
      return EXIT_FAILURE;
    XYZtoLL(intersection[0], &intersection_lon[0], &intersection_lat[0]);
    intersection_lon[0] /= YAC_RAD;
    intersection_lat[0] /= YAC_RAD;
    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,intersection_lat[0],20.0,30.0,
                   LAT_CIRCLE_EDGE, -20.0, 40.0,20.0, 40.0, intersection[1]))
      return EXIT_FAILURE;
    XYZtoLL(intersection[1], &intersection_lon[1], &intersection_lat[1]);
    intersection_lon[1] /= YAC_RAD;
    intersection_lat[1] /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-20.0,20.0, 0.0},
                          (double[]){30.0,30.0,intersection_lat[0]},
                          gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0,intersection_lon[1],
                                     -intersection_lon[1]},
                          (double[]){intersection_lat[0],40.0,40.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
    printf ("\n 12.23) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -20     20
    */

    /* cell b

         /\     70
        /  \
       /    \
      /      \
      --------  59
   -21   0   21
    */

  {
    double intersection[4][3], intersection_lon[4], intersection_lat[4];

    if (!intersect(GREAT_CIRCLE_EDGE, -21.0, 59.0,21.0, 59.0,
                   LON_CIRCLE_EDGE,    20.0,40.0,20.0,60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,70.0,21.0,59.0,
                   LON_CIRCLE_EDGE,  20.0,40.0,20.0,60.0, intersection[1]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,70.0,21.0,59.0,
                   LAT_CIRCLE_EDGE, -20.0, 60.0,20.0, 60.0, intersection[2]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE,-21.0, 59.0,21.0, 59.0,
                   LAT_CIRCLE_EDGE,    0.0,60.0,20.0,60.0, intersection[3]))
      return EXIT_FAILURE;
    for (int i = 0; i < 4; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-21.0, 0.0,21.0},
                          (double[]){59.0,70.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){20.0,20.0,intersection_lon[2],
                                     intersection_lon[3],-intersection_lon[3],
                                     -intersection_lon[2],-20.0,-20.0},
                          (double[]){intersection_lat[0],intersection_lat[1],
                                     60.0,60.0,60.0,60.0,intersection_lat[1],
                                     intersection_lat[0]},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 8),
        generate_cell_deg((double[]){20.0,20.0,intersection_lon[2],
                                     intersection_lon[3],-intersection_lon[3],
                                     -intersection_lon[2],-20.0,-20.0,
                                     -intersection_lon[3],intersection_lon[3]},
                          (double[]){intersection_lat[0],intersection_lat[1],
                                     60.0,60.0,60.0,60.0,intersection_lat[1],
                                     intersection_lat[0],60.0,60.0},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE}, 10),
        generate_cell_deg((double[]){20.0,20.0,intersection_lon[2],
                                     -intersection_lon[2],-20.0,-20.0,
                                     -intersection_lon[3],intersection_lon[3]},
                          (double[]){intersection_lat[0],intersection_lat[1],
                                     60.0,60.0,intersection_lat[1],
                                     intersection_lat[0],60.0,60.0},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 8)},
      // number of reference cells
      3);
  }

#ifdef VERBOSE
    printf ("\n 12.24) some mean tests \n\n");
#endif
    /* cell a

       +------+  60
       |      |
       |      |
       |      |
       +------+ 40
      -x      x      (x is choosen such that the upper corners of the
                      quadrangle touch the sloping edges of the triangle)
    */

    /* cell b

         /\     70
        /  \
       /    \
      /      \
      --------  59
   -21   0   21
    */

  {
    double intersection[3][3], intersection_lon[3], intersection_lat[3];

    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,70.0,21.0,59.0,
                   LAT_CIRCLE_EDGE,   0.0,60.0,25.0,60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE,-21.0,59.0,21.0,59.0,
                   LAT_CIRCLE_EDGE,    0.0,60.0,25.0,60.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }
    if (!intersect(GREAT_CIRCLE_EDGE, -21.0, 59.0,21.0, 59.0,
                   LON_CIRCLE_EDGE,   intersection_lon[0],60.0,intersection_lon[0],
                                 40.0, intersection[2]))
      return EXIT_FAILURE;
    XYZtoLL(intersection[2], &intersection_lon[2], &intersection_lat[2]);
    intersection_lon[2] /= YAC_RAD;
    intersection_lat[2] /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-intersection_lon[0],intersection_lon[0],
                                     intersection_lon[0],-intersection_lon[0]},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-21.0, 0.0,21.0},
                          (double[]){59.0,70.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon[2],intersection_lon[0],
                                     intersection_lon[1],-intersection_lon[1],
                                     -intersection_lon[0],-intersection_lon[2]},
                          (double[]){intersection_lat[2],60.0,60.0,60.0,60.0,
                                     intersection_lat[2]},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 6),
        generate_cell_deg((double[]){intersection_lon[2],intersection_lon[0],
                                     intersection_lon[1],-intersection_lon[1],
                                     -intersection_lon[0],-intersection_lon[2],
                                     -intersection_lon[1],intersection_lon[1]},
                          (double[]){intersection_lat[2],60.0,60.0,60.0,60.0,
                                     intersection_lat[2],60.0,60.0},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,LAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 8),
        generate_cell_deg((double[]){intersection_lon[2],intersection_lon[0],
                                     -intersection_lon[0],-intersection_lon[2],
                                     -intersection_lon[1],intersection_lon[1]},
                          (double[]){intersection_lat[2],60.0,60.0,
                                     intersection_lat[2],60.0,60.0},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,LAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 6)},
      // number of reference cells
      3);
  }

#ifdef VERBOSE
  printf ("\n 12.25) some mean tests \n\n");
#endif
  /* cell a

     +------+  60
     |      |
     |      |
     |      |
     +------+ 40
    -20     20
  */

  /* cell b

        --------  59
        \      /
         \    /
          \  /
           \/     50
     -21   0   21
  */

  {
    double intersection[3][3], intersection_lon[3], intersection_lat[3];

    if (!intersect(GREAT_CIRCLE_EDGE, 21.0,59.0,-21.0,59.0,
                   LAT_CIRCLE_EDGE,   0.0,60.0,20.0,60.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 21.0, 59.0,-21.0, 59.0,
                   LON_CIRCLE_EDGE,   20.0,40.0,20.0,60.0, intersection[1]))
      return EXIT_FAILURE;

    if (!intersect(GREAT_CIRCLE_EDGE, 0.0,50.0,21.0,59.0,
                   LON_CIRCLE_EDGE,  20.0,40.0,20.0,60.0, intersection[2]))
      return EXIT_FAILURE;
    for (int i = 0; i < 3; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-20.0,20.0,20.0,-20.0},
                          (double[]){ 40.0,40.0,60.0,60.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-21.0, 0.0,21.0},
                          (double[]){59.0,50.0,59.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon[0],20.0,20.0,0.0,-20.0,
                                     -20.0,-intersection_lon[0]},
                          (double[]){60.0,intersection_lat[1],
                                     intersection_lat[2],50.0,
                                     intersection_lat[2],intersection_lat[1],
                                     60.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE}, 7)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.26) some mean tests \n\n");
#endif
  /* cell a (circle around a pole)

         .--.
        /    \
        \    /
         '--'
  */

  /* cell b (circle around a pole; bigger than cell a)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  {

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5),
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 80.0,80.0,80.0,80.0,80.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.27) some mean tests \n\n");
#endif
  /* cell a (circle around a pole)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  /* cell b (circle around a pole; same size as cell a)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  {

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5),
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.28) some mean tests \n\n");
#endif
  /* cell a (circle around a pole)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  /* cell b (circle around different pole than cell a)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  {

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5),
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){-85.0,-85.0,-85.0,-85.0,-85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.29) some mean tests \n\n");
#endif
  /* cell a (circle around a pole; made of lat circle edges)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  /* cell b (circle around a pole; made of great circle edges;
             corners are on the lat circles)

        ____
       /    \
      /      \
     *        *
      \      /
       \____/
  */

  {

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5),
        generate_cell_deg((double[]){  0.0,60.0,120.0,180.0,240.0,300.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}, 6)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,60.0,120.0,180.0,240.0,300.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 6)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 12.30) some mean tests \n\n");
#endif
  /* cell a (circle around a pole; made of lat circle edges)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  /* cell b (circle around a pole; made of great circle edges;
             all great circles edges intersect with the lat circles twice)

        ____
       /    \
      /      \
     *        *
      \      /
       \____/
  */

  {

    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(LAT_CIRCLE_EDGE,   0.0,85.5,30.0,85.5,
                   GREAT_CIRCLE_EDGE, 0.0,85.0,60.0,85.0, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(LAT_CIRCLE_EDGE,  30.0,85.5,60.0,85.5,
                   GREAT_CIRCLE_EDGE, 0.0,85.0,60.0,85.0, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,72.0,144.0,216.0,288.0},
                          (double[]){ 85.5,85.5,85.5,85.5,85.5},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 5),
        generate_cell_deg((double[]){  0.0,60.0,120.0,180.0,240.0,300.0},
                          (double[]){ 85.0,85.0,85.0,85.0,85.0,85.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                          GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}, 6)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon[0]+60.0*0.0,
                                     intersection_lon[1]+60.0*0.0,
                                     intersection_lon[0]+60.0*1.0,
                                     intersection_lon[1]+60.0*1.0,
                                     intersection_lon[0]+60.0*2.0,
                                     intersection_lon[1]+60.0*2.0,
                                     intersection_lon[0]+60.0*3.0,
                                     intersection_lon[1]+60.0*3.0,
                                     intersection_lon[0]+60.0*4.0,
                                     intersection_lon[1]+60.0*4.0,
                                     intersection_lon[0]+60.0*5.0,
                                     intersection_lon[1]+60.0*5.0},
                          (double[]){85.5,85.5,85.5,85.5,85.5,85.5,85.5,85.5,
                                     85.5,85.5,85.5,85.5},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE},
                          12)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 13.) example taken from bug report by Uwe \n\n");
#endif

  /* 13.) source cell test data

     +------+ 5
     |      |
     |      |
     |      |
     +------+ 0
    355    360
  */

  /* target cell test data

     +------+   0
     |      |
     |      |
     |      |
     +------+ -10
    240    250
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){355.0,360.0,360.0,355.0},
                          (double[]){  0.0, 0.0, 5.0, 5.0}, gc_edges, 4),
        generate_cell_deg((double[]){240.0,250.0,250.0,240.0},
                          (double[]){-10.0,-10.0, 0.0, 0.0}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.1) example taken from icon toy \n\n");
#endif

  /* 14.1) source cell test data

     +------+ 2
     |      |
     |      |
     |      |
     +------+ 0
     0      2
  */

  /* target cell test data

     _______
    |      /  1
    |    /
    |  /
    |/       -1
    2      4
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 2.0,2.0,4.0},
                          (double[]){-1.0,1.0,1.0}, gc_edges, 3),
        generate_cell_deg((double[]){0.0,2.0,2.0,0.0},
                          (double[]){0.0,0.0,2.0,2.0}, latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.2) example taken from icon toy \n\n");
#endif

  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(LON_CIRCLE_EDGE,   1.0,90.0,1.0,80.0,
                   GREAT_CIRCLE_EDGE, 0.0,85.0,5.0,85.0, intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){  0.0,180.0, 5.0},
                          (double[]){85.0,90.0,85.0}, gc_edges, 3),
        generate_cell_deg((double[]){ 1.0,-5.0,-5.0,1.0},
                          (double[]){80.0,80.0,90.0,90.0}, latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){1.0,0.0,0.0},
                          (double[]){intersection_lat,85.0,90.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.3) example taken from icon toy \n\n");
#endif

  {

    double intersection[4][3], intersection_lon[4], intersection_lat[4];

    if (!intersect(LAT_CIRCLE_EDGE,   0.46633015951723494341,-0.88357293382212931387,
                                      0.49087385212340517437,-0.88357293382212931387,
                   GREAT_CIRCLE_EDGE, 0.47953965147452792817,-0.90109638077934106626,
                                      0.46519399935013672209,-0.88122392525212978054,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(LON_CIRCLE_EDGE,   0.46633015951723494341,-0.88357293382212931387,
                                      0.46633015951723494341,-0.85902924121595902740,
                   GREAT_CIRCLE_EDGE, 0.47953965147452792817,-0.90109638077934106626,
                                      0.46519399935013672209,-0.88122392525212978054,
                   intersection[1]))
      return EXIT_FAILURE;
    if (!intersect(LON_CIRCLE_EDGE,   0.46633015951723494341,-0.88357293382212931387,
                                      0.46633015951723494341,-0.85902924121595902740,
                   GREAT_CIRCLE_EDGE, 0.50100250912415544846,-0.88440597750634275531,
                                      0.46519399935013672209,-0.88122392525212978054,
                   intersection[2]))
      return EXIT_FAILURE;
    if (!intersect(LON_CIRCLE_EDGE,   0.49087385212340517437,-0.88357293382212931387,
                                      0.49087385212340517437,-0.85902924121595902740,
                   GREAT_CIRCLE_EDGE, 0.50100250912415544846,-0.88440597750634275531,
                                      0.46519399935013672209,-0.88122392525212978054,
                   intersection[3]))
      return EXIT_FAILURE;
    for (int i = 0; i < 4; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 0.46633015951723494341,
                                      0.49087385212340517437,
                                      0.49087385212340517437,
                                      0.46633015951723494341},
                          (double[]){-0.88357293382212931387,
                                     -0.88357293382212931387,
                                     -0.85902924121595902740,
                                     -0.85902924121595902740}, latlon_edges, 4),
        generate_cell_deg((double[]){ 0.47953965147452792817,
                                      0.50100250912415544846,
                                      0.46519399935013672209},
                          (double[]){-0.90109638077934106626,
                                     -0.88440597750634275531,
                                     -0.88122392525212978054}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.49087385212340517437,
                                     intersection_lon[0],
                                     0.46633015951723494341,
                                     0.46633015951723494341,
                                     0.49087385212340517437},
                          (double[]){-0.88357293382212931387,
                                     intersection_lat[0],
                                     intersection_lat[1],
                                     intersection_lat[2],
                                     intersection_lat[3]},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.4) example provided by Uwe \n\n");
#endif

  {

    double intersection[4][3], intersection_lon[4], intersection_lat[4];

    if (!intersect(GREAT_CIRCLE_EDGE, 180.25,-89.75,180.25,-89.50,
                   GREAT_CIRCLE_EDGE, 180.2,-89.6,180.3,-89.6, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 180.25,-89.75,180.25,-89.50,
                   GREAT_CIRCLE_EDGE, 180.3,-89.5,180.2,-89.5, intersection[1]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 180.25,-89.5,180.00,-89.5,
                   GREAT_CIRCLE_EDGE, 180.3,-89.5,180.2,-89.5, intersection[2]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 180.25,-89.5,180.00,-89.5,
                   GREAT_CIRCLE_EDGE, 180.2,-89.5,180.2,-89.6, intersection[3]))
      return EXIT_FAILURE;
    for (int i = 0; i < 4; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){180.2,180.3,180.3,180.2},
                          (double[]){-89.6,-89.6,-89.5,-89.5}, gc_edges, 4),
        generate_cell_deg((double[]){180.00,180.25,180.25,180.00},
                          (double[]){-89.75,-89.75,-89.50,-89.50}, gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){180.2,180.25,180.25,
                                     intersection_lon[2],180.20},
                          (double[]){-89.6,intersection_lat[0],
                                     intersection_lat[1],intersection_lat[2],
                                     intersection_lat[3]}, gc_edges, 5)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.5) example taken from icon toy \n\n");
#endif

  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(LAT_CIRCLE_EDGE,   1.20264093770234/YAC_RAD,
                                -0.76085447079128/YAC_RAD,
                                 1.22718463030851/YAC_RAD,
                                -0.76085447079128/YAC_RAD,
                   GREAT_CIRCLE_EDGE, 1.20552692320133/YAC_RAD,
                                -0.76075226373764/YAC_RAD,
                                 1.30774719967050/YAC_RAD,
                                -0.76075226373764/YAC_RAD, intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(LAT_CIRCLE_EDGE,   1.20264093770234/YAC_RAD,
                                -0.76085447079128/YAC_RAD,
                                 1.22718463030851/YAC_RAD,
                                -0.76085447079128/YAC_RAD,
                   GREAT_CIRCLE_EDGE, 1.20552692320133/YAC_RAD,
                                -0.76075226373764/YAC_RAD,
                                 1.25663706143592/YAC_RAD,
                                -0.82925269516665/YAC_RAD, intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells

      (struct grid_cell[2]){
        generate_cell_deg((double[]){ 1.20264093770234/YAC_RAD,
                                      1.22718463030851/YAC_RAD,
                                      1.22718463030851/YAC_RAD,
                                      1.20264093770234/YAC_RAD},
                          (double[]){-0.76085447079128/YAC_RAD,
                                     -0.76085447079128/YAC_RAD,
                                     -0.73631077818511/YAC_RAD,
                                     -0.73631077818511/YAC_RAD},
                          latlon_edges, 4),
        generate_cell_deg((double[]){ 1.20552692320133/YAC_RAD,
                                      1.25663706143592/YAC_RAD,
                                      1.30774719967050/YAC_RAD},
                          (double[]){-0.76075226373764/YAC_RAD,
                                     -0.82925269516665/YAC_RAD,
                                     -0.76075226373763/YAC_RAD}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){1.20552692320133/YAC_RAD,
                                     intersection_lon[0], intersection_lon[1]},
                          (double[]){-0.76075226373764/YAC_RAD,
                                     intersection_lat[0], intersection_lat[1]},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 3),
        generate_cell_deg((double[]){1.20552692320133/YAC_RAD,
                                     intersection_lon[0],
                                     1.22718463030851/YAC_RAD,
                                     intersection_lon[1]},
                          (double[]){-0.76075226373764/YAC_RAD,
                                     intersection_lat[0],
                                     -0.76085447079128/YAC_RAD,
                                     intersection_lat[1]},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 4)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 14.6) example in which all coordinates are 0 (occurred in a test"
          " of Uwe due to erroneous input data) \n\n");
#endif

  {
    test_clipping(
      // input cells

      (struct grid_cell[2]){
        generate_cell_deg((double[]){0, 0, 0, 0},
                          (double[]){0, 0, 0, 0}, latlon_edges, 4),
        generate_cell_deg((double[]){0, 0, 0},
                          (double[]){0, 0, 0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.7) example example take from ICON_toy \n\n");
#endif
  {
    test_clipping(
      // input cells

      (struct grid_cell[2]){
        generate_cell_deg((double[]){3.1415926372157803/YAC_RAD,
                                     0.9682712007396151/YAC_RAD,
                                     -1.8849555862661402/YAC_RAD},
                          (double[]){-1.5418135425481947/YAC_RAD,
                                     -1.5707963262077691/YAC_RAD,
                                     -1.5418135428532602/YAC_RAD},
                          gc_edges, 3),
        generate_cell_deg((double[]){1.2517283229146832/YAC_RAD,
                                     1.2762720155208536/YAC_RAD,
                                     1.2762720155208536/YAC_RAD,
                                     1.2517283229146832/YAC_RAD},
                          (double[]){-1.5707963267948966/YAC_RAD,
                                     -1.5707963267948966/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.8) example take from ICON_toy \n\n");
#endif
  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 0.9682712007396151/YAC_RAD,
                                -1.5707963262077691/YAC_RAD,
                                 3.1415926372157803/YAC_RAD,
                                -1.5418135425481947/YAC_RAD,
                   LAT_CIRCLE_EDGE,   3.1170489609836229/YAC_RAD,
                                -1.5462526341887264/YAC_RAD,
                                 3.1415926535897931/YAC_RAD,
                                -1.5462526341887264/YAC_RAD,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -1.8849555862661402/YAC_RAD,
                                 -1.5418135428532602/YAC_RAD,
                                  3.1415926372157803/YAC_RAD,
                                 -1.5418135425481947/YAC_RAD,
                   LON_CIRCLE_EDGE,    3.1415926535897931/YAC_RAD,
                                 -1.5462526341887264/YAC_RAD,
                                  3.1415926535897931/YAC_RAD,
                                  -1.5217089415825560/YAC_RAD,
                   intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells

      (struct grid_cell[]){
        generate_cell_deg((double[]){ 3.1415926372157803/YAC_RAD,
                                      0.9682712007396151/YAC_RAD,
                                     -1.8849555862661402/YAC_RAD},
                          (double[]){-1.5418135425481947/YAC_RAD,
                                     -1.5707963262077691/YAC_RAD,
                                     -1.5418135428532602/YAC_RAD}, gc_edges, 3),
        generate_cell_deg((double[]){ 3.1170489609836229/YAC_RAD,
                                      3.1415926535897931/YAC_RAD,
                                      3.1415926535897931/YAC_RAD,
                                      3.1170489609836229/YAC_RAD},
                          (double[]){-1.5462526341887264/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD,
                                     -1.5217089415825560/YAC_RAD,
                                     -1.5217089415825560/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 3.1415926535897931/YAC_RAD,
                                      intersection_lon[0],
                                      3.1415926372157803/YAC_RAD,
                                      intersection_lon[1]},
                          (double[]){-1.5462526341887264/YAC_RAD,
                                      intersection_lat[0],
                                     -1.5418135425481947/YAC_RAD,
                                      intersection_lat[1]},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 4),
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 14.9) example example take from ICON_toy \n\n");
#endif
  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, -1.2564300972590368/YAC_RAD,
                                  1.5418206553283160/YAC_RAD,
                                  1.2277260313871721/YAC_RAD,
                                  1.5707871564758291/YAC_RAD,
                   LON_CIRCLE_EDGE,    0.0000000000000000/YAC_RAD,
                                  1.5462526341887264/YAC_RAD,
                                  0.0000000000000000/YAC_RAD,
                                  1.5707963267948966/YAC_RAD,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -1.2564300972590368/YAC_RAD,
                                  1.5418206553283160/YAC_RAD,
                                  1.2277260313871721/YAC_RAD,
                                  1.5707871564758291/YAC_RAD,
                   LON_CIRCLE_EDGE,    0.0245436926061703/YAC_RAD,
                                  1.5462526341887264/YAC_RAD,
                                  0.0245436926061703/YAC_RAD,
                                  1.5707963267948966/YAC_RAD,
                   intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-1.2564300972590368/YAC_RAD,
                                      1.2277260313871721/YAC_RAD,
                                     -2.5134348157836111/YAC_RAD},
                          (double[]){ 1.5418206553283160/YAC_RAD,
                                      1.5707871564758291/YAC_RAD,
                                      1.5418209335124486/YAC_RAD}, gc_edges, 3),
        generate_cell_deg((double[]){ 0.0000000000000000/YAC_RAD,
                                      0.0245436926061703/YAC_RAD,
                                      0.0245436926061703/YAC_RAD,
                                      0.0000000000000000/YAC_RAD},
                          (double[]){ 1.5462526341887264/YAC_RAD,
                                      1.5462526341887264/YAC_RAD,
                                      1.5707963267948966/YAC_RAD,
                                      1.5707963267948966/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0,intersection_lon[0],intersection_lon[1]},
                          (double[]){90.0,intersection_lat[0],
                                     intersection_lat[1]},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.10) example example take from ICON_toy \n\n");
#endif
  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, -2.1798479726998332/YAC_RAD,
                                 -0.0311624518327736/YAC_RAD,
                                 -2.1991148570983681/YAC_RAD,
                                  0.0000000043633749/YAC_RAD,
                   LAT_CIRCLE_EDGE,    4.0742529726242633/YAC_RAD,
                                  0.0000000000000000/YAC_RAD,
                                  4.0987966652304335/YAC_RAD,
                                  0.0000000000000000/YAC_RAD,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -2.2190472361409284/YAC_RAD,
                                 -0.0338654325037920/YAC_RAD,
                                 -2.1991148570983681/YAC_RAD,
                                  0.0000000043633749/YAC_RAD,
                   LAT_CIRCLE_EDGE,    4.0742529726242633/YAC_RAD,
                                  0.0000000000000000/YAC_RAD,
                                  4.0987966652304335/YAC_RAD,
                                  0.0000000000000000/YAC_RAD,
                   intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-2.1798479726998332/YAC_RAD,
                                     -2.1991148570983681/YAC_RAD,
                                     -2.2190472361409284/YAC_RAD},
                          (double[]){-0.0311624518327736/YAC_RAD,
                                      0.0000000043633749/YAC_RAD,
                                     -0.0338654325037920/YAC_RAD},
                          gc_edges, 3),
        generate_cell_deg((double[]){ 4.0742529726242633/YAC_RAD,
                                      4.0987966652304335/YAC_RAD,
                                      4.0987966652304335/YAC_RAD,
                                      4.0742529726242633/YAC_RAD},
                          (double[]){ 0.0000000000000000/YAC_RAD,
                                      0.0000000000000000/YAC_RAD,
                                      0.0245436926061703/YAC_RAD,
                                      0.0245436926061703/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-2.1991148570983681/YAC_RAD,
                                     intersection_lon[0],
                                     intersection_lon[1]},
                          (double[]){0.0000000043633749/YAC_RAD,
                                     0.0000000000000000/YAC_RAD,
                                     0.0000000000000000/YAC_RAD},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.11) example example take from ICON_toy \n\n");
#endif
  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 3.1415926521557882/YAC_RAD,
                                -0.9805860393425995/YAC_RAD,
                                 3.1415926520689506/YAC_RAD,
                                -1.0172219678978514/YAC_RAD,
                   LAT_CIRCLE_EDGE,   3.1415926535897931/YAC_RAD,
                                -1.0062913968529805/YAC_RAD,
                                 3.1170489609836229/YAC_RAD,
                                -1.0062913968529805/YAC_RAD,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 3.1415926521557882/YAC_RAD,
                                -0.9805860393425995/YAC_RAD,
                                 3.1415926520689506/YAC_RAD,
                                -1.0172219678978514/YAC_RAD,
                   LAT_CIRCLE_EDGE,   3.1415926535897931/YAC_RAD,
                                -0.9817477042468103/YAC_RAD,
                                 3.1170489609836229/YAC_RAD,
                                -0.9817477042468103/YAC_RAD,
                   intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 3.1415926521557882/YAC_RAD,
                                      3.1415926520689506/YAC_RAD,
                                     -3.0774466652769172/YAC_RAD},
                          (double[]){-0.9805860393425995/YAC_RAD,
                                     -1.0172219678978514/YAC_RAD,
                                     -0.9979427097227050/YAC_RAD}, gc_edges, 3),
        generate_cell_deg((double[]){ 3.1170489609836229/YAC_RAD,
                                      3.1415926535897931/YAC_RAD,
                                      3.1415926535897931/YAC_RAD,
                                      3.1170489609836229/YAC_RAD},
                          (double[]){-1.0062913968529805/YAC_RAD,
                                     -1.0062913968529805/YAC_RAD,
                                     -0.9817477042468103/YAC_RAD,
                                     -0.9817477042468103/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 3.1415926535897931/YAC_RAD,
                                      3.1415926535897931/YAC_RAD,
                                     intersection_lon[0],intersection_lon[1]},
                          (double[]){-0.9817477042468103/YAC_RAD,
                                     -1.0062913968529805/YAC_RAD,
                                     -1.0062913968529805/YAC_RAD,
                                     -0.9817477042468103/YAC_RAD},
                          (enum yac_edge_type[]){
                            LON_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE}, 4),
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 14.12) example example take from ICON_toy \n\n");
#endif
  {

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 1.8849555814614660/YAC_RAD,
                                      0.9682712007396151/YAC_RAD,
                                      3.1415926372157803/YAC_RAD},
                          (double[]){-1.5418135420482968/YAC_RAD,
                                     -1.5707963262077691/YAC_RAD,
                                     -1.5418135425481947/YAC_RAD}, gc_edges, 3),
        generate_cell_deg((double[]){ 1.0799224746714915/YAC_RAD,
                                      1.1044661672776617/YAC_RAD,
                                      1.1044661672776617/YAC_RAD,
                                      1.0799224746714915/YAC_RAD},
                          (double[]){-1.5707963267948966/YAC_RAD,
                                     -1.5707963267948966/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.13) example example take from ICON_toy \n\n");
#endif
  {

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 3.1415926372157803/YAC_RAD,
                                      0.9682712007396151/YAC_RAD,
                                     -1.8849555862661402/YAC_RAD},
                          (double[]){-1.5418135425481947/YAC_RAD,
                                     -1.5707963262077691/YAC_RAD,
                                     -1.5418135428532602/YAC_RAD}, gc_edges, 3),
        generate_cell_deg((double[]){ 1.5707963267948966/YAC_RAD,
                                      1.5953400194010670/YAC_RAD,
                                      1.5953400194010670/YAC_RAD,
                                      1.5707963267948966/YAC_RAD},
                          (double[]){-1.5707963267948966/YAC_RAD,
                                     -1.5707963267948966/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 14.14) example example take from ICON_toy \n\n");
#endif
  {
    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 1.25664579522418/YAC_RAD,
                                 1.54180416508220/YAC_RAD,
                                 0.00031736492867/YAC_RAD,
                                 1.54181036404629/YAC_RAD,
                   LAT_CIRCLE_EDGE,   0.90811662642830/YAC_RAD,
                                 1.54625263418873/YAC_RAD,
                                 0.93266031903447/YAC_RAD,
                                 1.54625263418873/YAC_RAD,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 1.25664579522418/YAC_RAD,
                                 1.54180416508220/YAC_RAD,
                                 0.00031736492867/YAC_RAD,
                                 1.54181036404629/YAC_RAD,
                   LON_CIRCLE_EDGE,   0.93266031903447/YAC_RAD,
                                 1.52170894158256/YAC_RAD,
                                 0.93266031903447/YAC_RAD,
                                 1.54625263418873/YAC_RAD,
                   intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){1.25664579522418/YAC_RAD,
                                     1.22772603138717/YAC_RAD,
                                     0.00031736492867/YAC_RAD},
                          (double[]){1.54180416508220/YAC_RAD,
                                     1.57078715647583/YAC_RAD,
                                     1.54181036404629/YAC_RAD}, gc_edges, 3),
        generate_cell_deg((double[]){0.90811662642830/YAC_RAD,
                                     0.93266031903447/YAC_RAD,
                                     0.93266031903447/YAC_RAD,
                                     0.90811662642830/YAC_RAD},
                          (double[]){1.52170894158256/YAC_RAD,
                                     1.52170894158256/YAC_RAD,
                                     1.54625263418873/YAC_RAD,
                                     1.54625263418873/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.93266031903447/YAC_RAD,
                                     intersection_lon[0],
                                     intersection_lon[1]},
                          (double[]){1.54625263418873/YAC_RAD,
                                     intersection_lat[0],
                                     intersection_lat[1]},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 3),
        generate_cell_deg((double[]){0.93266031903447/YAC_RAD,
                                     0.90811662642830/YAC_RAD,
                                     intersection_lon[0],
                                     intersection_lon[1]},
                          (double[]){1.54625263418873/YAC_RAD,
                                     1.54625263418873/YAC_RAD,
                                     intersection_lat[0],
                                     intersection_lat[1]},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 4)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 14.15) example example take from ICON_toy \n\n");
#endif
  {
    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 4.9999997254688875/YAC_RAD,
                                      4.9999997254688875/YAC_RAD,
                                      5.0078539201557488/YAC_RAD,
                                      5.0078539201557488/YAC_RAD},
                          (double[]){-1.2878426515089207/YAC_RAD,
                                     -1.2946756570757916/YAC_RAD,
                                     -1.2946797849754812/YAC_RAD,
                                     -1.2878469125666649/YAC_RAD},
                          gc_edges, 4),
        generate_cell_deg((double[]){ 5.0069132916587327/YAC_RAD,
                                      5.0130492148102759/YAC_RAD,
                                      5.0130492148102759/YAC_RAD,
                                      5.0069132916587327/YAC_RAD},
                          (double[]){-1.2946797849754812/YAC_RAD,
                                     -1.2946797849754812/YAC_RAD,
                                     -1.2885438618239387/YAC_RAD,
                                     -1.2885438618239387/YAC_RAD},
                          latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 5.0069132916587327/YAC_RAD,
                                      5.0078539201557488/YAC_RAD,
                                      5.0078539201557488/YAC_RAD,
                                      5.0069132916587327/YAC_RAD},
                          (double[]){-1.2946797849754812/YAC_RAD,
                                     -1.2946797849754812/YAC_RAD,
                                     -1.2885438618239387/YAC_RAD,
                                     -1.2885438618239387/YAC_RAD},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 15) example provided by Uwe \n\n");
#endif
  {

    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, 2.5132741228718345/YAC_RAD,
                                 1.5494144670074232/YAC_RAD,
                                 3.7699111843077517/YAC_RAD,
                                 1.5494144670074232/YAC_RAD,
                   LAT_CIRCLE_EDGE,   3.2637657012293961/YAC_RAD,
                                 1.5533430342749532/YAC_RAD,
                                 3.2812189937493392/YAC_RAD,
                                 1.5533430342749532/YAC_RAD,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 2.5132741228718345/YAC_RAD,
                                 1.5494144670074232/YAC_RAD,
                                 3.7699111843077517/YAC_RAD,
                                 1.5494144670074232/YAC_RAD,
                   LON_CIRCLE_EDGE,   3.2812189937493392/YAC_RAD,
                                 1.5533430342749532/YAC_RAD,
                                 3.2812189937493392/YAC_RAD,
                                 1.5358897417550099/YAC_RAD,
                   intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){3.2637657012293961/YAC_RAD,
                                     3.2812189937493392/YAC_RAD,
                                     3.2812189937493392/YAC_RAD,
                                     3.2637657012293961/YAC_RAD},
                          (double[]){1.5358897417550099/YAC_RAD,
                                     1.5358897417550099/YAC_RAD,
                                     1.5533430342749532/YAC_RAD,
                                     1.5533430342749532/YAC_RAD},
                          latlon_edges, 4),
        generate_cell_deg((double[]){1.2566370614359175/YAC_RAD,
                                     2.5132741228718345/YAC_RAD,
                                     3.7699111843077517/YAC_RAD,
                                     5.0265482457436690/YAC_RAD,
                                     0.0000000000000000/YAC_RAD,
                                     0.0000000000000000/YAC_RAD},
                          (double[]){1.5494144670074232/YAC_RAD,
                                     1.5494144670074232/YAC_RAD,
                                     1.5494144670074232/YAC_RAD,
                                     1.5494144670074232/YAC_RAD,
                                     1.5494144670074232/YAC_RAD,
                                     1.5494144670074232/YAC_RAD},
                          gc_edges, 6)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){intersection_lon[0],
                                     3.2812189937493392/YAC_RAD,
                                     3.2812189937493392/YAC_RAD},
                          (double[]){1.5533430342749532/YAC_RAD,
                                     1.5533430342749532/YAC_RAD,
                                     intersection_lat[1]},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 16.1) example provided by Rene \n\n");
#endif
  {

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 1.8849563731522099/YAC_RAD,
                                      1.8465416851483734/YAC_RAD,
                                      1.8661276404324745/YAC_RAD},
                          (double[]){-0.4124872765902973/YAC_RAD,
                                     -0.4149531265147232/YAC_RAD,
                                     -0.4401304121445256/YAC_RAD},
                          gc_edges, 3),
        generate_cell_deg((double[]){ 1.8443992887189797/YAC_RAD,
                                      1.8849555921538759/YAC_RAD,
                                      1.9255118955887724/YAC_RAD},
                          (double[]){-0.4121554426487201/YAC_RAD,
                                     -0.4636476090008059/YAC_RAD,
                                     -0.4121554426487201/YAC_RAD},
                          gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 1.8849563731522099/YAC_RAD,
                                      1.8465416851483734/YAC_RAD,
                                      1.8661276404324745/YAC_RAD},
                          (double[]){-0.4124872765902973/YAC_RAD,
                                     -0.4149531265147232/YAC_RAD,
                                     -0.4401304121445256/YAC_RAD}, gc_edges, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 16.2) example provided by Rene \n\n");
#endif
  {
    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, -1.9037843097448526/YAC_RAD,
                                 -0.4401294427665800/YAC_RAD,
                                 -1.8849555925510342/YAC_RAD,
                                 -0.4636476094897370/YAC_RAD,
                   GREAT_CIRCLE_EDGE, -1.8849555921538759/YAC_RAD,
                                 -0.4636476090008059/YAC_RAD,
                                 -1.9528559331227824/YAC_RAD,
                                 -0.4822665528104840/YAC_RAD,
                   intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-1.9037843097448526/YAC_RAD,
                                     -1.8849555925510342/YAC_RAD,
                                     -1.8661268754099334/YAC_RAD},
                          (double[]){-0.4401294427665800/YAC_RAD,
                                     -0.4636476094897370/YAC_RAD,
                                     -0.4401294427950335/YAC_RAD},
                          gc_edges, 3),
        generate_cell_deg((double[]){-1.8849555921538759/YAC_RAD,
                                     -1.8849555921538759/YAC_RAD,
                                     -1.9528559331227824/YAC_RAD},
                          (double[]){-0.5268929705698896/YAC_RAD,
                                     -0.4636476090008059/YAC_RAD,
                                     -0.4822665528104840/YAC_RAD},
                          gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-1.8849555921538759/YAC_RAD,
                                     -1.8849555921538759/YAC_RAD,
                                     intersection_lon},
                          (double[]){-0.4636476094897370/YAC_RAD,
                                     -0.4636476090008059/YAC_RAD,
                                     intersection_lat}, gc_edges, 3),
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 16.3) example provided by Rene \n\n");
#endif
  {

    double intersection[2][3], intersection_lon[2], intersection_lat[2];

    if (!intersect(GREAT_CIRCLE_EDGE, -1.2566370614359172/YAC_RAD,
                                 -0.5535743588970451/YAC_RAD,
                                 -1.1715639552378114/YAC_RAD,
                                 -0.5519553785110131/YAC_RAD,
                   GREAT_CIRCLE_EDGE, -1.2566370573287293/YAC_RAD,
                                 -0.5535743562505872/YAC_RAD,
                                 -1.2346024154180213/YAC_RAD,
                                 -0.5882204740702072/YAC_RAD,
                   intersection[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, -1.2566370614359172/YAC_RAD,
                                 -0.5535743588970451/YAC_RAD,
                                 -1.1715639552378114/YAC_RAD,
                                 -0.5519553785110131/YAC_RAD,
                   GREAT_CIRCLE_EDGE, -1.2566370573287293/YAC_RAD,
                                 -0.5535743562505872/YAC_RAD,
                                 -1.2786717012447690/YAC_RAD,
                                 -0.5882204751704803/YAC_RAD,
                   intersection[1]))
      return EXIT_FAILURE;
    for (int i = 0; i < 2; ++i) {
      XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
      intersection_lon[i] /= YAC_RAD;
      intersection_lat[i] /= YAC_RAD;
    }

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-1.2346024154180213/YAC_RAD,
                                     -1.2566370573287293/YAC_RAD,
                                     -1.2786717012447690/YAC_RAD},
                          (double[]){-0.5882204740702072/YAC_RAD,
                                     -0.5535743562505872/YAC_RAD,
                                     -0.5882204751704803/YAC_RAD},
                          gc_edges, 3),
        generate_cell_deg((double[]){-1.2156135107669752/YAC_RAD,
                                     -1.2566370614359172/YAC_RAD,
                                     -1.1715639552378114/YAC_RAD},
                          (double[]){-0.4835340719879115/YAC_RAD,
                                     -0.5535743588970451/YAC_RAD,
                                     -0.5519553785110131/YAC_RAD},
                          gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0),
        generate_cell_deg((double[]){-1.2566370573287293/YAC_RAD,
                                     intersection_lon[0],
                                     intersection_lon[1]},
                          (double[]){-0.5535743562505872/YAC_RAD,
                                     intersection_lat[0],
                                     intersection_lat[1]}, gc_edges, 3)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 17.1) handling of 'triangle' in which all corners on on a gc \n\n");
#endif
  {
    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0, 1, 2},
                          (double[]){0, 0, 0}, gc_edges, 3),
        generate_cell_deg((double[]){0, 1, 2},
                          (double[]){0, 0, 0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0),
        generate_cell_deg((double[]){0, 1, 2},
                          (double[]){0, 0, 0}, gc_edges, 3)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 17.1) handling of 'triangle' in which all corners on on a gc \n\n");
#endif
  {
    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0, 1, 1},
                          (double[]){0, 0, 1}, gc_edges, 3),
        generate_cell_deg((double[]){0, 1, 2},
                          (double[]){0, 0, 0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0),
        generate_cell_deg((double[]){0, 1, 2},
                          (double[]){0, 0, 0}, gc_edges, 3)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 17.2)\n\n");
#endif
  {

    test_clipping(

      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 1.5462526341887264/YAC_RAD,
                                      1.5707963267948966/YAC_RAD,
                                      1.5707963267948966/YAC_RAD,
                                      1.5462526341887264/YAC_RAD},
                          (double[]){-1.5707963267948966/YAC_RAD,
                                     -1.5707963267948966/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD},
                          latlon_edges, 4),
        generate_cell_deg((double[]){ 1.5707963267948966/YAC_RAD,
                                      1.5707963267948966/YAC_RAD,
                                      1.2491662739534428/YAC_RAD,
                                      1.1072089649419674/YAC_RAD},
                          (double[]){-1.5462526341887277/YAC_RAD,
                                     -1.5339807878856395/YAC_RAD,
                                     -1.5319928455462422/YAC_RAD,
                                     -1.5433578475477794/YAC_RAD},
                          gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 17.3)\n\n");
#endif
  {

    double intersection[3], intersection_lon, intersection_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, 0/YAC_RAD,
                                -1.5339807878856395/YAC_RAD,
                                 0.32163005284145374/YAC_RAD,
                                -1.5319928455462422/YAC_RAD,
                   LON_CIRCLE_EDGE,   0.024543692606170259/YAC_RAD,
                                -1.5462526341887264/YAC_RAD,
                                 0.024543692606170259/YAC_RAD,
                                -1.5217089415825564/YAC_RAD,
                   intersection))
      return EXIT_FAILURE;
    XYZtoLL(intersection, &intersection_lon, &intersection_lat);
    intersection_lon /= YAC_RAD;
    intersection_lat /= YAC_RAD;

    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 0/YAC_RAD,
                                      0.024543692606170259/YAC_RAD,
                                      0.024543692606170259/YAC_RAD,
                                      0/YAC_RAD},
                          (double[]){-1.5462526341887264/YAC_RAD,
                                     -1.5462526341887264/YAC_RAD,
                                     -1.5217089415825564/YAC_RAD,
                                     -1.5217089415825564/YAC_RAD},
                          latlon_edges, 4),
        generate_cell_deg((double[]){ 0/YAC_RAD,
                                      0.32163005284145374/YAC_RAD,
                                      0.24480144176769322/YAC_RAD,
                                      0/YAC_RAD},
                          (double[]){-1.5339807878856395/YAC_RAD,
                                     -1.5319928455462422/YAC_RAD,
                                     -1.5202029837469972/YAC_RAD,
                                     -1.5217089415825564/YAC_RAD},
                          gc_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){ 0/YAC_RAD,
                                      0.024543692606170259/YAC_RAD,
                                      0.024543692606170259/YAC_RAD,
                                      0/YAC_RAD},
                          (double[]){-1.5217089415825564/YAC_RAD,
                                     -1.5217089415825564/YAC_RAD,
                                     intersection_lat,
                                     -1.5339807878856395/YAC_RAD},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE,LON_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,LON_CIRCLE_EDGE}, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 18.1) lon-lat cell that touches the pole \n\n");
#endif
  /* cell a (circle around a pole; made of lat circle edges)

       .-""-.
      /      \
     ;        ;
      \      /
       '-..-'
  */

  /* cell b (lon-lat cell whose upper bound is a zero lenght edge at the pole)

         /\
        /  \
       '-..-'
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-5.0, 5.0, 0.0},
                          (double[]){85.0, 85.0, 90.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE}, 3),
        generate_cell_deg((double[]){  0.0,60.0,120.0,180.0,240.0,300.0},
                          (double[]){ 84.0,84.0,84.0,84.0,84.0,84.0},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                            GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}, 6)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-5.0, 5.0, 0.0},
                          (double[]){85.0, 85.0, 90.0},
                          (enum yac_edge_type[]){
                            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE}, 3)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 19.1) two square that partially share two edges \n\n");
#endif
  /* 11.) source cell test data

    -1      1
     +------+  1
     |      |
     |      |
     |      |
     +------+ -1
  */
  /* target cell test data

     +------+  0.5
     |      |
     |      |
     |      |
     +------+ -0.5
    -1      1
  */

  {
    test_clipping(
      // input cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-1.0, 1.0, 1.0,-1.0},
                          (double[]){-1.0,-1.0, 1.0, 1.0}, latlon_edges, 4),
        generate_cell_deg((double[]){-1.0, 1.0, 1.0,-1.0},
                          (double[]){-0.5,-0.5, 0.5, 0.5}, latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-1.0, 1.0, 1.0,-1.0},
                          (double[]){-0.5,-0.5, 0.5, 0.5}, latlon_edges, 4)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 20.1) two polygons with two edges touch in a single point\n\n");
#endif
  /* 11.) source cell test data

    -10     10
     +------+ 90
     |      |
     |      |
     |      |
     +------+ 80 + dx
  */
  /* target cell test data

        --------  80
        \      /
         \    /
          \  /
           \/     0
      -10   0  10
  */

  {
    double touch_point[3], touch_lon, touch_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, -10.0, 80.0, 10.0, 80.0,
                   LON_CIRCLE_EDGE, 0.0, 70.0, 0.0, 90.0, touch_point))
      return EXIT_FAILURE;
    XYZtoLL(touch_point, &touch_lon, &touch_lat);
    touch_lon /= YAC_RAD;
    touch_lat /= YAC_RAD;

    test_clipping(
      // input cells

      (struct grid_cell[2]){
        generate_cell_deg((double[]){-10.0, 10.0, 10.0, -10.0},
                          (double[]){touch_lat, touch_lat, 90.0, 90.0},
                          latlon_edges, 4),
        generate_cell_deg((double[]){0.0, 10.0, -10.0},
                          (double[]){0.0, 80.0, 80.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 20.2) two polygons with two edges touch in a single point\n\n");
#endif
  /* 11.) source cell test data

    -15     15
     +------+ 80 + dx
     |      |
     |      |
     |      |
     +------+ 60
  */
  /* target cell test data

        --------  80
        \      /
         \    /
          \  /
           \/     70
      -10   0  10
  */

  {
    double touch_point[3], touch_lon, touch_lat;

    if (!intersect(GREAT_CIRCLE_EDGE, -10.0, 80.0, 10.0, 80.0,
                   LON_CIRCLE_EDGE, 0.0, 70.0, 0.0, 90.0, touch_point))
      return EXIT_FAILURE;
    XYZtoLL(touch_point, &touch_lon, &touch_lat);
    touch_lon /= YAC_RAD;
    touch_lat /= YAC_RAD;

    test_clipping(
      // input cells

      (struct grid_cell[2]){
        generate_cell_deg((double[]){-15.0, 15.0, 15.0, -15.0},
                          (double[]){60.0, 60.0, touch_lat, touch_lat},
                          latlon_edges, 4),
        generate_cell_deg((double[]){0.0, 10.0, -10.0},
                          (double[]){70.0, 80.0, 80.0}, gc_edges, 3)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){0.0, 10.0, -10.0},
                          (double[]){70.0, 80.0, 80.0}, gc_edges, 3),
        generate_cell_deg((double[]){0.0, 10.0, 0.0, -10.0},
                          (double[]){70.0, 80.0, touch_lat, 80.0}, gc_edges, 4)},
      // number of reference cells
      2);
  }

#ifdef VERBOSE
  printf ("\n 21.1) one cell is empty\n\n");
#endif
  /* 11.) source cell test data

     +------+
     |      |
     |      |
     |      |
     +------+
  */

  {
    test_clipping(
      // input cells

      (struct grid_cell[2]){
        generate_cell_deg((double[]){0.0, 1.0, 1.0, 0.0},
                          (double[]){0.0, 0.0, 1.0, 1.0},
                          latlon_edges, 4),
        generate_cell_deg((double[]){0.0}, (double[]){0.0}, gc_edges, 1)},
      // reference cells
      (struct grid_cell[]){generate_cell_deg(NULL, NULL, NULL, 0)},
      // number of reference cells
      1);
  }

#ifdef VERBOSE
  printf ("\n 22.1) star-shaped cell\n\n");
#endif
  /* source cell test data


             +
            / \
           /   \
          +     +
       -           -
    +                 +
       -           -
          +     +
           \   /
            \ /
             +
  */

  {
    test_clipping(
      // input cells

      (struct grid_cell[2]){
        generate_cell_deg((double[]){0.0,0.5,2.0,0.5,0.0,-0.5,-2.0,-0.5},
                          (double[]){2.0,0.5,0.0,-0.5,-2.0,-0.5,0.0,0.5},
                          gc_edges, 8),
        generate_cell_deg((double[]){-0.4,0.4,0.4,-0.4},
                          (double[]){0.4,0.4,-0.4,-0.4}, latlon_edges, 4)},
      // reference cells
      (struct grid_cell[]){
        generate_cell_deg((double[]){-0.4,0.4,0.4,-0.4},
                          (double[]){0.4,0.4,-0.4,-0.4}, latlon_edges, 4)},
      // number of reference cells
      1);
  }

//-----------------test yac_correct_weights------------------------------------

  {
    enum {num_weights = 5};
    double weights[num_weights] = {1,3,2,0.5,3.5};
    double ref_weights[num_weights] = {0.1,0.3,0.2,0.05,0.35};
    check_weight_correction(weights, num_weights, ref_weights);
  }

  {
    enum {num_weights = 1};
    double weights[num_weights] = {3};
    double ref_weights[num_weights] = {1};
    check_weight_correction(weights, num_weights, ref_weights);
  }

  {
    enum {num_weights = 4};
    double weights[num_weights] = {1,2,3,4};
    double ref_weights[num_weights] = {0.1,0.2,0.3,0.4};
    check_weight_correction(weights, num_weights, ref_weights);
  }

  yac_free_grid_cell(overlap_cell);
  yac_free_grid_cell(overlap_cell+1);

  return TEST_EXIT_CODE;
}

static int compare_cells(struct grid_cell a, struct grid_cell b) {

  double const tol = 1e-8;

  if ((yac_huiliers_area(a) < yac_angle_tol * yac_angle_tol) &&
      (yac_huiliers_area(b) < yac_angle_tol * yac_angle_tol)) return 0;

  if (a.num_corners != b.num_corners) return 1;

  if (a.num_corners == 0) return 0;

  for (int order = -1; order <= 1; order += 2) {
    for (int start = 0; start < (int)(a.num_corners); ++start) {

      int differences = 0;

      for (int i = 0; i < (int)(a.num_corners); ++i) {

        int j =
          ((int)(a.num_corners) + start + order * i) % (int)(a.num_corners);

        if (get_vector_angle(a.coordinates_xyz[i], b.coordinates_xyz[j]) > tol)
          ++differences;

        else if ((a.edge_type[i] == LAT_CIRCLE_EDGE) !=
                 (b.edge_type[
                    (((int)(a.num_corners))+j-(order<0))%
                    ((int)(a.num_corners))] == LAT_CIRCLE_EDGE))
          ++differences;
      }

      if (!differences) return 0;
    }
  }

  return 1;
}

static void test_clipping(struct grid_cell cells[2],
                          struct grid_cell * ref_cells,
                          unsigned num_ref_cells) {

  int order[2], start[2];
  double mem_dummy_[2][128][3];
  enum yac_edge_type edge_dummy[2][128];
  struct grid_cell overlap_cell,
                   test_cells[2] = {
        (struct grid_cell){.coordinates_xyz = mem_dummy_[0],
                           .edge_type = edge_dummy[0],
                           .num_corners = cells[0].num_corners},
        (struct grid_cell){.coordinates_xyz = mem_dummy_[1],
                           .edge_type = edge_dummy[1],
                           .num_corners = cells[1].num_corners}};

  yac_init_grid_cell(&overlap_cell);

  for (order[0] = -1; order[0] <= 1; order[0] += 2) {
    for (order[1] = -1; order[1] <= 1; order[1] += 2) {
      for (start[0] = 0; start[0] < (int)(cells[0].num_corners); ++start[0]) {
        for (start[1] = 0; start[1] < (int)(cells[1].num_corners); ++start[1]) {

          for (int k = 0; k <= 1; ++k) {
            for (int i = 0; i < ((int)(cells[k].num_corners)); ++i) {
              int j =
                (((int)(cells[k].num_corners))+start[k]+i*order[k])%
                ((int)(cells[k].num_corners));
              test_cells[k].coordinates_xyz[i][0] = cells[k].coordinates_xyz[j][0];
              test_cells[k].coordinates_xyz[i][1] = cells[k].coordinates_xyz[j][1];
              test_cells[k].coordinates_xyz[i][2] = cells[k].coordinates_xyz[j][2];
              j =
                (((int)(cells[k].num_corners)) + j - (order[k] < 0))%
                ((int)(cells[k].num_corners));
              test_cells[k].edge_type[i] = cells[k].edge_type[j];
            }
          }

          for (int k = 0; k <= 1; ++k) {
            yac_cell_clipping(1, test_cells + k, test_cells[k^1], &overlap_cell);
            int match = 0;
            for (int i = 0; i < (int)num_ref_cells; ++i)
              match |= !compare_cells(overlap_cell, ref_cells[i]);
            if (!match)
              PUT_ERR("ERROR: wrong clipping cell\n");
          }
        }
      }
    }
  }
  yac_free_grid_cell(&cells[0]);
  yac_free_grid_cell(&cells[1]);
  for (unsigned i = 0; i < num_ref_cells; ++i)
    yac_free_grid_cell(&ref_cells[i]);
  yac_free_grid_cell(&overlap_cell);
}

static void check_weight_correction(
  double * weights, size_t count, double * ref_weights) {

  double scales[] = {0.001, 0.01, 0.1, 1.0, 10.0, 100.0};
  size_t num_tests = sizeof(scales) / sizeof(scales[0]);

  for (size_t i = 0; i < num_tests; ++i) {

    double temp_weights[count];
    for (size_t j = 0; j < count; ++j)
      temp_weights[j] = weights[j] * scales[i];

    yac_correct_weights(count, temp_weights);

    double weight_diff = 1.0;
    for (size_t j = 0; j < count; ++j) {
      weight_diff -= temp_weights[j];
      if (fabs(ref_weights[j] - temp_weights[j]) > 1e-9)
        PUT_ERR("wrong weight");
    }
    if (fabs(weight_diff) > 1e-12) PUT_ERR("wrong weight sum");
  }
}
