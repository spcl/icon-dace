/**
 * @file test_partial_areas.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include "utils.h"
#include "area.h"
#include "tests.h"
#include "geometry.h"
#include "test_common.h"

static void check_partial_areas_deg(
  double * tgt_lon, double * tgt_lat,
  enum yac_edge_type * tgt_edges, int tgt_count,
  double ** src_lons, double ** src_lats,
  enum yac_edge_type ** src_edges, int * src_counts, size_t nSourceCells);
static void check_partial_areas_rad(
  double * tgt_lon, double * tgt_lat,
  enum yac_edge_type * tgt_edges, int tgt_count,
  double ** src_lons, double ** src_lats,
  enum yac_edge_type ** src_edges, int * src_counts, size_t nSourceCells);

int main (void) {

  enum yac_edge_type gc_edges[] = {GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                                   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                                   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                                   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE};
  enum yac_edge_type latlon_edges[] = {LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
                                       LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE};

  /* 4 source elements overlapping a target cell

    Target Cell

              0.0 (lon)  [1]
              0.7 (lat)
              / \
             /   \
   -0.5 [2] /     \  0.5 [0]
    0.0     \     /  0.0
             \   /
              \ /
              0.0  [3]
             -0.5
  */
  /* source cell test data

         0.0          1.0 (lon)
         1.8          1.8 (lat)
   +-------+--------+
   |       |        |
   |   1   |    0   |
   |       |        |
   |       |        |  1.0
   +-------+--------+  0.6
   |     0.0        |
   |     0.6        |
   |                |
   |   2   |    3   |
   |       |        |
   +-------+--------+  1.0
                      -0.6
  */
  {
    enum {nSourceCells = 4};
    check_partial_areas_deg(
      (double[]){0.5, 0.0, -0.5,  0.0},
      (double[]){0.0, 0.7,  0.0, -0.5}, gc_edges, 4,
      (double*[nSourceCells]){(double[]){1.0, 0.0, 0.0, 1.0},
                              (double[]){0.0, -1.0, -1.0, 0.0},
                              (double[]){0.0, -1.0, -1.0,  0.0},
                              (double[]){1.0, 0.0,  0.0,  1.0}},
      (double*[nSourceCells]){(double[]){1.8, 1.8, 0.6, 0.6},
                              (double[]){1.8,  1.8,  0.6, 0.6},
                              (double[]){0.6,  0.6, -0.6, -0.6},
                              (double[]){0.6, 0.6, -0.6, -0.6}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){4, 4, 4, 4}, nSourceCells);
  }

  /* Next test is over the pole

     Target cell located over the pole o with
     vertices x and source cells 0 to 3 arranged as


 135.0         90.0       45.0 (lon)
               89.0       88.5 (lat)
         ------------------
        |        |        |
        |        |        |
        |       1x        |
        |        |    0   |
-180.0   ---x----o----x----  0.0        pole @ -180.0; 90.0
        |   2    |        |
        |        |        |
        |       3x        |
        |        |        |
        |        |        |
-135.0   ------------------ -45.0
               -90.0         -0.6

  */
  {
    enum {nSourceCells = 4};
    check_partial_areas_deg(
      (double[]){ 0.0, 90.0, -180.0, -90.0},
      (double[]){89.5, 89.5,   89.5,  89.5}, gc_edges, 4,
      (double*[nSourceCells]){(double[]){45.0, 90.0, -180.0,  0.0},
                              (double[]){90.0, 135.0, -180.0, -180.0},
                              (double[]){-180.0, -180.0, -135.0, -90.0},
                              (double[]){ 0.0, -180.0, -90.0, -45.0}},
      (double*[nSourceCells]){(double[]){88.5, 89.0, 90.0, 89.0},
                              (double[]){89.0, 88.5, 89.0, 90.0},
                              (double[]){90.0, 89.0, 88.5, 89.0},
                              (double[]){89.0, 90.0, 89.0, 88.5}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){4, 4, 4, 4}, nSourceCells);
  }

  // checking overlap between identical cells
  {
    enum {nSourceCells = 1};
    check_partial_areas_deg(
      (double[]){-180.0, -150.0, -150.0, -180.0},
      (double[]){ -90.0,  -90.0,  -60.0,  -60.0}, gc_edges, 4,
      (double*[nSourceCells]){(double[]){-180.0, -150.0, -150.0, -180.0}},
      (double*[nSourceCells]){(double[]){ -90.0,  -90.0,  -60.0,  -60.0}},
      (enum yac_edge_type*[nSourceCells]){gc_edges},
      (int[nSourceCells]){4}, nSourceCells);
  }
  {
    enum {nSourceCells = 1};
    check_partial_areas_deg(
      (double[]){-180.0, -150.0, -150.0, -180.0},
      (double[]){  60.0,   60.0,   90.0,   90.0}, gc_edges, 4,
      (double*[nSourceCells]){(double[]){-180.0, -150.0, -150.0, -180.0}},
      (double*[nSourceCells]){(double[]){  60.0,   60.0,   90.0,   90.0}},
      (enum yac_edge_type*[nSourceCells]){gc_edges},
      (int[nSourceCells]){4}, nSourceCells);
  }

  // one out of three source cells completely overlaps with the target cell
  {
    enum {nSourceCells = 3};
    check_partial_areas_deg(
      (double[]){25.0, 30.0, 30.0, 25.0},
      (double[]){75.0, 75.0, 80.0, 80.0}, gc_edges, 4,
      (double*[nSourceCells]){(double[]){30.0, 35.0, 35.0, 30.0},
                              (double[]){25.0, 30.0, 30.0, 25.0},
                              (double[]){20.0, 25.0, 25.0, 20.0}},
      (double*[nSourceCells]){(double[]){75.0, 75.0, 80.0, 80.0},
                              (double[]){75.0, 75.0, 80.0, 80.0},
                              (double[]){75.0, 75.0, 80.0, 80.0}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){4, 4, 4}, nSourceCells);
  }

  // two cell are very close to each other
  {
    struct grid_cell TargetCell =
      generate_cell_rad((double[]){0.098707396293595859,
                                   0.117826295355522680,
                                   0.098174496441371703},
                        (double[]){0.333033372840672190,
                                   0.345484775276814930,
                                   0.353335609258583650}, gc_edges, 3);
    struct grid_cell SourceCell =
      generate_cell_rad((double[]){0.073631077818510776,
                                   0.098174770424681035,
                                   0.098174770424681035,
                                   0.073631077818510776},
                        (double[]){0.343611696486383620,
                                   0.343611696486383620,
                                   0.368155389092553910,
                                   0.368155389092553910}, gc_edges, 4);

    double overlap_area;
    yac_compute_overlap_areas(
      1, &SourceCell, TargetCell, &overlap_area);

    if (fabs(overlap_area) > 1e-9)
      PUT_ERR("ERROR in yac_compute_overlap_areas");

    yac_free_grid_cell(&SourceCell);
    yac_free_grid_cell(&TargetCell);
  }

  // overlapping lon-lat cells
  {
    enum {nSourceCells = 4};
    check_partial_areas_deg(
      (double[]){-0.5,  0.5, 0.5, -0.5},
      (double[]){-0.5, -0.5, 0.5,  0.5}, latlon_edges, 4,
      (double*[nSourceCells]){(double[]){-1,  0, 0, -1},
                              (double[]){ 0,  1, 1, 0},
                              (double[]){-1, 0, 0, -1},
                              (double[]){0, 1, 1, 0}},
      (double*[nSourceCells]){(double[]){-1, -1, 0,  0},
                              (double[]){-1, -1, 0,  0},
                              (double[]){ 0, 0, 1,  1},
                              (double[]){ 0, 0, 1,  1}},
      (enum yac_edge_type*[nSourceCells])
        {latlon_edges, latlon_edges, latlon_edges, latlon_edges},
      (int[nSourceCells]){4, 4, 4, 4}, nSourceCells);
  }

  // overlapping lon-lat cells
  {
    enum {nSourceCells = 1};
    check_partial_areas_deg(
      (double[]){-0.5,  0.5, 0.5, -0.5},
      (double[]){-0.5, -0.5, 0.5,  0.5}, latlon_edges, 4,
      (double*[nSourceCells]){(double[]){-1,  1, 1, -1}},
      (double*[nSourceCells]){(double[]){-1, -1, 1,  1}},
      (enum yac_edge_type*[nSourceCells])
        {latlon_edges},
      (int[nSourceCells]){4}, nSourceCells);
  }

  // lon-lat cells overlapping with gc cell
  {
    enum {nSourceCells = 1};
    check_partial_areas_deg(
      (double[]){-0.5,  0.5, 0.0},
      (double[]){-0.5, -0.5, 0.5}, gc_edges, 3,
      (double*[nSourceCells]){(double[]){-1,  1, 1, -1}},
      (double*[nSourceCells]){(double[]){-1, -1, 1,  1}},
      (enum yac_edge_type*[nSourceCells])
        {latlon_edges},
      (int[nSourceCells]){4}, nSourceCells);
  }

  // test case provided by a user
  {
    enum {nSourceCells = 4};
    check_partial_areas_rad(
      (double[]){2.5034566458293663,
                 2.5280003384355365,
                 2.5280003384355365,
                 2.5034566458293663},
      (double[]){0.5890486225480862,
                 0.5890486225480862,
                 0.6135923151542565,
                 0.6135923151542565}, latlon_edges, 4,
      (double*[nSourceCells]){(double[]){2.5132750564642796,
                                         2.5132751216401110,
                                         2.5530451220055963},
                              (double[]){2.5530451220055963,
                                         2.5132751216401110,
                                         2.5514531161521616},
                              (double[]){2.4735050741253155,
                                         2.5132751216401110,
                                         2.5132750564642796},
                              (double[]){2.4750971948235070,
                                         2.5132751216401110,
                                         2.4735050741253155}},
      (double*[nSourceCells]){(double[]){0.6227096831849960,
                                         0.5888365165503746,
                                         0.6044728875052012},
                              (double[]){0.6044728875052012,
                                         0.5888365165503746,
                                         0.5706926293017953},
                              (double[]){0.6044725603601205,
                                         0.5888365165503746,
                                         0.6227096831849960},
                              (double[]){0.5706923461529829,
                                         0.5888365165503746,
                                         0.6044725603601205}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){3, 3, 3, 3}, nSourceCells);
  }

  // test case provided by a user
  {
    enum {nSourceCells = 3};
    check_partial_areas_rad(
      (double[]){0.6626797003665970,
                 0.6872233929727672,
                 0.6872233929727672,
                 0.6626797003665970},
      (double[]){1.2271846303085130,
                 1.2271846303085130,
                 1.2517283229146832,
                 1.2517283229146832}, latlon_edges, 4,
      (double*[nSourceCells]){(double[]){0.6910166804158336,
                                         0.6283648515564006,
                                         0.7428780447168533},
                              (double[]){0.5657152198049494,
                                         0.6283648515564006,
                                         0.6910166804158336},
                              (double[]){0.6283674880262698,
                                         0.5657152198049494,
                                         0.6910166804158336}},
      (double*[nSourceCells]){(double[]){1.2515918327919406,
                                         1.2207207913058662,
                                         1.2190404189581385},
                              (double[]){1.2515932566235379,
                                         1.2207207913058662,
                                         1.2515918327919406},
                              (double[]){1.2831450892755811,
                                         1.2515932566235379,
                                         1.2515918327919406}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){3, 3, 3}, nSourceCells);
  }

  // test case provided by a user
  {
    enum {nSourceCells = 7};
    check_partial_areas_rad(
      (double[]){-2.1798479726998332,
                 -2.1991148570983681,
                 -2.2190472361409284},
      (double[]){-0.0311624518327736,
                  0.0000000043633749,
                 -0.0338654325037920}, gc_edges, 3,
      (double*[nSourceCells]){(double[]){ 4.0987966652304335,
                                          4.1233403578366037,
                                          4.1233403578366037,
                                          4.0987966652304335},
                              (double[]){ 4.0742529726242633,
                                          4.0987966652304335,
                                          4.0987966652304335,
                                          4.0742529726242633},
                              (double[]){ 4.0742529726242633,
                                          4.0987966652304335,
                                          4.0987966652304335,
                                          4.0742529726242633},
                              (double[]){ 4.0987966652304335,
                                          4.1233403578366037,
                                          4.1233403578366037,
                                          4.0987966652304335},
                              (double[]){ 4.0742529726242633,
                                          4.0987966652304335,
                                          4.0987966652304335,
                                          4.0742529726242633},
                              (double[]){ 4.0497092800180932,
                                          4.0742529726242633,
                                          4.0742529726242633,
                                          4.0497092800180932},
                              (double[]){ 4.0497092800180932,
                                          4.0742529726242633,
                                          4.0742529726242633,
                                          4.0497092800180932}},
      (double*[nSourceCells]){(double[]){-0.0490873852123405,
                                         -0.0490873852123405,
                                         -0.0245436926061703,
                                         -0.0245436926061703},
                              (double[]){-0.0490873852123405,
                                         -0.0490873852123405,
                                         -0.0245436926061703,
                                         -0.0245436926061703},
                              (double[]){ 0.0000000000000000,
                                          0.0000000000000000,
                                          0.0245436926061703,
                                          0.0245436926061703},
                              (double[]){-0.0245436926061703,
                                         -0.0245436926061703,
                                          0.0000000000000000,
                                          0.0000000000000000},
                              (double[]){-0.0245436926061703,
                                         -0.0245436926061703,
                                          0.0000000000000000,
                                          0.0000000000000000},
                              (double[]){-0.0245436926061703,
                                         -0.0245436926061703,
                                          0.0000000000000000,
                                          0.0000000000000000},
                              (double[]){-0.0490873852123405,
                                         -0.0490873852123405,
                                         -0.0245436926061703,
                                         -0.0245436926061703}},
      (enum yac_edge_type*[nSourceCells])
        {latlon_edges, latlon_edges, latlon_edges, latlon_edges,
         latlon_edges, latlon_edges, latlon_edges},
      (int[nSourceCells]){4, 4, 4, 4, 4, 4, 4}, nSourceCells);
  }

  // test case provided by a user
  {
    enum {nSourceCells = 4};
    check_partial_areas_rad(
      (double[]){ 3.4361169648638361,
                  3.4606606574700067,
                  3.4606606574700067,
                  3.4361169648638361},
      (double[]){-0.0245436926061703,
                 -0.0245436926061703,
                  0.0000000000000000,
                  0.0000000000000000}, latlon_edges, 4,
      (double*[nSourceCells]){(double[]){-2.8075010104343816,
                                         -2.8274333887952685,
                                         -2.8467002736213392},
                              (double[]){-2.8081665040027728,
                                         -2.8274333887952685,
                                         -2.7882305020025413},
                              (double[]){-2.7882305020025413,
                                         -2.8274333887952685,
                                         -2.8075010104343816},
                              (double[]){-2.8467002736213392,
                                         -2.8274333887952685,
                                         -2.8666362753717674}},
      (double*[nSourceCells]){(double[]){-0.0338654341705874,
                                          0.0000000011622542,
                                         -0.0311624540855508},
                              (double[]){ 0.0311624566805337,
                                          0.0000000011622542,
                                         -0.0026744729717723},
                              (double[]){-0.0026744729717723,
                                          0.0000000011622542,
                                         -0.0338654341705874},
                              (double[]){-0.0311624540855508,
                                          0.0000000011622542,
                                          0.0026744751289656}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){3, 3, 3, 3}, nSourceCells);
  }

  // test case provided by a user
  // (the overlap between the first source cell and the target cell
  //  is very small, which is why it is potentially dropped by YAC)
  {
    enum {nSourceCells = 3};
    check_partial_areas_rad(
      (double[]){ 3.1170489609836229,
                  3.1415926535897931,
                  3.1415926535897931,
                  3.1170489609836229},
      (double[]){-1.0062913968529805,
                 -1.0062913968529805,
                 -0.9817477042468103,
                 -0.9817477042468103}, latlon_edges, 4,
      (double*[nSourceCells]){(double[]){ 3.1415926521557882,
                                          3.1415926520689506,
                                         -3.0774466652769172},
                              (double[]){ 3.0774466624490269,
                                          3.1415926520689506,
                                          3.1415926521557882},
                              (double[]){ 3.0808931330453300,
                                          3.0774466624490269,
                                          3.1415926521557882}},
      (double*[nSourceCells]){(double[]){-0.9805860393425995,
                                         -1.0172219678978514,
                                         -0.9979427097227050},
                              (double[]){-0.9979427096430752,
                                         -1.0172219678978514,
                                         -0.9805860393425995},
                              (double[]){-0.9613498550843829,
                                         -0.9979427096430752,
                                         -0.9805860393425995}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){3, 3, 3}, nSourceCells);
  }

  // overlap of a triangle with a lot of very narrow lon lat cells
  {
    enum {nSourceCells = 198};
    double ** src_lons = xmalloc(2 * nSourceCells * sizeof(*src_lons));
    double ** src_lats = src_lons + nSourceCells;
    enum yac_edge_type ** src_edges =
      xmalloc(nSourceCells * sizeof(*src_edges));
    int * src_counts = xmalloc(nSourceCells * sizeof(*src_counts));

    src_lons[0] = xmalloc(2 * 4 * nSourceCells * sizeof(**src_lons));
    double dx = ((49.2 - 22.8) /((double)nSourceCells));
    for (size_t i = 0; i < nSourceCells; ++i) {
      src_lons[i] = src_lons[0] + i * 4;
      src_lats[i] = src_lons[0] + i * 4 + 4 * nSourceCells;
      src_lons[i][0] = 22.8000 + dx*((double)(i+0));
      src_lons[i][1] = 22.8000 + dx*((double)(i+1));
      src_lons[i][2] = 22.8000 + dx*((double)(i+1));
      src_lons[i][3] = 22.8000 + dx*((double)(i+0));
      src_lats[i][0] = 89.7333;
      src_lats[i][1] = 89.7333;
      src_lats[i][2] = 89.8667;
      src_lats[i][3] = 89.8667;
      src_edges[i] = latlon_edges;
      src_counts[i] = 4;
    }

     check_partial_areas_deg(
      (double[]){22.8313, 49.1687, 36.0000},
      (double[]){89.7418, 89.7418, 89.8335}, gc_edges, 3,
      src_lons, src_lats, src_edges, src_counts, nSourceCells);

    free(src_counts);
    free(src_edges);
    free(src_lons[0]);
    free(src_lons);
  }

  // test case provided by a user
  {
    enum {nSourceCells = 4};
    check_partial_areas_rad(
      (double[]){ 5.0069132916587327,
                  5.0130492148102759,
                  5.0130492148102759,
                  5.0069132916587327},
      (double[]){-1.2946797849754812,
                 -1.2946797849754812,
                 -1.2885438618239387,
                 -1.2885438618239387}, latlon_edges, 4,
      (double*[nSourceCells]){(double[]){ 4.9999997254688875,
                                          4.9999997254688875,
                                          5.0078539201557488,
                                          5.0078539201557488},
                              (double[]){ 5.0078539201557488,
                                          5.0078539201557488,
                                          5.0157075822103927,
                                          5.0157075822103927},
                              (double[]){ 5.0078539201557488,
                                          5.0078539201557488,
                                          5.0157075822103927,
                                          5.0157075822103927},
                              (double[]){ 4.9999997254688875,
                                          4.9999997254688875,
                                          5.0078539201557488,
                                          5.0078539201557488}},
      (double*[nSourceCells]){(double[]){-1.2878426515089207,
                                         -1.2946756570757916,
                                         -1.2946797849754812,
                                         -1.2878469125666649},
                              (double[]){-1.2946797849754812,
                                         -1.3015127905423520,
                                         -1.3015170516000960,
                                         -1.2946840460332254},
                              (double[]){-1.2878469125666649,
                                         -1.2946797849754812,
                                         -1.2946840460332254,
                                         -1.2878513067824635},
                              (double[]){-1.2946756570757916,
                                         -1.3015085294846078,
                                         -1.3015127905423520,
                                         -1.2946797849754812}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){4, 4, 4, 4}, nSourceCells);
  }

  // test case provided by a user
  {
    enum {nSourceCells = 2};
    check_partial_areas_rad(
      (double[]){-2.8460897445189417,
                 -2.8659852121839959,
                 -2.8853055611946927},
      (double[]){-0.0961602003314924,
                 -0.0622820530745658,
                 -0.0933150814318732}, gc_edges, 3,
      (double*[nSourceCells]){(double[]){-2.8655528147976832,
                                         -2.9038195735214831,
                                         -2.8254194450162888},
                              (double[]){-2.9440830021466091,
                                         -2.9038195735214831,
                                         -2.8655528147976832}},
      (double*[nSourceCells]){(double[]){-0.0615856127966505,
                                         -0.1228534567110898,
                                         -0.1292806814915037},
                              (double[]){-0.0558979744839555,
                                         -0.1228534567110898,
                                         -0.0615856127966505}},
      (enum yac_edge_type*[nSourceCells]){gc_edges, gc_edges},
      (int[nSourceCells]){3, 3}, nSourceCells);
  }

  // test case provided by a user
  {
    enum {nSourceCells = 6};
    check_partial_areas_rad(
      (double[]){3.0127596408133406,
                 2.9702182441873108,
                 2.9943996418061438},
      (double[]){0.5498596133493898,
                 0.5469967117356207,
                 0.5140000458840724}, gc_edges, 3,
      (double*[nSourceCells]){(double[]){2.9198084624304621,
                                         2.9720969737232750,
                                         3.0078074923914948},
                              (double[]){2.8889398582596519,
                                         2.9720969737232750,
                                         2.9198084624304621},
                              (double[]){3.0187962726049782,
                                         2.9720969737232750,
                                         2.9378118571213374},
                              (double[]){2.9378118571213374,
                                         2.9720969737232750,
                                         2.8889398582596519},
                              (double[]){3.0565195473916882,
                                         2.9720969737232750,
                                         3.0187962726049782},
                              (double[]){3.0078074923914948,
                                         2.9720969737232750,
                                         3.0565195473916882}},
      (double*[nSourceCells]){(double[]){0.6123031304013200,
                                         0.5471403770103002,
                                         0.6192726969308253},
                              (double[]){0.5392518296263771,
                                         0.5471403770103002,
                                         0.6123031304013200},
                              (double[]){0.4805223860987507,
                                         0.5471403770103002,
                                         0.4745836829538209},
                              (double[]){0.4745836829538209,
                                         0.5471403770103002,
                                         0.5392518296263771},
                              (double[]){0.5519553785110132,
                                         0.5471403770103002,
                                         0.4805223860987507},
                              (double[]){0.6192726969308253,
                                         0.5471403770103002,
                                         0.5519553785110132}},
      (enum yac_edge_type*[nSourceCells])
        {gc_edges, gc_edges, gc_edges, gc_edges, gc_edges, gc_edges},
      (int[nSourceCells]){3, 3, 3, 3, 3, 3}, nSourceCells);
  }

  { // test the barycenter returned by yac_compute_overlap_info for a
    // target cell with more than 4 edges
    struct grid_cell TargetCell =
      generate_cell_deg((double[]){-0.5,  0.5, 0.75, 0.5, -0.5, -0.75},
                        (double[]){-0.5, -0.5, 0.0 , 0.5,  0.5,  0.0},
                        gc_edges, 6);
    struct grid_cell SourceCells[] = {
      generate_cell_deg((double[]){0.3, -0.3,  0.0},
                        (double[]){0.3,  0.3, -0.3}, gc_edges, 3),
      generate_cell_deg((double[]){-0.4,  0.4, 0.4, -0.4},
                        (double[]){-0.4, -0.4, 0.4,  0.4}, gc_edges, 4),
      generate_cell_deg((double[]){-0.6,  0.6, 0.8, 0.6, -0.6, -0.8},
                        (double[]){-0.6, -0.6, 0.0, 0.6,  0.6,  0.0},
                        gc_edges, 6)};
    enum {nSourceCells = sizeof(SourceCells)/sizeof(SourceCells[0])};

    double ref_area[nSourceCells];
    double ref_barycenter[nSourceCells][3];
    for (int i = 0; i < nSourceCells; ++i)
      for (int j = 0; j < 3; ++j)
        ref_barycenter[i][j] = 0.0;
    {
      ref_area[0] =
        yac_huiliers_area_info(SourceCells[0], ref_barycenter[0], 1.0);
      ref_area[1] = yac_huiliers_area(SourceCells[1]);
      ref_area[2] = yac_huiliers_area(TargetCell);
      LLtoXYZ(0, 0, ref_barycenter[1]);
      LLtoXYZ(0, 0, ref_barycenter[2]);
    }

    double overlap_areas[nSourceCells];
    double overlap_barycenters[nSourceCells][3];
    yac_compute_overlap_info(
      nSourceCells, SourceCells, TargetCell,
      overlap_areas, &overlap_barycenters[0]);

    for (int i = 0; i < nSourceCells; ++i) {
      double area_tolerance = yac_huiliers_area(SourceCells[i]) * 1e-6;
      if (fabs(ref_area[i] - overlap_areas[i]) > area_tolerance)
        PUT_ERR("ERROR in yac_compute_overlap_info (area)");
      if (overlap_areas[i] > area_tolerance) {
        normalise_vector(ref_barycenter[i]);
        normalise_vector(overlap_barycenters[i]);
        if(get_vector_angle(
             ref_barycenter[i], overlap_barycenters[i]) > 1e-9)
          PUT_ERR("ERROR in yac_compute_overlap_info (barycenter)");
      }
    }

    for (size_t i = 0; i < nSourceCells; ++i)
      yac_free_grid_cell(&SourceCells[i]);
    yac_free_grid_cell(&TargetCell);
  }

  { // test the area and barycenter returned by
    // yac_compute_overlap_info for a triangle target cell
    
    struct grid_cell TargetCells[] =
       // normal triangle
      {generate_cell_deg((double[]){-0.5, 0.5, -0.5},
                         (double[]){-0.5, 0.5,  0.5}, gc_edges, 3),
       // triangle with a duplicated vertex to simulate a degenerated cell
       generate_cell_deg((double[]){-0.5, 0.5, -0.5, -0.5},
                         (double[]){-0.5, 0.5,  0.5,  0.5}, gc_edges, 4)};
    enum {nTargetCells = sizeof(TargetCells)/sizeof(TargetCells[0])};

    struct grid_cell SourceCells[] = {
      generate_cell_deg((double[]){-0.6, 0.7, -0.6},
                        (double[]){-0.7, 0.6,  0.6}, gc_edges, 3),
      generate_cell_deg((double[]){-0.4, 0.3, -0.4},
                        (double[]){-0.3, 0.4,  0.4}, gc_edges, 3),
      generate_cell_deg((double[]){-0.5, 0.5, -0.5},
                        (double[]){-0.5, 0.5,  0.5}, gc_edges, 3),
      generate_cell_deg((double[]){0.0, 1.0, 1.0, 0.0},
                        (double[]){0.0, 0.0, 1.0, 1.0}, latlon_edges, 4),
      generate_cell_deg((double[]){-1.0,  0.0, 0.0, -1.0},
                        (double[]){-1.0, -1.0, 1.0,  1.0}, latlon_edges, 4),
      generate_cell_deg((double[]){1.0, 2.0, 2.0, 1.0},
                        (double[]){1.0, 1.0, 2.0, 2.0}, latlon_edges, 4)};
    enum {nSourceCells = sizeof(SourceCells)/sizeof(SourceCells[0])};

    double ref_area[nSourceCells];
    double ref_barycenter[nSourceCells][3];
    for (int i = 0; i < nSourceCells; ++i)
      for (int j = 0; j < 3; ++j)
        ref_barycenter[i][j] = 0.0;
    {
      ref_area[0] = yac_huiliers_area_info(TargetCells[0], ref_barycenter[0], 1.0);
      ref_area[1] = yac_huiliers_area_info(SourceCells[1], ref_barycenter[1], 1.0);
      ref_area[2] = yac_huiliers_area_info(TargetCells[0], ref_barycenter[2], 1.0);
      {
        double intersection[3], intersection_lon, intersection_lat;
        if (!intersect(GREAT_CIRCLE_EDGE, -0.5, 0.5, 0.5, 0.5,
                       GREAT_CIRCLE_EDGE,  0.0, 0.0, 0.0, 1.0, intersection))
          return EXIT_FAILURE;
        XYZtoLL(intersection, &intersection_lon, &intersection_lat);
        intersection_lon /= YAC_RAD;
        intersection_lat /= YAC_RAD;
        struct grid_cell overlap_cell =
          generate_cell_deg((double[]){0.0, 0.5, intersection_lon},
                            (double[]){0.0, 0.5, intersection_lat},
                            gc_edges, 3);
        ref_area[3] =
          yac_huiliers_area_info(overlap_cell, ref_barycenter[3], 1.0);
        yac_free_grid_cell(&overlap_cell);
        overlap_cell =
          generate_cell_deg((double[]){-0.5, -0.5, 0.0, intersection_lon},
                            (double[]){ 0.5, -0.5, 0.0, intersection_lat},
                            gc_edges, 4);
        ref_area[4] =
          yac_huiliers_area_info(overlap_cell, ref_barycenter[4], 1.0);
        yac_free_grid_cell(&overlap_cell);
      }
      ref_area[5] = 0.0;
      ref_barycenter[5][0] = 0.0;
      ref_barycenter[5][1] = 0.0;
      ref_barycenter[5][2] = 0.0;
    }

    for (int i = 0; i < nTargetCells; ++i) {

      double overlap_areas[nSourceCells];
      double overlap_barycenters[nSourceCells][3];
      yac_compute_overlap_info(
        nSourceCells, SourceCells, TargetCells[i],
        overlap_areas, &overlap_barycenters[0]);

      for (int i = 0; i < nSourceCells; ++i) {
        double area_tolerance = yac_huiliers_area(SourceCells[i]) * 1e-6;
        if (fabs(ref_area[i] - overlap_areas[i]) > area_tolerance)
          PUT_ERR("ERROR in yac_compute_overlap_info (area)");
        if (overlap_areas[i] > area_tolerance) {
          normalise_vector(ref_barycenter[i]);
          normalise_vector(overlap_barycenters[i]);
          if(get_vector_angle(
               ref_barycenter[i], overlap_barycenters[i]) > 1e-9)
            PUT_ERR("ERROR in yac_compute_overlap_info (barycenter)");
        }
      }
    }
    for (int i = 0; i < nSourceCells; ++i)
      yac_free_grid_cell(&SourceCells[i]);
    for (int i = 0; i < nTargetCells; ++i)
      yac_free_grid_cell(&TargetCells[i]);
  }

  return TEST_EXIT_CODE;
}

static void check_partial_areas(
  double * tgt_lon, double * tgt_lat,
  enum yac_edge_type * tgt_edges, int tgt_count,
  double ** src_lons, double ** src_lats,
  enum yac_edge_type ** src_edges, int * src_counts, size_t nSourceCells,
  struct grid_cell (*generate_cell)(
    double*, double*, enum yac_edge_type*,size_t)) {

  /* For the test at the Equator we cannot get better than 1.0e-11, while
     for the Pole test we are somewhat better.

    tolerance test for the deviation from 1 of  the sum of weights: 1.0e-11

  */
  double const epsilon = 1.0e-10; // relative precision

  struct grid_cell TargetCell =
    generate_cell(tgt_lon, tgt_lat, tgt_edges, tgt_count);

  struct grid_cell * SourceCells =
    xmalloc(nSourceCells * sizeof(*SourceCells));
  for (size_t i = 0; i < nSourceCells; ++i)
    SourceCells[i] =
      generate_cell(
        src_lons[i], src_lats[i], src_edges[i], src_counts[i]);

  double tgt_area = yac_huiliers_area(TargetCell);
  double area_tolerance = MAX(tgt_area * 1e-6, 1e-10);

  double partial_areas[nSourceCells];
  yac_compute_overlap_areas(
    nSourceCells, SourceCells, TargetCell, partial_areas);

  double partial_areas_sum = 0.0;
  double weights[nSourceCells];
  for (size_t i = 0; i < nSourceCells; ++i) {
    partial_areas_sum += partial_areas[i];
    weights[i] = partial_areas[i] / tgt_area;
  }

  // check partial areas
  if (fabs(partial_areas_sum - tgt_area) > area_tolerance)
    PUT_ERR("ERROR in sum of partial areas\n");

  // test yac_correct_weights
  yac_correct_weights(nSourceCells, weights);
  double weights_sum = 0.0;
  for (size_t i = 0; i < nSourceCells; ++i) weights_sum += weights[i];
  if ( fabs(weights_sum-1.0) > epsilon )
    PUT_ERR("ERROR in sum of corrected weights\n");

  for (size_t i = 0; i < nSourceCells; ++i)
    yac_free_grid_cell(&SourceCells[i]);
  yac_free_grid_cell(&TargetCell);
  free(SourceCells);
}

static void check_partial_areas_deg(
  double * tgt_lon, double * tgt_lat,
  enum yac_edge_type * tgt_edges, int tgt_count,
  double ** src_lons, double ** src_lats,
  enum yac_edge_type ** src_edges, int * src_counts, size_t nSourceCells) {

  check_partial_areas(
    tgt_lon, tgt_lat, tgt_edges, tgt_count,
    src_lons, src_lats, src_edges, src_counts, nSourceCells,
    generate_cell_deg);
}

static void check_partial_areas_rad(
  double * tgt_lon, double * tgt_lat,
  enum yac_edge_type * tgt_edges, int tgt_count,
  double ** src_lons, double ** src_lats,
  enum yac_edge_type ** src_edges, int * src_counts, size_t nSourceCells) {

  check_partial_areas(
    tgt_lon, tgt_lat, tgt_edges, tgt_count,
    src_lons, src_lats, src_edges, src_counts, nSourceCells,
    generate_cell_rad);
}
