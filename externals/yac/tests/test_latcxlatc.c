/**
 * @file test_gcxgc.c
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

#include "test_cxc.h"

static void test_latcxlatc(double lon_a, double lat_a, double lon_b, double lat_b,
                           double lon_c, double lat_c, double lon_d, double lat_d,
                           double lon_ref_p, double lat_ref_p,
                           double lon_ref_q, double lat_ref_q, int ref_ret_val);

int main (void) {

   unsigned const p_between_ab = 1 << 0;
   unsigned const q_between_ab = 1 << 1;
   unsigned const p_between_cd = 1 << 2;
   unsigned const q_between_cd = 1 << 3;
   unsigned const circles_are_identically = 1 << 4;
   double lon_middle_ab, lat_middle_ab, lon_middle_cd, lat_middle_cd;

   // two circles of latitude on the same plane that intersect

   for (int i = 0; i < 2; ++i) {

      double lon_coords[2][2] = {{-10, 5}, {-5, 10}};
      test_latcxlatc(lon_coords[i][0],   30, // point a
                     lon_coords[i][1],   30, // point b
                     lon_coords[i^1][0], 30, // point c
                     lon_coords[i^1][1], 30, // point d
                     -5,                 30, // reference point p
                      5,                 30, // reference point q
                     p_between_ab + p_between_cd + q_between_ab + q_between_cd +
                     circles_are_identically); // reference return value
   }

   for (int i = 0; i < 2; ++i) {

      double lon_coords[2][2] = {{20, 25}, {30, 10}};
      test_latcxlatc(lon_coords[i][0],   -15, // point a
                     lon_coords[i][1],   -15, // point b
                     lon_coords[i^1][0], -15, // point c
                     lon_coords[i^1][1], -15, // point d
                     20,                 -15, // reference point p
                     25,                 -15, // reference point q
                     p_between_ab + p_between_cd + q_between_ab + q_between_cd +
                     circles_are_identically); // reference return value
   }

   // two circles of latitude on the same plane that do not intersect

   for (int i = 0; i < 2; ++i) {

      double lon_coords[2][2] = {{-80, -70}, {20, 30}};
      get_edge_middle_point(
        LAT_CIRCLE_EDGE, lon_coords[i][0], 10, lon_coords[i][1], 10,
        &lon_middle_ab, &lat_middle_ab);
      get_edge_middle_point(
        LAT_CIRCLE_EDGE, lon_coords[i^1][0], 10, lon_coords[i^1][1], 10,
        &lon_middle_cd, &lat_middle_cd);
      test_latcxlatc(lon_coords[i][0],   10, // point a
                     lon_coords[i][1],   10, // point b
                     lon_coords[i^1][0], 10, // point c
                     lon_coords[i^1][1], 10, // point d
                     lon_middle_ab,      10, // reference point p
                     lon_middle_cd,      10, // reference point q
                     p_between_ab + q_between_cd +
                     circles_are_identically); // reference return value
   }

   // two circles of latitude on different planes that do not intersect

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2] = {0, 5};
      test_latcxlatc(20, lat_coords[i],   // point a
                     30, lat_coords[i],   // point b
                     20, lat_coords[i^1], // point c
                     30, lat_coords[i^1], // point d
                     -1, -1,              // reference point p
                     -1, -1,              // reference point q
                     -1); // reference return value

      test_latcxlatc(-80, lat_coords[i],   // point a
                     -70, lat_coords[i],   // point b
                      20, lat_coords[i^1], // point c
                      30, lat_coords[i^1], // point d
                      -1, -1,              // reference point p
                      -1, -1,              // reference point q
                     -1); // reference return value
   }

   // two circles of latitude on the pole

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2] = {90, -90};
      test_latcxlatc(  0, lat_coords[i],           // point a
                      10, lat_coords[i],           // point b
                     180, lat_coords[i],           // point c
                     190, lat_coords[i],           // point d
                       0, lat_coords[i],           // reference point p
                       0, lat_coords[i^1],         // reference point q
                     p_between_ab + p_between_cd); // reference return value
   }

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2] = {90, -90};
      test_latcxlatc( 0, lat_coords[i],            // point a
                     10, lat_coords[i],            // point b
                      5, lat_coords[i],            // point c
                     15, lat_coords[i],            // point d
                      0, lat_coords[i],            // reference point p
                      0, lat_coords[i^1],          // reference point q
                     p_between_ab + p_between_cd); // reference return value
   }

   // two circles of latitude on different planes with length zero

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2] = {-10, 10};
      test_latcxlatc(45    , lat_coords[i],   // point a
                     45    , lat_coords[i],   // point b
                     90    , lat_coords[i^1], // point c
                     90    , lat_coords[i^1], // point d
                     -1    , -1,              // reference point p
                     -1    , -1,              // reference point q
                     -1); // reference return value
   }

   // two circles of latitude on the same plane that do not intersect
   // with length zero

   test_latcxlatc(45,   10, // point a
                  46,   10, // point b
                  90,   10, // point c
                  91,   10, // point d
                  45.5, 10, // reference point p
                  90.5, 10, // reference point q
                  p_between_ab + q_between_cd +
                  circles_are_identically); // reference return value

   // two circles of latitude on the same plane that do intersect
   // with length zero

   test_latcxlatc( 90,  10, // point a
                   90,  10, // point b
                   90,  10, // point c
                   90,  10, // point d
                   90,  10, // reference point p
                  -90, -10, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // two circles of latitude on the same plane, one has the length zero and
   // intersects with the other

   test_latcxlatc(  0,  10, // point a
                    0,  10, // point b
                   -5,  10, // point c
                    5,  10, // point d
                    0,  10, // reference point p
                  180,  10, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // two circles of latitude on the same plane, one has the length zero and
   // does not intersects with the other

   test_latcxlatc(180,  10, // point a
                  180,  10, // point b
                   -5,  10, // point c
                    5,  10, // point d
                  180,  10, // reference point p
                    0,  10, // reference point q
                  p_between_ab + q_between_cd); // reference return value

   test_latcxlatc( 90,  10, // point a
                   90,  10, // point b
                   -5,  10, // point c
                    5,  10, // point d
                   90,  10, // reference point p
                  -90,  10, // reference point q
                  p_between_ab); // reference return value

   return TEST_EXIT_CODE;
}

static void test_latcxlatc(double lon_a, double lat_a, double lon_b, double lat_b,
                           double lon_c, double lat_c, double lon_d, double lat_d,
                           double lon_ref_p, double lat_ref_p,
                           double lon_ref_q, double lat_ref_q, int ref_ret_val) {

   test_cxc(lon_a, lat_a, lon_b, lat_b, lon_c, lat_c, lon_d, lat_d,
            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            lon_ref_p, lat_ref_p, lon_ref_q, lat_ref_q, ref_ret_val);
}

