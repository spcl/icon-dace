/**
 * @file test_pxgc.c
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

#include "test_cxc.h"

static void test_gcxgc(double lon_a, double lat_a, double lon_b, double lat_b,
                      double lon_c, double lat_c, double lon_d, double lat_d,
                      double lon_ref_p, double lat_ref_p,
                      double lon_ref_q, double lat_ref_q, int ref_ret_val);

int main (void) {

   enum {
      p_between_ab = 1,
      q_between_ab = 2,
      p_between_cd = 4,
      q_between_cd = 8,
      circles_are_identically = 16
   };

   double lon_middle_ab, lat_middle_ab, lon_middle_cd, lat_middle_cd;

   // simple example

   test_gcxgc(-10,   0, // point a
               10,   0, // point b
                0,  10, // point c
                0, -10, // point d
                0,   0, // reference point p
              180,   0, // reference point q
              p_between_ab + p_between_cd); // reference return value

   // d at north pole

   test_gcxgc(-10,   0, // point a
               10,   0, // point b
                0, -20, // point c
                0,  90, // point d
                0,   0, // reference point p
              180,   0, // reference point q
              p_between_ab + p_between_cd); // reference return value

   // edge_a == edge_b

   test_gcxgc(-10,   0, // point a
               10,   0, // point b
              -10,   0, // point c
               10,   0, // point d
              -10,   0, // reference point p
               10,   0, // reference point q
              p_between_ab + p_between_cd + q_between_ab + q_between_cd +
              circles_are_identically); // reference return value

   // both circles on the same plane but no intersection

   get_edge_middle_point(
     GREAT_CIRCLE_EDGE, -45, 45, 0, 0, &lon_middle_ab, &lat_middle_ab);
   get_edge_middle_point(
     GREAT_CIRCLE_EDGE, 225, 45, 135, -45, &lon_middle_cd, &lat_middle_cd);
   test_gcxgc(-45,  45,                     // point a
                0,   0,                     // point b
              225,  45,                     // point c
              135, -45,                     // point d
              lon_middle_ab, lat_middle_ab, // reference point p
              lon_middle_cd, lat_middle_cd, // reference point q
              p_between_ab + q_between_cd +
              circles_are_identically); // reference return value

   // normal case + some cyclic stuff

   test_gcxgc(  45,   3, // point a
                45, -10, // point b
               -45,   0, // point c
              -135, -45, // point d
                45,  45, // reference point p
              -135, -45, // reference point q
              q_between_cd); // reference return value

   // normal case + some heavier cyclic stuff

   test_gcxgc(  67.5,  40, // point a
                67.5, -13, // point b
                 405,   0, // point c
                 450,   0, // point d
              -112.5,   0, // reference point p
                67.5,   0, // reference point q
              q_between_ab + q_between_cd); // reference return value


   // normal case + some heavier cyclic stuff + shifted end of edge b

   test_gcxgc(  67.5,  40, // point a
                67.5, -13, // point b
                  45,   0, // point c
                 450,   0, // point d
              -112.5,   0, // reference point p
                67.5,   0, // reference point q
              q_between_ab + q_between_cd); // reference return value

  // one of the two edge has length zero and is not on the other plane

   test_gcxgc( 10,  10, // point a
               10,  10, // point b
              -10, -10, // point c
               10,  11, // point d
               -1,  -1, // reference point p
               -1,  -1, // reference point q
              -1); // reference return value

   test_gcxgc(-10, -10, // point a
               10,  11, // point b
               10,  10, // point c
               10,  10, // point d
               -1,  -1, // reference point p
               -1,  -1, // reference point q
              -1); // reference return value

  // one of the two edge has length zero, is not on the other edge, but on the other plane

   test_gcxgc( 10,   10, // point a
               10,   10, // point b
              -10,  -10, // point c
                0,    0, // point d
               10,   10, // reference point p
               190, -10, // reference point q
              p_between_ab); // reference return value

   test_gcxgc(-10,  -10, // point a
                0,    0, // point b
               10,   10, // point c
               10,   10, // point d
               10,   10, // reference point p
               190, -10, // reference point q
              p_between_cd); // reference return value

   test_gcxgc( 10,   10, // point a
               10,   10, // point b
              190,  -10, // point c
              170,   10, // point d
               10,   10, // reference point p
               190, -10, // reference point q
              p_between_ab + q_between_cd); // reference return value

   test_gcxgc(190,  -10, // point a
              170,   10, // point b
               10,   10, // point c
               10,   10, // point d
               10,   10, // reference point p
               190, -10, // reference point q
              p_between_cd + q_between_ab); // reference return value

  // one of the two edge has length zero and touches the other edge

   test_gcxgc( 10,  10, // point a
               10,  10, // point b
              -10, -10, // point c
               10,  10, // point d
               10,   10, // reference point p
               190, -10, // reference point q
              p_between_ab + p_between_cd); // reference return value

   test_gcxgc( 10,  10, // point a
               10,  10, // point b
              -10, -10, // point c
               10,  10, // point d
               10,   10, // reference point p
               190, -10, // reference point q
              p_between_ab + p_between_cd); // reference return value

  // one of the two edge has length zero and is on the other edge

   test_gcxgc(  0,   0, // point a
                0,   0, // point b
              -10, -10, // point c
               10,  10, // point d
                0,   0, // reference point p
              180,   0, // reference point q
              p_between_ab + p_between_cd); // reference return value

   test_gcxgc(-10, -10, // point a
               10,  10, // point b
                0,   0, // point c
                0,   0, // point d
                0,   0, // reference point p
              180,   0, // reference point q
              p_between_ab + p_between_cd); // reference return value

   // both edges have length zero and do not match

   test_gcxgc(-10, -10, // point a
              -10, -10, // point b
               10,  10, // point c
               10,  10, // point d
               -1,  -1, // reference point p
               -1,  -1, // reference point q
              -1); // reference return value

   test_gcxgc( 10,  10, // point a
               10,  10, // point b
              -10, -10, // point c
              -10, -10, // point d
               -1,  -1, // reference point p
               -1,  -1, // reference point q
              -1); // reference return value

   // both edges have length zero and do not match, but are directly opposite
   // of each other

   test_gcxgc(-10, -10, // point a
              -10, -10, // point b
              170,  10, // point c
              170,  10, // point d
               -1,  -1, // reference point p
               -1,  -1, // reference point q
              -1); // reference return value

   test_gcxgc(170,  10, // point a
              170,  10, // point b
              -10, -10, // point c
              -10, -10, // point d
               -1,  -1, // reference point p
               -1,  -1, // reference point q
              -1); // reference return value

   // both edges have length zero and match

   test_gcxgc(-10, -10, // point a
              -10, -10, // point b
              -10, -10, // point c
              -10, -10, // point d
              -10, -10, // reference point p
              170,  10, // reference point q
              p_between_ab + p_between_cd); // reference return value

   return TEST_EXIT_CODE;
}

static void test_gcxgc(double lon_a, double lat_a, double lon_b, double lat_b,
                       double lon_c, double lat_c, double lon_d, double lat_d,
                       double lon_ref_p, double lat_ref_p,
                       double lon_ref_q, double lat_ref_q, int ref_ret_val) {

   test_cxc(lon_a, lat_a, lon_b, lat_b, lon_c, lat_c, lon_d, lat_d,
            GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
            lon_ref_p, lat_ref_p, lon_ref_q, lat_ref_q, ref_ret_val);
}

