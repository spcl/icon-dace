// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "test_cxc.h"
#include "geometry.h"

static void test_gcxlatc(double lon_a, double lat_a, double lon_b, double lat_b,
                         double lon_c, double lat_c, double lon_d, double lat_d,
                         double lon_ref_p, double lat_ref_p,
                         double lon_ref_q, double lat_ref_q, int ref_ret_val);

int main (void) {

   enum {
      no_intersection = 0,
      p_between_ab = 1,
      q_between_ab = 2,
      p_between_cd = 4,
      q_between_cd = 8,
      circles_are_identically = 16
   };

   // simple example

   test_gcxlatc(-10, 10, // point a
                 10, 10, // point b
                -10,  0, // point c
                 10,  0, // point d
                 90,  0, // reference point p
                -90,  0, // reference point q
                no_intersection); // reference return value

   // simple example

   test_gcxlatc(-10,  10, // point a
                 10, -10, // point b
                -10,   0, // point c
                 10,   0, // point d
                180,   0, // reference point p
                  0,   0, // reference point q
                q_between_ab + q_between_cd); // reference return value

   // example from test program

   test_gcxlatc(  0,  90, // point a
                140,  60, // point b
                 90, -45, // point c
                 39, -45, // point d
                140, -45, // reference point p
                -40, -45, // reference point q
                no_intersection); // reference return value

   // both circles do not intersect

   test_gcxlatc(-10, 10, // point a
                 10,  9, // point b
                -10, 70, // point c
                 10, 70, // point d
                -10, 70, // reference point p
                 10, 70, // reference point q
                -1);     // reference return value

   test_gcxlatc(-10,  9, // point a
                 10, 10, // point b
                -10, 70, // point c
                 10, 70, // point d
                -10, 70, // reference point p
                 10, 70, // reference point q
                -1);     // reference return value

   // the great circle only touches the circle of latitude

   test_gcxlatc(-90,  0, // point a
                  0, 45, // point b
                -10, 45, // point c
                 10, 45, // point d
                  0, 45, // reference point p
                  0, 45, // reference point q
                p_between_ab + p_between_cd +
                q_between_ab + q_between_cd); // reference return value

   // the circle of latitude is on the pole and does not intersect the great circle

   test_gcxlatc(-90,  0, // point a
                  0, 45, // point b
                 10, 90, // point c
                 20, 90, // point d
                 10, 90, // reference point p
                 20, 90, // reference point q
                -1);     // reference return value

   // the circle of latitude is on the pole and does intersect the great circle

   test_gcxlatc(15      ,  85, // point a
                180 + 15,  85, // point b
                10      ,  90, // point c
                20      ,  90, // point d
                       0,  90, // reference point p
                       0, -90, // reference point q
                p_between_ab + p_between_cd); // reference return value

   // the circle of latitude is on the pole and does intersect the great circle

   test_gcxlatc(30,  90, // point a
                30,  70, // point b
                10,  90, // point c
                20,  90, // point d
                 0,  90, // reference point p
                 0, -90, // reference point q
                p_between_ab + p_between_cd); // reference return value

   // the great circle edge has a length of zero and has no intersection with
   // the circle of latitude

   test_gcxlatc(30, 70, // point a
                30, 70, // point b
                10, 60, // point c
                20, 60, // point d
                -1, -1, // reference point p
                -1, -1, // reference point q
                -1);    // reference return value

   // the great circle edge has a length of zero and an intersection with
   // the circle of latitude
   test_gcxlatc( 15, 10, // point a
                 15, 10, // point b
                 10, 10, // point c
                 20, 10, // point d
                 15, 10, // reference point p
                195, 10, // reference point q
                p_between_ab + p_between_cd); // reference return value

   // the latitude circle edge has a length of zero and has no intersection with
   // the great circle edge

   test_gcxlatc(-5, -5, // point a
                 5,  5, // point b
                10, 89, // point c
                10, 89, // point d
                -1, -1, // reference point p
                -1, -1, // reference point q
                -1);    // reference return value

   // the latitude circle edge has a length of zero and has no intersection with
   // the great circle edge

   test_gcxlatc(-5, -5, // point a
                 5,  5, // point b
                10,  0, // point c
                10,  0, // point d
                -1, -1, // reference point p
                -1, -1, // reference point q
                -1);    // reference return value

   // the latitude circle edge has a length of zero and has an intersection with
   // the great circle edge

   test_gcxlatc(   -5, -5, // point a
                    5,  5, // point b
                    0,  0, // point c
                    0,  0, // point d
                    0,  0, // reference point p
                0+180,  0, // reference point q
                p_between_ab + p_between_cd); // reference return value

   {
      // example taken from a toy...caused a crash

      double a[3] = {-0.18125814883147034, -0.5578552203226167, -0.80990310324198245};
      double b[3] = {-0.18125814875183416, -0.557855220392383, -0.80990310321175063};
      double c[3] = {-0.18673822182292288, -0.54232717509597339, -0.8191520442889918};
      double d[3] = {-0.1772448664054937, -0.5455036073850148, -0.8191520442889918};
      double p[3], q[3];

      yac_intersect_vec(YAC_GREAT_CIRCLE_EDGE, a, b, YAC_LAT_CIRCLE_EDGE, c, d, p, q);
   }

   {
      // example taken from a toy...caused a crash

      double a[3] = {0.049067638730680319, 0.0012045429209767452, -0.99879473161693588};
      double b[3] = {0.049067674327417689, 0, -0.99879545620517241};
      double c[3] = {0.049067674327418126, 0, -0.99879545620517241};
      double d[3] = {0.049052896061289972, 0.0012041810087570048, -0.99879545620517241};
      double p[3], q[3];

      yac_intersect_vec(YAC_GREAT_CIRCLE_EDGE, a, b, YAC_LAT_CIRCLE_EDGE, c, d, p, q);
   }

   // both circles are on the equator and they do intersect

   test_gcxlatc(-10, 0, // point a
                 10, 0, // point b
                -10, 0, // point c
                 10, 0, // point d
                -10, 0, // reference point p
                 10, 0, // reference point q
                p_between_ab + p_between_cd +
                q_between_ab + q_between_cd +
                circles_are_identically); // reference return value

   test_gcxlatc(  0, 0, // point a
                 10, 0, // point b
                -10, 0, // point c
                 10, 0, // point d
                  0, 0, // reference point p
                 10, 0, // reference point q
                p_between_ab + p_between_cd +
                q_between_ab + q_between_cd +
                circles_are_identically); // reference return value

   test_gcxlatc( -5, 0, // point a
                  5, 0, // point b
                -10, 0, // point c
                 10, 0, // point d
                 -5, 0, // reference point p
                  5, 0, // reference point q
                p_between_ab + p_between_cd +
                q_between_ab + q_between_cd +
                circles_are_identically); // reference return value

   test_gcxlatc(-10, 0, // point a
                 10, 0, // point b
                 -5, 0, // point c
                  5, 0, // point d
                 -5, 0, // reference point p
                  5, 0, // reference point q
                p_between_ab + p_between_cd +
                q_between_ab + q_between_cd +
                circles_are_identically); // reference return value

   // both circles are on the equator and they do not intersect

   test_gcxlatc(170, 0, // point a
                190, 0, // point b
                -10, 0, // point c
                 10, 0, // point d
                180, 0, // reference point p
                  0, 0, // reference point q
                p_between_ab + q_between_cd +
                circles_are_identically); // reference return value

   return TEST_EXIT_CODE;
}

static void test_gcxlatc(double lon_a, double lat_a, double lon_b, double lat_b,
                         double lon_c, double lat_c, double lon_d, double lat_d,
                         double lon_ref_p, double lat_ref_p,
                         double lon_ref_q, double lat_ref_q, int ref_ret_val) {

   test_cxc(lon_a, lat_a, lon_b, lat_b, lon_c, lat_c, lon_d, lat_d,
            YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
            lon_ref_p, lat_ref_p, lon_ref_q, lat_ref_q, ref_ret_val);
}

