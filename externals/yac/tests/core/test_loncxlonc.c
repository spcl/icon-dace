// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "test_cxc.h"

static void test_loncxlonc(double lon_a, double lat_a, double lon_b, double lat_b,
                           double lon_c, double lat_c, double lon_d, double lat_d,
                           double lon_ref_p, double lat_ref_p,
                           double lon_ref_q, double lat_ref_q, int ref_ret_val);

int main (void) {

   int const p_between_ab = 1 << 0;
   int const q_between_ab = 1 << 1;
   int const p_between_cd = 1 << 2;
   int const q_between_cd = 1 << 3;
   int const circles_are_identically = 1 << 4;

   // two circles of longitude on the same plane that intersect

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{-10, 5}, {-5, 10}};
      test_loncxlonc(0,  lat_coords[i][0],   // point a
                     0,  lat_coords[i][1],   // point b
                     0,  lat_coords[i^1][0], // point c
                     0,  lat_coords[i^1][1], // point d
                     0, -5,                 // reference point p
                     0,  5,                 // reference point q
                     p_between_ab + p_between_cd + q_between_ab + q_between_cd +
                     circles_are_identically); // reference return value
   }

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{20, 25}, {30, 10}};
      test_loncxlonc(10,  lat_coords[i][0],   // point a
                     10,  lat_coords[i][1],   // point b
                     10,  lat_coords[i^1][0], // point c
                     10,  lat_coords[i^1][1], // point d
                     10, 20,                 // reference point p
                     10, 25,                 // reference point q
                     p_between_ab + p_between_cd + q_between_ab + q_between_cd +
                     circles_are_identically); // reference return value
   }

   // two circles of longitude on the same plane that do not intersect

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{-80, -70}, {20, 30}};
      test_loncxlonc(30, lat_coords[i][0],         // point a
                     30, lat_coords[i][1],         // point b
                     30, lat_coords[i^1][0],       // point c
                     30, lat_coords[i^1][1],       // point d
                     30, (lat_coords[i][0] +
                          lat_coords[i][1])*0.5,   // reference point p
                     30, (lat_coords[i^1][0] +
                          lat_coords[i^1][1])*0.5, // reference point q
                     p_between_ab + q_between_cd +
                     circles_are_identically);     // reference return value
   }

   // two circles of longitude on different planes that do not intersect

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{-80, -70}, {20, 30}};
      test_loncxlonc(-30,  lat_coords[i][0],   // point a
                     -30,  lat_coords[i][1],   // point b
                      30,  lat_coords[i^1][0], // point c
                      30,  lat_coords[i^1][1], // point d
                       0,  90,                 // reference point p
                       0, -90,                 // reference point q
                      0); // reference return value
   }

   // two circles of longitude on different planes that do not intersect
   // and one edge goes across the north pole

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{80, 85}, {70, 75}};
      test_loncxlonc(-45    ,  lat_coords[i][0],   // point a
                     -45+180,  lat_coords[i][1],   // point b
                      50    ,  lat_coords[i^1][0], // point c
                      50    ,  lat_coords[i^1][1], // point d
                      0     ,  90,                 // reference point p
                      0     , -90,                 // reference point q
                      p_between_ab); // reference return value
   }

   // two circles of longitude on different planes that do not intersect
   // and one circle touches the north pole

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{90, 80}, {-40, -30}};
      int ret_val[2] = {p_between_ab, p_between_cd};
      test_loncxlonc(60, lat_coords[i][0],   // point a
                     60, lat_coords[i][1],   // point b
                     90, lat_coords[i^1][0], // point c
                     90, lat_coords[i^1][1], // point d
                      0,  90,                // reference point p
                      0, -90,                // reference point q
                     ret_val[i]); // reference return value
   }

   // two circles of longitude on different planes that intersect at the
   // north pole

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{75, 80},{60, 85}};
      test_loncxlonc(60    ,  lat_coords[i][0],   // point a
                     60-180,  lat_coords[i][1],   // point b
                     90    ,  lat_coords[i^1][0], // point c
                     90-180,  lat_coords[i^1][1], // point d
                      0    ,  90,                 // reference point p
                      0    , -90,                 // reference point q
                      p_between_ab + p_between_cd); // reference return value
   }

   // two circles of longitude on different planes that intersect at the
   // south pole

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{-75, -80},{-60, -85}};
      test_loncxlonc(60    , lat_coords[i][0],   // point a
                     60-180, lat_coords[i][1],   // point b
                     90    , lat_coords[i^1][0], // point c
                     90-180, lat_coords[i^1][1], // point d
                      0    ,  90,                // reference point p
                      0    , -90,                // reference point q
                      q_between_ab + q_between_cd); // reference return value
   }

   for (int i = 0; i < 2; ++i) {

      double lat_coords[2][2] = {{-75, -80},{-60, -90}};
      test_loncxlonc(60    , lat_coords[i][0],   // point a
                     60-180, lat_coords[i][1],   // point b
                     90    , lat_coords[i^1][0], // point c
                     90-180, lat_coords[i^1][1], // point d
                      0    ,  90,                // reference point p
                      0    , -90,                // reference point q
                      q_between_ab + q_between_cd); // reference return value
   }

   // two circles of longitude on different planes that do not intersect
   // one of the circles has length zero

   test_loncxlonc( 0, -10, // point a
                   0,  10, // point b
                  10,   0, // point c
                  10,   0, // point d
                  -1,  -1, // reference point p
                  -1,  -1, // reference point q
                  -1);     // reference return value

   test_loncxlonc(10, -10, // point a
                  10,  10, // point b
                   0,   0, // point c
                   0,   0, // point d
                  -1,  -1, // reference point p
                  -1,  -1, // reference point q
                  -1);     // reference return value

   test_loncxlonc( 0,  10, // point a
                   0,  10, // point b
                  10,   0, // point c
                  10,   0, // point d
                  -1,  -1, // reference point p
                  -1,  -1, // reference point q
                  -1);     // reference return value

   test_loncxlonc( 0,   0, // point a
                   0,   0, // point b
                  10,  10, // point c
                  10,  10, // point d
                  -1,  -1, // reference point p
                  -1,  -1, // reference point q
                  -1);     // reference return value

   // two circles of longitude on the same planes but do not intersect
   // both circles have length zero

   test_loncxlonc( 0, 10, // point a
                   0, 10, // point b
                   0,  0, // point c
                   0,  0, // point d
                  -1, -1, // reference point p
                  -1, -1, // reference point q
                  -1);    // reference return value

   // two circles of longitude on the same planes
   // both circles have length zero and are at the same location

   test_loncxlonc(  0, 0, // point a
                    0, 0, // point b
                    0, 0, // point c
                    0, 0, // point d
                    0, 0, // reference point p
                  180, 0, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // two circles of longitude, both have length zero and are on the poles

   for (int i = 0; i < 8; ++i) {
      double lon_coords[2] = {0 , 90};
      double lat_coords[2] = {-90, 90};

      test_loncxlonc(lon_coords[i>>2], lat_coords[i&1],          // point a
                     lon_coords[i>>2], lat_coords[i&1],          // point b
                     lon_coords[(i>>2)^1], lat_coords[(i>>1)&1], // point c
                     lon_coords[(i>>2)^1], lat_coords[(i>>1)&1], // point d
                      0,  lat_coords[i&1],                       // reference point p
                      0, -lat_coords[i&1],                       // reference point q
                     (((i&1)==((i>>1)&1))?
                       (p_between_ab+p_between_cd):
                       (-1))); // reference return value
   }

   // two circles of longitude, one has length zero and is on the pole

   for (int j = 0; j < 2; ++j) {
      double lat_coords[2][2] = {{90, 90}, {-10, 10}};
      int ref_ret_val[2] = {p_between_ab, p_between_cd};
      for (int i = 0; i < 2; ++i) {
         double lon_coords[2] = {0, 10};
         test_loncxlonc(lon_coords[i],   lat_coords[j][0],   // point a
                        lon_coords[i],   lat_coords[j][1],   // point b
                        lon_coords[i^1], lat_coords[j^1][0], // point c
                        lon_coords[i^1], lat_coords[j^1][1], // point d
                        0,                90,                // reference point p
                        0,               -90,                // reference point q
                        ref_ret_val[j]); // reference return value
      }
   }

   // two circles that intersect at the pole

   for (int i = 0; i < 4; ++i) {
      double lat_coords[2][2] = {{-89, 89},{-90, 90}};
      int ret_val[2] = {p_between_ab+q_between_cd,
                        p_between_ab+p_between_cd};
      test_loncxlonc(  0,  lat_coords[0][i&1], // point a
                     180,  lat_coords[0][i&1], // point b
                      45, lat_coords[0][i>>1], // point c
                     225, lat_coords[0][i>>1], // point d
                       0,  lat_coords[1][i&1], // reference point p
                       0, -lat_coords[1][i&1], // reference point q
                     ret_val[(i&1)==(i>>1)]);
   }

   for (int i = 0; i < 4; ++i) {
      double lat_coords[2][2] = {{-85, 85}, {-90, 90}};
      int ret_val[2] = {p_between_ab + q_between_cd,
                        p_between_ab + p_between_cd};
      test_loncxlonc(  0,  lat_coords[0][i&1], // point a
                     180,  lat_coords[0][i&1], // point b
                      45, lat_coords[0][i>>1], // point c
                     225, lat_coords[0][i>>1], // point d
                       0, lat_coords[1][i&1],  // reference point p
                       0, lat_coords[1][(i&1)^1], // reference point q
                     ret_val[(i&1)==(i>>1)]);
   }

   // two edges on identical circles that cross different poles

   for (int j = 0; j < 4; ++j) {
      double lon_coords[2] = {0, 180};
      for (int i = 0; i < 2; ++i) {
         double lat_coords[2][2] = {{-85, 85}, {-90, 90}};
         test_loncxlonc(lon_coords[j&1],      lat_coords[0][i&1],     // point a
                        lon_coords[(j&1)^1],  lat_coords[0][i&1],     // point b
                        lon_coords[j>>1],     lat_coords[0][(i&1)^1], // point c
                        lon_coords[(j>>1)^1], lat_coords[0][(i&1)^1], // point d
                        0,                    lat_coords[1][i&1],     // reference point p
                        0,                    lat_coords[1][(i&1)^1], // reference point q
                        p_between_ab + q_between_cd + circles_are_identically);
      }
   }

   // two edges on different circles that touch at the poles

   for (int j = 0; j < 2; ++j) {
      double lon_coords[2] = {0, 60};
      for (int i = 0; i < 2; ++i) {
         double lat_coords[2][2] = {{-85, -90}, {85, 90}};
         test_loncxlonc(lon_coords[j],   lat_coords[i][0],   // point a
                        lon_coords[j],   lat_coords[i][1],   // point b
                        lon_coords[j^1], lat_coords[i][0],   // point c
                        lon_coords[j^1], lat_coords[i][1],   // point d
                        0,               lat_coords[i][1],   // reference point p
                        0,               lat_coords[i^1][1], // reference point q
                        p_between_ab + p_between_cd);
      }
   }

   // two edges on identical circles that touch at the poles

   for (int i = 0; i < 2; ++i) {
      double lat_coords[2][2] = {{-85, -90}, {85, 90}};
      test_loncxlonc(60,     lat_coords[i][0], // point a
                     60,     lat_coords[i][1], // point b
                     60+180, lat_coords[i][0], // point c
                     60+180, lat_coords[i][1], // point d
                     0,      lat_coords[i][1], // reference point p
                     0,      lat_coords[i][1], // reference point q
                     p_between_ab + p_between_cd +
                     q_between_ab + q_between_cd +
                     circles_are_identically);
   }

   // two edges on identical circles that overlap on the poles

   for (int i = 0; i < 2; ++i) {
      double lat_coords[2][2] = {{-85, -80}, {85, 80}};
      test_loncxlonc(60,     lat_coords[i][0], // point a
                     60+180, lat_coords[i][1], // point b
                     60+180, lat_coords[i][0], // point c
                     60,     lat_coords[i][1], // point d
                     60,     lat_coords[i][0], // reference point p
                     60+180, lat_coords[i][0], // reference point q
                     p_between_ab + p_between_cd +
                     q_between_ab + q_between_cd + circles_are_identically);
   }

   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
      double lat_coords[2][2] = {{-85, -80}, {85, 80}};
      test_loncxlonc(60,     lat_coords[i][j],   // point a
                     60+180, lat_coords[i][j],   // point b
                     60+180, lat_coords[i][j^1], // point c
                     60,     lat_coords[i][j^1], // point d
                     60,     lat_coords[i][0],   // reference point p
                     60+180, lat_coords[i][0],   // reference point q
                     p_between_ab + p_between_cd +
                     q_between_ab + q_between_cd + circles_are_identically);
      }
   }

   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
      double lat_coords[2][3] = {{-80, -85, -90}, {80, 85, 90}};
      test_loncxlonc(60, lat_coords[i][j],   // point a
                     60, lat_coords[i][2],   // point b
                     60, lat_coords[i][j^1], // point c
                     60, lat_coords[i][2],   // point d
                     60, lat_coords[i][1],   // reference point p
                     60, lat_coords[i][2],   // reference point q
                     p_between_ab + p_between_cd +
                     q_between_ab + q_between_cd + circles_are_identically);
      }
   }

   for (int i = 0; i < 2; ++i) {
   double lat_coords[2][2] = {{-85, -90}, {85, 90}};
   test_loncxlonc(60, lat_coords[i][0],   // point a
                  60, lat_coords[i][1],   // point b
                  60, lat_coords[i][0], // point c
                  60, lat_coords[i][1], // point d
                  60, lat_coords[i][0],   // reference point p
                  60, lat_coords[i][1],   // reference point q
                  p_between_ab + p_between_cd +
                  q_between_ab + q_between_cd + circles_are_identically);
   }

   return TEST_EXIT_CODE;
}

static void test_loncxlonc(double lon_a, double lat_a, double lon_b, double lat_b,
                           double lon_c, double lat_c, double lon_d, double lat_d,
                           double lon_ref_p, double lat_ref_p,
                           double lon_ref_q, double lat_ref_q, int ref_ret_val) {

  test_cxc(lon_a, lat_a, lon_b, lat_b, lon_c, lat_c, lon_d, lat_d,
           YAC_LON_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE,
           lon_ref_p, lat_ref_p, lon_ref_q, lat_ref_q, ref_ret_val);
}
