// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "test_cxc.h"

static void test_loncxlatc(double lon_a, double lat_a, double lon_b, double lat_b,
                           double lon_c, double lat_c, double lon_d, double lat_d,
                           double lon_ref_p, double lat_ref_p,
                           double lon_ref_q, double lat_ref_q, int ref_ret_val);

int main (void) {

   unsigned const p_between_ab = 1 << 0;
   unsigned const q_between_ab = 1 << 1;
   unsigned const p_between_cd = 1 << 2;

   // circles that intersect

   test_loncxlatc(  0,  10, // point a
                    0, -10, // point b
                   10,   0, // point c
                  -10,   0, // point d
                    0,   0, // reference point p
                  180,   0, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // circles that intersect

   test_loncxlatc(65,     35, // point a
                  65,     45, // point b
                  60,     40, // point c
                  70,     40, // point d
                  65,     40, // reference point p
                  65+180, 40, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // circles that intersect but the edges not

   test_loncxlatc(-30,     -60, // point a
                  -30,     -70, // point b
                  -31,     -65, // point c
                  -40,     -65, // point d
                  -30,     -65, // reference point p
                  -30+180, -65, // reference point q
                  p_between_ab); // reference return value

   // circles that intersect and the edges touch

   test_loncxlatc(-30,     -60, // point a
                  -30,     -70, // point b
                  -30,     -65, // point c
                  -40,     -65, // point d
                  -30,     -65, // reference point p
                  -30+180, -65, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // circles that intersect but the edges not

   test_loncxlatc(-45,     80, // point a
                  -45,     90, // point b
                  -40,     70, // point c
                  -50,     70, // point d
                  -45,     70, // reference point p
                  -45+180, 70, // reference point q
                  p_between_cd); // reference return value

   // circles that intersect and the edges touch

   test_loncxlatc(-45,     70, // point a
                  -45,     90, // point b
                  -40,     70, // point c
                  -50,     70, // point d
                  -45,     70, // reference point p
                  -45+180, 70, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // circles that intersect at the north pole and the lon
   // edge crosses the north pole

   test_loncxlatc(-45,     80, // point a
                  -45+180, 80, // point b
                   20,     90, // point c
                   30,     90, // point d
                    0,     90, // reference point p
                    0,    -90, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // circles that intersect at the south pole and the lon
   // edge crosses the south pole

   test_loncxlatc(-45,     -80, // point a
                  -45+180, -80, // point b
                   20,     -90, // point c
                   30,     -90, // point d
                    0,     -90, // reference point p
                    0,      90, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // circles that intersect at the north pole and the lon
   // edge does not cross the north pole

   test_loncxlatc(-45,  80, // point a
                  -45,  85, // point b
                   20,  90, // point c
                   30,  90, // point d
                    0,  90, // reference point p
                    0, -90, // reference point q
                  p_between_cd); // reference return value

   // circles that intersect at the south pole and the lon
   // edge does not cross the south pole

   test_loncxlatc(-45, -80, // point a
                  -45, -85, // point b
                   20, -90, // point c
                   30, -90, // point d
                    0, -90, // reference point p
                    0,  90, // reference point q
                  p_between_cd); // reference return value

   // circles that intersect at the south pole and the lon
   // edge touches the south pole

   test_loncxlatc(-45, -80, // point a
                  -45, -90, // point b
                   20, -90, // point c
                   30, -90, // point d
                    0, -90, // reference point p
                    0,  90, // reference point q
                  p_between_ab + p_between_cd); // reference return value


   // circles that intersect close to the south pole and the lon
   // edge intersects the lat circle twice

   test_loncxlatc(  0, -80, // point a
                  180, -80, // point b
                   50, -89, // point c
                   60, -89, // point d
                    0, -89, // reference point p
                  180, -89, // reference point q
                  p_between_ab + q_between_ab); // reference return value

   test_loncxlatc(  0, -80, // point a
                  180, -80, // point b
                  -10, -89, // point c
                   10, -89, // point d
                    0, -89, // reference point p
                  180, -89, // reference point q
                  p_between_ab + p_between_cd + q_between_ab); // reference return value

   // longitude circle consists of two identical points and is not on
   // the latitude circle

   test_loncxlatc( 0,  0, // point a
                   0,  0, // point b
                   5, 10, // point c
                   5, 10, // point d
                  -1, -1, // reference point p
                  -1, -1, // reference point q
                  -1);    // reference return value

   test_loncxlatc( 0,  0, // point a
                   0,  0, // point b
                   5, 10, // point c
                  10, 10, // point d
                  -1, -1, // reference point p
                  -1, -1, // reference point q
                  -1);    // reference return value

   // longitude circle consists of two identical points and is on the
   // latitude circle

   test_loncxlatc(  0, 0, // point a
                    0, 0, // point b
                    5, 0, // point c
                   -5, 0, // point d
                    0, 0, // reference point p
                  180, 0, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   // latitude circle consists of two identical points and is not on
   // the longitude circle

   test_loncxlatc(  0, -5, // point a
                    0,  5, // point b
                    5,  0, // point c
                    5,  0, // point d
                   -1, -1, // reference point p
                   -1, -1, // reference point q
                   -1); // reference return value

   test_loncxlatc( 0, -5, // point a
                   0,  5, // point b
                   5, 10, // point c
                   5, 10, // point d
                  -1, -1, // reference point p
                  -1, -1, // reference point q
                  -1); // reference return value

   // longitude circle consists of two identical points and is on the
   // latitude circle

   test_loncxlatc(  0, -5, // point a
                    0,  5, // point b
                    0,  0, // point c
                    0,  0, // point d
                    0,  0, // reference point p
                  180,  0, // reference point q
                  p_between_ab + p_between_cd); // reference return value

   test_loncxlatc(  0,  90, // point a
                  140,  60, // point b
                   90, -45, // point c
                   39, -45, // point d
                  140, -45, // reference point p
                  -40, -45, // reference point q
                    0); // reference return value

   // both edges have length zero and are on a pole

   for (int i = 0; i < 4; ++i) {

      double lat_coords[2] = {90, -90};
      int ret_val[2] = {-1, p_between_ab + p_between_cd};
      test_loncxlatc( 20,  lat_coords[i&1],   // point a
                     200,  lat_coords[i&1],   // point b
                      10,  lat_coords[i>>1],  // point c
                      20,  lat_coords[i>>1],  // point d
                       0,  lat_coords[i>>1],  // reference point p
                       0, -lat_coords[i>>1],  // reference point q
                     ret_val[(i&1)==(i>>1)]); // reference return value
   }

   {

      double a[3] = {-0.01728256081754178,
                     -0.0024289055230921597,
                      0.99984769515639127};
      double b[3] = {-0.034559857199638548,
                     -0.0048570711780326408,
                      0.99939082701909576};
      double c[3] = {-0.017298289573378096,
                     -0.0023142316840916847,
                      0.99984769515639127};
      double d[3] = { 0.01732231884620395,
                      0.0021269133133512588,
                      0.99984769515639127};

      double a_lon, a_lat, b_lon, b_lat, c_lon, c_lat, d_lon, d_lat;

      XYZtoLL(a, &a_lon, &a_lat);
      XYZtoLL(b, &b_lon, &b_lat);
      XYZtoLL(c, &c_lon, &c_lat);
      XYZtoLL(d, &d_lon, &d_lat);

      test_loncxlatc(a_lon/YAC_RAD, a_lat/YAC_RAD, // point a
                     b_lon/YAC_RAD, b_lat/YAC_RAD, // point b
                     c_lon/YAC_RAD, c_lat/YAC_RAD, // point c
                     d_lon/YAC_RAD, d_lat/YAC_RAD, // point d
                     a_lon/YAC_RAD, c_lat/YAC_RAD, // reference point p
                     a_lon/YAC_RAD+180.0, c_lat/YAC_RAD, // reference point q
                     p_between_ab + p_between_cd); // reference return value
   }

   return TEST_EXIT_CODE;
}

static void test_loncxlatc(double lon_a, double lat_a, double lon_b, double lat_b,
                           double lon_c, double lat_c, double lon_d, double lat_d,
                           double lon_ref_p, double lat_ref_p,
                           double lon_ref_q, double lat_ref_q, int ref_ret_val) {

   test_cxc(lon_a, lat_a, lon_b, lat_b, lon_c, lat_c, lon_d, lat_d,
            YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
            lon_ref_p, lat_ref_p, lon_ref_q, lat_ref_q, ref_ret_val);
}

