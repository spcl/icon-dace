// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tests.h"
#include "geometry.h"
#include "grid_cell.h"
#include "test_common.h"

double const tol = 1.0e-10;

static void check_latlon_cell(double * coordinates_x, double * coordinates_y);
static void check_gc_triangle(double * coordinates_x, double * coordinates_y);
static void check_gc_quad(double * coordinates_x, double * coordinates_y);
// tests whether all corners of cell are within the bounding circle
static void test_circle(struct yac_grid_cell cell, struct bounding_circle circle);
static void test_cirumscribe_circle(struct yac_grid_cell cell, struct bounding_circle circle);
// tests whether a point is within the bounding circle
static unsigned point_in_circle(double point[3], struct bounding_circle circle);

int main (void) {

   { // test regular cell

      double coordinates_x[] = {-1.0, 1.0, 1.0, -1.0};
      double coordinates_y[] = {-1.0, -1.0, 1.0, 1.0};
      check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell

      double coordinates_x[] = {40.0, 45.0, 45.0, 40.0};
      double coordinates_y[] = {20.0, 20.0, 25.0, 25.0};
      check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell

      double coordinates_x[] = { 175.0, -175.0, -175.0, 175.0};
      double coordinates_y[] = {-5.0, -5.0, 5.0, 5.0};
      check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell


      double coordinates_x[] = {30.0, 40.0, 40.0, 30.0};
      double coordinates_y[] = {80.0, 80.0, 85.0, 85.0};
      check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell

      double coordinates_x[] = {30.0, 40.0, 40.0, 30.0};
      double coordinates_y[] = {80.0, 80.0, 90.0, 90.0};
      check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {-5.0, 5.0, 0.0};
      double coordinates_y[] = {-5.0, 5.0, -5.0};
      check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {0.0, 120.0, -120.0};
      double coordinates_y[] = {85.0, 85.0, 85.0};
      check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {0.0, 120.0, -120.0};
      double coordinates_y[] = {-85.0, -85.0, -85.0};
      check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {-5.0, 5.0, 1.0};
      double coordinates_y[] = {0.0, 0.0, 1.0};
      check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {0.0, 170.0, 260.0};
      double coordinates_y[] = {85.0, 85.0, 89.0};
      check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test great circle quad

      double coordinates_x[] = {-1.0, 1.0, 1.0, -1.0};
      double coordinates_y[] = {-1.0, -1.0, 1.0, 1.0};
      check_gc_quad(coordinates_x, coordinates_y);
   }

   { // test great circle quad

      double coordinates_x[] = {0.0, 90.0, 180.0, 270.0};
      double coordinates_y[] = {85.0, 85.0, 85.0, 85.0};
      check_gc_quad(coordinates_x, coordinates_y);
   }

   { // test great circle quad

      double coordinates_x[] = {0.0, 90.0, 180.0, 270.0};
      double coordinates_y[] = {-85.0, -85.0, -85.0, -85.0};
      check_gc_quad(coordinates_x, coordinates_y);
   }

   { // test great circle quad


      double coordinates_x[] = {0.0, 10.0, 0.0, -10.0};
      double coordinates_y[] = {-10.0, 0.0, 10.0, 0.0};
      check_gc_quad(coordinates_x, coordinates_y);
   }

   return TEST_EXIT_CODE;
}

static void check_gc_triangle(double * coordinates_x, double * coordinates_y) {

   enum yac_edge_type edges[] = {
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};

   for (int order = -1; order <= 1; order += 2) {

      for (int start = 0; start < 3; ++start) {

         double temp_coordinates_x[3];
         double temp_coordinates_y[3];
         double coords[3][3];

         for (int i = 0; i < 3; ++i) {
            temp_coordinates_x[i] = coordinates_x[(3+i*order+start)%3];
            temp_coordinates_y[i] = coordinates_y[(3+i*order+start)%3];
            LLtoXYZ(temp_coordinates_x[i]*YAC_RAD,
                    temp_coordinates_y[i]*YAC_RAD, coords[i]);
         }

         struct yac_grid_cell cell =
           generate_cell_deg(temp_coordinates_x, temp_coordinates_y, edges, 3);

         struct bounding_circle bnd_circle[3];

         yac_get_cell_bounding_circle(cell, bnd_circle+0);
         yac_get_cell_bounding_circle_unstruct_triangle(
            coords[0], coords[1], coords[2], bnd_circle+1);
         yac_get_cell_circumscribe_circle_unstruct_triangle(
            coords[0], coords[1], coords[2], bnd_circle+2);

         // if (bnd_circle[1].inc_angle > bnd_circle[0].inc_angle)
         if (compare_angles(
               bnd_circle[1].inc_angle, bnd_circle[0].inc_angle) > 0)
            PUT_ERR("get_cell_bounding_circle_unstruct_triangle did not "
                    "compute smallest bounding circle\n");

         test_circle(cell, bnd_circle[0]);
         test_circle(cell, bnd_circle[1]);
         test_circle(cell, bnd_circle[2]);
         test_cirumscribe_circle(cell, bnd_circle[2]);

         yac_free_grid_cell(&cell);
      }
   }
}

static void check_gc_quad(double * coordinates_x, double * coordinates_y) {

   enum yac_edge_type edges[] = {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
                                 YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};

   for (unsigned i = 0; i < 4; ++i) {
     coordinates_x[i] *= YAC_RAD;
     coordinates_y[i] *= YAC_RAD;
   }

   for (int order = -1; order <= 1; order += 2) {

      for (int start = 0; start < 4; ++start) {

         double temp_coordinates_x[4];
         double temp_coordinates_y[4];

         for (int i = 0; i < 4; ++i) {
            temp_coordinates_x[i] = coordinates_x[(4+i*order+start)%4];
            temp_coordinates_y[i] = coordinates_y[(4+i*order+start)%4];
         }

         struct yac_grid_cell cell =
           generate_cell_deg(temp_coordinates_x, temp_coordinates_y, edges, 4);

         struct bounding_circle bnd_circle;

         yac_get_cell_bounding_circle(cell, &bnd_circle);

         test_circle(cell, bnd_circle);

         yac_free_grid_cell(&cell);
      }
   }
}

static void check_latlon_cell(double * coordinates_x, double * coordinates_y) {

   for (int order = -1; order <= 1; order += 2) {

      for (int start = 0; start < 4; ++start) {

         double temp_coordinates_x[4];
         double temp_coordinates_y[4];
         enum yac_edge_type edges[4];
         double coords[4][3];

         for (int i = 0; i < 4; ++i) {
            temp_coordinates_x[i] = coordinates_x[(4+i*order+start)%4];
            temp_coordinates_y[i] = coordinates_y[(4+i*order+start)%4];
            LLtoXYZ(temp_coordinates_x[i]*YAC_RAD,
                    temp_coordinates_y[i]*YAC_RAD, coords[i]);
         }


         for (int i = 0; i < 4; ++i) {
            int temp[2] =
             {fabs(temp_coordinates_x[i] - temp_coordinates_x[(i+3)%4]) > 0.0,
              fabs(temp_coordinates_y[i] - temp_coordinates_y[(i+3)%4]) > 0.0};
            YAC_ASSERT(temp[0] != temp[1], "internal error")
            edges[i] = (temp[0])?(YAC_LON_CIRCLE_EDGE):(YAC_LAT_CIRCLE_EDGE);
         }

         struct yac_grid_cell cell =
           generate_cell_deg(temp_coordinates_x, temp_coordinates_y, edges, 4);

         struct bounding_circle bnd_circle[3];

         yac_get_cell_bounding_circle(cell, bnd_circle+0);
         yac_get_cell_bounding_circle_reg_quad(coords[0], coords[1],
                                               coords[2], bnd_circle+1);
         yac_get_cell_circumscribe_circle_reg_quad(coords[0], coords[1],
                                                   coords[2], bnd_circle+2);

         // if (bnd_circle[1].inc_angle > bnd_circle[0].inc_angle + tol)
         if (compare_angles(
               bnd_circle[1].inc_angle,
               sum_angles_no_check(bnd_circle[0].inc_angle, SIN_COS_TOL)) > 0)
            PUT_ERR("get_cell_bounding_circle_reg_quad did not "
                    "compute smallest bounding circle\n");

         test_circle(cell, bnd_circle[0]);
         test_circle(cell, bnd_circle[1]);
         test_circle(cell, bnd_circle[2]);
         test_cirumscribe_circle(cell, bnd_circle[2]);

         yac_free_grid_cell(&cell);
      }
   }
}

static void test_circle(struct yac_grid_cell cell, struct bounding_circle circle) {

   for (size_t i = 0; i < cell.num_corners; ++i)
      if (!point_in_circle(cell.coordinates_xyz[i], circle))
         PUT_ERR("point is not in bounding circle\n");
}

static void test_cirumscribe_circle(
  struct yac_grid_cell cell, struct bounding_circle circle) {

   for (size_t i = 0; i < cell.num_corners; ++i) {

      struct sin_cos_angle inc_angle =
        sum_angles_no_check(circle.inc_angle, SIN_COS_TOL);
      struct sin_cos_angle angle =
        get_vector_angle_2(circle.base_vector, cell.coordinates_xyz[i]);
      struct sin_cos_angle diff_angle;
      int neg = sub_angles(inc_angle, angle, &diff_angle);

      if (neg) PUT_ERR("not all points are on the circumscribe circle\n");
   }
}

static unsigned point_in_circle(double point[3], struct bounding_circle circle) {

   // double const tol = 1.0e-12;
   // double const pi = 3.14159265358979323846;

   struct sin_cos_angle inc_angle =
      sum_angles_no_check(
        *(struct sin_cos_angle*)&(circle.inc_angle), SIN_COS_TOL);

   // if (circle.inc_angle + tol >= pi) return 1 == 1;
   if (compare_angles(inc_angle, SIN_COS_M_PI) >= 0) return 1 == 1;

   return
      compare_angles(
         get_vector_angle_2(circle.base_vector, point), inc_angle) <= 0;
}

