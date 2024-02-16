/**
 * @file test_circle.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
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

#include "clipping.h"
#include "geometry.h"
#include "test_common.h"
#include "tests.h"

static void check_compare_circles(struct yac_circle * circle_a, int a_is_lat,
                                  struct yac_circle * circle_b, int b_is_lat,
                                  int are_identical);

int main (void) {

  struct point_2d {double lon, lat;};

  { // test yac_circle_compare and yac_circle_contains_north_pole
    struct {
      struct point_2d a, b;
      enum yac_edge_type edge_type;
      int edge_ordering;
      int circle_idx;
      int contains_north_pole;
    } test_data[] =
    {{.a = {.lon = 45.0, .lat = 0.0}, // equator clockwise
      .b = {.lon = -45.0, .lat = 0.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .circle_idx = 0,
      .contains_north_pole = 1},
     {.a = {.lon = -45.0, .lat = 0.0}, // equator clockwise
      .b = {.lon = 45.0, .lat = 0.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .circle_idx = 0,
      .contains_north_pole = 1},
     {.a = {.lon = 45.0, .lat = 0.0}, // equator counter clockwise
      .b = {.lon = -45.0, .lat = 0.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .circle_idx = 1,
      .contains_north_pole = 0},
     {.a = {.lon = -45.0, .lat = 0.0}, // equator counter clockwise
      .b = {.lon = 45.0, .lat = 0.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .circle_idx = 1,
      .contains_north_pole = 0},
     {.a = {.lon = 0.0, .lat = 88.0}, // Prime meridian
      .b = {.lon = 180.0, .lat = 88.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .circle_idx = 2,
      .contains_north_pole = 0},
     {.a = {.lon = 180.0, .lat = 88.0}, // Prime meridian
      .b = {.lon = 0.0, .lat = 88.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .circle_idx = 2,
      .contains_north_pole = 0},
     {.a = {.lon = 0.0, .lat = -1.0}, // Prime meridian (other direction)
      .b = {.lon = 0.0, .lat = 1.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .circle_idx = 3,
      .contains_north_pole = 0},
     {.a = {.lon = 10.0, .lat = 10.0}, // lon circle
      .b = {.lon = 10.0, .lat = 20.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = -1,
      .circle_idx = 4,
      .contains_north_pole = 0},
     {.a = {.lon = 10.0, .lat = 10.0}, // lon circle
      .b = {.lon = 10.0, .lat = 20.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = 1,
      .circle_idx = 5,
      .contains_north_pole = 0},
     {.a = {.lon = 45.0, .lat = 89.0}, // lon circle
      .b = {.lon = 45.0+180.0, .lat = 85.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = 1,
      .circle_idx = 6,
      .contains_north_pole = 0},
     {.a = {.lon = 45.0, .lat = 89.0}, // lon circle
      .b = {.lon = 45.0+180.0, .lat = 85.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = -1,
      .circle_idx = 7,
      .contains_north_pole = 0},
     {.a = {.lon = 10.0, .lat = 70.0}, // lat circle
      .b = {.lon = 20.0, .lat = 70.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .circle_idx = 8,
      .contains_north_pole = 1},
     {.a = {.lon = 10.0, .lat = 70.0}, // lat circle
      .b = {.lon = 20.0, .lat = 70.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .circle_idx = 9,
      .contains_north_pole = 0},
     {.a = {.lon = -5.0, .lat = 15.0}, // lat circle
      .b = {.lon = 5.0, .lat = 15.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .circle_idx = 10,
      .contains_north_pole = 1}};
    enum {NUM_TESTS = sizeof(test_data) / sizeof(test_data[0])};

    struct yac_circle test_circles[NUM_TESTS];

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      double a[3], b[3];
      LLtoXYZ_deg(test_data[i].a.lon, test_data[i].a.lat, a);
      LLtoXYZ_deg(test_data[i].b.lon, test_data[i].b.lat, b);
      yac_circle_generate(
        a, b, test_data[i].edge_type, test_data[i].edge_ordering,
        &test_circles[i]);
    }
    for (size_t i = 0; i < NUM_TESTS; ++i) {
      if (yac_circle_contains_north_pole(&test_circles[i]) !=
          test_data[i].contains_north_pole)
        PUT_ERR("error in yac_circle_contains_north_pole");
      for (size_t j = 0; j < NUM_TESTS; ++j) {
        check_compare_circles(
          &test_circles[i], test_data[i].edge_type == LAT_CIRCLE_EDGE,
          &test_circles[j], test_data[j].edge_type == LAT_CIRCLE_EDGE,
          test_data[i].circle_idx == test_data[j].circle_idx);
      }
    }
  }

  { // test yac_circle_point_is_inside
    struct {
      struct point_2d a, b, point;
      enum yac_edge_type edge_type;
      int edge_ordering;
      int is_inside;
    } test_data[] =
    {{.a = {.lon = -45.0, .lat = 0.0},
      .b = {.lon = 45.0, .lat = 0.0},
      .point = {.lon = 0.0, .lat = 90.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .is_inside = 1},
     {.a = {.lon = -45.0, .lat = 0.0},
      .b = {.lon = 45.0, .lat = 0.0},
      .point = {.lon = 0.0, .lat = 90.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .is_inside = 0},
     {.a = {.lon = -45.0, .lat = 0.0},
      .b = {.lon = 45.0, .lat = 0.0},
      .point = {.lon = 0.0, .lat = 0.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .is_inside = 2},
     {.a = {.lon = -45.0, .lat = 0.0},
      .b = {.lon = 45.0, .lat = 0.0},
      .point = {.lon = 1.0, .lat = 1.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .is_inside = 1},
     {.a = {.lon = 15.0, .lat = 25.0},
      .b = {.lon = 15.0, .lat = 30.0},
      .point = {.lon = 10.0, .lat = 10.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = 1,
      .is_inside = 1},
     {.a = {.lon = 15.0, .lat = 25.0},
      .b = {.lon = 15.0, .lat = 30.0},
      .point = {.lon = 10.0, .lat = 10.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = -1,
      .is_inside = 0},
     {.a = {.lon = 15.0, .lat = 25.0},
      .b = {.lon = 15.0, .lat = 30.0},
      .point = {.lon = 10.0, .lat = 90.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = -1,
      .is_inside = 2},
     {.a = {.lon = 15.0, .lat = 25.0},
      .b = {.lon = 15.0, .lat = 30.0},
      .point = {.lon = 15.0, .lat = 60.0},
      .edge_type = LON_CIRCLE_EDGE,
      .edge_ordering = -1,
      .is_inside = 2},
     {.a = {.lon = 75.0, .lat = 75.0},
      .b = {.lon = 80.0, .lat = 75.0},
      .point = {.lon = 15.0, .lat = 80.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .is_inside = 1},
     {.a = {.lon = 75.0, .lat = 75.0},
      .b = {.lon = 80.0, .lat = 75.0},
      .point = {.lon = 15.0, .lat = 90.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .is_inside = 1},
     {.a = {.lon = 75.0, .lat = 75.0},
      .b = {.lon = 80.0, .lat = 75.0},
      .point = {.lon = 15.0, .lat = 80.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .is_inside = 0},
     {.a = {.lon = 75.0, .lat = 75.0},
      .b = {.lon = 80.0, .lat = 75.0},
      .point = {.lon = 15.0, .lat = 90.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = -1,
      .is_inside = 0},
     {.a = {.lon = 75.0, .lat = 75.0},
      .b = {.lon = 80.0, .lat = 75.0},
      .point = {.lon = 15.0, .lat = 75.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .edge_ordering = 1,
      .is_inside = 2}};
    enum {NUM_TESTS = sizeof(test_data) / sizeof(test_data[0])};

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      double point[3];
      LLtoXYZ_deg(test_data[i].point.lon,
                  test_data[i].point.lat, point);
      struct yac_circle test_circle;
      double a[3], b[3];
      LLtoXYZ_deg(test_data[i].a.lon, test_data[i].a.lat, a);
      LLtoXYZ_deg(test_data[i].b.lon, test_data[i].b.lat, b);
      yac_circle_generate(
        a, b, test_data[i].edge_type, test_data[i].edge_ordering,
        &test_circle);
      if (yac_circle_point_is_inside(point, &test_circle) !=
          test_data[i].is_inside)
        PUT_ERR("error in yac_circle_point_is_inside");
    }
  }

  { // test yac_circle_compare_distances
    struct {
      struct point_2d a, b, point_a, point_b;
      enum yac_edge_type edge_type;
      int compare;
    } test_data[] =
    {{.a = {.lon = -1.0, .lat = -1.0},
      .b = {.lon = 1.0, .lat = 1.0},
      .point_a = {.lon = -1.0, .lat = 1.0},
      .point_b = {.lon = -2.0, .lat = 2.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .compare = -1},
     {.a = {.lon = -1.0, .lat = -1.0},
      .b = {.lon = 1.0, .lat = 1.0},
      .point_a = {.lon = -2.0, .lat = 2.0},
      .point_b = {.lon = -1.0, .lat = 1.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .compare = 1},
     {.a = {.lon = -1.0, .lat = -1.0},
      .b = {.lon = 1.0, .lat = 1.0},
      .point_a = {.lon = -2.0, .lat = 2.0},
      .point_b = {.lon = 2.0, .lat = -2.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .compare = 0},
     {.a = {.lon = -1.0, .lat = -1.0},
      .b = {.lon = 1.0, .lat = 1.0},
      .point_a = {.lon = 0.0, .lat = 0.0},
      .point_b = {.lon = 2.0, .lat = -2.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .compare = -1},
     {.a = {.lon = -1.0, .lat = -1.0},
      .b = {.lon = 1.0, .lat = 1.0},
      .point_a = {.lon = 2.0, .lat = -2.0},
      .point_b = {.lon = 0.0, .lat = 0.0},
      .edge_type = GREAT_CIRCLE_EDGE,
      .compare = 1},
     {.a = {.lon = 2.0, .lat = -1.0},
      .b = {.lon = 2.0, .lat = 1.0},
      .point_a = {.lon = 1.0, .lat = 1.0},
      .point_b = {.lon = 4.0, .lat = -1.0},
      .edge_type = LON_CIRCLE_EDGE,
      .compare = -1},
     {.a = {.lon = 2.0, .lat = -1.0},
      .b = {.lon = 2.0, .lat = 1.0},
      .point_a = {.lon = 4.0, .lat = -1.0},
      .point_b = {.lon = 1.0, .lat = 1.0},
      .edge_type = LON_CIRCLE_EDGE,
      .compare = 1},
     {.a = {.lon = 2.0, .lat = -1.0},
      .b = {.lon = 2.0, .lat = 1.0},
      .point_a = {.lon = 3.0, .lat = -1.0},
      .point_b = {.lon = 1.0, .lat = 1.0},
      .edge_type = LON_CIRCLE_EDGE,
      .compare = 0},
     {.a = {.lon = -1.0, .lat = 60.0},
      .b = {.lon = 1.0, .lat = 60.0},
      .point_a = {.lon = 7.0, .lat = 59.0},
      .point_b = {.lon = -13.0, .lat = 62.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .compare = -1},
     {.a = {.lon = -1.0, .lat = 60.0},
      .b = {.lon = 1.0, .lat = 60.0},
      .point_a = {.lon = -13.0, .lat = 62.0},
      .point_b = {.lon = 7.0, .lat = 59.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .compare = 1},
     {.a = {.lon = -1.0, .lat = 60.0},
      .b = {.lon = 1.0, .lat = 60.0},
      .point_a = {.lon = 7.0, .lat = 59.0},
      .point_b = {.lon = -13.0, .lat = 61.0},
      .edge_type = LAT_CIRCLE_EDGE,
      .compare = 0},
     {.a = {.lon = 45.0, .lat = 45.0},
      .b = {.lon = 45.0, .lat = 45.0},
      .point_a = {.lon = 46.0, .lat = 46.0},
      .point_b = {.lon = 47.0, .lat = 47.0},
      .edge_type = POINT,
      .compare = -1},
     {.a = {.lon = 45.0, .lat = 45.0},
      .b = {.lon = 45.0, .lat = 45.0},
      .point_a = {.lon = 47.0, .lat = 47.0},
      .point_b = {.lon = 46.0, .lat = 46.0},
      .edge_type = POINT,
      .compare = 1},
     {.a = {.lon = 45.0, .lat = 45.0},
      .b = {.lon = 45.0, .lat = 45.0},
      .point_a = {.lon = 40.0, .lat = 40.0},
      .point_b = {.lon = -45.0, .lat = -45.0},
      .edge_type = POINT,
      .compare = -1},
     {.a = {.lon = 45.0, .lat = 45.0},
      .b = {.lon = 45.0, .lat = 45.0},
      .point_a = {.lon = -45.0, .lat = -45.0},
      .point_b = {.lon = 40.0, .lat = 40.0},
      .edge_type = POINT,
      .compare = 1},
     {.a = {.lon = 45.0, .lat = 45.0},
      .b = {.lon = 45.0, .lat = 45.0},
      .point_a = {.lon = 45.0, .lat = 90.0},
      .point_b = {.lon = 45.0, .lat = 0.0},
      .edge_type = POINT,
      .compare = 0}};
    enum {NUM_TESTS = sizeof(test_data) / sizeof(test_data[0])};

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      double point_a[3], point_b[3];
      LLtoXYZ_deg(test_data[i].point_a.lon,
                  test_data[i].point_a.lat, point_a);
      LLtoXYZ_deg(test_data[i].point_b.lon,
                  test_data[i].point_b.lat, point_b);
      struct yac_circle test_circle;
      double a[3], b[3];
      LLtoXYZ_deg(test_data[i].a.lon, test_data[i].a.lat, a);
      LLtoXYZ_deg(test_data[i].b.lon, test_data[i].b.lat, b);
      yac_circle_generate(
        a, b, test_data[i].edge_type, 1, &test_circle);
      if (yac_circle_compare_distances(point_a, point_b, &test_circle) !=
          test_data[i].compare)
        PUT_ERR("error in yac_circle_compare_distances");
    }
  }

  { // basic tests for yac_circle_intersect
    struct {
      struct {
        struct point_2d points[2];
        enum yac_edge_type edge_type;
      } circles[2];
      struct point_2d p, q;
      int ret_value;
    } test_data[] =
      {{.circles = {{.points = {{.lon = -1.0, .lat = -1.0},
                                {.lon =  1.0, .lat =  1.0}},
                    .edge_type = GREAT_CIRCLE_EDGE},
                    {.points = {{.lon = -1.0, .lat =  1.0},
                                {.lon =  1.0, .lat = -1.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 0.0},
        .q = {.lon = 180.0, .lat = 0.0},
        .ret_value = 2},
       {.circles = {{.points = {{.lon = 0.0, .lat = -1.0},
                                {.lon = 0.0, .lat =  1.0}},
                    .edge_type = GREAT_CIRCLE_EDGE},
                    {.points = {{.lon = 180.0, .lat = -2.0},
                                {.lon = 180.0, .lat =  2.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 0.0},
        .q = {.lon = 0.0, .lat = 0.0},
        .ret_value = -1},
       {.circles = {{.points = {{.lon = 0.0, .lat = -1.0},
                                {.lon = 0.0, .lat =  1.0}},
                    .edge_type = GREAT_CIRCLE_EDGE},
                    {.points = {{.lon = 6.0, .lat = 45.0},
                                {.lon = 9.0, .lat = 45.0}},
                    .edge_type = LAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 45.0},
        .q = {.lon = 180.0, .lat = 45.0},
        .ret_value = 2},
       {.circles = {{.points = {{.lon = 0.0, .lat = -1.0},
                                {.lon = 0.0, .lat =  1.0}},
                    .edge_type = GREAT_CIRCLE_EDGE},
                    {.points = {{.lon = 6.0, .lat = 90.0},
                                {.lon = 9.0, .lat = 90.0}},
                    .edge_type = LAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 90.0},
        .q = {.lon = 180.0, .lat = -90.0},
        .ret_value = 1},
       {.circles = {{.points = {{.lon = 1.0, .lat = 0.0},
                                {.lon = 2.0, .lat = 0.0}},
                    .edge_type = GREAT_CIRCLE_EDGE},
                    {.points = {{.lon = 6.0, .lat = 45.0},
                                {.lon = 9.0, .lat = 45.0}},
                    .edge_type = LAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 0.0},
        .q = {.lon = 0.0, .lat = 0.0},
        .ret_value = 0},
       {.circles = {{.points = {{.lon = 0.0, .lat = -1.0},
                                {.lon = 0.0, .lat =  1.0}},
                    .edge_type = LON_CIRCLE_EDGE},
                    {.points = {{.lon = 17.0, .lat = -1.0},
                                {.lon = 17.0, .lat =  1.0}},
                    .edge_type = LON_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat =  90.0},
        .q = {.lon = 0.0, .lat = -90.0},
        .ret_value = 2},
       {.circles = {{.points = {{.lon = 0.0, .lat = -1.0},
                                {.lon = 0.0, .lat =  1.0}},
                    .edge_type = LON_CIRCLE_EDGE},
                    {.points = {{.lon = 0.0, .lat = 45.0},
                                {.lon = 0.0, .lat = 50.0}},
                    .edge_type = LON_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat =  90.0},
        .q = {.lon = 0.0, .lat = -90.0},
        .ret_value = -1},
       {.circles = {{.points = {{.lon = 0.0, .lat = -1.0},
                                {.lon = 0.0, .lat =  1.0}},
                    .edge_type = LON_CIRCLE_EDGE},
                    {.points = {{.lon = -1.0, .lat = 45.0},
                                {.lon =  1.0, .lat = 45.0}},
                    .edge_type = LAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 45.0},
        .q = {.lon = 180.0, .lat = 45.0},
        .ret_value = 2},
       {.circles = {{.points = {{.lon = 0.0, .lat = -1.0},
                                {.lon = 0.0, .lat =  1.0}},
                    .edge_type = LON_CIRCLE_EDGE},
                    {.points = {{.lon = -1.0, .lat = -90.0},
                                {.lon =  1.0, .lat = -90.0}},
                    .edge_type = LAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = -90.0},
        .q = {.lon = 180.0, .lat = 90.0},
        .ret_value = 1},
       {.circles = {{.points = {{.lon = 0.0, .lat = 0.0},
                                {.lon = 0.0, .lat = 0.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}, // point
                    {.points = {{.lon = -1.0, .lat = 0.0},
                                {.lon =  1.0, .lat = 0.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 0.0},
        .q = {.lon = 180.0, .lat = 0.0},
        .ret_value = 1},
       {.circles = {{.points = {{.lon = 0.0, .lat = 1.0},
                                {.lon = 0.0, .lat = 1.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}, // point
                    {.points = {{.lon = -1.0, .lat = 0.0},
                                {.lon =  1.0, .lat = 0.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 0.0},
        .q = {.lon = 180.0, .lat = 0.0},
        .ret_value = 0},
       {.circles = {{.points = {{.lon = 0.0, .lat = 45.0},
                                {.lon = 0.0, .lat = 45.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}, // point
                    {.points = {{.lon = -1.0, .lat = 45.0},
                                {.lon =  1.0, .lat = 45.0}},
                    .edge_type = LAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 45.0},
        .q = {.lon = 180.0, .lat = 45.0},
        .ret_value = 1},
       {.circles = {{.points = {{.lon = 0.0, .lat = 0.0},
                                {.lon = 0.0, .lat = 0.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}, // point
                    {.points = {{.lon = -1.0, .lat = 45.0},
                                {.lon =  1.0, .lat = 45.0}},
                    .edge_type = LAT_CIRCLE_EDGE}},
        .p = {.lon = 0.0, .lat = 45.0},
        .q = {.lon = 180.0, .lat = 45.0},
        .ret_value = 0},
       {.circles = {{.points = {{.lon = 5.0, .lat = 5.0},
                                {.lon = 5.0, .lat = 5.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}, // point
                    {.points = {{.lon = 5.0, .lat = 5.0},
                                {.lon = 5.0, .lat = 5.0}},
                    .edge_type = GREAT_CIRCLE_EDGE}}, // point
        .p = {.lon = 5.0, .lat = 5.0},
        .q = {.lon = 185.0, .lat = -5.0},
        .ret_value = 1}};
    enum {NUM_TESTS = sizeof(test_data) / sizeof(test_data[0])};

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      struct yac_circle test_circles[2];
      for (int j = 0; j < 2; ++j) {
        double points_3d[2][3];
        for (int k = 0; k < 2; ++k)
          LLtoXYZ_deg(
            test_data[i].circles[j].points[k].lon,
            test_data[i].circles[j].points[k].lat, points_3d[k]);
        yac_circle_generate(
          points_3d[0], points_3d[1], test_data[i].circles[j].edge_type,
          1, &test_circles[j]);
      }
      double ref_p[3], ref_q[3], p[3], q[3];
      LLtoXYZ_deg(test_data[i].p.lon, test_data[i].p.lat, ref_p);
      LLtoXYZ_deg(test_data[i].q.lon, test_data[i].q.lat, ref_q);
      for (int j = 0; j < 2; ++j) {
        if (yac_circle_intersect(
              test_circles[j], test_circles[j^1], p, q) !=
            test_data[i].ret_value)
          PUT_ERR("error in yac_circle_intersect");
        if (test_data[i].ret_value > 0)
          if (!(((get_vector_angle(p, ref_p) < yac_angle_tol) &&
                (get_vector_angle(q, ref_q) < yac_angle_tol)) ||
                ((get_vector_angle(q, ref_p) < yac_angle_tol) &&
                (get_vector_angle(p, ref_q) < yac_angle_tol))))
            PUT_ERR("error in yac_circle_intersect");
      }
    }
  }

  return TEST_EXIT_CODE;
}

static void check_compare_circles(struct yac_circle * circle_a, int a_is_lat,
                                  struct yac_circle * circle_b, int b_is_lat,
                                  int are_identical) {

  int ret_a = yac_circle_compare(&circle_a, &circle_b);
  int ret_b = yac_circle_compare(&circle_b, &circle_a);
  if ((are_identical && (ret_a || ret_b)) ||
      (!are_identical &&
       ((!ret_a || !ret_b) ||
        (ret_a != - ret_b) ||
        ((a_is_lat != b_is_lat) &&
        ((ret_a > ret_b) ^ a_is_lat)))))
    PUT_ERR("error in yac_circle_compare");
}
