// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "basic_grid.h"
#include "clipping.h"
#include "geometry.h"
#include "tests.h"
#include "area.h"
#include "test_common.h"

static void test_lat_clipping(struct yac_grid_cell cells,
                              double deg_lat_bounds[2],
                              struct yac_grid_cell * ref_cells,
                              unsigned num_ref_cells);

static int compare_cells(struct yac_grid_cell a, struct yac_grid_cell b);

static double get_intersection(
  enum yac_edge_type edge_type, double points[2][2], double lat) {

  double intersection[3], intersection_lon, intersection_lat;

  if (!intersect(edge_type, points[0][0], points[0][1],
                            points[1][0], points[1][1],
                 YAC_LAT_CIRCLE_EDGE, points[0][0], lat,
                                  points[1][0], lat, intersection))
    exit(EXIT_FAILURE);
  XYZtoLL(intersection, &intersection_lon, &intersection_lat);
  return intersection_lon / YAC_RAD;
}

int main (void) {

  enum yac_edge_type gc_edges[] = {
    YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
    YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
    YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
    YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};
  enum yac_edge_type latlon_edges[] = {
    YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE};

  double temp_lon[16];

  // a gaussian cell outside the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){-5,-5,5,5}, latlon_edges, 4),
    (double[2]){10, 20}, NULL, 0);

  // a gaussian cell outside the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){-5,-5,10,10}, latlon_edges, 4),
    (double[2]){10, 20}, NULL, 0);

  // a gaussian cell outside the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){20,20,30,30}, latlon_edges, 4),
    (double[2]){10, 20}, NULL, 0);

  // a gaussian cell outside the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){80,80,90,90}, latlon_edges, 4),
    (double[2]){10, 20}, NULL, 0);

  // a gaussian cell outside the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){-80,-80,-90,-90}, latlon_edges, 4),
    (double[2]){10, 20}, NULL, 0);

  // a gaussian cell outside the lat bound (one bound is at a pole)
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){80, 90}, NULL, 0);

  // a gaussian cell outside the lat bound (both bounds are identical)
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){0, 0}, NULL, 0);

  // a gaussian cell outside the lat bound (both bounds are identical)
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){-10,-10,-20,-20}, latlon_edges, 4),
    (double[2]){0, 0}, NULL, 0);

  // a gaussian cell outside the lat bound (both bounds are identical)
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){-5,-5,5,5}, latlon_edges, 4),
    (double[2]){0, 0}, NULL, 0);

  // a gaussian cell overlapping with the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){-10,-10,10,10}, latlon_edges, 4),
    (double[2]){-5, 5},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){-5,5,5,-5},
       (double[4]){-5,-5,5,5}, latlon_edges, 4)}, 1);

  // a gaussian cell overlapping with the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){12, 15},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){12,12,15,15}, latlon_edges, 4)}, 1);

  // a gaussian cell overlapping with the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){15, 25},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){15,15,20,20}, latlon_edges, 4)}, 1);

  // a gaussian cell overlapping with the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){15, 5},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){15,15,10,10}, latlon_edges, 4)}, 1);

  // a gaussian cell overlapping with the lat bound
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){5, 25},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){10,10,20,20}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){90, 80}, NULL, 0);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){90,90,70,70}, latlon_edges, 4),
    (double[2]){90, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){90,90,80,80}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){90,90,85,85}, latlon_edges, 4),
    (double[2]){90, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){90,90,85,85}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){88,88,83,83}, latlon_edges, 4),
    (double[2]){90, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){88,88,83,83}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){88,88,83,83}, latlon_edges, 4),
    (double[2]){90, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){88,88,83,83}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){88,88,70,70}, latlon_edges, 4),
    (double[2]){90, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){88,88,80,80}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){10,10,20,20}, latlon_edges, 4),
    (double[2]){-90, -80}, NULL, 0);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-90,-90,-70,-70}, latlon_edges, 4),
    (double[2]){-90, -80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-90,-90,-80,-80}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-90,-90,-85,-85}, latlon_edges, 4),
    (double[2]){-90, -80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-90,-90,-85,-85}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-88,-88,-83,-83}, latlon_edges, 4),
    (double[2]){-90, -80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-88,-88,-83,-83}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-88,-88,-83,-83}, latlon_edges, 4),
    (double[2]){-90, -80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-88,-88,-83,-83}, latlon_edges, 4)}, 1);

  // a gaussian cell
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-88,-88,-70,-70}, latlon_edges, 4),
    (double[2]){-90, -80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[4]){40,45,45,40},
       (double[4]){-88,-88,-80,-80}, latlon_edges, 4)}, 1);

  // a gaussian cell that touches the north pole and partially overlaps
  // with the bounds
  test_lat_clipping(
    generate_cell_deg(
       (double[4]){0,5,5,0},
       (double[4]){90,90,85,85}, latlon_edges, 4),
    (double[2]){89.9, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
        (double[4]){0,5,5,0},
        (double[4]){89.9,89.9,85,85}, latlon_edges, 4)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[3]){10,20,15},
       (double[3]){30,30,40}, gc_edges, 3),
    (double[2]){35, 50},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){15,
                  get_intersection(
                     YAC_GREAT_CIRCLE_EDGE,
                     (double[2][2]){{15, 40}, {10, 30}}, 35),
                  get_intersection(
                     YAC_GREAT_CIRCLE_EDGE,
                     (double[2][2]){{15, 40}, {20, 30}}, 35)},
       (double[]){40,35,35},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE}, 3)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[3]){10,20,15},
       (double[3]){30,30,40}, gc_edges, 3),
    (double[2]){35, 25},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){20, 10,
                  get_intersection(
                     YAC_GREAT_CIRCLE_EDGE,
                     (double[2][2]){{15, 40}, {10, 30}}, 35),
                  get_intersection(
                     YAC_GREAT_CIRCLE_EDGE,
                     (double[2][2]){{15, 40}, {20, 30}}, 35)},
       (double[]){30,30,35,35},
       (enum yac_edge_type[]){YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
                              YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE}, 4)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){10,20,15},
       (double[]){30,30,40}, gc_edges, 3),
    (double[2]){45, 25},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){10,20,15},
       (double[]){30,30,40}, gc_edges, 3)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[3]){10,20,15},
       (double[3]){30,30,40}, gc_edges, 3),
    (double[2]){32, 38},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{15, 40}, {20, 30}}, 32),
                  get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{15, 40}, {20, 30}}, 38),
                  get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{15, 40}, {10, 30}}, 38),
                  get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{15, 40}, {10, 30}}, 32)},
       (double[]){32,38,38,32},
       (enum yac_edge_type[]){YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
                              YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE}, 4)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[3]){-20,20,0},
       (double[3]){80,80,85}, gc_edges, 3),
    (double[2]){80.001, 89},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){0,
                  get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{0, 85}, {-20, 80}}, 80.001),
                  -fabs(get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{-20, 80}, {20, 80}}, 80.001)),
                  fabs(get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{-20, 80}, {20, 80}}, 80.001)),
                  get_intersection(
                    YAC_GREAT_CIRCLE_EDGE,
                    (double[2][2]){{0, 85}, {20, 80}}, 80.001),},
       (double[]){85,80.001,80.001,80.001,80.001},
       (enum yac_edge_type[]){YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
                              YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE,
                              YAC_GREAT_CIRCLE_EDGE}, 5)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){10,20,15},
       (double[]){82,82,88}, gc_edges, 3),
    (double[2]){90, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){10,20,15},
       (double[]){82,82,88}, gc_edges, 3)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){10,20,15},
       (double[]){-82,-82,-88}, gc_edges, 3),
    (double[2]){90, 80}, NULL, 0);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){10,20,15},
       (double[]){-82,-82,-88}, gc_edges, 3),
    (double[2]){-90, -80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){10,20,15},
       (double[]){-82,-82,-88}, gc_edges, 3)}, 1);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){10,20,15},
       (double[]){82,82,88}, gc_edges, 3),
    (double[2]){-90, -80}, NULL, 0);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){10,20,15},
       (double[]){82,82,88}, gc_edges, 3),
    (double[2]){90, 90}, NULL, 0);

  // a triangle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){10,20,15},
       (double[]){82,82,88}, gc_edges, 3),
    (double[2]){-90, -90}, NULL, 0);

  // a triangle cell covering the pole
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){80,80,80}, gc_edges, 3),
    (double[2]){90, 70},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){0,120,240},
       (double[]){80,80,80}, gc_edges, 3)}, 1);

  // a triangle cell covering the pole
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){80,80,80}, gc_edges, 3),
    (double[2]){60, 70}, NULL, 0);

  // a triangle cell covering the pole
  temp_lon[0] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {120, 70}}, 80);
  temp_lon[1] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {240, 70}}, 80);
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){85,70,70}, gc_edges, 3),
    (double[2]){90, 80},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){0,temp_lon[0],temp_lon[1]},
       (double[]){85,80,80},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE}, 3)}, 1);

  // a triangle cell covering the pole
  temp_lon[0] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {120, 70}}, 72);
  temp_lon[1] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{120, 70}, {240, 70}}, 72);
  temp_lon[1] = (temp_lon[1] < 0)?(temp_lon[1] + 360):(temp_lon[1]);
  temp_lon[1] = (temp_lon[1] < 180)?(temp_lon[1]):(360 - temp_lon[1]);
  temp_lon[2] = 360 - temp_lon[1];
  temp_lon[3] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {240, 70}}, 72);
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){85,70,70}, gc_edges, 3),
    (double[2]){90, 72},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){0,temp_lon[0],temp_lon[1],temp_lon[2],temp_lon[3]},
       (double[]){85,72,72,72,72},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
         YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE}, 5)}, 1);

  // a triangle cell covering the pole
  temp_lon[0] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {120, 70}}, 80);
  temp_lon[1] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {240, 70}}, 80);
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){85,70,70}, gc_edges, 3),
    (double[2]){80, 65},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){temp_lon[0],120,240,temp_lon[1]},
       (double[]){80,70,70,80},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
         YAC_LAT_CIRCLE_EDGE}, 4)}, 1);

  // a triangle cell covering the pole
  temp_lon[0] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {120, 70}}, 84);
  temp_lon[1] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {120, 70}}, 80);
  temp_lon[2] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {240, 70}}, 80);
  temp_lon[3] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {240, 70}}, 84);
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){85,70,70}, gc_edges, 3),
    (double[2]){80, 84},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){temp_lon[0],temp_lon[1],temp_lon[2],temp_lon[3]},
       (double[]){84,80,80,84},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
         YAC_LAT_CIRCLE_EDGE}, 4)}, 1);

  // a triangle cell covering the pole
  temp_lon[0] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {120, 70}}, 80);
  temp_lon[1] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {120, 70}}, 72);
  temp_lon[2] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{240, 70}, {120, 70}}, 72);
  temp_lon[2] = (temp_lon[2] < 0)?(temp_lon[2] + 360):(temp_lon[2]);
  temp_lon[2] = (temp_lon[2] < 180)?(temp_lon[2]):(360 - temp_lon[2]);
  temp_lon[3] = 360 - temp_lon[2];
  temp_lon[4] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {240, 70}}, 72);
  temp_lon[5] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 85}, {240, 70}}, 80);
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){85,70,70}, gc_edges, 3),
    (double[2]){80, 72},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){temp_lon[0],temp_lon[1],temp_lon[2],
                  temp_lon[3],temp_lon[4],temp_lon[5]},
       (double[]){80,72,72,72,72,80},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
         YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE}, 6)}, 1);

  // a triangle cell touching the pole
  temp_lon[0] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 90}, {120, 80}}, 85);
  temp_lon[1] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 90}, {240, 80}}, 85);
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,120,240},
       (double[]){90,80,80}, gc_edges, 3),
    (double[2]){90, 85},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){0,temp_lon[0],temp_lon[1]},
       (double[]){90,85,85},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE}, 3)}, 1);

  // a triangle cell touching the pole
  temp_lon[0] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{0, 0}, {10, 5}}, 2.5);
  temp_lon[1] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{10, -5}, {10, 5}}, 2.5);
  temp_lon[2] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{10, -5}, {10, 5}}, -2.5);
  temp_lon[3] =
    get_intersection(YAC_GREAT_CIRCLE_EDGE, (double[2][2]){{10, -5}, {0, 0}}, -2.5);
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,10,10},
       (double[]){0,5,-5}, gc_edges, 3),
    (double[2]){2.5, -2.5},
    (struct yac_grid_cell[]){
      generate_cell_deg(
       (double[]){0,temp_lon[0],temp_lon[1],temp_lon[2],temp_lon[3]},
       (double[]){0,2.5,2.5,-2.5,-2.5},
       (enum yac_edge_type[]){
         YAC_GREAT_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
         YAC_LAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE}, 5)}, 1);

  // empty great circle cell
  test_lat_clipping(
    generate_cell_deg(
       (double[]){0,0,0},
       (double[]){-5,0,5}, gc_edges, 3),
    (double[2]){2.5, -2.5}, NULL, 0);

  return TEST_EXIT_CODE;
}

static int compare_cells(struct yac_grid_cell a, struct yac_grid_cell b) {

  // Trigonometric functions can be very inaccurate depending on the compiler
  // and the optimsation level. Therefore, we use angle_tol instead of
  // yac_angle_tol in this routine.
  double const angle_tol = 1e-6;

  if ((yac_huiliers_area(a) < angle_tol * angle_tol) &&
      (yac_huiliers_area(b) < angle_tol * angle_tol)) return 0;

  if (a.num_corners != b.num_corners) return 1;

  if (a.num_corners == 0) return 0;

  for (int order = -1; order <= 1; order += 2) {
    for (int start = 0; start < (int)(a.num_corners); ++start) {

      int differences = 0;

      for (int i = 0; i < (int)(a.num_corners); ++i) {

        int j =
          ((int)(a.num_corners) + start + order * i) % (int)(a.num_corners);

        if (get_vector_angle(
              a.coordinates_xyz[i], b.coordinates_xyz[j]) > angle_tol)
          ++differences;

        else if ((a.edge_type[i] == YAC_LAT_CIRCLE_EDGE) !=
                 (b.edge_type[
                    (((int)(a.num_corners))+j-(order<0))%
                    ((int)(a.num_corners))] == YAC_LAT_CIRCLE_EDGE))
          ++differences;
      }

      if (!differences) return 0;
    }
  }

  return 1;
}

static void test_lat_clipping(struct yac_grid_cell cell,
                              double deg_lat_bounds[2],
                              struct yac_grid_cell * ref_cells,
                              unsigned num_ref_cells) {

  double const area_tol = 1e-8 * yac_huiliers_area(cell);

  double rad_lat_bounds[2] =
    {deg_lat_bounds[0] * YAC_RAD, deg_lat_bounds[1] * YAC_RAD};

  double mem_dummy_[128][3];
  enum yac_edge_type edge_dummy[128];
  struct yac_grid_cell overlap_cell,
                   test_cell =
        (struct yac_grid_cell){.coordinates_xyz = mem_dummy_,
                           .edge_type = edge_dummy,
                           .num_corners = cell.num_corners};

  yac_init_grid_cell(&overlap_cell);

  for (int pole = 0; pole < 2; ++pole) {
    for (int order = -1; order <= 1; order += 2) {
      for (int start = 0; start < (int)(cell.num_corners); ++start) {
        for (int i = 0; i < ((int)(cell.num_corners)); ++i) {
          int j = (((int)(cell.num_corners))+start+i*order)%
                  ((int)(cell.num_corners));
          test_cell.coordinates_xyz[i][0] = cell.coordinates_xyz[j][0];
          test_cell.coordinates_xyz[i][1] = cell.coordinates_xyz[j][1];
          test_cell.coordinates_xyz[i][2] = cell.coordinates_xyz[j][2];
          j = (((int)(cell.num_corners)) + j - (order < 0))%
              ((int)(cell.num_corners));
          test_cell.edge_type[i] = cell.edge_type[j];
        }

        for (int bound_tol_a = -1; bound_tol_a <= 1; ++bound_tol_a) {
          for (int bound_tol_b = -1; bound_tol_b <= 1; ++bound_tol_b) {
            for (int bound_order = 0; bound_order <= 1; ++bound_order) {
              double test_rad_lat_bounds[2] = {
                MIN(M_PI_2,
                  MAX(-M_PI_2,
                      rad_lat_bounds[bound_order] +
                      (double)bound_tol_a * yac_angle_low_tol)),
                MIN(M_PI_2,
                  MAX(-M_PI_2,
                    rad_lat_bounds[bound_order^1] +
                    (double)bound_tol_b * yac_angle_low_tol))};
              yac_cell_lat_clipping(
                1, &test_cell,test_rad_lat_bounds, &overlap_cell);
              if (num_ref_cells == 0) {
                if ((overlap_cell.num_corners != 0) &&
                    (yac_huiliers_area(overlap_cell) > area_tol))
                  PUT_ERR("ERROR: wrong clipping cell\n");
              } else {
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
    }

    // switch to the other pole
    for (int i = 0; i < (int)(cell.num_corners); ++i)
      cell.coordinates_xyz[i][2] *= -1.0;
    rad_lat_bounds[0] *= -1.0, rad_lat_bounds[1] *= -1.0;
    for (unsigned i = 0; i < num_ref_cells; ++i)
      for (int j = 0; j < (int)(ref_cells[i].num_corners); ++j)
        ref_cells[i].coordinates_xyz[j][2] *= -1.0;
  }

  yac_free_grid_cell(&cell);
  for (unsigned i = 0; i < num_ref_cells; ++i)
    yac_free_grid_cell(&ref_cells[i]);
  yac_free_grid_cell(&overlap_cell);
}
