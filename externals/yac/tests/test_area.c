/**
 * @file test_area.c
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
#include "clipping.h"
#include "area.h"
#include "tests.h"
#include "geometry.h"
#include "test_common.h"

//#define VERBOSE

static int compare_area(double a, double b);
static void test_area(struct grid_cell cell, double ref_area,
                      double (*area_func)(struct grid_cell));

int main (void) {

  enum yac_edge_type gc_edges[5] =
    {GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
     GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE};
  enum yac_edge_type lat_edges[4] =
    {LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE};

  double area;

#ifdef VERBOSE
  puts ("----- calculating area -----\n");
#endif

  {
    /* set the test data over the Equator

               0.0 (lon)
               0.5 (lat)
               / \
              /   \
    -0.5     /     \  0.5
     0.0     \     /  0.0
              \   /
               \ /
               0.0
              -0.5         
    */

    struct grid_cell Quad =
      generate_cell_deg((double[4]){0.5, 0.0, -0.5, 0.0},
                        (double[4]){0.0, 0.5, 0.0, -0.5}, gc_edges, 4);

    area = yac_pole_area ( Quad );
#ifdef VERBOSE
    printf ( "pole     area of quad   over Equator    is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Quad );
#ifdef VERBOSE
    printf ( "huiliers area of quad   over Equator    is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif
    yac_free_grid_cell(&Quad);
  }

  {
  /* Triangle, half of the above Quad */

    struct grid_cell Tri =
      generate_cell_deg((double[3]){0.5, 0.0, -0.5},
                        (double[3]){0.0, 0.5, 0.0}, gc_edges, 3);

    area = yac_triangle_area ( Tri );
#ifdef VERBOSE
    printf ( "ICON     area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_pole_area ( Tri );
#ifdef VERBOSE
    printf ( "pole     area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Tri );
#ifdef VERBOSE
    printf ( "huiliers area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif
    yac_free_grid_cell(&Tri);
  }

  {
    /* set the test data over the Equator

    -0.5               0.5 (lon)
    0.5  -----------  0.5 (lat)
         |         |
         |         |
         |         |
         |         |
         |         |
    -0.5  -----------  0.5
    -0.5              -0.5

    */

    struct grid_cell Quad =
      generate_cell_deg((double[4]){-0.5, 0.5, 0.5, -0.5},
                        (double[4]){-0.5, -0.5, 0.5, 0.5}, gc_edges, 4);

    area = yac_pole_area ( Quad );
#ifdef VERBOSE
    printf ( "pole     area of quad   over Equator    is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Quad );
#ifdef VERBOSE
    printf ( "huiliers area of quad   over Equator    is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    double ref_area = yac_pole_area(Quad);

    test_area(Quad, ref_area, yac_pole_area);
    test_area(Quad, ref_area, yac_huiliers_area);

    yac_free_grid_cell(&Quad);
  }

  {
    /* Half of the above Quad */

    struct grid_cell Tri =
      generate_cell_deg((double[3]){-0.5, 0.5, 0.5},
                        (double[3]){-0.5, -0.5, 0.5}, gc_edges, 3);

    area = yac_triangle_area ( Tri );
#ifdef VERBOSE
    printf ( "ICON     area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_pole_area ( Tri );
#ifdef VERBOSE
    printf ( "pole     area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Tri );
#ifdef VERBOSE
    printf ( "huiliers area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    yac_free_grid_cell(&Tri);
  }

  {

    struct grid_cell Tri =
      generate_cell_deg((double[3]){150.0, 180.0, 180.0},
                        (double[3]){-60.0, -60.0, -90.0}, gc_edges, 3);

    area = yac_triangle_area ( Tri );
#ifdef VERBOSE
    printf ( "ICON     area of triangle over South Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_pole_area ( Tri );
#ifdef VERBOSE
    printf ( "pole     area of triangle over South Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Tri );
#ifdef VERBOSE
    printf ( "huiliers area of triangle over South Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    yac_free_grid_cell(&Tri);
  }

  {
    /* set the test data over the North Pole

    90.0                0.0 (lon)
    89.5  -----------  89.5 (lat)
          |         |
          |         |
          |    x    |
          |         |
          |         |
    180.0  ----------- 270.0
    89.5               89.5

    */

    struct grid_cell Quad =
      generate_cell_deg((double[4]){0.0, 90.0, -180.0, -90.0},
                        (double[4]){89.5, 89.5, 89.5, 89.5}, gc_edges, 4);

    area = yac_pole_area ( Quad );
#ifdef VERBOSE
    printf ( "pole   area of quad     over North Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Quad );
#ifdef VERBOSE
    printf ( "huiliers area of quad   over North Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    yac_free_grid_cell(&Quad);
  }

  {
    struct grid_cell Quad =
      generate_cell_deg((double[4]){-180.0, -150.0, -150.0, -180.0},
                        (double[4]){60.0, 60.0, 90.0, 90.0}, gc_edges, 4);

    area = yac_pole_area ( Quad );
#ifdef VERBOSE
    printf ( "pole     area of quad   near North Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Quad );
#ifdef VERBOSE
    printf ( "huiliers area of quad   near North Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    yac_free_grid_cell(&Quad);
  }

  {
    struct grid_cell Quad =
      generate_cell_deg((double[4]){-180.0, -150.0, -150.0, -180.0},
                        (double[4]){-90.0, -90.0, -60.0, -60.0}, gc_edges, 4);

#ifdef VERBOSE
    printf ( "ICON     area of quad   near South Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_pole_area ( Quad );
#ifdef VERBOSE
    printf ( "pole     area of quad   near South Pole is %f sqr Km \n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Quad );
#ifdef VERBOSE
    printf ( "huiliers area of quad   near South Pole is %f sqr Km \n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    yac_free_grid_cell(&Quad);
  }

  {
    struct grid_cell Quin =
      generate_cell_deg((double[5]){-180.0, -150.0, 0.0, -150.0, -180.0},
                        (double[5]){60.0, 60.0, 90.0, 90.0, 90.0}, gc_edges, 5);

    area = yac_pole_area ( Quin );
#ifdef VERBOSE
    printf ( "pole     area of quin   near North Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Quin );
#ifdef VERBOSE
    printf ( "huiliers area of quin   near North Pole is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    yac_free_grid_cell(&Quin);
  }

  {
    /* Triangle that contains an edge of length zero */

    struct grid_cell Tri =
      generate_cell_deg((double[3]){-0.5,  0.5, -0.5},
                        (double[3]){-0.5, -0.5, -0.5}, gc_edges, 3);

    area = yac_triangle_area ( Tri );
#ifdef VERBOSE
    printf ( "ICON     area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_pole_area ( Tri );
#ifdef VERBOSE
    printf ( "pole     area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area ( Tri );
#ifdef VERBOSE
    printf ( "huiliers area of triangle over Equator  is %f sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    if ((area < 0.0) || (area > 1e-9))
      PUT_ERR("ERROR in yac_huiliers_area for zero size triangle");

    yac_free_grid_cell(&Tri);
  }

  {
    struct grid_cell Quad =
      generate_cell_deg((double[4]){0, 1, 1, 0},
                        (double[4]){0, 0, 1, 1}, gc_edges, 4);
    struct grid_cell Tri =
      generate_cell_deg((double[3]){1, 1, 1},
                        (double[3]){0, 0, 1}, gc_edges, 3);

    double base_area = yac_huiliers_area(Quad);

    yac_free_grid_cell(&Quad);
    yac_free_grid_cell(&Tri);

    /* A dx of 1e-8 deg corresponds to a dx of ~1 mm on the Earth surface */

    double dx[] = {1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,
                   0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    for (size_t i = 0; i < sizeof(dx) / sizeof(dx[0]); ++i) {


      struct grid_cell Quad =
        generate_cell_deg((double[4]){0, 1.0-dx[i], 1, 0},
                          (double[4]){0, 0, 1, 1}, gc_edges, 4);
      struct grid_cell Tri =
        generate_cell_deg((double[3]){1.0-dx[i], 1, 1},
                          (double[3]){0, 0, 1}, gc_edges, 3);

      double partial_area = yac_huiliers_area(Quad);
      partial_area += yac_huiliers_area(Tri);

      yac_free_grid_cell(&Quad);
      yac_free_grid_cell(&Tri);

      if (compare_area(base_area, partial_area))
        PUT_ERR("yac_huiliers_area computed wrong area\n");

#ifdef VERBOSE
      printf ("accuracy check for (base_area - partial_area) %8.1e sqr km for dx %.1e\n",
              base_area - partial_area, dx[i]);
#endif
    }
  }

  {
    struct grid_cell Tri =
      generate_cell_deg(
        (double[3]){0.09817477/YAC_RAD, 0.09817477/YAC_RAD, 0.098174496/YAC_RAD},
        (double[3]){0.353325247/YAC_RAD, 0.353335501/YAC_RAD, 0.353335609/YAC_RAD},
        gc_edges, 3);

    area = yac_pole_area(Tri);

    if (area < 0)
       PUT_ERR("pole_area computed wrong area\n");

#ifdef VERBOSE
    printf ("pole   accuracy check of very small area %.2e sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    area = yac_huiliers_area(Tri);

    if (area < 0)
       PUT_ERR("huiliers_area computed wrong area\n");

#ifdef VERBOSE
    printf ("huiliers accuracy check of very small area %.2e sqr km\n", area * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS );
#endif

    yac_free_grid_cell(&Tri);
  }

  {
    double edge_length_deg = 1.0e-8; // ~1.11 mm

#ifdef VERBOSE
     printf ("\n accuracy check for a unit square\n");
     printf (" edge length [m] |        plain       |       simple       |        pole        |      huiliers\n");
     printf (" ----------------+--------------------+--------------------+--------------------+--------------------\n");
#endif
     for ( unsigned i = 1; i < 11; ++i ) {

        double edge_length_m =
          edge_length_deg * YAC_RAD * YAC_EARTH_RADIUS * 1.0e3;

        struct grid_cell Quad =
          generate_cell_deg((double[4]){0, edge_length_deg, edge_length_deg, 0},
                            (double[4]){0, 0, edge_length_deg, edge_length_deg},
                            gc_edges, 4);

        double t_area = edge_length_m                * edge_length_m;
        double p_area = yac_planar_3dcell_area(Quad) * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS * 1.0e6;
        double s_area = yac_pole_area(Quad)          * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS * 1.0e6;
        double h_area = yac_huiliers_area(Quad)      * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS * 1.0e6;

#ifndef VERBOSE
        if (0)
#endif
          printf (" %.8e m|%.8e sqr m|%.8e sqr m|%.8e sqr m|%.8e sqr m\n",
                  edge_length_m, p_area, t_area, s_area, h_area);

        edge_length_deg *= 10.0;

        yac_free_grid_cell(&Quad);
     }

  }

  { // cell with two edges (one great circle and one circle of latitude)

    for (int y = -90; y <= 90; y += 5) {

      struct grid_cell cell =
        generate_cell_deg((double[2]){0.0, 20.0},
                          (double[2]){(double)y, (double)y},
                          (enum yac_edge_type[]){
                            GREAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 2);

      // compute reference area
      double ref_area;

      {
        double coordinates_x[3] = {0.0, 20.0, 0.0};
        double coordinates_y[3] = {(double)y, (double)y, 0};

        if ((y == -90) || (y == 90) || (y == 0)) {

          ref_area = 0.0;

        } else {

          if (y < 0) {
            ref_area = (coordinates_x[1]*YAC_RAD - coordinates_x[0]*YAC_RAD) *
                        (1 - cos(M_PI_2 + coordinates_y[0]*YAC_RAD));
            coordinates_y[2] = -90.0;
          } else {
            ref_area = (coordinates_x[1]*YAC_RAD - coordinates_x[0]*YAC_RAD) *
                        (1 - cos(M_PI_2 - coordinates_y[0]*YAC_RAD));
            coordinates_y[2] = 90.0;
          }

          struct grid_cell tmp_cell =
            generate_cell_deg(coordinates_x, coordinates_y, gc_edges, 3);
          // ref_area *= YAC_EARTH_RADIUS * YAC_EARTH_RADIUS;
          ref_area -= yac_huiliers_area(tmp_cell);

          yac_free_grid_cell(&tmp_cell);
        }
      }

      test_area(cell, ref_area, yac_pole_area);
      test_area(cell, ref_area, yac_huiliers_area);

      yac_free_grid_cell(&cell);
    }
  }

  { // cell with lon lat edges and one edge on the equator

    double ref_area = sin(0.5 * YAC_RAD) * (0.5 * YAC_RAD);
    // ref_area *= YAC_EARTH_RADIUS * YAC_EARTH_RADIUS;

    double coordinates_x[2][4] = {{-0.5, 0.0, 0.0, -0.5},
                                  { 0.0, 0.5, 0.5,  0.0}};
    double coordinates_y[2][4] = {{-0.5, -0.5, 0.0, 0.0},
                                  { 0.0,  0.0, 0.5, 0.5}};

    for (int x = 0; x < 2; ++x) {
      for (int y = 0; y < 2; ++y) {
        struct grid_cell cell =
          generate_cell_deg(coordinates_x[x], coordinates_y[y], lat_edges, 4);

        test_area(cell, ref_area, yac_pole_area);
        test_area(cell, ref_area, yac_huiliers_area);

        yac_free_grid_cell(&cell);
      }
    }
  }

  { // cell with lon lat edges (both lat edges on the different side of the
    // equator)
    double ref_area = (2.0 * sin(0.5 * YAC_RAD)) * (1.0 * YAC_RAD);
    // ref_area *= YAC_EARTH_RADIUS * YAC_EARTH_RADIUS;

    struct grid_cell cell =
      generate_cell_deg((double[4]){-0.5, 0.5, 0.5, -0.5},
                        (double[4]){-0.5, -0.5, 0.5, 0.5}, lat_edges, 4);

    test_area(cell, ref_area, yac_pole_area);
    test_area(cell, ref_area, yac_huiliers_area);

    yac_free_grid_cell(&cell);
  }

  { // cell with lon lat edges (both lat edges on the same side of the equator)
    double ref_area = sin(1.0 * YAC_RAD);
    ref_area -= sin(0.5 * YAC_RAD);
    ref_area *= 1.0 * YAC_RAD;
    // ref_area *= YAC_EARTH_RADIUS * YAC_EARTH_RADIUS;
    
    for (int i = 0; i < 2; ++i) {

      double coordinates_x[4] = {-0.5, 0.5, 0.5, -0.5};
      double coordinates_y[2][4] = {{ 0.5, 0.5, 1.0,  1.0},
                                    {-0.5,-0.5,-1.0, -1.0}};
      enum yac_edge_type edge_types[4] =
        {LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE};

      struct grid_cell cell =
        generate_cell_deg(coordinates_x, coordinates_y[i], edge_types, 4);

      test_area(cell, ref_area, yac_pole_area);
      test_area(cell, ref_area, yac_huiliers_area);

      yac_free_grid_cell(&cell);
    }
  }

  { // concave cell

    double intersections[2][3];

    if (!intersect(GREAT_CIRCLE_EDGE, 20.0, 59.0, 0.0, 80.0,
                   LAT_CIRCLE_EDGE, 0.0, 60.0, 20.0, 60.0,
                   intersections[0]))
      return EXIT_FAILURE;
    if (!intersect(GREAT_CIRCLE_EDGE, 20.0, 59.0, -20.0, 59.0,
                   LAT_CIRCLE_EDGE, 0.0, 60.0, 20.0, 60.0,
                   intersections[1]))
      return EXIT_FAILURE;

    double intersections_lon[2], intersections_lat[2];
    XYZtoLL(intersections[0], &(intersections_lon[0]), &(intersections_lat[0]));
    XYZtoLL(intersections[1], &(intersections_lon[1]), &(intersections_lat[1]));

    // this cell is convex
    struct grid_cell partial_ref_cell =
      generate_cell_deg(
        (double[]){20.0, intersections_lon[0]/YAC_RAD,
                   intersections_lon[1]/YAC_RAD},
        (double[]){59.0, 60.0, 60.0},
        (enum yac_edge_type[]){GREAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}, 3);

    double partial_area = yac_pole_area(partial_ref_cell);

    if (compare_area(partial_area,
                     yac_huiliers_area(partial_ref_cell)))
      PUT_ERR("pole_area != huiliers_area for convex cell");

    double ref_area = 2.0 * partial_area;

    struct grid_cell cells[2] =
      {generate_cell_deg(
         (double[]){20.0,
                     intersections_lon[0]/YAC_RAD,
                     intersections_lon[1]/YAC_RAD,
                    -intersections_lon[1]/YAC_RAD,
                    -intersections_lon[0]/YAC_RAD, -20.0},
         (double[]){59.0, 60.0, 60.0, 60.0, 60.0, 59.0},
         (enum yac_edge_type[]){
           GREAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
           LAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}, 6),
       generate_cell_deg(
         (double[]){20.0,
                     intersections_lon[0]/YAC_RAD,
                    -intersections_lon[0]/YAC_RAD,
                    -20.0,
                    -intersections_lon[1]/YAC_RAD,
                     intersections_lon[1]/YAC_RAD},
         (double[]){59.0, 60.0, 60.0, 59.0, 60.0, 60.0},
         (enum yac_edge_type[]){
           GREAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
           GREAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}, 6)};

    for (int i = 0; i < 2; ++i) {
      test_area(cells[i], ref_area, yac_pole_area);
      test_area(cells[i], ref_area, yac_huiliers_area);
    }
    yac_free_grid_cell(&(cells[0]));
    yac_free_grid_cell(&(cells[1]));
    yac_free_grid_cell(&partial_ref_cell);
  }

  { // case extracted from ICON example, which requires a high accuracy in the
    // area computation
    struct grid_cell cells[] = {
      generate_cell_deg((double[]){-151.1489316535362661,
                                   -151.7087167856307133,
                                   -150.6174180645230933},
                        (double[]){  10.6333746562670566,
                                     11.5836995038666437,
                                     11.6243355661022782}, gc_edges, 3),
      generate_cell_deg((double[]){-152.8132296086113513,
                                   -153.8786613849353557,
                                   -151.6929566531634066},
                        (double[]){  13.4512385235805798,
                                     11.4650630285209161,
                                     11.5568484830228009}, gc_edges, 3),
      generate_cell_deg((double[]){-152.8132296086113513,
                                   -151.6929566531634066,
                                   -150.6265735326831532},
                        (double[]){  13.4512385235805798,
                                     11.5568484830228009,
                                     13.5446528932721275}, gc_edges, 3),
      generate_cell_deg((double[]){-153.8786613849353557,
                                   -151.6929566531634066,
                                   -152.7627337461449031},
                        (double[]){  11.4650630285209161,
                                     11.5568484830228009,
                                      9.5827293336264265}, gc_edges, 3),
      generate_cell_deg((double[]){-151.6929566531634066,
                                   -152.7627337461449031,
                                   -150.5840076321059371},
                        (double[]){  11.5568484830228009,
                                      9.5827293336264265,
                                      9.6520960946755832}, gc_edges, 3),
      generate_cell_deg((double[]){-151.6929566531634066,
                                   -150.6265735326831532,
                                   -149.5139953075027108},
                        (double[]){  11.5568484830228009,
                                     13.5446528932721275,
                                     11.6284420125700851}, gc_edges, 3),
      generate_cell_deg((double[]){-151.6929566531634066,
                                   -149.5139953075027108,
                                   -150.5840076321059371},
                        (double[]){  11.5568484830228009,
                                     11.6284420125700851,
                                      9.6520960946755832}, gc_edges, 3)};
    size_t num_cells = sizeof(cells) / sizeof(cells[0]);

    struct grid_cell overlap_cells[sizeof(cells) / sizeof(cells[0]) - 1];
    for (size_t i = 0; i < num_cells-1; ++i)
      yac_init_grid_cell(&overlap_cells[i]);

    yac_cell_clipping(num_cells - 1, &cells[1], cells[0], &overlap_cells[0]);

    double tgt_area = yac_huiliers_area(cells[0]);
    double area_tol = tgt_area * 1e-8;
    double diff_area = tgt_area;

    for (size_t i = 0; i < num_cells-1; ++i)
      diff_area -= yac_huiliers_area(overlap_cells[i]);

    if (fabs(diff_area) > area_tol)
      PUT_ERR("yac_huiliers_area is too inaccurate");

    for (size_t i = 0; i < num_cells; ++i)
      yac_free_grid_cell(&cells[i]);
    for (size_t i = 0; i < num_cells-1; ++i)
      yac_free_grid_cell(&overlap_cells[i]);
  }

  { // check for a bug found by Uwe Schulzweida

    struct grid_cell Tri =
      generate_cell_deg((double[3]){0.5, 0.5, -0.5},
                        (double[3]){0.0, 0.0, 0.0}, gc_edges, 3);

    if (yac_triangle_area(Tri) != 0.0)
      PUT_ERR("ERROR in yac_triangle_area");
    yac_free_grid_cell(&Tri);

    Tri = generate_cell_deg((double[3]){0.5, -0.5, 0.5},
                            (double[3]){0.0, 0.0, 0.0}, gc_edges, 3);

    if (yac_triangle_area(Tri) != 0.0)
      PUT_ERR("ERROR in yac_triangle_area");
    yac_free_grid_cell(&Tri);

    Tri = generate_cell_deg((double[3]){-0.5, 0.5, 0.5},
                            (double[3]){0.0, 0.0, 0.0}, gc_edges, 3);

    if (yac_triangle_area(Tri) != 0.0)
      PUT_ERR("ERROR in yac_triangle_area");
    yac_free_grid_cell(&Tri);
  }

  return TEST_EXIT_CODE;
}

static void test_area(struct grid_cell cell, double ref_area,
                      double (*area_func)(struct grid_cell)) {

  double mem_dummy_[128][3];
  enum yac_edge_type edge_dummy[128];
  struct grid_cell test_cell = {.coordinates_xyz = mem_dummy_,
                                .edge_type = edge_dummy,
                                .num_corners = cell.num_corners};

  for (int order = -1; order <= 1; order += 2) {
    for (int start = 0; start < (int)(cell.num_corners); ++start) {

      for (int i = 0; i < (int)(cell.num_corners); ++i) {
        int j = ((int)(cell.num_corners)+start+i*order)%(int)(cell.num_corners);
        test_cell.coordinates_xyz[i][0] = cell.coordinates_xyz[j][0];
        test_cell.coordinates_xyz[i][1] = cell.coordinates_xyz[j][1];
        test_cell.coordinates_xyz[i][2] = cell.coordinates_xyz[j][2];
        j = ((int)(cell.num_corners) + j - (order < 0))%(int)(cell.num_corners);
        test_cell.edge_type[i] = cell.edge_type[j];
      }

      double area = area_func(test_cell);

      if (compare_area(area, ref_area)) PUT_ERR("ERROR: wrong area\n");
    }
  }
}

static int compare_area(double a, double b) {

  double diff = fabs(a - b);
  double tol = MAX(MIN(a,b)*1e-6, YAC_AREA_TOL);

  if (diff < tol) return 0;
  else return (a > b) - (a < b);
}
