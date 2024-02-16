/**
 * @file test_compute_overlap_area.c
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

#include "tests.h"
#include "clipping.h"
#include "geometry.h"
#include "area.h"
#include "test_common.h"
#include "vtk_output.h"

static enum yac_edge_type gc_edges[] =
  {GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
   GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE};

static void check_overlap(
  double * lon_a, double * lat_a, int count_a,
  double * lon_b, double * lat_b, int count_b, double ref_area,
  double const * ref_barycenter);
static void check_overlap_polygons(
  struct grid_cell * Polygons, double ref_area,
  double const * ref_barycenter);

//#define VERBOSE

int main (void) {

  {
    /* intersection between concave pentagon with a triangle

           10.0 (lon)   11.0 (lon)
           32.0 (lat)   32.0 (lat)
                  x--------x
                 /          \
                /    10.5    \
               /     31.0     \
   9.5 (lon)  /       x        \   11.5
  30.0 (lat) x                  x  30.0


   with


                    10.5 (lon)
                    31.5 (lat)
                      x
                     / \
                    /   \
           10.0    /     \   11.0
           29.5   x-------x  29.5
    */

    // the reference area any barycenter
    double ref_area;
    double ref_barycenter[3] = {0, 0, 0};
    {
      double intersection[2][3], intersection_lon[2], intersection_lat[2];

      if (!intersect(GREAT_CIRCLE_EDGE, 11.5, 30.0, 10.5, 31.0,
                     GREAT_CIRCLE_EDGE, 11.0, 29.5, 10.5, 31.5,
                     intersection[0]) ||
          !intersect(GREAT_CIRCLE_EDGE, 9.5, 30.0, 10.5, 31.0,
                     GREAT_CIRCLE_EDGE, 10.0, 29.5, 10.5, 31.5,
                     intersection[1]))
        return EXIT_FAILURE;
      for (int i = 0; i < 2; ++i) {
        XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
        intersection_lon[i] /= YAC_RAD;
        intersection_lat[i] /= YAC_RAD;
      }

      struct grid_cell Polygon =
        generate_cell_deg(
          (double[]){10.5, intersection_lon[0], 10.5, intersection_lon[1]},
          (double[]){31.0, intersection_lat[0], 31.5, intersection_lat[1]},
          gc_edges, 4);
      ref_area = yac_huiliers_area_info(Polygon, ref_barycenter, 1.0);
      normalise_vector(ref_barycenter);
      yac_free_grid_cell(&Polygon);
    }

    check_overlap(
      (double[]){10.0, 11.0, 11.5, 10.5,  9.5},
      (double[]){32.0, 32.0, 30.0, 31.0, 30.0}, 5,
      (double[]){10.5, 11.0, 10.0},
      (double[]){31.5, 29.5, 29.5}, 3,
      ref_area, ref_barycenter);
  }

  {
    /* intersection between square and a triangle

                -5.0      5.0
                  x--------x  5.0
                  |        |
                  |        |
                  |        |
                  x--------x -5.0

   with

                -4.0 0.0 4.0
                      x      6.0
                     / \
                    /   \
                   /     \
                  x-------x -4.0
    */
    // the reference area any barycenter
    double ref_area;
    double ref_barycenter[3] = {0, 0, 0};
    {
      double intersection[2][3], intersection_lon[2], intersection_lat[2];

      if (!intersect(GREAT_CIRCLE_EDGE, 0, 6, 4, -4,
                     GREAT_CIRCLE_EDGE, -5, 5, 5, 5,
                     intersection[0]) ||
          !intersect(GREAT_CIRCLE_EDGE, 0, 6, -4, -4,
                     GREAT_CIRCLE_EDGE, -5, 5, 5, 5,
                     intersection[1]))
        return EXIT_FAILURE;
      for (int i = 0; i < 2; ++i) {
        XYZtoLL(intersection[i], &intersection_lon[i], &intersection_lat[i]);
        intersection_lon[i] /= YAC_RAD;
        intersection_lat[i] /= YAC_RAD;
      }

      struct grid_cell Polygon =
        generate_cell_deg(
          (double[]){-4, 4, intersection_lon[0], intersection_lon[1]},
          (double[]){-4, -4, intersection_lat[0], intersection_lat[1]}, gc_edges, 4);
      ref_area = yac_huiliers_area_info(Polygon, ref_barycenter, 1.0);
      normalise_vector(ref_barycenter);
      yac_free_grid_cell(&Polygon);
    }

    check_overlap(
      (double[]){-5,  5, 5, -5}, (double[]){-5, -5, 5,  5}, 4,
      (double[]){-4, 4, 0}, (double[]){-4, -4, 6}, 3,
      ref_area, ref_barycenter);
  }

  {
    /* set the test cells (the target node is one of the corners of the
                           concave polygon; the center of both shapes
                           is the north pole)


                   +------+
                   |      |
                   |  /\  |
                   | /  \ |
                   +/    \+

   with

                   +------+
                   |      |
                   |      |
                   |      |
                   +------+

    */

    // the reference area and barycenter
    double ref_area = 0.0;
    double ref_barycenter[] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
      double tri_lon[3][3] = {{0,0,90},{90,90,180},{180,180,270}};
      double tri_lat[3] = {90,89,89};
      struct grid_cell Tri =
        generate_cell_deg(tri_lon[i], tri_lat, gc_edges, 3);
      ref_area += yac_huiliers_area_info(Tri, ref_barycenter, 1.0);
      yac_free_grid_cell(&Tri);
    }
    normalise_vector(ref_barycenter);

    check_overlap(
      (double[]){0,0,90,180,270}, (double[]){90,88,88,88,88}, 5,
      (double[]){0,90,180,270}, (double[]){89,89,89,89}, 4,
      ref_area, ref_barycenter);
  }

  {
    /*

      +---+        +---+
      |   |        |   |
      |   |  +--+  |   |
      |   |  |  |  |   |
      |   |  +--+  |   |
      |   |        |   |
      |   +--------+   |
      |                |
      +----------------+

    */

    // the reference area and barycenter
    double ref_area = 0.0;
    double * ref_barycenter = NULL;

    check_overlap(
      (double[]){2,3,3,-3,-3,-2,-2,2}, (double[]){3,3,-3,-3,3,3,-2,-2}, 8,
      (double[]){1,1,-1,-1}, (double[]){1,-1,-1,1}, 4,
      ref_area, ref_barycenter);
  }

  {
    /*

    +------+------+
    |......|......|
    |..+---+---+..|
    |..|       |..|
    |..| +---+ |..|
    |..| |...| |..|
    |..| +---+ |..|
    |..|       |..|
    |..+-------+..|
    |.............|
    +-------------+

    */

    // the reference area and barycenter
    double ref_area = 0.0;
    double * ref_barycenter = NULL;

    check_overlap(
      (double[]){0,3,3,-3,-3,0,0,-2,-2,2,2,0},
      (double[]){3,3,-3,-3,3,3,2,2,-2,-2,2,2}, 12,
      (double[]){1,1,-1,-1}, (double[]){1,-1,-1,1}, 4,
      ref_area, ref_barycenter);
  }

  {
    /*

    +                   +
    |\                 /|
    |..\             /..|
    |....\         /....|
    |.....+       +.....|
    |.....|       |.....|
    |.....|   +   |.....|
    |.....| /   \ |.....|
    |...../-------\.....|
    |.../...........\...|
    |..+-------------+..|
    |...................|
    +-------------------+

    */

    // the reference area and barycenter
    double ref_area;
    double ref_barycenter[3] = {0.0, 0.0, 0.0};
    {
      struct grid_cell Tri_a =
        generate_cell_deg(
          (double[]){90,180,90},
          (double[]){90,87.5,87.5}, gc_edges, 3);
      struct grid_cell Tri_b =
        generate_cell_deg(
          (double[]){90,180,90},
          (double[]){90,88,88}, gc_edges, 3);
      ref_area = yac_huiliers_area_info(Tri_a, ref_barycenter, 1.0) +
                 yac_huiliers_area_info(Tri_b, ref_barycenter, -1.0);
      if (ref_area < 0.0) {
        ref_area = -ref_area;
        ref_barycenter[0] = -ref_barycenter[0];
        ref_barycenter[1] = -ref_barycenter[1];
        ref_barycenter[2] = -ref_barycenter[2];
      }
      normalise_vector(ref_barycenter);
      yac_free_grid_cell(&Tri_b);
      yac_free_grid_cell(&Tri_a);
    }

    check_overlap(
      (double[]){0,90,180,270,270,180,90,0},
      (double[]){88,88,88,88,87,87,87,87}, 8,
      (double[]){90,180,90},
      (double[]){90,87.5,87.5}, 3,
      ref_area, ref_barycenter);
  }

  {
    /*

    +---------+---------+
    |.........|.........|
    |.........|.........|
    |.........|.........|
    |.....+---+---+.....|
    |.....|       |.....|
    |.....|   +   |.....|
    |.....| /   \ |.....|
    |...../-------\.....|
    |.../...........\...|
    |..+-------------+..|
    |...................|
    +-------------------+

    */

    // the reference area and barycenter
    double ref_area = 0.0;
    double ref_barycenter[3] = {0.0, 0.0, 0.0};
    {
      struct grid_cell Quad =
        generate_cell_deg(
          (double[]){90,180,180,90},
          (double[]){87.5,87.5,89,89}, gc_edges, 4);
      ref_area = yac_huiliers_area_info(Quad, ref_barycenter, 1.0);
      normalise_vector(ref_barycenter);
      yac_free_grid_cell(&Quad);
    }

    check_overlap(
      (double[]){315,270,180,90,0,315,
                 315,0,90,180,270,315},
      (double[]){89,89,89,89,89,89,
                 86,86,86,86,86,86}, 12,
      (double[]){90,90,180},
      (double[]){90,87.5,87.5}, 3,
      ref_area, ref_barycenter);
  }

  {
    /*

    +---------+---------+
    |.........|.........|
    |.........|......+..|
    |.........|...../|..|
    |.....+---+---/..|..|
    |.....|     / |..|..|
    +.....+   /   +..|..+
    |.....| /     |..|..|
    |...../---+---+..|..|
    |.../............|..|
    |..+-------------+..|
    |...................|
    +---------+---------+

    */

    // the reference area and barycenter
    for (int i = 0; i < 4; ++i) {

      double tri_rotation = (double)i * 90.0;

      double ref_area = 0.0;
      double ref_barycenter[3] = {0.0, 0.0, 0.0};
      {
        struct grid_cell Poly =
          generate_cell_deg(
            (double[]){0+tri_rotation,90+tri_rotation,180+tri_rotation,
                       180+tri_rotation,135+tri_rotation,90+tri_rotation,
                       45+tri_rotation,0+tri_rotation},
            (double[]){87.5,87.5,87.5,89,89,89,89,89}, gc_edges, 8);
        ref_area = yac_huiliers_area_info(Poly, ref_barycenter, 1.0);
        normalise_vector(ref_barycenter);
        yac_free_grid_cell(&Poly);
      }

      check_overlap(
        (double[]){315,270,225,180,135,90,45,0,315,
                   315,0,45,90,135,180,225,270,315},
        (double[]){89,89,89,89,89,89,89,89,89,
                   86,86,86,86,86,86,86,86,86}, 18,
        (double[]){0+tri_rotation,90+tri_rotation,180+tri_rotation},
        (double[]){87.5,87.5,87.5}, 3,
        ref_area, ref_barycenter);
    }
  }

  {
    /*

    +---------+---------+
    |.........|.........|
    |..+------|------+..|
    |..|......|......|..|
    |..|..+---+---+..|..|
    |..|..|       |..|..|
    +..|..+       +..|..+
    |..|..|       |..|..|
    |..|..+---+---+..|..|
    |..|.............|..|
    |..+-------------+..|
    |...................|
    +---------+---------+

    */

    // the reference area and barycenter
    double ref_area = 0.0;
    double ref_barycenter[3];
    {
      struct grid_cell Square =
        generate_cell_deg(
          (double[]){0,90,180,270},
          (double[]){87.5,87.5,87.5,87.5}, gc_edges, 4);
      struct grid_cell Polygon =
        generate_cell_deg(
          (double[]){0,45,90,135,180,225,270,315},
          (double[]){89,89,89,89,89,89,89,89}, gc_edges, 8);
      ref_area = yac_huiliers_area(Square) - yac_huiliers_area(Polygon);
      yac_free_grid_cell(&Polygon);
      yac_free_grid_cell(&Square);
      LLtoXYZ_deg(90, 90, ref_barycenter);
    }

    check_overlap(
      (double[]){315,270,225,180,135,90,45,0,315,
                 315,0,45,90,135,180,225,270,315},
      (double[]){89,89,89,89,89,89,89,89,89,
                 86,86,86,86,86,86,86,86,86}, 18,
      (double[]){0,90,180,270},
      (double[]){87.5,87.5,87.5,87.5}, 4,
      ref_area, ref_barycenter);
  }

  {
    /*

    +---------+---------+
    |.........|.........|
    |.........|.........|
    |.........|.........|
    |.....+---+---+.....|
    |.....|       |.....|
    +.....+       +.....+
    |.....|       |.....|
    |.....+---+---+.....|
    |...................|
    |...................|
    |...................|
    +---------+---------+

       +------+------+
       |.............|
       |.............|
       |.....+++.....|
       +.....+ +.....+
       |.....+++.....|
       |......|......|
       |......|......|
       +------+------+

    */

    // the reference area and barycenter
    double ref_area;
    double ref_barycenter[3];
    {
      struct grid_cell Polygon =
        generate_cell_deg(
          (double[]){315,270,225,180,135,90,45,0,315,
                     315,0,45,90,135,180,225,270,315},
          (double[]){88,88,88,88,88,88,88,88,88,
                     87,87,87,87,87,87,87,87,87}, gc_edges, 18);
      ref_area = yac_huiliers_area(Polygon);
      yac_free_grid_cell(&Polygon);
      LLtoXYZ_deg(90, 90, ref_barycenter);
    }

    check_overlap(
      (double[]){135,90,45,0,315,270,225,180,135,
                 135,180,225,270,315,0,45,90,135},
      (double[]){88,88,88,88,88,88,88,88,88,
                 86,86,86,86,86,86,86,86,86}, 18,
      (double[]){315,270,225,180,135,90,45,0,315,
                 315,0,45,90,135,180,225,270,315},
      (double[]){89,89,89,89,89,89,89,89,89,
                 87,87,87,87,87,87,87,87,87}, 18,
      ref_area, ref_barycenter);
  }

  { // test reproducing an issue found in toy_multi example
    struct grid_cell cells[2] =
      {{.coordinates_xyz =
          (double[3][3])
            {{-0.11609962857526578,
               0.1297445660881493,
              -0.98472697932740894},
             {-0.1873171704582404,
               0.13116089479097318,
              -0.97350351685504954},
             {-0.16162842294547278,
               0.064722709158966607,
              -0.98472697932740894}},
        .edge_type =
          (enum yac_edge_type[3])
            {GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE},
        .num_corners = 3,
        .array_size = 3},
       {.coordinates_xyz =
          (double[4][3])
            {{-0.10626911810438046,
               0.10117661128547845,
              -0.98917650996478101},
             {-0.10872011017156734,
               0.098538164069449888,
              -0.98917650996478101},
             {-0.12667440386975459,
               0.11481098733454023,
              -0.98527764238894122},
             {-0.12381864923052141,
               0.11788515390505606,
              -0.98527764238894122}},
        .edge_type =
          (enum yac_edge_type[4])
            {LAT_CIRCLE_EDGE,
             LON_CIRCLE_EDGE,
             LAT_CIRCLE_EDGE,
             LON_CIRCLE_EDGE},
        .num_corners = 4,
        .array_size = 4}};
    struct grid_cell overlap_buffer =
      {.coordinates_xyz = NULL,
       .edge_type = NULL,
       .num_corners = 0,
       .array_size = 0};
    yac_cell_clipping(1, &cells[0], cells[1], &overlap_buffer);
    double ref_overlap_areas = yac_huiliers_area(overlap_buffer);
    double ref_overlap_barycenters[3] = {0.0, 0.0, 0.0};
    check_overlap_polygons(cells, ref_overlap_areas, ref_overlap_barycenters);
    yac_free_grid_cell(&overlap_buffer);
  }

  { // test reproducing an issue found in toy_multi example
    struct grid_cell cells[2] =
      {{.coordinates_xyz =
          (double[4][3])
            {{0.12226322617998441,
              0.0060064071448744901,
              0.99247953459870997},
             {0.12207899817177321,
              0.0090050879009709785,
              0.99247953459870997},
             {0.097751558641606923,
              0.0072105881536315003,
              0.99518472667219682},
             {0.097899074391390048,
              0.0048094731201957178,
              0.99518472667219682}},
        .edge_type =
          (enum yac_edge_type[4])
            {LAT_CIRCLE_EDGE,
             LON_CIRCLE_EDGE,
             LAT_CIRCLE_EDGE,
             LON_CIRCLE_EDGE},
        .num_corners = 4,
        .array_size = 4},
       {.coordinates_xyz =
          (double[4][3])
            {{ 0.13173419827179747,
              -0.0024316641355669943,
               0.99128209305687476},
             { 0.12602575200455376,
               0.0081957056315718202,
               0.99199311501687715},
             { 0.11521731465958367,
               0.0019341931455107773,
               0.99333842636812875},
             { 0.12147807203778695,
              -0.0083897691028181811,
               0.99255865810962707}},
        .edge_type =
          (enum yac_edge_type[4])
            {GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE},
        .num_corners = 4,
        .array_size = 4}};
    struct grid_cell overlap_buffer =
      {.coordinates_xyz = NULL,
       .edge_type = NULL,
       .num_corners = 0,
       .array_size = 0};
    yac_cell_clipping(1, &cells[0], cells[1], &overlap_buffer);
    double ref_overlap_areas = yac_huiliers_area(overlap_buffer);
    double ref_overlap_barycenters[3] = {0.0, 0.0, 0.0};
    yac_free_grid_cell(&overlap_buffer);
    check_overlap_polygons(cells, ref_overlap_areas, ref_overlap_barycenters);
  }

  { // test reproducing an issue found in toy_multi example
    struct grid_cell cells[2] =
      {{.coordinates_xyz =
          (double[4][3])
            {{-0.67635734277981696,
              -0.34044967630401168,
               0.65317284295377664},
             {-0.66779858328675834,
              -0.35694577933893479,
               0.65317284295377664},
             {-0.65346055329338182,
              -0.34928194263987855,
               0.67155895484701855},
             {-0.66183555116521386,
              -0.33314002067992043,
               0.67155895484701855}},
        .edge_type =
          (enum yac_edge_type[4])
            {LAT_CIRCLE_EDGE,
             LON_CIRCLE_EDGE,
             LAT_CIRCLE_EDGE,
             LON_CIRCLE_EDGE},
        .num_corners = 4,
        .array_size = 4},
       {.coordinates_xyz =
          (double[4][3])
            {{-0.659061024045835,
              -0.36232186404812883,
               0.65906102404583489},
             {-0.66524659712730871,
              -0.33896007142593126,
               0.66524659712730871},
              {-0.64171186251713952,
               -0.34818611452313608,
                0.68335372623412638},
              {-0.63542071129861122,
               -0.37199386619312252,
                0.67665433063526625}},
        .edge_type =
          (enum yac_edge_type[4])
            {GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE},
        .num_corners = 4,
        .array_size = 4}};
    struct grid_cell overlap_buffer =
      {.coordinates_xyz = NULL,
       .edge_type = NULL,
       .num_corners = 0,
       .array_size = 0};
    yac_cell_clipping(1, &cells[0], cells[1], &overlap_buffer);
    double ref_overlap_areas = yac_huiliers_area(overlap_buffer);
    double ref_overlap_barycenters[3] = {0.0, 0.0, 0.0};
    yac_free_grid_cell(&overlap_buffer);
    check_overlap_polygons(cells, ref_overlap_areas, ref_overlap_barycenters);
  }

  { // testing issue found with FESOM grid

    /*

              +
              / \
            /   \
            +     +
        -           -
      +                 +
        -           -
            +     +
            \   /
              \ /
              +
    */

    // the reference area and barycenter
    double ref_area = 0.0;
    double ref_barycenter[3];
    {
      struct grid_cell Polygon =
        generate_cell_deg(
          (double[]){0.0,0.5,2.0,0.5,0.0,-0.5,-2.0,-0.5},
          (double[]){2.0,0.5,0.0,-0.5,-2.0,-0.5,0.0,0.5}, gc_edges, 8);
      ref_area = yac_huiliers_area(Polygon);
      yac_free_grid_cell(&Polygon);
      LLtoXYZ_deg(0, 0, ref_barycenter);
    }
    check_overlap(
      (double[]){0.0,0.5,2.0,0.5,0.0,-0.5,-2.0,-0.5},
      (double[]){2.0,0.5,0.0,-0.5,-2.0,-0.5,0.0,0.5}, 8,
      (double[]){-3.0,3.0,3.0,-3.0},
      (double[]){3.0,3.0,-3.0,-3.0}, 4,
      ref_area, ref_barycenter);
  }

  { // testing issue found with FESOM grid

    /*

              +
              / \
            /   \
            +     +
        -           -
      +                 +
        -           -
            +     +
            \   /
              \ /
              +
    */

    // the reference area and barycenter
    double ref_area = 0.0;
    double ref_barycenter[3];
    {
      struct grid_cell Polygon =
        generate_cell_deg(
          (double[]){-0.5,0.5,0.5,-0.5},
          (double[]){0.5,0.5,-0.5,-0.5}, gc_edges, 4);
      ref_area = yac_huiliers_area(Polygon);
      yac_free_grid_cell(&Polygon);
      LLtoXYZ_deg(0, 0, ref_barycenter);
    }
    check_overlap(
      (double[]){0.0,0.5,2.0,0.5,0.0,-0.5,-2.0,-0.5},
      (double[]){2.0,0.5,0.0,-0.5,-2.0,-0.5,0.0,0.5}, 8,
      (double[]){-0.5,0.5,0.5,-0.5},
      (double[]){0.5,0.5,-0.5,-0.5}, 4,
      ref_area, ref_barycenter);
  }

  {
    struct grid_cell src_cell =
        {.coordinates_xyz =
          (double[4][3])
            {{5.06731853971386536628e-01,
              5.06731853971386536628e-01,
              6.97456562333055085645e-01},
            {4.56966311668627278575e-01,
              6.28960169645094047119e-01,
              6.28960169645094047119e-01},
            {5.77350269189625842081e-01,
              5.77350269189625731059e-01,
              5.77350269189625842081e-01},
            {6.28960169645094047119e-01,
              4.56966311668627389597e-01,
              6.28960169645094047119e-01}},
        .edge_type =
          (enum yac_edge_type[4])
            {GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE},
        .num_corners = 4,
        .array_size = 4};
    struct grid_cell tgt_cell =
        {.coordinates_xyz =
          (double[3][3])
            {{5.26518818746262273756e-01,
              4.92120108892697694092e-01,
              6.93250122199397744716e-01},
            {5.28392100138220688343e-01,
              4.99031313570251600087e-01,
              6.86854814781020395209e-01},
            {5.20491669253990818511e-01,
              4.99252122636297701597e-01,
              6.92701768642426274347e-01}},
        .edge_type =
          (enum yac_edge_type[3])
            {GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE,
             GREAT_CIRCLE_EDGE},
        .num_corners = 3,
        .array_size = 3};

    double areas[2];
    double barycenter[3];
    yac_compute_overlap_info(
      1, &src_cell, tgt_cell, &areas[0], NULL);
    yac_compute_overlap_info(
      1, &src_cell, tgt_cell, &areas[1], &barycenter);

    if (fabs(areas[0] - areas[1]) > (areas[0] * 1e-6))
      PUT_ERR("error in yac_compute_overlap_info");
  }

  return TEST_EXIT_CODE;
}

static void check_overlap_polygons(
  struct grid_cell * Polygons, double ref_area,
  double const * ref_barycenter) {

  double area_tolerance = MAX(ref_area * 1e-6, 1e-10);

  for (int test_order = 0; test_order < 2; ++test_order) {
    double area;
    if (ref_barycenter) {
      double barycenter[3] = {0.0,0.0,0.0};
      yac_compute_overlap_info(
        1, &Polygons[test_order], Polygons[test_order^1],
        &area, &barycenter);
      if (fabs(ref_area-area) > area_tolerance)
        PUT_ERR("error in yac_compute_overlap_info "
                "area difference is too large.\n");
      if ((barycenter[0] != barycenter[0]) ||
          (barycenter[1] != barycenter[1]) ||
          (barycenter[2] != barycenter[2]))
        PUT_ERR("error in yac_compute_overlap_info "
                "barycenter contains nan's.\n");
      if ((ref_barycenter[0] != 0.0) &&
          (ref_barycenter[1] != 0.0) &&
          (ref_barycenter[2] != 0.0)) {
        double angle_diff = get_vector_angle(ref_barycenter, barycenter);
        if (fabs(angle_diff) > 1e-9)
          PUT_ERR("error in yac_compute_overlap_info "
                  "barycenter difference is too large.\n");
      }
    }
    yac_compute_overlap_areas(
      1, &Polygons[test_order], Polygons[test_order^1], &area);
    if (fabs(ref_area-area) > area_tolerance)
      PUT_ERR("error in yac_compute_overlap_areas "
              "area difference is too large.\n");
  }
}

static void check_overlap(
  double * lon_a, double * lat_a, int count_a,
  double * lon_b, double * lat_b, int count_b, double ref_area,
  double const * ref_barycenter) {

  for (int order_a = -1; order_a <=1; order_a += 2) {
    for (int order_b = -1; order_b <=1; order_b += 2) {
      for (int start_idx_a = 0; start_idx_a < count_a; ++start_idx_a) {
        for (int start_idx_b = 0; start_idx_b < count_b; ++start_idx_b) {

          double curr_lon_a[count_a], curr_lat_a[count_a];
          double curr_lon_b[count_b], curr_lat_b[count_b];

          for (int i = 0; i < count_a; ++i) {
            int idx = (start_idx_a + order_a * i + count_a) % count_a;
            curr_lon_a[i] = lon_a[idx];
            curr_lat_a[i] = lat_a[idx];
          }
          for (int i = 0; i < count_b; ++i) {
            int idx = (start_idx_b + order_b * i + count_b) % count_b;
            curr_lon_b[i] = lon_b[idx];
            curr_lat_b[i] = lat_b[idx];
          }

          struct grid_cell Polygons[2] = {
            generate_cell_deg(curr_lon_a, curr_lat_a, gc_edges, count_a),
            generate_cell_deg(curr_lon_b, curr_lat_b, gc_edges, count_b)};

          check_overlap_polygons(Polygons, ref_area, ref_barycenter);

          for (int i = 0; i < 2; ++i) yac_free_grid_cell(&Polygons[i]);
        }
      }
    }
  }
}
