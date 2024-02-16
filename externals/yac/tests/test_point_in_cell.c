/**
 * @file test_point_in_cell.c
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

#include "tests.h"
#include "grid.h"
#include "geometry.h"
#include "test_common.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

int main (void) {

   enum yac_edge_type gc_edges[8] =
     {GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
      GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
      GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
      GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE};
   enum yac_edge_type latlon_edges[4] =
     {LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE};

   double point_3d[3];
   struct grid_cell cell;

   { // simple inside test
      LLtoXYZ(0, 0, point_3d);
      cell = generate_cell_deg((double[]){-1,1,1,-1},
                               (double[]){-1,-1,1,1}, gc_edges, 4);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (simple inside test)\n")
      yac_free_grid_cell(&cell);
   }

   { // simple outside test
      LLtoXYZ_deg(0, 1.5, point_3d);
      cell = generate_cell_deg((double[]){-1,1,1,-1},
                               (double[]){-1,-1,1,1}, gc_edges, 4);

      if (yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (simple outside test)\n")
      yac_free_grid_cell(&cell);
   }

   { // test point on edge
      LLtoXYZ_deg(0, 1, point_3d);
      cell = generate_cell_deg((double[]){-1,1,1,-1},
                               (double[]){-1,-1,1,1}, gc_edges, 4);
      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (test point on edge)\n")
      yac_free_grid_cell(&cell);
   }

   { // test near the pole
      LLtoXYZ_deg(22.5, 85.1, point_3d);
      cell = generate_cell_deg((double[]){0,45,45,0},
                               (double[]){85,85,88,88}, latlon_edges, 4);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (lon-lat cell near the pole)\n")
      yac_free_grid_cell(&cell);
   }

   { // test near the pole
      LLtoXYZ_deg(22.5, 85.1, point_3d);
      cell = generate_cell_deg((double[]){0,45,45,0},
                               (double[]){85,85,88,88}, gc_edges, 4);

      if (yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (great circle cell near the pole)\n")
      yac_free_grid_cell(&cell);
   }

   { // simple test at the pole
      LLtoXYZ_deg(45, 88.3, point_3d);
      cell = generate_cell_deg((double[]){0,90,180,270},
                               (double[]){88,88,88,88},
                               (enum yac_edge_type[]){
                                  LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
                                  LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE}, 4);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (lon-lat cell on the pole)\n")
      yac_free_grid_cell(&cell);
   }

   { // simple test at the pole
      LLtoXYZ_deg(45, 88.3, point_3d);
      cell = generate_cell_deg((double[]){0,90,180,270},
                               (double[]){88,88,88,88}, gc_edges, 4);

      if (yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (great circle cell on the pole)\n")
      yac_free_grid_cell(&cell);
   }

   { // point and cell at pole
      for (int i = -1; i <= 1; i += 2) {
         for (int j = -1; j <= 1; j += 2) {
            LLtoXYZ_deg(45, 90.0 * (double)i, point_3d);
            double lat_coords[5] = {88.0, 88.0, 88.0, 88.0, 88.0};
            for (int k = 0; k < 5; ++k) lat_coords[k] *= (double)j;
            cell = generate_cell_deg((double[]){0,60,120,180,240},
                                     lat_coords, gc_edges, 5);

            if (yac_point_in_cell(point_3d, cell) ^ (i == j))
               PUT_ERR("error in point_in_cell (great circle cell on the pole)\n")
            yac_free_grid_cell(&cell);
         }
      }
   }

   { // test point on cell corner and edge
      LLtoXYZ_deg(1, 1, point_3d);
      cell = generate_cell_deg((double[]){0,1,1,0},
                               (double[]){0,0,1,1}, latlon_edges, 4);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (point on corner of lon-lat cell)\n")

      LLtoXYZ_deg(0.5, 0, point_3d);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (point on edge of lon-lat cell)\n")
      yac_free_grid_cell(&cell);
   }

   { // test point on cell corner
      LLtoXYZ_deg(1, 1, point_3d);
      cell = generate_cell_deg((double[]){0,1,0,-1},
                               (double[]){0,1,2,1}, gc_edges, 4);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (point on corner of great circle cell)\n")

      LLtoXYZ_deg(0, 2, point_3d);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (point on corner of great circle cell)\n")

      LLtoXYZ_deg(0, 1, point_3d);

      if (!yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (point in cell)\n")

      LLtoXYZ_deg(0, 3, point_3d);

      if (yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell (point in cell)\n")
      yac_free_grid_cell(&cell);
   }

   { // test point outside of cell
      LLtoXYZ(1.5, 0.5, point_3d);
      cell = generate_cell_deg((double[]){0,1,1,0},
                               (double[]){0,0,1,1}, latlon_edges, 4);

      if (yac_point_in_cell(point_3d, cell))
         PUT_ERR("error in point_in_cell\n")
      yac_free_grid_cell(&cell);
   }

   { // test point on longitude edge of cell


      for (unsigned i = 0; i < 2; ++i) {
         for (int j = 0; j < 3; ++j) {

            enum yac_edge_type edge_type[2][3] =
               {{GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE},
                {LON_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}};

            double point_lat[3] = {40, 45, 50};

            LLtoXYZ_deg(0.0, point_lat[j], point_3d);
            cell = generate_cell_deg((double[]){0,0,5},
                                     (double[]){40,50,45}, edge_type[i], 3);

            if (!yac_point_in_cell(point_3d, cell))
               PUT_ERR("error in point_in_cell\n")
            yac_free_grid_cell(&cell);
         }
      }
   }

   { // test point on latitude edge of cell

      for (unsigned i = 0; i < 2; ++i) {

         LLtoXYZ_deg(5.0, 10, point_3d);
         double y_coords[2][3] = {{10,10,15}, {10,10,5}};
         cell = generate_cell_deg((double[]){0,10,0}, y_coords[i],
                                  (enum yac_edge_type[])
                                    {LAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
                                     GREAT_CIRCLE_EDGE}, 3);

         if (!yac_point_in_cell(point_3d, cell))
            PUT_ERR("error in point_in_cell\n")
         yac_free_grid_cell(&cell);
      }
   }

   { // test point equator edge of cell

      for (unsigned i = 0; i < 2; ++i) {

         LLtoXYZ_deg(5.0, 0, point_3d);
         double y_coords[2][3] = {{0,0,5}, {0,0,-5}};
         cell = generate_cell_deg((double[]){0,10,0},
                                  y_coords[i], gc_edges, 3);

         if (!yac_point_in_cell(point_3d, cell))
            PUT_ERR("error in point_in_cell\n")
         yac_free_grid_cell(&cell);
      }
   }

   { // test point outside of cell


      for (unsigned i = 0; i < 2; ++i) {

         LLtoXYZ_deg(0.0, 11, point_3d);
         enum yac_edge_type edge_type[2][3] =
            {{LON_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE},
             {GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}};
         cell = generate_cell_deg((double[]){0,0,5},
                                  (double[]){5,10,20},
                                  edge_type[i], 3);

         if (yac_point_in_cell(point_3d, cell))
            PUT_ERR("error in point_in_cell\n")
         yac_free_grid_cell(&cell);
      }
   }

   { // test point inside of concave cell

      for (unsigned i = 0; i < 2; ++i) {

         LLtoXYZ_deg(0.0, 7.0, point_3d);
         enum yac_edge_type edge_type[2][4] =
            {{LON_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
              GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE},
             {GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
              GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}};
         cell = generate_cell_deg((double[]){0, 0, -5, 5},
                                  (double[]){0, 5, 10, 10},
                                  edge_type[i], 4);

         if (!yac_point_in_cell(point_3d, cell))
            PUT_ERR("error in point_in_cell\n")
         yac_free_grid_cell(&cell);
      }
   }

   { // test point outside of concave cell

      for (unsigned i = 0; i < 2; ++i) {

         LLtoXYZ_deg(-2.5, 6.0, point_3d);
         enum yac_edge_type edge_type[2][4] =
            {{LON_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
              GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE},
             {GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE,
              GREAT_CIRCLE_EDGE, GREAT_CIRCLE_EDGE}};
         cell = generate_cell_deg((double[]){0, 0, -5, 5},
                                  (double[]){0, 5, 10, 10},
                                  edge_type[i], 4);

         if (yac_point_in_cell(point_3d, cell))
            PUT_ERR("error in point_in_cell\n")
         yac_free_grid_cell(&cell);
      }
   }

   { // from Renes clipping tests

      double x_coords[] = { 0.0, 0.5, 0.0,-0.5};
      double y_coords[] = {-0.5, 0.0, 0.5, 0.0};
      cell = generate_cell_deg(x_coords, y_coords, gc_edges, 4);

      for (unsigned i = 0; i < 4; ++i) {
         LLtoXYZ_deg(x_coords[i], y_coords[i], point_3d);
         if (!yac_point_in_cell(point_3d, cell))
            PUT_ERR("error in point_in_cell\n")
      }
      yac_free_grid_cell(&cell);
   }

   { // example from toy

      cell = generate_cell_rad((double[]){-2.4523357834544868,
                                          -2.5132741228718340,
                                          -2.5742124622891822},
                               (double[]){-0.8955349680874743,
                                          -0.8292526951666540,
                                          -0.8955349680874748},
                               gc_edges, 3);

      LLtoXYZ(-2.5132741229357842, -0.8511807007280388, point_3d);
      if (!yac_point_in_cell(point_3d, cell))
            PUT_ERR("error in point_in_cell\n");
      yac_free_grid_cell(&cell);
   }

   { // example from scrip toy

      cell = generate_cell_3d((double[][3]){{-0.98109801397208807,
                                             -1.4373158006630782e-09,
                                             -0.19351146472502487},
                                            {-0.98296571213463857,
                                             -0.016290939681434788,
                                             -0.18306560040580719},
                                            {-0.98550150905720235,
                                             -0.016447047624685553,
                                             -0.16886761166786302},
                                            {-0.98739170499602125,
                                             -1.9584538881393463e-09,
                                             -0.15829599143708675},
                                            {-0.98550150905854683,
                                             0.016447044383749162,
                                             -0.16886761197567171},
                                            {-0.98296571214686901,
                                             0.016290936203869018,
                                             -0.18306560064960375}},
                               gc_edges, 6);

      if (!yac_point_in_cell((double[]){-0.98437791422272558,
                                        -1.2055152618042087e-16,
                                        -0.17606851504603632}, cell))
            PUT_ERR("error in point_in_cell\n");
      yac_free_grid_cell(&cell);
   }

   { // test degenerated triangles (with one very short edge)

      double lon[3] = {0,0 + yac_angle_tol,0};
      double lat[3] = {0,0,1};
      struct {
        double lon, lat, xyz[3];
        int is_inside;
      } test_points[] =
        {{.lon = lon[0],              .lat = lat[0], .is_inside = 1},
         {.lon = lon[1],              .lat = lat[1], .is_inside = 1},
         {.lon = lon[2],              .lat = lat[2], .is_inside = 1},
         {.lon = yac_angle_tol * 0.5, .lat = 0,      .is_inside = 1},
         {.lon = yac_angle_tol * 0.5, .lat = 0.5,    .is_inside = 1},
         {.lon = -0.01,               .lat = 0,      .is_inside = 0},
         {.lon = 0.01,                .lat = 0,      .is_inside = 0},
         {.lon = 0,                   .lat = 1.001,  .is_inside = 0},
         {.lon = -0.01,               .lat = 0.5,    .is_inside = 0}};
      enum{num_test_points = sizeof(test_points) / sizeof(test_points[0])};
      for (size_t i = 0; i < num_test_points; ++i)
        LLtoXYZ_deg(test_points[i].lon, test_points[i].lat, test_points[i].xyz);

      // bounding circle that contains all points and the cell
      struct bounding_circle bnd_circle;
      LLtoXYZ_deg(0.0, 0.5, bnd_circle.base_vector);
      bnd_circle.inc_angle =
         sin_cos_angle_new(sin(1.0 * YAC_RAD), cos(1.0 * YAC_RAD));
      bnd_circle.sq_crd = DBL_MAX;

      for (int direction = -1; direction <= 1; direction += 2) {

        for (int start = 0; start < 3; ++start) {

          double temp_lon[3], temp_lat[3];
          for (int edge_idx = start, i = 0; i < 3; ++i, edge_idx += direction) {
            temp_lon[i] = lon[(edge_idx + 3) % 3];
            temp_lat[i] = lat[(edge_idx + 3) % 3];
          }
          cell = generate_cell_deg(temp_lon, temp_lat, gc_edges, 3);

          for (size_t i = 0; i < num_test_points; ++i) {
            if (yac_point_in_cell(test_points[i].xyz, cell) !=
                test_points[i].is_inside)
               PUT_ERR("error in point_in_cell\n")
            if (yac_point_in_cell2(test_points[i].xyz, cell, bnd_circle) !=
                test_points[i].is_inside)
               PUT_ERR("error in point_in_cell2\n")
          }

          yac_free_grid_cell(&cell);
        }
      }
   }

   { // test degenerated triangles (with three very short edges)

      double lon[3] = {0,0 + yac_angle_tol,0};
      double lat[3] = {0,0,0 + yac_angle_tol};
      struct {
        double lon, lat, xyz[3];
        int is_inside;
      } test_points[] =
        {{.lon = lon[0],              .lat = lat[0], .is_inside = 1},
         {.lon = lon[1],              .lat = lat[1], .is_inside = 1},
         {.lon = lon[2],              .lat = lat[2], .is_inside = 1},
         {.lon = yac_angle_tol * 0.5, .lat = 0,      .is_inside = 1},
         {.lon = 0,                   .lat = yac_angle_tol * 0.5,
                                                     .is_inside = 1},
         {.lon = -0.001,              .lat = 0,      .is_inside = 0},
         {.lon = 0.001,               .lat = 0,      .is_inside = 0},
         {.lon = 0,                   .lat = 0.001,  .is_inside = 0},
         {.lon = 0,                   .lat = -0.001, .is_inside = 0}};
      size_t num_test_points = sizeof(test_points) / sizeof(test_points[0]);
      for (size_t i = 0; i < num_test_points; ++i)
        LLtoXYZ_deg(test_points[i].lon, test_points[i].lat, test_points[i].xyz);

      // bounding circle that contains all points and the cell
      struct bounding_circle bnd_circle;
      LLtoXYZ_deg(0.0, 0.0, bnd_circle.base_vector);
      bnd_circle.inc_angle =
         sin_cos_angle_new(sin(0.5 * YAC_RAD), cos(0.5 * YAC_RAD));
      bnd_circle.sq_crd = DBL_MAX;

      for (int direction = -1; direction <= 1; direction += 2) {

        for (int start = 0; start < 3; ++start) {

          double temp_lon[3], temp_lat[3];
          for (int edge_idx = start, i = 0; i < 3; ++i, edge_idx += direction) {
            temp_lon[i] = lon[(edge_idx + 3) % 3];
            temp_lat[i] = lat[(edge_idx + 3) % 3];
          }
          cell = generate_cell_deg(temp_lon, temp_lat, gc_edges, 3);

          for (size_t i = 0; i < num_test_points; ++i) {
            if (yac_point_in_cell(test_points[i].xyz, cell) !=
                test_points[i].is_inside)
               PUT_ERR("error in point_in_cell\n")
            if (yac_point_in_cell2(test_points[i].xyz, cell, bnd_circle) !=
                test_points[i].is_inside)
               PUT_ERR("error in point_in_cell2\n")
          }

          yac_free_grid_cell(&cell);
        }
      }
   }

   { // example from HD

      cell = generate_cell_3d((double[][3]){{0.39177807766083478,
                                             0.48669359266172357,
                                             0.78079400915120067},
                                            {0.39106981008544894,
                                             0.48726288480996138,
                                             0.78079400915120067},
                                            {0.39178019297633498,
                                             0.48814800354785914,
                                             0.77988448312789571},
                                            {0.39248974712806789,
                                             0.48757767727376555,
                                             0.77988448312789571}},
                               latlon_edges, 4);

      if (!yac_point_in_cell((double[]){0.39178019303094858,
                                        0.48814800355996663,
                                        0.77988448309288172}, cell))
            PUT_ERR("error in point_in_cell\n");
      yac_free_grid_cell(&cell);
   }

   return TEST_EXIT_CODE;
}

