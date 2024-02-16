/**
 * @file test_quicksort.c
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

#include <stdlib.h>
#include <stdio.h>
#include "dist_grid_utils.h"
#include "tests.h"

static void check_basic_grid_data(
  struct basic_grid_data grid_data, struct basic_grid_data ref_grid_data);

int main (void) {

  { // test grid without halo
    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double coordinates_y[] = {-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
    size_t const num_cells[2] = {9, 8};
    size_t local_start[2] = {3, 3};
    size_t local_count[2] = {3, 2};
    int with_halo = 0;

    for (size_t i = 0; i < sizeof(coordinates_x)/sizeof(coordinates_x[0]); ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i < sizeof(coordinates_y)/sizeof(coordinates_y[0]); ++i)
      coordinates_y[i] *= YAC_RAD;

    double ref_vertex_coordinates[12][3];
    yac_int ref_cell_ids[6] = {30, 31, 32,
                               39, 40, 41};
    yac_int ref_vertex_ids[12] = {33, 34, 35, 36,
                                  43, 44, 45, 46,
                                  53, 54, 55, 56};
    yac_int ref_edge_ids[17] = {63, 64, 65, 66, 67, 68, 70,
                                82, 83, 84, 85, 86, 87, 89,
                                101, 103, 105};
    size_t ref_num_cells = 6;
    size_t ref_num_vertices = 12;
    size_t ref_num_edges = 17;
    int ref_core_cell_mask[6] = {1, 1, 1,
                                 1, 1, 1};
    int ref_core_vertex_mask[12] = {1, 1, 1, 1,
                                    1, 1, 1, 1,
                                    1, 1, 1, 1};
    int ref_core_edge_mask[17] = {1, 1, 1, 1, 1, 1, 1,
                                  1, 1, 1, 1, 1, 1, 1,
                                  1, 1, 1};
    int ref_num_vertices_per_cell[6] = {4, 4, 4,
                                        4, 4, 4};
    int ref_num_cells_per_vertex[12] = {1, 2, 2, 1,
                                        2, 4, 4, 2,
                                        1, 2, 2, 1};
    size_t ref_cell_to_vertex[6][4] = {{0,1,5,4}, {1,2,6,5}, {2,3,7,6},
                                       {4,5,9,8}, {5,6,10,9}, {6,7,11,10}};
    size_t ref_cell_to_vertex_offsets[6] = {0,4,8,12,16,20};
    size_t ref_cell_to_edge[6][4] = {{0,3,7,1}, {2,5,9,3}, {4,6,11,5},
                                     {7,10,14,8}, {9,12,15,10}, {11,13,16,12}};
    size_t ref_cell_to_edge_offsets[6] = {0,4,8,12,16,20};
    size_t ref_vertex_to_cell[] =
      {0, 0,1, 1,2, 2, 0,3, 0,1,3,4, 1,2,4,5, 2,5, 3, 3,4, 4,5, 5};
    size_t ref_vertex_to_cell_offsets[12] = {0,1,3,5,6,8,12,16,18,19,21,23};
    size_t ref_edge_to_vertex[17][2] = {{0,1}, {0,4}, {1,2}, {1,5}, {2,3}, {2,6}, {3,7},
                                        {4,5}, {4,8}, {5,6}, {5,9}, {6,7}, {6,10}, {7,11},
                                        {8,9}, {9,10}, {10,11}};
    enum yac_edge_type ref_edge_type[17] = {
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE};

    for (size_t i = 0, k = 0; i < local_count[1] + 1; ++i)
      for (size_t j = 0; j < local_count[0] + 1; ++j, ++k)
        LLtoXYZ(coordinates_x[local_start[0]+j],
                coordinates_y[local_start[1]+i],
                &(ref_vertex_coordinates[k][0]));

    struct basic_grid_data ref_grid_data = {
      .vertex_coordinates     = ref_vertex_coordinates,
      .cell_ids               = ref_cell_ids,
      .vertex_ids             = ref_vertex_ids,
      .edge_ids               = ref_edge_ids,
      .num_cells              = ref_num_cells,
      .num_vertices           = ref_num_vertices,
      .num_edges              = ref_num_edges,
      .core_cell_mask         = ref_core_cell_mask,
      .core_vertex_mask       = ref_core_vertex_mask,
      .core_edge_mask         = ref_core_edge_mask,
      .num_vertices_per_cell  = ref_num_vertices_per_cell,
      .num_cells_per_vertex   = ref_num_cells_per_vertex,
      .cell_to_vertex         = &(ref_cell_to_vertex[0][0]),
      .cell_to_vertex_offsets = ref_cell_to_vertex_offsets,
      .cell_to_edge           = &(ref_cell_to_edge[0][0]),
      .cell_to_edge_offsets   = ref_cell_to_edge_offsets,
      .vertex_to_cell         = ref_vertex_to_cell,
      .vertex_to_cell_offsets = ref_vertex_to_cell_offsets,
      .edge_to_vertex         = (size_t_2_pointer)&(ref_edge_to_vertex[0][0]),
      .edge_type              = ref_edge_type,
      .num_total_cells        = ref_num_cells,
      .num_total_vertices     = ref_num_vertices,
      .num_total_edges        = ref_num_edges
    };

    // generate basic grid data

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(coordinates_x, coordinates_y,
                                         num_cells, local_start, local_count,
                                         with_halo);

    // check basic grid data
    check_basic_grid_data(grid_data, ref_grid_data);
  }

  { // test grid with halo
    double coordinates_x[] = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
    size_t const num_cells[2] = {6, 6};
    size_t local_start[2] = {2, 2};
    size_t local_count[2] = {2, 2};
    int with_halo = 1;

    for (size_t i = 0; i < sizeof(coordinates_x)/sizeof(coordinates_x[0]); ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i < sizeof(coordinates_y)/sizeof(coordinates_y[0]); ++i)
      coordinates_y[i] *= YAC_RAD;

    double ref_vertex_coordinates[25][3];
    yac_int ref_cell_ids[16] = {7,8,9,10, 13,14,15,16, 19,20,21,22, 25,26,27,28};
    yac_int ref_vertex_ids[25] = {8,9,10,11,12, 15,16,17,18,19, 22,23,24,25,26,
                                  29,30,31,32,33, 36,37,38,39,40};
    yac_int ref_edge_ids[40] = {15,16,17,18,19,20,21,22,24,
                                28,29,30,31,32,33,34,35,37,
                                41,42,43,44,45,46,47,48,50,
                                54,55,56,57,58,59,60,61,63,
                                67,69,71,73};
    size_t ref_num_cells = 16;
    size_t ref_num_vertices = 25;
    size_t ref_num_edges = 40;
    int ref_core_cell_mask[16] = {0,0,0,0,
                                  0,1,1,0,
                                  0,1,1,0,
                                  0,0,0,0};
    int ref_core_vertex_mask[25] = {0,0,0,0,0,
                                    0,1,1,1,0,
                                    0,1,1,1,0,
                                    0,1,1,1,0,
                                    0,0,0,0,0};
    int ref_core_edge_mask[40] = {0,0,0,0,0,0,0,0,0,
                                  0,0,1,1,1,1,0,1,0,
                                  0,0,1,1,1,1,0,1,0,
                                  0,0,1,0,1,0,0,0,0,
                                  0,0,0,0};
    int ref_num_vertices_per_cell[16] = {4,4,4,4,
                                         4,4,4,4,
                                         4,4,4,4,
                                         4,4,4,4};
    int ref_num_cells_per_vertex[25] = {1,2,2,2,1,
                                        2,4,4,4,2,
                                        2,4,4,4,2,
                                        2,4,4,4,2,
                                        1,2,2,2,1};
    size_t ref_cell_to_vertex[16][4] =
      {{0,1,6,5}, {1,2,7,6}, {2,3,8,7}, {3,4,9,8},
       {5,6,11,10}, {6,7,12,11}, {7,8,13,12}, {8,9,14,13},
       {10,11,16,15}, {11,12,17,16}, {12,13,18,17}, {13,14,19,18},
       {15,16,21,20}, {16,17,22,21}, {17,18,23,22}, {18,19,24,23}};
    size_t ref_cell_to_vertex_offsets[16] = {0,4,8,12,
                                             16,20,24,28,
                                             32,36,40,44,
                                             48,52,56,60};
    size_t ref_cell_to_edge[16][4] =
      {{0,3,9,1}, {2,5,11,3}, {4,7,13,5}, {6,8,15,7},
       {9,12,18,10}, {11,14,20,12}, {13,16,22,14}, {15,17,24,16},
       {18,21,27,19}, {20,23,29,21}, {22,25,31,23}, {24,26,33,25},
       {27,30,36,28}, {29,32,37,30}, {31,34,38,32}, {33,35,39,34}};
    size_t ref_cell_to_edge_offsets[16] = {0,4,8,12,
                                           16,20,24,28,
                                           32,36,40,44,
                                           48,52,56,60};
    size_t ref_vertex_to_cell[] =
      {0, 0,1, 1,2, 2,3, 3,
       0,4, 0,1,4,5, 1,2,5,6, 2,3,6,7, 3,7,
       4,8, 4,5,8,9, 5,6,9,10, 6,7,10,11, 7,11,
       8,12, 8,9,12,13, 9,10,13,14, 10,11,14,15, 11,15,
       12, 12,13, 13,14, 14,15, 15};
    size_t ref_vertex_to_cell_offsets[25] = {0,
                                             1,3,5,7,8,
                                             10,14,18,22,24,
                                             26,30,34,38,40,
                                             42,46,50,54,56,
                                             57,59,61,63};
    size_t ref_edge_to_vertex[40][2] = {{0,1}, {0,5}, {1,2}, {1,6}, {2,3}, {2,7}, {3,4}, {3,8}, {4,9},
                                        {5,6}, {5,10}, {6,7}, {6,11}, {7,8}, {7,12}, {8,9}, {8,13}, {9,14},
                                        {10,11}, {10,15}, {11,12}, {11,16}, {12,13}, {12,17}, {13,14}, {13,18}, {14,19},
                                        {15,16}, {15,20}, {16,17}, {16,21}, {17,18},{17,22}, {18,19}, {18,23}, {19,24},
                                        {20,21}, {21,22}, {22,23}, {23,24}};
    enum yac_edge_type ref_edge_type[40] = {
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE};

#ifndef __NEC__
    for (size_t i = 0, k = 0; i < local_count[1] + 1 + 2; ++i)
      for (size_t j = 0; j < local_count[0] + 1 + 2; ++j, ++k)
#else
// NEC compiler optimizes the following loops out when
// CFLAGS='-O2 -finline-functions' if they are specified
// as above
    for (size_t i = 0, k = 0; i < local_count[1] + 3; ++i)
      for (size_t j = 0; j < local_count[0] + 3; ++j, ++k)
#endif
        LLtoXYZ(coordinates_x[local_start[0]-1+j],
                coordinates_y[local_start[1]-1+i],
                &(ref_vertex_coordinates[k][0]));

    struct basic_grid_data ref_grid_data = {
      .vertex_coordinates     = ref_vertex_coordinates,
      .cell_ids               = ref_cell_ids,
      .vertex_ids             = ref_vertex_ids,
      .edge_ids               = ref_edge_ids,
      .num_cells              = ref_num_cells,
      .num_vertices           = ref_num_vertices,
      .num_edges              = ref_num_edges,
      .core_cell_mask         = ref_core_cell_mask,
      .core_vertex_mask       = ref_core_vertex_mask,
      .core_edge_mask         = ref_core_edge_mask,
      .num_vertices_per_cell  = ref_num_vertices_per_cell,
      .num_cells_per_vertex   = ref_num_cells_per_vertex,
      .cell_to_vertex         = &(ref_cell_to_vertex[0][0]),
      .cell_to_vertex_offsets = ref_cell_to_vertex_offsets,
      .cell_to_edge           = &(ref_cell_to_edge[0][0]),
      .cell_to_edge_offsets   = ref_cell_to_edge_offsets,
      .vertex_to_cell         = ref_vertex_to_cell,
      .vertex_to_cell_offsets = ref_vertex_to_cell_offsets,
      .edge_to_vertex         = (size_t_2_pointer)&(ref_edge_to_vertex[0][0]),
      .edge_type              = ref_edge_type,
      .num_total_cells        = ref_num_cells,
      .num_total_vertices     = ref_num_vertices,
      .num_total_edges        = ref_num_edges
    };

    // generate basic grid data

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(coordinates_x, coordinates_y,
                                         num_cells, local_start, local_count,
                                         with_halo);

    // check basic grid data
    check_basic_grid_data(grid_data, ref_grid_data);
  }

  { // test grid, which is at the boundary of the global grid, with halo
    double coordinates_x[] = {-2.0, -1.0, 0.0, 1.0, 2.0};
    double coordinates_y[] = {-2.0, -1.0, 0.0, 1.0, 2.0};
    size_t const num_cells[2] = {4, 4};
    size_t local_start[2] = {0, 0};
    size_t local_count[2] = {2, 2};
    int with_halo = 1;

    for (size_t i = 0; i < sizeof(coordinates_x)/sizeof(coordinates_x[0]); ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i < sizeof(coordinates_y)/sizeof(coordinates_y[0]); ++i)
      coordinates_y[i] *= YAC_RAD;

    double ref_vertex_coordinates[16][3];
    yac_int ref_cell_ids[9] = {0,1,2, 4,5,6, 8,9,10};
    yac_int ref_vertex_ids[16] = {0,1,2,3, 5,6,7,8, 10,11,12,13, 15,16,17,18};
    yac_int ref_edge_ids[24] = {0,1,2,3,4,5,7,
                                9,10,11,12,13,14,16,
                                18,19,20,21,22,23,25,
                                27,29,31};
    size_t ref_num_cells = 9;
    size_t ref_num_vertices = 16;
    size_t ref_num_edges = 24;
    int ref_core_cell_mask[9] = {1,1,0,
                                 1,1,0,
                                 0,0,0};
    int ref_core_vertex_mask[16] = {1,1,1,0,
                                    1,1,1,0,
                                    1,1,1,0,
                                    0,0,0,0};
    int ref_core_edge_mask[24] = {1,1,1,1,0,1,0,
                                  1,1,1,1,0,1,0,
                                  1,0,1,0,0,0,0,
                                  0,0,0};
    int ref_num_vertices_per_cell[9] = {4,4,4,
                                        4,4,4,
                                        4,4,4};
    int ref_num_cells_per_vertex[16] = {1,2,2,1,
                                        2,4,4,2,
                                        2,4,4,2,
                                        1,2,2,1};
    size_t ref_cell_to_vertex[9][4] =
      {{0,1,5,4}, {1,2,6,5}, {2,3,7,6},
       {4,5,9,8}, {5,6,10,9}, {6,7,11,10},
       {8,9,13,12}, {9,10,14,13}, {10,11,15,14}};
    size_t ref_cell_to_vertex_offsets[9] = {0,4,8, 12,16,20, 24,28,32};
    size_t ref_cell_to_edge[9][4] =
      {{0,3,7,1}, {2,5,9,3}, {4,6,11,5},
       {7,10,14,8}, {9,12,16,10}, {11,13,18,12},
       {14,17,21,15}, {16,19,22,17}, {18,20,23,19}};
    size_t ref_cell_to_edge_offsets[9] = {0,4,8, 12,16,20, 24,28,32};
    size_t ref_vertex_to_cell[] =
      {0, 0,1, 1,2, 2,
       0,3, 0,1,3,4, 1,2,4,5, 2,5,
       3,6, 3,4,6,7, 4,5,7,8, 5,8,
       6, 6,7, 7,8, 8};
    size_t ref_vertex_to_cell_offsets[16] = {0,1,3,5,
                                             6,8,12,16,
                                             18,20,24,28,
                                             30,31,33,35};
    size_t ref_edge_to_vertex[24][2] = {{0,1}, {0,4}, {1,2}, {1,5}, {2,3}, {2,6}, {3,7},
                                        {4,5}, {4,8}, {5,6}, {5,9}, {6,7}, {6,10}, {7,11},
                                        {8,9}, {8,12}, {9,10}, {9,13}, {10,11}, {10,14}, {11,15},
                                        {12,13}, {13,14}, {14,15}};

    enum yac_edge_type ref_edge_type[24] = {
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
        LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE};

#ifndef __NEC__
    for (size_t i = 0, k = 0; i < local_count[1] + 1 + 1; ++i)
      for (size_t j = 0; j < local_count[0] + 1 + 1; ++j, ++k)
#else
// NEC compiler optimizes the following loops out when
// CFLAGS='-O2 -finline-functions' if they are specified
// as above
    for (size_t i = 0, k = 0; i < local_count[1] + 2; ++i)
      for (size_t j = 0; j < local_count[0] + 2; ++j, ++k)
#endif
        LLtoXYZ(coordinates_x[local_start[0]+j],
                coordinates_y[local_start[1]+i],
                &(ref_vertex_coordinates[k][0]));

    struct basic_grid_data ref_grid_data = {
      .vertex_coordinates     = ref_vertex_coordinates,
      .cell_ids               = ref_cell_ids,
      .vertex_ids             = ref_vertex_ids,
      .edge_ids               = ref_edge_ids,
      .num_cells              = ref_num_cells,
      .num_vertices           = ref_num_vertices,
      .num_edges              = ref_num_edges,
      .core_cell_mask         = ref_core_cell_mask,
      .core_vertex_mask       = ref_core_vertex_mask,
      .core_edge_mask         = ref_core_edge_mask,
      .num_vertices_per_cell  = ref_num_vertices_per_cell,
      .num_cells_per_vertex   = ref_num_cells_per_vertex,
      .cell_to_vertex         = &(ref_cell_to_vertex[0][0]),
      .cell_to_vertex_offsets = ref_cell_to_vertex_offsets,
      .cell_to_edge           = &(ref_cell_to_edge[0][0]),
      .cell_to_edge_offsets   = ref_cell_to_edge_offsets,
      .vertex_to_cell         = ref_vertex_to_cell,
      .vertex_to_cell_offsets = ref_vertex_to_cell_offsets,
      .edge_to_vertex         = (size_t_2_pointer)&(ref_edge_to_vertex[0][0]),
      .edge_type              = ref_edge_type,
      .num_total_cells        = ref_num_cells,
      .num_total_vertices     = ref_num_vertices,
      .num_total_edges        = ref_num_edges
    };

    // generate basic grid data

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(coordinates_x, coordinates_y,
                                         num_cells, local_start, local_count,
                                         with_halo);

    // check basic grid data
    check_basic_grid_data(grid_data, ref_grid_data);
  }

   return TEST_EXIT_CODE;
}


static void check_basic_grid_data(
  struct basic_grid_data grid_data, struct basic_grid_data ref_grid_data) {

  if (grid_data.num_cells != ref_grid_data.num_cells)
    PUT_ERR("wrong num_cells");
  if (grid_data.num_vertices != ref_grid_data.num_vertices)
    PUT_ERR("wrong num_vertices");
  if (grid_data.num_edges != ref_grid_data.num_edges)
    PUT_ERR("wrong num_edges");

  for (size_t i = 0; i < ref_grid_data.num_vertices; ++i)
    for (int j = 0; j < 3; ++j)
      if (fabs(grid_data.vertex_coordinates[i][j] -
               ref_grid_data.vertex_coordinates[i][j]) > 1e-9)
        PUT_ERR("wrong grid_data.vertex_coordinates");

  for (size_t i = 0; i < ref_grid_data.num_cells; ++i) {
    if (grid_data.cell_ids[i] != ref_grid_data.cell_ids[i])
      PUT_ERR("wrong cell_ids");
    if (grid_data.core_cell_mask[i] != ref_grid_data.core_cell_mask[i])
      PUT_ERR("wrong core_cell_mask");
    if (grid_data.num_vertices_per_cell[i] != ref_grid_data.num_vertices_per_cell[i])
      PUT_ERR("wrong num_vertices_per_cell");
    if (grid_data.cell_to_vertex_offsets[i] != ref_grid_data.cell_to_vertex_offsets[i])
      PUT_ERR("wrong cell_to_vertex_offsets");
    for (size_t j = 0; j < ref_grid_data.num_vertices_per_cell[i]; ++j)
      if (grid_data.cell_to_vertex[grid_data.cell_to_vertex_offsets[i]+j] !=
          ref_grid_data.cell_to_vertex[ref_grid_data.cell_to_vertex_offsets[i]+j])
        PUT_ERR("wrong cell_to_vertex");
    if (grid_data.cell_to_edge_offsets[i] != ref_grid_data.cell_to_edge_offsets[i])
      PUT_ERR("wrong cell_to_edge_offsets");
    for (size_t j = 0; j < ref_grid_data.num_vertices_per_cell[i]; ++j)
      if (grid_data.cell_to_edge[grid_data.cell_to_edge_offsets[i]+j] !=
          ref_grid_data.cell_to_edge[ref_grid_data.cell_to_edge_offsets[i]+j])
        PUT_ERR("wrong cell_to_edge");

  }
  for (size_t i = 0; i < ref_grid_data.num_vertices; ++i) {
    if (grid_data.vertex_ids[i] != ref_grid_data.vertex_ids[i])
      PUT_ERR("wrong vertex_ids");
    if (grid_data.core_vertex_mask[i] != ref_grid_data.core_vertex_mask[i])
      PUT_ERR("wrong core_vertex_mask");
    if (grid_data.num_cells_per_vertex[i] != ref_grid_data.num_cells_per_vertex[i])
      PUT_ERR("wrong num_cells_per_vertex");
    if (grid_data.vertex_to_cell_offsets[i] != ref_grid_data.vertex_to_cell_offsets[i])
      PUT_ERR("wrong vertex_to_cell_offsets");
    for (size_t j = 0; j < ref_grid_data.num_cells_per_vertex[i]; ++j)
      if (grid_data.vertex_to_cell[grid_data.vertex_to_cell_offsets[i]+j] !=
          ref_grid_data.vertex_to_cell[ref_grid_data.vertex_to_cell_offsets[i]+j])
        PUT_ERR("wrong vertex_to_cell");
  }
  for (size_t i = 0; i < ref_grid_data.num_edges; ++i) {
    if (grid_data.edge_ids[i] != ref_grid_data.edge_ids[i])
      PUT_ERR("wrong edge_ids");
    if (grid_data.core_edge_mask[i] != ref_grid_data.core_edge_mask[i])
      PUT_ERR("wrong core_edge_mask");
    if (grid_data.edge_type[i] != ref_grid_data.edge_type[i])
      PUT_ERR("wrong edge_type");
    if (!(((grid_data.edge_to_vertex[i][0] == ref_grid_data.edge_to_vertex[i][0]) &&
           (grid_data.edge_to_vertex[i][1] == ref_grid_data.edge_to_vertex[i][1])) ||
          ((grid_data.edge_to_vertex[i][0] == ref_grid_data.edge_to_vertex[i][1]) &&
           (grid_data.edge_to_vertex[i][1] == ref_grid_data.edge_to_vertex[i][0]))))
      PUT_ERR("wrong edge_to_vertex");
  }

  if (grid_data.num_total_cells != ref_grid_data.num_total_cells)
    PUT_ERR("wrong num_total_cells");
  if (grid_data.num_total_vertices != ref_grid_data.num_total_vertices)
    PUT_ERR("wrong num_total_vertices");
  if (grid_data.num_total_edges != ref_grid_data.num_total_edges)
    PUT_ERR("wrong num_total_edges");
}
