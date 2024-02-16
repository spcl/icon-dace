/**
 * @file test_dist_grid_pair.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "dist_grid.h"
#include "yac_mpi.h"
#include "dist_grid_utils.h"
#include "io_utils.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

static void get_basic_grid_data(
  char * filename, size_t * num_cells, size_t * num_vertices,
  size_t * num_edges);

static void check_indices(
  const_yac_int_pointer ids, size_t * indices,
  size_t count, size_t global_count);

static void check_global_ids(yac_int const * ids, size_t count,
                             int * mask, int global_count,
                             MPI_Comm comm);

static void check_unmasked_global_ids(yac_int const * ids, size_t count,
                                      int * mask, int global_count,
                                      int * global_core_mask, MPI_Comm comm);

static void check_dist_owner(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  enum yac_location location, size_t global_size, MPI_Comm comm);

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  set_even_io_rank_list(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  {
    char * filenames[] ={"icon_grid_R02B02.nc", "icon_grid_R02B03_G.nc"};
    struct basic_grid_data grid_data[2];
    char const * grid_names[2] = {"grid_1", "grid2"};

    for (int i = 0; i < 2; ++i)
      grid_data[i] =
        read_icon_grid_information_parallel2(
          filenames[i], MPI_COMM_WORLD);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[0], grid_data[0]),
       yac_basic_grid_new(grid_names[1], grid_data[1])};

    struct {
      yac_int * global_ids;
      int * core_mask;
    } global_ids[2][3];

    for (int i = 0; i < 2; ++i) {
      global_ids[i][0].global_ids = grid_data[i].cell_ids;
      global_ids[i][0].core_mask  = grid_data[i].core_cell_mask;
      global_ids[i][1].global_ids = grid_data[i].vertex_ids;
      global_ids[i][1].core_mask  = grid_data[i].core_vertex_mask;
      global_ids[i][2].global_ids = grid_data[i].edge_ids;
      global_ids[i][2].core_mask  = grid_data[i].core_edge_mask;
    }

    for (int k = 0; k < 8; ++k) {

      for (int i = 0; i < 2; ++i) {
        struct basic_grid_data * curr_basic_grid_data =
          yac_basic_grid_get_data(grids[i]);
        if (k & (1 << 0)) {
          curr_basic_grid_data->cell_ids = NULL;
          curr_basic_grid_data->core_cell_mask = NULL;
        } else {
          curr_basic_grid_data->cell_ids = global_ids[i][0].global_ids;
          curr_basic_grid_data->core_cell_mask = global_ids[i][0].core_mask;
        }
        if (k & (1 << 1)) {
          curr_basic_grid_data->vertex_ids = NULL;
          curr_basic_grid_data->core_vertex_mask = NULL;
        } else {
          curr_basic_grid_data->vertex_ids = global_ids[i][1].global_ids;
          curr_basic_grid_data->core_vertex_mask = global_ids[i][1].core_mask;
        }
        if (k & (1 << 2)) {
          curr_basic_grid_data->edge_ids = NULL;
          curr_basic_grid_data->core_edge_mask = NULL;
        } else {
          curr_basic_grid_data->edge_ids = global_ids[i][2].global_ids;
          curr_basic_grid_data->core_edge_mask = global_ids[i][2].core_mask;
        }
      }

      struct dist_grid_pair * grid_pair =
        yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

      for (int i = 0; i < 2; ++i) {

        // test yac_dist_grid_pair_get_dist_grid
        struct dist_grid * dist_grid =
          yac_dist_grid_pair_get_dist_grid(grid_pair, grid_names[i]);

        enum yac_location locations[] = {CELL, CORNER, EDGE};
        size_t counts[3];

        get_basic_grid_data(filenames[i], &counts[0], &counts[1], &counts[2]);

        const struct const_basic_grid_data * grid_data =
          yac_dist_grid_get_basic_grid_data(dist_grid);

        const_yac_int_pointer ids[] =
          {grid_data->ids[CELL], grid_data->ids[CORNER], grid_data->ids[EDGE]};

        for (size_t j = 0; j < 3; ++j) {

          // test yac_dist_grid_get_local_count
          uint64_t local_count =
            (uint64_t)yac_dist_grid_get_local_count(dist_grid, locations[j]);

          uint64_t global_count;
          yac_mpi_call(
            MPI_Allreduce(&local_count, &global_count, 1, MPI_UINT64_T, MPI_SUM,
                          MPI_COMM_WORLD), MPI_COMM_WORLD);
          if ((uint64_t)(counts[j]) != global_count)
            PUT_ERR("error in yac_dist_grid_get_local_count\n");

          // test yac_dist_grid_get_local_unmasked_points
          size_t * indices, count;
          yac_dist_grid_get_local_unmasked_points(
            dist_grid, (struct interp_field){.location = locations[j],
                                             .coordinates_idx = SIZE_MAX,
                                             .masks_idx = SIZE_MAX},
            &indices, &count);

          if (count != local_count)
            PUT_ERR("error in yac_dist_grid_get_local_unmasked_points\n");

          check_indices(ids[j], indices, count, counts[j]);

          free(indices);
        }
      }

      yac_dist_grid_pair_delete(grid_pair);

      for (int i = 0; i < 2; ++i) {
        struct basic_grid_data * curr_basic_grid_data =
          yac_basic_grid_get_data(grids[i]);
        if (k & (1 << 0)) {
          free(curr_basic_grid_data->cell_ids);
          free(curr_basic_grid_data->core_cell_mask);
          curr_basic_grid_data->cell_ids = NULL;
          curr_basic_grid_data->core_cell_mask = NULL;
        }
        if (k & (1 << 1)) {
          free(curr_basic_grid_data->vertex_ids);
          free(curr_basic_grid_data->core_vertex_mask);
          curr_basic_grid_data->vertex_ids = NULL;
          curr_basic_grid_data->core_vertex_mask = NULL;
        }
        if (k & (1 << 2)) {
          free(curr_basic_grid_data->edge_ids);
          free(curr_basic_grid_data->core_edge_mask);
          curr_basic_grid_data->edge_ids = NULL;
          curr_basic_grid_data->core_edge_mask = NULL;
        }
      }
    }

    for (int i = 0; i < 2; ++i) {
      yac_basic_grid_delete(grids[i]);
      for (int j = 0; j < 3; ++j) {
        free(global_ids[i][j].global_ids);
        free(global_ids[i][j].core_mask);
      }
    }
  }

  if (comm_size >= 4) { // test with artificial grids

    // we only need 4 processes
    int do_test = comm_rank < 4;
    MPI_Comm comm;
    yac_mpi_call(
      MPI_Comm_split(MPI_COMM_WORLD, do_test, 0, &comm), MPI_COMM_WORLD);

    if (do_test) {

      int comm_rank, comm_size;
      yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
      yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

      int is_tgt = comm_rank >= 2;

      double coordinates_x[2][5] = {{0.0,1.0,2.0,3.0,4.0}, {0.5,1.5,2.5,3.5,-1.0}};
      double coordinates_y[2][4] = {{0.0,1.0,2.0,3.0}, {0.5,1.5,2.5,-1.0}};
      size_t num_cells[2][2] = {{4,3},{3,2}};
      size_t local_start[4][2] = {{0,0},{2,0},{0,0},{2,0}};
      size_t local_count[4][2] = {{2,3},{2,3},{2,2},{1,2}};
      int core_vertex_mask[2][20] = {{0,0,0,0,0,
                                      0,1,1,1,0,
                                      0,1,1,1,0,
                                      0,0,0,0,0},
                                     {0,0,0,0,
                                      0,1,1,0,
                                      0,0,0,0}};
      int core_edge_mask[2][31] = {{0,0,0,1,0,1,0,1,0,
                                    1,0,1,1,1,1,1,1,0,
                                    1,0,1,1,1,1,1,1,0,
                                    0,0,0,0},
                                    {0,0,0,1,0,1,0,
                                     1,0,1,1,1,1,0,
                                     0,0,0}};
      int with_halo = 0;
      for (int i = 0; i < 2; ++i){
        for (int j = 0; j < 5; ++j) coordinates_x[i][j] *= YAC_RAD;
        for (int j = 0; j < 4; ++j) coordinates_y[i][j] *= YAC_RAD;
      }

      struct basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg2d(
          coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
          local_start[comm_rank], local_count[comm_rank], with_halo);
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        grid_data.core_vertex_mask[i] =
          core_vertex_mask[is_tgt][grid_data.vertex_ids[i]];
      for (size_t i = 0; i < grid_data.num_edges; ++i)
        grid_data.core_edge_mask[i] =
          core_edge_mask[is_tgt][grid_data.edge_ids[i]];
      char const * grid_names[2] = {"src_grid", "tgt_grid"};
      struct yac_basic_grid * grids[2] =
        {yac_basic_grid_new(grid_names[is_tgt], grid_data),
         yac_basic_grid_empty_new(grid_names[is_tgt^1])};


      struct dist_grid_pair * grid_pair =
        yac_dist_grid_pair_new(grids[is_tgt], grids[is_tgt^1], comm);

      // check dist owner information
      {
        size_t num_points[2][3] = {{12,20,31}, {6,12,17}};
        enum yac_location locations[3] = {CELL, CORNER, EDGE};
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 3; ++j)
            check_dist_owner(
              grid_pair, grid_names[i], locations[j], num_points[i][j], comm);
      }

      // check basic grid data of local part of the grid
      {
        int ref_num_cells[2] = {12, 6};
        int ref_num_vertices[2] = {20, 12};
        int ref_num_edges[2] = {31, 17};
        int ref_cell_to_vertex[2][12][4] =
          {{{0,1,6,5},{1,2,7,6},{2,3,8,7},{3,4,9,8},
            {5,6,11,10},{6,7,12,11},{7,8,13,12},{8,9,14,13},
            {10,11,16,15},{11,12,17,16},{12,13,18,17},{13,14,19,18}},
           {{0,1,5,4},{1,2,6,5},{2,3,7,6},
            {4,5,9,8},{5,6,10,9},{6,7,11,10}}};
        int ref_cell_to_edge[2][12][4] =
          {{{0,3,9,1},{2,5,11,3},{4,7,13,5},{6,8,15,7},
            {9,12,18,10},{11,14,20,12},{13,16,22,14},{15,17,24,16},
            {18,21,27,19},{20,23,28,21},{22,25,29,23},{24,26,30,25}},
           {{0,3,7,1},{2,5,9,3},{4,6,11,5},
            {7,10,14,8},{9,12,15,10},{11,13,16,12}}};
        enum yac_edge_type ref_edge_type[2][31] =
          {{LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE},
           {LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE}};

        int cell_mask[12];
        int vertex_mask[20];
        int edge_mask[31];

        for (int i = 0; i < 2; ++i) {

          struct const_basic_grid_data * grid_data =
            yac_dist_grid_get_basic_grid_data(
              yac_dist_grid_pair_get_dist_grid(grid_pair, grid_names[i]));

          check_global_ids(grid_data->ids[CELL], grid_data->count[CELL],
                           cell_mask, ref_num_cells[i], comm);
          check_global_ids(grid_data->ids[CORNER], grid_data->count[CORNER],
                           vertex_mask, ref_num_vertices[i], comm);
          check_global_ids(grid_data->ids[EDGE], grid_data->count[EDGE],
                           edge_mask, ref_num_edges[i], comm);

          for (size_t j = 0; j < grid_data->count[CELL]; ++j) {

            int curr_cell_id = (int)(grid_data->ids[CELL][j]);

            if (grid_data->num_vertices_per_cell[j] != 4)
              PUT_ERR("error in num_vertices_per_cell");

            for (int k = 0; k < 4; ++k) {
              if (ref_cell_to_vertex[i][curr_cell_id][k] !=
                  grid_data->ids[CORNER][
                    grid_data->cell_to_vertex[
                      grid_data->cell_to_vertex_offsets[j]+k]])
                PUT_ERR("error in cell_to_vertex");

              if (ref_cell_to_edge[i][curr_cell_id][k] !=
                  grid_data->ids[EDGE][
                    grid_data->cell_to_edge[
                      grid_data->cell_to_edge_offsets[j]+k]])
                PUT_ERR("error in cell_to_edge");
            }
          }

          for (size_t j = 0; j < grid_data->count[EDGE]; ++j)
            if (ref_edge_type[i][grid_data->ids[EDGE][j]] !=
                grid_data->edge_type[j])
              PUT_ERR("error in edge_type");

          // check yac_dist_grid_get_local_unmasked_points
          {
            size_t * local_vertex_indices, local_vertex_count;
            yac_dist_grid_get_local_unmasked_points(
              yac_dist_grid_pair_get_dist_grid(grid_pair, grid_names[i]),
              (struct interp_field){.location = CORNER,
                                    .coordinates_idx = SIZE_MAX,
                                    .masks_idx = SIZE_MAX},
              &local_vertex_indices, &local_vertex_count);

            yac_int global_vertex_indices[20];

            for (size_t j = 0; j < local_vertex_count; ++j)
              global_vertex_indices[j] =
                grid_data->ids[CORNER][local_vertex_indices[j]];

            check_unmasked_global_ids(
              global_vertex_indices, local_vertex_count,
              vertex_mask, ref_num_vertices[i], core_vertex_mask[i],
              comm);

            free(local_vertex_indices);
          }

          // check yac_dist_grid_global_to_local
          {
            yac_int global_edge_ids[ref_num_edges[i]];
            size_t local_edge_ids[ref_num_edges[i]];
            size_t edge_count = 0;
            for (int j = 0; j < ref_num_edges[i]; ++j)
              if (core_edge_mask[i][j])
                global_edge_ids[edge_count++] = j;

            yac_dist_grid_global_to_local(
              yac_dist_grid_pair_get_dist_grid(grid_pair, grid_names[i]), EDGE,
              global_edge_ids, edge_count, local_edge_ids);

            struct const_basic_grid_data * grid_data =
              yac_dist_grid_get_basic_grid_data(
                yac_dist_grid_pair_get_dist_grid(grid_pair, grid_names[i]));

            for (int j = 0; j < edge_count; ++j) {

              if (grid_data->ids[EDGE][local_edge_ids[j]] != global_edge_ids[j])
                PUT_ERR("ERROR in yac_dist_grid_global_to_local");

              if (ref_edge_type[i][global_edge_ids[j]] !=
                  grid_data->edge_type[local_edge_ids[j]])
                PUT_ERR("ERROR in yac_dist_grid_global_to_local");
            }
          }
        }
      }

      yac_dist_grid_pair_delete(grid_pair);
      yac_basic_grid_delete(grids[1]);
      yac_basic_grid_delete(grids[0]);
    }

    yac_mpi_call(MPI_Comm_free(&comm), MPI_COMM_WORLD);

  } else {
    PUT_ERR("insufficient number of processes");
  }

  if (comm_size >= 3) { // test with degenerated artificial grids

    // we only need 3 processes
    int do_test = comm_rank < 3;
    MPI_Comm comm;
    yac_mpi_call(
      MPI_Comm_split(MPI_COMM_WORLD, do_test, 0, &comm), MPI_COMM_WORLD);

    if (do_test) {

      int comm_rank, comm_size;
      yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
      yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

      int is_tgt = comm_rank == 2;

      struct basic_grid_data grid_data;
      if (is_tgt) {
        size_t nbr_vertices[2] = {5, 2};
        int cyclic[2] = {0,0};
        double lon_vertices[] = {-2.0,-1.0,0.0,1.0,2.0};
        double lat_vertices[] = {-0.5,0.5};
        grid_data =
          yac_generate_basic_grid_data_reg_2d_deg(
            nbr_vertices, cyclic, lon_vertices, lat_vertices);
      } else {
        size_t nbr_vertices[2] = {7,8};
        size_t nbr_cells[2] = {4,4};
        int num_vertices_per_cell[2][4] = {{5,3,3,0}, {3,6,3,0}};
        double x_vertices[2][8] = {{-1.5,-0.5,0.5,-2.0,0.0,-0.75,0.5},
                                   {-0.5,0.5,1.5,0.0,2.0,-0.75,0.5,1.5}};
        double y_vertices[2][8] = {{-0.5,-0.5,-0.5,0.0,0.0,0.5,0.5},
                                   {-0.5,-0.5,-0.5,0.0,0.0,0.5,0.5,0.5}};
        int cell_to_vertex[2][12] = {{0,1,4,5,3, 1,2,4, 4,6,5},
                                     {0,1,3, 1,2,4,7,6,3, 3,6,5}};
        static int core_cell_mask[2][4] = {{1,1,0,1},{0,1,1,1}};
        static int core_vertex_mask[2][8] = {{1,1,0,1,1,1,0},{0,1,1,1,1,0,1,1}};
        grid_data =
          yac_generate_basic_grid_data_unstruct_deg(
            nbr_vertices[comm_rank], nbr_cells[comm_rank],
            num_vertices_per_cell[comm_rank],
            x_vertices[comm_rank], y_vertices[comm_rank],
            cell_to_vertex[comm_rank]);

        grid_data.core_cell_mask = TO_POINTER(core_cell_mask[comm_rank]);
        grid_data.core_vertex_mask = TO_POINTER(core_vertex_mask[comm_rank]);
      }

      char const * grid_names[2] = {"src_grid", "tgt_grid"};
      struct yac_basic_grid * grids[2] =
        {yac_basic_grid_new(grid_names[is_tgt], grid_data),
         yac_basic_grid_empty_new(grid_names[is_tgt^1])};


      struct dist_grid_pair * grid_pair =
        yac_dist_grid_pair_new(grids[is_tgt], grids[is_tgt^1], comm);

      // check dist owner information
      {
        size_t num_points[2][3] = {{5,10,13}, {4,10,13}};
        enum yac_location locations[3] = {CELL, CORNER, EDGE};
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 3; ++j)
            check_dist_owner(
              grid_pair, grid_names[i], locations[j], num_points[i][j], comm);
      }

      yac_dist_grid_pair_delete(grid_pair);
      yac_basic_grid_delete(grids[1]);
      yac_basic_grid_delete(grids[0]);
    }

    yac_mpi_call(MPI_Comm_free(&comm), MPI_COMM_WORLD);

  } else {
    PUT_ERR("insufficient number of processes");
  }

  // test on a single process
  {

    // we only need 1 processes
    int do_test = comm_rank < 1;
    MPI_Comm comm;
    yac_mpi_call(
      MPI_Comm_split(MPI_COMM_WORLD, do_test, 0, &comm), MPI_COMM_WORLD);

    if (do_test) {

      int comm_rank, comm_size;
      yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
      yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

      double coordinates_x[2][5] = {{0.0,1.0,2.0,3.0,4.0}, {0.5,1.5,2.5,3.5,-1.0}};
      double coordinates_y[2][4] = {{0.0,1.0,2.0,3.0}, {0.5,1.5,2.5,-1.0}};
      size_t num_cells[2][2] = {{4,3},{3,2}};
      size_t local_start[2][2] = {{0,0},{0,0}};
      size_t local_count[2][2] = {{4,3},{3,2}};
      int with_halo = 0;
      for (int i = 0; i < 2; ++i){
        for (int j = 0; j < 5; ++j) coordinates_x[i][j] *= YAC_RAD;
        for (int j = 0; j < 4; ++j) coordinates_y[i][j] *= YAC_RAD;
      }

      struct basic_grid_data grid_data[2] =
        {yac_generate_basic_grid_data_reg2d(
           coordinates_x[0], coordinates_y[0], num_cells[0],
           local_start[0], local_count[0], with_halo),
         yac_generate_basic_grid_data_reg2d(
           coordinates_x[1], coordinates_y[1], num_cells[1],
           local_start[1], local_count[1], with_halo)};
      char const * grid_names[2] = {"src_grid", "tgt_grid"};
      struct yac_basic_grid * grids[2] =
        {yac_basic_grid_new(grid_names[0], grid_data[0]),
         yac_basic_grid_new(grid_names[1], grid_data[1])};


      struct dist_grid_pair * grid_pair =
        yac_dist_grid_pair_new(grids[0], grids[1], comm);

      // check basic grid data of local part of the grid
      {
        int ref_num_cells[2] = {12, 6};
        int ref_num_vertices[2] = {20, 12};
        int ref_num_edges[2] = {31, 17};
        int ref_cell_to_vertex[2][12][4] =
          {{{0,1,6,5},{1,2,7,6},{2,3,8,7},{3,4,9,8},
            {5,6,11,10},{6,7,12,11},{7,8,13,12},{8,9,14,13},
            {10,11,16,15},{11,12,17,16},{12,13,18,17},{13,14,19,18}},
           {{0,1,5,4},{1,2,6,5},{2,3,7,6},
            {4,5,9,8},{5,6,10,9},{6,7,11,10}}};
        int ref_cell_to_edge[2][12][4] =
          {{{0,3,9,1},{2,5,11,3},{4,7,13,5},{6,8,15,7},
            {9,12,18,10},{11,14,20,12},{13,16,22,14},{15,17,24,16},
            {18,21,27,19},{20,23,28,21},{22,25,29,23},{24,26,30,25}},
           {{0,3,7,1},{2,5,9,3},{4,6,11,5},
            {7,10,14,8},{9,12,15,10},{11,13,16,12}}};
        enum yac_edge_type ref_edge_type[2][31] =
          {{LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE},
           {LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LON_CIRCLE_EDGE, LON_CIRCLE_EDGE, LAT_CIRCLE_EDGE, LAT_CIRCLE_EDGE,
            LAT_CIRCLE_EDGE}};

         int cell_mask[12];
         int vertex_mask[20];
         int edge_mask[31];

         for (int i = 0; i < 2; ++i) {

           struct const_basic_grid_data * grid_data =
             yac_dist_grid_get_basic_grid_data(
               yac_dist_grid_pair_get_dist_grid(grid_pair, grid_names[i]));

           check_global_ids(grid_data->ids[CELL], grid_data->count[CELL],
                            cell_mask, ref_num_cells[i], comm);
           check_global_ids(grid_data->ids[CORNER], grid_data->count[CORNER],
                            vertex_mask, ref_num_vertices[i], comm);
           check_global_ids(grid_data->ids[EDGE], grid_data->count[EDGE],
                            edge_mask, ref_num_edges[i], comm);

           for (size_t j = 0; j < grid_data->count[CELL]; ++j) {

             int curr_cell_id = (int)(grid_data->ids[CELL][j]);

             if (grid_data->num_vertices_per_cell[j] != 4)
               PUT_ERR("error in num_vertices_per_cell");

             for (int k = 0; k < 4; ++k) {
               if (ref_cell_to_vertex[i][curr_cell_id][k] !=
                   grid_data->ids[CORNER][
                     grid_data->cell_to_vertex[
                       grid_data->cell_to_vertex_offsets[j]+k]])
                 PUT_ERR("error in cell_to_vertex");

               if (ref_cell_to_edge[i][curr_cell_id][k] !=
                   grid_data->ids[EDGE][
                     grid_data->cell_to_edge[
                       grid_data->cell_to_edge_offsets[j]+k]])
                 PUT_ERR("error in cell_to_edge");
             }
           }

           for (size_t j = 0; j < grid_data->count[EDGE]; ++j)
             if (ref_edge_type[i][grid_data->ids[EDGE][j]] !=
                 grid_data->edge_type[j])
               PUT_ERR("error in edge_type");
         }
      }

      yac_dist_grid_pair_delete(grid_pair);
      yac_basic_grid_delete(grids[1]);
      yac_basic_grid_delete(grids[0]);
    }

    yac_mpi_call(MPI_Comm_free(&comm), MPI_COMM_WORLD);
  }

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void get_basic_grid_data(
  char * filename, size_t * num_cells, size_t * num_vertices,
  size_t * num_edges) {

  int ncid;

  yac_nc_open(filename, NC_NOWRITE, &ncid);

  int dimid;
  HANDLE_ERROR(nc_inq_dimid(ncid, "cell", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, num_cells));
  HANDLE_ERROR(nc_inq_dimid(ncid, "vertex", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, num_vertices));
  HANDLE_ERROR(nc_inq_dimid(ncid, "edge", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, num_edges));
}

static void check_indices(
  const_yac_int_pointer ids, size_t * indices,
  size_t count, size_t global_count) {

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  size_t local_start =
    (global_count * (size_t)comm_rank + (size_t)comm_size - 1) /
    (size_t)comm_size;
  size_t next_local_start =
    (global_count * (size_t)(comm_rank + 1) + (size_t)comm_size - 1) /
    (size_t)comm_size;
  size_t local_size = next_local_start - local_start;

  Xt_int * src_global_ids = xmalloc(count * sizeof(*indices));
  for (size_t i = 0; i < count; ++i)
    src_global_ids[i] = (Xt_int)(ids[indices[i]]);

  struct Xt_stripe dst_stripe =
    {.start = (Xt_int)local_start, .stride = 1, .nstrides = (int)local_size};

  Xt_idxlist src_idxlist = xt_idxvec_new(src_global_ids, (int)count);
  Xt_idxlist dst_idxlist = xt_idxstripes_new(&dst_stripe, 1);
  Xt_xmap xmap = xt_xmap_dist_dir_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
  Xt_redist redist = xt_redist_p2p_new(xmap, yac_int_dt);

  Xt_int * collected_ids = xmalloc(local_size * sizeof(*collected_ids));

  xt_redist_s_exchange1(redist, src_global_ids, collected_ids);

  for (size_t i = 0; i < local_size; ++i)
    if (collected_ids[i] != (Xt_int)(i + local_start)) PUT_ERR("missing id");

  free(collected_ids);
  xt_redist_delete(redist);
  xt_xmap_delete(xmap);
  xt_idxlist_delete(dst_idxlist);
  xt_idxlist_delete(src_idxlist);
  free(src_global_ids);
}

static void check_global_ids(yac_int const * ids, size_t count,
                             int * mask, int global_count,
                             MPI_Comm comm) {

  for (int i = 0; i < global_count; ++i) mask[i] = 0;

  for (size_t i = 0; i < count; ++i) mask[ids[i]] = 1;

  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, mask, global_count, MPI_INT, MPI_MAX, comm), comm);

  for (int i = 0; i < global_count; ++i) if (!mask[i]) PUT_ERR("missing id");
}

static void check_unmasked_global_ids(yac_int const * ids, size_t count,
                                      int * mask, int global_count,
                                      int * global_core_mask, MPI_Comm comm) {

  for (int i = 0; i < global_count; ++i) mask[i] = 0;

  for (size_t i = 0; i < count; ++i) mask[ids[i]] = 1;

  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, mask, global_count, MPI_INT, MPI_SUM, comm), comm);

  for (int i = 0; i < global_count; ++i) {
    if (mask[i] > 1) PUT_ERR("multiple owners");
    if (mask[i] != global_core_mask[i])
      PUT_ERR("masked data is owned by a process");
  }
}

static void check_dist_owner(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  enum yac_location location, size_t global_size, MPI_Comm comm) {

  struct const_basic_grid_data * grid_data =
    yac_dist_grid_get_basic_grid_data(
      yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name));

  size_t local_size = grid_data->count[location];
  yac_int const * local_ids = grid_data->ids[location];

  size_t points[local_size];
  for (size_t i = 0; i < local_size; ++i) points[i] = i;

  // request dist owners for all points in the local part of the dist grid
  int local_ranks[local_size];
  yac_dist_grid_pair_determine_dist_owner(
    grid_pair, grid_name, points, local_size, location, local_ranks);

  // combine the result from all processes
  int global_ranks[global_size];
  for (size_t i = 0; i < global_size; ++i) global_ranks[i] = INT_MAX;
  for (size_t i = 0; i < local_size; ++i)
    global_ranks[local_ids[i]] = local_ranks[i];
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, global_ranks, (int)global_size, MPI_INT, MPI_MIN, comm),
    comm);

  // every point should be available on at least one process
  for (size_t i = 0; i < global_size; ++i)
    if (global_ranks[i] == INT_MAX)
      PUT_ERR("ERROR in yac_dist_grid_pair_determine_dist_owner");

  // for any point the dist owner should be the same on all ranks
  for (size_t i = 0; i < local_size; ++i)
    if (global_ranks[local_ids[i]] != local_ranks[i])
      PUT_ERR("ERROR in yac_dist_grid_pair_determine_dist_owner");

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  // all points should be available on the dist owners
  for (size_t i = 0; i < global_size; ++i) {
    if (global_ranks[i] == comm_rank) {
      int found = 0;
      for(size_t j = 0; (j < local_size) && !found; ++j)
        if (local_ids[j] == (size_t)i) found = 1;
      if (!found)
        PUT_ERR("ERROR in yac_dist_grid_pair_determine_dist_owner");
    }
  }
}
