/**
 * @file test_interp_grid_parallel.c
 *
 * @copyright Copyright  (C)  2019 DKRZ, MPI-M
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
#include <string.h>
#include <float.h>

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "interp_grid.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"

static int check_global_ids(
  yac_int * global_ids, int count, int ref_global_count, MPI_Comm comm);
static int compare_global_ids(
  yac_int * global_ids, yac_int * ref_global_ids, size_t count);

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  set_even_io_rank_list(MPI_COMM_WORLD);

  { // general test using ICON grid files
    char * filenames[] ={"icon_grid_R02B02.nc", "icon_grid_R02B03_G.nc"};
    struct basic_grid_data grid_data[2];
    char const * grid_names[2] = {"R02B02", "R02B03"};

    for (int i = 0; i < 2; ++i)
      grid_data[i] =
        read_icon_grid_information_parallel2(
          filenames[i], MPI_COMM_WORLD);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[0], grid_data[0]),
       yac_basic_grid_new(grid_names[1], grid_data[1])};

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    if (yac_interp_grid_get_num_src_fields(interp_grid) != 1)
      PUT_ERR("error in yac_interp_grid_get_num_src_fields");

    if (yac_interp_grid_get_src_field_location(interp_grid, 0) !=
        src_fields[0].location)
      PUT_ERR("error in yac_interp_grid_get_src_field_location");

    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);

    for (int i = 0; i < 2; ++i) yac_basic_grid_delete(grids[i]);
  }

  if (comm_size >= 4) { // test with artificial grids

    // we only need 4 processes
    int do_test = comm_rank < 4;
    MPI_Comm comm;
    yac_mpi_call(
      MPI_Comm_split(MPI_COMM_WORLD, do_test, 0, &comm), MPI_COMM_WORLD);
    int comm_rank, comm_size;
    yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
    yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

    if (do_test) {

      {

        int is_tgt = comm_rank >= 2;

        double coordinates_x[2][5] = {{0.0,1.0,2.0,3.0,4.0}, {0.5,1.5,2.5,3.5,-1.0}};
        double coordinates_y[2][4] = {{0.0,1.0,2.0,3.0}, {0.5,1.5,2.5,-1.0}};
        size_t num_cells[2][2] = {{4,3},{3,2}};
        size_t local_start[4][2] = {{0,0},{2,0},{0,0},{2,0}};
        size_t local_count[4][2] = {{2,3},{2,3},{2,2},{1,2}};
        int with_halo = 0;
        for (int i = 0; i < 2; ++i){
          for (int j = 0; j < 5; ++j) coordinates_x[i][j] *= YAC_RAD;
          for (int j = 0; j < 4; ++j) coordinates_y[i][j] *= YAC_RAD;
        }

        struct basic_grid_data grid_data =
          yac_generate_basic_grid_data_reg2d(
            coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
            local_start[comm_rank], local_count[comm_rank], with_halo);

        char const * grid_names[2] = {"src_grid", "tgt_grid"};
        struct yac_basic_grid * grids[] =
          {yac_basic_grid_new(grid_names[is_tgt], grid_data),
          yac_basic_grid_empty_new(grid_names[is_tgt^1])};

        struct dist_grid_pair * grid_pair =
          yac_dist_grid_pair_new(grids[0], grids[1], comm);

        struct interp_field src_fields[] =
          {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
          {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
          {.location = EDGE, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
        size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
        struct interp_field tgt_field =
          {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

        struct interp_grid * interp_grid =
          yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                              num_src_fields, src_fields, tgt_field);

        // check

        size_t * tgt_indices, tgt_count;
        yac_interp_grid_get_tgt_points(interp_grid, &tgt_indices, &tgt_count);

        yac_int tgt_global_ids[12];
        yac_interp_grid_get_tgt_global_ids(
          interp_grid, tgt_indices, tgt_count, tgt_global_ids);

        double tgt_coordinates[12][3];
        yac_interp_grid_get_tgt_coordinates(
          interp_grid, tgt_indices, tgt_count, tgt_coordinates);

        {
          int count_ = (int)tgt_count, counts[4], displs[4];
          yac_mpi_call(
            MPI_Allgather(&count_, 1, MPI_INT, counts, 1, MPI_INT, comm), comm);
          for (int i = 0, accu = 0; i < comm_size; accu += counts[i++])
            displs[i] = accu;

          int global_count = 0;
          for (int i = 0; i < comm_size; ++i) global_count += counts[i];
          if (global_count != 12)
            PUT_ERR("error in yac_interp_grid_get_tgt_points");

          yac_int all_tgt_global_ids[12];
          yac_mpi_call(MPI_Gatherv(tgt_global_ids, count_, yac_int_dt,
                                  all_tgt_global_ids, counts, displs, yac_int_dt,
                                  0, comm), comm);
          double all_tgt_coordinates[12][3];
          for (int i = 0; i < comm_size; ++i) counts[i] *= 3;
          for (int i = 0; i < comm_size; ++i) displs[i] *= 3;
          yac_mpi_call(MPI_Gatherv(tgt_coordinates, 3*count_, MPI_DOUBLE,
                                  all_tgt_coordinates, counts, displs, MPI_DOUBLE,
                                  0, comm), comm);

          if (comm_rank == 0) {
            double ref_tgt_coords[12][3];
            for (size_t i = 0, k = 0; i < num_cells[1][1]+1; ++i)
              for (size_t j = 0; j < num_cells[1][0]+1; ++j, ++k)
                LLtoXYZ(
                  coordinates_x[1][j], coordinates_y[1][i],
                  &(ref_tgt_coords[k][0]));
            for (yac_int i = 0; i < 12; ++i) {
              if ((all_tgt_global_ids[i] >= 12) || (all_tgt_global_ids[i] < 0))
                PUT_ERR("error in yac_interp_grid_get_tgt_points");
              if (get_vector_angle(ref_tgt_coords[all_tgt_global_ids[i]],
                                  all_tgt_coordinates[i]) > yac_angle_tol)
                PUT_ERR("error in yac_interp_grid_get_tgt_points");
            }
          }
        }

        size_t ref_global_count[3] = {12, 20, 31};
        for (size_t src_field_idx = 0; src_field_idx < 3; ++src_field_idx) {

          if (yac_interp_grid_get_src_field_location(
                interp_grid, src_field_idx) != src_fields[src_field_idx].location)
            PUT_ERR("error in yac_interp_grid_get_src_field_location");

          size_t * src_indices, src_count;
          yac_interp_grid_get_src_points(
            interp_grid, src_field_idx, &src_indices, &src_count);

          yac_int src_global_ids[31];
          yac_interp_grid_get_src_global_ids(
            interp_grid, src_indices, src_count, src_field_idx, src_global_ids);

          if (check_global_ids(
                src_global_ids, src_count, ref_global_count[src_field_idx], comm))
            PUT_ERR("error in yac_interp_grid_get_src_points or "
                    "yac_interp_grid_get_src_global_ids");
          free(src_indices);
        }

        {
          double search_coords[15][3];
          for (size_t i = 0, k = 0; i < num_cells[1][1]+1; ++i)
            for (size_t j = 0; j < num_cells[1][0]+1; ++j, ++k)
              LLtoXYZ(
                coordinates_x[1][j], coordinates_y[1][i], &(search_coords[k][0]));
          LLtoXYZ(4.5*YAC_RAD, 0.5*YAC_RAD, &(search_coords[12][0]));
          LLtoXYZ(4.5*YAC_RAD, 1.5*YAC_RAD, &(search_coords[13][0]));
          LLtoXYZ(4.5*YAC_RAD, 2.5*YAC_RAD, &(search_coords[14][0]));

          // test yac_interp_grid_do_points_search
          size_t src_cells[15];
          yac_interp_grid_do_points_search(
            interp_grid, search_coords, 15, src_cells);

          struct grid_cell cell;
          yac_init_grid_cell(&cell);

          struct const_basic_grid_data * src_grid_data =
            yac_interp_grid_get_basic_grid_data_src(interp_grid);

          // check search results
          for (size_t i = 0; i < 15; ++i) {

            size_t curr_src_cell = src_cells[i];
            if (curr_src_cell != SIZE_MAX) {
              yac_const_basic_grid_data_get_grid_cell(
                src_grid_data, curr_src_cell, &cell);
              if (!yac_point_in_cell(search_coords[i], cell))
                PUT_ERR("error in yac_interp_grid_do_points_search");
            }
          }

          yac_free_grid_cell(&cell);
        }

        {
          double a[3], b[3];
          struct bounding_circle bnd_circles[5];
          size_t count = sizeof(bnd_circles) / sizeof(bnd_circles[0]);

          LLtoXYZ_deg(0.0, 0.0, a);
          LLtoXYZ_deg(0.4, 0.4, b);
          LLtoXYZ_deg(0.0, 0.0, bnd_circles[0].base_vector);
          bnd_circles[0].inc_angle = get_vector_angle_2(a, b);
          bnd_circles[0].sq_crd = DBL_MAX;

          LLtoXYZ_deg(0.0, 0.0, a);
          LLtoXYZ_deg(0.0, 0.9, b);
          LLtoXYZ_deg(0.0, 0.0, bnd_circles[1].base_vector);
          bnd_circles[1].inc_angle = get_vector_angle_2(a, b);
          bnd_circles[1].sq_crd = DBL_MAX;

          LLtoXYZ_deg(-1.0, -1.0, a);
          LLtoXYZ_deg(-0.9, -0.0, b);
          LLtoXYZ_deg(-1.0, -1.0, bnd_circles[2].base_vector);
          bnd_circles[2].inc_angle = get_vector_angle_2(a, b);
          bnd_circles[2].sq_crd = DBL_MAX;

          LLtoXYZ_deg(2.0, 1.5, a);
          LLtoXYZ_deg(4.0, 3.0, b);
          LLtoXYZ_deg(2.0, 1.5, bnd_circles[3].base_vector);
          bnd_circles[3].inc_angle = get_vector_angle_2(a, b);
          bnd_circles[3].sq_crd = DBL_MAX;

          LLtoXYZ_deg(2.0, 1.5, a);
          LLtoXYZ_deg(2.1, 1.5, b);
          LLtoXYZ_deg(2.0, 1.5, bnd_circles[4].base_vector);
          bnd_circles[4].inc_angle = get_vector_angle_2(a, b);
          bnd_circles[4].sq_crd = DBL_MAX;

          size_t ref_num_src_per_bnd_circle[] = {1, 3, 0, 12, 2};
          yac_int ref_src_cell_global_ids[] =
            {0, 0,1,4,   0,1,2,3,4,5,6,7,8,9,10,11, 5,6};

          size_t * cells, num_cells_per_bnd_circle[count];

          yac_interp_grid_do_bnd_circle_search_src(
            interp_grid, bnd_circles, count, 0,
            &cells, num_cells_per_bnd_circle);

          for (size_t i = 0, displ = 0; i < count; ++i) {
            if (num_cells_per_bnd_circle[i] !=  ref_num_src_per_bnd_circle[i]) {
              PUT_ERR("error in yac_interp_grid_do_bnd_circle_search_src");
            } else {
              yac_int src_cell_global_ids[20];
              yac_interp_grid_get_src_global_ids(
                interp_grid, cells + displ, num_cells_per_bnd_circle[i], 0,
                src_cell_global_ids);
              if (compare_global_ids(
                    src_cell_global_ids, ref_src_cell_global_ids + displ,
                    num_cells_per_bnd_circle[i]))
                PUT_ERR("error in yac_interp_grid_do_bnd_circle_search_src");
            }
            displ += num_cells_per_bnd_circle[i];
          }
          free(cells);
        }

        free(tgt_indices);
        yac_interp_grid_delete(interp_grid);
        yac_dist_grid_pair_delete(grid_pair);
        yac_basic_grid_delete(grids[1]);
        yac_basic_grid_delete(grids[0]);
      }

      { // test yac_interp_grid_relocate_src_tgt_pairs_orig

        for (int to_tgt_owner = 0; to_tgt_owner <= 1; ++to_tgt_owner) {

          int is_tgt = comm_rank >= 2;

          double coordinates_x[2][5] = {{0.0,1.0,2.0,3.0,4.0},
                                        {0.5,1.5,2.5,3.5,-1.0}};
          double coordinates_y[2][4] = {{0.0,1.0,2.0,3.0},
                                        {0.5,1.5,2.5,-1.0}};
          size_t num_cells[2][2] = {{4,3},{3,2}};
          size_t local_start[4][2] = {{0,0},{2,0},{0,0},{2,0}};
          size_t local_count[4][2] = {{2,3},{2,3},{2,2},{1,2}};
          int with_halo = 0;
          for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 5; ++j) coordinates_x[i][j] *= YAC_RAD;
            for (int j = 0; j < 4; ++j) coordinates_y[i][j] *= YAC_RAD;
          }

          struct basic_grid_data grid_data =
            yac_generate_basic_grid_data_reg2d(
              coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
              local_start[comm_rank], local_count[comm_rank], with_halo);

          char const * grid_names[2] = {"src_grid", "tgt_grid"};
          struct yac_basic_grid * grids[] =
            {yac_basic_grid_new(grid_names[is_tgt], grid_data),
            yac_basic_grid_empty_new(grid_names[is_tgt^1])};

          if (!is_tgt) {

            coordinate_pointer src_cell_coordinates =
              malloc(grid_data.num_cells * sizeof(*src_cell_coordinates));
            for (size_t i = 0; i < grid_data.num_cells; ++i) {
              for (int j = 0; j < 3; ++j) src_cell_coordinates[i][j] = 0.0;
              size_t * cell_vertices =
                grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
              for (int j = 0; j < grid_data.num_vertices_per_cell[i]; ++j)
                for (int k = 0; k < 3; ++k)
                  src_cell_coordinates[i][k] +=
                    grid_data.vertex_coordinates[cell_vertices[j]][k];
              normalise_vector(src_cell_coordinates[i]);
            }

            yac_basic_grid_add_coordinates(
              grids[0], CELL, src_cell_coordinates, grid_data.num_cells);
            free(src_cell_coordinates);
          }

          struct dist_grid_pair * grid_pair =
            yac_dist_grid_pair_new(grids[0], grids[1], comm);

          struct interp_field src_fields[] =
            {{.location = CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX}};
          size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
          struct interp_field tgt_field =
            {.location = CORNER, .coordinates_idx = SIZE_MAX,
             .masks_idx = SIZE_MAX};

          struct interp_grid * interp_grid =
            yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                                num_src_fields, src_fields, tgt_field);

          // get local tgt indices
          size_t * tgt_points, tgt_count;
          yac_interp_grid_get_tgt_points(interp_grid, &tgt_points, &tgt_count);

          // get local global tgt ids
          yac_int tgt_global_ids[12];
          yac_interp_grid_get_tgt_global_ids(
            interp_grid, tgt_points, tgt_count, tgt_global_ids);

          // get local tgt coordinates
          coordinate_pointer tgt_coordinates =
            malloc(tgt_count * sizeof(*tgt_coordinates));
          yac_interp_grid_get_tgt_coordinates(
            interp_grid, tgt_points, tgt_count, tgt_coordinates);

          // get matching source cells for each local tgt point
          size_t * src_points = malloc(tgt_count * sizeof(*src_points));
          yac_interp_grid_do_nnn_search_src(
            interp_grid, tgt_coordinates, tgt_count, 1, src_points);
          free(tgt_coordinates);

          // generate dummy weights
          double * weights = malloc(tgt_count * sizeof(*weights));
          for (size_t i = 0; i < tgt_count; ++i)
            weights[i] = (double)tgt_global_ids[i];

          // relocate target and matching source points to original owners
          yac_interp_grid_relocate_src_tgt_pairs_orig(
            interp_grid, to_tgt_owner, CELL, &src_points,
            &tgt_points, &weights, &tgt_count);

          // sort src points, tgt points and weights
          yac_int src_global_ids[12];
          yac_interp_grid_get_src_global_ids(
            interp_grid, src_points, tgt_count, 0, src_global_ids);
          yac_interp_grid_get_tgt_global_ids(
            interp_grid, tgt_points, tgt_count, tgt_global_ids);
          yac_quicksort_index_yac_int_yac_int_double(
            tgt_global_ids, tgt_count, src_global_ids, weights);

          // check results
          size_t ref_num_points[2][4] = {{6,6,0,0},{0,0,9,3}};
          yac_int ref_global_ids[2][4][9] =
            {{{0,1,4,5,8,9},{2,3,6,7,10,11},{-1},{-1}},
             {{-1},{-1},{0,1,2,4,5,6,8,9,10},{3,7,11}}};
          if (tgt_count != ref_num_points[to_tgt_owner][comm_rank])
            PUT_ERR(
              "error in yac_interp_grid_relocate_src_tgt_pairs_orig");
          for (size_t i = 0; i < tgt_count; ++i) {
            if (src_global_ids[i] != ref_global_ids[to_tgt_owner][comm_rank][i])
              PUT_ERR(
                "error in yac_interp_grid_relocate_src_tgt_pairs_orig");
            if (tgt_global_ids[i] != ref_global_ids[to_tgt_owner][comm_rank][i])
              PUT_ERR(
                "error in yac_interp_grid_relocate_src_tgt_pairs_orig");
            if (weights[i] != (double)ref_global_ids[to_tgt_owner][comm_rank][i])
              PUT_ERR(
                "error in yac_interp_grid_relocate_src_tgt_pairs_orig");
          }

          free(tgt_points);
          free(weights);
          free(src_points);
          yac_interp_grid_delete(interp_grid);
          yac_dist_grid_pair_delete(grid_pair);
          yac_basic_grid_delete(grids[1]);
          yac_basic_grid_delete(grids[0]);
        }
      }
    }

    yac_mpi_call(MPI_Comm_free(&comm), MPI_COMM_WORLD);

  } else {
    PUT_ERR("insufficient number of processes");
  }

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static int check_global_ids(
  yac_int * global_ids, int count, int ref_global_count, MPI_Comm comm) {

  int error_count = 0;

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int * counts = xmalloc((size_t)comm_size * sizeof(*counts));
  int * displs = xmalloc((size_t)comm_size * sizeof(*displs));
  yac_mpi_call(
    MPI_Allgather(&count, 1, MPI_INT, counts, 1, MPI_INT, comm), comm);
  for (int i = 0, accu = 0; i < comm_size; accu += counts[i++])
    displs[i] = accu;

  int global_count = 0;
  for (int i = 0; i < comm_size; ++i) global_count += counts[i];

  yac_int * all_global_ids =
    (comm_rank == 0)?
      xmalloc((size_t)global_count * sizeof(*all_global_ids)):
      NULL;
  yac_mpi_call(
    MPI_Gatherv(global_ids, count, yac_int_dt, all_global_ids,
                counts, displs, yac_int_dt, 0, comm), comm);

  if (comm_rank == 0) {
    int * flag = xcalloc((size_t)ref_global_count, sizeof(*flag));
    int flag_count = 0;
    for (int i = 0; i < global_count; ++i) {
      if ((all_global_ids[i] >= ref_global_count) ||
          (all_global_ids[i] < 0)) {
        ++error_count;
      } else {
        if (!flag[all_global_ids[i]]) ++flag_count;
        flag[all_global_ids[i]] = 1;
      }
    }
    if (flag_count != ref_global_count) ++error_count;
    free(flag);
  }

  free(all_global_ids);
  free(displs);
  free(counts);

  return error_count;
}

static int compare_global_ids(
  yac_int * global_ids, yac_int * ref_global_ids, size_t count) {

  size_t match_count = 0;
  for (size_t i = 0; i < count; ++i)
    for (size_t j = 0; j < count; ++j)
      if (global_ids[i] == ref_global_ids[j]) ++match_count;
  return match_count != count;
}
