// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "interp_method.h"
#include "interp_method_callback.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

struct grid_info {
  size_t * cell_to_vertex;
  yac_int * vertex_ids;
  yac_int * cell_ids;
  size_t num_grid_cells;
};

static void test1_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);
static void test2_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);
static void test3_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);
static void test4_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);
static char const * grid_names[2] = {"src_grid", "tgt_grid"};

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 3) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  { // test setting and getting callbacks

    struct {
      yac_func_compute_weights func;
      char const * key;
      int user_data;
    } callbacks[] =
      {{.func = test1_callback, .key = "test1", .user_data = 1},
       {.func = test2_callback, .key = "test2", .user_data = 2},
       {.func = test3_callback, .key = "test3", .user_data = 3}};
    enum {NUM_CALLBACKS = sizeof(callbacks) / sizeof(callbacks[0])};

    for (size_t i = 0; i < NUM_CALLBACKS; ++i)
      yac_interp_method_callback_add_compute_weights_callback(
        callbacks[i].func, (void*)&(callbacks[i].user_data), callbacks[i].key);
    // add the second key again
    yac_interp_method_callback_add_compute_weights_callback(
      callbacks[1].func, (void*)&(callbacks[1].user_data), callbacks[1].key);

    for (size_t i = 0; i < NUM_CALLBACKS; ++i) {

      yac_func_compute_weights func;
      void * user_data;

      yac_interp_method_callback_get_compute_weights_callback(
        callbacks[i].key, &func, &user_data);

      if (func != callbacks[i].func)
        PUT_ERR(
          "ERROR in yac_interp_method_callback_get_compute_weights_callback");
      if (user_data != (void*)&(callbacks[i].user_data))
        PUT_ERR(
          "ERROR in yac_interp_method_callback_get_compute_weights_callback");
    }
  }

  { // test with two source processes, one target process and two source fields
    // the global source grid is a 4x2 grid:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    //---------------
    // setup
    //---------------

    int is_tgt = comm_rank == 2;
    double coordinates_x[2][5] = {{0.0,1.0,2.0,3.0}, {-0.5,0.5,1.5,2.5,3.5}};
    double coordinates_y[2][4] = {{0.0,1.0,2.0}, {-0.5,0.5,1.5,2.5}};
    size_t const num_cells[2][2] = {{3,2}, {4,3}};
    size_t local_start[2][2][2] = {{{0,0},{1,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{1,2},{2,2}}, {{4,3}}};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[is_tgt][i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[is_tgt][i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
        local_start[is_tgt][(is_tgt)?0:comm_rank],
        local_count[is_tgt][(is_tgt)?0:comm_rank], with_halo);

    yac_int * cell_ids =
      xmalloc(grid_data.num_cells * sizeof(*cell_ids));
    memcpy(cell_ids, grid_data.cell_ids,
           grid_data.num_cells * sizeof(*cell_ids));
    size_t * cell_to_vertex =
      xmalloc(4 * grid_data.num_cells * sizeof(*cell_to_vertex));
    memcpy(cell_to_vertex, grid_data.cell_to_vertex,
           4 * grid_data.num_cells * sizeof(*cell_to_vertex));
    yac_int * vertex_ids =
      xmalloc(grid_data.num_vertices * sizeof(*vertex_ids));
    memcpy(vertex_ids, grid_data.vertex_ids,
           grid_data.num_vertices * sizeof(*vertex_ids));

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX},
       {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct grid_info grid_info = {
      .cell_to_vertex = cell_to_vertex,
      .vertex_ids = vertex_ids,
      .cell_ids = cell_ids,
      .num_grid_cells = grid_data.num_cells};
    yac_func_compute_weights test_functions[] = {test1_callback,
                                                 test2_callback,
                                                 test4_callback};
    enum {
      NUM_TEST_FUNCTIONS = sizeof(test_functions) / sizeof(test_functions[0])};
    void *user_data[NUM_TEST_FUNCTIONS] =
      {(void*)&num_src_fields, (void*)&grid_info, (void*)&grid_info};
    double ref_tgt_data[NUM_TEST_FUNCTIONS][20] =
      {{-1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, -1.0, -1.0, -1.0, -1.0},
       {-1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, 0.5*0+0.125*(0+1+4+5), 0.5*1+0.125*(1+2+5+6), 0.5*2+0.125*(2+3+6+7), -1.0,
        -1.0, 0.5*3+0.125*(4+5+8+9), 0.5*4+0.125*(5+6+9+10), 0.5*5+0.125*(6+7+10+11), -1.0,
        -1.0, -1.0, -1.0, -1.0, -1.0},
       {-1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, -1.0, 0.5*1+0.125*(1+2+5+6), -1.0, -1.0,
        -1.0, 0.5*3+0.125*(4+5+8+9), -1.0, 0.5*5+0.125*(6+7+10+11), -1.0,
        -1.0, -1.0, -1.0, -1.0, -1.0}};
    double collection_factor[] = {0.0, 10.0, 10.0};

    for (size_t test_idx = 0; test_idx < NUM_TEST_FUNCTIONS; ++test_idx) {

      struct interp_method * method_stack[] =
        {yac_interp_method_callback_new(
           test_functions[test_idx], user_data[test_idx]),
         yac_interp_method_fixed_new(-1.0), NULL};

      struct yac_interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);

      enum yac_interp_weights_reorder_type reorder_type[2] =
        {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

      for (size_t collection_size = 1; collection_size <= 16;
           collection_size *= 2) {

        for (size_t reorder_idx = 0; reorder_idx < 2; ++reorder_idx) {

          struct yac_interpolation * interpolation =
            yac_interp_weights_get_interpolation(
              weights, reorder_type[reorder_idx], collection_size,
              YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

          // check generated interpolation
          {
            double *** src_data = NULL;
            double ** tgt_data = NULL;

            if (is_tgt) {
              tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
              for (size_t collection_idx = 0; collection_idx < collection_size;
                   ++collection_idx) {
                tgt_data[collection_idx] =
                  xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
                for (size_t i = 0; i < grid_data.num_vertices; ++i)
                  tgt_data[collection_idx][i] = -1337.0;
              }
            } else {
              src_data = xmalloc(collection_size * sizeof(*src_data));
              for (size_t collection_idx = 0; collection_idx < collection_size;
                   ++collection_idx) {
                src_data[collection_idx] = xmalloc(2 * sizeof(**src_data));
                src_data[collection_idx][0] =
                  xmalloc(grid_data.num_vertices * sizeof(***src_data));
                for (size_t i = 0; i < grid_data.num_vertices; ++i)
                  src_data[collection_idx][0][i] =
                    (double)(vertex_ids[i] + 10 * collection_idx);
                src_data[collection_idx][1] =
                  xmalloc(grid_data.num_cells * sizeof(***src_data));
                for (size_t i = 0; i < grid_data.num_cells; ++i)
                  src_data[collection_idx][1][i] =
                    (double)(cell_ids[i] + 10 * collection_idx);
              }
            }

            yac_interpolation_execute(interpolation, src_data, tgt_data);

            if (is_tgt) {

              // check results
              for (size_t collection_idx = 0; collection_idx < collection_size;
                   ++collection_idx) {
                for (size_t i = 0; i < grid_data.num_vertices; ++i)
                  if (ref_tgt_data[test_idx][vertex_ids[i]] == -1.0) {
                    if (tgt_data[collection_idx][i] != -1.0)
                      PUT_ERR("wrong results");
                  } else {
                    if (fabs(tgt_data[collection_idx][i] -
                             (ref_tgt_data[test_idx][vertex_ids[i]] +
                              (double)collection_idx *
                              collection_factor[test_idx])) > 1e-9)
                    PUT_ERR("wrong results");
                  }
              }

              for (size_t collection_idx = 0; collection_idx < collection_size;
                   ++collection_idx) free(tgt_data[collection_idx]);
              free(tgt_data);
            } else {
              for (size_t collection_idx = 0; collection_idx < collection_size;
                   ++collection_idx) {
                for (size_t src_field_idx = 0; src_field_idx < num_src_fields;
                     ++src_field_idx)
                  free(src_data[collection_idx][src_field_idx]);
                free(src_data[collection_idx]);
              }
              free(src_data);
            }
          }
          yac_interpolation_delete(interpolation);
        }
      }

      yac_interp_weights_delete(weights);
      yac_interp_method_delete(method_stack);
    }
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
    free(cell_to_vertex);
    free(vertex_ids);
    free(cell_ids);
  }

  {
    // All three ranks are source processes. Rank 0 is also a target process.
    //---------------
    // setup
    //---------------

    double src_coordinates_x[51];
    double src_coordinates_y[51];
    size_t const src_num_cells[2] = {50,50};
    size_t src_local_start[3][2] = {{20,10},{20,0},{0,0}};
    size_t src_local_count[3][2] = {{30,40},{30,10},{20,50}};
    int with_halo = 0;
    for (size_t i = 0; i <= src_num_cells[0]; ++i)
      src_coordinates_x[i] = (double)i * YAC_RAD;
    for (size_t i = 0; i <= src_num_cells[1]; ++i)
      src_coordinates_y[i] = (double)i * YAC_RAD;

    struct yac_basic_grid_data src_grid_data =
      yac_generate_basic_grid_data_reg2d(
        src_coordinates_x, src_coordinates_y, src_num_cells,
        src_local_start[comm_rank], src_local_count[comm_rank], with_halo);

    double tgt_coordinates_x[50];
    double tgt_coordinates_y[50];
    size_t const tgt_num_cells[2] = {49,49};
    size_t tgt_local_start[2] = {0,0};
    size_t tgt_local_count[2] = {49,49};
    for (size_t i = 0; i <= tgt_num_cells[0]; ++i)
      tgt_coordinates_x[i] = (0.5 + (double)i) * YAC_RAD;
    for (size_t i = 0; i <= tgt_num_cells[1]; ++i)
      tgt_coordinates_y[i] = (0.5 + (double)i) * YAC_RAD;

    struct yac_basic_grid_data tgt_grid_data;
    if (comm_rank == 0)
      tgt_grid_data =
        yac_generate_basic_grid_data_reg2d(
          tgt_coordinates_x, tgt_coordinates_y, tgt_num_cells,
          tgt_local_start, tgt_local_count, with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(grid_names[0], src_grid_data);
    struct yac_basic_grid * tgt_grid = NULL;
    if (comm_rank == 0)
      tgt_grid = yac_basic_grid_new(grid_names[1], tgt_grid_data);
    else
      tgt_grid = yac_basic_grid_empty_new(grid_names[1]);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_callback_new(test3_callback, NULL),
       yac_interp_method_fixed_new(-1.0), NULL};

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, YAC_MAPPING_ON_SRC, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

    // check generated interpolation
    {
      double *** src_data = NULL;
      double ** tgt_data = NULL;

      if (comm_rank == 0) {
        tgt_data = xmalloc(1 * sizeof(*tgt_data));
        tgt_data[0] = xmalloc(tgt_grid_data.num_vertices * sizeof(**tgt_data));
        for (size_t i = 0; i < tgt_grid_data.num_vertices; ++i)
          tgt_data[0][i] = -1.0;
      }

      src_data = xmalloc(1 * sizeof(*src_data));
      src_data[0] = xmalloc(1 * sizeof(**src_data));
      src_data[0][0] =
        xmalloc(src_grid_data.num_cells * sizeof(***src_data));
      for (size_t i = 0; i < src_grid_data.num_cells; ++i)
        src_data[0][0][i] = (double)(src_grid_data.cell_ids[i]);

      yac_interpolation_execute(interpolation, src_data, tgt_data);

      if (comm_rank == 0) {

        // check results
        for (size_t i = 0; i < tgt_grid_data.num_vertices; ++i)
          if (tgt_data[0][i] != (double)(tgt_grid_data.vertex_ids[i]))
            PUT_ERR("wrong results");

        free(tgt_data[0]);
        free(tgt_data);
      }

      free(src_data[0][0]);
      free(src_data[0]);
      free(src_data);
    }
    yac_interpolation_delete(interpolation);
    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);
  }

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void test1_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  UNUSED(tgt_coords);
  UNUSED(src_cell_id);
  UNUSED(src_cell_idx);

  size_t num_src_fields = *(size_t*)user_data;

  for (size_t i = 0; i < num_src_fields; ++i) {
    global_results_points[i] = NULL;
    result_weights[i] = NULL;
    result_count[i] = 0;
  }
}

static void test2_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  UNUSED(tgt_coords);

  struct grid_info * info = (struct grid_info *)user_data;

  // consistency check
  if ((info->num_grid_cells <= src_cell_idx) ||
      (info->cell_ids[src_cell_idx] != src_cell_id)) {
    fputs("ERROR(test2_callback): inconsistent data\n", stderr);
    exit(EXIT_FAILURE);
  }

  static double vertex_weights[4] = {0.125, 0.125, 0.125, 0.125};
  static double cell_weight[1] = {0.5};

  static int vertex_result_points[4];
  static int cell_result_points[1];

  for (size_t i = 0; i < 4; ++i)
    vertex_result_points[i] =
      (int)(info->vertex_ids[info->cell_to_vertex[4*src_cell_idx + i]]);
  cell_result_points[0] = (int)src_cell_id;

  global_results_points[0] = vertex_result_points;
  global_results_points[1] = cell_result_points;
  result_weights[0] = vertex_weights;
  result_weights[1] = cell_weight;
  result_count[0] = 4;
  result_count[1] = 1;
}

static void test3_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  UNUSED(tgt_coords);
  UNUSED(src_cell_id);
  UNUSED(src_cell_idx);
  UNUSED(user_data);

  static int result_points[1];
  static double weight[1] = {1.0};

  result_points[0] = (int)src_cell_id;

  global_results_points[0] = result_points;
  result_weights[0] = weight;
  result_count[0] = 1;
}

static void test4_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  UNUSED(tgt_coords);

  struct grid_info * info = (struct grid_info *)user_data;

  // consistency check
  if ((info->num_grid_cells <= src_cell_idx) ||
      (info->cell_ids[src_cell_idx] != src_cell_id)) {
    fputs("ERROR(test2_callback): inconsistent data\n", stderr);
    exit(EXIT_FAILURE);
  }

  static double vertex_weights[4] = {0.125, 0.125, 0.125, 0.125};
  static double cell_weight[1] = {0.5};

  static int vertex_result_points[4];
  static int cell_result_points[1];

  for (size_t i = 0; i < 4; ++i)
    vertex_result_points[i] =
      (int)(info->vertex_ids[info->cell_to_vertex[4*src_cell_idx + i]]);
  cell_result_points[0] = (int)src_cell_id;

  if (src_cell_id&1) {
    global_results_points[0] = vertex_result_points;
    global_results_points[1] = cell_result_points;
    result_weights[0] = vertex_weights;
    result_weights[1] = cell_weight;
    result_count[0] = 4;
    result_count[1] = 1;
  } else {
    global_results_points[0] = NULL;
    global_results_points[1] = NULL;
    result_weights[0] = NULL;
    result_weights[1] = NULL;
    result_count[0] = 0;
    result_count[1] = 0;
  }
}
