/**
 * @file test_dynamic_config_c.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Nils-Arne Dreier <dreier@dkrz.de>
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

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "tests.h"
#include "yac_interface.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

static void gen_weight_file(void);
static void comp1_main(void);
static void comp2_main(void);
static int gen_field_avg(int yac_id);
static int gen_field_nnn(int yac_id);
static int gen_field_conserv(int yac_id);
static int gen_field_spmap(int yac_id);
static int gen_field_hcsbb(int yac_id);
static int gen_field_user_file(int yac_id);
static int gen_field_fixed(int yac_id);
static int gen_field_creep(int yac_id);
static int gen_field_check(int yac_id);
static int gen_field_callback(int yac_id);

int n_cells[2] = {360, 180};
int n_corners[2] = {361, 181};
int cyclic[2] = {0, 0};
double * data, * corner_lon, * corner_lat, * cell_lon, * cell_lat;
int timestep_counter;
char const * weight_file_name = "test_dynamic_config_c.nc";
char const * weight_file_name_dummy = "test_dynamic_config_dummy_c.nc";
int comp_id, point_id;

void compute_weights_callback(
  double const tgt_coords[3],
  int src_cell_id, size_t src_cell_idx,
  int const **global_results_points, double **result_weights,
  size_t *result_count, void *user_data) {
  (void)tgt_coords;
  (void)src_cell_id;
  (void)src_cell_idx;
  *global_results_points = NULL;
  *result_weights = NULL;
  *result_count = 0;
  (void)user_data;
}

int main (void) {

  MPI_Init(NULL, NULL);

  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);

  if (global_size != 4) {
    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }


  data = malloc((size_t)(n_cells[0] * n_cells[1]) * sizeof(*data));
  corner_lon = malloc((size_t)(n_corners[0]) * sizeof(*corner_lon));
  corner_lat = malloc((size_t)(n_corners[1]) * sizeof(*corner_lat));
  cell_lon = malloc((size_t)(n_cells[0]) * sizeof(*cell_lon));
  cell_lat = malloc((size_t)(n_cells[1]) * sizeof(*cell_lat));

  double const delta_lon = 1.0, delta_lat=1.0;
  for (int i = 0; i < n_corners[0]; ++i)
    corner_lon[i] = ((double)i * delta_lon - 180.0) * YAC_RAD;
  for (int i = 0; i < n_corners[1]; ++i)
    corner_lat[i] = ((double)i * delta_lat - 90.0) * YAC_RAD;
  for (int i = 0; i < n_cells[0]; ++i)
    cell_lon[i] = (((double)i + 0.5) * delta_lon - 180.0) * YAC_RAD;
  for (int i = 0; i < n_cells[1]; ++i)
    cell_lat[i] = (((double)i + 0.5) * delta_lat - 90.0) * YAC_RAD;

  yac_cadd_compute_weights_callback(
    compute_weights_callback, NULL, "compute_weights_callback");
  if (global_rank == 0) gen_weight_file();

  timestep_counter = 0;

  if ((global_rank % 2) == 0) comp1_main();
  else                        comp2_main();

  if (global_rank == 0) {
    unlink(weight_file_name);
    unlink(weight_file_name_dummy);
  }

  free(cell_lat);
  free(cell_lon);
  free(corner_lat);
  free(corner_lon);
  free(data);

  if (timestep_counter != 11) PUT_ERR("wrong final timestep_counter");

  MPI_Finalize();
  return TEST_EXIT_CODE;
}

static void comp1_main(void) {

  int grid_id;

  yac_cinit();
  // this component defines the datetime
  yac_cdef_datetime("2000-01-01T00:00:00", "2000-01-01T00:00:10");
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  // test yac_cpredef_comp:
  yac_cpredef_comp("comp1", &comp_id);
  yac_cdef_comps(NULL, 0, NULL);
  yac_cdef_grid_reg2d(
    "grid1", n_corners, cyclic, corner_lon, corner_lat, &grid_id);
  yac_cdef_points_reg2d(
    grid_id, n_cells, YAC_LOCATION_CELL, cell_lon, cell_lat, &point_id);

  int field_1, field_2, field_3, field_4, field_out;
  yac_cdef_field(
    "A", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_1);
  yac_cdef_field(
    "B", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_2);
  yac_cdef_field(
    "C", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_3);
  yac_cdef_field(
    "D", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_4);

  yac_cdef_field(
    "OUT", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_out);

  int avg_interp;
  yac_cget_interp_stack_config(&avg_interp);
  yac_cadd_interp_stack_config_average(avg_interp, YAC_AVG_ARITHMETIC, 1);
  yac_cdef_couple(
    "comp2", "grid2", "4", "comp1", "grid1", "D",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_ACCUMULATE, avg_interp,
    0, 0);
  yac_cfree_interp_stack_config(avg_interp);

  yac_cenddef();

  int info = YAC_ACTION_NONE;
  while (info != YAC_ACTION_PUT_FOR_RESTART) {
    int ierror;
    timestep_counter++;
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = 42.0;
    yac_cput_(field_1, 1, data, &info, &ierror);
    printf("time of field_1: %s\n", yac_cget_field_datetime(field_1));
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = 43.0;
    yac_cput_(field_2, 1, data, &info, &ierror);
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = 44.0;
    yac_cput_(field_3, 1, data, &info, &ierror);
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = -1.0;
    yac_cget(field_4, 1, &data, &info, &ierror);
    if (fabs(data[0] - 40.0) > 1.0e-14) PUT_ERR("error in data for field_4");
    yac_cput_(field_out, 1, data, &info, &ierror);
  }

  yac_cfinalize();
}

static void comp2_main(void) {

  int yac_id, grid_id, avg_interp;
  yac_cinit_instance(&yac_id);
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  // test yac_cpredef_comp_instance
  yac_cpredef_comp_instance(yac_id, "comp2", &comp_id);
  yac_cdef_comps_instance(yac_id, NULL, 0, NULL);

  yac_csync_def_instance(yac_id);

  // define fields after sync_def to tests the sync mechanism
  yac_cdef_grid_reg2d(
    "grid2", n_corners, cyclic, corner_lon, corner_lat, &grid_id);
  yac_cdef_points_reg2d(
    grid_id, n_cells, YAC_LOCATION_CELL, cell_lon, cell_lat, &point_id);

  int field_1, field_2, field_3, field_4;
  yac_cdef_field(
    "1", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_1);
  yac_cdef_field(
    "2", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_2);
  yac_cdef_field(
    "3", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_3);
  yac_cdef_field(
    "4", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_4);

  yac_cget_interp_stack_config(&avg_interp);
  yac_cadd_interp_stack_config_average(avg_interp, YAC_AVG_ARITHMETIC, 1);

  yac_cdef_couple_instance(
    yac_id, "comp1", "grid1", "A", "comp2", "grid2", "1",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_ACCUMULATE, avg_interp,
    0, 0);
  yac_cdef_couple_instance(
    yac_id, "comp1", "grid1", "B", "comp2", "grid2", "2",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_ACCUMULATE, avg_interp,
    0, 0);
  yac_cdef_couple_instance(
    yac_id, "comp1", "grid1", "C", "comp2", "grid2", "3",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_ACCUMULATE, avg_interp,
    0, 0);

  yac_cfree_interp_stack_config(avg_interp);

  // check all interpolation methods
  int field_in[] =
    {gen_field_avg(yac_id),
     gen_field_nnn(yac_id),
     gen_field_conserv(yac_id),
     gen_field_spmap(yac_id),
     gen_field_hcsbb(yac_id),
     gen_field_user_file(yac_id),
     gen_field_fixed(yac_id),
     gen_field_creep(yac_id),
     gen_field_check(yac_id),
     gen_field_callback(yac_id)};

  yac_cenddef_instance(yac_id);

  int info = YAC_ACTION_NONE;
  while(info != YAC_ACTION_GET_FOR_RESTART) {
    int ierror;
    timestep_counter++;
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = -1.0;
    yac_cget(field_1, 1, &data, &info, &ierror);
    printf("time of field_1: %s\n", yac_cget_field_datetime(field_1));
    if (fabs(data[0] - 42.0) > 1.0e-14) PUT_ERR("error in data for field_1");
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = -1.0;
    yac_cget(field_2, 1, &data, &info, &ierror);
    if (fabs(data[0] - 43.0) > 1.0e-14) PUT_ERR("error in data for field_2");
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = -1.0;
    yac_cget(field_3, 1, &data, &info, &ierror);
    if (fabs(data[0] - 44.0) > 1.0e-14) PUT_ERR("error in data for field_3");
    for (int i = 0; i < n_cells[0] * n_cells[1]; ++i) data[i] = 40.0;
    yac_cput_(field_4, 1, data, &info, &ierror);
    for (int i = 0; i < sizeof(field_in) / sizeof(field_in[0]); ++i)
      yac_cget(field_in[i], 1, &data, &info, &ierror);
  }

  yac_cfinalize_instance(yac_id);
}

static int gen_field(int yac_id, char const * field_name, int interp_stack) {
  int field_id;
  yac_cdef_field(
    field_name, comp_id, &point_id, 1, 1, "1",
    YAC_TIME_UNIT_SECOND, &field_id);
  yac_cdef_couple_instance(
    yac_id, "comp1", "grid1", "OUT", "comp2", "grid2", field_name,
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, interp_stack,
    0, 0);
  yac_cfree_interp_stack_config(interp_stack);
  return field_id;
}

static int gen_field_avg(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_average(
    interp_stack, YAC_AVG_ARITHMETIC, 1);
  return gen_field(yac_id, "avg", interp_stack);
}

static int gen_field_nnn(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_nnn(
    interp_stack, YAC_NNN_AVG, 3, 1.0);
  return gen_field(yac_id, "nn", interp_stack);
}

static int gen_field_conserv(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_conservative(
    interp_stack, 1, 1, 1, YAC_CONSERV_DESTAREA);
  return gen_field(yac_id, "conserv", interp_stack);
}

static int gen_field_spmap(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_spmap(
    interp_stack, 0.0, 0.0, YAC_SPMAP_AVG);
  return gen_field(yac_id, "spmap", interp_stack);
}

static int gen_field_hcsbb(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_hcsbb(interp_stack);
  return gen_field(yac_id, "hcsbb", interp_stack);
}

static int gen_field_user_file(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_user_file(
      interp_stack, weight_file_name, "grid1", "grid2");
  return gen_field(yac_id, "user_file", interp_stack);
}

static int gen_field_fixed(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_fixed(interp_stack, -1.0);
  return gen_field(yac_id, "fixed", interp_stack);
}

static int gen_field_creep(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_creep(interp_stack, 0);
  return gen_field(yac_id, "creep", interp_stack);
}

static int gen_field_check(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_check(interp_stack, NULL, NULL);
  return gen_field(yac_id, "check", interp_stack);
}

static int gen_field_callback(int yac_id) {
  int interp_stack;
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_user_callback(
    interp_stack, "compute_weights_callback");
  return gen_field(yac_id, "callback", interp_stack);
}

static void gen_weight_file(void) {

  int yac_id, comp_id, grid_id_out, grid_id_in, point_id_out, point_id_in,
      field_id_out, field_id_in, interp_stack, ext_couple_config;
  int mapping_side = 0;
  double scale_factor = 10.0;
  double scale_summand = 0.5;
  yac_cinit_comm_instance(MPI_COMM_SELF, &yac_id);
  yac_cdef_datetime_instance(yac_id, "2000-01-01T00:00:00", "2000-01-01T00:00:10");
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_comp_instance(yac_id, "comp", &comp_id);
  yac_cdef_grid_reg2d(
    "grid1", n_corners, cyclic, corner_lon, corner_lat, &grid_id_out);
  yac_cdef_grid_reg2d(
    "grid2", n_corners, cyclic, corner_lon, corner_lat, &grid_id_in);
  yac_cdef_points_reg2d(
    grid_id_out, n_cells, YAC_LOCATION_CELL, cell_lon, cell_lat, &point_id_out);
  yac_cdef_points_reg2d(
    grid_id_in, n_cells, YAC_LOCATION_CELL, cell_lon, cell_lat, &point_id_in);
  yac_cdef_field(
    "out", comp_id, &point_id_out, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_id_out);
  yac_cdef_field(
    "in", comp_id, &point_id_in, 1, 1, "1", YAC_TIME_UNIT_SECOND, &field_id_in);
  yac_cget_interp_stack_config(&interp_stack);
  yac_cadd_interp_stack_config_average(interp_stack, YAC_AVG_ARITHMETIC, 1);
  yac_cget_ext_couple_config(&ext_couple_config);
  yac_cset_ext_couple_config_weight_file(
    ext_couple_config, weight_file_name);
  yac_cset_ext_couple_config_mapping_side(
    ext_couple_config, mapping_side);
  yac_cset_ext_couple_config_scale_factor(
    ext_couple_config, scale_factor);
  yac_cset_ext_couple_config_scale_summand(
    ext_couple_config, scale_summand);
  yac_cdef_couple_custom_instance(
    yac_id, "comp", "grid1", "out", "comp", "grid2", "in",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, interp_stack,
    0, 0, ext_couple_config);
  char const * temp_weight_file_name;
  yac_cget_ext_couple_config_weight_file(
    ext_couple_config, &temp_weight_file_name);
  if (strcmp(weight_file_name, temp_weight_file_name))
    PUT_ERR("error in yac_cget_ext_couple_config_weight_file");
  int temp_mapping_side;
  yac_cget_ext_couple_config_mapping_side(
    ext_couple_config, &temp_mapping_side);
  if (mapping_side != temp_mapping_side)
    PUT_ERR("error in yac_cget_ext_couple_config_mapping_side");
  double temp_scale_factor;
  yac_cget_ext_couple_config_scale_factor(
    ext_couple_config, &temp_scale_factor);
  if (scale_factor != temp_scale_factor)
    PUT_ERR("error in yac_cget_ext_couple_config_scale_factor");
  double temp_scale_summand;
  yac_cget_ext_couple_config_scale_summand(
    ext_couple_config, &temp_scale_summand);
  if (scale_summand != temp_scale_summand)
    PUT_ERR("error in yac_cget_ext_couple_config_scale_summand");
  yac_cfree_ext_couple_config(ext_couple_config);
  yac_cfree_interp_stack_config(interp_stack);
  yac_cenddef_instance(yac_id);
  yac_cfinalize_instance(yac_id);
}
