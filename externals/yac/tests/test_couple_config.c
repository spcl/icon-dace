/**
 * @file test_couple_config.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Teresa Holfeld <teresa.holfeld@mpimet.mpg.de>
 *         Moritz Hanke <hanke@dkrz.de>
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "tests.h"
#include "utils.h"
#include "config_yaml.h"
#include "instance.h"
#include "yac_mpi.h"
#include "event.h"

static void check_couple_config(struct yac_couple_config * couple_config);
static struct yac_couple_config * generate_couple_config_from_YAML_parallel(
  char const * config_filename);
static struct yac_couple_config * generate_couple_config_from_YAML(
  char const * config_filename, int parse_flags);
static void write_couple_config_to_YAML(
  struct yac_couple_config * couple_config, char const * config_filename,
  int emit_flags);

int main(void) {

  yac_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  int size, rank;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  if (size != 2) {
    fputs("ERROR wrong number of processes (has to be 2)", stderr);
    exit(EXIT_FAILURE);
  }

  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

  struct yac_couple_config * couple_config =
    generate_couple_config_from_YAML_parallel("couple_config_test.yaml");

  // check the emitting of coupling configurations (on rank zero only)
  if (rank == 0) {

    struct {
      int emit_flags;
      int parse_flags;
      char const * ext;
    } formats[] =
    {{.emit_flags = YAC_YAML_EMITTER_DEFAULT,
      .parse_flags = YAC_YAML_PARSER_DEFAULT,
      .ext = "yaml"},
     {.emit_flags =  YAC_YAML_EMITTER_JSON,
      .parse_flags = YAC_YAML_PARSER_JSON_FORCE,
      .ext = "json"}};
    enum {NUM_FORMATS = sizeof(formats) / sizeof(formats[0])};

    for (int i = 0; i < NUM_FORMATS; ++i) {
      char config_filename[32];
      snprintf(config_filename, sizeof(config_filename),
               "temp_couple_config_test.%s", formats[i].ext);
      write_couple_config_to_YAML(
        couple_config, config_filename, formats[i].emit_flags);

      check_couple_config(
        generate_couple_config_from_YAML(
          config_filename, formats[i].parse_flags));

      unlink(config_filename);
    }
  }

  check_couple_config(couple_config);

  { // reproduces a bug in routine dist_merge
    struct yac_couple_config * couple_config = yac_couple_config_new();
    struct yac_interp_stack_config * interp_stack_config =
      yac_interp_stack_config_new();
    yac_interp_stack_config_add_fixed(interp_stack_config, -1.0);
    yac_couple_config_def_couple(
      couple_config, "comp_a", "grid_a", "field_a",
      "comp_b", "grid_b", "field_b", "1", YAC_REDUCTION_TIME_NONE,
      interp_stack_config, 0, 0, NULL, 0, 1.0, 0.0, 0, NULL, NULL);
    yac_interp_stack_config_delete(interp_stack_config);
    yac_couple_config_add_component(
      couple_config, (rank == 0)?"comp_c":"comp_d");
    yac_couple_config_sync(couple_config, MPI_COMM_WORLD);
    yac_couple_config_delete(couple_config);
  }

  yac_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
}

static void check_couple_config(struct yac_couple_config * couple_config) {

  if (yac_couple_config_get_redirect_stdout(couple_config) != 1)
    PUT_ERR("ERROR in yac_couple_config_get_redirect_stdout\n");

  if (yac_couple_config_get_num_components(couple_config) != 4)
    PUT_ERR("ERROR in yac_couple_config_get_num_components\n");

  if (!yac_couple_config_component_name_is_valid(couple_config, "ICON-O"))
    PUT_ERR("ERROR in yac_couple_config_component_name_is_valid\n");
  if (!yac_couple_config_component_name_is_valid(couple_config, "ICON-A"))
    PUT_ERR("ERROR in yac_couple_config_component_name_is_valid\n");
  if (!yac_couple_config_component_name_is_valid(couple_config, "DUMMY"))
    PUT_ERR("ERROR in yac_couple_config_component_name_is_valid\n");
  if (!yac_couple_config_component_name_is_valid(couple_config, "DUMMY_2"))
    PUT_ERR("ERROR in yac_couple_config_component_name_is_valid\n");
  if (yac_couple_config_component_name_is_valid(couple_config, "INVALID"))
    PUT_ERR("ERROR in yac_couple_config_component_name_is_valid\n");

  if (!yac_couple_config_contains_grid_name(couple_config, "grid1"))
    PUT_ERR("ERROR in yac_couple_config_contains_grid_name\n");
  if (!yac_couple_config_contains_grid_name(couple_config, "grid2"))
    PUT_ERR("ERROR in yac_couple_config_contains_grid_name\n");
  if (!yac_couple_config_contains_grid_name(couple_config, "grid3"))
    PUT_ERR("ERROR in yac_couple_config_contains_grid_name\n");
  if (!yac_couple_config_contains_grid_name(couple_config, "grid4"))
    PUT_ERR("ERROR in yac_couple_config_contains_grid_name\n");
  if (yac_couple_config_contains_grid_name(couple_config, "grid5"))
    PUT_ERR("ERROR in yac_couple_config_contains_grid_name\n");

  if (strcmp("grid1_meta",
             yac_couple_config_grid_get_metadata(couple_config, "grid1")))
    PUT_ERR("ERROR in yac_couple_config_grid_get_metadata\n");

  if (strcmp("grid2_meta",
             yac_couple_config_grid_get_metadata(couple_config, "grid2")))
    PUT_ERR("ERROR in yac_couple_config_grid_get_metadata\n");

  char const * start_datetime =
    yac_couple_config_get_start_datetime(couple_config);
  if (!start_datetime || strcmp(start_datetime, "2008-03-09T16:05:07"))
    PUT_ERR("ERROR in yac_couple_config_get_start_datetime\n");

  char const * end_datetime =
    yac_couple_config_get_end_datetime(couple_config);
  if (!end_datetime || strcmp(end_datetime, "2008-03-10T16:05:07"))
    PUT_ERR("ERROR in yac_couple_config_get_end_datetime\n");

  yac_couple_config_set_datetime(
    couple_config, "2009-03-09T16:05:07", "2009-03-10T16:05:07");

  if (strcmp(yac_couple_config_get_start_datetime(couple_config),
             "2009-03-09T16:05:07"))
    PUT_ERR("ERROR in yac_couple_config_get_start_datetime\n");

  if (strcmp(yac_couple_config_get_end_datetime(couple_config),
             "2009-03-10T16:05:07"))
    PUT_ERR("ERROR in yac_couple_config_get_end_datetime\n");

  size_t ref_num_couples = 2;
  if (ref_num_couples !=
      yac_couple_config_get_num_couples(couple_config))
    PUT_ERR("ERROR in yac_couple_config_get_num_couples\n");

  char const * ref_component_names[2][2] =
    {{"ICON-A", "ICON-O"}, {"ICON-O", "DUMMY"}};
  size_t ref_num_field_couples[2] = {5, 2};
  char const * ref_field_grid_names[2][5][2] =
    {{{"grid1","grid3"},
      {"grid3","grid1"},
      {"grid1","grid3"},
      {"grid3","grid2"},
      {"grid2","grid3"}},
     {{"grid2","grid4"},
      {"grid4","grid2"}}};
  char const * ref_field_component_names[2][5][2] =
    {{{"ICON-O","ICON-A"},
      {"ICON-A","ICON-O"},
      {"ICON-O","ICON-A"},
      {"ICON-A","ICON-O"},
      {"ICON-O","ICON-A"}},
     {{"ICON-O", "DUMMY"},
      {"DUMMY", "ICON-O"}}};
  double ref_frac_mask_fallback_value[2][5] =
    {{YAC_FRAC_MASK_NO_VALUE,1.0,0.0,0.0,YAC_FRAC_MASK_NO_VALUE},
     {YAC_FRAC_MASK_NO_VALUE, YAC_FRAC_MASK_NO_VALUE}};
  double ref_scale_factor[2][5] =  {{1.0,10.0, 1.0, 0.5,1.0}, {1.0,9.0/5.0}};
  double ref_scale_summand[2][5] = {{0.0, 0.0,-1.0,-0.5,0.0}, {0.0,32.0}};
  size_t ref_field_collection_size[2][5] = {{3,3,4,4,5},{5, 2}};
  char const * ref_field_name[2][5] =
    {{"sea_surface_temperature",
      "wind_speed",
      "water_flux_into_sea_water",
      "grid_eastward_wind",
      "grid_northward_wind"},
     {"grid_northward_wind",
      "manual_field"}};
  int ref_mapping_on_source[2][5] = {{1,1,0,0,1},{1, 0}};
  char const * ref_coupling_period[2][5] =
    {{"10","20","30","40","50"},{"50", "60"}};
  enum yac_reduction_type ref_coupling_period_operation[2][5] =
    {{TIME_NONE, TIME_ACCUMULATE, TIME_AVERAGE, TIME_MINIMUM, TIME_MAXIMUM},
     {TIME_MAXIMUM, TIME_NONE}};
  char const * ref_timestep[2][5][2] =
    {{{"1","10"},{"2","20"},{"3","30"},{"4","40"},{"5","50"}},
     {{"5","50"}, {"6", "60"}}};
  int ref_lag[2][5][2] = {{{0,4},{1,3},{2,2},{3,1},{4,0}},{{4,0},{0,0}}};
  int ref_enforce_write_weight_file[2][5] = {{1,0,1,0,1},{1,0}};
  char const * ref_weight_file_name[2][5] =
    {{"weights1.nc",
      "weights2.nc",
      "weights3.nc",
      "weights4.nc",
      "weights5.nc"},
     {"weights6.nc",
      "weights7.nc"}};
  struct yac_interp_stack_config * ref_interp_stack_config[2][5];
  size_t ref_num_src_mask_names[2][5] = {{1,3,0,0,0}, {0,2}};
  char const * const * ref_src_mask_names[2][5] =
    {{(char const*[]){"src_sst_mask"},
      (char const*[]){"src_wind_mask1", "src_wind_mask2", "src_wind_mask3"},
      NULL, NULL, NULL},
     {NULL, (char const*[]){"src_mask1", "src_mask2"}}};
  char const * ref_tgt_mask_name[2][5] =
    {{"tgt_sst_mask", NULL, NULL, NULL, NULL}, {NULL, "tgt_mask"}};

  ref_interp_stack_config[0][0] = yac_interp_stack_config_new();
  yac_interp_stack_config_add_nnn(
    ref_interp_stack_config[0][0], NNN_AVG, 16, YAC_DEFAULT_GAUSS_SCALE);
  yac_interp_stack_config_add_average(
    ref_interp_stack_config[0][0], AVG_ARITHMETIC, 0);
  yac_interp_stack_config_add_conservative(
    ref_interp_stack_config[0][0], 1, 0, 0, CONSERV_DESTAREA);
  yac_interp_stack_config_add_hcsbb(
    ref_interp_stack_config[0][0]);
  yac_interp_stack_config_add_user_file(
    ref_interp_stack_config[0][0], "weights.nc", "grid1", "grid3");
  yac_interp_stack_config_add_fixed(
    ref_interp_stack_config[0][0], -1.0);

  ref_interp_stack_config[0][1] = yac_interp_stack_config_new();
  yac_interp_stack_config_add_average(
    ref_interp_stack_config[0][1], AVG_ARITHMETIC, 0);
  yac_interp_stack_config_add_nnn(
    ref_interp_stack_config[0][1], NNN_DIST, 2, YAC_DEFAULT_GAUSS_SCALE);
  yac_interp_stack_config_add_conservative(
    ref_interp_stack_config[0][1], 2, 1, 1, CONSERV_FRACAREA);

  ref_interp_stack_config[0][2] = yac_interp_stack_config_new();
  yac_interp_stack_config_add_check(
    ref_interp_stack_config[0][2], "", "");
  yac_interp_stack_config_add_check(
    ref_interp_stack_config[0][2], "check_constructor", "");
  yac_interp_stack_config_add_check(
    ref_interp_stack_config[0][2], "", "check_do_search");
  yac_interp_stack_config_add_check(
    ref_interp_stack_config[0][2], "check_constructor", "check_do_search");
  yac_interp_stack_config_add_nnn(
    ref_interp_stack_config[0][2], NNN_RBF, 4, YAC_DEFAULT_RBF_SCALE);

  ref_interp_stack_config[0][3] = yac_interp_stack_config_new();
  yac_interp_stack_config_add_nnn(
    ref_interp_stack_config[0][3], NNN_GAUSS, 8, 0.2);
  yac_interp_stack_config_add_spmap(
    ref_interp_stack_config[0][3], 5.0 * YAC_RAD,
    YAC_DEFAULT_MAX_SEARCH_DISTANCE, SPMAP_DIST);

  ref_interp_stack_config[0][4] = yac_interp_stack_config_new();
  yac_interp_stack_config_add_creep(
    ref_interp_stack_config[0][4], 5);
  yac_interp_stack_config_add_user_callback(
    ref_interp_stack_config[0][4], "compute_weights");
  yac_interp_stack_config_add_fixed(
    ref_interp_stack_config[0][4], -2.0);

  ref_interp_stack_config[1][0] = yac_interp_stack_config_new();
  yac_interp_stack_config_add_creep(
    ref_interp_stack_config[1][0], -1);
  yac_interp_stack_config_add_user_callback(
    ref_interp_stack_config[1][0], "compute_weights");
  yac_interp_stack_config_add_fixed(
    ref_interp_stack_config[1][0], -1.0);

  ref_interp_stack_config[1][1] = yac_interp_stack_config_new();
  yac_interp_stack_config_add_average(
    ref_interp_stack_config[1][1], AVG_ARITHMETIC, 0);
  yac_interp_stack_config_add_fixed(
    ref_interp_stack_config[1][1], -1.0);

  for (size_t ref_couple_idx = 0; ref_couple_idx < ref_num_couples;
       ++ref_couple_idx) {

    // find matching couple
    size_t couple_idx = SIZE_MAX;
    for (size_t i = 0; (i < ref_num_couples) && (couple_idx == SIZE_MAX);
         ++i) {

      char const * couple_comp_names[2];
      yac_couple_config_get_couple_component_names(
        couple_config, i, couple_comp_names);
      if ((!strcmp(
              couple_comp_names[0],
              ref_component_names[ref_couple_idx][0]) &&
           !strcmp(
              couple_comp_names[1],
              ref_component_names[ref_couple_idx][1])) ||
          (!strcmp(
              couple_comp_names[0],
              ref_component_names[ref_couple_idx][1]) &&
           !strcmp(
              couple_comp_names[1],
              ref_component_names[ref_couple_idx][0])))
        couple_idx = i;
    }

    if (couple_idx == SIZE_MAX) {
      PUT_ERR("ERROR: no matching couple found\n");
      continue;
    }

    size_t num_field_couples =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

    if (ref_num_field_couples[ref_couple_idx] !=
        yac_couple_config_get_num_couple_fields(couple_config, couple_idx))
      PUT_ERR("ERROR in yac_couple_config_get_num_couple_fields\n");

    for (size_t ref_field_couple_idx = 0;
         ref_field_couple_idx < num_field_couples; ++ref_field_couple_idx) {

      // find matching field couple
      size_t field_couple_idx = SIZE_MAX;
      char const * src_component_name, * tgt_component_name;
      char const * src_grid_name, * tgt_grid_name;
      char const * src_field_name, * tgt_field_name;
      char const * ref_src_component_name = ref_field_component_names[ref_couple_idx][ref_field_couple_idx][0];
      char const * ref_tgt_component_name = ref_field_component_names[ref_couple_idx][ref_field_couple_idx][1];
      char const * ref_src_grid_name = ref_field_grid_names[ref_couple_idx][ref_field_couple_idx][0];
      char const * ref_tgt_grid_name = ref_field_grid_names[ref_couple_idx][ref_field_couple_idx][1];
      char const * ref_src_field_name = ref_field_name[ref_couple_idx][ref_field_couple_idx];
      char const * ref_tgt_field_name = ref_field_name[ref_couple_idx][ref_field_couple_idx];
      for (size_t i = 0;
           (i < num_field_couples) && (field_couple_idx == SIZE_MAX); ++i) {

        yac_couple_config_get_field_grid_names(
          couple_config, couple_idx, i, &src_grid_name, &tgt_grid_name);
        yac_couple_config_get_field_couple_component_names(
          couple_config, couple_idx, i, &src_component_name, &tgt_component_name);
        yac_couple_config_get_field_names(
          couple_config, couple_idx, i, &src_field_name, &tgt_field_name);

        if (!strcmp(src_component_name, ref_src_component_name) &&
            !strcmp(tgt_component_name, ref_tgt_component_name) &&
            !strcmp(src_grid_name, ref_src_grid_name) &&
            !strcmp(tgt_grid_name, ref_tgt_grid_name) &&
            !strcmp(src_field_name, ref_src_field_name) &&
            !strcmp(tgt_field_name, ref_tgt_field_name))
          field_couple_idx = i;
      }
      if (field_couple_idx == SIZE_MAX) {
          PUT_ERR("ERROR: no matching field couple found\n");
          continue;
      }

      if (yac_couple_config_get_frac_mask_fallback_value(
            couple_config, src_component_name,
            src_grid_name, src_field_name) !=
          ref_frac_mask_fallback_value[ref_couple_idx][ref_field_couple_idx])
        PUT_ERR("ERROR in yac_couple_config_get_frac_mask_fallback_value\n");

      if (yac_couple_config_get_collection_size(
            couple_config, src_component_name,
            src_grid_name, src_field_name) !=
          ref_field_collection_size[ref_couple_idx][ref_field_couple_idx])
        PUT_ERR("ERROR in yac_couple_config_get_collection_size\n");

      if (yac_couple_config_get_field_role(
            couple_config, src_component_name,
            src_grid_name, src_field_name) != YAC_EXCHANGE_TYPE_SOURCE)
        PUT_ERR("ERROR in yac_couple_config_get_field_role\n");

      if (yac_couple_config_get_field_role(
            couple_config, tgt_component_name,
            tgt_grid_name, tgt_field_name) != YAC_EXCHANGE_TYPE_TARGET)
        PUT_ERR("ERROR in yac_couple_config_get_field_role\n");

      if (yac_couple_config_mapping_on_source(
            couple_config, couple_idx, field_couple_idx) !=
          ref_mapping_on_source[ref_couple_idx][ref_field_couple_idx])
        PUT_ERR("ERROR in yac_couple_config_mapping_on_source\n");

      if (strcmp(
            yac_couple_config_get_coupling_period(
              couple_config, couple_idx, field_couple_idx),
            yac_time_to_ISO(
              ref_coupling_period[ref_couple_idx][ref_field_couple_idx],
              YAC_TIME_UNIT_SECOND)))
        PUT_ERR("ERROR in yac_couple_config_get_coupling_period\n");

      if (yac_couple_config_get_coupling_period_operation(
            couple_config, couple_idx, field_couple_idx) !=
          ref_coupling_period_operation[ref_couple_idx][ref_field_couple_idx])
        PUT_ERR("ERROR in yac_couple_config_get_coupling_period_operation\n");

      char const * src_timestep =
        yac_couple_config_get_source_timestep(
          couple_config, couple_idx, field_couple_idx);
      int src_lag =
        yac_couple_config_get_source_lag(
          couple_config, couple_idx, field_couple_idx);

      if (strcmp(
            src_timestep,
            yac_time_to_ISO(
              ref_timestep[ref_couple_idx][ref_field_couple_idx][0],
              YAC_TIME_UNIT_SECOND)) ||
          (src_lag != ref_lag[ref_couple_idx][ref_field_couple_idx][0]))
        PUT_ERR("ERROR in yac_couple_config_get_source_timestep\n");

      char const * tgt_timestep =
        yac_couple_config_get_target_timestep(
          couple_config, couple_idx, field_couple_idx);
      int tgt_lag =
        yac_couple_config_get_target_lag(
          couple_config, couple_idx, field_couple_idx);

      if (strcmp(
            tgt_timestep,
            yac_time_to_ISO(
              ref_timestep[ref_couple_idx][ref_field_couple_idx][1],
              YAC_TIME_UNIT_SECOND)) ||
          (tgt_lag != ref_lag[ref_couple_idx][ref_field_couple_idx][1]))
        PUT_ERR("ERROR in yac_couple_config_get_target_timestep\n");

      if (yac_couple_config_get_scale_factor(
            couple_config, couple_idx, field_couple_idx) !=
          ref_scale_factor[ref_couple_idx][ref_field_couple_idx])
        PUT_ERR("ERROR in yac_couple_config_get_scale_factor");

      if (yac_couple_config_get_scale_summand(
            couple_config, couple_idx, field_couple_idx) !=
          ref_scale_summand[ref_couple_idx][ref_field_couple_idx])
        PUT_ERR("ERROR in yac_couple_config_get_scale_summand");

      if (yac_couple_config_enforce_write_weight_file(
            couple_config, couple_idx, field_couple_idx) !=
          ref_enforce_write_weight_file[ref_couple_idx][ref_field_couple_idx])
        PUT_ERR("ERROR in yac_couple_config_enforce_write_weight_file\n");

      if (ref_enforce_write_weight_file[ref_couple_idx][ref_field_couple_idx]) {
        char const * weight_file_name =
          yac_couple_config_get_weight_file_name(
            couple_config, couple_idx, field_couple_idx);
        if (strcmp(
              weight_file_name, ref_weight_file_name[ref_couple_idx][ref_field_couple_idx]))
          PUT_ERR("ERROR in yac_couple_config_get_weight_file_name\n");
      }

      if (yac_interp_stack_config_compare(
             ref_interp_stack_config[ref_couple_idx][ref_field_couple_idx],
             yac_couple_config_get_interp_stack(
                couple_config, couple_idx, field_couple_idx)))
        PUT_ERR("ERROR in yac_interp_stack_config");

      yac_interp_stack_config_delete(
        ref_interp_stack_config[ref_couple_idx][ref_field_couple_idx]);

      char const * const * src_mask_names;
      size_t num_src_mask_names;
      yac_couple_config_get_src_mask_names(
        couple_config, couple_idx, field_couple_idx,
        &src_mask_names, &num_src_mask_names);
      if (ref_num_src_mask_names[ref_couple_idx][ref_field_couple_idx] !=
          num_src_mask_names)
        PUT_ERR("ERROR in yac_couple_config_get_src_mask_names");
      for (size_t i = 0; i < num_src_mask_names; ++i)
        if (strcmp(ref_src_mask_names[ref_couple_idx][ref_field_couple_idx][i],
                   src_mask_names[i]))
          PUT_ERR("ERROR in yac_couple_config_get_src_mask_names");

      char const * tgt_mask_name =
        yac_couple_config_get_tgt_mask_name(
          couple_config, couple_idx, field_couple_idx);
      if ((ref_tgt_mask_name[ref_couple_idx][ref_field_couple_idx] == NULL) !=
          (tgt_mask_name == NULL))
        PUT_ERR("ERROR in yac_couple_config_get_tgt_mask_name");
      if (tgt_mask_name != NULL)
        if (strcmp(ref_tgt_mask_name[ref_couple_idx][ref_field_couple_idx],
                   tgt_mask_name))
          PUT_ERR("ERROR in yac_couple_config_get_tgt_mask_name");
    }
  }

  yac_couple_config_delete(couple_config);
}

static struct yac_couple_config * generate_couple_config_from_YAML_parallel(
  char const * config_filename) {

  int size, rank;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  struct yac_couple_config * couple_config = yac_couple_config_new();

  yac_couple_config_add_grid(couple_config, "grid1");
  yac_couple_config_add_grid(couple_config, "grid2");
  yac_couple_config_add_grid(couple_config, "grid3");
  yac_couple_config_add_grid(couple_config, "grid4");

  if (rank == 0)
    yac_couple_config_grid_set_metadata(
      couple_config, "grid1", "grid1_meta");
  if (rank == 1)
    yac_couple_config_grid_set_metadata(
      couple_config, "grid2", "grid2_meta");

  yac_couple_config_add_component(couple_config, "ICON-O");
  yac_couple_config_add_component(couple_config, "ICON-A");
  yac_couple_config_add_component(couple_config, "DUMMY");
  yac_couple_config_add_component(couple_config, "DUMMY_2");

  // only rank 0 provides a valid collection size
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "sea_surface_temperature",
    yac_time_to_ISO("1", YAC_TIME_UNIT_SECOND), (rank == 0)?3:SIZE_MAX);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "sea_surface_temperature",
    yac_time_to_ISO("10", YAC_TIME_UNIT_SECOND), (rank == 0)?3:SIZE_MAX);

  // only rank 1 provides a valid collection size and
  // fractional mask fallback value
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "wind_speed",
    yac_time_to_ISO("2", YAC_TIME_UNIT_SECOND), (rank == 1)?3:SIZE_MAX);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "wind_speed",
    yac_time_to_ISO("20", YAC_TIME_UNIT_SECOND), (rank == 1)?3:SIZE_MAX);
  if (rank == 1) {
    yac_couple_config_field_enable_frac_mask(
      couple_config, "ICON-A", "grid3", "wind_speed", 1.0);
    yac_couple_config_field_enable_frac_mask(
      couple_config, "ICON-O", "grid1", "wind_speed", 1.0);
  }

  // only rank 0 provides a valid timestep
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water",
    (rank == 0)?yac_time_to_ISO("3", YAC_TIME_UNIT_SECOND):NULL, 4);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water",
    (rank == 0)?yac_time_to_ISO("30", YAC_TIME_UNIT_SECOND):NULL, 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water", 0.0);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water", 0.0);

  // only rank 1 provides a valid timestep
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind",
    (rank == 1)?yac_time_to_ISO("4", YAC_TIME_UNIT_SECOND):NULL, 4);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_eastward_wind",
    (rank == 1)?yac_time_to_ISO("40", YAC_TIME_UNIT_SECOND):NULL, 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind", 0.0);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid2", "grid_eastward_wind", 0.0);

  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_northward_wind",
    yac_time_to_ISO("5", YAC_TIME_UNIT_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_northward_wind",
    yac_time_to_ISO("50", YAC_TIME_UNIT_SECOND), 5);

  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_northward_wind",
    yac_time_to_ISO("5", YAC_TIME_UNIT_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "grid_northward_wind",
    yac_time_to_ISO("50", YAC_TIME_UNIT_SECOND), 5);

  // add couple by hand
  if (rank == 0) {
    yac_couple_config_component_add_field(
      couple_config, "ICON-O", "grid2", "manual_field",
      yac_time_to_ISO("60", YAC_TIME_UNIT_SECOND), 2);
    struct yac_interp_stack_config * interp_stack =
      yac_interp_stack_config_new();
    yac_interp_stack_config_add_average(interp_stack, AVG_ARITHMETIC, 0);
    yac_interp_stack_config_add_fixed(interp_stack, -1.0);
    yac_couple_config_def_couple(
      couple_config, "DUMMY", "grid4", "manual_field",
      "ICON-O", "grid2", "manual_field",
      yac_time_to_ISO("60", YAC_TIME_UNIT_SECOND), TIME_NONE,
      interp_stack, 0, 0, NULL, 0, 9.0/5.0, 32.0,
      2, (char const *[]){"src_mask1", "src_mask2"},
      "tgt_mask");
    yac_interp_stack_config_delete(interp_stack);
  } else if (rank == 1) {
    yac_couple_config_component_add_field(
      couple_config, "DUMMY", "grid4", "manual_field",
      yac_time_to_ISO("6", YAC_TIME_UNIT_SECOND), 2);
  }

  // rank zero reads in couplings from YAML configuration file
  if (rank == 0)
    yac_yaml_read_coupling(
      couple_config, config_filename, YAC_YAML_PARSER_DEFAULT);

  // synchronise coupling configuration across all processes
  yac_couple_config_sync(couple_config, MPI_COMM_WORLD);

  yac_couple_config_set_redirect_stdout(couple_config, 1);

  return couple_config;
}

static void write_couple_config_to_YAML(
  struct yac_couple_config * couple_config, char const * config_filename,
  int emit_flags) {

  FILE * yaml_file = fopen(config_filename, "w");

  char * str_couple_config =
    yac_yaml_emit_coupling(couple_config, emit_flags);

  fputs(str_couple_config, yaml_file);
  free(str_couple_config);
  fclose(yaml_file);
}

static struct yac_couple_config * generate_couple_config_from_YAML(
  char const * config_filename, int parse_flags) {

  struct yac_couple_config * couple_config = yac_couple_config_new();
  yac_yaml_read_coupling(
    couple_config, config_filename, parse_flags);

  yac_couple_config_set_redirect_stdout(couple_config, 1);

  // add stuff not included in the configuration file
  yac_couple_config_add_component(couple_config, "DUMMY_2");
  yac_couple_config_grid_set_metadata(
    couple_config, "grid1", "grid1_meta");
  yac_couple_config_grid_set_metadata(
    couple_config, "grid2", "grid2_meta");
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "sea_surface_temperature",
    yac_time_to_ISO("1", YAC_TIME_UNIT_SECOND), 3);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "sea_surface_temperature",
    yac_time_to_ISO("10", YAC_TIME_UNIT_SECOND), 3);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "wind_speed",
    yac_time_to_ISO("2", YAC_TIME_UNIT_SECOND), 3);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "wind_speed", 1.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "wind_speed",
    yac_time_to_ISO("20", YAC_TIME_UNIT_SECOND), 3);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid1", "wind_speed", 1.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water",
    yac_time_to_ISO("3", YAC_TIME_UNIT_SECOND), 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water", 0.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water",
    yac_time_to_ISO("30", YAC_TIME_UNIT_SECOND), 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water", 0.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind",
    yac_time_to_ISO("4", YAC_TIME_UNIT_SECOND), 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind", 0.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_eastward_wind",
    yac_time_to_ISO("40", YAC_TIME_UNIT_SECOND), 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid2", "grid_eastward_wind", 0.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_northward_wind",
    yac_time_to_ISO("5", YAC_TIME_UNIT_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_northward_wind",
    yac_time_to_ISO("50", YAC_TIME_UNIT_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "grid_northward_wind",
    yac_time_to_ISO("50", YAC_TIME_UNIT_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "manual_field",
    yac_time_to_ISO("60", YAC_TIME_UNIT_SECOND), 2);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "manual_field",
    yac_time_to_ISO("6", YAC_TIME_UNIT_SECOND), 2);

  return couple_config;
}
