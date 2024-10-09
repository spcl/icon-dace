// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <mpi.h>

#include "tests.h"
#include "yac.h"
#include "config_yaml.h"
#include "instance.h"
#include "event.h"
#include "geometry.h"
#include "io_utils.h"

char const * ref_start_datetime ="2008-03-09T16:05:07.000";
char const * ref_end_datetime = "2008-03-10T16:05:07.000";

static void check_couple_config(struct yac_couple_config * couple_config);
static void check_couple_config_no_delete(struct yac_couple_config * couple_config);
static struct yac_couple_config * generate_couple_config_from_YAML_parallel(
  char const * config_filename, int parse_flags);
static struct yac_couple_config * generate_couple_config_from_YAML(
  char const * config_filename, int parse_flags);
static void write_couple_config_to_YAML(
  struct yac_couple_config * couple_config, char const * config_filename,
  int emit_flags);
static int compare_string(char const * a, char const * b);

#define DEF_INTERP_STACK(IDX, ...) \
  { \
    struct yac_interp_stack_config * interp_stack = \
      ((ref_field_couple_data[IDX].interp_stack_config = \
          yac_interp_stack_config_new())); \
    __VA_ARGS__ \
  }
#define ADD_INTERP(NAME, ...) \
  yac_interp_stack_config_add_ ## NAME( \
    interp_stack, __VA_ARGS__);
#define ADD_INTERP_NO_PARAM(NAME) \
  yac_interp_stack_config_add_ ## NAME(interp_stack);

int main(int argc, char** argv) {

  MPI_Init(NULL, NULL);

  int size, rank;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  if (size != 2) {
    fputs("ERROR wrong number of processes (has to be 2)", stderr);
    exit(EXIT_FAILURE);
  }

  if (argc != 2) {
    PUT_ERR("ERROR: missing config file directory");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

  char const * filename = "couple_config_test.yaml";
  char * config_filename = malloc(strlen(argv[1]) + 32);

  struct yac_couple_config * couple_config =
    generate_couple_config_from_YAML_parallel(
      strcat(strcpy(config_filename, argv[1]), filename),
      YAC_YAML_PARSER_DEFAULT);

  check_couple_config_no_delete(couple_config);

  // check the emitting of coupling configurations (on rank zero only)
  {

    struct {
      int emit_flags;
      int parse_flags;
      enum yac_text_filetype filetype;
      char const * ext;
    } formats[] =
    {{.emit_flags = YAC_YAML_EMITTER_DEFAULT,
      .parse_flags = YAC_YAML_PARSER_DEFAULT,
      .filetype = YAC_TEXT_FILETYPE_YAML,
      .ext = "yaml"},
     {.emit_flags =  YAC_YAML_EMITTER_JSON,
      .parse_flags = YAC_YAML_PARSER_JSON_FORCE,
      .filetype = YAC_TEXT_FILETYPE_JSON,
      .ext = "json"}};
    enum {NUM_FORMATS = sizeof(formats) / sizeof(formats[0])};

    for (int i = 0; i < NUM_FORMATS; ++i) {

      sprintf(
        config_filename, "temp_couple_config_test.%s", formats[i].ext);

      if (rank == 0) {
        write_couple_config_to_YAML(
          couple_config, config_filename, formats[i].emit_flags);

        check_couple_config(
          generate_couple_config_from_YAML(
            config_filename, formats[i].parse_flags));

        unlink(config_filename);
      }

      MPI_Barrier(MPI_COMM_WORLD);

      char const * debug_global_config_file_def_comp = "test_couple_config_def_comp.yaml";
      char const * debug_global_config_file_sync_def = "test_couple_config_sync_def.json";
      char const * debug_global_config_file_enddef = "test_couple_config_enddef.yaml";

      if (yac_file_exists(config_filename))
        PUT_ERR("ERROR config file should not exist");
      if (yac_file_exists(debug_global_config_file_def_comp))
        PUT_ERR("ERROR config file should not exist");
      if (yac_file_exists(debug_global_config_file_sync_def))
        PUT_ERR("ERROR config file should not exist");
      if (yac_file_exists(debug_global_config_file_enddef))
        PUT_ERR("ERROR config file should not exist");

      char test_reference[32];
      sprintf(test_reference, "test_ref_%s", formats[i].ext);
      int include_definitions = 0;
      yac_couple_config_set_config_output_filename(
        couple_config, config_filename, formats[i].filetype, test_reference,
        include_definitions);
      yac_couple_config_set_config_output_filename(
        couple_config, debug_global_config_file_def_comp,
        YAC_TEXT_FILETYPE_YAML, YAC_INSTANCE_CONFIG_OUTPUT_REF_COMP,
        include_definitions);
      yac_couple_config_sync(couple_config, MPI_COMM_WORLD, test_reference);

      if (rank == 0) {

        if (!yac_file_exists(config_filename))
          PUT_ERR("ERROR config file should not exist");
        if (yac_file_exists(debug_global_config_file_def_comp))
          PUT_ERR("ERROR config file should not exist");
        if (yac_file_exists(debug_global_config_file_sync_def))
          PUT_ERR("ERROR config file should not exist");
        if (yac_file_exists(debug_global_config_file_enddef))
          PUT_ERR("ERROR config file should not exist");

        check_couple_config(
          generate_couple_config_from_YAML(
            config_filename, YAC_YAML_PARSER_JSON_AUTO));

        unlink(config_filename);
      }

      yac_couple_config_sync(
        couple_config, MPI_COMM_WORLD, YAC_INSTANCE_CONFIG_OUTPUT_REF_SYNC);

      if (rank == 0) {

        if (yac_file_exists(config_filename))
            PUT_ERR("ERROR config file should not exist");
        if (yac_file_exists(debug_global_config_file_def_comp))
          PUT_ERR("ERROR config file should not exist");
        if (!yac_file_exists(debug_global_config_file_sync_def))
          PUT_ERR("ERROR config file should not exist");
        if (yac_file_exists(debug_global_config_file_enddef))
          PUT_ERR("ERROR config file should not exist");

        check_couple_config(
          generate_couple_config_from_YAML(
            debug_global_config_file_sync_def, YAC_TEXT_FILETYPE_JSON));

        unlink(debug_global_config_file_sync_def);
      }

      yac_couple_config_sync(
        couple_config, MPI_COMM_WORLD, YAC_INSTANCE_CONFIG_OUTPUT_REF_ENDDEF);

      if (rank == 0) {

        if (yac_file_exists(config_filename))
          PUT_ERR("ERROR config file should not exist");
        if (yac_file_exists(debug_global_config_file_def_comp))
          PUT_ERR("ERROR config file should not exist");
        if (yac_file_exists(debug_global_config_file_sync_def))
          PUT_ERR("ERROR config file should not exist");
        if (!yac_file_exists(debug_global_config_file_enddef))
          PUT_ERR("ERROR config file should not exist");

        check_couple_config(
          generate_couple_config_from_YAML(
            debug_global_config_file_enddef, YAC_TEXT_FILETYPE_YAML));

        unlink(debug_global_config_file_enddef);
      }
    }
  }

  free(config_filename);

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
    yac_couple_config_sync(couple_config, MPI_COMM_WORLD, NULL);
    yac_couple_config_delete(couple_config);
  }

  {
    char const * str_interp_stack_config =
      "- conservative:\n"
      "    enforced_conservation: false\n"
      "    normalisation: fracarea\n"
      "    partial_coverage: false\n"
      "- fixed:\n"
      "    user_value: -999.0\n";
    struct yac_interp_stack_config * interp_stack_config =
      yac_yaml_parse_interp_stack_config_string(
        str_interp_stack_config, YAC_YAML_PARSER_DEFAULT);

    struct yac_interp_stack_config * ref_interp_stack_config =
      yac_interp_stack_config_new();
    yac_interp_stack_config_add_conservative(
      ref_interp_stack_config, 1, 0, 0, YAC_INTERP_CONSERV_FRACAREA);
    yac_interp_stack_config_add_fixed(
      ref_interp_stack_config, -999.0);

    if (yac_interp_stack_config_compare(
          ref_interp_stack_config, interp_stack_config))
      PUT_ERR("ERROR in yac_yaml_parse_interp_stack_config_string");

    yac_interp_stack_config_delete(interp_stack_config);
    yac_interp_stack_config_delete(ref_interp_stack_config);
  }

  { // reproduces a bug in routine merge_components
    struct yac_couple_config * couple_config = yac_couple_config_new();
    struct yac_interp_stack_config * interp_stack_config =
      yac_interp_stack_config_new();
    yac_interp_stack_config_add_fixed(interp_stack_config, -1.0);
    yac_couple_config_add_component(
      couple_config, (rank == 0)?"comp_a":"comp_b");
    yac_couple_config_add_component(
      couple_config, (rank == 0)?"comp_b":"comp_a");
    yac_couple_config_def_couple(
      couple_config, "comp_a", "grid_a", "field_a",
      "comp_b", "grid_b", "field_b", "1", YAC_REDUCTION_TIME_NONE,
      interp_stack_config, 0, 0, NULL, 0, 1.0, 0.0, 0, NULL, NULL);
    yac_interp_stack_config_delete(interp_stack_config);
    yac_couple_config_sync(couple_config, MPI_COMM_WORLD, NULL);
    if (yac_couple_config_get_num_couples(couple_config) != 1)
      PUT_ERR("ERROR in yac_couple_config_get_num_couples\n");
    yac_couple_config_delete(couple_config);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void check_couple_config_no_delete(struct yac_couple_config * couple_config) {

  //---------------------------------------------------------------------------
  // setting up reference data
  //---------------------------------------------------------------------------

  struct {
    char const * name;
    char const * metadata;
  } ref_comp_data[] =
  {{.name = "ICON-O",
    .metadata = "a lot of water"},
   {.name = "ICON-A",
    .metadata = "a lot of hot air"},
   {.name = "DUMMY",
    .metadata = NULL},
   {.name = "DUMMY_2",
    .metadata = NULL},
   {.name = "src_comp",
    .metadata = NULL},
   {.name = "tgt_comp",
    .metadata = NULL}};
  enum {REF_NUM_COMPS = sizeof(ref_comp_data) / sizeof(ref_comp_data[0])};

  struct {
    char const * name;
    char const * output_filename;
    char const * metadata;
  } ref_grid_data[] =
  {{.name = "grid1", .output_filename = "debug_grid.nc", .metadata = "grid1_meta"},
   {.name = "grid2", .output_filename = NULL, .metadata = "grid2_meta"},
   {.name = "grid3", .output_filename = NULL, .metadata = NULL},
   {.name = "grid4", .output_filename = "debug_grid.nc", .metadata = NULL},
   {.name = "src_grid_1", .output_filename = NULL, .metadata = NULL},
   {.name = "src_grid_2", .output_filename = NULL, .metadata = NULL},
   {.name = "tgt_grid_1", .output_filename = NULL, .metadata = NULL},
   {.name = "tgt_grid_2", .output_filename = NULL, .metadata = NULL}};
  enum {REF_NUM_GRIDS = sizeof(ref_grid_data) / sizeof(ref_grid_data[0])};

  struct {
    char const * name;
    char const * metadata;
    size_t ref_comp_data_idx;
    size_t ref_grid_data_idx;
    size_t ref_couple_data_idx;
    double frac_mask_value;
    char const * timestep;
    size_t collection_size;
    int role;
  } ref_field_data[] =
  {{.name = "sea_surface_temperature",
    .metadata = "T in C",
    .ref_comp_data_idx = 0,
    .ref_grid_data_idx = 0,
    .ref_couple_data_idx = 0,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "1",
    .collection_size = 3,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "sea_surface_temperature",
    .metadata = "T in K",
    .ref_comp_data_idx = 1,
    .ref_grid_data_idx = 2,
    .ref_couple_data_idx = 0,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "10",
    .collection_size = 3,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "wind_speed",
    .metadata = "v in km/h",
    .ref_comp_data_idx = 0,
    .ref_grid_data_idx = 0,
    .ref_couple_data_idx = 0,
    .frac_mask_value = 1.0,
    .timestep = "20",
    .collection_size = 3,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "wind_speed",
    .metadata = "v in m/s",
    .ref_comp_data_idx = 1,
    .ref_grid_data_idx = 2,
    .ref_couple_data_idx = 0,
    .frac_mask_value = NAN,
    .timestep = "2",
    .collection_size = 3,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "water_flux_into_sea_water",
    .metadata = NULL,
    .ref_comp_data_idx = 0,
    .ref_grid_data_idx = 0,
    .ref_couple_data_idx = 0,
    .frac_mask_value = 0.0,
    .timestep = "3",
    .collection_size = 4,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "water_flux_into_sea_water",
    .metadata = NULL,
    .ref_comp_data_idx = 1,
    .ref_grid_data_idx = 2,
    .ref_couple_data_idx = 0,
    .frac_mask_value = 0.0,
    .timestep = "30",
    .collection_size = 4,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "grid_eastward_wind",
    .metadata = NULL,
    .ref_comp_data_idx = 0,
    .ref_grid_data_idx = 1,
    .ref_couple_data_idx = 0,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "40",
    .collection_size = 4,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "grid_eastward_wind",
    .metadata = NULL,
    .ref_comp_data_idx = 1,
    .ref_grid_data_idx = 2,
    .ref_couple_data_idx = 0,
    .frac_mask_value = 0.0,
    .timestep = "4",
    .collection_size = 4,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "grid_northward_wind",
    .metadata = NULL,
    .ref_comp_data_idx = 0,
    .ref_grid_data_idx = 1,
    .ref_couple_data_idx = 0,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "5",
    .collection_size = 5,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "grid_northward_wind",
    .metadata = NULL,
    .ref_comp_data_idx = 1,
    .ref_grid_data_idx = 2,
    .ref_couple_data_idx = 0,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "50",
    .collection_size = 5,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "grid_northward_wind",
    .metadata = NULL,
    .ref_comp_data_idx = 2,
    .ref_grid_data_idx = 3,
    .ref_couple_data_idx = 1,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "50",
    .collection_size = 5,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "manual_field",
    .metadata = NULL,
    .ref_comp_data_idx = 0,
    .ref_grid_data_idx = 1,
    .ref_couple_data_idx = 1,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "60",
    .collection_size = 2,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "manual_field",
    .metadata = NULL,
    .ref_comp_data_idx = 2,
    .ref_grid_data_idx = 3,
    .ref_couple_data_idx = 1,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "6",
    .collection_size = 2,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "manual_field_uncoupled_a",
    .metadata = NULL,
    .ref_comp_data_idx = 2,
    .ref_grid_data_idx = 3,
    .ref_couple_data_idx = 1,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "6",
    .collection_size = 2,
    .role = YAC_EXCHANGE_TYPE_NONE},
   {.name = "manual_field_uncoupled_b",
    .metadata = NULL,
    .ref_comp_data_idx = 2,
    .ref_grid_data_idx = 3,
    .ref_couple_data_idx = 1,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "6",
    .collection_size = SIZE_MAX,
    .role = YAC_EXCHANGE_TYPE_NONE},
   {.name = "manual_field_uncoupled_c",
    .metadata = NULL,
    .ref_comp_data_idx = 2,
    .ref_grid_data_idx = 3,
    .ref_couple_data_idx = 1,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = NULL,
    .collection_size = 2,
    .role = YAC_EXCHANGE_TYPE_NONE},
   {.name = "multi_grid_field",
    .metadata = NULL,
    .ref_comp_data_idx = 4,
    .ref_grid_data_idx = 4,
    .ref_couple_data_idx = 2,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "1",
    .collection_size = 1,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "multi_grid_field",
    .metadata = NULL,
    .ref_comp_data_idx = 4,
    .ref_grid_data_idx = 5,
    .ref_couple_data_idx = 2,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "1",
    .collection_size = 1,
    .role = YAC_EXCHANGE_TYPE_SOURCE},
   {.name = "multi_grid_field",
    .metadata = NULL,
    .ref_comp_data_idx = 5,
    .ref_grid_data_idx = 6,
    .ref_couple_data_idx = 2,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "1",
    .collection_size = 1,
    .role = YAC_EXCHANGE_TYPE_TARGET},
   {.name = "multi_grid_field",
    .metadata = NULL,
    .ref_comp_data_idx = 5,
    .ref_grid_data_idx = 7,
    .ref_couple_data_idx = 2,
    .frac_mask_value = YAC_FRAC_MASK_NO_VALUE,
    .timestep = "1",
    .collection_size = 1,
    .role = YAC_EXCHANGE_TYPE_TARGET}};
  enum {REF_NUM_FIELDS = sizeof(ref_field_data) / sizeof(ref_field_data[0])};

  struct {
    size_t ref_comp_data_idxs[2];
  } ref_couple_data[] =
  {{.ref_comp_data_idxs = {0, 1}},
   {.ref_comp_data_idxs = {0, 2}},
   {.ref_comp_data_idxs = {4, 5}}};
  enum {REF_NUM_COUPLES = sizeof(ref_couple_data) / sizeof(ref_couple_data[0])};

  struct {
    size_t ref_src_field_data_idx;
    size_t ref_tgt_field_data_idx;
    struct yac_interp_stack_config * interp_stack_config;
    int mapping_on_source;
    int src_lag;
    int tgt_lag;
    char const * src_timestep;
    char const * tgt_timestep;
    char const * coupling_period;
    enum yac_reduction_type coupling_period_operation;
    int enforce_write_weight_file;
    char const * weight_file_name;
    double scale_factor;
    double scale_summand;
    size_t num_src_mask_names;
    char const * const * src_mask_names;
    char const * tgt_mask_name;
  } ref_field_couple_data[] =
  {{.ref_src_field_data_idx = 0,
    .ref_tgt_field_data_idx = 1,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 0,
    .tgt_lag = 4,
    .src_timestep = "1",
    .tgt_timestep = "10",
    .coupling_period = "10",
    .coupling_period_operation = TIME_NONE,
    .enforce_write_weight_file = 1,
    .weight_file_name = "weights1.nc",
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 1,
    .src_mask_names = (char const*[]){"src_sst_mask"},
    .tgt_mask_name = "tgt_sst_mask"},
   {.ref_src_field_data_idx = 3,
    .ref_tgt_field_data_idx = 2,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 1,
    .tgt_lag = 3,
    .src_timestep = "2",
    .tgt_timestep = "20",
    .coupling_period = "20",
    .coupling_period_operation = TIME_ACCUMULATE,
    .enforce_write_weight_file = 0,
    .weight_file_name = "\0",
    .scale_factor = 10.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 3,
    .src_mask_names = (char const*[]){"src_wind_mask1", "src_wind_mask2", "src_wind_mask3"},
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 4,
    .ref_tgt_field_data_idx = 5,
    .interp_stack_config = NULL,
    .mapping_on_source = 0,
    .src_lag = 2,
    .tgt_lag = 2,
    .src_timestep = "3",
    .tgt_timestep = "30",
    .coupling_period = "30",
    .coupling_period_operation = TIME_AVERAGE,
    .enforce_write_weight_file = 1,
    .weight_file_name = "weights3.nc",
    .scale_factor = 1.0,
    .scale_summand = -1.0,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 7,
    .ref_tgt_field_data_idx = 6,
    .interp_stack_config = NULL,
    .mapping_on_source = 0,
    .src_lag = 3,
    .tgt_lag = 1,
    .src_timestep = "4",
    .tgt_timestep = "40",
    .coupling_period = "40",
    .coupling_period_operation = TIME_MINIMUM,
    .enforce_write_weight_file = 0,
    .weight_file_name = "\0",
    .scale_factor = 0.5,
    .scale_summand = -0.5,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 8,
    .ref_tgt_field_data_idx = 9,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 4,
    .tgt_lag = 0,
    .src_timestep = "5",
    .tgt_timestep = "50",
    .coupling_period = "50",
    .coupling_period_operation = TIME_MAXIMUM,
    .enforce_write_weight_file = 1,
    .weight_file_name = "weights5.nc",
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 8,
    .ref_tgt_field_data_idx = 10,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 4,
    .tgt_lag = 0,
    .src_timestep = "5",
    .tgt_timestep = "50",
    .coupling_period = "50",
    .coupling_period_operation = TIME_MAXIMUM,
    .enforce_write_weight_file = 1,
    .weight_file_name = "weights6.nc",
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 12,
    .ref_tgt_field_data_idx = 11,
    .interp_stack_config = NULL,
    .mapping_on_source = 0,
    .src_lag = 0,
    .tgt_lag = 0,
    .src_timestep = "6",
    .tgt_timestep = "60",
    .coupling_period = "60",
    .coupling_period_operation = TIME_NONE,
    .enforce_write_weight_file = 0,
    .weight_file_name = "\0",
    .scale_factor = 9.0/5.0,
    .scale_summand = 32.0,
    .num_src_mask_names = 2,
    .src_mask_names = (char const*[]){"src_mask1", "src_mask2"},
    .tgt_mask_name = "tgt_mask"},
   {.ref_src_field_data_idx = 16,
    .ref_tgt_field_data_idx = 18,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 0,
    .tgt_lag = 0,
    .src_timestep = "1",
    .tgt_timestep = "1",
    .coupling_period = "1",
    .coupling_period_operation = TIME_NONE,
    .enforce_write_weight_file = 0,
    .weight_file_name = "\0",
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 16,
    .ref_tgt_field_data_idx = 19,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 0,
    .tgt_lag = 0,
    .src_timestep = "1",
    .tgt_timestep = "1",
    .coupling_period = "1",
    .coupling_period_operation = TIME_NONE,
    .enforce_write_weight_file = 0,
    .weight_file_name = "\0",
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 17,
    .ref_tgt_field_data_idx = 18,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 0,
    .tgt_lag = 0,
    .src_timestep = "1",
    .tgt_timestep = "1",
    .coupling_period = "1",
    .coupling_period_operation = TIME_NONE,
    .enforce_write_weight_file = 0,
    .weight_file_name = "\0",
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL},
   {.ref_src_field_data_idx = 17,
    .ref_tgt_field_data_idx = 19,
    .interp_stack_config = NULL,
    .mapping_on_source = 1,
    .src_lag = 0,
    .tgt_lag = 0,
    .src_timestep = "1",
    .tgt_timestep = "1",
    .coupling_period = "1",
    .coupling_period_operation = TIME_NONE,
    .enforce_write_weight_file = 0,
    .weight_file_name = "\0",
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .num_src_mask_names = 0,
    .src_mask_names = NULL,
    .tgt_mask_name = NULL}};
  enum {
    REF_NUM_FIELD_COUPLES =
      sizeof(ref_field_couple_data) / sizeof(ref_field_couple_data[0])};

  DEF_INTERP_STACK(0,
    ADD_INTERP(nnn,
      YAC_INTERP_NNN_AVG, 16,
      YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
      YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT)
    ADD_INTERP(average, YAC_INTERP_AVG_ARITHMETIC, 0)
    ADD_INTERP(conservative, 1, 0, 0, YAC_INTERP_CONSERV_DESTAREA)
    ADD_INTERP_NO_PARAM(hcsbb)
    ADD_INTERP(user_file, "weights.nc")
    ADD_INTERP(fixed, -1.0))

  DEF_INTERP_STACK(1,
    ADD_INTERP(average, YAC_INTERP_AVG_DIST, 1)
    ADD_INTERP(nnn,
      YAC_INTERP_NNN_DIST, 2, YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
      YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT)
    ADD_INTERP(conservative, 2, 1, 1, YAC_INTERP_CONSERV_FRACAREA))

  DEF_INTERP_STACK(2,
    ADD_INTERP(check, "", "")
    ADD_INTERP(check, "check_constructor", "")
    ADD_INTERP(check, "", "check_do_search")
    ADD_INTERP(check, "check_constructor", "check_do_search")
    ADD_INTERP(nnn,
      YAC_INTERP_NNN_RBF, 4, YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT,
      YAC_INTERP_RBF_SCALE_DEFAULT))

  DEF_INTERP_STACK(3,
    ADD_INTERP(nnn, YAC_INTERP_NNN_GAUSS, 8, M_PI_2, 0.2)
    ADD_INTERP(spmap,
      5.0 * YAC_RAD, YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
      YAC_INTERP_SPMAP_DIST, YAC_INTERP_SPMAP_NONE,
      YAC_INTERP_SPMAP_SRC_SPHERE_RADIUS_DEFAULT, 2.0)
    ADD_INTERP(ncc, YAC_INTERP_NCC_DIST, 1)
    ADD_INTERP(nnn,
      YAC_INTERP_NNN_RBF, 4, M_PI_2, YAC_INTERP_RBF_SCALE_DEFAULT))

  DEF_INTERP_STACK(4,
    ADD_INTERP(creep, 5)
    ADD_INTERP(user_callback, "compute_weights")
    ADD_INTERP(fixed, -2.0))

  DEF_INTERP_STACK(5,
    ADD_INTERP(creep, -1)
    ADD_INTERP(user_callback, "compute_weights")
    ADD_INTERP(fixed, -1.0))

  DEF_INTERP_STACK(6,
    ADD_INTERP(average, YAC_INTERP_AVG_ARITHMETIC, 0)
    ADD_INTERP(fixed, -1.0))

  DEF_INTERP_STACK(7,
    ADD_INTERP(fixed, -1.0))
  DEF_INTERP_STACK(8,
    ADD_INTERP(fixed, -1.0))
  DEF_INTERP_STACK(9,
    ADD_INTERP(fixed, -1.0))
  DEF_INTERP_STACK(10,
    ADD_INTERP(fixed, -1.0))

  int ref_missing_definition_is_fatal =
    !MISSING_DEFINITION_IS_FATAL_DEFAULT_VALUE;

  //---------------------------------------------------------------------------
  // testing date stuff
  //---------------------------------------------------------------------------

  char * start_datetime =
    yac_couple_config_get_start_datetime(couple_config);
  if (compare_string(start_datetime, ref_start_datetime))
    PUT_ERR("ERROR in yac_couple_config_get_start_datetime\n");
  free(start_datetime);

  char * end_datetime =
    yac_couple_config_get_end_datetime(couple_config);
  if (compare_string(end_datetime, ref_end_datetime))
    PUT_ERR("ERROR in yac_couple_config_get_end_datetime\n");
  free(end_datetime);

  if (yac_couple_config_get_num_components(couple_config) != REF_NUM_COMPS)
    PUT_ERR("ERROR in yac_couple_config_get_num_components");

  //---------------------------------------------------------------------------
  // testing component stuff
  //---------------------------------------------------------------------------

  for (size_t i = 0; i < REF_NUM_COMPS; ++i) {

    if (!yac_couple_config_component_name_is_valid(
           couple_config, ref_comp_data[i].name))
      PUT_ERR("ERROR in yac_couple_config_component_name_is_valid");

    if (compare_string(
          yac_couple_config_component_get_metadata(
            couple_config, ref_comp_data[i].name),
            ref_comp_data[i].metadata))
      PUT_ERR("ERROR in yac_couple_config_component_get_metadata");

    size_t comp_idx =
      yac_couple_config_get_component_idx(
        couple_config, ref_comp_data[i].name);
    if (compare_string(
          yac_couple_config_get_component_name(couple_config, comp_idx),
          ref_comp_data[i].name))
      PUT_ERR("ERROR in yac_couple_config_get_component_name");

    // count number of associated fields
    size_t ref_num_comp_fields = 0;
    for (size_t j = 0; j < REF_NUM_FIELDS; ++j)
      if (ref_field_data[j].ref_comp_data_idx == i)
        ++ref_num_comp_fields;

    if (yac_couple_config_get_num_fields(couple_config, comp_idx) !=
        ref_num_comp_fields)
      PUT_ERR("ERROR in yac_couple_config_get_num_fields");
  }

  if (yac_couple_config_component_name_is_valid(couple_config, "INVALID"))
    PUT_ERR("ERROR in yac_couple_config_component_name_is_valid");

  //---------------------------------------------------------------------------
  // testing grid stuff
  //---------------------------------------------------------------------------

  if (yac_couple_config_get_num_grids(couple_config) != REF_NUM_GRIDS)
    PUT_ERR("ERROR in yac_couple_config_get_num_grids");

  for (size_t i = 0; i < REF_NUM_GRIDS; ++i) {

    if (!yac_couple_config_contains_grid_name(
           couple_config, ref_grid_data[i].name))
      PUT_ERR("ERROR in yac_couple_config_contains_grid_name\n");

    if (compare_string(
          yac_couple_config_grid_get_output_filename(
            couple_config, ref_grid_data[i].name),
          ref_grid_data[i].output_filename))
      PUT_ERR("ERROR in yac_couple_config_grid_get_output_filename");

    if (compare_string(
          yac_couple_config_grid_get_metadata(
            couple_config, ref_grid_data[i].name), ref_grid_data[i].metadata))
      PUT_ERR("ERROR in yac_couple_config_grid_get_metadata");

    size_t grid_idx =
      yac_couple_config_get_grid_idx(
        couple_config, ref_grid_data[i].name);
    if (compare_string(
          yac_couple_config_get_grid_name(couple_config, grid_idx),
          ref_grid_data[i].name))
      PUT_ERR("ERROR in yac_couple_config_get_grid_name");
  }

  if (yac_couple_config_contains_grid_name(couple_config, "grid5"))
    PUT_ERR("ERROR in yac_couple_config_contains_grid_name\n");

  //---------------------------------------------------------------------------
  // testing field stuff
  //---------------------------------------------------------------------------

  for (size_t i = 0; i < REF_NUM_FIELDS; ++i) {

    char const * field_comp_name =
      ref_comp_data[ref_field_data[i].ref_comp_data_idx].name;
    char const * field_grid_name =
      ref_grid_data[ref_field_data[i].ref_grid_data_idx].name;
    char const * field_name = ref_field_data[i].name;

    if (compare_string(
          yac_couple_config_field_get_metadata(
            couple_config, field_comp_name, field_grid_name, field_name),
          ref_field_data[i].metadata))
      PUT_ERR("ERROR in yac_couple_config_field_get_metadata");

    // use memcmp to compare the values because they can be nan
    double frac_mask_fallback_value =
      yac_couple_config_get_frac_mask_fallback_value(
        couple_config, field_comp_name, field_grid_name, field_name);
    if (memcmp(
          &frac_mask_fallback_value, &ref_field_data[i].frac_mask_value,
          sizeof(frac_mask_fallback_value)))
      PUT_ERR("ERROR in yac_couple_config_get_frac_mask_fallback_value");

    size_t comp_idx =
      yac_couple_config_get_component_idx(
        couple_config, ref_comp_data[ref_field_data[i].ref_comp_data_idx].name);
    size_t grid_idx =
      yac_couple_config_get_grid_idx(
        couple_config, ref_grid_data[ref_field_data[i].ref_grid_data_idx].name);
    size_t field_idx =
      yac_couple_config_get_field_idx(
        couple_config, comp_idx, grid_idx, ref_field_data[i].name);

    if (compare_string(
          yac_couple_config_get_field_grid_name(
            couple_config, comp_idx, field_idx), field_grid_name))
      PUT_ERR("ERROR in yac_couple_config_get_field_grid_name");

    if (compare_string(
          yac_couple_config_get_field_name(
            couple_config, comp_idx, field_idx), field_name))
      PUT_ERR("ERROR in yac_couple_config_get_field_name");

    // timestep can only be queried if it is valid
    if (ref_field_data[i].timestep != NULL)
      if (compare_string(
            yac_couple_config_get_field_timestep(
              couple_config, field_comp_name, field_grid_name, field_name),
            yac_time_to_ISO(ref_field_data[i].timestep, C_SECOND)))
        PUT_ERR("ERROR in yac_couple_config_get_field_timestep");

    // collection size can only be queried if it is valid
    if (ref_field_data[i].collection_size != SIZE_MAX)
      if (yac_couple_config_get_field_collection_size(
            couple_config, field_comp_name, field_grid_name, field_name) !=
          ref_field_data[i].collection_size)
        PUT_ERR("ERROR in yac_couple_config_get_field_collection_size");

    if (yac_couple_config_get_field_role(
          couple_config, field_comp_name, field_grid_name, field_name) !=
        ref_field_data[i].role)
      PUT_ERR("ERROR in yac_couple_config_get_field_role");

    if (yac_couple_config_field_is_valid(couple_config, comp_idx, field_idx) !=
        ((ref_field_data[i].timestep != NULL) &&
         (ref_field_data[i].collection_size != SIZE_MAX)))
      PUT_ERR("ERROR in yac_couple_config_field_is_valid");
  }

  //---------------------------------------------------------------------------
  // testing couple stuff
  //---------------------------------------------------------------------------

  if (yac_couple_config_get_num_couples(couple_config) != REF_NUM_COUPLES)
    PUT_ERR("ERROR in yac_couple_config_get_num_couples\n");

  for (size_t i = 0; i < REF_NUM_COUPLES; ++i) {

    char const * ref_comp_names[2];
    ref_comp_names[0] = ref_comp_data[ref_couple_data[i].ref_comp_data_idxs[0]].name;
    ref_comp_names[1] = ref_comp_data[ref_couple_data[i].ref_comp_data_idxs[1]].name;

    // find matching couple
    size_t couple_idx = SIZE_MAX;
    for (size_t j = 0; (j < REF_NUM_COUPLES) && (couple_idx == SIZE_MAX); ++j) {

      char const * couple_comp_names[2];
      yac_couple_config_get_couple_component_names(
        couple_config, j, couple_comp_names);
      if ((!compare_string(couple_comp_names[0], ref_comp_names[0]) &&
           !compare_string(couple_comp_names[1], ref_comp_names[1])) ||
          (!compare_string(couple_comp_names[0], ref_comp_names[1]) &&
           !compare_string(couple_comp_names[1], ref_comp_names[0])))
        couple_idx = j;
    }

    if (couple_idx == SIZE_MAX) {
      PUT_ERR("ERROR: no matching couple found\n");
      continue;
    }

    size_t num_couple_fields = 0;
    for (size_t j = 0; j < REF_NUM_FIELD_COUPLES; ++j)
      if (((ref_field_data[ref_field_couple_data[j].ref_src_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[i].ref_comp_data_idxs[0]) &&
           (ref_field_data[ref_field_couple_data[j].ref_tgt_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[i].ref_comp_data_idxs[1])) ||
          ((ref_field_data[ref_field_couple_data[j].ref_src_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[i].ref_comp_data_idxs[1]) &&
           (ref_field_data[ref_field_couple_data[j].ref_tgt_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[i].ref_comp_data_idxs[0])))
        ++num_couple_fields;

    if (yac_couple_config_get_num_couple_fields(couple_config, couple_idx) !=
        num_couple_fields)
      PUT_ERR("ERROR in yac_couple_config_get_num_couple_fields");
  }

  //-------------------------------------------------------------------------
  // testing field couple stuff
  //-------------------------------------------------------------------------

  for (size_t i = 0; i < REF_NUM_FIELD_COUPLES; ++i) {

    // find reference couple data
    size_t ref_couple_data_idx;
    for (ref_couple_data_idx = 0; ref_couple_data_idx < REF_NUM_COUPLES; ++ref_couple_data_idx)
      if (((ref_field_data[ref_field_couple_data[i].ref_src_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[ref_couple_data_idx].ref_comp_data_idxs[0]) &&
           (ref_field_data[ref_field_couple_data[i].ref_tgt_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[ref_couple_data_idx].ref_comp_data_idxs[1])) ||
          ((ref_field_data[ref_field_couple_data[i].ref_src_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[ref_couple_data_idx].ref_comp_data_idxs[1]) &&
           (ref_field_data[ref_field_couple_data[i].ref_tgt_field_data_idx].ref_comp_data_idx ==
            ref_couple_data[ref_couple_data_idx].ref_comp_data_idxs[0]))) break;

    char const * ref_comp_names[2];
    ref_comp_names[0] =
      ref_comp_data[ref_couple_data[ref_couple_data_idx].ref_comp_data_idxs[0]].name;
    ref_comp_names[1] =
      ref_comp_data[ref_couple_data[ref_couple_data_idx].ref_comp_data_idxs[1]].name;

    // find matching couple
    size_t couple_idx = SIZE_MAX;
    for (size_t j = 0; (j < REF_NUM_COUPLES) && (couple_idx == SIZE_MAX); ++j) {

      char const * couple_comp_names[2];
      yac_couple_config_get_couple_component_names(
        couple_config, j, couple_comp_names);
      if ((!compare_string(couple_comp_names[0], ref_comp_names[0]) &&
           !compare_string(couple_comp_names[1], ref_comp_names[1])) ||
          (!compare_string(couple_comp_names[0], ref_comp_names[1]) &&
           !compare_string(couple_comp_names[1], ref_comp_names[0])))
        couple_idx = j;
    }

    if (couple_idx == SIZE_MAX) {
      PUT_ERR("ERROR: no matching couple found\n");
      continue;
    }

    size_t num_couples =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

    char const * src_comp_name, * tgt_comp_name;
    char const * src_grid_name, * tgt_grid_name;
    char const * src_field_name, * tgt_field_name;

    // find matching field couple index
    size_t field_couple_idx = SIZE_MAX;
    for (size_t j = 0; (j < num_couples) && (field_couple_idx == SIZE_MAX); ++j) {

      yac_couple_config_get_field_couple_component_names(
        couple_config, couple_idx, j, &src_comp_name, &tgt_comp_name);
      yac_couple_config_get_field_grid_names(
        couple_config, couple_idx, j, &src_grid_name, &tgt_grid_name);
      yac_couple_config_get_field_names(
        couple_config, couple_idx, j, &src_field_name, &tgt_field_name);

      if (!compare_string(
             src_comp_name,
               ref_comp_data[
                 ref_field_data[
                   ref_field_couple_data[i].
                    ref_src_field_data_idx].ref_comp_data_idx].name) &&
          !compare_string(
             tgt_comp_name,
               ref_comp_data[
                 ref_field_data[
                   ref_field_couple_data[i].
                    ref_tgt_field_data_idx].ref_comp_data_idx].name) &&
          !compare_string(
             src_grid_name,
               ref_grid_data[
                 ref_field_data[
                   ref_field_couple_data[i].
                    ref_src_field_data_idx].ref_grid_data_idx].name) &&
          !compare_string(
             tgt_grid_name,
               ref_grid_data[
                 ref_field_data[
                   ref_field_couple_data[i].
                    ref_tgt_field_data_idx].ref_grid_data_idx].name) &&
          !compare_string(
             src_field_name,
               ref_field_data[
                 ref_field_couple_data[i].ref_src_field_data_idx].name) &&
          !compare_string(
             tgt_field_name,
               ref_field_data[
                 ref_field_couple_data[i].ref_tgt_field_data_idx].name))
        field_couple_idx = j;
    }

    if (field_couple_idx == SIZE_MAX) {
      PUT_ERR("ERROR: no matching field couple found\n");
      continue;
    }

    if (yac_interp_stack_config_compare(
            ref_field_couple_data[i].interp_stack_config,
            yac_couple_config_get_interp_stack(
              couple_config, couple_idx, field_couple_idx)))
      PUT_ERR("ERROR in yac_interp_stack_config");

    yac_interp_stack_config_delete(
      ref_field_couple_data[i].interp_stack_config);

    if (yac_couple_config_mapping_on_source(
          couple_config, couple_idx, field_couple_idx) !=
        ref_field_couple_data[i].mapping_on_source)
      PUT_ERR("ERROR in yac_couple_config_mapping_on_source");

    if (yac_couple_config_get_source_lag(
          couple_config, couple_idx, field_couple_idx) !=
        ref_field_couple_data[i].src_lag)
      PUT_ERR("ERROR in yac_couple_config_get_source_lag");

    if (yac_couple_config_get_target_lag(
          couple_config, couple_idx, field_couple_idx) !=
        ref_field_couple_data[i].tgt_lag)
      PUT_ERR("ERROR in yac_couple_config_get_target_lag");

    if (compare_string(
          yac_couple_config_get_source_timestep(
            couple_config, couple_idx, field_couple_idx),
          yac_time_to_ISO(ref_field_couple_data[i].src_timestep, C_SECOND)))
      PUT_ERR("ERROR in yac_couple_config_get_source_timestep");

    if (compare_string(
          yac_couple_config_get_target_timestep(
            couple_config, couple_idx, field_couple_idx),
          yac_time_to_ISO(ref_field_couple_data[i].tgt_timestep, C_SECOND)))
      PUT_ERR("ERROR in yac_couple_config_get_target_timestep");

    if (compare_string(
          yac_couple_config_get_coupling_period(
            couple_config, couple_idx, field_couple_idx),
          yac_time_to_ISO(ref_field_couple_data[i].coupling_period, C_SECOND)))
      PUT_ERR("ERROR in yac_couple_config_get_coupling_period");

    if (yac_couple_config_get_coupling_period_operation(
          couple_config, couple_idx, field_couple_idx) !=
        ref_field_couple_data[i].coupling_period_operation)
      PUT_ERR("ERROR in yac_couple_config_get_coupling_period_operation");

    if (yac_couple_config_enforce_write_weight_file(
          couple_config, couple_idx, field_couple_idx) !=
        ref_field_couple_data[i].enforce_write_weight_file)
      PUT_ERR("ERROR in yac_couple_config_enforce_write_weight_file");

    if (compare_string(
          yac_couple_config_get_weight_file_name(
            couple_config, couple_idx, field_couple_idx),
          ref_field_couple_data[i].weight_file_name))
      PUT_ERR("ERROR in yac_couple_config_get_weight_file_name");

    if (yac_couple_config_get_scale_factor(
          couple_config, couple_idx, field_couple_idx) !=
        ref_field_couple_data[i].scale_factor)
      PUT_ERR("ERROR in yac_couple_config_get_scale_factor");

    if (yac_couple_config_get_scale_summand(
          couple_config, couple_idx, field_couple_idx) !=
        ref_field_couple_data[i].scale_summand)
      PUT_ERR("ERROR in yac_couple_config_get_scale_summand");

    char const * const * src_mask_names;
    size_t num_src_mask_names;
    yac_couple_config_get_src_mask_names(
      couple_config, couple_idx, field_couple_idx,
      &src_mask_names, &num_src_mask_names);
    if (ref_field_couple_data[i].num_src_mask_names != num_src_mask_names)
      PUT_ERR("ERROR in yac_couple_config_get_src_mask_names");
    if (ref_field_couple_data[i].num_src_mask_names == num_src_mask_names)
      for (size_t j = 0; j < num_src_mask_names; ++j)
        if (compare_string(
              ref_field_couple_data[i].src_mask_names[j], src_mask_names[j]))
          PUT_ERR("ERROR in yac_couple_config_get_src_mask_names");

    char const * tgt_mask_name =
      yac_couple_config_get_tgt_mask_name(
        couple_config, couple_idx, field_couple_idx);
    if (compare_string(
          ref_field_couple_data[i].tgt_mask_name, tgt_mask_name))
      PUT_ERR("ERROR in yac_couple_config_get_tgt_mask_name");
  }

  //-------------------------------------------------------------------------
  // testing missing_definition_is_fatal-flag
  //-------------------------------------------------------------------------

  if (yac_couple_config_get_missing_definition_is_fatal(couple_config) !=
      ref_missing_definition_is_fatal)
    PUT_ERR("ERROR in yac_couple_config_get_missing_definition_is_fatal");
}

static void check_couple_config(struct yac_couple_config * couple_config) {

  check_couple_config_no_delete(couple_config);

  yac_couple_config_set_datetime(
    couple_config, "2008-03-09T16:05:07", "2008-03-10T16:05:07");

  char * start_datetime =
    yac_couple_config_get_start_datetime(couple_config);
  if (compare_string(start_datetime, ref_start_datetime))
    PUT_ERR("ERROR in yac_couple_config_get_start_datetime\n");
  free(start_datetime);

  char * end_datetime =
    yac_couple_config_get_end_datetime(couple_config);
  if (compare_string(end_datetime, ref_end_datetime))
    PUT_ERR("ERROR in yac_couple_config_get_end_datetime\n");
  free(end_datetime);

  yac_couple_config_delete(couple_config);
}

static struct yac_couple_config * generate_couple_config_from_YAML_parallel(
  char const * config_filename, int parse_flags) {

  int size, rank;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  struct yac_couple_config * couple_config = yac_couple_config_new();

  // same datetime but with additional ".000" at the end
  if (rank == 0)
    yac_couple_config_set_datetime(
      couple_config, "2008-03-09T16:05:07.000", "2008-03-10T16:05:07.000");
  if (rank == 1)
    yac_couple_config_set_datetime(
      couple_config, "2008-03-09T16:05:07", "2008-03-10T16:05:07");

  yac_couple_config_add_grid(couple_config, "grid1");
  yac_couple_config_add_grid(couple_config, "grid2");
  yac_couple_config_add_grid(couple_config, "grid3");
  yac_couple_config_add_grid(couple_config, "grid4");
  yac_couple_config_add_grid(couple_config, "src_grid_1");
  yac_couple_config_add_grid(couple_config, "src_grid_2");
  yac_couple_config_add_grid(couple_config, "tgt_grid_1");
  yac_couple_config_add_grid(couple_config, "tgt_grid_2");

  if (rank == 0)
    yac_couple_config_grid_set_metadata(
      couple_config, "grid1", "grid1_meta");
  if (rank == 1)
    yac_couple_config_grid_set_metadata(
      couple_config, "grid2", "grid2_meta");
  if (rank == 0)
    yac_couple_config_grid_set_output_filename(
      couple_config, "grid4", "debug_grid.nc");

  yac_couple_config_add_component(couple_config, "ICON-O");
  yac_couple_config_add_component(couple_config, "ICON-A");
  yac_couple_config_add_component(couple_config, "DUMMY");
  yac_couple_config_add_component(couple_config, "DUMMY_2");
  yac_couple_config_add_component(couple_config, "src_comp");
  yac_couple_config_add_component(couple_config, "tgt_comp");

  // add metadata to some components
  if (rank == 0)
    yac_couple_config_component_set_metadata(
      couple_config, "ICON-O", "a lot of water");
  if (rank == 1)
    yac_couple_config_component_set_metadata(
      couple_config, "ICON-A", "a lot of hot air");

  // only rank 0 provides a valid collection size and metadata
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "sea_surface_temperature",
    yac_time_to_ISO("1", C_SECOND), (rank == 0)?3:SIZE_MAX);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "sea_surface_temperature",
    yac_time_to_ISO("10", C_SECOND), (rank == 0)?3:SIZE_MAX);
  if (rank == 0) {
    yac_couple_config_field_set_metadata(
      couple_config, "ICON-O", "grid1", "sea_surface_temperature", "T in C");
    yac_couple_config_field_set_metadata(
      couple_config, "ICON-A", "grid3", "sea_surface_temperature", "T in K");
  }

  // only rank 1 provides a valid collection size,
  // fractional mask fallback value, and metadata
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "wind_speed",
    yac_time_to_ISO("2", C_SECOND), (rank == 1)?3:SIZE_MAX);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "wind_speed",
    yac_time_to_ISO("20", C_SECOND), (rank == 1)?3:SIZE_MAX);
  if (rank == 1) {
    yac_couple_config_field_enable_frac_mask(
      couple_config, "ICON-A", "grid3", "wind_speed", NAN);
    yac_couple_config_field_enable_frac_mask(
      couple_config, "ICON-A", "grid3", "wind_speed", NAN);
    yac_couple_config_field_enable_frac_mask(
      couple_config, "ICON-O", "grid1", "wind_speed", 1.0);
    yac_couple_config_field_set_metadata(
      couple_config, "ICON-A", "grid3", "wind_speed", "v in m/s");
    yac_couple_config_field_set_metadata(
      couple_config, "ICON-O", "grid1", "wind_speed", "v in km/h");
  }

  // only rank 0 provides a valid timestep
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water",
    (rank == 0)?yac_time_to_ISO("3", C_SECOND):NULL, 4);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water",
    (rank == 0)?yac_time_to_ISO("30", C_SECOND):NULL, 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water", 0.0);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water", 0.0);

  // only rank 1 provides a valid timestep
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind",
    (rank == 1)?yac_time_to_ISO("4", C_SECOND):NULL, 4);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_eastward_wind",
    (rank == 1)?yac_time_to_ISO("40", C_SECOND):NULL, 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind", 0.0);

  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_northward_wind",
    yac_time_to_ISO("5", C_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_northward_wind",
    yac_time_to_ISO("50", C_SECOND), 5);

  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_northward_wind",
    yac_time_to_ISO("5", C_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "grid_northward_wind",
    yac_time_to_ISO("50", C_SECOND), 5);

  yac_couple_config_component_add_field(
    couple_config, "src_comp", "src_grid_1", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);
  yac_couple_config_component_add_field(
    couple_config, "src_comp", "src_grid_2", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);
  yac_couple_config_component_add_field(
    couple_config, "tgt_comp", "tgt_grid_1", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);
  yac_couple_config_component_add_field(
    couple_config, "tgt_comp", "tgt_grid_2", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);

  // add couple by hand
  if (rank == 0) {
    yac_couple_config_component_add_field(
      couple_config, "ICON-O", "grid2", "manual_field",
      yac_time_to_ISO("60", C_SECOND), 2);
    struct yac_interp_stack_config * interp_stack =
      yac_interp_stack_config_new();
    yac_interp_stack_config_add_average(interp_stack, YAC_INTERP_AVG_ARITHMETIC, 0);
    yac_interp_stack_config_add_fixed(interp_stack, -1.0);
    yac_couple_config_def_couple(
      couple_config, "DUMMY", "grid4", "manual_field",
      "ICON-O", "grid2", "manual_field",
      yac_time_to_ISO("60", C_SECOND), TIME_NONE,
      interp_stack, 0, 0, NULL, 0, 9.0/5.0, 32.0,
      2, (char const *[]){"src_mask1", "src_mask2"},
      "tgt_mask");
    yac_interp_stack_config_delete(interp_stack);
  } else if (rank == 1) {
    yac_couple_config_component_add_field(
      couple_config, "DUMMY", "grid4", "manual_field",
      yac_time_to_ISO("6", C_SECOND), 2);
  }

  // rank 1 adds some uncoupled fields to component DUMMY
  if (rank == 1) {
    yac_couple_config_component_add_field(
      couple_config, "DUMMY", "grid4", "manual_field_uncoupled_a",
      yac_time_to_ISO("6", C_SECOND), 2);
    yac_couple_config_component_add_field(
      couple_config, "DUMMY", "grid4", "manual_field_uncoupled_b",
      yac_time_to_ISO("6", C_SECOND), SIZE_MAX);
    yac_couple_config_component_add_field(
      couple_config, "DUMMY", "grid4", "manual_field_uncoupled_c",
      NULL, 2);
  }

  // rank zero reads in couplings from YAML configuration file
  if (rank == 0) {
    // reading the file twice should not cause any issues
    yac_yaml_read_coupling(
      couple_config, config_filename, parse_flags);
    yac_yaml_read_coupling(
      couple_config, config_filename, parse_flags);
  }

  if (rank == 1)
    yac_couple_config_set_missing_definition_is_fatal(
      couple_config, !MISSING_DEFINITION_IS_FATAL_DEFAULT_VALUE);

  // synchronise coupling configuration across all processes
  yac_couple_config_sync(couple_config, MPI_COMM_WORLD, NULL);

  return couple_config;
}

static void write_couple_config_to_YAML(
  struct yac_couple_config * couple_config, char const * config_filename,
  int emit_flags) {

  FILE * yaml_file = fopen(config_filename, "w");

  int include_definitions = 0;
  char * str_couple_config =
    yac_yaml_emit_coupling(couple_config, emit_flags, include_definitions);

  fputs(str_couple_config, yaml_file);
  free(str_couple_config);
  fclose(yaml_file);
}

static struct yac_couple_config * generate_couple_config_from_YAML(
  char const * config_filename, int parse_flags) {

  struct yac_couple_config * couple_config = yac_couple_config_new();
  yac_yaml_read_coupling(
    couple_config, config_filename, parse_flags);

  // add stuff not included in the configuration file
  yac_couple_config_add_component(couple_config, "DUMMY_2");
  yac_couple_config_component_set_metadata(
    couple_config, "ICON-O", "a lot of water");
  yac_couple_config_component_set_metadata(
    couple_config, "ICON-A", "a lot of hot air");
  yac_couple_config_grid_set_metadata(
    couple_config, "grid1", "grid1_meta");
  yac_couple_config_grid_set_metadata(
    couple_config, "grid2", "grid2_meta");
  yac_couple_config_grid_set_output_filename(
    couple_config, "grid4", "debug_grid.nc");
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "sea_surface_temperature",
    yac_time_to_ISO("1", C_SECOND), 3);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "sea_surface_temperature",
    yac_time_to_ISO("10", C_SECOND), 3);
  yac_couple_config_field_set_metadata(
    couple_config, "ICON-O", "grid1", "sea_surface_temperature", "T in C");
  yac_couple_config_field_set_metadata(
    couple_config, "ICON-A", "grid3", "sea_surface_temperature", "T in K");
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "wind_speed",
    yac_time_to_ISO("2", C_SECOND), 3);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "wind_speed", NAN);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "wind_speed",
    yac_time_to_ISO("20", C_SECOND), 3);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid1", "wind_speed", 1.0);
  yac_couple_config_field_set_metadata(
    couple_config, "ICON-A", "grid3", "wind_speed", "v in m/s");
  yac_couple_config_field_set_metadata(
    couple_config, "ICON-O", "grid1", "wind_speed", "v in km/h");
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water",
    yac_time_to_ISO("3", C_SECOND), 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-O", "grid1", "water_flux_into_sea_water", 0.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water",
    yac_time_to_ISO("30", C_SECOND), 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "water_flux_into_sea_water", 0.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind",
    yac_time_to_ISO("4", C_SECOND), 4);
  yac_couple_config_field_enable_frac_mask(
    couple_config, "ICON-A", "grid3", "grid_eastward_wind", 0.0);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_eastward_wind",
    yac_time_to_ISO("40", C_SECOND), 4);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "grid_northward_wind",
    yac_time_to_ISO("5", C_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "ICON-A", "grid3", "grid_northward_wind",
    yac_time_to_ISO("50", C_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "grid_northward_wind",
    yac_time_to_ISO("50", C_SECOND), 5);
  yac_couple_config_component_add_field(
    couple_config, "ICON-O", "grid2", "manual_field",
    yac_time_to_ISO("60", C_SECOND), 2);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "manual_field",
    yac_time_to_ISO("6", C_SECOND), 2);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "manual_field_uncoupled_a",
    yac_time_to_ISO("6", C_SECOND), 2);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "manual_field_uncoupled_b",
    yac_time_to_ISO("6", C_SECOND), SIZE_MAX);
  yac_couple_config_component_add_field(
    couple_config, "DUMMY", "grid4", "manual_field_uncoupled_c",
    NULL, 2);
  yac_couple_config_component_add_field(
    couple_config, "src_comp", "src_grid_1", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);
  yac_couple_config_component_add_field(
    couple_config, "src_comp", "src_grid_2", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);
  yac_couple_config_component_add_field(
    couple_config, "tgt_comp", "tgt_grid_1", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);
  yac_couple_config_component_add_field(
    couple_config, "tgt_comp", "tgt_grid_2", "multi_grid_field",
    yac_time_to_ISO("1", C_SECOND), 1);

  return couple_config;
}

static int compare_string(char const * a, char const * b) {

  int ret = (a == NULL) - (b == NULL);
  if (ret || (a == NULL)) return ret;
  return strcmp((void*)a,(void*)b);
}
