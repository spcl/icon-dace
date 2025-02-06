// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// libfyaml is not very clean...so we have to suppress some warnings

#if defined(__NVCOMPILER)
#  pragma diag_suppress unsigned_compare_with_zero
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpedantic"
#  pragma GCC diagnostic ignored "-Wall"
#  pragma GCC diagnostic ignored "-Wextra"
#endif
#include <libfyaml.h>
#if defined(__NVCOMPILER)
#  pragma diag_default unsigned_compare_with_zero
#elif defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "yac.h"
#include "utils_mci.h"
#include "config_yaml.h"
#include "mtime_calendar.h"
#include "geometry.h"
#include "io_utils.h"
#include "interp_method_avg.h"
#include "interp_method_callback.h"
#include "interp_method_check.h"
#include "interp_method_conserv.h"
#include "interp_method_creep.h"
#include "interp_method_file.h"
#include "interp_method_fixed.h"
#include "interp_method_hcsbb.h"
#include "interp_method_ncc.h"
#include "interp_method_nnn.h"
#include "interp_method_spmap.h"
#include "interp_stack_config.h"
#include "instance.h"
#include "fields.h"

typedef struct fy_document * fy_document_t;
typedef struct fy_node * fy_node_t;
typedef struct fy_node_pair * fy_node_pair_t;

enum {
  EMITTER_DEFAULT = FYECF_DEFAULT,
  EMITTER_JSON = FYECF_MODE_JSON,
  PARSER_DEFAULT = 0,
  PARSER_JSON_AUTO = FYPCF_JSON_AUTO,
  PARSER_JSON_FORCE = FYPCF_JSON_FORCE,
};

int const  YAC_YAML_EMITTER_DEFAULT = EMITTER_DEFAULT;
int const  YAC_YAML_EMITTER_JSON = EMITTER_JSON;
int const  YAC_YAML_PARSER_DEFAULT = PARSER_DEFAULT;
int const  YAC_YAML_PARSER_JSON_AUTO = PARSER_JSON_AUTO;
int const  YAC_YAML_PARSER_JSON_FORCE = PARSER_JSON_FORCE;

char const * yac_time_to_ISO(
  char const * time, enum yac_time_unit_type time_unit);

struct field_couple_buffer {
  struct {
    char const * comp_name;
    struct {
      char const ** name;
      size_t count;
    } grid;
    int lag;
  } src, tgt;
  char const * coupling_period;
  enum yac_reduction_type time_reduction;
  struct yac_interp_stack_config * interp_stack;
  char const * weight_file_name;
  int mapping_on_source;
  double scale_factor, scale_summand;

  struct field_couple_field_names {
    char const * src, * tgt;
  } * field_names;
  size_t num_field_names;

  char const ** src_mask_names;
  size_t num_src_mask_names;
  char const * tgt_mask_name;
};

enum yaml_base_key_types {
  START_DATE,
  END_DATE,
  CALENDAR,
  TIMESTEP_UNIT,
  COUPLING,
  DEBUG,
};

enum yaml_couple_key_types {
  SOURCE_COMPONENT,
  SOURCE_GRID,
  TARGET_COMPONENT,
  TARGET_GRID,
  FIELD,
  COUPLING_PERIOD,
  TIME_REDUCTION,
  SOURCE_LAG,
  TARGET_LAG,
  WEIGHT_FILE_NAME,
  MAPPING_SIDE,
  SCALE_FACTOR,
  SCALE_SUMMAND,
  INTERPOLATION,
  SOURCE_MASK_NAME,
  SOURCE_MASK_NAMES,
  TARGET_MASK_NAME,
};

enum yaml_debug_key_types {
  GLOBAL_CONFIG,
  GLOBAL_DEFS,
  OUTPUT_GRIDS,
  MISSING_DEF
};

enum yaml_debug_sync_loc_key_types {
  SYNC_LOC_DEF_COMP,
  SYNC_LOC_SYNC_DEF,
  SYNC_LOC_ENDDEF,
  SYNC_LOC_COUNT, // has to be the last entry
};

struct debug_config_file_buffer {
  struct debug_config_file {
    char const * name;
    enum yac_text_filetype type;
    int include_definitions;
  } config_file[SYNC_LOC_COUNT];
  char const * sync_loc_ref[SYNC_LOC_COUNT];
};

enum yaml_debug_output_grid_key_types{
  OUTPUT_GRID_GRID_NAME,
  OUTPUT_GRID_FILE_NAME
};

#define CAST_NAME_TYPE_PAIRS(...) (struct yac_name_type_pair[]) {__VA_ARGS__}
#define COUNT_NAME_TYPE_PAIRS(...) \
  sizeof(CAST_NAME_TYPE_PAIRS(__VA_ARGS__)) / sizeof(struct yac_name_type_pair)
#define DEF_NAME_TYPE_PAIR(NAME, TYPE) {.name = #NAME, .type = (int)(TYPE)}
#define DEF_NAME_TYPE_PAIRS(NAME, ...) \
  static const struct yac_name_type_pair NAME [] = {__VA_ARGS__}; \
  static const size_t num_ ## NAME = COUNT_NAME_TYPE_PAIRS(__VA_ARGS__);

DEF_NAME_TYPE_PAIRS(
  yaml_base_keys,
  DEF_NAME_TYPE_PAIR(start_date,    START_DATE),
  DEF_NAME_TYPE_PAIR(end_date,      END_DATE),
  DEF_NAME_TYPE_PAIR(calendar,      CALENDAR),
  DEF_NAME_TYPE_PAIR(timestep_unit, TIMESTEP_UNIT),
  DEF_NAME_TYPE_PAIR(coupling,      COUPLING),
  DEF_NAME_TYPE_PAIR(debug,         DEBUG))

DEF_NAME_TYPE_PAIRS(
  yaml_couple_keys,
  DEF_NAME_TYPE_PAIR(src_component,    SOURCE_COMPONENT),
  DEF_NAME_TYPE_PAIR(src_grid,         SOURCE_GRID),
  DEF_NAME_TYPE_PAIR(tgt_component,    TARGET_COMPONENT),
  DEF_NAME_TYPE_PAIR(tgt_grid,         TARGET_GRID),
  DEF_NAME_TYPE_PAIR(field,            FIELD),
  DEF_NAME_TYPE_PAIR(coupling_period,  COUPLING_PERIOD),
  DEF_NAME_TYPE_PAIR(time_reduction,   TIME_REDUCTION),
  DEF_NAME_TYPE_PAIR(src_lag,          SOURCE_LAG),
  DEF_NAME_TYPE_PAIR(tgt_lag,          TARGET_LAG),
  DEF_NAME_TYPE_PAIR(weight_file_name, WEIGHT_FILE_NAME),
  DEF_NAME_TYPE_PAIR(mapping_side,     MAPPING_SIDE),
  DEF_NAME_TYPE_PAIR(scale_factor,     SCALE_FACTOR),
  DEF_NAME_TYPE_PAIR(scale_summand,    SCALE_SUMMAND),
  DEF_NAME_TYPE_PAIR(interpolation,    INTERPOLATION),
  DEF_NAME_TYPE_PAIR(src_mask_name,    SOURCE_MASK_NAME),
  DEF_NAME_TYPE_PAIR(src_mask_names,   SOURCE_MASK_NAMES),
  DEF_NAME_TYPE_PAIR(tgt_mask_name,    TARGET_MASK_NAME))

DEF_NAME_TYPE_PAIRS(
  yaml_debug_sync_loc_keys,
  DEF_NAME_TYPE_PAIR(def_comp, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(DEF_COMP, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(def_comps, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(DEF_COMPS, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(sync_def, SYNC_LOC_SYNC_DEF),
  DEF_NAME_TYPE_PAIR(SYNC_DEF, SYNC_LOC_SYNC_DEF),
  DEF_NAME_TYPE_PAIR(enddef, SYNC_LOC_ENDDEF),
  DEF_NAME_TYPE_PAIR(ENDDEF, SYNC_LOC_ENDDEF))

DEF_NAME_TYPE_PAIRS(
  yaml_debug_output_grid_keys,
  DEF_NAME_TYPE_PAIR(grid_name, OUTPUT_GRID_GRID_NAME),
  DEF_NAME_TYPE_PAIR(file_name, OUTPUT_GRID_FILE_NAME))

DEF_NAME_TYPE_PAIRS(
  bool_names,
  DEF_NAME_TYPE_PAIR(true,  true),
  DEF_NAME_TYPE_PAIR(TRUE,  true),
  DEF_NAME_TYPE_PAIR(yes,   true),
  DEF_NAME_TYPE_PAIR(YES,   true),
  DEF_NAME_TYPE_PAIR(false, false),
  DEF_NAME_TYPE_PAIR(FALSE, false),
  DEF_NAME_TYPE_PAIR(no,    false),
  DEF_NAME_TYPE_PAIR(NO,    false))

DEF_NAME_TYPE_PAIRS(
  timestep_units,
  DEF_NAME_TYPE_PAIR(millisecond, C_MILLISECOND),
  DEF_NAME_TYPE_PAIR(second,      C_SECOND),
  DEF_NAME_TYPE_PAIR(minute,      C_MINUTE),
  DEF_NAME_TYPE_PAIR(hour,        C_HOUR),
  DEF_NAME_TYPE_PAIR(day,         C_DAY),
  DEF_NAME_TYPE_PAIR(month,       C_MONTH),
  DEF_NAME_TYPE_PAIR(year,        C_YEAR),
  DEF_NAME_TYPE_PAIR(ISO_format,  C_ISO_FORMAT))

DEF_NAME_TYPE_PAIRS(
  time_operations,
  DEF_NAME_TYPE_PAIR(accumulate, TIME_ACCUMULATE),
  DEF_NAME_TYPE_PAIR(average,    TIME_AVERAGE),
  DEF_NAME_TYPE_PAIR(minimum,    TIME_MINIMUM),
  DEF_NAME_TYPE_PAIR(maximum,    TIME_MAXIMUM),
  DEF_NAME_TYPE_PAIR(none,       TIME_NONE))

DEF_NAME_TYPE_PAIRS(
  calendar_types,
  DEF_NAME_TYPE_PAIR(proleptic-gregorian, PROLEPTIC_GREGORIAN),
  DEF_NAME_TYPE_PAIR(360d, YEAR_OF_360_DAYS),
  DEF_NAME_TYPE_PAIR(365d, YEAR_OF_365_DAYS))

DEF_NAME_TYPE_PAIRS(
  mapping_sides,
  DEF_NAME_TYPE_PAIR(source, 1),
  DEF_NAME_TYPE_PAIR(target, 0))

DEF_NAME_TYPE_PAIRS(
  yaml_debug_keys,
  DEF_NAME_TYPE_PAIR(global_config, GLOBAL_CONFIG),
  DEF_NAME_TYPE_PAIR(output_grids, OUTPUT_GRIDS),
  DEF_NAME_TYPE_PAIR(missing_definition_is_fatal, MISSING_DEF))

DEF_NAME_TYPE_PAIRS(
  config_filetypes,
  DEF_NAME_TYPE_PAIR(yaml, YAC_TEXT_FILETYPE_YAML),
  DEF_NAME_TYPE_PAIR(YAML, YAC_TEXT_FILETYPE_YAML),
  DEF_NAME_TYPE_PAIR(json, YAC_TEXT_FILETYPE_JSON),
  DEF_NAME_TYPE_PAIR(JSON, YAC_TEXT_FILETYPE_JSON))

DEF_NAME_TYPE_PAIRS(
  role_types,
  DEF_NAME_TYPE_PAIR(target, TARGET),
  DEF_NAME_TYPE_PAIR(source, SOURCE),
  DEF_NAME_TYPE_PAIR(nothing, NOTHING))

typedef
#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 23 || __NVCOMPILER_MAJOR__ == 24 && __NVCOMPILER_MINOR__ <= 3)
// Older versions of NVHPC have serious problems with unions that contain
// pointers that are not the first members: versions 23.7 and older fail with
// 'Internal compiler error. unhandled type', the newer versions produce code
// that fails at the runtime with
// 'Segmentation fault: address not mapped to object'.
struct
#else
union
#endif
interp_method_parameter_value {
  int enum_value;
  int int_value;
  double dble_value;
  int bool_value;
  char const * str_value;
} interp_method_parameter_value;

enum interp_method_parameter_value_type{
  ENUM_PARAM,
  INT_PARAM,
  DBLE_PARAM,
  BOOL_PARAM,
  STR_PARAM,
  DEG_PARAM,
};

struct interp_method_parameter {

  char const * name;
  enum interp_method_parameter_value_type type;

#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 23 || __NVCOMPILER_MAJOR__ == 24 && __NVCOMPILER_MINOR__ <= 3)
  // Older versions of NVHPC generate invalid LLVM IR for unions that contain
  // structures with doubles: an attempt to initialize a member of type double
  // that is not the first member of the respective structure will lead to a
  // type mismatch.
  struct
#else
  union
#endif
  {
    struct { // enum
      struct yac_name_type_pair const * valid_values;
      size_t num_valid_values;
    } enum_param;
    struct { // integer
      int valid_min, valid_max;
    } int_param;
    struct { // double
      double valid_min, valid_max;
    } dble_param;
    struct { // string
      size_t max_str_len;
    } str_param;
    struct { // bool
      int dummy;
    } bool_param;
  } data;
  interp_method_parameter_value const default_value;
};

#define DEF_INTERP_METHOD_ADD_FUNC(NAME, FUNC) \
  static void add_interp_method_ ## NAME ( \
    struct yac_interp_stack_config * interp_stack, \
    interp_method_parameter_value * parameters, \
    char const * yaml_filename) { \
    char const * routine_name = "add_interp_method_" #NAME; \
    (void)parameters; \
    (void)yaml_filename; \
    (void)routine_name; \
    {FUNC} }
#define DEF_INTERP_METHOD_GET_FUNC(NAME, FUNC) \
  static void get_interp_method_ ## NAME ( \
    union yac_interp_stack_config_entry const * interp_stack_entry, \
    interp_method_parameter_value * parameter_values) { \
    char const * routine_name = "get_interp_method_" #NAME; \
    (void)interp_stack_entry; \
    (void)parameter_values; \
    (void)routine_name; \
    {FUNC} }

#define DEF_ENUM_PARAM(NAME, DEFAULT, ...) \
  {.name = #NAME, \
   .type = ENUM_PARAM, \
   .data.enum_param = \
     {.valid_values = CAST_NAME_TYPE_PAIRS(__VA_ARGS__), \
      .num_valid_values = COUNT_NAME_TYPE_PAIRS(__VA_ARGS__)}, \
   .default_value.enum_value = (int)(DEFAULT)}
#define DEF_INT_PARAM(NAME, DEFAULT, VALID_MIN, VALID_MAX) \
  {.name = #NAME, \
   .type = INT_PARAM, \
   .data.int_param = \
     {.valid_min = (int)(VALID_MIN), \
      .valid_max = (int)(VALID_MAX)}, \
   .default_value.int_value = (int)(DEFAULT)}
#define DEF_DBLE_PARAM(NAME, DEFAULT, VALID_MIN, VALID_MAX) \
  {.name = #NAME, \
   .type = DBLE_PARAM, \
   .data.dble_param = \
     {.valid_min = (double)(VALID_MIN), \
      .valid_max = (double)(VALID_MAX)}, \
   .default_value.dble_value = (double)(DEFAULT)}
#define DEF_DEG_PARAM(NAME, DEFAULT, VALID_MIN, VALID_MAX) \
  {.name = #NAME, \
   .type = DEG_PARAM, \
   .data.dble_param = \
     {.valid_min = (double)(VALID_MIN), \
      .valid_max = (double)(VALID_MAX)}, \
   .default_value.dble_value = (double)(DEFAULT)}
#define DEF_BOOL_PARAM(NAME, DEFAULT) \
  {.name = #NAME, \
   .type = BOOL_PARAM, \
   .default_value.bool_value = (int)(DEFAULT)}
#define DEF_STR_PARAM(NAME, DEFAULT, MAX_STR_LEN) \
  {.name = #NAME, \
   .type = STR_PARAM, \
   .data.str_param.max_str_len = (MAX_STR_LEN), \
   .default_value.str_value = (DEFAULT)}
#define DEF_INTERP_METHOD_PARAM(NAME, ...) \
  static struct interp_method_parameter \
    interp_method_parameters_ ## NAME [] = {__VA_ARGS__};

#define YAML_ASSERT(CHECK, MSG) \
  YAC_ASSERT_F((CHECK), \
    "ERROR(%s): " MSG " in YAML configuration file \"%s\"", \
    routine_name, yaml_filename)
#define YAML_ASSERT_F(CHECK, MSG, ...) \
  YAC_ASSERT_F((CHECK), \
    "ERROR(%s): " MSG " in YAML configuration file \"%s\"", \
    routine_name, __VA_ARGS__, yaml_filename)

#define DEF_INTERP_STACK_ADD(NAME, ...) \
  yac_interp_stack_config_add_ ## NAME ( \
    interp_stack, __VA_ARGS__);

#define DEF_INTERP_STACK_ADD_NO_PARAM(NAME) \
  yac_interp_stack_config_add_ ## NAME ( \
    interp_stack);

#define DEF_INTERP_STACK_GET(NAME, ...) \
  yac_interp_stack_config_entry_get_ ## NAME ( \
    interp_stack_entry, __VA_ARGS__);

#define DEF_INTERP_STACK_GET_NO_PARAM(NAME) \
  {}

#define DEF_INTERP_METHOD(NAME, FUNC_ADD, FUNC_GET, ...) \
  DEF_INTERP_METHOD_ADD_FUNC(NAME, FUNC_ADD) \
  DEF_INTERP_METHOD_GET_FUNC(NAME, FUNC_GET) \
  DEF_INTERP_METHOD_PARAM(NAME, __VA_ARGS__)

// interpolation method average
DEF_INTERP_METHOD(average,
  DEF_INTERP_STACK_ADD(average,
    (enum yac_interp_avg_weight_type)parameters[0].enum_value,
    parameters[1].bool_value),
  enum yac_interp_avg_weight_type reduction_type;
  DEF_INTERP_STACK_GET(average,
    &reduction_type, &parameter_values[1].bool_value)
  parameter_values[0].enum_value = (int)reduction_type;,
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_AVG_WEIGHT_TYPE_DEFAULT,
    DEF_NAME_TYPE_PAIR(distance_weighted, YAC_INTERP_AVG_DIST),
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_AVG_ARITHMETIC),
    DEF_NAME_TYPE_PAIR(barycentric_coordinate, YAC_INTERP_AVG_BARY)),
  DEF_BOOL_PARAM(
    partial_coverage, YAC_INTERP_AVG_PARTIAL_COVERAGE_DEFAULT))

// interpolation method nearest corner cells
DEF_INTERP_METHOD(ncc,
  DEF_INTERP_STACK_ADD(ncc,
    (enum yac_interp_ncc_weight_type)parameters[0].enum_value,
    parameters[1].bool_value),
  enum yac_interp_ncc_weight_type weight_type;
  DEF_INTERP_STACK_GET(ncc,
    &weight_type, &parameter_values[1].bool_value)
  parameter_values[0].enum_value = (int)weight_type;,
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_NCC_WEIGHT_TYPE_DEFAULT,
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_NCC_AVG),
    DEF_NAME_TYPE_PAIR(distance_weighted, YAC_INTERP_NCC_DIST)),
  DEF_BOOL_PARAM(
    partial_coverage, YAC_INTERP_NCC_PARTIAL_COVERAGE_DEFAULT))

// interpolation method n-nearest neighbor
DEF_INTERP_METHOD(nnn,
  DEF_INTERP_STACK_ADD(nnn,
    (enum yac_interp_nnn_weight_type)parameters[0].enum_value,
    (size_t)parameters[1].int_value,
    parameters[2].dble_value,
    parameters[3].dble_value),
  enum yac_interp_nnn_weight_type type;
  size_t n;
  DEF_INTERP_STACK_GET(nnn,
    &type, &n, &parameter_values[2].dble_value,
    &parameter_values[3].dble_value)
  parameter_values[0].enum_value = (int)type;
  parameter_values[1].int_value = (int)n;
  parameter_values[2].dble_value /= YAC_RAD;,
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_NNN_WEIGHTED_DEFAULT,
    DEF_NAME_TYPE_PAIR(distance_weighted,  YAC_INTERP_NNN_DIST),
    DEF_NAME_TYPE_PAIR(gauss_weighted,     YAC_INTERP_NNN_GAUSS),
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_NNN_AVG),
    DEF_NAME_TYPE_PAIR(zero, YAC_INTERP_NNN_ZERO)),
  DEF_INT_PARAM(n, YAC_INTERP_NNN_N_DEFAULT, 1, INT_MAX),
  DEF_DEG_PARAM(
    max_search_distance, YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
    0.0, 179.9999),
  DEF_DBLE_PARAM(
    gauss_scale, YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT, -DBL_MAX, DBL_MAX))

// interpolation method conservative
DEF_INTERP_METHOD(conservative,
  DEF_INTERP_STACK_ADD(conservative,
    parameters[0].int_value,
    parameters[1].bool_value,
    parameters[2].bool_value,
    (enum yac_interp_method_conserv_normalisation)parameters[3].enum_value),
  enum yac_interp_method_conserv_normalisation normalisation;
  DEF_INTERP_STACK_GET(conservative,
    &parameter_values[0].int_value,
    &parameter_values[1].bool_value,
    &parameter_values[2].bool_value, &normalisation)
  parameter_values[3].enum_value = (int)normalisation;,
  DEF_INT_PARAM(order, YAC_INTERP_CONSERV_ORDER_DEFAULT, 1, 2),
  DEF_BOOL_PARAM(
    enforced_conservation, YAC_INTERP_CONSERV_ENFORCED_CONSERV_DEFAULT),
  DEF_BOOL_PARAM(
    partial_coverage, YAC_INTERP_CONSERV_PARTIAL_COVERAGE_DEFAULT),
  DEF_ENUM_PARAM(
    normalisation, YAC_INTERP_CONSERV_NORMALISATION_DEFAULT,
    DEF_NAME_TYPE_PAIR(fracarea, YAC_INTERP_CONSERV_FRACAREA),
    DEF_NAME_TYPE_PAIR(destarea, YAC_INTERP_CONSERV_DESTAREA)))

// interpolation method source to target mapping
DEF_INTERP_METHOD(source_to_target_map,
  DEF_INTERP_STACK_ADD(spmap,
    parameters[0].dble_value,
    parameters[1].dble_value,
    (enum yac_interp_spmap_weight_type)parameters[2].enum_value,
    (enum yac_interp_spmap_scale_type)parameters[3].enum_value,
    parameters[4].dble_value,
    parameters[5].dble_value),
  enum yac_interp_spmap_weight_type weight_type;
  enum yac_interp_spmap_scale_type scale_type;
  DEF_INTERP_STACK_GET(spmap,
    &parameter_values[0].dble_value,
    &parameter_values[1].dble_value,
    &weight_type, &scale_type,
    &parameter_values[4].dble_value,
    &parameter_values[5].dble_value)
  parameter_values[0].dble_value /= YAC_RAD;
  parameter_values[1].dble_value /= YAC_RAD;
  parameter_values[2].enum_value = (int)weight_type;
  parameter_values[3].enum_value = (int)scale_type;,
  DEF_DEG_PARAM(
    spread_distance, YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT, 0.0, 89.9999),
  DEF_DEG_PARAM(
    max_search_distance,
    YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT, 0.0, 179.9999),
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
    DEF_NAME_TYPE_PAIR(distance_weighted, YAC_INTERP_SPMAP_DIST),
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_SPMAP_AVG)),
  DEF_ENUM_PARAM(
    scale, YAC_INTERP_SPMAP_SCALE_DEFAULT,
    DEF_NAME_TYPE_PAIR(none, YAC_INTERP_SPMAP_NONE),
    DEF_NAME_TYPE_PAIR(srcarea, YAC_INTERP_SPMAP_SRCAREA),
    DEF_NAME_TYPE_PAIR(invtgtarea, YAC_INTERP_SPMAP_INVTGTAREA),
    DEF_NAME_TYPE_PAIR(fracarea, YAC_INTERP_SPMAP_FRACAREA)),
  DEF_DBLE_PARAM(
    src_sphere_radius,
    YAC_INTERP_SPMAP_SRC_SPHERE_RADIUS_DEFAULT, DBL_MIN, DBL_MAX),
  DEF_DBLE_PARAM(
    tgt_sphere_radius,
    YAC_INTERP_SPMAP_TGT_SPHERE_RADIUS_DEFAULT, DBL_MIN, DBL_MAX))

// interpolation method fixed
DEF_INTERP_METHOD(fixed,
  YAML_ASSERT(
    parameters[0].dble_value != DBL_MAX,
    "parameter 'user_value' of interpolation method 'fixed' is unset");
  DEF_INTERP_STACK_ADD(fixed, parameters[0].dble_value),
  DEF_INTERP_STACK_GET(fixed,
    &parameter_values[0].dble_value),
  DEF_DBLE_PARAM(
    user_value, YAC_INTERP_FIXED_VALUE_DEFAULT, -DBL_MAX, DBL_MAX))

// interpolation method user file
DEF_INTERP_METHOD(user_file,
  YAML_ASSERT(
    parameters[0].str_value,
    "parameter \"filename\" of interpolation method \"user file\" is unset")
  DEF_INTERP_STACK_ADD(user_file,
    (char*)(parameters[0].str_value)),
  DEF_INTERP_STACK_GET(user_file,
    &parameter_values[0].str_value),
  DEF_STR_PARAM(
    filename, YAC_INTERP_FILE_WEIGHT_FILE_NAME_DEFAULT,
    YAC_MAX_FILE_NAME_LENGTH))

// interpolation method check
DEF_INTERP_METHOD(check,
  DEF_INTERP_STACK_ADD(check,
    (char*)(parameters[0].str_value),
    (char*)(parameters[1].str_value)),
  DEF_INTERP_STACK_GET(check,
    &parameter_values[0].str_value,
    &parameter_values[1].str_value),
  DEF_STR_PARAM(
    constructor_key, YAC_INTERP_CHECK_CONSTRUCTOR_KEY_DEFAULT,
    YAC_MAX_ROUTINE_NAME_LENGTH),
  DEF_STR_PARAM(
    do_search_key, YAC_INTERP_CHECK_DO_SEARCH_KEY_DEFAULT,
    YAC_MAX_ROUTINE_NAME_LENGTH))

// interpolation method Bernstein Bezier
DEF_INTERP_METHOD(bernstein_bezier,
  DEF_INTERP_STACK_ADD_NO_PARAM(hcsbb),
  DEF_INTERP_STACK_GET_NO_PARAM(hcsbb),
  DEF_BOOL_PARAM(dummy, 0))

// interpolation method radial basis function
DEF_INTERP_METHOD(rbf,
  DEF_INTERP_STACK_ADD(nnn,
    (enum yac_interp_nnn_weight_type)YAC_INTERP_NNN_RBF,
    (size_t)parameters[0].int_value,
    parameters[1].dble_value,
    parameters[2].dble_value),
  enum yac_interp_nnn_weight_type type;
  size_t n;
  DEF_INTERP_STACK_GET(nnn,
    &type, &n, &parameter_values[1].dble_value,
    &parameter_values[2].dble_value)
  YAC_ASSERT(
    type == YAC_INTERP_NNN_RBF,
    "ERROR(get_interp_method_nnn): n-nearest-neighbor type missmatch");
  parameter_values[0].int_value = (int)n;
  parameter_values[1].dble_value /= YAC_RAD;
  parameter_values[3].enum_value = (int)0;,
  DEF_INT_PARAM(n, YAC_INTERP_RBF_N_DEFAULT, 1, INT_MAX),
  DEF_DEG_PARAM(
    max_search_distance, YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT,
    0.0, 179.9999),
  DEF_DBLE_PARAM(rbf_scale, YAC_INTERP_RBF_SCALE_DEFAULT, -DBL_MAX, DBL_MAX),
  DEF_ENUM_PARAM(
    rbf_kernel, YAC_INTERP_RBF_KERNEL_DEFAULT,
    DEF_NAME_TYPE_PAIR(gauss_kernel, 0)))

// interpolation method creep
DEF_INTERP_METHOD(creep,
  DEF_INTERP_STACK_ADD(creep, parameters[0].int_value),
  DEF_INTERP_STACK_GET(creep,
    &parameter_values[0].int_value),
  DEF_INT_PARAM(
    creep_distance, YAC_INTERP_CREEP_DISTANCE_DEFAULT, -1, INT_MAX))

// interpolation method user_callback
DEF_INTERP_METHOD(user_callback,
  YAML_ASSERT(
    parameters[0].str_value,
    "parameter \"func_compute_weights\" "
    "of interpolation method \"user callback\" is unset")
  DEF_INTERP_STACK_ADD(user_callback, (char*)(parameters[0].str_value)),
  DEF_INTERP_STACK_GET(user_callback,
    &parameter_values[0].str_value),
  DEF_STR_PARAM(
    func_compute_weights, YAC_INTERP_CALLBACK_COMPUTE_WEIGHTS_KEY_DEFAULT,
    YAC_MAX_ROUTINE_NAME_LENGTH))


#define ADD_INTERPOLATION(NAME, TYPE) \
  {.name = #NAME , \
   .type = TYPE , \
   .add_interpolation = add_interp_method_ ## NAME , \
   .get_interpolation = get_interp_method_ ## NAME , \
   .parameters = interp_method_parameters_ ## NAME , \
   .num_parameters = \
     sizeof(interp_method_parameters_ ## NAME ) / \
     sizeof(interp_method_parameters_ ## NAME [0])}
#define ADD_INTERPOLATION_NO_PARAM(NAME, TYPE) \
  {.name = #NAME , \
   .type = TYPE , \
   .add_interpolation = add_interp_method_ ## NAME , \
   .get_interpolation = get_interp_method_ ## NAME , \
   .parameters = interp_method_parameters_ ## NAME , \
   .num_parameters = 0}

struct yac_interpolation_method {
  char const * name;
  enum yac_interpolation_list type;
  void(*add_interpolation)(
    struct yac_interp_stack_config * interp_stack,
    interp_method_parameter_value * parameters,
    char const * yaml_filename);
  void(*get_interpolation)(
    union yac_interp_stack_config_entry const * interp_stack_entry,
    interp_method_parameter_value * parameter_values);
  struct interp_method_parameter const * parameters;
  size_t num_parameters;
} const interpolation_methods[] =
  {ADD_INTERPOLATION(average, YAC_AVERAGE),
   ADD_INTERPOLATION(ncc, YAC_NEAREST_CORNER_CELLS),
   ADD_INTERPOLATION(nnn, YAC_N_NEAREST_NEIGHBOR),
   ADD_INTERPOLATION(conservative, YAC_CONSERVATIVE),
   ADD_INTERPOLATION(source_to_target_map, YAC_SOURCE_TO_TARGET_MAP),
   ADD_INTERPOLATION(fixed, YAC_FIXED_VALUE),
   ADD_INTERPOLATION(user_file, YAC_USER_FILE),
   ADD_INTERPOLATION(check, YAC_CHECK),
   ADD_INTERPOLATION_NO_PARAM(bernstein_bezier, YAC_BERNSTEIN_BEZIER),
   ADD_INTERPOLATION(rbf, YAC_RADIAL_BASIS_FUNCTION),
   ADD_INTERPOLATION(creep, YAC_CREEP),
   ADD_INTERPOLATION(user_callback, YAC_USER_CALLBACK)};
enum {
  NUM_INTERPOLATION_METHODS =
    sizeof(interpolation_methods)/sizeof(interpolation_methods[0]),
};

static char const * yaml_parse_string_value(
  fy_node_t value_node, char const * name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_string_value";

  YAML_ASSERT_F(
    value_node && fy_node_is_scalar(value_node),
    "unsupported node type for \"%s\" (the node is expected to be scalar)",
    name);

  return fy_node_get_scalar0(value_node);
}

static calendarType yaml_parse_calendar_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_calendar_value";

  char const * calendar_name =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int calendar_type =
    yac_name_type_pair_get_type(
      calendar_types, num_calendar_types, calendar_name);

  YAML_ASSERT_F(
    calendar_type != INT_MAX,
    "\"%s\" is not a valid calendar name", calendar_name);

  return (calendarType)calendar_type;
}

static char const * yaml_parse_timestep_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename,
  enum yac_time_unit_type time_unit) {
  char const * routine_name = "yaml_parse_timestep_value";

  YAML_ASSERT(
    time_unit != TIME_UNIT_UNDEFINED, "time unit is not yet defined");
  YAML_ASSERT_F(
    value_node && fy_node_is_scalar(value_node),
    "unsupported node type for \"%s\" (the node is expected to be scalar)",
    key_name);
  char const * timestep =
    yaml_parse_string_value(value_node, key_name, yaml_filename);
  char const * timestep_iso =
    yac_time_to_ISO(timestep, time_unit);

  YAML_ASSERT_F(
    timestep_iso, "valid to convert timestep \"%s\" to ISO 8601 format",
    timestep);

  return strdup(timestep_iso);
}

static enum yac_time_unit_type yaml_parse_timestep_unit_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_timestep_unit_value";

  char const * timestep_unit_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int timestep_unit =
    yac_name_type_pair_get_type(
      timestep_units, num_timestep_units, timestep_unit_str);

  YAML_ASSERT_F(
    timestep_unit != INT_MAX,
    "\"%s\" is not a valid time step unit", timestep_unit_str);

  return (enum yac_time_unit_type)timestep_unit;
}

static enum yac_reduction_type yaml_parse_time_reduction_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_time_reduction_value";

  char const * time_reduction_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int time_reduction =
    yac_name_type_pair_get_type(
      time_operations, num_time_operations, time_reduction_str);

  YAML_ASSERT_F(
    time_reduction != INT_MAX,
    "\"%s\" is not a valid time reduction type in", time_reduction_str);

  return (enum yac_reduction_type)time_reduction;
}

static int yaml_parse_integer_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_integer_value";

  char const * integer_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  char * endptr;
  long int long_value = strtol(integer_str, &endptr, 10);

  YAML_ASSERT_F(
    (endptr != integer_str) && (*endptr == '\0') &&
    (long_value >= INT_MIN) && (long_value <= INT_MAX),
    "\"%s\" is not a valid integer value", integer_str);

  return (int)long_value;
}

static double yaml_parse_double_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_double_value";

  char const * double_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  char * endptr;
  double dble_value = strtod(double_str, &endptr);

  YAML_ASSERT_F(
    (endptr != double_str) && (*endptr == '\0'),
    "\"%s\" is not a valid double value", double_str);

  return dble_value;
}

static int yaml_parse_enum_value(
  struct yac_name_type_pair const * valid_values, size_t num_valid_values,
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_enum_value";

  char const * value_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int value =
    yac_name_type_pair_get_type(
      valid_values, num_valid_values, value_str);

  YAML_ASSERT_F(
    value != INT_MAX,
    "\"%s\" is not a valid enum value for \"%s\" ", value_str, key_name);

  return value;
}

static void yaml_parse_base_interp_method_node(
  char const ** interpolation_type_str, fy_node_t * parameter_node,
  fy_node_t interp_method_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_base_interp_method_node";

  YAML_ASSERT(
    fy_node_is_scalar(interp_method_node) ||
    fy_node_is_mapping(interp_method_node),
    "unsupported interpolation method node type "
    "(interpolation methods are expected to be defined as either scalar "
    "or maps)");

  switch(fy_node_get_type(interp_method_node)) {
    default:
    case (FYNT_SCALAR):
      *interpolation_type_str =
        yaml_parse_string_value(
          interp_method_node, "interpolation method name", yaml_filename);
      *parameter_node = NULL;
      break;

    case (FYNT_MAPPING):

      YAML_ASSERT(
        fy_node_mapping_item_count(interp_method_node) == 1,
        "base interpolation method node is only allowed to have one pair ");

      fy_node_pair_t base_interp_method_pair =
        fy_node_mapping_get_by_index(interp_method_node, 0);

      fy_node_t base_interp_method_key_node =
        fy_node_pair_key(base_interp_method_pair);

      *interpolation_type_str =
        yaml_parse_string_value(
          base_interp_method_key_node, "interpolation method name",
          yaml_filename);

      fy_node_t base_interp_method_value_node =
        fy_node_pair_value(base_interp_method_pair);

      YAML_ASSERT_F(
        fy_node_is_mapping(base_interp_method_value_node),
        "unsupported base interpolation method value node type "
        "for interpolation method \"%s\" "
        "(interpolation method parameters are expected to be "
        "defined as maps)", *interpolation_type_str);

      *parameter_node = base_interp_method_value_node;
      break;
  }
}

static interp_method_parameter_value
  yaml_parse_interp_method_parameter_value(
  struct interp_method_parameter const * parameter,
  fy_node_t value_node, char const * interpolation_name,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method_parameter_value";

  YAML_ASSERT_F(
    (parameter->type == ENUM_PARAM) ||
    (parameter->type == INT_PARAM) ||
    (parameter->type == DBLE_PARAM) ||
    (parameter->type == DEG_PARAM) ||
    (parameter->type == BOOL_PARAM) ||
    (parameter->type == STR_PARAM),
    "unsupported parameter type for interpolation method \"%s\"",
    interpolation_name);

  interp_method_parameter_value parameter_value;

  switch(parameter->type) {
    default:
    case (ENUM_PARAM):
      parameter_value.enum_value =
        yaml_parse_enum_value(
          parameter->data.enum_param.valid_values,
          parameter->data.enum_param.num_valid_values,
          value_node, "interpolation method enum parameter value",
          yaml_filename);
      break;
    case (INT_PARAM):
      parameter_value.int_value =
        yaml_parse_integer_value(
          value_node, "interpolation method integer parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        (parameter_value.int_value >= parameter->data.int_param.valid_min) &&
        (parameter_value.int_value <= parameter->data.int_param.valid_max),
        "\"%d\" is not a valid integer parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(valid range: %d <= value <= %d)",
        parameter_value.int_value, parameter->name, interpolation_name,
        parameter->data.int_param.valid_min,
        parameter->data.int_param.valid_max);
      break;
    case (DBLE_PARAM):
      parameter_value.dble_value =
        yaml_parse_double_value(
          value_node, "interpolation method double parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        (parameter_value.dble_value >= parameter->data.dble_param.valid_min) &&
        (parameter_value.dble_value <= parameter->data.dble_param.valid_max),
        "\"%lf\" is not a valid double parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(valid range: %e <= value <= %e)",
        parameter_value.dble_value, parameter->name, interpolation_name,
        parameter->data.dble_param.valid_min,
        parameter->data.dble_param.valid_max);
      break;
    case (DEG_PARAM):
      parameter_value.dble_value =
        yaml_parse_double_value(
          value_node, "interpolation method degree parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        (parameter_value.dble_value >= parameter->data.dble_param.valid_min) &&
        (parameter_value.dble_value <= parameter->data.dble_param.valid_max),
        "\"%lf\" is not a valid degree parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(valid range: %e <= value <= %e)",
        parameter_value.dble_value, parameter->name, interpolation_name,
        parameter->data.dble_param.valid_min,
        parameter->data.dble_param.valid_max);
      parameter_value.dble_value *= YAC_RAD;
      break;
    case (BOOL_PARAM):
      parameter_value.bool_value =
        yaml_parse_enum_value(
          bool_names, num_bool_names, value_node,
          "interpolation method bool parameter value", yaml_filename);
      break;
    case (STR_PARAM):
      parameter_value.str_value =
        yaml_parse_string_value(
          value_node, "interpolation method string parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        strlen(parameter_value.str_value) <
          parameter->data.str_param.max_str_len,
        "\"%s\" is not a valid string parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(maximum string length: %d)",
        parameter_value.str_value, parameter->name, interpolation_name,
        (int)(parameter->data.str_param.max_str_len - 1));
      break;
  };

  return parameter_value;
}

static void yaml_parse_interp_method_parameter(
  interp_method_parameter_value * parameter_values,
  struct interp_method_parameter const * parameters, size_t num_parameters,
  fy_node_pair_t parameter_pair, char const * interpolation_name,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method_parameter";

  char const * parameter_name =
    yaml_parse_string_value(
      fy_node_pair_key(parameter_pair),
      "interpolation method parameter name", yaml_filename);

  int found_flag = 0;
  for (size_t i = 0; (i < num_parameters) && !found_flag; ++i)
    if ((found_flag = !strcmp(parameter_name, parameters[i].name)))
      parameter_values[i] =
        yaml_parse_interp_method_parameter_value(
          &parameters[i], fy_node_pair_value(parameter_pair),
          interpolation_name, yaml_filename);

  YAML_ASSERT_F(
    found_flag,
    "\"%s\" is not a valid parameter for interpolation method \"%s\"",
    parameter_name, interpolation_name);
}

static void yaml_parse_interp_method_parameters(
  interp_method_parameter_value * parameter_values,
  struct yac_interpolation_method const * interp_method,
  fy_node_t parameter_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method_parameters";

  // set default parameter values
  for (size_t i = 0; i < interp_method->num_parameters; ++i)
    parameter_values[i] = interp_method->parameters[i].default_value;

  // if there are no parameters
  if (!parameter_node) return;

  YAML_ASSERT_F(
    fy_node_is_mapping(parameter_node),
    "unsupported interpolation method parameter node type "
    "for interpolation method \"%s\" "
    "(interpolation method parameters are expected to be defined as maps)",
    interp_method->name);

  // parse parameter
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(parameter_node, &iter)))
    yaml_parse_interp_method_parameter(
      parameter_values, interp_method->parameters,
      interp_method->num_parameters, pair, interp_method->name,
      yaml_filename);
}

static void yaml_parse_interp_method(
  struct yac_interp_stack_config * interp_stack,
  fy_node_t interp_method_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method";

  char const * interpolation_type_str;
  fy_node_t parameter_node;

  yaml_parse_base_interp_method_node(
    &interpolation_type_str, &parameter_node,
    interp_method_node, yaml_filename);

  struct yac_interpolation_method const * interp_method = NULL;

  for (int i = 0; (i < NUM_INTERPOLATION_METHODS) && (!interp_method); ++i)
    if (!strcmp(interpolation_type_str, interpolation_methods[i].name))
      interp_method = &interpolation_methods[i];

  YAML_ASSERT_F(
    interp_method,
    "\"%s\" is not a valid interpolation method",
    interpolation_type_str);

  interp_method_parameter_value
    interp_method_paramter_values[interp_method->num_parameters];
  yaml_parse_interp_method_parameters(
    interp_method_paramter_values, interp_method,
    parameter_node, yaml_filename);
  interp_method->add_interpolation(
    interp_stack, interp_method_paramter_values, yaml_filename);
}

static void yaml_parse_interp_stack_value(
  struct yac_interp_stack_config * interp_stack,
  fy_node_t interp_stack_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_stack_value";

  YAML_ASSERT(
    fy_node_is_sequence(interp_stack_node),
    "unsupported interpolation stack node type"
    "(interpolation stacks are expected to be defined as a sequence)");

  // parse couplings
  void * iter = NULL;
  fy_node_t interp_stack_item;
  while ((interp_stack_item =
            fy_node_sequence_iterate(interp_stack_node, &iter)))
    yaml_parse_interp_method(
      interp_stack, interp_stack_item, yaml_filename);
}

static struct field_couple_field_names yaml_parse_field_name(
  fy_node_t field_node, const char * yaml_filename) {
  char const * routine_name = "yaml_parse_field_name";

  YAML_ASSERT(
    fy_node_is_scalar(field_node) ||
    fy_node_is_mapping(field_node),
    "unsupported field name node type "
    "(field name is either scalars or a map)");

  struct field_couple_field_names field_name;

  // if the node contains one name for both source and target field
  if (fy_node_is_scalar(field_node)) {

    field_name.src =
      ((field_name.tgt =
          yaml_parse_string_value(field_node, "field name", yaml_filename)));

  // if the node contains different names for the source and target field
  } else {
    field_name.src =
      fy_node_mapping_lookup_scalar0_by_simple_key(
        field_node, "src", (size_t)-1);
    field_name.tgt =
      fy_node_mapping_lookup_scalar0_by_simple_key(
        field_node, "tgt", (size_t)-1);

    YAML_ASSERT(
      field_name.src && field_name.tgt &&
      (fy_node_mapping_item_count(field_node) == 2),
      "invalid field name mapping node "
      "(field name mapping node has to contain two maps "
      "with the keys \"src\" and \"tgt\")")
  }

  return field_name;
}

static void yaml_parse_string_sequence(
  char const *** values, size_t * num_values,
  fy_node_t values_node, char const * sequence_name,
  const char * yaml_filename) {
  char const * routine_name = "yaml_parse_string_sequence";

  YAML_ASSERT_F(
    (*values == NULL) && (*num_values == 0),
    "values have already been set for sequence \"%s\"",
    sequence_name);

  // if the field node contains multiple fields
  if (fy_node_is_sequence(values_node)) {

    *num_values = (size_t)fy_node_sequence_item_count(values_node);
    *values = xmalloc(*num_values * sizeof(**values));
    for (size_t value_idx = 0; value_idx < *num_values; ++value_idx)
      (*values)[value_idx] =
        yaml_parse_string_value(
          fy_node_sequence_get_by_index(values_node, value_idx),
          sequence_name, yaml_filename);
  } else {
    *num_values = 1;
    *values = xmalloc(sizeof(**values));
    **values =
      yaml_parse_string_value(
        values_node, sequence_name, yaml_filename);
  }
}

static void yaml_parse_field_names(
  struct field_couple_field_names ** field_names,
  size_t * num_field_names, fy_node_t fields_node,
  const char * yaml_filename) {

  // if the field node contains multiple fields
  if (fy_node_is_sequence(fields_node)) {

    size_t start_idx = *num_field_names;
    *num_field_names += (size_t)fy_node_sequence_item_count(fields_node);
    *field_names =
      xrealloc(*field_names, *num_field_names * sizeof(**field_names));
    for (size_t i = start_idx; i < *num_field_names; ++i)
      (*field_names)[i] =
        yaml_parse_field_name(
          fy_node_sequence_get_by_index(fields_node, i), yaml_filename);
  } else {
    ++*num_field_names;
    *field_names =
      xrealloc(*field_names, *num_field_names * sizeof(**field_names));
    (*field_names)[*num_field_names-1] =
      yaml_parse_field_name(fields_node, yaml_filename);
  }
}

static void yaml_parse_couple_map_pair(
  struct field_couple_buffer * field_buffer,
  fy_node_pair_t couple_pair, const char * yaml_filename,
  enum yac_time_unit_type time_unit) {

  enum yaml_couple_key_types couple_key_type =
    (enum yaml_couple_key_types)
      yaml_parse_enum_value(
        yaml_couple_keys, num_yaml_couple_keys,
        fy_node_pair_key(couple_pair),
        "couple configuration parameter name", yaml_filename);
  char const * couple_key_name =
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, couple_key_type);

  fy_node_t value_node = fy_node_pair_value(couple_pair);

  switch (couple_key_type) {
    default:
    case (SOURCE_COMPONENT):
      field_buffer->src.comp_name =
        yaml_parse_string_value(value_node, couple_key_name, yaml_filename);
      break;
    case (SOURCE_GRID):
      yaml_parse_string_sequence(
        &field_buffer->src.grid.name, &field_buffer->src.grid.count,
        value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_COMPONENT):
      field_buffer->tgt.comp_name =
        yaml_parse_string_value(value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_GRID):
      yaml_parse_string_sequence(
        &field_buffer->tgt.grid.name, &field_buffer->tgt.grid.count,
        value_node, couple_key_name, yaml_filename);
      break;
    case (FIELD):
      yaml_parse_field_names(
        &(field_buffer->field_names),
        &(field_buffer->num_field_names),
        value_node, yaml_filename);
      break;
    case (COUPLING_PERIOD):
      field_buffer->coupling_period =
        yaml_parse_timestep_value(
          value_node, couple_key_name, yaml_filename, time_unit);
      break;
    case (TIME_REDUCTION):
      field_buffer->time_reduction =
        yaml_parse_time_reduction_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (SOURCE_LAG):
      field_buffer->src.lag =
        yaml_parse_integer_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_LAG):
      field_buffer->tgt.lag =
        yaml_parse_integer_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (WEIGHT_FILE_NAME):
      field_buffer->weight_file_name =
        yaml_parse_string_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (MAPPING_SIDE):
      field_buffer->mapping_on_source =
        yaml_parse_enum_value(
          mapping_sides, num_mapping_sides,
          value_node, couple_key_name, yaml_filename);
      break;
    case (SCALE_FACTOR):
      field_buffer->scale_factor =
        yaml_parse_double_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (SCALE_SUMMAND):
      field_buffer->scale_summand =
        yaml_parse_double_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (INTERPOLATION):
      yaml_parse_interp_stack_value(
        field_buffer->interp_stack, value_node, yaml_filename);
      break;
    case (SOURCE_MASK_NAME):
    case (SOURCE_MASK_NAMES):
      yaml_parse_string_sequence(
        &(field_buffer->src_mask_names),
        &(field_buffer->num_src_mask_names),
        value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_MASK_NAME):
      field_buffer->tgt_mask_name =
        yaml_parse_string_value(
          value_node, couple_key_name, yaml_filename);
      break;
  }
}

static void yaml_parse_couple(
  struct yac_couple_config * couple_config, fy_node_t couple_node,
  char const * yaml_filename, enum yac_time_unit_type time_unit) {
  char const * routine_name = "yaml_parse_couple";

  YAML_ASSERT(
    fy_node_is_mapping(couple_node),
    "unsupported couple node type "
    "(couples are expected to be defined as a mapping)");

  // initialise field configuration buffer with default values
  struct field_couple_buffer field_buffer = {
    .src.comp_name = NULL,
    .src.grid.name = NULL,
    .src.grid.count = 0,
    .tgt.comp_name = NULL,
    .tgt.grid.name = NULL,
    .tgt.grid.count = 0,
    .field_names = NULL,
    .num_field_names = 0,
    .coupling_period = NULL,
    .time_reduction = TIME_NONE,
    .interp_stack = yac_interp_stack_config_new(),
    .src.lag = 0,
    .tgt.lag = 0,
    .weight_file_name = NULL,
    .mapping_on_source = 1,
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .src_mask_names = NULL,
    .num_src_mask_names = 0,
    .tgt_mask_name = NULL};

  // parse couple
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(couple_node, &iter)))
    yaml_parse_couple_map_pair(
      &field_buffer, pair, yaml_filename, time_unit);

  YAML_ASSERT(
    field_buffer.src.comp_name, "missing source component name");
  YAML_ASSERT_F(
    field_buffer.src.grid.count > 0,
    "missing source grid name (component \"%s\")",
    field_buffer.src.comp_name);
  for (size_t i = 0; i < field_buffer.src.grid.count; ++i)
    YAML_ASSERT_F(
      (field_buffer.src.grid.name[i] != NULL) &&
      (field_buffer.src.grid.name[i][0] != '\0'),
      "invalid source grid name (component \"%s\" grid idx %zu)",
      field_buffer.src.comp_name, i);
  YAML_ASSERT(
    field_buffer.tgt.comp_name,
    "missing target component name");
  YAML_ASSERT_F(
    field_buffer.tgt.grid.count > 0,
    "missing target grid name (component \"%s\")",
    field_buffer.tgt.comp_name);
  for (size_t i = 0; i < field_buffer.tgt.grid.count; ++i)
    YAML_ASSERT_F(
      (field_buffer.tgt.grid.name[i] != NULL) &&
      (field_buffer.tgt.grid.name[i][0] != '\0'),
      "invalid target grid name (component \"%s\" grid idx %zu)",
      field_buffer.tgt.comp_name, i);
  YAML_ASSERT_F(
    field_buffer.num_field_names > 0,
    "missing field names "
    "(source component \"%s\" source grid \"%s\" "
    "target component \"%s\" target grid \"%s\")",
    field_buffer.src.comp_name, field_buffer.src.grid.name[0],
    field_buffer.tgt.comp_name, field_buffer.tgt.grid.name[0]);
  YAML_ASSERT_F(
    field_buffer.coupling_period,
    "missing coupling period "
    "(source component \"%s\" source grid \"%s\" "
    "target component \"%s\" target grid \"%s\")",
    field_buffer.src.comp_name, field_buffer.src.grid.name[0],
    field_buffer.tgt.comp_name, field_buffer.tgt.grid.name[0]);

  for (size_t i = 0; i < field_buffer.num_field_names; ++i)
    for (size_t j = 0; j < field_buffer.src.grid.count; ++j)
      for (size_t k = 0; k < field_buffer.tgt.grid.count; ++k)
        yac_couple_config_def_couple(
          couple_config,
          field_buffer.src.comp_name,
          field_buffer.src.grid.name[j],
          field_buffer.field_names[i].src,
          field_buffer.tgt.comp_name,
          field_buffer.tgt.grid.name[k],
          field_buffer.field_names[i].tgt,
          field_buffer.coupling_period,
          field_buffer.time_reduction,
          field_buffer.interp_stack,
          field_buffer.src.lag,
          field_buffer.tgt.lag,
          field_buffer.weight_file_name,
          field_buffer.mapping_on_source,
          field_buffer.scale_factor,
          field_buffer.scale_summand,
          field_buffer.num_src_mask_names,
          field_buffer.src_mask_names,
          field_buffer.tgt_mask_name);

  // cleanup
  free((void*)field_buffer.src.grid.name);
  free((void*)field_buffer.tgt.grid.name);
  free((void*)field_buffer.field_names);
  free((void*)field_buffer.coupling_period);
  yac_interp_stack_config_delete(field_buffer.interp_stack);
  free((void*)field_buffer.src_mask_names);
}

static void yaml_parse_coupling(
  struct yac_couple_config * couple_config,
  fy_node_t coupling_node, char const * yaml_filename,
  enum yac_time_unit_type time_unit) {
  char const * routine_name = "yaml_parse_coupling";

  // check if the coupling node is empty -> nothing to be read
  if (!coupling_node) return;

  YAML_ASSERT(
    fy_node_is_sequence(coupling_node),
    "unsupported coupling node type "
    "(couplings are expected to be defined as a sequence)");

  // parse couplings
  void * iter = NULL;
  fy_node_t couple_node;
  while ((couple_node = fy_node_sequence_iterate(coupling_node, &iter)))
    yaml_parse_couple(couple_config, couple_node, yaml_filename, time_unit);
}

static struct debug_config_file yaml_parse_config_file_value(
  fy_node_t config_file_node, char const * file_type_name,
  const char * yaml_filename) {
  char const * routine_name = "yaml_parse_config_file_value";

  YAML_ASSERT_F(
    fy_node_is_scalar(config_file_node) ||
    fy_node_is_mapping(config_file_node),
    "unsupported config file node type "
    "(%s is either scalar or a map)", file_type_name);

  struct debug_config_file config_file =
    {.name = NULL, .type = YAC_TEXT_FILETYPE_YAML, .include_definitions = 0};

  char * str_buffer = xmalloc(strlen(file_type_name) + 32);

  // if the node contains only the filename
  if (fy_node_is_scalar(config_file_node)) {

    config_file.name =
      yaml_parse_string_value(
        config_file_node,
        strcat(strcpy(str_buffer, file_type_name), " name"),
        yaml_filename);

  // if the node contains the name and the type
  } else {

    fy_node_t filename_node =
      fy_node_mapping_lookup_by_string(
        config_file_node, "filename", (size_t)-1);
    fy_node_t filetype_node =
      fy_node_mapping_lookup_by_string(
        config_file_node, "filetype", (size_t)-1);
    fy_node_t include_definitions_node =
      fy_node_mapping_lookup_by_string(
        config_file_node, "include_definitions", (size_t)-1);

    YAML_ASSERT_F(
      filename_node,
      "invalid %s mapping node "
      "(global config file mapping node has to include a map "
      "with the keys \"filename\")", file_type_name)

    config_file.name =
      yaml_parse_string_value(
        filename_node, strcat(strcpy(str_buffer, file_type_name), " name"),
        yaml_filename);
    config_file.type =
      (filetype_node)?
        (enum yac_text_filetype)
          yaml_parse_enum_value(
            config_filetypes, num_config_filetypes,
            filetype_node, strcat(strcpy(str_buffer, file_type_name), " type"),
            yaml_filename):YAC_TEXT_FILETYPE_YAML;
    config_file.include_definitions =
      (include_definitions_node)?
        yaml_parse_enum_value(
          bool_names, num_bool_names,
          include_definitions_node,
          strcat(strcpy(str_buffer, file_type_name), " include definitions"),
          yaml_filename):0;
  }

  YAML_ASSERT_F(
    config_file.name, "missing filename for %s", file_type_name);

  free(str_buffer);

  return config_file;
}

static void yaml_parse_debug_config_file_map_pair(
  struct debug_config_file_buffer * config_file_buffer,
  fy_node_pair_t config_file_pair, char const * config_file_type_name,
  const char * yaml_filename) {

  enum yaml_debug_sync_loc_key_types
    debug_sync_loc_key_type =
      (enum yaml_debug_sync_loc_key_types)
        yaml_parse_enum_value(
          yaml_debug_sync_loc_keys,
          num_yaml_debug_sync_loc_keys,
          fy_node_pair_key(config_file_pair),
          "config synchronisation location parameter name",
          yaml_filename);
  char const * debug_sync_loc_key_name =
    yac_name_type_pair_get_name(
      yaml_debug_sync_loc_keys, num_yaml_debug_sync_loc_keys,
      debug_sync_loc_key_type);

  fy_node_t value_node = fy_node_pair_value(config_file_pair);

  char config_file_type_name_sync[
    strlen(config_file_type_name) + strlen(debug_sync_loc_key_name) + 8];
  sprintf(
    config_file_type_name_sync, "%s (%s)",
    config_file_type_name, debug_sync_loc_key_name);

  config_file_buffer->config_file[debug_sync_loc_key_type] =
    yaml_parse_config_file_value(
      value_node, config_file_type_name_sync, yaml_filename);
}

static struct debug_config_file_buffer yaml_parse_debug_config_file_buffer(
  fy_node_t config_file_node, char const * config_file_type_name,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug_config_file_buffer";

  YAML_ASSERT_F(
    fy_node_is_mapping(config_file_node),
    "unsupported %s node type "
    "(%s is expected to be defined as a mapping)",
    config_file_type_name, config_file_type_name);

  struct debug_config_file_buffer config_file_buffer;
  config_file_buffer.sync_loc_ref[SYNC_LOC_DEF_COMP] =
    YAC_INSTANCE_CONFIG_OUTPUT_REF_COMP;
  config_file_buffer.sync_loc_ref[SYNC_LOC_SYNC_DEF] =
    YAC_INSTANCE_CONFIG_OUTPUT_REF_SYNC;
  config_file_buffer.sync_loc_ref[SYNC_LOC_ENDDEF] =
    YAC_INSTANCE_CONFIG_OUTPUT_REF_ENDDEF;
  for (int i = 0; i < SYNC_LOC_COUNT; ++i) {
    config_file_buffer.config_file[i].name = NULL;
    config_file_buffer.config_file[i].type = YAC_TEXT_FILETYPE_YAML;
    config_file_buffer.config_file[i].include_definitions = 0;
    YAML_ASSERT_F(
      config_file_buffer.sync_loc_ref[i] != NULL,
      "invalid unsupported synchronisation location (%d) for %s",
      i, config_file_type_name);
  }

  // parse couplings
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(config_file_node, &iter)))
    yaml_parse_debug_config_file_map_pair(
      &config_file_buffer, pair, config_file_type_name, yaml_filename);

  return config_file_buffer;
}

static void yaml_parse_debug_global_config(
  struct yac_couple_config * couple_config, fy_node_t global_config_node,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug_global_config";

  char const * config_file_type_name = "debug global config file";

  YAML_ASSERT_F(
    fy_node_is_mapping(global_config_node),
    "unsupported %s node type "
    "(%s is expected to be defined as a mapping)",
    config_file_type_name, config_file_type_name);

  struct debug_config_file_buffer global_config_buffer =
    yaml_parse_debug_config_file_buffer(
      global_config_node, config_file_type_name, yaml_filename);


  for (int i = 0; i < SYNC_LOC_COUNT; ++i)
    if (global_config_buffer.config_file[i].name != NULL)
      yac_couple_config_set_config_output_filename(
        couple_config, global_config_buffer.config_file[i].name,
        global_config_buffer.config_file[i].type,
        global_config_buffer.sync_loc_ref[i],
        global_config_buffer.config_file[i].include_definitions);
}

static void yaml_parse_output_grid_pair(
  char const ** grid_name, char const ** file_name,
  fy_node_pair_t output_grid_pair, char const * yaml_filename) {

  enum yaml_debug_output_grid_key_types output_grid_key_type =
    (enum yaml_debug_output_grid_key_types)
      yaml_parse_enum_value(
        yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
        fy_node_pair_key(output_grid_pair),
        "output grid parameter name", yaml_filename);
  char const * debug_output_key_name =
    yac_name_type_pair_get_name(
      yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
      output_grid_key_type);

  fy_node_t value_node = fy_node_pair_value(output_grid_pair);

  switch(output_grid_key_type) {
    default:
    case(OUTPUT_GRID_GRID_NAME): {
      *grid_name =
        yaml_parse_string_value(value_node, debug_output_key_name, yaml_filename);
      break;
    }
    case(OUTPUT_GRID_FILE_NAME): {
      *file_name =
        yaml_parse_string_value(value_node, debug_output_key_name, yaml_filename);
      break;
    }
  };
}

static void yaml_parse_output_grid(
  struct yac_couple_config * couple_config, fy_node_t output_grid_node,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_output_grid";

  YAML_ASSERT(
    fy_node_is_mapping(output_grid_node),
    "unsupported output grid node type "
    "(output grids are expected to be defined as a mapping)");

  char const * grid_name = NULL;
  char const * file_name = NULL;

  // parse output grid
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(output_grid_node, &iter)))
    yaml_parse_output_grid_pair(&grid_name, &file_name, pair, yaml_filename);

  YAML_ASSERT(grid_name, "missing grid name");
  YAML_ASSERT_F(file_name, "missing file name for grid \"%s\"", grid_name);

  yac_couple_config_add_grid(couple_config, grid_name);
  yac_couple_config_grid_set_output_filename(
    couple_config, grid_name, file_name);
}

static void yaml_parse_debug_output_grids(
  struct yac_couple_config * couple_config, fy_node_t output_grids_node,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug_output_grids";

  YAML_ASSERT(
    fy_node_is_sequence(output_grids_node),
    "unsupported debug output grids node type "
    "(debug output grids is expected to be defined as a sequence)");

  // parse output grids
  void * iter = NULL;
  fy_node_t output_grid_node;
  while ((output_grid_node =
            fy_node_sequence_iterate(output_grids_node, &iter)))
    yaml_parse_output_grid(
      couple_config, output_grid_node, yaml_filename);
}

static void yaml_parse_debug_map_pair(
  struct yac_couple_config * couple_config, fy_node_pair_t debug_pair,
  const char * yaml_filename) {

  enum yaml_debug_key_types debug_key_type =
    (enum yaml_debug_key_types)
      yaml_parse_enum_value(
        yaml_debug_keys, num_yaml_debug_keys,
        fy_node_pair_key(debug_pair),
        "debug configuration parameter name", yaml_filename);

  fy_node_t value_node = fy_node_pair_value(debug_pair);

  switch (debug_key_type) {
    default:
    case(GLOBAL_CONFIG):
      yaml_parse_debug_global_config(
        couple_config, value_node, yaml_filename);
      break;
    case(OUTPUT_GRIDS):
      yaml_parse_debug_output_grids(
        couple_config, value_node, yaml_filename);
      break;
    case(MISSING_DEF):
      yac_couple_config_set_missing_definition_is_fatal(
        couple_config,
        yaml_parse_enum_value(
          bool_names, num_bool_names, value_node,
          "\"missing definition is fatal\" bool value", yaml_filename));
    break;
  };
}

static void yaml_parse_debug(
  struct yac_couple_config * couple_config,
  fy_node_t debug_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug";

  // check if the debug node is empty -> nothing to be read
  if (!debug_node) return;

  YAML_ASSERT(
    fy_node_is_mapping(debug_node),
    "unsupported debug node type "
    "(debug is expected to be defined as a mapping)");

  // parse couplings
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(debug_node, &iter)))
    yaml_parse_debug_map_pair(couple_config, pair, yaml_filename);
}

static void yaml_parse_base_map_pair(
  struct yac_couple_config * couple_config, fy_node_pair_t base_pair,
  const char * yaml_filename, enum yac_time_unit_type * time_unit,
  char const ** start_datetime, char const ** end_datetime) {
  char const * routine_name = "yaml_parse_base_map_pair";

  fy_node_t key_node = fy_node_pair_key(base_pair);

  YAML_ASSERT(
    fy_node_is_scalar(key_node),
    "unsupported key node type "
    "(key nodes are expected to be scalar nodes)");

  char const * base_key_name = fy_node_get_scalar0(key_node);
  int base_key_type =
    yac_name_type_pair_get_type(
      yaml_base_keys, num_yaml_base_keys, base_key_name);

  YAML_ASSERT_F(
    (base_key_type == INT_MAX) || // for unknown keys, which are skipped
    (base_key_type == START_DATE) ||
    (base_key_type == END_DATE) ||
    (base_key_type == CALENDAR) ||
    (base_key_type == TIMESTEP_UNIT) ||
    (base_key_type == COUPLING) ||
    (base_key_type == DEBUG),
    "unsupported base key name \"%s\"",
    base_key_name);

  fy_node_t value_node = fy_node_pair_value(base_pair);

  switch (base_key_type) {
    case (START_DATE):
      *start_datetime =
        yaml_parse_string_value(
          value_node, base_key_name, yaml_filename);
      break;
    case (END_DATE):
      *end_datetime =
        yaml_parse_string_value(
          value_node, base_key_name, yaml_filename);
      break;
    case (CALENDAR):
      yac_cdef_calendar(
        (int)yaml_parse_calendar_value(
          value_node, base_key_name, yaml_filename));
      break;
    case (TIMESTEP_UNIT): {
      enum yac_time_unit_type time_unit_ =
        yaml_parse_timestep_unit_value(
          value_node, base_key_name, yaml_filename);
      YAML_ASSERT(
        (*time_unit == TIME_UNIT_UNDEFINED) || (*time_unit == time_unit_),
        "inconsistent redefinition of time unit")
      *time_unit = time_unit_;
      break;
    }
    case (COUPLING):
      YAML_ASSERT(
        (*time_unit != TIME_UNIT_UNDEFINED),
        "time unit has to be defined before the couplings")
      yaml_parse_coupling(
        couple_config, value_node, yaml_filename, *time_unit);
      break;
    case (DEBUG):
      yaml_parse_debug(
        couple_config, value_node, yaml_filename);
    default:
      // nothing to be done
      break;
  }
}

static void yaml_parse_document(
  struct yac_couple_config * couple_config, fy_document_t document,
  const char * yaml_filename) {
  char const * routine_name = "yaml_parse_document";

  // get root node of document
  fy_node_t root_node = fy_document_root(document);

  // if the configuration file is empty
  if (!root_node) return;

  YAML_ASSERT(
    fy_node_is_mapping(root_node),
    "unsupported root node type (root node is expected to be a mapping node)");

  char const * start_datetime = NULL;
  char const * end_datetime = NULL;

  // parse base root mappings
  enum yac_time_unit_type time_unit = TIME_UNIT_UNDEFINED;
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(root_node, &iter)))
    yaml_parse_base_map_pair(
      couple_config, pair, yaml_filename, &time_unit,
      &start_datetime, &end_datetime);

  if ((start_datetime != NULL) || (end_datetime != NULL)) {

    YAC_ASSERT(
      yac_cget_calendar() != YAC_CALENDAR_NOT_SET,
      "ERROR(yaml_parse_document): "
      "cannot set start/end datetime because calendar has not yet been set");

    yac_couple_config_set_datetime(
      couple_config, start_datetime, end_datetime);
  }
}

void yac_yaml_read_coupling(
  struct yac_couple_config * couple_config, const char * yaml_filename,
  int parse_flags) {
  char const * routine_name = "yac_yaml_read_coupling";

  // check whether the yaml configuration file exists
  YAC_ASSERT_F(
    yac_file_exists(yaml_filename),
    "ERROR(%s): YAML configuration file could not be found \"%s\"",
    routine_name, yaml_filename);

  // open yaml configuration file
  FILE * config_file = xfopen(yaml_filename, "r");
  YAC_ASSERT_F(
    config_file,
    "ERROR(%s): could not open YAML configuration file \"%s\"",
    routine_name, yaml_filename);

  // parse yaml configuration file into document
  struct fy_parse_cfg parse_config =
    {.search_path = ".",
     .userdata = NULL,
     .diag = NULL};
  parse_config.flags =
    (enum fy_parse_cfg_flags)parse_flags;

  fy_document_t document =
    fy_document_build_from_fp(&parse_config, config_file);
  YAC_ASSERT_F(
    document,
    "ERROR(%s): could not parse YAML configuration file \"%s\"",
    routine_name, yaml_filename);

  // resolve anchors and merge keys
  YAML_ASSERT(
    !fy_document_resolve(document),
    "could not resolve anchors and merge keys");

  // parse document into couple configuration
  yaml_parse_document(couple_config, document, yaml_filename);

  // cleanup
  fy_document_destroy(document);
  xfclose(config_file);

  return;
}

struct yac_interp_stack_config *
  yac_yaml_parse_interp_stack_config_string(
    char const * str_interp_stack_config, int parse_flags) {
  char const * routine_name =
    "yac_yaml_parse_interp_stack_config_string";
  char const * yaml_filename =
    "user provided interp stack config string";

  YAML_ASSERT(
    str_interp_stack_config != NULL, "interpolation stack string is NULL");

  // parse string into document
  struct fy_parse_cfg parse_config =
    {.search_path = ".",
     .userdata = NULL,
     .diag = NULL};
  parse_config.flags =
    (enum fy_parse_cfg_flags)parse_flags;

  fy_document_t document =
    fy_document_build_from_string(
      &parse_config, str_interp_stack_config, (size_t)-1);
  YAML_ASSERT(document, "failed parsing");

  // resolve anchors and merge keys
  YAML_ASSERT(
    !fy_document_resolve(document),
    "could not resolve anchors and merge keys");

  fy_node_t interp_stack_config_node = fy_document_root(document);
  YAML_ASSERT(interp_stack_config_node, "invalid root node");

  struct yac_interp_stack_config * interp_stack_config =
    yac_interp_stack_config_new();

  yaml_parse_interp_stack_value(
    interp_stack_config, interp_stack_config_node,
    yaml_filename);

  // cleanup
  fy_document_destroy(document);

  return interp_stack_config;
}

static fy_node_t yac_yaml_create_scalar(
  fy_document_t document, char const * value) {

  // return NULL if value is empty
  if (!value) return (fy_node_t)NULL;

  fy_node_t scalar_node =
    fy_node_create_scalar_copy(document, value, strlen(value));
  YAC_ASSERT(
    scalar_node,
    "ERROR(yac_yaml_create_scalar): failed to create scalar node");

  return scalar_node;
}

static fy_node_t yac_yaml_create_scalar_int(
  fy_document_t document, int value) {

  char str_value[16];
  int value_size = snprintf(str_value, sizeof(str_value), "%d", value);
  YAC_ASSERT_F(
    (value_size >= 0) && ((size_t)value_size < sizeof(str_value)),
    "ERROR(yac_yaml_create_scalar_int): "
    "could not write \"%d\" to string buffer of size %zu",
    value, sizeof(str_value));

  return yac_yaml_create_scalar(document, str_value);
}

static fy_node_t yac_yaml_create_scalar_dble(
  fy_document_t document, double value) {

  char str_value[32];
  int value_size = snprintf(str_value, sizeof(str_value), "%g", value);
  YAC_ASSERT_F(
    (value_size >= 0) && ((size_t)value_size < sizeof(str_value)),
    "ERROR(yac_yaml_create_scalar_dble): "
    "could not write \"%g\" to string buffer of size %zu",
    value, sizeof(str_value));

  return yac_yaml_create_scalar(document, str_value);
}

static fy_node_t yac_yaml_create_sequence_scalar(
  fy_document_t document, char const * const * values, size_t num_values) {

  // return NULL if sequence is empty
  if (num_values == 0) return (fy_node_t)NULL;

  YAC_ASSERT(
    values, "ERROR(yac_yaml_create_sequence_scalar): no values provided");

  // create sequence node
  fy_node_t sequence_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    sequence_node, "ERROR(yac_yaml_create_sequence_scalar): "
    "failed to create sequence node");

  for (size_t value_idx = 0; value_idx < num_values; ++value_idx) {
    char const * value = values[value_idx];
    YAC_ASSERT_F(
      value, "ERROR(yac_yaml_create_sequence_scalar): "
      "invalid value at idx %zu", value_idx);
    int appending_failed =
      fy_node_sequence_append(
        sequence_node, yac_yaml_create_scalar(document, value));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_sequence_scalar): "
      "failed to append interpolation node");
  }

  return sequence_node;
}

static void yac_yaml_map_append(
  fy_node_t map, char const * key, fy_node_t value) {

  // if the value node is empty, return
  if (!value) return;

  YAC_ASSERT(
    key, "ERROR(yac_yaml_map_append): NULL key is not supported");

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document,
    "ERROR(yac_yaml_map_append): failed to get document from node");

  // set key and value for root node
  int appending_failed =
    fy_node_mapping_append(
      map, yac_yaml_create_scalar(document, key), value);
  YAC_ASSERT(
    !appending_failed,
    "ERROR(yac_yaml_map_append): failed to append mapping node pair");
}

static void yac_yaml_map_append_scalar(
  fy_node_t map, char const * key, char const * value) {

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_map_append_scalar): "
    "failed to get document from node");

  yac_yaml_map_append(map, key, yac_yaml_create_scalar(document, value));
}

static void yac_yaml_map_append_scalar_int(
  fy_node_t map, char const * key, int value) {

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_map_append_scalar_int): "
    "failed to get document from node");

  yac_yaml_map_append(
    map, key, yac_yaml_create_scalar_int(document, value));
}

static void yac_yaml_map_append_scalar_dble(
  fy_node_t map, char const * key, double value) {

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_map_append_scalar_dble): "
    "failed to get document from node");

  yac_yaml_map_append(
    map, key, yac_yaml_create_scalar_dble(document, value));
}

static fy_node_t yac_yaml_create_field_name_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  const char * src_field_name;
  const char * tgt_field_name;
  yac_couple_config_get_field_names(
    couple_config, couple_idx, field_couple_idx,
    &src_field_name, &tgt_field_name);

  // if both names are identical
  if (!strcmp(src_field_name, tgt_field_name))
    return yac_yaml_create_scalar(document, src_field_name);

  // create field name node
  fy_node_t field_name_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    field_name_node, "ERROR(yac_yaml_create_field_name_node): "
    "failed to create mapping node");

  // add source field name
  yac_yaml_map_append_scalar(field_name_node, "src", src_field_name);

  // add target field name
  yac_yaml_map_append_scalar(field_name_node, "tgt", tgt_field_name);

  return field_name_node;
}

static void yac_yaml_map_append_parameter(
  fy_node_t map,
  struct interp_method_parameter const * parameter,
  interp_method_parameter_value value) {

  YAC_ASSERT(
    (parameter->type == ENUM_PARAM) ||
    (parameter->type == INT_PARAM) ||
    (parameter->type == DBLE_PARAM) ||
    (parameter->type == BOOL_PARAM) ||
    (parameter->type == STR_PARAM) ||
    (parameter->type == DEG_PARAM),
    "ERROR(yac_yaml_map_append_parameter): unsupported parameter type");

  char const * key = parameter->name;

  switch (parameter->type) {
    default:
    case (ENUM_PARAM): {
      yac_yaml_map_append_scalar(
        map, key,
        yac_name_type_pair_get_name(
          parameter->data.enum_param.valid_values,
          parameter->data.enum_param.num_valid_values,
          value.enum_value));
      break;
    }
    case (INT_PARAM): {
      yac_yaml_map_append_scalar_int(map, key, value.int_value);
      break;
    }
    case (DEG_PARAM):
    case (DBLE_PARAM): {
      yac_yaml_map_append_scalar_dble(map, key, value.dble_value);
      break;
    }
    case (BOOL_PARAM): {
      yac_yaml_map_append_scalar(
        map, key,
        yac_name_type_pair_get_name(
          bool_names, num_bool_names, value.bool_value));
      break;
    }
    case (STR_PARAM): {
      yac_yaml_map_append_scalar(map, key, value.str_value);
      break;
    }
  };
}

static int compare_parameter_values(
  interp_method_parameter_value const * a,
  interp_method_parameter_value const * b,
  enum interp_method_parameter_value_type type) {

  YAC_ASSERT(
    (type == ENUM_PARAM) ||
    (type == INT_PARAM) ||
    (type == DBLE_PARAM) ||
    (type == BOOL_PARAM) ||
    (type == STR_PARAM) ||
    (type == DEG_PARAM),
    "ERROR(compare_parameter_values): "
    "invalid interpolation method parameter value type");

  switch (type) {
    default:
    case (ENUM_PARAM):
      return
        (a->enum_value > b->enum_value) - (a->enum_value < b->enum_value);
    case (INT_PARAM):
      return
        (a->int_value > b->int_value) - (a->int_value < b->int_value);
    case (DEG_PARAM):
    case (DBLE_PARAM):
      return
        (a->dble_value > b->dble_value) - (a->dble_value < b->dble_value);
    case (BOOL_PARAM):
      return
        (a->bool_value > b->bool_value) - (a->bool_value < b->bool_value);
    case (STR_PARAM):
      if ((a->str_value != NULL) && (b->str_value != NULL))
        return strcmp(a->str_value, b->str_value);
      else
        return (a->str_value != NULL) - (b->str_value != NULL);
  }
}

static fy_node_t yac_yaml_create_interpolation_node(
  fy_document_t document,
  union yac_interp_stack_config_entry const * interp_stack_entry) {

  enum yac_interpolation_list interp_type =
    yac_interp_stack_config_entry_get_type(interp_stack_entry);

  struct yac_interpolation_method const * interp_method = NULL;
  for (size_t i = 0;
       (i < NUM_INTERPOLATION_METHODS) && !interp_method; ++i)
    if (interpolation_methods[i].type == interp_type)
      interp_method = interpolation_methods + i;

  struct interp_method_parameter const * parameters =
    interp_method->parameters;
  size_t num_parameters = interp_method->num_parameters;

  interp_method_parameter_value
    parameter_values[MAX(num_parameters,1)];
  for (size_t param_idx = 0; param_idx < num_parameters; ++param_idx)
    parameter_values[param_idx] = parameters[param_idx].default_value;

  interp_method->get_interpolation(interp_stack_entry, parameter_values);

  // create parameter node
  fy_node_t parameter_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    parameter_node, "ERROR(yac_yaml_create_interpolation_node): "
    "failed to create mapping node");

  size_t num_non_default_parameters = 0;
  for (size_t param_idx = 0; param_idx < num_parameters; ++param_idx) {
    if (compare_parameter_values(
          &parameters[param_idx].default_value,
          parameter_values + param_idx, parameters[param_idx].type)) {
      ++num_non_default_parameters;
      yac_yaml_map_append_parameter(
        parameter_node, parameters + param_idx,
        parameter_values[param_idx]);
    }
  }

  fy_node_t interpolation_node;
  if (num_non_default_parameters) {

    // create interpolation node
    interpolation_node = fy_node_create_mapping(document);
    YAC_ASSERT(
      interpolation_node, "ERROR(yac_yaml_create_interpolation_node): "
      "failed to create mapping node");

    yac_yaml_map_append(
      interpolation_node, interp_method->name, parameter_node);
  } else {

    // free parameter node
    fy_node_free(parameter_node);

    interpolation_node =
      yac_yaml_create_scalar(document, interp_method->name);
  }

  return interpolation_node;
}

static fy_node_t yac_yaml_create_interpolation_stack_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  struct yac_interp_stack_config * interp_stack =
    yac_couple_config_get_interp_stack(
      couple_config, couple_idx, field_couple_idx);

  YAC_ASSERT(
    interp_stack, "ERROR(yac_yaml_create_interpolation_stack_node): "
    "invalid interpolation stack");

  // create interpolation stack node
  fy_node_t interp_stack_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    interp_stack_node, "ERROR(yac_yaml_create_interpolation_stack_node): "
    "failed to create sequence node");

  size_t interp_stack_size =
    yac_interp_stack_config_get_size(interp_stack);
  YAC_ASSERT(
    interp_stack_size, "ERROR(yac_yaml_create_interpolation_stack_node): "
    "invalid interpolation stack size");

  for (size_t interp_stack_idx = 0; interp_stack_idx < interp_stack_size;
       ++interp_stack_idx) {
    int appending_failed =
      fy_node_sequence_append(
        interp_stack_node,
        yac_yaml_create_interpolation_node(
          document,
          yac_interp_stack_config_get_entry(
            interp_stack, interp_stack_idx)));
    YAC_ASSERT(
      !appending_failed,
      "ERROR(yac_yaml_create_interpolation_stack_node): "
      "failed to append interpolation node");
  }

  return interp_stack_node;
}

static void yac_yaml_append_couple_field_nodes(
  fy_node_t coupling_node, struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  fy_document_t document = fy_node_document(coupling_node);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_append_couple_field_nodes): "
    "failed to get document from node");

  // create couple node
  fy_node_t field_couple_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    coupling_node, "ERROR(yac_yaml_append_couple_field_nodes): "
    "failed to create mapping node");

  // get component names
  char const * src_component_name;
  char const * tgt_component_name;
  yac_couple_config_get_field_couple_component_names(
    couple_config, couple_idx, field_couple_idx,
    &src_component_name, &tgt_component_name);

  // add source component name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SOURCE_COMPONENT),
    src_component_name);

  // add target component name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TARGET_COMPONENT),
    tgt_component_name);

  // get grid names
  char const * src_grid_name;
  char const * tgt_grid_name;
  yac_couple_config_get_field_grid_names(
    couple_config, couple_idx, field_couple_idx,
    &src_grid_name, &tgt_grid_name);

  // add source grid name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SOURCE_GRID),
    src_grid_name);

  // add target grid name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TARGET_GRID),
    tgt_grid_name);

  // add field names
  yac_yaml_map_append(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, FIELD),
    yac_yaml_create_field_name_node(
      document, couple_config, couple_idx, field_couple_idx));

  // add coupling period
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, COUPLING_PERIOD),
    yac_couple_config_get_coupling_period(
      couple_config, couple_idx, field_couple_idx));

  // add time reduction
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TIME_REDUCTION),
    yac_name_type_pair_get_name(
      time_operations, num_time_operations,
       yac_couple_config_get_coupling_period_operation(
         couple_config, couple_idx, field_couple_idx)));

  // add source lag
  yac_yaml_map_append_scalar_int(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SOURCE_LAG),
    yac_couple_config_get_source_lag(
      couple_config, couple_idx, field_couple_idx));

  // add target lag
  yac_yaml_map_append_scalar_int(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TARGET_LAG),
    yac_couple_config_get_target_lag(
      couple_config, couple_idx, field_couple_idx));

  // add weight file name
  if (yac_couple_config_enforce_write_weight_file(
        couple_config, couple_idx, field_couple_idx))
      yac_yaml_map_append_scalar(
        field_couple_node,
        yac_name_type_pair_get_name(
          yaml_couple_keys, num_yaml_couple_keys, WEIGHT_FILE_NAME),
        yac_couple_config_get_weight_file_name(
          couple_config, couple_idx, field_couple_idx));

  // add mapping side
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, MAPPING_SIDE),
    yac_name_type_pair_get_name(
      mapping_sides, num_mapping_sides,
       yac_couple_config_mapping_on_source(
         couple_config, couple_idx, field_couple_idx)));

  // add scale factor
  yac_yaml_map_append_scalar_dble(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SCALE_FACTOR),
    yac_couple_config_get_scale_factor(
      couple_config, couple_idx, field_couple_idx));

  // add scale summand
  yac_yaml_map_append_scalar_dble(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SCALE_SUMMAND),
    yac_couple_config_get_scale_summand(
      couple_config, couple_idx, field_couple_idx));

  // add interpolation
  yac_yaml_map_append(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, INTERPOLATION),
    yac_yaml_create_interpolation_stack_node(
      document, couple_config, couple_idx, field_couple_idx));

  // add source mask names
  char const * const * src_mask_names;
  size_t num_src_mask_names;
  yac_couple_config_get_src_mask_names(
    couple_config, couple_idx, field_couple_idx,
    &src_mask_names, &num_src_mask_names);
  if (num_src_mask_names == 1)
    yac_yaml_map_append_scalar(
      field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, SOURCE_MASK_NAME),
        src_mask_names[0]);
    else if (num_src_mask_names > 1)
      yac_yaml_map_append(
        field_couple_node,
        yac_name_type_pair_get_name(
          yaml_couple_keys, num_yaml_couple_keys, SOURCE_MASK_NAMES),
          yac_yaml_create_sequence_scalar(
            document, src_mask_names, num_src_mask_names));

  // add target mask name
  char const * tgt_mask_name =
    yac_couple_config_get_tgt_mask_name(
      couple_config, couple_idx, field_couple_idx);
  if (tgt_mask_name)
    yac_yaml_map_append_scalar(
      field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, TARGET_MASK_NAME),
      tgt_mask_name);

  int appending_failed =
    fy_node_sequence_append(coupling_node, field_couple_node);
  YAC_ASSERT(
    !appending_failed,
    "ERROR(yac_yaml_append_couple_field_nodes): "
    "failed to append field couple node");
}

static void yac_yaml_append_couple_nodes(
  fy_node_t coupling_node, struct yac_couple_config * couple_config,
  size_t couple_idx) {

  size_t num_couple_fields =
    yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

  for (size_t field_couple_idx = 0;
       field_couple_idx < num_couple_fields; ++field_couple_idx)
    yac_yaml_append_couple_field_nodes(
      coupling_node, couple_config, couple_idx, field_couple_idx);
}

static fy_node_t yac_yaml_create_output_grid_node(
  fy_document_t document, char const * grid_name, char const * file_name) {

  fy_node_t output_grid_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    output_grid_node, "ERROR(yac_yaml_create_output_grid_node): "
    "failed to create mapping node");

  yac_yaml_map_append_scalar(
    output_grid_node,
    yac_name_type_pair_get_name(
      yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
      OUTPUT_GRID_GRID_NAME),
    grid_name);
  yac_yaml_map_append_scalar(
    output_grid_node,
    yac_name_type_pair_get_name(
      yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
      OUTPUT_GRID_FILE_NAME),
    file_name);

  return output_grid_node;
}

static fy_node_t yac_yaml_create_output_grids_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // count the number of output grids
  size_t num_output_grids = 0;
  size_t num_grids =
    yac_couple_config_get_num_grids(couple_config);
  for (size_t grid_idx = 0; grid_idx < num_grids; ++grid_idx)
    if (yac_couple_config_grid_get_output_filename(
          couple_config,
          yac_couple_config_get_grid_name(couple_config, grid_idx)) != NULL)
      ++num_output_grids;

  fy_node_t output_grids_node = NULL;

  if (num_output_grids > 0) {

    // create output grids node
    output_grids_node = fy_node_create_sequence(document);
    YAC_ASSERT(
      output_grids_node, "ERROR(yac_yaml_create_output_grids_node): "
      "failed to create sequence node");

    // for all output grids
    for (size_t grid_idx = 0; grid_idx < num_grids; ++grid_idx) {

      char const * grid_name =
        yac_couple_config_get_grid_name(couple_config, grid_idx);
      char const * file_name =
        yac_couple_config_grid_get_output_filename(couple_config, grid_name);

      if (file_name != NULL) {
        int appending_failed =
          fy_node_sequence_append(
            output_grids_node,
            yac_yaml_create_output_grid_node(document, grid_name, file_name));
        YAC_ASSERT(
          !appending_failed, "ERROR(yac_yaml_create_output_grids_node): "
          "failed to append output grid node");
      }
    }
  }

  return output_grids_node;
}

static fy_node_t yac_yaml_create_debug_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // create debug node
  fy_node_t debug_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    debug_node, "ERROR(yac_yaml_create_debug_node): "
    "failed to create mapping node");

  // add output grids node
  yac_yaml_map_append(
    debug_node,
    yac_name_type_pair_get_name(
      yaml_debug_keys, num_yaml_debug_keys, OUTPUT_GRIDS),
    yac_yaml_create_output_grids_node(document, couple_config));

  // add "missing_definition_is_fatal" node
  yac_yaml_map_append_scalar(
    debug_node,
    yac_name_type_pair_get_name(
      yaml_debug_keys, num_yaml_debug_keys, MISSING_DEF),
    yac_name_type_pair_get_name(
      bool_names, num_bool_names,
      yac_couple_config_get_missing_definition_is_fatal(couple_config)));

  return debug_node;
}

static fy_node_t yac_yaml_create_coupling_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  size_t num_couples =
    yac_couple_config_get_num_couples(couple_config);

  if (!num_couples) return NULL;

  // create coupling node
  fy_node_t coupling_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    coupling_node, "ERROR(yac_yaml_create_coupling_node): "
    "failed to create sequence node");

  // for all couples
  for (size_t couple_idx = 0; couple_idx < num_couples; ++couple_idx)
    yac_yaml_append_couple_nodes(
      coupling_node, couple_config, couple_idx);

  return coupling_node;
}

static fy_node_t yac_yaml_create_field_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t component_idx, size_t field_idx) {

  fy_node_t field_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    field_node, "ERROR(yac_yaml_create_field_node): "
    "failed to create mapping node");

  // get field parameters
  char const * component_name =
    yac_couple_config_get_component_name(
      couple_config, component_idx);
  char const * field_name =
    yac_couple_config_get_field_name(
      couple_config, component_idx, field_idx);
  char const * grid_name =
    yac_couple_config_get_field_grid_name(
      couple_config, component_idx, field_idx);
  char const * metadata =
    yac_couple_config_field_get_metadata(
      couple_config, component_name, grid_name, field_name);
  size_t collection_size =
    yac_couple_config_get_field_collection_size(
      couple_config, component_name, grid_name, field_name);
  char const * timestep =
    yac_couple_config_get_field_timestep(
      couple_config, component_name, grid_name, field_name);
  int role =
    yac_couple_config_get_field_role(
      couple_config, component_name, grid_name, field_name);
  double frac_mask_fallback_value =
    yac_couple_config_get_frac_mask_fallback_value(
      couple_config, component_name, grid_name, field_name);

  // add field name
  yac_yaml_map_append_scalar(field_node, "name", field_name);

  // add grid name
  yac_yaml_map_append_scalar(field_node, "grid_name", grid_name);

  // add metadata
  if (metadata)
    yac_yaml_map_append_scalar(field_node, "metadata", metadata);

  // add collection_size
  if (collection_size != SIZE_MAX)
    yac_yaml_map_append_scalar_int(
      field_node, "collection_size", (int)collection_size);

  // add timestep
  if (timestep)
    yac_yaml_map_append_scalar(field_node, "timestep", timestep);

  // add role
  yac_yaml_map_append_scalar(
    field_node, "role",
    yac_name_type_pair_get_name(
      role_types, num_role_types, (enum yac_field_exchange_type)role));

  // add fractional fallback value
  if (frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE)
    yac_yaml_map_append_scalar_dble(
      field_node, "frac_mask_fallback_value", frac_mask_fallback_value);

  return field_node;
}

static fy_node_t yac_yaml_create_fields_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t component_idx) {

  size_t num_fields =
    yac_couple_config_get_num_fields(couple_config, component_idx);

  if (num_fields == 0) return (fy_node_t)NULL;

  fy_node_t fields_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    fields_node, "ERROR(yac_yaml_create_fields_node): "
    "failed to create sequence node");

  for (size_t field_idx = 0; field_idx < num_fields; ++field_idx) {

    int appending_failed =
      fy_node_sequence_append(
        fields_node,
        yac_yaml_create_field_node(
          document, couple_config, component_idx, field_idx));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_fields_node): "
      "failed to append field node");
  }

  return fields_node;
}

static fy_node_t yac_yaml_create_component_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t component_idx) {

  fy_node_t component_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    component_node, "ERROR(yac_yaml_create_component_node): "
    "failed to create mapping node");

  // get component name, component metadata, and fields
  char const * component_name =
    yac_couple_config_get_component_name(couple_config, component_idx);
  char const * metadata =
    yac_couple_config_component_get_metadata(couple_config, component_name);
  fy_node_t fields_node =
    yac_yaml_create_fields_node(
      document, couple_config, component_idx);

  // add component name
  yac_yaml_map_append_scalar(component_node, "name", component_name);

  // add metadata
  if (metadata)
    yac_yaml_map_append_scalar(component_node, "metadata", metadata);

  // add fields
  yac_yaml_map_append(component_node, "fields", fields_node);

  return component_node;
}

static fy_node_t yac_yaml_create_components_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // get number of components in coupling configuration
  size_t num_components = yac_couple_config_get_num_components(couple_config);

  if (num_components == 0) return (fy_node_t)NULL;

  // create sequence node
  fy_node_t components_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    components_node, "ERROR(yac_yaml_create_components_node): "
    "failed to create sequence node");

  for (size_t component_idx = 0; component_idx < num_components;
      ++component_idx) {

    int appending_failed =
      fy_node_sequence_append(
        components_node,
        yac_yaml_create_component_node(document, couple_config, component_idx));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_components_node): "
      "failed to append component node");
  }

  return components_node;
}

static fy_node_t yac_yaml_create_grid_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t grid_idx) {

  fy_node_t grid_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    grid_node, "ERROR(yac_yaml_create_grid_node): "
    "failed to create mapping node");

  // get grid name and metadata
  char const * grid_name =
    yac_couple_config_get_grid_name(couple_config, grid_idx);
  char const * metadata =
    yac_couple_config_grid_get_metadata(couple_config, grid_name);

  // add grid name
  yac_yaml_map_append_scalar(grid_node, "name", grid_name);

  // add metadata
  if (metadata)
    yac_yaml_map_append_scalar(grid_node, "metadata", metadata);

  return grid_node;
}

static fy_node_t yac_yaml_create_grids_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // get number of grids in coupling configuration
  size_t num_grids = yac_couple_config_get_num_grids(couple_config);

  if (num_grids == 0) return (fy_node_t)NULL;

  // create sequence node
  fy_node_t grids_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    grids_node, "ERROR(yac_yaml_create_grids_node): "
    "failed to create sequence node");

  for (size_t grids_idx = 0; grids_idx < num_grids; ++grids_idx) {

    int appending_failed =
      fy_node_sequence_append(
        grids_node,
        yac_yaml_create_grid_node(document, couple_config, grids_idx));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_grids_node): "
      "failed to append grid node");
  }

  return grids_node;
}

static fy_node_t yac_yaml_create_definitions_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // create definition
  fy_node_t definition_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    definition_node,
    "ERROR(yac_yaml_create_definitions_node): "
    "failed to create mapping node");

  // add components
  yac_yaml_map_append(
    definition_node,
    "components", yac_yaml_create_components_node(document, couple_config));

  // add grids
  yac_yaml_map_append(
    definition_node, "grids", yac_yaml_create_grids_node(document, couple_config));

  return definition_node;
}

static fy_node_t yac_yaml_create_couple_config_nodes(
  fy_document_t document, struct yac_couple_config * couple_config,
  int include_definitions) {

  // create root node
  fy_node_t root_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    root_node,
    "ERROR(yac_yaml_create_couple_config_nodes): "
    "failed to create mapping node");

  // add debug
  yac_yaml_map_append(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, DEBUG),
    yac_yaml_create_debug_node(document, couple_config));

  // add user definitions (components, grids, and fields)
  if (include_definitions)
    yac_yaml_map_append(
      root_node, "definitions",
      yac_yaml_create_definitions_node(document, couple_config));

  // add start datetime
  char * start_datetime = yac_couple_config_get_start_datetime(couple_config);
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, START_DATE), start_datetime);
  free(start_datetime);

  // add end datetime
  char * end_datetime = yac_couple_config_get_end_datetime(couple_config);
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, END_DATE), end_datetime);
  free(end_datetime);

  // add calendar
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, CALENDAR),
    yac_name_type_pair_get_name(
      calendar_types, num_calendar_types, getCalendarType()));

  // add timestep unit
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, TIMESTEP_UNIT),
    yac_name_type_pair_get_name(
      timestep_units, num_timestep_units, C_ISO_FORMAT));

  // add couplings
  yac_yaml_map_append(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, COUPLING),
    yac_yaml_create_coupling_node(document, couple_config));

  return root_node;
}

char * yac_yaml_emit_coupling(
  struct yac_couple_config * couple_config, int emit_flags,
  int include_definitions) {

  // create an empty document
  fy_document_t document = fy_document_create(NULL);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_emit): failed to create document");

  // create nodes from coupling configuration
  fy_node_t root_node =
    yac_yaml_create_couple_config_nodes(
      document, couple_config, include_definitions);

  // set root node of the document
  int setting_root_failed =
    fy_document_set_root(document, root_node);
  YAC_ASSERT(
    !setting_root_failed,
    "ERROR(yac_yaml_emit): failed to add root node to document");

  // emit document to string
  char * str_document =
    fy_emit_document_to_string(
      document, (enum fy_emitter_cfg_flags)emit_flags);

  YAC_ASSERT(
    str_document, "ERROR(yac_yaml_emit): failed to emit document to string");

  // destroy document
  fy_document_destroy(document);

  return str_document;
}
