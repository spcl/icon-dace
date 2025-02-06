// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <math.h>
#include "tests.h"
#include "interp_stack_config.h"

#define ARGS(...) __VA_ARGS__
#define _GET_NTH_ARG(_1, _2, _3, _4, _5, _6, N, ...) N
#define EXPAND(x) x
#define FOREACH(name, ...) \
  { \
    enum {NUM_ ## name = sizeof( name ) / sizeof( name [0])}; \
    int name ## _idx[2]; \
    for (name ## _idx[0] = 0; name ## _idx[0] < NUM_ ## name; \
         ++ name ## _idx[0]) { \
      for (name ## _idx[1] = 0; name ## _idx[1] < NUM_ ## name; \
           ++ name ## _idx[1]) { \
        configs_differ += (name ## _idx[0]) != (name ## _idx[1]); \
        {__VA_ARGS__} \
        configs_differ -= (name ## _idx[0]) != (name ## _idx[1]); \
      } \
    } \
  }
#define FOREACH_ENUM(name, values, ...) \
  { \
    enum yac_ ## name name [] = {values}; \
    FOREACH(name, __VA_ARGS__) \
  }
#define FOREACH_TYPE(name, type, values, ...) \
  { \
    type name[] = {values}; \
    FOREACH(name, __VA_ARGS__) \
  }
#define FOREACH_INT(name, values, ...) \
  FOREACH_TYPE(name, int, ARGS(values), __VA_ARGS__)
#define FOREACH_DBLE(name, values, ...) \
  FOREACH_TYPE(name, double, ARGS(values), __VA_ARGS__)
#define FOREACH_BOOL(name, ...) FOREACH_INT(name, ARGS(0, 1), __VA_ARGS__)
#define FOREACH_STRING(name, values, ...) \
  FOREACH_TYPE(name, ARGS(char const *), ARGS(values), __VA_ARGS__)
#define _CHECK_STACKS(interp_name, config) \
  { \
    int config_idx; \
    struct yac_interp_stack_config * a = yac_interp_stack_config_new(); \
    struct yac_interp_stack_config * b = yac_interp_stack_config_new(); \
    config_idx = 0, yac_interp_stack_config_add_ ## interp_name ( a, config ); \
    config_idx = 1, yac_interp_stack_config_add_ ## interp_name ( b, config ); \
    check_compare_stacks(a, b, configs_differ); \
  }
#define _CONFIG_ARGS1(arg_name) arg_name[arg_name ## _idx[config_idx]]
#define _CONFIG_ARGS2(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS1(__VA_ARGS__)
#define _CONFIG_ARGS3(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS2(__VA_ARGS__)
#define _CONFIG_ARGS4(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS3(__VA_ARGS__)
#define _CONFIG_ARGS5(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS4(__VA_ARGS__)
#define _CONFIG_ARGS6(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS5(__VA_ARGS__)
#define CHECK_STACKS(interp_name, ... ) \
  _CHECK_STACKS(interp_name, \
    EXPAND(_GET_NTH_ARG(__VA_ARGS__, _CONFIG_ARGS6, \
                                     _CONFIG_ARGS5, \
                                     _CONFIG_ARGS4, \
                                     _CONFIG_ARGS3, \
                                     _CONFIG_ARGS2, \
                                     _CONFIG_ARGS1)(__VA_ARGS__)))

static void check_compare_stacks(
  struct yac_interp_stack_config * a, struct yac_interp_stack_config * b,
  int configs_differ);

int main (void) {

  int configs_differ = 0;

  { // stack with different sizes
    struct yac_interp_stack_config * a = yac_interp_stack_config_new();
    struct yac_interp_stack_config * b = yac_interp_stack_config_new();
    yac_interp_stack_config_add_average(a, YAC_INTERP_AVG_ARITHMETIC, 1);
    yac_interp_stack_config_add_fixed(a, -1.0);
    yac_interp_stack_config_add_average(a, YAC_INTERP_AVG_ARITHMETIC, 1);
    check_compare_stacks(a, b, 1);
  }

  { // compare empty config
    struct yac_interp_stack_config * a = yac_interp_stack_config_new();
    struct yac_interp_stack_config * b = yac_interp_stack_config_new();
    check_compare_stacks(a, b, 0);
  }

  // compare average config
  FOREACH_ENUM(
    interp_avg_weight_type,
    ARGS(YAC_INTERP_AVG_ARITHMETIC, YAC_INTERP_AVG_DIST, YAC_INTERP_AVG_BARY),
    FOREACH_BOOL(
      partial_coverage,
      CHECK_STACKS(average, interp_avg_weight_type, partial_coverage)))

  // compare ncc config
  FOREACH_ENUM(
    interp_ncc_weight_type,
    ARGS(YAC_INTERP_NCC_AVG, YAC_INTERP_NCC_DIST),
    FOREACH_BOOL(
      partial_coverage,
      CHECK_STACKS(ncc, interp_ncc_weight_type, partial_coverage)))

  // compare nnn config
  //  for YAC_INTERP_NNN_AVG, YAC_INTERP_NNN_DIST, and YAC_INTERP_NNN_ZERO
  //  the scale parameter is being ignored
  FOREACH_ENUM(
    interp_nnn_weight_type,
    ARGS(YAC_INTERP_NNN_AVG, YAC_INTERP_NNN_DIST, YAC_INTERP_NNN_ZERO),
    FOREACH_INT(
      counts, ARGS(1,3,9),
      FOREACH_DBLE(
        max_search_distance, ARGS(0.0, M_PI_2),
        FOREACH_DBLE(
          scales, -1.0,
          CHECK_STACKS(
            nnn, interp_nnn_weight_type, counts,
            max_search_distance, scales)))))

  // compare nnn config
  //   for YAC_INTERP_NNN_GAUSS and YAC_INTERP_NNN_RBF the scale
  //   parameter is being interpreted
  FOREACH_ENUM(
    interp_nnn_weight_type,
    ARGS(YAC_INTERP_NNN_GAUSS, YAC_INTERP_NNN_RBF),
    FOREACH_INT(
      counts, ARGS(1,3,9),
      FOREACH_DBLE(
        max_search_distance, ARGS(0.0, M_PI_2),
        FOREACH_DBLE(
          scales, ARGS(0.5, 1.0),
          CHECK_STACKS(
            nnn, interp_nnn_weight_type, counts,
            max_search_distance, scales)))))

  // compare conservative config
  FOREACH_INT(
    order, ARGS(1,2),
    FOREACH_BOOL(
      enforced_conserv,
      FOREACH_BOOL(
        partial_coverage,
        FOREACH_ENUM(
          interp_method_conserv_normalisation,
          ARGS(YAC_INTERP_CONSERV_DESTAREA, YAC_INTERP_CONSERV_FRACAREA),
          CHECK_STACKS(conservative,
            order, enforced_conserv, partial_coverage,
            interp_method_conserv_normalisation)))))

  // compare source point mapping
  FOREACH_DBLE(
    spread_distance, ARGS(0.0, 1.0, 2.0),
    FOREACH_DBLE(
      max_search_distance, ARGS(0.0, 1.0, 2.0),
      FOREACH_ENUM(
        interp_spmap_weight_type,
        ARGS(YAC_INTERP_SPMAP_AVG, YAC_INTERP_SPMAP_DIST),
        FOREACH_ENUM(
          interp_spmap_scale_type,
          ARGS(
            YAC_INTERP_SPMAP_NONE, YAC_INTERP_SPMAP_SRCAREA,
            YAC_INTERP_SPMAP_INVTGTAREA, YAC_INTERP_SPMAP_FRACAREA),
          FOREACH_DBLE(
            src_sphere_radius, ARGS(1.0, 6371000.0),
            FOREACH_DBLE(
              tgt_sphere_radius, ARGS(1.0, 6371000.0),
              CHECK_STACKS(spmap,
                spread_distance, max_search_distance,
                interp_spmap_weight_type, interp_spmap_scale_type,
                src_sphere_radius, tgt_sphere_radius)))))))

  // compare user file
  FOREACH_STRING(
    filename, ARGS(
      "test_interp_stack_config_file_a.nc",
      "test_interp_stack_config_file_b.nc"),
    CHECK_STACKS(user_file, filename))

  // compare fixes
  FOREACH_DBLE(
    fixed_value, ARGS(-1.0, 0.0, 1.0),
    CHECK_STACKS(fixed, fixed_value))

  // compare check
  FOREACH_STRING(
    constructor_key, ARGS(NULL, "constructor_a", "constructor_b"),
    FOREACH_STRING(
      do_search_key, ARGS(NULL, "do_search_key_a", "do_search_key_b"),
      CHECK_STACKS(check, constructor_key, do_search_key)))

  // compare creep
  FOREACH_INT(
    creep_distance, ARGS(-1, 0, 1),
    CHECK_STACKS(creep, creep_distance))

  // compare user callback
  FOREACH_STRING(
    compute_weights_key, ARGS("compute_weights_a", "compute_weights_b"),
    CHECK_STACKS(user_callback, compute_weights_key))

  return TEST_EXIT_CODE;
}

static void check_compare_stacks(
  struct yac_interp_stack_config * a, struct yac_interp_stack_config * b,
  int configs_differ) {

  configs_differ = configs_differ != 0;

  if (yac_interp_stack_config_compare(a, a))
    PUT_ERR("error in yac_interp_stack_config_compare (a != a)")
  if (yac_interp_stack_config_compare(b, b))
    PUT_ERR("error in yac_interp_stack_config_compare (b != b)")

  int cmp_a = yac_interp_stack_config_compare(a, b);
  int cmp_b = yac_interp_stack_config_compare(b, a);

  if ((cmp_a != cmp_b) ^ configs_differ)
    PUT_ERR("error in yac_interp_stack_config_compare ((a > b) == (a < b))")
  if ((cmp_a != 0) ^ configs_differ)
    PUT_ERR("error in yac_interp_stack_config_compare ((a > b) == 0)")
  if ((cmp_b != 0) ^ configs_differ)
    PUT_ERR("error in yac_interp_stack_config_compare ((a > b) == 0)")

  yac_interp_stack_config_delete(b);
  yac_interp_stack_config_delete(a);
}
