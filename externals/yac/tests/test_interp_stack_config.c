/**
 * @file test_interp_stack_config.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include "tests.h"
#include "interp_stack_config.h"

#define ARGS(...) __VA_ARGS__
#define _GET_NTH_ARG(_1, _2, _3, _4, N, ...) N
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
#define CHECK_STACKS(interp_name, ... ) \
  _CHECK_STACKS(interp_name, \
    EXPAND(_GET_NTH_ARG(__VA_ARGS__, _CONFIG_ARGS4, \
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
    yac_interp_stack_config_add_average(a, AVG_ARITHMETIC, 1);
    yac_interp_stack_config_add_fixed(a, -1.0);
    yac_interp_stack_config_add_average(a, AVG_ARITHMETIC, 1);
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
    ARGS(AVG_ARITHMETIC, AVG_DIST, AVG_BARY),
    FOREACH_BOOL(
      partial_coverage,
      CHECK_STACKS(average, interp_avg_weight_type, partial_coverage)))

  // compare nnn config
  //  for NNN_AVG and NNN_DIST the scale parameter is being ignored
  FOREACH_ENUM(
    interp_nnn_weight_type,
    ARGS(NNN_AVG, NNN_DIST),
    FOREACH_INT(
      counts, ARGS(1,3,9),
      FOREACH_DBLE(
        scales, -1.0,
        CHECK_STACKS(nnn, interp_nnn_weight_type, counts, scales))))

  // compare nnn config
  //   for NNN_GAUSS and NNN_RBF the scale parameter is being interpreted
  FOREACH_ENUM(
    interp_nnn_weight_type,
    ARGS(NNN_GAUSS, NNN_RBF),
    FOREACH_INT(
      counts, ARGS(1,3,9),
      FOREACH_DBLE(
        scales, ARGS(0.5, 1.0),
        CHECK_STACKS(nnn, interp_nnn_weight_type, counts, scales))))

  // compare conservative config
  FOREACH_INT(
    order, ARGS(1,2),
    FOREACH_BOOL(
      enforced_conserv,
      FOREACH_BOOL(
        partial_coverage,
        FOREACH_ENUM(
          interp_method_conserv_normalisation,
          ARGS(CONSERV_DESTAREA, CONSERV_FRACAREA),
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
        ARGS(SPMAP_AVG, SPMAP_DIST),
        CHECK_STACKS(spmap,
          spread_distance, max_search_distance, interp_spmap_weight_type))))

  // compare user file
  FOREACH_STRING(
    filename, ARGS("file_a.nc", "file_b.nc"),
    FOREACH_STRING(
      src_grid_name, ARGS("src_grid1", "src_grid2"),
      FOREACH_STRING(
        tgt_grid_name, ARGS("tgt_grid1", "tgt_grid2"),
        CHECK_STACKS(user_file,
          filename, src_grid_name, tgt_grid_name))))

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
