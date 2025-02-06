// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_STACK_CONFIG_H
#define INTERP_STACK_CONFIG_H

#include "interp_method_avg.h"
#include "interp_method_check.h"
#include "interp_method_conserv.h"
#include "interp_method_creep.h"
#include "interp_method_ncc.h"
#include "interp_method_nnn.h"
#include "interp_method_fixed.h"
#include "interp_method_file.h"
#include "interp_method_hcsbb.h"
#include "interp_method_spmap.h"
#include "interp_method_callback.h"

// YAC PUBLIC HEADER START

#define YAC_MAX_ROUTINE_NAME_LENGTH (256)
#define YAC_MAX_FILE_NAME_LENGTH (512)

enum yac_interpolation_list {
   YAC_UNDEFINED             = 0,
   YAC_AVERAGE               = 1,
   YAC_N_NEAREST_NEIGHBOR    = 2,
   YAC_CONSERVATIVE          = 3,
   YAC_SOURCE_TO_TARGET_MAP  = 4,
   YAC_FIXED_VALUE           = 5,
   YAC_USER_FILE             = 6,
   YAC_CHECK                 = 7,
   YAC_BERNSTEIN_BEZIER      = 8,
   YAC_RADIAL_BASIS_FUNCTION = 9,
   YAC_CREEP                 = 10,
   YAC_USER_CALLBACK         = 11,
   YAC_NEAREST_CORNER_CELLS  = 12,
};

/** \example test_interp_stack_config.c
 * Tests for interpolation stack interface routines.
 */

struct yac_interp_stack_config;
union yac_interp_stack_config_entry;

struct yac_interp_stack_config * yac_interp_stack_config_new();
void yac_interp_stack_config_delete(
  struct yac_interp_stack_config * interp_stack_config);
struct yac_interp_stack_config * yac_interp_stack_config_copy(
  struct yac_interp_stack_config * interp_stack);

void yac_interp_stack_config_add_average(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_avg_weight_type reduction_type, int partial_coverage);
void yac_interp_stack_config_add_ncc(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_ncc_weight_type type, int partial_coverage);
void yac_interp_stack_config_add_nnn(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_nnn_weight_type type, size_t n,
  double max_search_distance, double scale);
void yac_interp_stack_config_add_conservative(
  struct yac_interp_stack_config * interp_stack_config,
  int order, int enforced_conserv, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation);
void yac_interp_stack_config_add_spmap(
  struct yac_interp_stack_config * interp_stack_config,
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  enum yac_interp_spmap_scale_type scale_type,
  double src_sphere_radius, double tgt_sphere_radius);
void yac_interp_stack_config_add_hcsbb(
  struct yac_interp_stack_config * interp_stack_config);
void yac_interp_stack_config_add_user_file(
  struct yac_interp_stack_config * interp_stack_config, char const * filename);
void yac_interp_stack_config_add_fixed(
  struct yac_interp_stack_config * interp_stack_config, double value);
void yac_interp_stack_config_add_check(
  struct yac_interp_stack_config * interp_stack_config,
  char const * constructor_key, char const * do_search_key);
void yac_interp_stack_config_add_creep(
  struct yac_interp_stack_config * interp_stack_config, int creep_distance);
void yac_interp_stack_config_add_user_callback(
  struct yac_interp_stack_config * interp_stack_config,
  char const * func_compute_weights_key);

int yac_interp_stack_config_compare(void const * a, void const * b);
struct interp_method ** yac_interp_stack_config_generate(
  struct yac_interp_stack_config * interp_stack);

size_t yac_interp_stack_config_get_pack_size(
  struct yac_interp_stack_config * interp_stack, MPI_Comm comm);
void yac_interp_stack_config_pack(
  struct yac_interp_stack_config * interp_stack,
  void * buffer, int buffer_size, int * position, MPI_Comm comm);
struct yac_interp_stack_config * yac_interp_stack_config_unpack(
  void * buffer, int buffer_size, int * position, MPI_Comm comm);

size_t yac_interp_stack_config_get_size(
  struct yac_interp_stack_config * interp_stack);
union yac_interp_stack_config_entry const *
  yac_interp_stack_config_get_entry(
    struct yac_interp_stack_config * interp_stack,
    size_t interp_stack_idx);

enum yac_interpolation_list yac_interp_stack_config_entry_get_type(
  union yac_interp_stack_config_entry const * interp_stack_entry);

void yac_interp_stack_config_entry_get_average(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_avg_weight_type * reduction_type,
  int * partial_coverage);
void yac_interp_stack_config_entry_get_ncc(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_ncc_weight_type * type, int * partial_coverage);
void yac_interp_stack_config_entry_get_nnn(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_nnn_weight_type * type, size_t * n,
  double * max_search_distance, double * scale);
void yac_interp_stack_config_entry_get_conservative(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * order, int * enforced_conserv, int * partial_coverage,
  enum yac_interp_method_conserv_normalisation * normalisation);
void yac_interp_stack_config_entry_get_spmap(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  double * spread_distance, double * max_search_distance,
  enum yac_interp_spmap_weight_type * weight_type,
  enum yac_interp_spmap_scale_type * scale_type,
  double * src_sphere_radius, double * tgt_sphere_radius);
void yac_interp_stack_config_entry_get_user_file(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** filename);
void yac_interp_stack_config_entry_get_fixed(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  double * value);
void yac_interp_stack_config_entry_get_check(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** constructor_key, char const ** do_search_key);
void yac_interp_stack_config_entry_get_creep(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * creep_distance);
void yac_interp_stack_config_entry_get_user_callback(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** func_compute_weights_key);

// YAC PUBLIC HEADER STOP

#endif // INTERP_STACK_CONFIG_H
