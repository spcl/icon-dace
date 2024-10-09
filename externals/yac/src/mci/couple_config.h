// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

/** \example test_couple_config.c
 * This example tests the generation of a coupling configuration.
 */

#ifndef COUPLE_CONFIG_H
#define COUPLE_CONFIG_H

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

enum yac_time_unit_type {
   C_MILLISECOND = 0,
   C_SECOND      = 1,
   C_MINUTE      = 2,
   C_HOUR        = 3,
   C_DAY         = 4,
   C_MONTH       = 5,
   C_YEAR        = 6,
   C_ISO_FORMAT  = 7,
   TIME_UNIT_UNDEFINED,
};

enum yac_reduction_type {
  TIME_NONE       = 0,
  TIME_ACCUMULATE = 1,
  TIME_AVERAGE    = 2,
  TIME_MINIMUM    = 3,
  TIME_MAXIMUM    = 4,
};

enum yac_text_filetype {
  YAC_TEXT_FILETYPE_YAML   = 0, //!< YAML format
  YAC_TEXT_FILETYPE_JSON   = 1, //!< JSON format
};

#define MISSING_DEFINITION_IS_FATAL_DEFAULT_VALUE (1)

struct yac_couple_config;

// general couple config routines
struct yac_couple_config * yac_couple_config_new();
void yac_couple_config_delete(struct yac_couple_config * couple_config);

// basic run parameters
char * yac_couple_config_get_start_datetime(
  struct yac_couple_config * couple_config);
char * yac_couple_config_get_end_datetime(
  struct yac_couple_config * couple_config);
void yac_couple_config_set_datetime(
  struct yac_couple_config * couple_config,
  char const * start, char const * end);

// components
size_t yac_couple_config_get_num_components(
  struct yac_couple_config * couple_config);
int yac_couple_config_component_name_is_valid(
  struct yac_couple_config * couple_config, char const * component_name);
size_t yac_couple_config_get_component_idx(
  struct yac_couple_config * couple_config, char const * component_name);
char const * yac_couple_config_get_component_name(
  struct yac_couple_config * couple_config, size_t component_idx);
void yac_couple_config_add_component(
  struct yac_couple_config * couple_config, char const * name);
void yac_couple_config_component_set_metadata(
  struct yac_couple_config * couple_config,
  char const * comp_name, const char* metadata);
const char * yac_couple_config_component_get_metadata(
  struct yac_couple_config * couple_config,
  char const * comp_name);
int yac_couple_config_component_name_is_valid(
  struct yac_couple_config * couple_config, char const * component_name);
size_t yac_couple_config_get_num_fields(
  struct yac_couple_config * couple_config, size_t component_idx);

// grids
size_t yac_couple_config_get_num_grids(
  struct yac_couple_config * couple_config);
int yac_couple_config_contains_grid_name(
  struct yac_couple_config * couple_config, char const * grid_name);
void yac_couple_config_add_grid(
  struct yac_couple_config * couple_config, char const * name);
void yac_couple_config_grid_set_metadata(
  struct yac_couple_config * couple_config,
  char const * grid_name, const char* metadata);
const char * yac_couple_config_grid_get_metadata(
  struct yac_couple_config * couple_config,
  char const * grid_name);
size_t yac_couple_config_get_grid_idx(
  struct yac_couple_config * couple_config, char const * grid_name);
char const * yac_couple_config_get_grid_name(
  struct yac_couple_config * couple_config, size_t grid_idx);
void yac_couple_config_grid_set_output_filename(
  struct yac_couple_config * couple_config, const char * grid_name,
  char const * output_filename);
const char* yac_couple_config_grid_get_output_filename(
  struct yac_couple_config * couple_config, const char * grid_name);

// fields
void yac_couple_config_component_add_field(
  struct yac_couple_config * couple_config,
    const char * component_name, const char * grid_name,
    const char * name, const char * timestep, size_t collection_size);
void yac_couple_config_field_set_metadata(
  struct yac_couple_config * couple_config,
  char const * comp_name, char const * grid_name, char const * field_name,
  const char* metadata);
const char * yac_couple_config_field_get_metadata(
  struct yac_couple_config * couple_config,
  char const * comp_name, char const * grid_name, char const * field_name);
size_t yac_couple_config_get_field_idx(
  struct yac_couple_config * couple_config, size_t component_idx,
    size_t grid_idx, char const * field_name);
void yac_couple_config_field_enable_frac_mask(
  struct yac_couple_config * couple_config,
  char const * comp_name, char const * grid_name, char const * field_name,
  double frac_mask_fallback_value);
double yac_couple_config_get_frac_mask_fallback_value(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name);
char const * yac_couple_config_get_field_grid_name(
  struct yac_couple_config * couple_config, size_t component_idx,
  size_t field_idx);
char const * yac_couple_config_get_field_name(
  struct yac_couple_config * couple_config, size_t component_idx,
  size_t field_idx);
char const * yac_couple_config_get_field_timestep(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name);
size_t yac_couple_config_get_field_collection_size(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name);
int yac_couple_config_get_field_role(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name);
int yac_couple_config_field_is_valid(
  struct yac_couple_config * couple_config,
  size_t component_idx, size_t field_idx);

// couples
void yac_couple_config_def_couple(
  struct yac_couple_config * couple_config,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_period, int time_reduction,
  struct yac_interp_stack_config * interp_stack,
  int src_lag, int tgt_lag,
  char const * weight_file_name, int mapping_on_source,
  double scale_factor, double scale_summand,
  size_t num_src_mask_names, char const * const * src_mask_names,
  char const * tgt_mask_name);
size_t yac_couple_config_get_num_couples(
  struct yac_couple_config * couple_config);
size_t yac_couple_config_get_num_couple_fields(
  struct yac_couple_config * couple_config, size_t couple_idx);
void yac_couple_config_get_couple_component_names(
  struct yac_couple_config * couple_config, size_t couple_idx,
  char const * couple_component_names[2]);
void yac_couple_config_get_field_couple_component_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const ** src_component_name, char const ** tgt_component_name);
void yac_couple_config_get_field_grid_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const ** src_grid_name, char const ** tgt_grid_name);
void yac_couple_config_get_field_names(
  struct yac_couple_config * couple_config,
    size_t couple_idx, size_t field_couple_idx,
    const char ** src_field_name, const char ** tgt_field_name);
struct yac_interp_stack_config * yac_couple_config_get_interp_stack(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
int yac_couple_config_mapping_on_source(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
int yac_couple_config_get_source_lag(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
int yac_couple_config_get_target_lag(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
char const * yac_couple_config_get_coupling_period(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
char const * yac_couple_config_get_source_timestep(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
char const * yac_couple_config_get_target_timestep(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
enum yac_reduction_type yac_couple_config_get_coupling_period_operation(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
int yac_couple_config_enforce_write_weight_file(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
char const * yac_couple_config_get_weight_file_name(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
double yac_couple_config_get_scale_factor(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
double yac_couple_config_get_scale_summand(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
void yac_couple_config_get_src_mask_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const * const ** mask_names, size_t * num_mask_names);
char const * yac_couple_config_get_tgt_mask_name(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);

/**
 * synchronises the coupling configuration across all processes in comm
 * @param[in] couple_config coupling configuration
 * @param[in] comm          MPI communicator
 * @param[in] output_ref    The coupling configuration will be written
 *                          to file after it has been synchronised between
 *                          all processes, if a filename and -type have been
 *                          set for the provided reference
 *                          (see \ref yac_couple_config_set_config_output_filename)
 */
void yac_couple_config_sync(
  struct yac_couple_config * couple_config, MPI_Comm comm,
  char const * output_ref);

// output configuration

/**
 * enables the writing of the synchronised coupling configuration to file by
 * \ref yac_couple_config_sync
 * @param[in] couple_config       coupling configuration
 * @param[in] filename            name of the output file
 * @param[in] filetype            type of the output file
 * @param[in] ref                 reference, which has to be provided to
 *                                \ref yac_couple_config_sync in order to select
 *                                the filename
 * @param[in] include_definitions include user definitions (components, grids,
 *                                and fields) in the output file
 */
void yac_couple_config_set_config_output_filename(
  struct yac_couple_config * couple_config,
  char const * filename, enum yac_text_filetype filetype, char const * ref,
  int include_definitions);

/**
 * returns whether YAC aborts if for a defined couple at least one
 * associated field was not defined by the user
 * @param[in] couple_config coupling configuration
 * @return `missing_definition_is_fatal` flag
 */
int yac_couple_config_get_missing_definition_is_fatal(
  struct yac_couple_config * couple_config);

/**
 * sets whether YAC aborts if for a defined couple at least one
 * associated field was not defined by the user
 * @param[in] couple_config               coupling configuration
 * @param[in] missing_definition_is_fatal `missing_definition_is_fatal` flag
 */
void yac_couple_config_set_missing_definition_is_fatal(
  struct yac_couple_config * couple_config,
  int missing_definition_is_fatal);

#endif // COUPLE_CONFIG_H
