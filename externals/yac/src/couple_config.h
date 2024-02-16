/**
 * @file couple_config.c
 * @brief Structs and interfaces to defined coupling configuration
 *
 * Stores reads in coupling configuration and sets up data structures
 * to do the coupling.
 *
 * @copyright Copyright  (C)  2021 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
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

/** \example test_couple_config.c
 * This example tests the generation of a coupling configuration.
 */

#ifndef COUPLE_CONFIG_H
#define COUPLE_CONFIG_H

#include "interp_stack_config.h"

#define MAX_ROUTINE_NAME_LENGTH (256)
#define MAX_FILE_NAME_LENGTH (512)

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

struct yac_couple_config;

// general couple config routines
struct yac_couple_config * yac_couple_config_new();
void yac_couple_config_delete(struct yac_couple_config * couple_config);

// basic run parameters
char const * yac_couple_config_get_start_datetime(
  struct yac_couple_config * couple_config);
char const * yac_couple_config_get_end_datetime(
  struct yac_couple_config * couple_config);
void yac_couple_config_set_datetime(
  struct yac_couple_config * couple_config,
  char const * start, char const * end);
int yac_couple_config_get_redirect_stdout(
  struct yac_couple_config * couple_config);
void yac_couple_config_set_redirect_stdout(
  struct yac_couple_config * couple_config, int redirect_stdout);

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

// grids
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
int yac_couple_config_get_field_role(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name);
int yac_couple_config_field_is_valid(
  struct yac_couple_config * couple_config,
  size_t component_idx, size_t field_idx);

// couples
size_t yac_couple_config_get_num_couples(
  struct yac_couple_config * couple_config);
size_t yac_couple_config_get_num_couple_fields(
  struct yac_couple_config * couple_config, size_t couple_idx);
void yac_couple_config_get_couple_component_names(
  struct yac_couple_config * couple_config, size_t couple_idx,
  char const * couple_component_names[2]);
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
int yac_couple_config_component_name_is_valid(
  struct yac_couple_config * couple_config, char const * component_name);
size_t yac_couple_config_get_num_components(
  struct yac_couple_config * couple_config);
size_t yac_couple_config_get_component_idx(
  struct yac_couple_config * couple_config, char const * component_name);
size_t yac_couple_config_get_num_grids(
  struct yac_couple_config * couple_config);
size_t yac_couple_config_get_num_fields(
  struct yac_couple_config * couple_config, size_t component_idx);
char const * yac_couple_config_get_component_name(
  struct yac_couple_config * couple_config, size_t component_idx);
struct yac_interp_stack_config * yac_couple_config_get_interp_stack(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx);
void yac_couple_config_get_field_grid_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const ** src_grid_name, char const ** tgt_grid_name);
void yac_couple_config_get_field_couple_component_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const ** src_component_name, char const ** tgt_component_name);
size_t yac_couple_config_get_collection_size(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name);
void yac_couple_config_get_field_names(
  struct yac_couple_config * couple_config,
    size_t couple_idx, size_t field_couple_idx,
    const char ** src_field_name, const char ** tgt_field_name);
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
char const * yac_couple_config_get_start_datetime(
  struct yac_couple_config * couple_config);
char const * yac_couple_config_get_end_datetime(
  struct yac_couple_config * couple_config);
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
 */
void yac_couple_config_sync(
  struct yac_couple_config * couple_config, MPI_Comm comm);

#endif // COUPLE_CONFIG_H
