/**
 * @file instance.c
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "utils.h"
#include "yac_interface.h"
#include "event.h"
#include "yac_mpi.h"
#include "fields.h"
#include "component.h"
#include "interpolation.h"
#include "ensure_array_size.h"
#include "geometry.h"
#include "dist_grid.h"
#include "interp_grid.h"
#include "interp_weights.h"
#include "interp_method.h"
#include "config_yaml.h"

enum yac_instance_phase {
  INSTANCE_DEFINITION = 0, // after yac_cinit
  INSTANCE_DEFINITION_COMP = 1, // after yac_cdef_comp
  INSTANCE_DEFINITION_SYNC = 2, // after yac_csync_def
  INSTANCE_EXCHANGE = 3, // after yac_cenddef
  INSTANCE_UNKNOWN = 4,
};

static const char * yac_instance_phase_str[] =
  {"definition phase",
   "definition phase (after component definition)",
   "definition phase (after synchronisation)",
   "exchange phase",
   "unknown phase"};

#define CHECK_PHASE(FUNC_NAME, REF_PHASE, NEW_PHASE) \
  { \
    enum yac_instance_phase ref_phase_ = (REF_PHASE); \
    YAC_ASSERT_F( \
      instance->phase == (ref_phase_), \
      "ERROR(%s): Invalid phase " \
      "(current phase: \"%s\" expected phase: \"%s\")", \
      #FUNC_NAME, yac_instance_phase_str[instance->phase], \
      yac_instance_phase_str[(ref_phase_)]); \
    instance->phase = (NEW_PHASE); \
  }
#define CHECK_MIN_PHASE(FUNC_NAME, MIN_REF_PHASE) \
  { \
    enum yac_instance_phase ref_min_phase_ = (MIN_REF_PHASE); \
    YAC_ASSERT_F( \
      instance->phase >= (ref_min_phase_), \
      "ERROR(%s): Invalid phase " \
      "(current phase: \"%s\" minimum expected phase: \"%s\")", \
      #FUNC_NAME, yac_instance_phase_str[instance->phase], \
      yac_instance_phase_str[(ref_min_phase_)]); \
  }
#define CHECK_MAX_PHASE(FUNC_NAME, MAX_REF_PHASE) \
  { \
    enum yac_instance_phase ref_max_phase_ = (MAX_REF_PHASE); \
    YAC_ASSERT_F( \
      instance->phase <= (ref_max_phase_), \
      "ERROR(%s): Invalid phase " \
      "(current phase: \"%s\" maximum expected phase: \"%s\")", \
      #FUNC_NAME, \
      yac_instance_phase_str[ \
        MIN(instance->phase,INSTANCE_UNKNOWN)], \
      yac_instance_phase_str[(ref_max_phase_)]); \
  }

struct yac_instance {

  struct yac_couple_config * couple_config;

  struct yac_component_config * comp_config;

  struct coupling_field ** cpl_fields;
  size_t num_cpl_fields;

  MPI_Comm comm;

  enum yac_instance_phase phase;
};

enum field_type {
  SRC = 1,
  TGT = 2
};

struct comp_grid_config {
  const char * grid_name;
  const char * comp_name;
};

struct comp_grid_pair_config {
  struct comp_grid_config config[2];
};

static struct field_config_event_data {
  char const * timestep;
  char const * coupling_period;
  int timelag;
  enum yac_reduction_type reduction_operation;
  char const * start_datetime;
  char const * end_datetime;
} empty_event_data =
  {.timestep = NULL,
   .coupling_period = NULL,
   .timelag = 0,
   .reduction_operation = TIME_NONE,
   .start_datetime = NULL,
   .end_datetime = NULL};

struct src_field_config {
  struct coupling_field * field;
  char const * name;
  int id;

  char const * const * mask_names;
  size_t num_mask_names;

  struct yac_interp_stack_config * interp_stack;
  const char * weight_file_name;
  struct field_config_event_data event_data;
};

struct tgt_field_config {
  struct coupling_field * field;
  char const * name;
  int id;
  char const * mask_name;
  struct field_config_event_data event_data;
};

struct field_config {

  char const * name;
  int id;

  struct comp_grid_pair_config comp_grid_pair;
  int src_comp_idx;

  struct src_field_config src_interp_config;
  struct tgt_field_config tgt_interp_config;

  enum interp_weights_reorder_type reorder_type;

  double frac_mask_fallback_value;

  double scale_factor, scale_summand;

  size_t collection_size;
};

static struct yac_basic_grid *
  get_basic_grid(const char * grid_name, char const * component_name,
                 size_t num_fields, struct coupling_field ** coupling_fields,
                 int * delete_flag) {

  *delete_flag = 0;

  for (size_t i = 0; i < num_fields; ++i) {

    struct coupling_field * field = coupling_fields[i];

    if (strcmp(component_name, yac_get_coupling_field_comp_name(field)))
      continue;

    struct yac_basic_grid * grid = yac_coupling_field_get_basic_grid(field);

    if (strcmp(grid_name, yac_basic_grid_get_name(grid)) == 0)
      return grid;
  }

  *delete_flag = 1;
  return yac_basic_grid_empty_new(grid_name);
}

static int contains_component(
  struct component ** components, size_t num_components,
  char const * component_name) {

  for (size_t i = 0; i < num_components; ++i)
    if (!strcmp(component_name, yac_component_get_name(components[i])))
      return 1;
  return 0;
}

static
int compare_comp_grid_config(const void * a, const void * b) {

  struct comp_grid_config * a_ = (struct comp_grid_config *)a,
                          * b_ = (struct comp_grid_config *)b;

  int ret;
  if ((ret = strcmp(a_->comp_name, b_->comp_name))) return ret;
  else return strcmp(a_->grid_name, b_->grid_name);
}

static struct coupling_field * get_coupling_field(
  char const * component_name, const char * field_name,
  const char * grid_name, size_t num_fields,
  struct coupling_field ** coupling_fields) {

  for (size_t i = 0; i < num_fields; ++i) {
    struct coupling_field * curr_field = coupling_fields[i];
    if (!strcmp(component_name, yac_get_coupling_field_comp_name(curr_field)) &&
        !strcmp(field_name, yac_get_coupling_field_name(curr_field)) &&
        !strcmp(grid_name,
                yac_basic_grid_get_name(
                  yac_coupling_field_get_basic_grid(curr_field))))
      return curr_field;
  }

  return NULL;
}

static struct src_field_config get_src_interp_config(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  char const * const * src_mask_names;
  size_t num_src_mask_names;
  yac_couple_config_get_src_mask_names(
    couple_config, couple_idx, field_couple_idx,
    &src_mask_names, &num_src_mask_names);

  return
    (struct src_field_config) {
    .field = NULL,
    .name = NULL,
    .id = INT_MAX,
    .mask_names = src_mask_names,
    .num_mask_names = num_src_mask_names,
    .interp_stack =
      yac_couple_config_get_interp_stack(
        couple_config, couple_idx, field_couple_idx),
    .weight_file_name =
      (yac_couple_config_enforce_write_weight_file(
         couple_config, couple_idx, field_couple_idx))?
        yac_couple_config_get_weight_file_name(
          couple_config, couple_idx, field_couple_idx):NULL,
    .event_data = empty_event_data,
  };
}

static struct tgt_field_config get_tgt_interp_config(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  return
    (struct tgt_field_config) {
    .field = NULL,
    .name = NULL,
    .id = INT_MAX,
    .mask_name =
      yac_couple_config_get_tgt_mask_name(
        couple_config, couple_idx, field_couple_idx),
    .event_data = empty_event_data,
  };
}

static int compare_field_config_fields(
  struct coupling_field * a,
  size_t num_mask_names_a, char const * const * mask_names_a,
  struct coupling_field * b,
  size_t num_mask_names_b, char const * const * mask_names_b) {

  if ((a == NULL) || (b == NULL)) return (a == NULL) - (b == NULL);

  size_t num_interp_fields[2] =
    {yac_coupling_field_get_num_interp_fields(a),
     yac_coupling_field_get_num_interp_fields(b)};

  YAC_ASSERT(
    (num_mask_names_a == 0) ||
    (num_mask_names_a == num_interp_fields[0]),
    "ERROR(compare_field_config_fields): inconsistent mask names defined")
  YAC_ASSERT(
    (num_mask_names_b == 0) ||
    (num_mask_names_b == num_interp_fields[1]),
    "ERROR(compare_field_config_fields): inconsistent mask names defined")

  // if both fields have different numbers of interpolation fields
  if (num_interp_fields[0] != num_interp_fields[1])
    return (int)num_interp_fields[0] - (int)num_interp_fields[1];

  struct yac_basic_grid * grid[2] =
    {yac_coupling_field_get_basic_grid(a),
     yac_coupling_field_get_basic_grid(b)};

  struct interp_field const * interp_fields[2] =
    {yac_coupling_field_get_interp_fields(a),
     yac_coupling_field_get_interp_fields(b)};

  for (size_t i = 0; i < num_interp_fields[0]; ++i) {

    int ret;

    enum yac_location location_a = interp_fields[0][i].location;
    enum yac_location location_b = interp_fields[1][i].location;

    if ((ret = (location_a > location_b) -
               (location_a < location_b))) return ret;

    coordinate_pointer coordinates_a =
      yac_basic_grid_get_field_coordinates(grid[0], interp_fields[0][i]);
    coordinate_pointer coordinates_b =
      yac_basic_grid_get_field_coordinates(grid[1], interp_fields[1][i]);

    if ((ret = (coordinates_a > coordinates_b) -
               (coordinates_a < coordinates_b))) return ret;

    char const * const * mask_names[2] =
      {num_mask_names_a?mask_names_a:NULL,
       num_mask_names_b?mask_names_b:NULL};
    int const * masks[2];

    for (int j = 0; j < 2; ++j) {
      struct interp_field temp_interp_field = interp_fields[j][i];
      if (mask_names[j] != NULL)
        temp_interp_field.masks_idx =
          yac_basic_grid_get_named_mask_idx(
            grid[j], temp_interp_field.location, mask_names[j][i]);

      masks[j] = yac_basic_grid_get_field_mask(grid[j], temp_interp_field);
    }
    if ((ret = (masks[0] > masks[1]) - (masks[0] < masks[1]))) return ret;
  }

  return 0;
}

static int compare_field_config_field_ids(const void * a, const void * b) {

  struct field_config * a_ = (struct field_config *)a,
                      * b_ = (struct field_config *)b;

  int ret;
  if ((ret = a_->src_interp_config.id - b_->src_interp_config.id)) return ret;
  return a_->tgt_interp_config.id - b_->tgt_interp_config.id;

}

static int
compare_field_config_interp_method(const void * a, const void * b) {

  struct field_config * a_ = (struct field_config *)a,
                      * b_ = (struct field_config *)b;

  int ret;
  if ((ret = (int)(a_->src_interp_config.weight_file_name == NULL) -
             (int)(b_->src_interp_config.weight_file_name == NULL))) return ret;
  if ((a_->src_interp_config.weight_file_name != NULL) &&
      (ret = strcmp(a_->src_interp_config.weight_file_name,
                    b_->src_interp_config.weight_file_name))) return ret;

  return
    yac_interp_stack_config_compare(
      a_->src_interp_config.interp_stack,
      b_->src_interp_config.interp_stack);
}

static
int compare_field_config(const void * a, const void * b) {

  struct field_config * a_ = (struct field_config *)a,
                      * b_ = (struct field_config *)b;

  int ret;
  if ((ret = strcmp(a_->comp_grid_pair.config[0].comp_name,
                    b_->comp_grid_pair.config[0].comp_name))) return ret;
  if ((ret = strcmp(a_->comp_grid_pair.config[1].comp_name,
                    b_->comp_grid_pair.config[1].comp_name))) return ret;
  if ((ret = strcmp(a_->comp_grid_pair.config[0].grid_name,
                    b_->comp_grid_pair.config[0].grid_name))) return ret;
  if ((ret = strcmp(a_->comp_grid_pair.config[1].grid_name,
                    b_->comp_grid_pair.config[1].grid_name))) return ret;
  if ((ret = compare_field_config_field_ids(a_, b_))) return ret;
  if ((ret = compare_field_config_interp_method(a_, b_))) return ret;
  if ((ret = (int)(a_->reorder_type) - (int)(b_->reorder_type))) return ret;
  if ((ret = (a_->scale_factor > b_->scale_factor) -
             (a_->scale_factor < b_->scale_factor))) return ret;
  if ((ret = (a_->scale_summand > b_->scale_summand) -
             (a_->scale_summand < b_->scale_summand))) return ret;

  // use memcmp to compare frac_mask_fallback_value, because they can be nan
  return
    memcmp(
      &(a_->frac_mask_fallback_value),
      &(b_->frac_mask_fallback_value), sizeof(double));
}

struct field_config_event_data get_event_data(
  struct yac_instance * instance, int couple_idx, int field_couple_idx,
  enum field_type field_type) {

  struct yac_couple_config * couple_config = instance->couple_config;

  enum yac_reduction_type reduction_operation =
    yac_couple_config_get_coupling_period_operation(
      couple_config, couple_idx, field_couple_idx);

  char const * coupling_period =
    yac_couple_config_get_coupling_period(
      couple_config, couple_idx, field_couple_idx);
  char const * timestep;
  int timelag;
  if ( field_type == SRC ){
    timestep =
      yac_couple_config_get_source_timestep(
        couple_config, couple_idx, field_couple_idx);
    timelag =
      yac_couple_config_get_source_lag(
        couple_config, couple_idx, field_couple_idx);
  }else{
    reduction_operation = TIME_NONE;
    timestep =
      yac_couple_config_get_target_timestep(
        couple_config, couple_idx, field_couple_idx);
    timelag =
      yac_couple_config_get_target_lag(
        couple_config, couple_idx, field_couple_idx);
  }

  return
    (struct field_config_event_data)
      {.timestep = timestep,
       .coupling_period = coupling_period,
       .timelag = timelag,
       .reduction_operation = reduction_operation,
       .start_datetime = yac_couple_config_get_start_datetime(couple_config),
       .end_datetime = yac_couple_config_get_end_datetime(couple_config)};
}

static struct event * generate_event(
  struct field_config_event_data event_data) {

  struct event * event = yac_event_new();
  yac_event_add(
    event, event_data.timestep, event_data.coupling_period,
    event_data.timelag, event_data.reduction_operation,
    event_data.start_datetime, event_data.end_datetime);
  return event;
}

static void get_field_configuration(
  struct yac_instance * instance,
  struct field_config ** field_configs_, size_t * count) {

  struct yac_couple_config * couple_config = instance->couple_config;
  MPI_Comm comm = instance->comm;
  size_t num_fields = instance->num_cpl_fields;
  struct coupling_field ** coupling_fields = instance->cpl_fields;

  size_t num_couples = yac_couple_config_get_num_couples(couple_config);
  size_t total_num_fields = 0;
  for (size_t couple_idx = 0; couple_idx < num_couples; ++couple_idx)
    total_num_fields +=
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

  // determine for which configuration source and target field
  // are available
  int (*fields_available)[2] =
    xmalloc(total_num_fields * sizeof(*fields_available));
  for (size_t couple_idx = 0, i = 0; couple_idx < num_couples;
       ++couple_idx) {
    size_t curr_num_fields =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);
    for (size_t field_couple_idx = 0; field_couple_idx < curr_num_fields;
         ++field_couple_idx, ++i) {
      char const * src_component_name;
      char const * tgt_component_name;
      char const * src_grid_name;
      char const * tgt_grid_name;
      const char * src_field_name;
      const char * tgt_field_name;
      yac_couple_config_get_field_couple_component_names(
        couple_config, couple_idx, field_couple_idx,
        &src_component_name, &tgt_component_name);
      yac_couple_config_get_field_grid_names(
        couple_config, couple_idx, field_couple_idx,
        &src_grid_name, &tgt_grid_name);
      yac_couple_config_get_field_names(
        couple_config, couple_idx, field_couple_idx,
          &src_field_name, &tgt_field_name);
      fields_available[i][0] =
        get_coupling_field(
          src_component_name, src_field_name, src_grid_name,
          num_fields, coupling_fields) != NULL;
      fields_available[i][1] =
        get_coupling_field(
          tgt_component_name, tgt_field_name, tgt_grid_name,
          num_fields, coupling_fields) != NULL;
    }
  }
  yac_mpi_call(
    MPI_Allreduce(MPI_IN_PLACE, fields_available,
                  2 * total_num_fields,
                  MPI_INT, MPI_MAX, comm), comm);
  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  size_t new_total_num_fields = 0;
  for (size_t couple_idx = 0, i = 0; couple_idx < num_couples;
       ++couple_idx) {
    size_t curr_num_fields =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);
    for (size_t field_couple_idx = 0; field_couple_idx < curr_num_fields;
         ++field_couple_idx, ++i) {
      if (fields_available[i][0] && fields_available[i][1]) {
        ++new_total_num_fields;
      } else if (comm_rank == 0) {
        char const * src_component_name;
        char const * tgt_component_name;
        char const * src_grid_name;
        char const * tgt_grid_name;
        const char * src_field_name;
        const char * tgt_field_name;
        yac_couple_config_get_field_couple_component_names(
          couple_config, couple_idx, field_couple_idx,
          &src_component_name, &tgt_component_name);
        yac_couple_config_get_field_grid_names(
          couple_config, couple_idx, field_couple_idx,
          &src_grid_name, &tgt_grid_name);
        yac_couple_config_get_field_names(
          couple_config, couple_idx, field_couple_idx,
            &src_field_name, &tgt_field_name);
        fprintf(stderr, "WARNING: couple defined for field: \n"
                        "  source (%s):\n"
                        "    component name: %s\n"
                        "    grid name:      %s\n"
                        "    field name:     %s\n"
                        "  target(%s):\n"
                        "    component name: %s\n"
                        "    grid name:      %s\n"
                        "    field name:     %s\n",
                (fields_available[i][0])?"defined":"not defined",
                src_component_name, src_grid_name, src_field_name,
                (fields_available[i][1])?"defined":"not defined",
                tgt_component_name, tgt_grid_name, tgt_field_name);
      }
    }
  }
  total_num_fields = new_total_num_fields;

  struct field_config * field_configs =
    xmalloc(total_num_fields * sizeof(*field_configs));

  size_t field_config_idx = 0;

  // extract component, grid configurations from all fields in
  // the coupling configuration
  for (size_t couple_idx = 0, i = 0; couple_idx < num_couples; ++couple_idx) {

    size_t curr_num_fields =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

    for (size_t field_couple_idx = 0; field_couple_idx < curr_num_fields;
         ++field_couple_idx, ++i) {

      if (!fields_available[i][0] || !fields_available[i][1]) continue;

      struct comp_grid_config src_config, tgt_config;

      yac_couple_config_get_field_couple_component_names(
        couple_config, couple_idx, field_couple_idx,
        &(src_config.comp_name), &(tgt_config.comp_name));

      yac_couple_config_get_field_grid_names(
        couple_config, couple_idx, field_couple_idx,
        &(src_config.grid_name), &(tgt_config.grid_name));

      int src_comp_idx =
        compare_comp_grid_config(&src_config, &tgt_config) > 0;

      const char * src_field_name;
      const char * tgt_field_name;
      yac_couple_config_get_field_names(
        couple_config, couple_idx, field_couple_idx,
          &src_field_name, &tgt_field_name);
      struct coupling_field * src_field =
        get_coupling_field(
          src_config.comp_name, src_field_name, src_config.grid_name,
          num_fields, coupling_fields);
      struct coupling_field * tgt_field =
        get_coupling_field(
          tgt_config.comp_name, tgt_field_name, tgt_config.grid_name,
          num_fields, coupling_fields);

      double frac_mask_fallback_value =
        yac_couple_config_get_frac_mask_fallback_value(
          couple_config, src_config.comp_name, src_config.grid_name,
          src_field_name);
      double scale_factor =
        yac_couple_config_get_scale_factor(
          couple_config, couple_idx, field_couple_idx);
      double scale_summand =
        yac_couple_config_get_scale_summand(
          couple_config, couple_idx, field_couple_idx);
      size_t collection_size =
        yac_couple_config_get_collection_size(
          couple_config, src_config.comp_name, src_config.grid_name,
          src_field_name);

      field_configs[field_config_idx].comp_grid_pair.config[src_comp_idx] =
        src_config;
      field_configs[field_config_idx].comp_grid_pair.config[src_comp_idx^1] =
        tgt_config;
      field_configs[field_config_idx].src_comp_idx = src_comp_idx;
      field_configs[field_config_idx].src_interp_config =
        get_src_interp_config(couple_config, couple_idx, field_couple_idx);
      field_configs[field_config_idx].tgt_interp_config =
        get_tgt_interp_config(couple_config, couple_idx, field_couple_idx);
      field_configs[field_config_idx].src_interp_config.field = src_field;
      field_configs[field_config_idx].tgt_interp_config.field = tgt_field;
      field_configs[field_config_idx].reorder_type =
        (yac_couple_config_mapping_on_source(
           couple_config, couple_idx, field_couple_idx))?
        (MAPPING_ON_SRC):(MAPPING_ON_TGT);
      field_configs[field_config_idx].src_interp_config.event_data =
        get_event_data(instance, couple_idx, field_couple_idx, SRC);
      field_configs[field_config_idx].tgt_interp_config.event_data =
        get_event_data(instance, couple_idx, field_couple_idx, TGT);
      field_configs[field_config_idx].frac_mask_fallback_value =
        frac_mask_fallback_value;
      field_configs[field_config_idx].scale_factor = scale_factor;
      field_configs[field_config_idx].scale_summand = scale_summand;
      field_configs[field_config_idx].collection_size = collection_size;
      field_configs[field_config_idx].src_interp_config.name = src_field_name;
      field_configs[field_config_idx].tgt_interp_config.name = tgt_field_name;
      ++field_config_idx;
    }
  }

  free(fields_available);

  { // determine unique field config ids, for each unique configuration
    int * field_src_field_config_compare_matrix =
      xmalloc((size_t)(total_num_fields * (total_num_fields - 1)) *
              sizeof(*field_src_field_config_compare_matrix));
    int * field_tgt_field_config_compare_matrix =
      field_src_field_config_compare_matrix +
      (size_t)(total_num_fields * (total_num_fields - 1)) / 2;

    for (unsigned i = 0, k = 0; i < total_num_fields; ++i) {
      for (unsigned j = i + 1; j < total_num_fields; ++j, ++k) {
        field_src_field_config_compare_matrix[k] =
          compare_field_config_fields(
            field_configs[i].src_interp_config.field,
            field_configs[i].src_interp_config.num_mask_names,
            field_configs[i].src_interp_config.mask_names,
            field_configs[j].src_interp_config.field,
            field_configs[j].src_interp_config.num_mask_names,
            field_configs[j].src_interp_config.mask_names) != 0;
        field_tgt_field_config_compare_matrix[k] =
          compare_field_config_fields(
            field_configs[i].tgt_interp_config.field,
            field_configs[i].tgt_interp_config.mask_name != NULL,
            &(field_configs[i].tgt_interp_config.mask_name),
            field_configs[j].tgt_interp_config.field,
            field_configs[j].tgt_interp_config.mask_name != NULL,
            &(field_configs[j].tgt_interp_config.mask_name)) != 0;
      }
    }

    yac_mpi_call(
      MPI_Allreduce(MPI_IN_PLACE, field_src_field_config_compare_matrix,
                    total_num_fields * (total_num_fields - 1),
                    MPI_INT, MPI_LOR, comm), comm);

    int * field_src_field_config_ids =
      xcalloc(2 * (size_t)total_num_fields,
              sizeof(*field_src_field_config_ids));
    int * field_tgt_field_config_ids =
      field_src_field_config_ids + (size_t)total_num_fields;

    for (unsigned i = 0, k = 0, id = 1; i < total_num_fields; ++i) {
      if (field_src_field_config_ids[i] != 0) {
        k += total_num_fields - i - 1;
        continue;
      } else {
        field_src_field_config_ids[i] = (int)id;
        for (unsigned j = i + 1; j < total_num_fields; ++j, ++k)
          if (!field_src_field_config_compare_matrix[k])
            field_src_field_config_ids[j] = (int)id;
        id++;
      }
    }
    for (unsigned i = 0, k = 0, id = 1; i < total_num_fields; ++i) {
      if (field_tgt_field_config_ids[i] != 0) {
        k += total_num_fields - i - 1;
        continue;
      } else {
        field_tgt_field_config_ids[i] = (int)id;
        for (unsigned j = i + 1; j < total_num_fields; ++j, ++k)
          if (!field_tgt_field_config_compare_matrix[k])
            field_tgt_field_config_ids[j] = (int)id;
        id++;
      }
    }
    free(field_src_field_config_compare_matrix);

    for (unsigned i = 0; i < total_num_fields; ++i) {
      field_configs[i].src_interp_config.id =
        field_src_field_config_ids[i];
      field_configs[i].tgt_interp_config.id =
        field_tgt_field_config_ids[i];
    }
    free(field_src_field_config_ids);
  }

  { // determine unique field config ids, for each unique configuration
    int * field_config_compare_matrix =
      xmalloc((size_t)(total_num_fields * (total_num_fields - 1)) / 2 *
              sizeof(*field_config_compare_matrix));
    for (unsigned i = 0, k = 0; i < total_num_fields; ++i)
      for (unsigned j = i + 1; j < total_num_fields; ++j, ++k)
        field_config_compare_matrix[k] =
          compare_field_config(
            &(field_configs[i]), &(field_configs[j])) != 0;

    yac_mpi_call(
      MPI_Allreduce(MPI_IN_PLACE, field_config_compare_matrix,
                    (total_num_fields * (total_num_fields - 1)) / 2,
                    MPI_INT, MPI_LOR, comm), comm);

    int * field_config_ids =
      xcalloc((size_t)total_num_fields, sizeof(*field_config_ids));

    for (unsigned i = 0, k = 0, id = 1; i < total_num_fields; ++i) {
      if (field_config_ids[i] != 0) {
        k += total_num_fields - i - 1;
        continue;
      } else {
        field_config_ids[i] = (int)id;
        for (unsigned j = i + 1; j < total_num_fields; ++j, ++k)
          if (!field_config_compare_matrix[k])
            field_config_ids[j] = (int)id;
        id++;
      }
    }
    free(field_config_compare_matrix);

    for (unsigned i = 0; i < total_num_fields; ++i)
      field_configs[i].id = field_config_ids[i];
    free(field_config_ids);
  }

  *field_configs_ = field_configs;
  *count = total_num_fields;
}

static
int compare_field_config_ids(const void * a, const void * b) {

  struct field_config * a_ = (struct field_config *)a,
                      * b_ = (struct field_config *)b;

  return a_->id - b_->id;
}

static struct interp_weights * generate_interp_weights(
  struct src_field_config src_interp_config,
  struct interp_grid * interp_grid) {

  struct interp_method ** method_stack =
    yac_interp_stack_config_generate(src_interp_config.interp_stack);
  struct interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);
  free(method_stack);

  return weights;
}

static void get_interp_fields_from_coupling_field(
  struct coupling_field * field, char const * const * mask_names,
  size_t num_mask_names, struct interp_field ** interp_fields,
  size_t * num_fields, MPI_Comm comm) {

  struct interp_field * interp_fields_;
  size_t num_fields_;

  if (field != NULL) {

    num_fields_ = yac_coupling_field_get_num_interp_fields(field);

    YAC_ASSERT_F(
      (num_mask_names == 0) || (num_fields_ == num_mask_names),
      "ERROR(get_interp_fields_from_coupling_field): "
      "missmatch in number of interpolation fields of coupling field \"%s\" "
      "and number of provided mask names (%zu != %zu)",
      yac_get_coupling_field_name(field), num_fields_, num_mask_names)

    uint64_t local_counts[2] =
      {(uint64_t)num_fields_, (uint64_t)num_mask_names};
    uint64_t global_counts[2];

    yac_mpi_call(
      MPI_Allreduce(
        local_counts, global_counts, 2, MPI_UINT64_T, MPI_MAX, comm), comm);

    YAC_ASSERT_F(
      local_counts[0] == global_counts[0],
      "ERROR(get_interp_fields_from_coupling_field): missmatch in number of"
      "local interpolation fields for coupling field \"%s\" and global "
      "(%zu != %zu)",
      yac_get_coupling_field_name(field), num_fields_, (size_t)(global_counts[0]))

    YAC_ASSERT_F(
      (num_mask_names != 0) || (global_counts[1] == 0),
      "ERROR(get_interp_fields_from_coupling_field): local process did "
      "not provide mask names for coupling field \"%s\" while others did",
      yac_get_coupling_field_name(field));

    // make a copy of the interpolation fields of the coupling field
    interp_fields_ = xmalloc(num_fields_ * sizeof(*interp_fields_));
    memcpy(
      interp_fields_, yac_coupling_field_get_interp_fields(field),
      num_fields_ * sizeof(*interp_fields_));

    // if mask names are provided, overwrite already existing masks
    struct yac_basic_grid * grid = yac_coupling_field_get_basic_grid(field);
    for (size_t i = 0; i < num_mask_names; ++i) {
      char const * mask_name = mask_names[i];
      YAC_ASSERT_F(
        mask_name != NULL,
        "ERROR(get_interp_fields_from_coupling_field): "
        "make_names[%zu] is NULL", i);
      interp_fields_[i].masks_idx =
        yac_basic_grid_get_named_mask_idx(
          grid, interp_fields_[i].location, mask_name);
    }

    uint64_t data[num_fields_][3];

    for (size_t i = 0; i < num_fields_; ++i) {
      data[i][0] = (uint64_t)interp_fields_[i].location;
      data[i][1] = (uint64_t)interp_fields_[i].coordinates_idx;
      data[i][2] = (uint64_t)interp_fields_[i].masks_idx;
    }

    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, data, 3 * (int)num_fields_,
        MPI_UINT64_T, MPI_MIN, comm), comm);

    for (size_t i = 0; i < num_fields_; ++i) {
      YAC_ASSERT(
        data[i][0] == (uint64_t)(interp_fields_[i].location),
        "ERROR(get_interp_fields_from_coupling_field): location mismatch")
      YAC_ASSERT(
        data[i][1] == (uint64_t)(interp_fields_[i].coordinates_idx),
        "ERROR(get_interp_fields_from_coupling_field): "
        "coordinates index mismatch")
      YAC_ASSERT(
        data[i][2] == (uint64_t)(interp_fields_[i].masks_idx),
        "ERROR(get_interp_fields_from_coupling_field): "
        "masks index mismatch")
    }

  } else {

    uint64_t zero_counts[2] = {0,0};
    uint64_t counts[2];

    yac_mpi_call(
      MPI_Allreduce(
        zero_counts, counts, 2,
        MPI_UINT64_T, MPI_MAX, comm), comm);

    num_fields_ = (size_t)(counts[0]);
    interp_fields_ = xmalloc(num_fields_ * sizeof(*interp_fields_));

    uint64_t data[num_fields_][3];

    for (size_t i = 0; i < num_fields_; ++i) {
      data[i][0] = (uint64_t)LOC_UNDEFINED;
      data[i][1] = (uint64_t)UINT64_MAX;
      data[i][2] = (uint64_t)UINT64_MAX;
    }

    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, data, 3 * (int)num_fields_,
        MPI_UINT64_T, MPI_MIN, comm), comm);

    for (size_t i = 0; i < num_fields_; ++i) {
      interp_fields_[i].location = (enum yac_location)data[i][0];
      interp_fields_[i].coordinates_idx = (size_t)data[i][1];
      interp_fields_[i].masks_idx = (size_t)data[i][2];
    }
  }

  *interp_fields = interp_fields_;
  *num_fields = num_fields_;
}

static void generate_interpolations(struct yac_instance * instance) {

  MPI_Comm comm = instance->comm;
  size_t num_components =
    yac_component_config_get_num_components(instance->comp_config);
  struct component ** components =
    yac_component_config_get_components(instance->comp_config);

  // get information about all fields
  struct field_config * field_configs;
  size_t field_count;
  get_field_configuration(
    instance, &field_configs, &field_count);

  // sort field configurations
  qsort(field_configs, field_count, sizeof(*field_configs),
        compare_field_config_ids);

  struct dist_grid_pair * dist_grid_pair = NULL;
  struct interp_grid * interp_grid = NULL;
  struct interp_weights * interp_weights = NULL;
  struct interpolation * interp = NULL;
  struct comp_grid_pair_config * prev_comp_grid_pair = NULL;
  MPI_Comm comp_pair_comm = MPI_COMM_NULL;
  int is_active = 0;

  // loop over all fields to build interpolations
  for (size_t i = 0; i < field_count; ++i) {

    struct field_config * curr_field_config = field_configs + i;
    struct comp_grid_pair_config * curr_comp_grid_pair =
      &(curr_field_config->comp_grid_pair);

    int build_flag = 0;
    // if the current comp pair is different from the previous one
    // => generate a component pair communicator
    if ((prev_comp_grid_pair == NULL) ||
        (strcmp(prev_comp_grid_pair->config[0].comp_name,
                curr_comp_grid_pair->config[0].comp_name) ||
         strcmp(prev_comp_grid_pair->config[1].comp_name,
                curr_comp_grid_pair->config[1].comp_name))) {

      build_flag = 1;

      if (comp_pair_comm != MPI_COMM_NULL)
        yac_mpi_call(MPI_Comm_free(&comp_pair_comm), comm);

      int comp_is_available[2] =
        {contains_component(
           components, num_components,
           curr_comp_grid_pair->config[0].comp_name),
         contains_component(
           components, num_components,
           curr_comp_grid_pair->config[1].comp_name)};

      is_active = comp_is_available[0] || comp_is_available[1];

      yac_mpi_call(MPI_Comm_split(comm, is_active, 0,
                                  &comp_pair_comm), comm);

      // determine whether both components are available in the setup
      yac_mpi_call(MPI_Allreduce(MPI_IN_PLACE, comp_is_available, 2,
                                 MPI_INT, MPI_LOR, comp_pair_comm), comm);
      is_active = comp_is_available[0] && comp_is_available[1];
    }

    // if the current configuration differs from the previous one and the local
    // process is involved in this configuration
    if (is_active &&
        (build_flag ||
         strcmp(prev_comp_grid_pair->config[0].grid_name,
                curr_comp_grid_pair->config[0].grid_name) ||
         strcmp(prev_comp_grid_pair->config[1].grid_name,
                curr_comp_grid_pair->config[1].grid_name))) {

      build_flag = 1;

      if (dist_grid_pair != NULL) yac_dist_grid_pair_delete(dist_grid_pair);

      char const * grid_names[2] =
        {curr_comp_grid_pair->config[0].grid_name,
         curr_comp_grid_pair->config[1].grid_name};
      char const * component_names[2] =
        {curr_comp_grid_pair->config[0].comp_name,
         curr_comp_grid_pair->config[1].comp_name};

      int delete_flags[2];
      struct yac_basic_grid * basic_grid[2] =
        {get_basic_grid(
           grid_names[0], component_names[0], instance->num_cpl_fields,
           instance->cpl_fields, &delete_flags[0]),
         get_basic_grid(
           grid_names[1], component_names[1], instance->num_cpl_fields,
           instance->cpl_fields, &delete_flags[1])};

      dist_grid_pair =
        yac_dist_grid_pair_new(basic_grid[0], basic_grid[1], comp_pair_comm);

      for (int i = 0; i < 2; ++i)
        if (delete_flags[i]) yac_basic_grid_delete(basic_grid[i]);
    }

    // if the current source or target field data differes from the previous
    // one
    if (is_active &&
        (build_flag ||
         compare_field_config_field_ids(
           curr_field_config - 1, curr_field_config))) {

      build_flag = 1;

      struct interp_field * src_fields;
      size_t num_src_fields;
      get_interp_fields_from_coupling_field(
          curr_field_config->src_interp_config.field,
          curr_field_config->src_interp_config.mask_names,
          curr_field_config->src_interp_config.num_mask_names,
          &src_fields, &num_src_fields, comp_pair_comm);
      struct interp_field * tgt_fields;
      size_t num_tgt_fields;
      get_interp_fields_from_coupling_field(
          curr_field_config->tgt_interp_config.field,
          &(curr_field_config->tgt_interp_config.mask_name),
          curr_field_config->tgt_interp_config.mask_name != NULL,
          &tgt_fields, &num_tgt_fields, comp_pair_comm);

      YAC_ASSERT(
        num_tgt_fields == 1,
        "ERROR(generate_interpolations): "
        "only one point set per target field supported")

      if (interp_grid != NULL) yac_interp_grid_delete(interp_grid);

      int src_comp_idx = field_configs[i].src_comp_idx;
      interp_grid = yac_interp_grid_new(
        dist_grid_pair,
        curr_comp_grid_pair->config[src_comp_idx].grid_name,
        curr_comp_grid_pair->config[src_comp_idx^1].grid_name,
        num_src_fields, src_fields, *tgt_fields);

      free(tgt_fields);
      free(src_fields);
    }

    // if the current interpolation method stack differes from the previous
    // configuration
    if (is_active &&
        (build_flag ||
         compare_field_config_interp_method(
           curr_field_config - 1, curr_field_config))) {

      build_flag = 1;

      if (interp_weights != NULL) yac_interp_weights_delete(interp_weights);

      // generate interp weights
      interp_weights = generate_interp_weights(
        curr_field_config->src_interp_config, interp_grid);
    }

    if (is_active &&
        (curr_field_config->src_interp_config.weight_file_name != NULL)) {

      int src_comp_idx = field_configs[i].src_comp_idx;

      yac_interp_weights_write_to_file(
        interp_weights,
        curr_field_config->src_interp_config.weight_file_name,
        curr_comp_grid_pair->config[src_comp_idx].grid_name,
        curr_comp_grid_pair->config[src_comp_idx^1].grid_name);
    }

    // if the current weight reorder method differs from the previous
    // configuration
    // (use memcmp to compare frac_mask_fallback_value, because they can be nan)
    if (is_active &&
        (build_flag ||
         (curr_field_config[-1].reorder_type !=
          curr_field_config->reorder_type) ||
         (curr_field_config[-1].collection_size !=
          curr_field_config->collection_size) ||
         (curr_field_config[-1].scale_factor !=
          curr_field_config->scale_factor) ||
         (curr_field_config[-1].scale_summand !=
          curr_field_config->scale_summand) ||
         memcmp(
           &(curr_field_config[-1].frac_mask_fallback_value),
           &(curr_field_config->frac_mask_fallback_value),
           sizeof(double)))) {

      if (interp != NULL) yac_interpolation_delete(interp);

      // generate interpolation
      interp =
        yac_interp_weights_get_interpolation(
          interp_weights, curr_field_config->reorder_type,
          curr_field_config->collection_size,
          curr_field_config->frac_mask_fallback_value,
          curr_field_config->scale_factor,
          curr_field_config->scale_summand);
    }

    if (is_active) {

      int is_source = curr_field_config->src_interp_config.field != NULL;
      int is_target = curr_field_config->tgt_interp_config.field != NULL;

      struct interpolation * interp_copy = yac_interpolation_copy(interp);

      if (is_source) {
        yac_set_coupling_field_put_op(
          curr_field_config->src_interp_config.field,
          generate_event(
            curr_field_config->src_interp_config.event_data),
          interp_copy);
        yac_interpolation_inc_ref_count(interp_copy);
      }

      if (is_target) {
        yac_set_coupling_field_get_op(
          curr_field_config->tgt_interp_config.field,
          generate_event(
            curr_field_config->tgt_interp_config.event_data),
          interp_copy);
        yac_interpolation_inc_ref_count(interp_copy);
      }

      yac_interpolation_delete(interp_copy);

    }

    prev_comp_grid_pair = curr_comp_grid_pair;
  }

  yac_interpolation_delete(interp);
  yac_interp_weights_delete(interp_weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(dist_grid_pair);
  if (comp_pair_comm != MPI_COMM_NULL)
    yac_mpi_call(MPI_Comm_free(&comp_pair_comm), comm);
  free(field_configs);
}

void yac_instance_sync_def(struct yac_instance * instance) {
  CHECK_PHASE(
    yac_instance_sync_def, INSTANCE_DEFINITION_COMP, INSTANCE_DEFINITION_SYNC);
  YAC_ASSERT(instance->comp_config,
    "ERROR(yac_instance_sync_def): no components have been defined");
  yac_couple_config_sync(instance->couple_config, instance->comm);
}

void yac_instance_setup(struct yac_instance * instance) {
  // if definitions have not yet been synced
  int requires_def_sync = (instance->phase == INSTANCE_DEFINITION_COMP);
  if (requires_def_sync)
    yac_instance_sync_def(instance);
  CHECK_PHASE(yac_instance_setup, INSTANCE_DEFINITION_SYNC, INSTANCE_EXCHANGE);

  YAC_ASSERT(
    instance->comp_config,
    "ERROR(yac_instance_setup): no components have been defined");

  // if not all precesses have synced the definitions in this function
  // -> sync again
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &requires_def_sync, 1, MPI_INT, MPI_MIN, instance->comm),
      instance->comm);
  if(requires_def_sync == 0)
    yac_couple_config_sync(instance->couple_config, instance->comm);

  generate_interpolations(instance);
}

char * yac_instance_setup_and_emit_config(
  struct yac_instance * instance, int emit_flags) {

  yac_instance_setup(instance);

  return yac_yaml_emit_coupling(instance->couple_config, emit_flags);
}

MPI_Comm yac_instance_get_comps_comm(
  struct yac_instance * instance,
  char const ** comp_names, size_t num_comp_names) {
  CHECK_MIN_PHASE(yac_instance_get_comps_comm, INSTANCE_DEFINITION_COMP);
  return
    yac_component_config_get_comps_comm(
      instance->comp_config, comp_names, num_comp_names);
}

struct yac_instance * yac_instance_new(MPI_Comm comm) {

  struct yac_instance * instance = xmalloc(1 * sizeof(*instance));

  instance->couple_config = yac_couple_config_new();

  instance->comp_config = NULL;

  instance->cpl_fields = NULL;
  instance->num_cpl_fields = 0;

  yac_mpi_call(MPI_Comm_split(comm, 0, 0, &(instance->comm)), comm);

  instance->phase = INSTANCE_DEFINITION;

  return instance;
}


void yac_instance_dummy_new(MPI_Comm comm) {

  MPI_Comm dummy_comm;
  yac_mpi_call(MPI_Comm_split(comm, MPI_UNDEFINED, 0, &dummy_comm), comm);
}

void yac_instance_delete(struct yac_instance * instance) {

  if (instance == NULL) return;

  yac_component_config_delete(instance->comp_config);

  for (size_t i = 0; i < instance->num_cpl_fields; ++i)
    yac_coupling_field_delete(instance->cpl_fields[i]);
  free(instance->cpl_fields);

  yac_couple_config_delete(instance->couple_config);

  yac_mpi_call(MPI_Comm_free(&(instance->comm)), MPI_COMM_WORLD);

  free(instance);
}

struct yac_couple_config * yac_instance_get_couple_config(
  struct yac_instance * instance) {

  return instance->couple_config;
}

void yac_instance_set_couple_config(
  struct yac_instance * instance,
    struct yac_couple_config * couple_config) {
  CHECK_MAX_PHASE("yac_instance_set_couple_config", INSTANCE_DEFINITION);
  if (instance->couple_config == couple_config) return;
  yac_couple_config_delete(instance->couple_config);
  instance->couple_config = couple_config;
}

void yac_instance_def_datetime(
  struct yac_instance * instance, const char * start_datetime,
  const char * end_datetime ) {
  CHECK_MAX_PHASE("yac_instance_def_datetime", INSTANCE_DEFINITION_COMP);
  yac_couple_config_set_datetime(
    instance->couple_config, start_datetime, end_datetime);
}

char * yac_instance_get_start_datetime(struct yac_instance * instance) {
  CHECK_MIN_PHASE("yac_instance_get_start_datetime", INSTANCE_DEFINITION_COMP);
  return (char*)yac_couple_config_get_start_datetime(instance->couple_config);
}

char * yac_instance_get_end_datetime(struct yac_instance * instance) {
  CHECK_MIN_PHASE("yac_instance_get_end_datetime", INSTANCE_DEFINITION_COMP);
  return (char*)yac_couple_config_get_end_datetime(instance->couple_config);
}

void yac_instance_def_components(
  struct yac_instance * instance,
  char const ** comp_names, size_t num_comps) {
  CHECK_PHASE(
    yac_instance_def_components, INSTANCE_DEFINITION, INSTANCE_DEFINITION_COMP);

  YAC_ASSERT(
    !instance->comp_config,
    "ERROR(yac_instance_def_components): components have already been defined")

  // add components to coupling configuration
  for (size_t i = 0; i < num_comps; ++i)
    yac_couple_config_add_component(instance->couple_config, comp_names[i]);

  // synchronise coupling configuration
  yac_couple_config_sync(instance->couple_config, instance->comm);

  instance->comp_config =
    yac_component_config_new(
      instance->couple_config, comp_names, num_comps, instance->comm);

  int redirect = yac_couple_config_get_redirect_stdout(instance->couple_config);

  if ( redirect && (num_comps > 0)) {
    MPI_Comm comm =
      yac_component_get_communicator(
        yac_component_config_get_components(instance->comp_config)[0]);
    int rank, size;
    yac_mpi_call(MPI_Comm_rank(comm, &rank), comm);
    yac_mpi_call(MPI_Comm_size(comm, &size), comm);
    int redirect_only_root = (redirect != 1);
    yac_redirstdout(comp_names[0], redirect_only_root, rank, size);
  }
}

int yac_instance_components_are_defined(
  struct yac_instance * instance) {
  return instance->phase > INSTANCE_DEFINITION_COMP;
}

struct component * yac_instance_get_component(
  struct yac_instance * instance, size_t component_index) {
  CHECK_MIN_PHASE(yac_instance_get_component, INSTANCE_DEFINITION_COMP);
  YAC_ASSERT(
    yac_component_config_get_num_components(instance->comp_config) >
     component_index,
    "ERROR(yac_instance_get_component): invalid component_index")

  return
    yac_component_config_get_components(instance->comp_config)[component_index];
}

struct coupling_field * yac_instance_add_field(
  struct yac_instance * instance, char const * field_name,
  size_t component_index, struct yac_basic_grid * grid,
    struct interp_field * interp_fields, size_t num_interp_fields,
    int collection_size, char const * timestep) {
  CHECK_MIN_PHASE(yac_instance_get_component, INSTANCE_DEFINITION_COMP);
  CHECK_MAX_PHASE(yac_instance_get_component, INSTANCE_DEFINITION_SYNC);

  struct component * component =
    yac_instance_get_component(instance, component_index);

  struct yac_couple_config * couple_config = instance->couple_config;
  char const * grid_name = yac_basic_grid_get_name(grid);
  char const * component_name = yac_component_get_name(component);
  if(!yac_couple_config_contains_grid_name(couple_config, grid_name))
    yac_couple_config_add_grid(couple_config, grid_name);

  YAC_ASSERT(
    field_name, "ERROR(yac_instance_add_field): "
    "\"NULL\" is not a valid field name")
  YAC_ASSERT(
    strlen(field_name) <= YAC_MAX_CHARLEN,
    "ERROR(yac_instance_add_field): field name is too long "
    "(maximum is YAC_MAX_CHARLEN)")
  YAC_ASSERT_F(
    (collection_size > 0) && (collection_size < INT_MAX),
    "ERROR(yac_instance_add_field): \"%d\" is not a valid collection size "
    "(component \"%s\" grid \"%s\" field \"%s\")",
    collection_size, component_name, grid_name, field_name)

  // add field to coupling configuration
  yac_couple_config_component_add_field(
    couple_config, component_name, grid_name, field_name,
    timestep, collection_size);

  // check whether the field is already defined
  for (size_t i = 0; i < instance->num_cpl_fields; ++i) {
    struct coupling_field * cpl_field = instance->cpl_fields[i];
    YAC_ASSERT_F(
      strcmp(yac_get_coupling_field_name(cpl_field), field_name) ||
      (yac_get_coupling_field_component(cpl_field) != component) ||
      (yac_coupling_field_get_basic_grid(cpl_field) != grid),
      "ERROR(yac_instance_add_field): "
      "field with the name \"%s\" has already been defined",
      field_name);
  }

  struct coupling_field * cpl_field =
    yac_coupling_field_new(
      field_name, component, grid, interp_fields, num_interp_fields,
        collection_size, timestep);

  instance->cpl_fields =
    xrealloc(
      instance->cpl_fields,
      (instance->num_cpl_fields + 1) * sizeof(*(instance->cpl_fields)));
  instance->cpl_fields[instance->num_cpl_fields] = cpl_field;
  instance->num_cpl_fields++;

  return cpl_field;
}

void yac_instance_def_couple(
  struct yac_instance * instance,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_period, int time_reduction,
  struct yac_interp_stack_config * interp_stack_config, int src_lag, int tgt_lag,
  const char* weight_file_name, int mapping_on_source,
  double scale_factor, double scale_summand, size_t num_src_mask_names,
  char const * const * src_mask_names, char const * tgt_mask_name) {

  CHECK_MIN_PHASE(yac_instance_def_couple, INSTANCE_DEFINITION);
  CHECK_MAX_PHASE(yac_instance_def_couple, INSTANCE_DEFINITION_SYNC);

  yac_couple_config_def_couple(
    instance->couple_config, src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name, coupling_period,
    time_reduction, interp_stack_config, src_lag, tgt_lag, weight_file_name,
    mapping_on_source, scale_factor, scale_summand, num_src_mask_names,
    src_mask_names, tgt_mask_name);
}

struct coupling_field* yac_instance_get_field(struct yac_instance * instance,
  const char * comp_name, const char* grid_name, const char * field_name){
  CHECK_MIN_PHASE("yac_instance_get_field", INSTANCE_DEFINITION_COMP);
  return get_coupling_field(comp_name, field_name,
    grid_name, instance->num_cpl_fields, instance->cpl_fields);
}
