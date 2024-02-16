/**
 * @file couple_config.c
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

#include "yac_mpi.h"
#include "couple_config.h"
#include "yac_interface.h"
#include "interp_method_avg.h"
#include "interp_method_check.h"
#include "interp_method_conserv.h"
#include "interp_method_creep.h"
#include "interp_method_nnn.h"
#include "interp_method_fixed.h"
#include "interp_method_file.h"
#include "interp_method_hcsbb.h"
#include "interp_method_spmap.h"
#include "interp_method_callback.h"
#include "mtime_calendar.h"

struct yac_couple_config_grid {
  char const * name;
  char* metadata;
};

struct yac_couple_config_field {
  char const * name;
  size_t grid_idx;
  char const * timestep;
  char * metadata;
  double frac_mask_fallback_value;
  size_t collection_size;
};

struct yac_couple_config_component {

  struct yac_couple_config_field * fields;
  size_t num_fields;

  char const * name;
  char * metadata;
};

struct yac_couple_config_field_couple {

  struct {
    size_t component_idx;
    size_t field_idx;
    int lag;
  } source, target;

  int mapping_on_source;

  struct yac_interp_stack_config * interp_stack;

  enum yac_reduction_type coupling_period_operation;
  char const * coupling_period;

  int enforce_write_weight_file;
  char const * weight_file_name;

  double scale_factor, scale_summand;

  char ** src_mask_names;
  size_t num_src_mask_names;
  char * tgt_mask_name;
};

struct yac_couple_config_couple {

  size_t component_indices[2];

  struct yac_couple_config_field_couple * field_couples;
  size_t num_field_couples;
};

struct yac_couple_config {

  char const * start_datetime;
  char const * end_datetime;

  struct yac_couple_config_grid * grids;
  size_t num_grids;

  struct yac_couple_config_component * components;
  size_t num_components;

  struct yac_couple_config_couple * couples;
  size_t num_couples;

  int redirect_stdout;
};

struct yac_couple_config * yac_couple_config_new() {

  struct yac_couple_config * couple_config =
    xmalloc(1 * sizeof(*couple_config));

  couple_config->start_datetime = NULL;
  couple_config->end_datetime = NULL;

  couple_config->grids = NULL;
  couple_config->num_grids = 0;

  couple_config->components = NULL;
  couple_config->num_components = 0;

  couple_config->couples = NULL;
  couple_config->num_couples = 0;

  couple_config->redirect_stdout = 0;

  return couple_config;
}

static char const * string_dup(char const * string) {
  return (string != NULL)?strdup(string):NULL;
}

static void yac_couple_config_field_free(void * field_){
  struct yac_couple_config_field * field = field_;
  free((void*)field->name);
  free((void*)field->timestep);
  free(field->metadata);
}


static void yac_couple_config_component_free(void * component_) {

  struct yac_couple_config_component * component = component_;

  for (size_t i = 0; i < component->num_fields; ++i)
    yac_couple_config_field_free(component->fields + i);
  free(component->fields);
  free((void*)(component->name));
  free(component->metadata);
}

static void yac_couple_config_field_couple_free(void * field_couple_) {

  struct yac_couple_config_field_couple * field_couple = field_couple_;
  yac_interp_stack_config_delete(field_couple->interp_stack);
  free((void*)field_couple->coupling_period);
  free((void*)field_couple->weight_file_name);
  for (size_t i = 0; i < field_couple->num_src_mask_names; ++i)
    free((void*)field_couple->src_mask_names[i]);
  free(field_couple->src_mask_names);
  free(field_couple->tgt_mask_name);
}

static void yac_couple_config_couple_free(void * couple_) {

  struct yac_couple_config_couple * couple = couple_;
  for (size_t i = 0; i < couple->num_field_couples; ++i)
    yac_couple_config_field_couple_free(couple->field_couples + i);
  free(couple->field_couples);
  couple->field_couples = NULL;
  couple->num_field_couples = 0;
}

static void yac_couple_config_grid_free(void * grid_) {

  struct yac_couple_config_grid * grid = grid_;

  free((void*)grid->name);
  free(grid->metadata);
}

static int yac_couple_config_grid_compare(void const * a, void const * b) {
  YAC_ASSERT(((struct yac_couple_config_grid const *)a)->name != NULL &&
             ((struct yac_couple_config_grid const *)b)->name != NULL,
    "ERROR(yac_couple_config_grid_compare): "
    "invalid name (NULL is not allowed)");
  return strcmp(((struct yac_couple_config_grid const *)a)->name,
                ((struct yac_couple_config_grid const *)b)->name);
}

static int yac_couple_config_field_compare(void const * a_, void const * b_) {
  struct yac_couple_config_field const * a = a_;
  struct yac_couple_config_field const * b = b_;
  int ret = (a->grid_idx > b->grid_idx) - (a->grid_idx < b->grid_idx);
  if (ret) return ret;
  YAC_ASSERT(
    a->name != NULL && b->name != NULL,
    "ERROR(yac_couple_config_field_compare): "
    "invalid name (NULL is not allowed)");
  return strcmp(a->name, b->name);
}

static int yac_couple_config_component_compare(
  void const * a_, void const * b_) {
  struct yac_couple_config_component const * a = a_;
  struct yac_couple_config_component const * b = b_;
  YAC_ASSERT(a->name != NULL && b->name != NULL,
    "ERROR(yac_couple_config_component_compare): "
    "invalid name (NULL is not allowed)");
  return strcmp(a->name, b->name);
}

static int yac_couple_config_field_couple_compare(
  void const * a_, void const * b_) {
  struct yac_couple_config_field_couple const * a = a_;
  struct yac_couple_config_field_couple const * b = b_;
  int ret;
  ret = (a->source.component_idx > b->source.component_idx) -
        (a->source.component_idx < b->source.component_idx);
  if (ret) return ret;
  ret = (a->target.component_idx > b->target.component_idx) -
        (a->target.component_idx < b->target.component_idx);
  if (ret) return ret;
  ret = (a->source.field_idx > b->source.field_idx) -
        (a->source.field_idx < b->source.field_idx);
  if (ret) return ret;
  return (a->target.field_idx > b->target.field_idx) -
         (a->target.field_idx < b->target.field_idx);
}

static int yac_couple_config_couple_compare_basic(
  void const * a_, void const * b_) {
  struct yac_couple_config_couple const * a = a_;
  struct yac_couple_config_couple const * b = b_;
  int ret;
  ret = (a->component_indices[0] > b->component_indices[0]) -
        (a->component_indices[0] < b->component_indices[0]);
  if (ret) return ret;
  return (a->component_indices[1] > b->component_indices[1]) -
         (a->component_indices[1] < b->component_indices[1]);
}

static void couple_config_sync_string(
  char const * string_name, char ** string, MPI_Comm comm) {

  int rank;
  MPI_Comm_rank(comm, &rank);
  struct {
    int len;
    int rank;
  } data_pair;
  size_t len = (*string != NULL)?(strlen(*string) + 1):0;
  YAC_ASSERT_F(
    len <= INT_MAX,
    "ERROR(couple_config_sync_string): \"%s\" too long", string_name);
  data_pair.len = (int)len;
  data_pair.rank = rank;

  // determine the broadcaster
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &data_pair, 1, MPI_2INT, MPI_MAXLOC, comm), comm);
  if (data_pair.len == 0) return;

  // broadcast string
  char const * string_bak = NULL;
  if (data_pair.rank != rank) {
    string_bak = *string;
    *string = xmalloc((size_t)data_pair.len * sizeof(**string));
  }
  yac_mpi_call(
    MPI_Bcast(*string, data_pair.len, MPI_CHAR, data_pair.rank, comm), comm);

  // check for consistency
  if (data_pair.rank != rank) {
    YAC_ASSERT_F(
      (string_bak == NULL) ||
      !strcmp(string_bak, *string),
      "ERROR(couple_config_sync_string): inconsistent \"%s\" definition "
      "(\"%s\" != \"%s\")", string_name, string_bak, *string);
    free((void*)string_bak);
  }
}

static void couple_config_sync_calendar(MPI_Comm comm) {

  int calendar = (int)getCalendarType();
  if (calendar == CALENDAR_NOT_SET) calendar = INT_MAX;

  // broadcast calendar
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &calendar, 1, MPI_INT, MPI_MIN, comm), comm);

  // if no process has defined a calendar
  if (calendar == INT_MAX) return;

  // set calendar (if local process has already defined a calendar,
  // the definition is checked for consistency)
  yac_cdef_calendar(calendar);
}

static void yac_couple_config_component_merge(
  void * a_, void * b_, MPI_Comm comm) {

  struct yac_couple_config_component * a = a_;
  struct yac_couple_config_component * b = b_;

  if (b) {
    a->num_fields = b->num_fields;
    a->fields = b->fields;
    b->num_fields = 0;
    b->fields = NULL;
    a->metadata = b->metadata;
    b->metadata = NULL;
  }

  couple_config_sync_string(
    "component metadata", (char **)&(a->metadata), comm);
}

static void yac_couple_config_grid_merge(
  void * a_, void * b_, MPI_Comm comm) {

  struct yac_couple_config_grid * a = a_;
  struct yac_couple_config_grid * b = b_;

  if (b!= NULL) {
    a->metadata = b->metadata;
    b->metadata = NULL;
  }

  couple_config_sync_string("grid metadata", (char **)&(a->metadata), comm);
}

static void yac_couple_config_field_merge(
  void * a_, void * b_, MPI_Comm comm) {

  struct yac_couple_config_field * a = a_;
  struct yac_couple_config_field * b = b_;

  enum {
    TIMESTEP_IDX,
    METADATA_IDX,
    FRAC_MASK_IDX,
    COLLECTION_SIZE_IDX,
    DATA_PAIR_COUNT
  };

  int rank;
  MPI_Comm_rank(comm, &rank);
  struct {
    int len;
    int rank;
  } data_pairs[DATA_PAIR_COUNT];
  size_t timestep_len =
    ((b != NULL) && (b->timestep != NULL))?(strlen(b->timestep) + 1):0;
  size_t metadata_len =
    ((b != NULL) && (b->metadata != NULL))?(strlen(b->metadata) + 1):0;
  YAC_ASSERT(
    timestep_len <= INT_MAX,
    "ERROR(yac_couple_config_field_merge): timestep string too long");
  YAC_ASSERT(
    metadata_len <= INT_MAX,
    "ERROR(yac_couple_config_field_merge): metadata string too long");
  YAC_ASSERT_F(
    (b == NULL) || (b->collection_size <= INT_MAX) ||
    (b->collection_size == SIZE_MAX),
    "ERROR(yac_couple_config_field_merge): invalid collection size \"%zu\"",
    b->collection_size);
  data_pairs[TIMESTEP_IDX].len = (int)timestep_len;
  data_pairs[TIMESTEP_IDX].rank = rank;
  data_pairs[METADATA_IDX].len = (int)metadata_len;
  data_pairs[METADATA_IDX].rank = rank;
  data_pairs[FRAC_MASK_IDX].len =
    (b != NULL) && (b->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE);
  data_pairs[FRAC_MASK_IDX].rank = rank;
  data_pairs[COLLECTION_SIZE_IDX].len =
    ((b == NULL) || (b->collection_size == SIZE_MAX))?
      -1:(int)b->collection_size;
  data_pairs[COLLECTION_SIZE_IDX].rank = rank;

  // determine the broadcaster
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, data_pairs, DATA_PAIR_COUNT, MPI_2INT,
      MPI_MAXLOC, comm), comm);

  // if at least one process has a valid timestep
  if (data_pairs[TIMESTEP_IDX].len > 0) {

    // broadcast timestep
    char * timestep_buffer = xmalloc((size_t)data_pairs[TIMESTEP_IDX].len);
    if (data_pairs[TIMESTEP_IDX].rank == rank)
      strcpy(timestep_buffer, b->timestep);
    yac_mpi_call(
      MPI_Bcast(
        timestep_buffer, data_pairs[TIMESTEP_IDX].len, MPI_CHAR,
        data_pairs[TIMESTEP_IDX].rank, comm), comm);

    // check for consistency
    YAC_ASSERT_F(
      (b == NULL) || (b->timestep == NULL) ||
      !strcmp(b->timestep, timestep_buffer),
      "ERROR(yac_couple_config_field_merge): "
      "inconsistent timestep definition (\"%s\" != \"%s\")",
      b->timestep, timestep_buffer);

      // update timestep
      free((void*)(a->timestep));
      a->timestep = timestep_buffer;
  }

  // if at least one process has a valid metadata
  if (data_pairs[METADATA_IDX].len > 0) {

    // broadcast metadata
    char * metadata_buffer = xmalloc((size_t)data_pairs[METADATA_IDX].len);
    if (data_pairs[METADATA_IDX].rank == rank)
      strcpy(metadata_buffer, b->metadata);
    yac_mpi_call(
      MPI_Bcast(
        metadata_buffer, data_pairs[METADATA_IDX].len, MPI_CHAR,
        data_pairs[METADATA_IDX].rank, comm), comm);

    // check for consistency
    YAC_ASSERT_F(
      (b == NULL) || (b->metadata == NULL) ||
      !strcmp(b->metadata, metadata_buffer),
      "ERROR(yac_couple_config_field_merge): "
      "inconsistent metadata definition (\"%s\" != \"%s\")",
      b->metadata, metadata_buffer);

    // update metadata
    free(a->metadata);
    a->metadata = metadata_buffer;
  }

  // if at least one process has a valid fractional mask fallback value
  if (data_pairs[FRAC_MASK_IDX].len != 0) {

    // broadcast fractional mask fallback value
    double frac_mask_fallback_value;
    if (data_pairs[FRAC_MASK_IDX].rank == rank)
      frac_mask_fallback_value = b->frac_mask_fallback_value;
    yac_mpi_call(
      MPI_Bcast(
        &frac_mask_fallback_value, 1, MPI_DOUBLE,
        data_pairs[FRAC_MASK_IDX].rank, comm), comm);

    // check for consistency
    // (use memcmp for comparison, because it can be nan)
    YAC_ASSERT_F(
      (b == NULL) || (b->frac_mask_fallback_value == YAC_FRAC_MASK_NO_VALUE) ||
      !memcmp(
        &b->frac_mask_fallback_value, &frac_mask_fallback_value,
        sizeof(double)),
      "ERROR(yac_couple_config_field_merge): "
      "inconsistent fractional mask fallback value definition "
      "(%lf != %lf)",
      b->frac_mask_fallback_value, frac_mask_fallback_value);

    // update fractional mask fallback value
    a->frac_mask_fallback_value = frac_mask_fallback_value;
  }


  // if at least one process has a valid collection size
  if (data_pairs[COLLECTION_SIZE_IDX].len > -1) {

    // check for consistency
    YAC_ASSERT_F(
      (b == NULL) || (b->collection_size == SIZE_MAX) ||
      (b->collection_size == (size_t)(data_pairs[COLLECTION_SIZE_IDX].len)),
      "ERROR(yac_couple_config_field_merge): "
      "inconsistent collection size definition (\"%zu\" != \"%d\")",
      b->collection_size, data_pairs[COLLECTION_SIZE_IDX].len);

    // update collection size
    a->collection_size = (size_t)(data_pairs[COLLECTION_SIZE_IDX].len);
  }
}

static void yac_couple_config_field_couple_merge(
  void * a_, void * b_, MPI_Comm comm) {

  if (b_ == NULL) return;

  struct yac_couple_config_field_couple * a = a_;
  struct yac_couple_config_field_couple * b = b_;

  YAC_ASSERT_F(
    (a->source.lag == b->source.lag),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent source lag (%d != %d)", a->source.lag, b->source.lag)
  YAC_ASSERT_F(
    (a->target.lag == b->target.lag),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent target lag (%d != %d)", a->target.lag, b->target.lag)
  YAC_ASSERT_F(
    (a->mapping_on_source == b->mapping_on_source),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent mapping side (%d != %d)",
    a->mapping_on_source, b->mapping_on_source)
  YAC_ASSERT(
    !yac_interp_stack_config_compare(a->interp_stack, b->interp_stack),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent interpolation stack")
  YAC_ASSERT_F(
    a->coupling_period_operation == b->coupling_period_operation,
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent coupling period operation (%d != %d)",
    (int)(a->coupling_period_operation), (int)(b->coupling_period_operation))
  YAC_ASSERT_F(
    !strcmp(a->coupling_period, b->coupling_period),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent coupling period (%s != %s)",
    a->coupling_period, b->coupling_period)
  YAC_ASSERT_F(
    a->enforce_write_weight_file == b->enforce_write_weight_file,
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent defintion of enforce_write_weight_file (%d != %d)",
    a->enforce_write_weight_file, b->enforce_write_weight_file)
  YAC_ASSERT_F(
    !a->enforce_write_weight_file ||
    !strcmp(a->weight_file_name, b->weight_file_name),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent weight_file_name (%s != %s)",
    a->weight_file_name, b->weight_file_name)
  YAC_ASSERT_F(
    a->scale_factor == b->scale_factor,
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent scale factor (%lf != %lf)",
    a->scale_factor, b->scale_factor)
  YAC_ASSERT_F(
    a->scale_summand == b->scale_summand,
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent scale summand (%lf != %lf)",
    a->scale_summand, b->scale_summand)
  YAC_ASSERT_F(
    a->num_src_mask_names == b->num_src_mask_names,
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent number of source mask names (%zu != %zu)",
    a->num_src_mask_names, b->num_src_mask_names)
  YAC_ASSERT_F(
    (a->src_mask_names != NULL) == (b->src_mask_names != NULL),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent availability of source mask names (%s != %s)",
    (a->src_mask_names != NULL)?"TRUE":"FALSE",
    (b->src_mask_names != NULL)?"TRUE":"FALSE")
  for (size_t i = 0; i < a->num_src_mask_names; ++i)
    YAC_ASSERT_F(
      !strcmp(a->src_mask_names[i], b->src_mask_names[i]),
      "ERROR(yac_couple_config_field_couple_merge): "
      "inconsistent source mask names at index %zu (\"%s\" != \"%s\")",
      i, a->src_mask_names[i], b->src_mask_names[i])
  YAC_ASSERT_F(
    (a->tgt_mask_name != NULL) == (b->tgt_mask_name != NULL),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent availability of target mask name (%s != %s)",
    (a->tgt_mask_name != NULL)?"TRUE":"FALSE",
    (b->tgt_mask_name != NULL)?"TRUE":"FALSE")
  YAC_ASSERT_F(
    (a->tgt_mask_name == NULL) ||
    !strcmp(a->tgt_mask_name, b->tgt_mask_name),
    "ERROR(yac_couple_config_field_couple_merge): "
    "inconsistent target mask name (\"%s\" != \"%s\")",
    a->tgt_mask_name, b->tgt_mask_name)
}

static void merge_field_couples(
  size_t * num_field_couples,
  struct yac_couple_config_field_couple ** field_couples, MPI_Comm comm);

static void yac_couple_config_couple_merge(
  void * a_, void * b_, MPI_Comm comm) {

  struct yac_couple_config_couple * a = a_;
  struct yac_couple_config_couple * b = b_;

  if (b) {
    a->num_field_couples = b->num_field_couples;
    a->field_couples = b->field_couples;
    b->num_field_couples = 0;
    b->field_couples = NULL;
  }

  // distribute and merge field couples
  merge_field_couples(&(a->num_field_couples), &(a->field_couples), comm);
}

void yac_couple_config_delete(struct yac_couple_config * couple_config) {

  free((void*)couple_config->start_datetime);
  free((void*)couple_config->end_datetime);

  for (size_t i = 0; i < couple_config->num_grids; ++i){
    yac_couple_config_grid_free(couple_config->grids + i);
  }
  free(couple_config->grids);

  for (size_t i = 0;
       i < couple_config->num_components; ++i)
    yac_couple_config_component_free(couple_config->components + i);
  free(couple_config->components);

  for (size_t couple_idx = 0; couple_idx < couple_config->num_couples;
       ++couple_idx)
    yac_couple_config_couple_free(couple_config->couples + couple_idx);
  free(couple_config->couples);
  free(couple_config);
}

static size_t yac_couple_config_add_grid_(
  struct yac_couple_config * couple_config, char const * name) {

  YAC_ASSERT(
    name != NULL,
    "ERROR(yac_couple_config_add_grid_): "
    "invalid name (NULL is not allowed)")

  for (size_t i = 0; i < couple_config->num_grids; ++i)
    if (!strcmp(couple_config->grids[i].name, name)) return i;

  size_t grid_idx = couple_config->num_grids;
  couple_config->num_grids++;
  couple_config->grids =
    xrealloc(
      couple_config->grids,
      couple_config->num_grids * sizeof(*(couple_config->grids)));
  couple_config->grids[grid_idx].name = string_dup(name);
  couple_config->grids[grid_idx].metadata = NULL;
  return grid_idx;
}

void yac_couple_config_add_grid(
  struct yac_couple_config * couple_config, char const * name) {

  yac_couple_config_add_grid_(couple_config, name);
}

static size_t yac_couple_config_add_component_(
  struct yac_couple_config * couple_config, char const * name) {

  YAC_ASSERT(
    name != NULL,
    "ERROR(yac_couple_config_add_component_): "
    "invalid name (NULL is not allowed)")

  for (size_t i = 0; i < couple_config->num_components; ++i)
    if (!strcmp(couple_config->components[i].name, name)) return i;

  size_t component_idx = couple_config->num_components;
  couple_config->num_components++;
  couple_config->components =
    xrealloc(
      couple_config->components,
      couple_config->num_components * sizeof(*(couple_config->components)));
  couple_config->components[component_idx].fields = NULL;
  couple_config->components[component_idx].num_fields = 0;
  couple_config->components[component_idx].name = strdup(name);
  couple_config->components[component_idx].metadata = NULL;
  return component_idx;
}

void yac_couple_config_add_component(
  struct yac_couple_config * couple_config, char const * name) {

  yac_couple_config_add_component_(couple_config, name);
}

void yac_couple_config_component_set_metadata(
  struct yac_couple_config * couple_config,
    char const * comp_name, const char* metadata) {
  size_t comp_idx = yac_couple_config_get_component_idx(couple_config, comp_name);
  free(couple_config->components[comp_idx].metadata);
  couple_config->components[comp_idx].metadata
    = metadata==NULL?NULL:strdup(metadata);
}

void yac_couple_config_grid_set_metadata(
  struct yac_couple_config * couple_config,
    char const * grid_name, const char* metadata) {
  if(!yac_couple_config_contains_grid_name(couple_config, grid_name))
    yac_couple_config_add_grid(couple_config, grid_name);
  size_t grid_idx = yac_couple_config_get_grid_idx(couple_config, grid_name);
  free(couple_config->grids[grid_idx].metadata);
  couple_config->grids[grid_idx].metadata
    = metadata==NULL?NULL:strdup(metadata);
}

void yac_couple_config_field_set_metadata(
  struct yac_couple_config * couple_config,
    const char* comp_name, const char * grid_name, const char* field_name,
    const char* metadata) {
  size_t comp_idx = yac_couple_config_get_component_idx(couple_config, comp_name);
  size_t grid_idx = yac_couple_config_get_grid_idx(couple_config, grid_name);
  size_t field_idx = yac_couple_config_get_field_idx(couple_config,
    comp_idx, grid_idx, field_name);
  free(couple_config->components[comp_idx].fields[field_idx].metadata);
  couple_config->components[comp_idx].fields[field_idx].metadata
    = metadata==NULL?NULL:strdup(metadata);
}

const char* yac_couple_config_component_get_metadata(
  struct yac_couple_config * couple_config,
    const char * comp_name) {
  size_t comp_idx = yac_couple_config_get_component_idx(couple_config, comp_name);
  return couple_config->components[comp_idx].metadata;
}

const char* yac_couple_config_grid_get_metadata(
  struct yac_couple_config * couple_config,
    const char * grid_name) {
  size_t grid_idx = yac_couple_config_get_grid_idx(couple_config, grid_name);
  return couple_config->grids[grid_idx].metadata;
}

const char* yac_couple_config_field_get_metadata(
  struct yac_couple_config * couple_config,
    const char* comp_name, const char * grid_name, const char* field_name) {
  size_t comp_idx = yac_couple_config_get_component_idx(couple_config, comp_name);
  size_t grid_idx = yac_couple_config_get_grid_idx(couple_config, grid_name);
  size_t field_idx = yac_couple_config_get_field_idx(couple_config,
    comp_idx, grid_idx, field_name);
  return couple_config->components[comp_idx].fields[field_idx].metadata;
}

static void check_component_idx(
  struct yac_couple_config * couple_config, size_t component_idx,
  char const * routine_name, int line) {

  YAC_ASSERT_F(
    component_idx < couple_config->num_components,
    "ERROR(%s:%d:%s): invalid component_idx", __FILE__, line, routine_name)
}

static void check_grid_idx(
  struct yac_couple_config * couple_config, size_t grid_idx,
  char const * routine_name, int line) {

  YAC_ASSERT_F(
    grid_idx < couple_config->num_grids,
    "ERROR(%s:%d:%s): invalid grid_idx", __FILE__, line, routine_name)
}

static size_t yac_couple_config_component_add_field_(
  struct yac_couple_config * couple_config,
  size_t comp_idx, size_t grid_idx, char const * name,
  char const * timestep, size_t collection_size) {

  check_component_idx(
    couple_config, comp_idx,
    "yac_couple_config_component_add_field_", __LINE__);
  check_grid_idx(
    couple_config, grid_idx,
    "yac_couple_config_component_add_field_", __LINE__);

  // check whether the field already exists
  struct yac_couple_config_component * component =
    couple_config->components + comp_idx;
  for (size_t i = 0; i < component->num_fields; i++) {
    if((strcmp(component->fields[i].name, name) == 0) &&
       (component->fields[i].grid_idx == grid_idx)) {
      // if no timestep is defined for the field
      if (!component->fields[i].timestep && timestep)
        component->fields[i].timestep = strdup(timestep);
      // if no collection size is defined for the field
      if (component->fields[i].collection_size == SIZE_MAX)
        component->fields[i].collection_size = collection_size;
      YAC_ASSERT_F(
        !timestep ||
        !strcmp(timestep, component->fields[i].timestep),
        "ERROR(yac_couple_config_component_add_field): "
        "inconsistent timestep definition (\"%s\" != \"%s\")",
        timestep, component->fields[i].timestep);
      YAC_ASSERT_F(
        collection_size == SIZE_MAX ||
        collection_size == component->fields[i].collection_size,
        "ERROR(yac_couple_config_component_add_field): "
        "inconsistent collection_size definition (%zu != %zu)",
        collection_size, component->fields[i].collection_size);
      return i;
    }
  }

  size_t field_idx = component->num_fields;
  component->num_fields++;

  component->fields =
    xrealloc(
      component->fields,
      component->num_fields *
      sizeof(*(component->fields)));
  struct yac_couple_config_field * field =
    component->fields + field_idx;

  field->name = strdup(name);
  field->grid_idx = grid_idx;
  field->timestep = timestep?strdup(timestep):NULL;
  field->metadata = NULL;
  field->frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
  field->collection_size = collection_size;
  return field_idx;
}

void yac_couple_config_component_add_field(
  struct yac_couple_config * couple_config, const char* component_name,
    const char* grid_name, const char* name, char const * timestep,
    size_t collection_size) {

  yac_couple_config_component_add_field_(
    couple_config,
    yac_couple_config_get_component_idx(couple_config, component_name),
    yac_couple_config_get_grid_idx(couple_config, grid_name),
    name, timestep, collection_size);
}

size_t yac_couple_config_get_num_couples(
  struct yac_couple_config * couple_config) {

  return couple_config->num_couples;
}

static void check_couple_idx(
  struct yac_couple_config * couple_config, size_t couple_idx,
  char const * routine_name, int line) {

  YAC_ASSERT_F(
    couple_idx < couple_config->num_couples,
    "ERROR(%s:%d:%s): invalid couple_idx", __FILE__, line, routine_name);
}

size_t yac_couple_config_get_num_couple_fields(
  struct yac_couple_config * couple_config, size_t couple_idx) {

  check_couple_idx(
    couple_config, couple_idx, "yac_couple_config_get_num_couple_fields",
    __LINE__);

  return couple_config->couples[couple_idx].num_field_couples;
}

void yac_couple_config_get_couple_component_names(
  struct yac_couple_config * couple_config, size_t couple_idx,
  char const * couple_component_names[2]) {

  check_couple_idx(
    couple_config, couple_idx, "yac_couple_config_get_couple_component_names",
    __LINE__);

  for (int i = 0; i < 2; ++i)
    couple_component_names[i] =
      couple_config->components[
        couple_config->couples[couple_idx].component_indices[i]].name;
}

int yac_couple_config_component_name_is_valid(
  struct yac_couple_config * couple_config, char const * component_name) {

  YAC_ASSERT(
    component_name,
    "ERROR(yac_couple_config_component_name_is_valid): component name is NULL")
  YAC_ASSERT(
    strlen(component_name) <= YAC_MAX_CHARLEN,
    "ERROR(yac_couple_config_component_name_is_valid): "
    "component name is too long (maximum is YAC_MAX_CHARLEN)")

  for (size_t component_idx = 0; component_idx < couple_config->num_components;
       ++component_idx)
    if (!strcmp(component_name, couple_config->components[component_idx].name))
      return 1;
  return 0;
}

size_t yac_couple_config_get_num_components(
  struct yac_couple_config * couple_config) {

  return couple_config->num_components;
}


size_t yac_couple_config_get_num_grids(
  struct yac_couple_config * couple_config) {

  return couple_config->num_grids;
}

size_t yac_couple_config_get_num_fields(
  struct yac_couple_config * couple_config, size_t component_idx) {
  check_component_idx(couple_config, component_idx,
    "yac_couple_config_get_num_fields", __LINE__);
  return couple_config->components[component_idx].num_fields;
}

size_t yac_couple_config_get_component_idx(
  struct yac_couple_config * couple_config, char const * component_name) {

  size_t component_idx = SIZE_MAX;
  for (size_t i = 0; (i < couple_config->num_components) &&
                     (component_idx == SIZE_MAX); ++i)
    if (!strcmp(couple_config->components[i].name, component_name))
      component_idx = i;

  YAC_ASSERT_F(
    component_idx != SIZE_MAX,
    "ERROR(yac_couple_config_get_component_idx): "
    "Component \"%s\" not found in coupling configuration",
    component_name);

  return component_idx;
}

size_t yac_couple_config_get_grid_idx(
  struct yac_couple_config * couple_config, char const * grid_name) {

  size_t grid_idx = SIZE_MAX;
  for (size_t i = 0;
       (i < couple_config->num_grids) && (grid_idx == SIZE_MAX); ++i)
    if (!strcmp(couple_config->grids[i].name, grid_name))
      grid_idx = i;

  YAC_ASSERT_F(
    grid_idx != SIZE_MAX,
    "ERROR(yac_couple_config_get_grid_idx): "
    "grid name \"%s\" not in list of grids", grid_name)

  return grid_idx;
}

size_t yac_couple_config_get_field_idx(
  struct yac_couple_config * couple_config, size_t component_idx,
    size_t grid_idx, char const * field_name) {
  check_component_idx(
    couple_config, component_idx,
    "yac_couple_config_get_component_name", __LINE__);

  size_t field_idx = SIZE_MAX;
  struct yac_couple_config_component * component =
    couple_config->components + component_idx;
  size_t nbr_fields = component->num_fields;
  for(size_t i=0;
      (i<nbr_fields) && (field_idx == SIZE_MAX); ++i)
    if((component->fields[i].grid_idx == grid_idx) &&
       !strcmp(
          field_name,
            component->fields[i].name))
      field_idx = i;

  YAC_ASSERT_F(
    field_idx != INT_MAX,
    "ERROR(yac_couple_config_get_field_idx): "
    "field not found "
    "(component_idx %zu grid_idx %zu field_name \"%s\"",
    component_idx, grid_idx, field_name);

  return field_idx;
}

char const * yac_couple_config_get_component_name(
  struct yac_couple_config * couple_config, size_t component_idx) {

  check_component_idx(
    couple_config, component_idx,
    "yac_couple_config_get_component_name", __LINE__);

  return couple_config->components[component_idx].name;
}

static void check_field_idx(
  struct yac_couple_config * couple_config, size_t component_idx,
  size_t field_idx, char const * routine_name, int line) {

  check_component_idx(
    couple_config, component_idx, routine_name, __LINE__);

  YAC_ASSERT_F(
    field_idx <
      couple_config->components[component_idx].num_fields,
    "ERROR(%s:%d:%s): invalid field_idx", __FILE__,
    line, routine_name)
}

char const * yac_couple_config_get_field_grid_name(
  struct yac_couple_config * couple_config, size_t component_idx,
  size_t field_idx) {

  check_field_idx(
    couple_config, component_idx, field_idx,
    "yac_couple_config_get_field_grid_name", __LINE__);

  return
    couple_config->grids[
      couple_config->
        components[component_idx].
        fields[field_idx].grid_idx].name;
}

char const * yac_couple_config_get_field_name(
  struct yac_couple_config * couple_config, size_t component_idx,
  size_t field_idx) {

  check_field_idx(
    couple_config, component_idx, field_idx,
    "yac_couple_config_get_field_name", __LINE__);

  return
      couple_config->
        components[component_idx].
        fields[field_idx].name;
}

char const * yac_couple_config_get_field_timestep(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name) {

  size_t component_idx =
    yac_couple_config_get_component_idx(couple_config, component_name);
  size_t grid_idx =
    yac_couple_config_get_grid_idx(couple_config, grid_name);
  size_t field_idx =
    yac_couple_config_get_field_idx(
      couple_config, component_idx, grid_idx, field_name);

  struct yac_couple_config_field * field =
    couple_config->components[component_idx].fields +
    field_idx;

  YAC_ASSERT_F(
    field->timestep,
    "ERROR(yac_couple_config_get_field_timestep): "
    "no valid timestep defined (component: \"%s\" field \"%s\")",
    couple_config->components[component_idx].name, field->name);

  return field->timestep;
}

int yac_couple_config_get_field_role(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name) {

  size_t component_idx =
    yac_couple_config_get_component_idx(couple_config, component_name);
  size_t grid_idx =
    yac_couple_config_get_grid_idx(couple_config, grid_name);
  size_t field_idx =
    yac_couple_config_get_field_idx(
      couple_config, component_idx, grid_idx, field_name);

  size_t nbr_couples = couple_config->num_couples;
  for(size_t couple_idx = 0; couple_idx<nbr_couples; ++couple_idx){
    struct yac_couple_config_couple * couple = couple_config->couples + couple_idx;
    if(couple->component_indices[0] != component_idx &&
       couple->component_indices[1] != component_idx)
      continue;
    size_t nbr_trans_couples = couple->num_field_couples;
    for(size_t trans_couple_idx = 0; trans_couple_idx < nbr_trans_couples;
        ++trans_couple_idx){
      struct yac_couple_config_field_couple * transcouple =
        couple->field_couples + trans_couple_idx;
      if(transcouple->source.component_idx == component_idx &&
        transcouple->source.field_idx == field_idx)
        return YAC_EXCHANGE_TYPE_SOURCE;
      if(transcouple->target.component_idx == component_idx &&
        transcouple->target.field_idx == field_idx)
        return YAC_EXCHANGE_TYPE_TARGET;
    }
  }
  return YAC_EXCHANGE_TYPE_NONE;
}

int yac_couple_config_field_is_valid(
  struct yac_couple_config * couple_config,
  size_t component_idx, size_t field_idx) {

  check_field_idx(
    couple_config, component_idx, field_idx,
      "yac_couple_config_field_is_valid", __LINE__);

  struct yac_couple_config_field * field =
    couple_config->components[component_idx].fields +
    field_idx;

  return (field->collection_size != SIZE_MAX) &&
         (field->timestep != NULL);
}

static void check_field_couple_idx(
  struct yac_couple_config * couple_config, size_t couple_idx,
  size_t field_couple_idx, char const * routine_name, int line) {

  check_couple_idx(couple_config, couple_idx, routine_name, __LINE__);

  YAC_ASSERT_F(
    field_couple_idx <
      couple_config->couples[couple_idx].num_field_couples,
    "ERROR(%s:%d:%s): invalid field_couple_idx",
    __FILE__, line, routine_name)
}

struct yac_interp_stack_config * yac_couple_config_get_interp_stack(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_interp_stack", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].interp_stack;
}

void yac_couple_config_get_field_grid_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const ** src_grid_name, char const ** tgt_grid_name) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_field_grid_names", __LINE__);

  size_t src_component_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].
        source.component_idx;
  size_t src_field_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].
        source.field_idx;

  size_t tgt_component_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].
        target.component_idx;
  size_t tgt_field_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].
        target.field_idx;

  *src_grid_name =
    couple_config->grids[
      couple_config->components[src_component_idx].
        fields[src_field_idx].grid_idx].name;
  *tgt_grid_name =
    couple_config->grids[
      couple_config->components[tgt_component_idx].
        fields[tgt_field_idx].grid_idx].name;
}

void yac_couple_config_field_enable_frac_mask(
  struct yac_couple_config * couple_config,
  char const * comp_name, char const * grid_name, char const * field_name,
  double frac_mask_fallback_value) {

  YAC_ASSERT_F(
    (frac_mask_fallback_value != YAC_FRAC_MASK_UNDEF) &&
    (frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE),
    "ERROR(yac_couple_config_field_enable_frac_mask): "
    "\"%lf\" is not a valid fractional mask fallback value "
    "(component: \"%s\" grid: \"%s\" field \"%s\")",
    frac_mask_fallback_value, comp_name, grid_name, field_name);

  size_t comp_idx = yac_couple_config_get_component_idx(couple_config, comp_name);
  size_t grid_idx = yac_couple_config_get_grid_idx(couple_config, grid_name);
  size_t field_idx = yac_couple_config_get_field_idx(couple_config,
    comp_idx, grid_idx, field_name);

  double old_frac_mask_fallback_value =
    couple_config->components[comp_idx].fields[field_idx].
      frac_mask_fallback_value;

  YAC_ASSERT_F(
    (old_frac_mask_fallback_value == YAC_FRAC_MASK_UNDEF) ||
    (old_frac_mask_fallback_value == YAC_FRAC_MASK_NO_VALUE) ||
    (old_frac_mask_fallback_value == frac_mask_fallback_value),
    "ERROR(yac_couple_config_field_enable_frac_mask): "
    "fractional mask fallback value was already set:\n"
    "\told value \"%lf\" new value \"%lf\" "
    "(component: \"%s\" grid: \"%s\" field \"%s\")",
    old_frac_mask_fallback_value, frac_mask_fallback_value,
    comp_name, grid_name, field_name);

  couple_config->components[comp_idx].fields[field_idx].
      frac_mask_fallback_value = frac_mask_fallback_value;
}

double yac_couple_config_get_frac_mask_fallback_value(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name) {

  size_t component_idx =
    yac_couple_config_get_component_idx(couple_config, component_name);
  size_t grid_idx =
    yac_couple_config_get_grid_idx(couple_config, grid_name);
  size_t field_idx =
    yac_couple_config_get_field_idx(
      couple_config, component_idx, grid_idx, field_name);

  struct yac_couple_config_field * field =
    couple_config->components[component_idx].fields +
    field_idx;

  YAC_ASSERT_F(
    field->frac_mask_fallback_value != YAC_FRAC_MASK_UNDEF,
    "ERROR(yac_couple_config_get_frac_mask_fallback_value): "
    "no valid fractional mask fallback value defined "
    "(component: \"%s\" grid: \"%s\" field \"%s\")",
    couple_config->components[component_idx].name, grid_name, field->name);

  return field->frac_mask_fallback_value;
}

size_t yac_couple_config_get_collection_size(
  struct yac_couple_config * couple_config,
  char const * component_name, char const * grid_name,
  char const * field_name) {

  size_t component_idx =
    yac_couple_config_get_component_idx(couple_config, component_name);
  size_t grid_idx =
    yac_couple_config_get_grid_idx(couple_config, grid_name);
  size_t field_idx =
    yac_couple_config_get_field_idx(
      couple_config, component_idx, grid_idx, field_name);

  struct yac_couple_config_field * field =
    couple_config->components[component_idx].fields +
    field_idx;

  YAC_ASSERT_F(
    field->collection_size != SIZE_MAX,
    "ERROR(yac_couple_config_get_collection_size): "
    "no valid collection size defined "
    "(component: \"%s\" grid: \"%s\" field \"%s\")",
    couple_config->components[component_idx].name, grid_name, field->name);

  return field->collection_size;
}

void yac_couple_config_get_field_couple_component_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const ** src_component_name, char const ** tgt_component_name) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_field_couple_component_names", __LINE__);

  *src_component_name =
    couple_config->components[
      couple_config->couples[couple_idx].
        field_couples[field_couple_idx].source.component_idx].name;
  *tgt_component_name =
    couple_config->components[
      couple_config->couples[couple_idx].
        field_couples[field_couple_idx].target.component_idx].name;
}

void yac_couple_config_get_field_names(
  struct yac_couple_config * couple_config,
    size_t couple_idx, size_t field_couple_idx,
    char const ** src_field_name, const char ** tgt_field_name) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_field_names", __LINE__);

  size_t src_component_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].source.component_idx;
  size_t tgt_component_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].target.component_idx;
  size_t src_field_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].source.field_idx;
  size_t tgt_field_idx =
    couple_config->couples[couple_idx].
      field_couples[field_couple_idx].target.field_idx;

  *src_field_name =
    couple_config->components[src_component_idx].
    fields[src_field_idx].name;
  *tgt_field_name =
      couple_config->components[tgt_component_idx].
        fields[tgt_field_idx].name;
}

int yac_couple_config_mapping_on_source(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_mapping_on_source", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          mapping_on_source;
}

int yac_couple_config_get_source_lag(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_source_lag", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          source.lag;
}

int yac_couple_config_get_target_lag(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_target_lag", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          target.lag;
}

char const * yac_couple_config_get_coupling_period(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_coupling_period", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          coupling_period;
}

char const * yac_couple_config_get_source_timestep(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_source_timestep", __LINE__);

  struct yac_couple_config_field_couple * tcouple =
    couple_config->couples[couple_idx].field_couples + field_couple_idx;
  char const * timestep =
    couple_config->components[tcouple->source.component_idx].
      fields[tcouple->source.field_idx].timestep;

  YAC_ASSERT_F(
    timestep,
    "ERROR(yac_couple_config_get_source_timestep): "
    "no valid timestep defined (component: \"%s\" field \"%s\")",
    couple_config->components[tcouple->source.component_idx].name,
    couple_config->components[tcouple->source.component_idx].
      fields[tcouple->source.field_idx].name);

  return timestep;
}

char const * yac_couple_config_get_target_timestep(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_target_timestep", __LINE__);

  struct yac_couple_config_field_couple * tcouple =
    couple_config->couples[couple_idx].field_couples + field_couple_idx;
  char const * timestep =
    couple_config->components[tcouple->target.component_idx].
      fields[tcouple->target.field_idx].timestep;

  YAC_ASSERT_F(
    timestep,
    "ERROR(yac_couple_config_get_target_timestep): "
    "no valid timestep defined (component: \"%s\" field \"%s\")",
    couple_config->components[tcouple->target.component_idx].name,
    couple_config->components[tcouple->target.component_idx].
      fields[tcouple->target.field_idx].name);

  return timestep;
}

enum yac_reduction_type yac_couple_config_get_coupling_period_operation(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_coupling_period_operation", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          coupling_period_operation;
}

void yac_couple_config_set_datetime(
  struct yac_couple_config * couple_config,
  char const * start, char const * end) {

  if ((start != NULL) && (strlen(start) > 0)) {

    // in case start points to couple_config->start_datetime,
    // we have to first strdup start, before freeing it
    char const * old_start = couple_config->start_datetime;
    couple_config->start_datetime = strdup(start);
    free((void*)old_start);
  }

  if ((end != NULL) && (strlen(end) > 0)) {

    // in case end points to couple_config->end_datetime,
    // we have to first strdup end, before freeing it
    char const * old_end = couple_config->end_datetime;
    couple_config->end_datetime = strdup(end);
    free((void*)old_end);
  }
}

char const * yac_couple_config_get_start_datetime(
  struct yac_couple_config * couple_config) {

  YAC_ASSERT(couple_config->start_datetime,
    "ERROR(yac_couple_config_get_start_datetime): "
    "start_datetime not yet defined");

  return couple_config->start_datetime;
}

char const * yac_couple_config_get_end_datetime(
  struct yac_couple_config * couple_config) {

  YAC_ASSERT(couple_config->end_datetime,
    "ERROR(yac_couple_config_get_start_datetime): "
    "start_datetime not yet defined");

  return couple_config->end_datetime;
}

char const * yac_couple_config_get_grid_name(
  struct yac_couple_config * couple_config, size_t grid_idx) {

  YAC_ASSERT_F(
    grid_idx < couple_config->num_grids,
      "ERROR(yac_couple_config_get_grid_name): "
      "Invalid grid idx %zu", grid_idx);

  return couple_config->grids[grid_idx].name;
}

int yac_couple_config_enforce_write_weight_file(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_enforce_write_weight_file", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          enforce_write_weight_file;
}

char const * yac_couple_config_get_weight_file_name(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_weight_file_name", __LINE__);

  static char dummy[] = "\0";
  char const * weight_file_name =
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          weight_file_name;

  return (weight_file_name != NULL)?weight_file_name:dummy;
}

double yac_couple_config_get_scale_factor(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_scale_factor", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          scale_factor;
}

double yac_couple_config_get_scale_summand(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_scale_summand", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          scale_summand;
}

void yac_couple_config_get_src_mask_names(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  char const * const ** mask_names, size_t * num_mask_names) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_src_mask_names", __LINE__);

  *mask_names =
    (char const * const *)(
      couple_config->
        couples[couple_idx].
          field_couples[field_couple_idx].
            src_mask_names);
  *num_mask_names =
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          num_src_mask_names;
}

char const * yac_couple_config_get_tgt_mask_name(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  check_field_couple_idx(
    couple_config, couple_idx, field_couple_idx,
    "yac_couple_config_get_tgt_mask_name", __LINE__);

  return
    couple_config->
      couples[couple_idx].
        field_couples[field_couple_idx].
          tgt_mask_name;
}

int yac_couple_config_get_redirect_stdout(
  struct yac_couple_config * couple_config) {

  return couple_config->redirect_stdout;
}

void yac_couple_config_set_redirect_stdout(
  struct yac_couple_config * couple_config, int redirect_stdout) {

  YAC_ASSERT(
    (redirect_stdout >= 0) && (redirect_stdout <= 2),
    "ERROR(yac_couple_config_set_redirect_stdout): "
    "invalid redirect_stdout out value (has to be 0 <= x <= 2)");

  couple_config->redirect_stdout = redirect_stdout;
}

int yac_couple_config_contains_grid_name(
  struct yac_couple_config * couple_config, char const * grid_name) {

  for (size_t grid_idx = 0; grid_idx < couple_config->num_grids;
       ++grid_idx)
    if (!strcmp(couple_config->grids[grid_idx].name, grid_name))
      return 1;
  return 0;
}

static size_t yac_couple_config_get_string_pack_size(
  char const * string, MPI_Comm comm) {

  int strlen_pack_size, string_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &strlen_pack_size), comm);

  if (string != NULL) {
    yac_mpi_call(
      MPI_Pack_size(
        (int)(strlen(string)), MPI_CHAR, comm, &string_pack_size), comm);
  } else {
    string_pack_size = 0;
  }

  return (size_t)strlen_pack_size + (size_t)string_pack_size;
}

static size_t yac_couple_config_get_grid_pack_size(
  void * grid_, MPI_Comm comm) {

  struct yac_couple_config_grid * grid = grid_;

  return yac_couple_config_get_string_pack_size(grid->name, comm);
}

static size_t yac_couple_config_get_grids_pack_size(
  size_t num_grids, void * grids_, MPI_Comm comm) {

  struct yac_couple_config_grid * grids = grids_;

  int num_grids_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &num_grids_pack_size), comm);

  size_t grids_pack_size = 0;
  for (size_t i = 0; i < num_grids; ++i)
    grids_pack_size +=
      yac_couple_config_get_grid_pack_size(grids + i, comm);

  return (size_t)num_grids_pack_size + grids_pack_size;
}

static size_t yac_couple_config_get_field_pack_size(
  struct yac_couple_config_field * field,
  MPI_Comm comm) {

  YAC_ASSERT(
    field->grid_idx <= INT_MAX,
    "ERROR(yac_couple_config_get_field_pack_size):"
    "grid_idx is too big")

  int int_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);
  int double_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_DOUBLE, comm, &double_pack_size), comm);
  size_t name_pack_size = yac_couple_config_get_string_pack_size(
    field->name, comm);
  size_t timestep_pack_size = yac_couple_config_get_string_pack_size(
    field->timestep, comm);
  size_t metadata_pack_size = yac_couple_config_get_string_pack_size(
    field->metadata, comm);

  return int_pack_size + // grid_idx
         int_pack_size + // collection_size
         double_pack_size + // frac_mask_fallback_value
         name_pack_size +
         timestep_pack_size +
         metadata_pack_size;
}

static size_t yac_couple_config_get_fields_pack_size(
  size_t num_fields, void * fields_, MPI_Comm comm) {

  struct yac_couple_config_field * fields = fields_;

  int num_fields_pack_size;
  yac_mpi_call(
    MPI_Pack_size(1, MPI_INT, comm, &num_fields_pack_size), comm);

  size_t fields_pack_size = 0;
  for (size_t i = 0; i < num_fields; ++i)
    fields_pack_size +=
      yac_couple_config_get_field_pack_size(
        fields + i, comm);

  return (size_t)num_fields_pack_size +
         fields_pack_size;
}

static size_t yac_couple_config_get_component_pack_size(
  struct yac_couple_config_component * component, MPI_Comm comm) {

  return
    yac_couple_config_get_string_pack_size(component->name, comm);
}

static size_t yac_couple_config_get_components_pack_size(
  size_t num_components, void * components_, MPI_Comm comm) {

  struct yac_couple_config_component * components = components_;

  int num_components_pack_size;
  yac_mpi_call(
    MPI_Pack_size(1, MPI_INT, comm, &num_components_pack_size), comm);

  size_t components_pack_size = 0;
  for (size_t i = 0; i < num_components; ++i)
    components_pack_size +=
      yac_couple_config_get_component_pack_size(components + i, comm);

  return (size_t)num_components_pack_size + components_pack_size;
}

static size_t yac_couple_config_get_field_couple_pack_size(
  struct yac_couple_config_field_couple * field_couple, MPI_Comm comm) {

  int ints_pack_size;
  yac_mpi_call(MPI_Pack_size(10, MPI_INT, comm, &ints_pack_size), comm);
  int doubles_pack_size;
  yac_mpi_call(MPI_Pack_size(2, MPI_DOUBLE, comm, &doubles_pack_size), comm);
  int src_mask_names_pack_size = 0;
  if (field_couple->num_src_mask_names > 0)
    for (size_t i = 0; i < field_couple->num_src_mask_names; ++i)
      src_mask_names_pack_size +=
        yac_couple_config_get_string_pack_size(
          field_couple->src_mask_names[i], comm);

  return
    (size_t)ints_pack_size + // source.component_idx
                             // source.field_idx
                             // source.lag
                             // target.component_idx
                             // target.field_idx
                             // target.lag
                             // mapping_on_source
                             // coupling_period_operation
                             // enforce_write_weight_file
                             // num_src_mask_names
    yac_interp_stack_config_get_pack_size(
      field_couple->interp_stack, comm) +
    yac_couple_config_get_string_pack_size(
      field_couple->coupling_period, comm) +
    yac_couple_config_get_string_pack_size(
      field_couple->weight_file_name, comm) +
    doubles_pack_size + // scale_factor
                        // scale_summand
    src_mask_names_pack_size +
    yac_couple_config_get_string_pack_size(
      field_couple->tgt_mask_name, comm);
}

static size_t yac_couple_config_get_field_couples_pack_size(
  size_t num_field_couples, void * field_couples_, MPI_Comm comm) {

  struct yac_couple_config_field_couple * field_couples = field_couples_;

  int num_field_couples_pack_size;
  yac_mpi_call(
    MPI_Pack_size(1, MPI_INT, comm, &num_field_couples_pack_size), comm);

  size_t field_couples_pack_size = 0;
  for (size_t i = 0; i < num_field_couples; ++i)
    field_couples_pack_size +=
      yac_couple_config_get_field_couple_pack_size(field_couples + i, comm);

  return (size_t)num_field_couples_pack_size +
         field_couples_pack_size;
}

static size_t yac_couple_config_get_couple_pack_size_basic(
  struct yac_couple_config_couple * couple, MPI_Comm comm) {

  int component_indices_pack_size;
  yac_mpi_call(
    MPI_Pack_size(2, MPI_INT, comm, &component_indices_pack_size), comm);

  return (size_t)component_indices_pack_size;
}

static size_t yac_couple_config_get_couples_pack_size_basic(
  size_t num_couples, void * couples_, MPI_Comm comm) {

  struct yac_couple_config_couple * couples = couples_;

  int num_couples_pack_size;
  yac_mpi_call(
    MPI_Pack_size(1, MPI_INT, comm, &num_couples_pack_size), comm);

  size_t couples_pack_size = 0;
  for (size_t i = 0; i < num_couples; ++i)
    couples_pack_size +=
      yac_couple_config_get_couple_pack_size_basic(couples + i, comm);

  return (size_t)num_couples_pack_size + couples_pack_size;
}

static void yac_couple_config_pack_string(
  char const * string, void * buffer, int buffer_size, int * position,
  MPI_Comm comm) {

  size_t len = (string == NULL)?0:strlen(string);

  YAC_ASSERT(
    len <= INT_MAX, "ERROR(yac_couple_config_pack_string): string too long")

  int len_int = (int)len;

  yac_mpi_call(
    MPI_Pack(
      &len_int, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  if (len > 0)
    yac_mpi_call(
      MPI_Pack(
        string, len_int, MPI_CHAR, buffer, buffer_size, position, comm),
      comm);
}

static void yac_couple_config_pack_grid(
  struct yac_couple_config_grid * grid,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  yac_couple_config_pack_string(
    grid->name, buffer, buffer_size, position, comm);
}

static void yac_couple_config_pack_grids(
  size_t num_grids, void * grids_,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  struct yac_couple_config_grid * grids = grids_;

  YAC_ASSERT(
    num_grids <= INT_MAX,
    "ERROR(yac_couple_config_pack_grids): num_grids bigger than INT_MAX")

  int num_grids_int = (int)num_grids;
  yac_mpi_call(
    MPI_Pack(
      &num_grids_int, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  for (size_t i = 0; i < num_grids; ++i)
    yac_couple_config_pack_grid(
      grids + i, buffer, buffer_size, position, comm);
}

static void yac_couple_config_pack_field(
  struct yac_couple_config_field * field,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  YAC_ASSERT(
    field->grid_idx <= INT_MAX,
    "ERROR(yac_couple_config_pack_field): grid_idx is too big")
  YAC_ASSERT(
    (field->collection_size < INT_MAX) ||
    (field->collection_size == SIZE_MAX),
    "ERROR(yac_couple_config_pack_field): invalid collection size")

  int grid_idx = (int)field->grid_idx;
  double frac_mask_fallback_value = field->frac_mask_fallback_value;
  int collection_size =
    (field->collection_size == SIZE_MAX)?
      INT_MAX:(int)(field->collection_size);

  yac_mpi_call(
    MPI_Pack(
      &grid_idx, 1, MPI_INT, buffer, buffer_size, position, comm), comm);
  yac_mpi_call(
    MPI_Pack(
      &frac_mask_fallback_value, 1, MPI_DOUBLE, buffer, buffer_size,
      position, comm), comm);
  yac_mpi_call(
    MPI_Pack(
      &collection_size, 1, MPI_INT, buffer, buffer_size, position, comm), comm);
  yac_couple_config_pack_string(field->name, buffer,
    buffer_size, position, comm);
  yac_couple_config_pack_string(field->timestep, buffer,
    buffer_size, position, comm);
  yac_couple_config_pack_string(field->metadata, buffer,
    buffer_size, position, comm);
}

static void yac_couple_config_pack_fields(
  size_t num_fields, void * fields_, void * buffer, int buffer_size,
  int * position, MPI_Comm comm) {

  struct yac_couple_config_field * fields = fields_;

  YAC_ASSERT(
    num_fields <= INT_MAX,
    "ERROR(yac_couple_config_pack_fields): "
    "num_fields bigger than INT_MAX")

  int num_fields_int = (int)num_fields;
  yac_mpi_call(
    MPI_Pack(
      &num_fields_int, 1, MPI_INT,
      buffer, buffer_size, position, comm), comm);

  for (size_t i = 0; i < num_fields; ++i)
    yac_couple_config_pack_field(
      fields + i, buffer, buffer_size, position, comm);
}

static void yac_couple_config_pack_component(
  struct yac_couple_config_component * component,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  yac_couple_config_pack_string(
    component->name, buffer, buffer_size, position, comm);
}

static void yac_couple_config_pack_components(
  size_t num_components, void * components_,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  struct yac_couple_config_component * components = components_;

  YAC_ASSERT(
    num_components <= INT_MAX,
    "ERROR(yac_couple_config_pack_components): "
    "num_components bigger than INT_MAX")

  int num_components_int = (int)num_components;
  yac_mpi_call(
    MPI_Pack(
      &num_components_int, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  for (size_t i = 0; i < num_components; ++i)
    yac_couple_config_pack_component(
      components + i, buffer, buffer_size, position, comm);
}

static void yac_couple_config_pack_field_couple(
  struct yac_couple_config_field_couple * field_couple,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  YAC_ASSERT(
    field_couple->source.component_idx <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "source.component_idx bigger than INT_MAX")
  YAC_ASSERT(
    field_couple->source.field_idx <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "source.field_idx bigger than INT_MAX")
  YAC_ASSERT(
    field_couple->target.component_idx <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "target.component_idx bigger than INT_MAX")
  YAC_ASSERT(
    field_couple->target.field_idx <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "target.field_idx bigger than INT_MAX")
  YAC_ASSERT(
    field_couple->mapping_on_source <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "mapping_on_source bigger than INT_MAX")
  YAC_ASSERT(
    field_couple->coupling_period_operation <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "coupling_period_operation bigger than INT_MAX")
  YAC_ASSERT(
    field_couple->enforce_write_weight_file <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "enforce_write_weight_file bigger than INT_MAX")
  YAC_ASSERT(
    field_couple->num_src_mask_names <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couple): "
    "num_src_mask_names bigger than INT_MAX")

  int ints[10] = {
    field_couple->source.component_idx,
    field_couple->source.field_idx,
    field_couple->source.lag,
    field_couple->target.component_idx,
    field_couple->target.field_idx,
    field_couple->target.lag,
    field_couple->mapping_on_source,
    field_couple->coupling_period_operation,
    field_couple->enforce_write_weight_file,
    field_couple->num_src_mask_names};

  yac_mpi_call(
    MPI_Pack(ints, 10, MPI_INT, buffer, buffer_size, position, comm), comm);

  yac_couple_config_pack_string(
    field_couple->coupling_period, buffer, buffer_size, position, comm);
  yac_couple_config_pack_string(
    field_couple->weight_file_name, buffer, buffer_size, position, comm);

  double doubles[2] = {
    field_couple->scale_factor,
    field_couple->scale_summand};

  yac_mpi_call(
    MPI_Pack(doubles, 2, MPI_DOUBLE, buffer, buffer_size, position, comm),
    comm);

  yac_interp_stack_config_pack(
    field_couple->interp_stack, buffer, buffer_size, position, comm);

  if (field_couple->num_src_mask_names > 0)
    for (size_t i = 0; i < field_couple->num_src_mask_names; ++i)
      yac_couple_config_pack_string(
        field_couple->src_mask_names[i], buffer, buffer_size, position, comm);

  yac_couple_config_pack_string(
    field_couple->tgt_mask_name, buffer, buffer_size, position, comm);
}

static void yac_couple_config_pack_field_couples(
  size_t num_field_couples, void * field_couples_,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  struct yac_couple_config_field_couple * field_couples = field_couples_;

  YAC_ASSERT(
    num_field_couples <= INT_MAX,
    "ERROR(yac_couple_config_pack_field_couples): "
    "num_field_couples bigger than INT_MAX")

  int num_field_couples_int = (int)num_field_couples;
  yac_mpi_call(
    MPI_Pack(
      &num_field_couples_int, 1, MPI_INT,
      buffer, buffer_size, position, comm), comm);

  for (size_t i = 0; i < num_field_couples; ++i)
    yac_couple_config_pack_field_couple(
      field_couples + i, buffer, buffer_size, position, comm);
}

static void yac_couple_config_pack_couple_basic(
  struct yac_couple_config_couple * couple,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  YAC_ASSERT(
    (couple->component_indices[0] <= INT_MAX) &&
    (couple->component_indices[1] <= INT_MAX),
    "ERROR(yac_couple_config_pack_couple_basic): "
    "component_indices bigger than INT_MAX")

  int component_indices[2] = {(int)(couple->component_indices[0]),
                              (int)(couple->component_indices[1])};

  yac_mpi_call(
    MPI_Pack(
      component_indices, 2, MPI_INT, buffer, buffer_size, position, comm),
    comm);
}

static void yac_couple_config_pack_couples_basic(
  size_t num_couples, void * couples_,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  struct yac_couple_config_couple * couples = couples_;

  YAC_ASSERT(
    num_couples <= INT_MAX,
    "ERROR(yac_couple_config_pack_couples_basic): "
    "num_couples bigger than INT_MAX")

  int num_couples_int = (int)num_couples;
  yac_mpi_call(
    MPI_Pack(
      &num_couples_int, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  for (size_t i = 0; i < num_couples; ++i)
    yac_couple_config_pack_couple_basic(
      couples + i, buffer, buffer_size, position, comm);
}

static char * yac_couple_config_unpack_string(
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int string_len;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &string_len, 1, MPI_INT, comm), comm);

  if (string_len <= 0) return NULL;

  char * string = xmalloc(((size_t)string_len + 1) * sizeof(*string));
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, string, string_len, MPI_CHAR, comm), comm);
  string[string_len] = '\0';
  return string;
}

static void yac_couple_config_unpack_grid(
  void * buffer, int buffer_size, int * position,
  struct yac_couple_config_grid * grid, MPI_Comm comm) {

  grid->name =
    yac_couple_config_unpack_string(buffer, buffer_size, position, comm);
  grid->metadata = NULL;
}

static void yac_couple_config_unpack_grids(
  void * buffer, int buffer_size, int * position,
  size_t * num_grids, void * grids_, MPI_Comm comm) {

  struct yac_couple_config_grid ** grids = grids_;

  int num_grids_int;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &num_grids_int, 1, MPI_INT, comm), comm);

  YAC_ASSERT(
    num_grids_int >= 0,
    "ERROR(yac_couple_config_unpack_grids): invalid number of grids")

  *num_grids = (size_t)num_grids_int;
  *grids = xmalloc(*num_grids * sizeof(**grids));

  for (size_t i = 0; i < *num_grids; ++i)
    yac_couple_config_unpack_grid(
      buffer, buffer_size, position, *grids + i, comm);
}

static void yac_couple_config_unpack_field(
  void * buffer, int buffer_size, int * position,
  struct yac_couple_config_field * field, MPI_Comm comm) {

  int grid_idx;
  double frac_mask_fallback_value;
  int collection_size;

  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &grid_idx, 1, MPI_INT, comm), comm);
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &frac_mask_fallback_value, 1,
      MPI_DOUBLE, comm), comm);
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &collection_size, 1, MPI_INT, comm),
    comm);

  YAC_ASSERT(
    grid_idx >= 0,
    "ERROR(yac_couple_config_unpack_field): invalid number of grid_idx")
  YAC_ASSERT(
    collection_size >= 0,
    "ERROR(yac_couple_config_unpack_field): invalid collection_size")

  field->grid_idx = (size_t)grid_idx;
  field->frac_mask_fallback_value = frac_mask_fallback_value;
  field->collection_size =
    (collection_size == INT_MAX)?SIZE_MAX:((int)collection_size);
  field->name = yac_couple_config_unpack_string(buffer,
    buffer_size, position, comm);
  field->timestep = yac_couple_config_unpack_string(buffer,
    buffer_size, position, comm);
  field->metadata = yac_couple_config_unpack_string(buffer,
    buffer_size, position, comm);
}

static void yac_couple_config_unpack_fields(
  void * buffer, int buffer_size, int * position,
  size_t * num_fields, void * fields_, MPI_Comm comm) {

  struct yac_couple_config_field ** fields = fields_;

  int num_fields_int;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &num_fields_int, 1,
      MPI_INT, comm), comm);

  YAC_ASSERT(
    num_fields_int >= 0,
    "ERROR(yac_couple_config_unpack_fields): "
    "invalid number of fields")

  *num_fields = (size_t)num_fields_int;
  *fields =
    xmalloc(*num_fields * sizeof(**fields));

  for (size_t i = 0; i < *num_fields; ++i)
    yac_couple_config_unpack_field(
      buffer, buffer_size, position, *fields + i, comm);
}

static void yac_couple_config_unpack_component(
  void * buffer, int buffer_size, int * position,
  struct yac_couple_config_component * component, MPI_Comm comm) {

  component->name =
    yac_couple_config_unpack_string(buffer, buffer_size, position, comm);
  component->metadata = NULL;
  component->num_fields = 0;
  component->fields = NULL;
}

static void yac_couple_config_unpack_components(
  void * buffer, int buffer_size, int * position,
  size_t * num_components, void * components_, MPI_Comm comm) {

  struct yac_couple_config_component ** components = components_;

  int num_components_int;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &num_components_int, 1, MPI_INT, comm),
    comm);

  YAC_ASSERT(
    num_components_int >= 0,
    "ERROR(yac_couple_config_unpack_components): "
    "invalid number of components")

  *num_components = (size_t)num_components_int;
  *components = xmalloc(*num_components * sizeof(**components));

  for (size_t i = 0; i < *num_components; ++i)
    yac_couple_config_unpack_component(
      buffer, buffer_size, position, *components + i, comm);
}

static void yac_couple_config_unpack_field_couple(
  void * buffer, int buffer_size, int * position,
  struct yac_couple_config_field_couple * field_couple,
  MPI_Comm comm) {

  int ints[10];
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, ints, 10, MPI_INT, comm), comm);

  YAC_ASSERT(
    ints[0] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "invalid source.component_idx")
  YAC_ASSERT(
    ints[1] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "invalid source.field_idx")
  YAC_ASSERT(
    ints[3] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "target.component_idx bigger than INT_MAX")
  YAC_ASSERT(
    ints[4] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "invalid target.field_idx")
  YAC_ASSERT(
    ints[6] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "invalid mapping_on_source")
  YAC_ASSERT(
    ints[7] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "invalid coupling_period_operation")
  YAC_ASSERT(
    ints[8] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "invalid enforce_write_weight_file")
  YAC_ASSERT(
    ints[9] >= 0,
    "ERROR(yac_couple_config_unpack_field_couple): "
    "invalid num_src_mask_names")

  field_couple->source.component_idx = (size_t)(ints[0]);
  field_couple->source.field_idx = (size_t)(ints[1]);
  field_couple->source.lag = ints[2];
  field_couple->target.component_idx = (size_t)(ints[3]);
  field_couple->target.field_idx = (size_t)(ints[4]);
  field_couple->target.lag = ints[5];
  field_couple->mapping_on_source = ints[6];
  field_couple->coupling_period_operation =
    (enum yac_reduction_type)(ints[7]);
  field_couple->enforce_write_weight_file = ints[8];
  field_couple->num_src_mask_names = (size_t)(ints[9]);

  field_couple->coupling_period =
    yac_couple_config_unpack_string(buffer, buffer_size, position, comm);
  field_couple->weight_file_name =
    yac_couple_config_unpack_string(buffer, buffer_size, position, comm);

  double doubles[2];
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, doubles, 2, MPI_DOUBLE, comm), comm);

  field_couple->scale_factor = doubles[0];
  field_couple->scale_summand = doubles[1];

  field_couple->interp_stack =
    yac_interp_stack_config_unpack(buffer, buffer_size, position, comm);

  if (field_couple->num_src_mask_names > 0) {
    field_couple->src_mask_names =
      xmalloc(
        field_couple->num_src_mask_names *
        sizeof(*(field_couple->src_mask_names)));
    for (size_t i = 0; i < field_couple->num_src_mask_names; ++i)
      field_couple->src_mask_names[i] =
        yac_couple_config_unpack_string(
          buffer, buffer_size, position, comm);
  } else {
    field_couple->src_mask_names = NULL;
  }

  field_couple->tgt_mask_name =
    yac_couple_config_unpack_string(
      buffer, buffer_size, position, comm);
}

static void yac_couple_config_unpack_field_couples(
  void * buffer, int buffer_size, int * position,
  size_t * num_field_couples, void * field_couples_, MPI_Comm comm) {

  struct yac_couple_config_field_couple ** field_couples = field_couples_;

  int num_field_couples_int;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &num_field_couples_int, 1,
      MPI_INT, comm), comm);

  YAC_ASSERT(
    num_field_couples_int >= 0,
    "ERROR(yac_couple_config_unpack_field_couples): "
    "invalid number of field_couples")

  *num_field_couples = (size_t)num_field_couples_int;
  *field_couples =
    xmalloc(*num_field_couples * sizeof(**field_couples));

  for (size_t i = 0; i < *num_field_couples; ++i)
    yac_couple_config_unpack_field_couple(
      buffer, buffer_size, position, *field_couples + i, comm);
}

static void yac_couple_config_unpack_couple_basic(
  void * buffer, int buffer_size, int * position,
  struct yac_couple_config_couple * couple, MPI_Comm comm) {

  int component_indices[2];
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, component_indices, 2, MPI_INT, comm),
    comm);

  YAC_ASSERT(
    (component_indices[0] >= 0) && (component_indices[1] >= 0),
    "ERROR(yac_couple_config_unpack_couple_basic): invalid component indices")

  couple->component_indices[0] = (size_t)(component_indices[0]);
  couple->component_indices[1] = (size_t)(component_indices[1]);
  couple->num_field_couples = 0;
  couple->field_couples = NULL;
}

static void yac_couple_config_unpack_couples_basic(
  void * buffer, int buffer_size, int * position,
  size_t * num_couples, void * couples_, MPI_Comm comm) {

  struct yac_couple_config_couple ** couples = couples_;

  int num_couples_int;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &num_couples_int, 1, MPI_INT, comm),
    comm);

  YAC_ASSERT(
    num_couples_int >= 0,
    "ERROR(yac_couple_config_unpack_couples_basic): "
    "invalid number of couples")

  *num_couples = (size_t)num_couples_int;
  *couples = xmalloc(*num_couples * sizeof(**couples));

  for (size_t i = 0; i < *num_couples; ++i)
    yac_couple_config_unpack_couple_basic(
      buffer, buffer_size, position, *couples + i, comm);
}

void yac_couple_config_def_couple(
  struct yac_couple_config * couple_config,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_period, int time_reduction,
  struct yac_interp_stack_config * interp_stack,
  int src_lag, int tgt_lag,
  const char* weight_file_name, int mapping_on_source,
  double scale_factor, double scale_summand,
  size_t num_src_mask_names, char const * const * src_mask_names,
  char const * tgt_mask_name) {

  YAC_ASSERT(src_comp_name && src_comp_name[0] != '\0',
    "ERROR(yac_couple_config_def_couple): invalid parameter: src_comp_name");
  YAC_ASSERT(src_grid_name && src_grid_name[0] != '\0',
    "ERROR(yac_couple_config_def_couple): invalid parameter: src_grid_name");
  YAC_ASSERT(src_field_name && src_field_name[0] != '\0',
    "ERROR(yac_couple_config_def_couple): invalid parameter: src_field_name");
  YAC_ASSERT(tgt_comp_name && tgt_comp_name[0] != '\0',
    "ERROR(yac_couple_config_def_couple): invalid parameter: tgt_comp_name");
  YAC_ASSERT(tgt_grid_name && tgt_grid_name[0] != '\0',
    "ERROR(yac_couple_config_def_couple): invalid parameter: tgt_grid_name");
  YAC_ASSERT(tgt_field_name && tgt_field_name[0] != '\0',
    "ERROR(yac_couple_config_def_couple): invalid parameter: tgt_field_name");
  YAC_ASSERT(coupling_period && coupling_period[0] != '\0',
    "ERROR(yac_couple_config_def_couple): invalid parameter: coupling_period");
  YAC_ASSERT(
    (time_reduction == TIME_NONE) ||
    (time_reduction == TIME_ACCUMULATE) ||
    (time_reduction == TIME_AVERAGE) ||
    (time_reduction == TIME_MINIMUM) ||
    (time_reduction == TIME_MAXIMUM),
    "ERROR(yac_couple_config_def_couple): invalid parameter: time_reduction");
  YAC_ASSERT_F(isnormal(scale_factor),
    "ERROR(yac_couple_config_def_couple): \"%lf\" is not a valid scale factor",
    scale_factor);
  YAC_ASSERT_F(isnormal(scale_summand) || (scale_summand == 0.0),
    "ERROR(yac_couple_config_def_couple): \"%lf\" is not a valid scale summand",
    scale_summand);

  // get component indices
  size_t src_comp_idx =
    yac_couple_config_add_component_(couple_config, src_comp_name);
  size_t tgt_comp_idx =
    yac_couple_config_add_component_(couple_config, tgt_comp_name);
  size_t src_grid_idx =
    yac_couple_config_add_grid_(couple_config, src_grid_name);
  size_t tgt_grid_idx =
    yac_couple_config_add_grid_(couple_config, tgt_grid_name);

  // check if couple exists
  size_t component_indices[2];
  if(src_comp_idx < tgt_comp_idx){
    component_indices[0] = src_comp_idx;
    component_indices[1] = tgt_comp_idx;
  }else{
    component_indices[0] = tgt_comp_idx;
    component_indices[1] = src_comp_idx;
  }
  struct yac_couple_config_couple * couple = NULL;
  for(size_t i = 0; (i < couple_config->num_couples) && !couple; ++i)
    if(couple_config->couples[i].component_indices[0] == component_indices[0] &&
       couple_config->couples[i].component_indices[1] == component_indices[1])
      couple = &couple_config->couples[i];

  // create if couple does not exists
  if(!couple){
    couple_config->couples =
      xrealloc(
        couple_config->couples,
        (couple_config->num_couples + 1) * sizeof(*couple_config->couples));
    couple = &couple_config->couples[couple_config->num_couples];
    couple_config->num_couples++;
    couple->component_indices[0] = component_indices[0];
    couple->component_indices[1] = component_indices[1];
    couple->num_field_couples = 0;
    couple->field_couples = NULL;
  }

  // get field indices
  size_t src_field_idx =
    yac_couple_config_component_add_field_(
      couple_config, src_comp_idx, src_grid_idx,
      src_field_name, NULL, SIZE_MAX);
  size_t tgt_field_idx =
    yac_couple_config_component_add_field_(
      couple_config, tgt_comp_idx, tgt_grid_idx,
      tgt_field_name, NULL, SIZE_MAX);

  // check if field_couple exists
  for(size_t i = 0; i < couple->num_field_couples; ++i)
    YAC_ASSERT_F(
      !(couple->field_couples[i].source.component_idx == src_comp_idx &&
        couple->field_couples[i].target.component_idx == tgt_comp_idx &&
        couple->field_couples[i].source.field_idx == src_field_idx &&
        couple->field_couples[i].target.field_idx == tgt_field_idx),
      "ERROR(yac_couple_config_def_couple): Coupling is already defined \n"
      "source (comp_name: \"%s\" grid_name: \"%s\" field_name: \"%s\")\n"
      "target (comp_name: \"%s\" grid_name: \"%s\" field_name: \"%s\")\n",
      src_comp_name, src_grid_name, src_field_name,
      tgt_comp_name, tgt_grid_name, tgt_field_name);

  couple->field_couples = xrealloc(couple->field_couples,
    (couple->num_field_couples + 1) * sizeof(*couple->field_couples));
  struct yac_couple_config_field_couple * field_couple =
    couple->field_couples + couple->num_field_couples;
  couple->num_field_couples++;
  field_couple->source.component_idx = src_comp_idx;
  field_couple->source.field_idx = src_field_idx;
  field_couple->source.lag = src_lag;
  field_couple->target.component_idx = tgt_comp_idx;
  field_couple->target.field_idx = tgt_field_idx;
  field_couple->target.lag = tgt_lag;
  field_couple->coupling_period = strdup(coupling_period);
  field_couple->coupling_period_operation =
    (enum yac_reduction_type)time_reduction;
  field_couple->interp_stack = yac_interp_stack_config_copy(interp_stack);
  field_couple->mapping_on_source = mapping_on_source;
  field_couple->weight_file_name =
    weight_file_name?strdup(weight_file_name):NULL;
  field_couple->scale_factor = scale_factor;
  field_couple->scale_summand = scale_summand;
  field_couple->enforce_write_weight_file = weight_file_name != NULL;
  field_couple->num_src_mask_names = num_src_mask_names;
  if (num_src_mask_names > 0) {
    field_couple->src_mask_names =
      xmalloc(num_src_mask_names * sizeof(*src_mask_names));
    for (size_t i = 0; i < num_src_mask_names; ++i)
      field_couple->src_mask_names[i] = strdup(src_mask_names[i]);
  } else {
    field_couple->src_mask_names = NULL;
  }
  field_couple->tgt_mask_name =
    (tgt_mask_name != NULL)?strdup(tgt_mask_name):NULL;
}

static void couple_config_sync_time(
  struct yac_couple_config * couple_config, MPI_Comm comm) {

  couple_config_sync_calendar(comm);
  couple_config_sync_string(
    "start_datetime", (char**)&(couple_config->start_datetime), comm);
  couple_config_sync_string(
    "end_datetime", (char**)&(couple_config->end_datetime), comm);
}

struct dist_merge_vtable {
  size_t(*get_pack_size)(size_t, void*, MPI_Comm);
  void(*pack)(size_t, void*, void *, int, int *, MPI_Comm);
  void(*unpack)(void *, int, int *, size_t *, void *, MPI_Comm);
  int(*compare)(void const*, void const*);
  void(*merge)(void*, void *, MPI_Comm);
  void(*free_data)(void*);
} dist_merge_vtable_grid =
    {.get_pack_size = yac_couple_config_get_grids_pack_size,
     .pack = yac_couple_config_pack_grids,
     .unpack = yac_couple_config_unpack_grids,
     .compare = yac_couple_config_grid_compare,
     .merge = yac_couple_config_grid_merge,
     .free_data = yac_couple_config_grid_free},
  dist_merge_vtable_component =
    {.get_pack_size = yac_couple_config_get_components_pack_size,
     .pack = yac_couple_config_pack_components,
     .unpack = yac_couple_config_unpack_components,
     .compare = yac_couple_config_component_compare,
     .merge = yac_couple_config_component_merge,
     .free_data = yac_couple_config_component_free},
  dist_merge_vtable_field =
    {.get_pack_size = yac_couple_config_get_fields_pack_size,
     .pack = yac_couple_config_pack_fields,
     .unpack = yac_couple_config_unpack_fields,
     .compare = yac_couple_config_field_compare,
     .merge = yac_couple_config_field_merge,
     .free_data = yac_couple_config_field_free},
  dist_merge_vtable_field_couple =
    {.get_pack_size = yac_couple_config_get_field_couples_pack_size,
     .pack = yac_couple_config_pack_field_couples,
     .unpack = yac_couple_config_unpack_field_couples,
     .compare = yac_couple_config_field_couple_compare,
     .merge = yac_couple_config_field_couple_merge,
     .free_data = yac_couple_config_field_couple_free},
  dist_merge_vtable_couple =
    {.get_pack_size = yac_couple_config_get_couples_pack_size_basic,
     .pack = yac_couple_config_pack_couples_basic,
     .unpack = yac_couple_config_unpack_couples_basic,
     .compare = yac_couple_config_couple_compare_basic,
     .merge = yac_couple_config_couple_merge,
     .free_data = yac_couple_config_couple_free};

static void dist_merge(size_t* len, void** arr, size_t element_size, MPI_Comm comm,
  struct dist_merge_vtable* vtable, size_t** idx_old_to_new){
  int rank;
  MPI_Comm_rank(comm, &rank);

  // initialise
  unsigned char * input = *arr;
  size_t input_len = *len;
  *idx_old_to_new = xmalloc(input_len * sizeof(size_t));
  size_t* idx = xmalloc(input_len * sizeof(size_t));
  for(size_t i = 0;i<input_len;++i) idx[i] = i;
  unsigned char * arr_new = NULL;
  size_t len_new = 0;

  // sort input data
  yac_qsort_index(input, input_len, element_size, vtable->compare, idx);

  void * buffer = NULL;
  while(1){

    // determine pack size of local data
    size_t pack_size =
      (input_len > 0)?vtable->get_pack_size(input_len, input, comm):0;
    YAC_ASSERT(
      pack_size <= LONG_MAX, "ERROR(dist_merge): packing size too big");

    // determine rank with most amount of data
    struct {
      long pack_size;
      int rank;
    } data_pair = {.pack_size = (long)pack_size, .rank = rank};
    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, &data_pair, 1, MPI_LONG_INT, MPI_MAXLOC, comm), comm);

    // if there is no more data to exchange
    if (data_pair.pack_size == 0) break;

    // pack data into buffer buffer
    pack_size = (size_t)data_pair.pack_size;
    if (!buffer) buffer = xmalloc(pack_size);
    int position = 0;
    if(data_pair.rank == rank)
      vtable->pack(input_len, input, buffer, pack_size, &position, comm);

    // broadcast and unpack data
    yac_mpi_call(
      MPI_Bcast(buffer, pack_size, MPI_PACKED, data_pair.rank, comm), comm);
    void * recved = NULL;
    size_t num_recved;
    position = 0;
    vtable->unpack(buffer, pack_size, &position, &num_recved, &recved, comm);

    // add received elements to list
    arr_new = xrealloc(arr_new, (len_new + num_recved)*element_size);
    memcpy(
      arr_new + len_new*element_size, recved, num_recved*element_size);
    free(recved);

    // for all received elements
    size_t input_idx, input_len_new, i = len_new;
    len_new += num_recved;
    for(input_idx = 0, input_len_new = 0; i < len_new; ++i) {
      void* recved_element = arr_new + i*element_size;

      // search for matching element in input list
      int cmp = 0;
      void * input_element = NULL;
      while ((input_idx < input_len) &&
             (((cmp =
                  vtable->compare(
                    ((input_element = input + input_idx * element_size)),
                    recved_element))) < 0)) {

        if (input_idx != input_len_new) {
          memcpy(
            input + input_len_new * element_size,
            input_element, element_size);
          idx[input_len_new] = idx[input_idx];
        }
        ++input_len_new;
        ++input_idx;
      }

      // if the end of the local input list was reached
      if (input_idx == input_len) break;

      // if a matching element was found in the input list
      if (!cmp) {

        // merge input list element with received element
        if (vtable->merge)
          vtable->merge(recved_element, input_element, comm);
        // remove element
        vtable->free_data(input_element);
        // update index
        (*idx_old_to_new)[idx[input_idx]] = i;
        // upate input idx
        input_idx++;

      // if no matching element was found in the input list
      } else if (vtable->merge) {
        // since the merge operation can potentially be collective, we have
        // to call it, even if no matching element was found
        vtable->merge(recved_element, NULL, comm);
      }
    }
    // process remaining received elements
    if (vtable->merge)
      for(; i < len_new; ++i)
        vtable->merge(arr_new + i*element_size, NULL, comm);
    // pack remaining elements in input list
    for (; input_idx < input_len; ++input_idx, ++input_len_new) {
      if (input_idx != input_len_new) {
        memcpy(
          input + input_len_new * element_size,
          input + input_idx * element_size, element_size);
        idx[input_len_new] = idx[input_idx];
      }
    }
    // update length of input list
    input_len = input_len_new;
  }
  free(buffer);
  free(input);
  free(idx);
  *arr = arr_new;
  *len = len_new;
}

static void merge_grids(
  struct yac_couple_config * couple_config, MPI_Comm comm) {

  size_t* old_to_new_idx;
  void * p_grids = couple_config->grids;
  dist_merge(
    &couple_config->num_grids, &p_grids,
      sizeof(couple_config->grids[0]),
      comm, &dist_merge_vtable_grid, &old_to_new_idx);
  couple_config->grids = p_grids;

  // set new grid_idx in fields
  for(size_t comp_idx = 0; comp_idx < couple_config->num_components;
      ++comp_idx) {
    struct yac_couple_config_component * component =
      couple_config->components + comp_idx;
    for(size_t field_idx = 0; field_idx < component->num_fields; ++field_idx)
      component->fields[field_idx].grid_idx =
        old_to_new_idx[component->fields[field_idx].grid_idx];
  }

  free(old_to_new_idx);
}

static void merge_fields(
  struct yac_couple_config * couple_config, size_t comp_idx, MPI_Comm comm) {

  struct yac_couple_config_component * component =
    couple_config->components + comp_idx;

  size_t* old_to_new_idx;
  void * p_fields = component->fields;
  dist_merge(
    &component->num_fields, &p_fields,
    sizeof(component->fields[0]),
    comm, &dist_merge_vtable_field, &old_to_new_idx);
  component->fields = p_fields;

  // set new field_idx in all field_couples
  for(size_t couple_idx = 0; couple_idx < couple_config->num_couples;
      ++couple_idx) {
    struct yac_couple_config_couple * couple =
      couple_config->couples + couple_idx;
    for(size_t field_couple_idx = 0;
        field_couple_idx < couple->num_field_couples; ++field_couple_idx){
      struct yac_couple_config_field_couple * field_couple =
        couple->field_couples + field_couple_idx;
      if(field_couple->source.component_idx == comp_idx)
        field_couple->source.field_idx =
          old_to_new_idx[field_couple->source.field_idx];
      if(field_couple->target.component_idx == comp_idx)
        field_couple->target.field_idx =
          old_to_new_idx[field_couple->target.field_idx];
    }
  }

  free(old_to_new_idx);
}

static void merge_components(
  struct yac_couple_config * couple_config, MPI_Comm comm) {

  // distribute and merge basic component information while keeping the
  // individual fields
  size_t* old_to_new_idx;
  void * p_components = couple_config->components;
  dist_merge(
    &couple_config->num_components, &p_components,
      sizeof(couple_config->components[0]),
      comm, &dist_merge_vtable_component, &old_to_new_idx);
  couple_config->components = p_components;

  // set new component_idx in couples
  for(size_t couple_idx = 0; couple_idx < couple_config->num_couples;
      ++couple_idx) {
    struct yac_couple_config_couple * couple =
      couple_config->couples + couple_idx;
    couple->component_indices[0] = old_to_new_idx[couple->component_indices[0]];
    couple->component_indices[1] = old_to_new_idx[couple->component_indices[1]];
    for(size_t field_couple_idx = 0;
        field_couple_idx < couple->num_field_couples; ++field_couple_idx) {
      couple->field_couples[field_couple_idx].source.component_idx =
        old_to_new_idx[
          couple->field_couples[field_couple_idx].source.component_idx];
      couple->field_couples[field_couple_idx].target.component_idx =
        old_to_new_idx[
          couple->field_couples[field_couple_idx].target.component_idx];
    }
  }
  free(old_to_new_idx);

  // merge the fields of each component
  for (size_t comp_idx = 0; comp_idx < couple_config->num_components; ++comp_idx)
    merge_fields(couple_config, comp_idx, comm);
}

static void merge_field_couples(
  size_t * num_field_couples,
  struct yac_couple_config_field_couple ** field_couples, MPI_Comm comm) {

  size_t* old_to_new_idx;
  void * p_field_couples = *field_couples;
  dist_merge(
    num_field_couples, &p_field_couples, sizeof(**field_couples),
    comm, &dist_merge_vtable_field_couple, &old_to_new_idx);
  free(old_to_new_idx);
  *field_couples = p_field_couples;
}

static void merge_couples(
  struct yac_couple_config * couple_config, MPI_Comm comm) {

  // distribute and merge basic couple information while keeping the
  // individual field couples
  size_t* old_to_new_idx;
  void * p_couples = couple_config->couples;
  dist_merge(
    &couple_config->num_couples, &p_couples,
    sizeof(couple_config->couples[0]),
    comm, &dist_merge_vtable_couple, &old_to_new_idx);
  couple_config->couples = p_couples;
  free(old_to_new_idx);
}

void yac_couple_config_sync(
  struct yac_couple_config * couple_config, MPI_Comm comm){

  // sync time stuff
  couple_config_sync_time(couple_config, comm);

  merge_grids(couple_config, comm);
  merge_components(couple_config, comm);
  merge_couples(couple_config, comm);
}
