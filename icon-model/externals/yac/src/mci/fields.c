// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <string.h>
#include "event.h"
#include "utils_mci.h"
#include "utils_common.h"
#include "fields.h"

struct coupling_field {

  // general field related information

  char * name;

  char * component_name;

  struct yac_basic_grid * grid;

  struct yac_interp_field * interp_fields;
  size_t num_interp_fields;
  size_t collection_size; //!< number of vertical levels or bundles
  char* timestep;

  char current_datetime[MAX_DATETIME_STR_LEN]; // updated only in yac_coupling_field_get_datetime

  // data exchange

  enum yac_field_exchange_type exchange_type;

  union {
    struct {
      struct {
        struct event * event;
        struct yac_interpolation * interpolation;
        int time_accumulation_count;
        double *** send_field_acc;
        double *** send_frac_mask_acc;
      } * puts;
      unsigned num_puts;
      int ** mask;
    } put;
    struct {
      struct event * event;
      struct yac_interpolation * interpolation;
      int * mask;
    } get;
  } exchange_data;
};

struct coupling_field *
yac_coupling_field_new(char const * field_name, char const * component_name,
  struct yac_basic_grid * grid, struct yac_interp_field * interp_fields,
  unsigned num_interp_fields, size_t collection_size, const char* timestep) {

  struct coupling_field * cpl_field = xmalloc(1 * sizeof(*cpl_field));

  cpl_field->name              = xmalloc((strlen(field_name) + 1));
  strcpy(cpl_field->name, field_name);
  cpl_field->component_name    = component_name?strdup(component_name):NULL;
  cpl_field->grid              = grid;
  cpl_field->interp_fields     = NULL;
  cpl_field->num_interp_fields = (size_t)num_interp_fields;
  cpl_field->collection_size   = collection_size;
  cpl_field->timestep          = strdup(timestep);
  cpl_field->exchange_type     = NOTHING;

  if (num_interp_fields != 0) {

    YAC_ASSERT(
      interp_fields != NULL,
      "ERROR: num_interp_fields > 0 but not interpolation field provided")

    cpl_field->interp_fields =
      xmalloc((size_t)num_interp_fields * sizeof(*(cpl_field->interp_fields)));
    memcpy(cpl_field->interp_fields, interp_fields,
           (size_t)num_interp_fields * sizeof(*(cpl_field->interp_fields)));

  }

  return cpl_field;
}

unsigned yac_get_coupling_field_collection_size(struct coupling_field * field) {

   return field->collection_size;
}

const char* yac_get_coupling_field_timestep(struct coupling_field * field) {

   return field->timestep;
}

size_t yac_coupling_field_get_num_interp_fields(struct coupling_field * field) {

  return field->num_interp_fields;
}

enum yac_location
yac_get_coupling_field_get_interp_field_location(
  struct coupling_field * field, size_t interp_field_idx) {

  return field->interp_fields[interp_field_idx].location;
}

const char * yac_get_coupling_field_name(struct coupling_field * field) {

  return field->name;
}

struct yac_basic_grid * yac_coupling_field_get_basic_grid(
  struct coupling_field * field) {

  return field->grid;
}

char const * yac_get_coupling_field_comp_name(struct coupling_field * field) {

  return field->component_name;
}

enum yac_field_exchange_type
yac_get_coupling_field_exchange_type(struct coupling_field * field) {

  return field->exchange_type;
}

static size_t get_interp_field_size(
  struct yac_interp_field field, struct yac_basic_grid * grid) {

  return yac_basic_grid_get_data_size(grid, field.location);
}

size_t yac_coupling_field_get_data_size(
  struct coupling_field * field, enum yac_location location) {

  return yac_basic_grid_get_data_size(field->grid, location);
}

static void init_put_op_acc(
  struct coupling_field * field, double init_value, double *** acc) {

  size_t num_interp_fields = field->num_interp_fields;
  size_t collection_size = field->collection_size;

  for (size_t h = 0; h < collection_size; ++h) {
    for (size_t j = 0; j < num_interp_fields; j++) {
      size_t num_points =
        get_interp_field_size(field->interp_fields[j], field->grid);
      YAC_OMP_PARALLEL
      {
        YAC_OMP_FOR
        for (size_t m = 0; m < num_points; m++ ) {
          acc[h][j][m] = init_value;
        }
      }
    }
  }
}

static void get_put_op_acc(
  struct coupling_field * field, unsigned put_idx,
  char const * routine_name, double **** acc) {

  YAC_ASSERT_F(
    field->exchange_type == SOURCE,
    "ERROR(%s): wrong field exchange type", routine_name)

  YAC_ASSERT_F(
    put_idx < field->exchange_data.put.num_puts,
    "ERROR(%s): put_idx is invalid", routine_name)

  if (*acc != NULL) return;

  size_t num_interp_fields = field->num_interp_fields;
  size_t collection_size = field->collection_size;

  /* Allocate memory for the accumulation */

  *acc = xmalloc(collection_size * sizeof(**acc));

  for (size_t h = 0; h < collection_size; ++h) {

    (*acc)[h] = xmalloc(num_interp_fields * sizeof(***acc));

    for (size_t j = 0; j < num_interp_fields; j++)
      (*acc)[h][j] =
        xmalloc(
          get_interp_field_size(
            field->interp_fields[j], field->grid) * sizeof(****acc));
  }

  /* Initialise memory for the accumulation */
  init_put_op_acc(field, 0.0, *acc);
}

double *** yac_get_coupling_field_put_op_send_field_acc(
  struct coupling_field * field, unsigned put_idx) {

  get_put_op_acc(
    field, put_idx, "yac_get_coupling_field_put_op_send_field_acc",
    &(field->exchange_data.put.puts[put_idx].send_field_acc));

  return field->exchange_data.put.puts[put_idx].send_field_acc;
}

double *** yac_get_coupling_field_put_op_send_frac_mask_acc(
  struct coupling_field * field, unsigned put_idx) {

  get_put_op_acc(
    field, put_idx, "yac_get_coupling_field_put_op_send_frac_mask_acc",
    &(field->exchange_data.put.puts[put_idx].send_frac_mask_acc));

  return field->exchange_data.put.puts[put_idx].send_frac_mask_acc;
}

static int check_core_mask(
  struct coupling_field * field, size_t field_idx) {

  int const * core_mask =
    yac_basic_grid_get_core_mask(
      field->grid, field->interp_fields[field_idx].location);

  int flag = 0;

  // if the core mask is defined
  if (core_mask != NULL) {

    size_t field_size =
      get_interp_field_size(field->interp_fields[field_idx], field->grid);

    // check if the core mask contains values other than 1
    for (size_t i = 0; (i < field_size) && !flag; ++i)
      flag = core_mask[i] != 1;
  }

  return flag;
}

static int check_field_mask(
  struct coupling_field * field, size_t field_idx) {

  int const * field_mask =
    yac_basic_grid_get_field_mask(
      field->grid, field->interp_fields[field_idx]);

  int flag = 0;

  // if a field mask is defined
  if (field_mask != NULL) {

    size_t field_size =
      get_interp_field_size(field->interp_fields[field_idx], field->grid);

    // check if the core mask contains values other than 1
    for (size_t i = 0; (i < field_size) && !flag; ++i)
      flag = field_mask[i] != 1;
  }

  return flag;
}

int ** yac_get_coupling_field_put_mask(struct coupling_field * field) {

  YAC_ASSERT(
    field->exchange_type == SOURCE,
    "ERROR(yac_get_coupling_field_put_mask): "
    "wrong field exchange type")

  if (field->exchange_data.put.mask != NULL)
    return field->exchange_data.put.mask;

  size_t num_interp_fields = field->num_interp_fields;

  // check whether for this put operation a field or core mask has been defined
  int has_mask = 0;
  for (size_t i = 0; (i < num_interp_fields) && !has_mask; ++i)
    has_mask = check_core_mask(field, i) || check_field_mask(field, i);

  int ** mask = NULL;

  if (has_mask) {

    mask = xmalloc(num_interp_fields * sizeof(*mask));

    // generate mask from core and field mask
    for (size_t i = 0; i < num_interp_fields; ++i) {
      int const * core_mask =
        yac_basic_grid_get_core_mask(
          field->grid, field->interp_fields[i].location);
      int const * field_mask =
        yac_basic_grid_get_field_mask(
          field->grid, field->interp_fields[i]);
      size_t field_size =
        get_interp_field_size(field->interp_fields[i], field->grid);
      mask[i] = xmalloc(field_size * sizeof(**mask));
      for (size_t j = 0; j < field_size; ++j)
        mask[i][j] = ((core_mask == NULL) || (core_mask[j] != 0)) &&
                     ((field_mask == NULL) || (field_mask[j] != 0));
    }
  }

  field->exchange_data.put.mask = mask;

  return mask;
}

int * yac_get_coupling_field_get_mask(struct coupling_field * field) {

  YAC_ASSERT(
    field->exchange_type == TARGET,
    "ERROR(yac_get_coupling_field_put_mask): "
    "wrong field exchange type")

  if (field->exchange_data.get.mask != NULL)
    return field->exchange_data.get.mask;

  // check whether for this put operation a field or core mask has been defined
  int has_mask = check_core_mask(field, 0) || check_field_mask(field, 0);

  int * mask = NULL;

  if (has_mask) {

    // generate mask from core and field mask
    int const * core_mask =
      yac_basic_grid_get_core_mask(
        field->grid, field->interp_fields[0].location);
    int const * field_mask =
      yac_basic_grid_get_field_mask(
        field->grid, field->interp_fields[0]);
    size_t field_size =
      get_interp_field_size(field->interp_fields[0], field->grid);
    mask = xmalloc(field_size * sizeof(*mask));
    for (size_t i = 0; i < field_size; ++i)
      mask[i] = ((core_mask == NULL) || (core_mask[i] != 0)) &&
                 ((field_mask == NULL) || (field_mask[i] != 0));
  }

  field->exchange_data.get.mask = mask;

  return mask;
}

void yac_init_coupling_field_put_op_send_field_acc(
  struct coupling_field * field, unsigned put_idx, double init_value) {

  YAC_ASSERT(
    field->exchange_type == SOURCE,
    "ERROR(yac_init_coupling_field_put_op_send_field_acc): "
    "wrong field exchange type")

  YAC_ASSERT(
    put_idx < field->exchange_data.put.num_puts,
    "ERROR(yac_init_coupling_field_put_op_send_field_acc): "
    "put_idx is invalid")

  double *** send_field_acc =
    (field->exchange_data.put.puts[put_idx].send_field_acc == NULL)?
    yac_get_coupling_field_put_op_send_field_acc(field, put_idx):
    field->exchange_data.put.puts[put_idx].send_field_acc;

  init_put_op_acc(field, init_value, send_field_acc);
}

void yac_init_coupling_field_put_op_send_frac_mask_acc(
  struct coupling_field * field, unsigned put_idx, double init_value) {

  YAC_ASSERT(
    field->exchange_type == SOURCE,
    "ERROR(yac_init_coupling_field_put_op_send_frac_mask_acc): "
    "wrong field exchange type")

  YAC_ASSERT(
    put_idx < field->exchange_data.put.num_puts,
    "ERROR(yac_init_coupling_field_put_op_send_frac_mask_acc): "
    "put_idx is invalid")

  double *** send_frac_mask_acc =
    (field->exchange_data.put.puts[put_idx].send_frac_mask_acc == NULL)?
    yac_get_coupling_field_put_op_send_frac_mask_acc(field, put_idx):
    field->exchange_data.put.puts[put_idx].send_frac_mask_acc;

  init_put_op_acc(field, init_value, send_frac_mask_acc);
}

struct event * yac_get_coupling_field_put_op_event(
  struct coupling_field * field, unsigned put_idx) {

  YAC_ASSERT(
    field->exchange_type == SOURCE, "ERROR: wrong field exchange type")

  YAC_ASSERT(
    put_idx < field->exchange_data.put.num_puts, "ERROR: put_idx is invalid")

  return field->exchange_data.put.puts[put_idx].event;
}

struct yac_interpolation * yac_get_coupling_field_put_op_interpolation(
  struct coupling_field * field, unsigned put_idx) {

  YAC_ASSERT(
    field->exchange_type == SOURCE, "ERROR: wrong field exchange type")

  YAC_ASSERT(
    put_idx < field->exchange_data.put.num_puts, "ERROR: put_idx is invalid")

  return field->exchange_data.put.puts[put_idx].interpolation;
}

int yac_get_coupling_field_put_op_time_accumulation_count(
  struct coupling_field * field, unsigned put_idx) {

  YAC_ASSERT(
    field->exchange_type == SOURCE, "ERROR: wrong field exchange type")

  YAC_ASSERT(
    put_idx < field->exchange_data.put.num_puts, "ERROR: put_idx is invalid")

  return field->exchange_data.put.puts[put_idx].time_accumulation_count;
}

void
yac_set_coupling_field_put_op_time_accumulation_count(struct coupling_field * field,
                                                      unsigned put_idx, int count) {

  YAC_ASSERT(
    field->exchange_type == SOURCE, "ERROR: wrong field exchange type")

  YAC_ASSERT(
    put_idx < field->exchange_data.put.num_puts, "ERROR: put_idx is invalid")

  field->exchange_data.put.puts[put_idx].time_accumulation_count = count;
}

unsigned yac_get_coupling_field_num_puts(struct coupling_field * field) {

  YAC_ASSERT(
    field->exchange_type == SOURCE, "ERROR: wrong field exchange type")

  return field->exchange_data.put.num_puts;
}

struct event * yac_get_coupling_field_get_op_event(
  struct coupling_field * field) {

  YAC_ASSERT(
    field->exchange_type == TARGET, "ERROR: wrong field exchange type")

  return field->exchange_data.get.event;
}

struct yac_interpolation *
yac_get_coupling_field_get_op_interpolation(struct coupling_field * field) {

  YAC_ASSERT(
    field->exchange_type == TARGET, "ERROR: wrong field exchange type")

  return field->exchange_data.get.interpolation;
}

void yac_set_coupling_field_put_op(
  struct coupling_field * field, struct event * event,
  struct yac_interpolation * interpolation) {

  YAC_ASSERT_F(
    field->exchange_type != TARGET,
    "ERROR: a get operation has already been set (field \"%s\")",
    field->name)

  if (field->exchange_type == NOTHING) {
    field->exchange_data.put.puts = NULL;
    field->exchange_data.put.num_puts = 0;
    field->exchange_type = SOURCE;
    field->exchange_data.get.mask = NULL;
  }

  unsigned num_puts = field->exchange_data.put.num_puts;
  field->exchange_data.put.num_puts++;

  field->exchange_data.put.puts =
    xrealloc(field->exchange_data.put.puts, (num_puts + 1) *
             sizeof(*(field->exchange_data.put.puts)));

  field->exchange_data.put.puts[num_puts].event = event;
  field->exchange_data.put.puts[num_puts].interpolation = interpolation;
  field->exchange_data.put.puts[num_puts].time_accumulation_count = 0;
  field->exchange_data.put.puts[num_puts].send_field_acc = NULL;
  field->exchange_data.put.puts[num_puts].send_frac_mask_acc = NULL;
}

void yac_set_coupling_field_get_op(
  struct coupling_field * field, struct event * event,
  struct yac_interpolation * interpolation) {

  YAC_ASSERT_F(
    field->exchange_type == NOTHING,
    "ERROR: a put or get operation has already been set (field \"%s\")",
    field->name)

  field->exchange_type = TARGET;
  field->exchange_data.get.event = event;
  field->exchange_data.get.interpolation = interpolation;
  field->exchange_data.put.mask = NULL;
}

struct yac_interp_field const * yac_coupling_field_get_interp_fields(
  struct coupling_field * cpl_field) {

  return cpl_field->interp_fields;
}

char* yac_coupling_field_get_datetime(
  struct coupling_field * cpl_field) {
    YAC_ASSERT_F(
    (cpl_field->exchange_type == SOURCE) ||
    (cpl_field->exchange_type == TARGET),
      "ERROR(yac_coupling_field_get_datetime): "
      "datetime is not defined for non-exchanged fields "
      "(field \"%s\")", cpl_field->name);
  struct event * event =
    (cpl_field->exchange_type == SOURCE)?
      cpl_field->exchange_data.put.puts[0].event:
      cpl_field->exchange_data.get.event;
  return
    yac_get_event_current_datetime(
      event, cpl_field->current_datetime);
}

void yac_coupling_field_delete(struct coupling_field * cpl_field) {

  size_t collection_size = cpl_field->collection_size;
  size_t num_interp_fields = cpl_field->num_interp_fields;

  if (cpl_field->exchange_type == SOURCE) {

    for (unsigned put_idx = 0; put_idx < cpl_field->exchange_data.put.num_puts;
         ++put_idx) {

      struct event * event = cpl_field->exchange_data.put.puts[put_idx].event;
      struct yac_interpolation * interpolation =
        cpl_field->exchange_data.put.puts[put_idx].interpolation;
      double ***send_field_acc =
        cpl_field->exchange_data.put.puts[put_idx].send_field_acc;
      double ***send_frac_mask_acc =
        cpl_field->exchange_data.put.puts[put_idx].send_frac_mask_acc;

      if (send_field_acc != NULL) {
        for (size_t h = 0; h < collection_size; ++h) {
          for (size_t j = 0; j < num_interp_fields; ++j)
            free(send_field_acc[h][j]);
          free(send_field_acc[h]);
        }
        free(send_field_acc);
      }

      if (send_frac_mask_acc != NULL) {
        for (size_t h = 0; h < collection_size; ++h) {
          for (size_t j = 0; j < num_interp_fields; ++j)
            free(send_frac_mask_acc[h][j]);
          free(send_frac_mask_acc[h]);
        }
        free(send_frac_mask_acc);
      }

      yac_event_delete(event);
      yac_interpolation_delete(interpolation);
    }

    free(cpl_field->exchange_data.put.puts);
    if (cpl_field->exchange_data.put.mask) {
      for (size_t i = 0; i < num_interp_fields; ++i)
        free(cpl_field->exchange_data.put.mask[i]);
      free(cpl_field->exchange_data.put.mask);
    }

  } else if (cpl_field->exchange_type == TARGET) {

    yac_event_delete(cpl_field->exchange_data.get.event);
    yac_interpolation_delete(cpl_field->exchange_data.get.interpolation);
    if (cpl_field->exchange_data.get.mask)
      free(cpl_field->exchange_data.get.mask);
  }

  free(cpl_field->name);
  free(cpl_field->component_name);
  free(cpl_field->interp_fields);
  free(cpl_field->timestep);
  free(cpl_field);
}
