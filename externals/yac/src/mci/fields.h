// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIELDS_H
#define FIELDS_H

#include "component.h"

// general coupling field struct

enum yac_field_exchange_type {
  NOTHING = 0,
  SOURCE  = 1,
  TARGET  = 2,
};

// forward declaration

struct coupling_field;
struct yac_interpolation;

// generation of a field

/**
 * Constructs a coupling field
 *
 * @param[in] field_name        name of the coupling field
 * @param[in] component_name    component name
 * @param[in] grid              grid
 * @param[in] interp_fields     interpolation fields
 * @param[in] num_interp_fields number of entries in interp_fields
  * @param[in] collection_size
  * @param[in] timestep
 */
struct coupling_field *
yac_coupling_field_new(char const * field_name, char const * component_name,
  struct yac_basic_grid * grid, struct yac_interp_field * interp_fields,
  unsigned num_interp_fields, size_t collection_size, const char* timestep);

/**
 * gets the number of vertical levels or bundles for a field.
 *
 * @param[in] field
 * @return collection size of the field
 */
unsigned yac_get_coupling_field_collection_size(
  struct coupling_field * field);

/**
 * gets the timestep or bundles for a field.
 *
 * @param[in] field
 * @return timestep
 */
const char* yac_get_coupling_field_timestep(
  struct coupling_field * field);

/**
 * gets the name of the coupling field
 * @param[in] field
 * @return name of the coupling field
 */
const char * yac_get_coupling_field_name(struct coupling_field * field);

/**
 * gets the component name of the coupling field
 * @param[in] field
 * @return component name
 */
char const * yac_get_coupling_field_comp_name(struct coupling_field * field);

/**
 * gets the grid_data of the coupling field
 * @param[in] field
 * @return basic grid data of the coupling field
 */
struct yac_basic_grid * yac_coupling_field_get_basic_grid(
  struct coupling_field * field);

/**
 * get interpolation field data location
 * @param[in] field
 * @param[in] interp_field_idx index of interpolation field
 * @return data location of requested interpolation field
 */
enum yac_location
yac_get_coupling_field_get_interp_field_location(
  struct coupling_field * field, size_t interp_field_idx);

/**
 * gets the interpolation fields of the coupling field
 * @param[in] field
 * @return interpolation fields
 * @remarks dimensions: [point set idx]
 */
struct yac_interp_field const * yac_coupling_field_get_interp_fields(
  struct coupling_field * field);

/**
 * gets the number of interpolation fields of the coupling field
 * @param[in] field
 * @return number of interpolation fields
 */
size_t yac_coupling_field_get_num_interp_fields(struct coupling_field * field);

/**
 * gets the number of points in the grids associated with the coupling field
 * @param[in] field
 * @param[in] location
 * @return number of grid points
 */
size_t yac_coupling_field_get_data_size(
  struct coupling_field * field, enum yac_location location);

/**
 * gets the exchange type of the coupling field
 * @param[in] field
 * @return exchange type of the coupling field
 */
enum yac_field_exchange_type
yac_get_coupling_field_exchange_type(struct coupling_field * field);

/**
 * gets the number of put operations of the coupling field
 * @param[in] field
 * @return number of put operations of the coupling field
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 */
unsigned yac_get_coupling_field_num_puts(struct coupling_field * field);

/**
 * gets the event of a specified put operation
 * @param[in] field
 * @param[in] put_idx index of the put operation
 * @return event of specified put operation
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @see get_coupling_field_num_puts
 */
struct event * yac_get_coupling_field_put_op_event(
  struct coupling_field * field, unsigned put_idx);

/**
 * gets the field interpolation of a specified put operation
 * @param[in] field
 * @param[in] put_idx index of the put operation
 * @return field interpolation of specified put operation
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @see get_coupling_field_num_puts
 */
struct yac_interpolation *
yac_get_coupling_field_put_op_interpolation(struct coupling_field * field,
                                            unsigned put_idx);

/**
 * gets the send field accumulator of a specified put operation
 * @param[in] field
 * @param[in] put_idx index of the put operation
 * @return send field accumulator of specified put operation
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @remarks when this routine is called the first time for a put operation
 *          it will be allocated and initialised with zeros
 * @see get_coupling_field_num_puts
 */
double ***
yac_get_coupling_field_put_op_send_field_acc(struct coupling_field * field,
                                             unsigned put_idx);

/**
 * gets the fractional mask accumulator of a specified put operation
 * @param[in] field
 * @param[in] put_idx index of the put operation
 * @return fractional mask accumulator of specified put operation
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @remarks when this routine is called the first time for a put operation
 *          it will be allocated and initialised with zeros
 * @see get_coupling_field_num_puts
 */
double ***
yac_get_coupling_field_put_op_send_frac_mask_acc(struct coupling_field * field,
                                                 unsigned put_idx);

/**
 * initialises the send field accumulator of a specified put operation
 * @param[in] field
 * @param[in] put_idx    index of the put operation
 * @param[in] init_value value that is to be used for the initialisation of the
 *                       send field accumulator
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @see get_coupling_field_num_puts
 */
void
yac_init_coupling_field_put_op_send_field_acc(struct coupling_field * field,
                                              unsigned put_idx,
                                              double init_value);

/**
 * initialises the fractional mask accumulator of a specified put operation
 * @param[in] field
 * @param[in] put_idx    index of the put operation
 * @param[in] init_value value that is to be used for the initialisation of the
 *                       fractional mask accumulator
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @see get_coupling_field_num_puts
 */
void
yac_init_coupling_field_put_op_send_frac_mask_acc(struct coupling_field * field,
                                                  unsigned put_idx,
                                                  double init_value);

/**
 * gets the time accumulation count of a specified put operation
 * @param[in] field
 * @param[in] put_idx    index of the put operation
 * @return time accumulation count of specified put operation
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @see get_coupling_field_num_puts
 */
int
yac_get_coupling_field_put_op_time_accumulation_count(struct coupling_field * field,
                                                      unsigned put_idx);

/**
 * sets the time accumulation count of a specified put operation
 * @param[in] field
 * @param[in] put_idx index of the put operation
 * @param[in] count   new time accumulation count
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than SOURCE yields
 *          undefined results
 * @see get_coupling_field_num_puts
 */
void
yac_set_coupling_field_put_op_time_accumulation_count(struct coupling_field * field,
                                                      unsigned put_idx,
                                                      int count);

/**
 * gets a combined core and field mask
 * @param[in] field
 * @return combined core and field mask
 */
int ** yac_get_coupling_field_put_mask(struct coupling_field * field);

/**
 * gets a combined core and field mask
 * @param[in] field
 * @return combined core and field mask
 */
int * yac_get_coupling_field_get_mask(struct coupling_field * field);

/**
 * gets the event of the get operation
 * @param[in] field
 * @return event of the get operation
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than TARGET yields
 *          undefined results
 */
struct event * yac_get_coupling_field_get_op_event(
  struct coupling_field * field);

/**
 * gets the field interpolation of the get operation
 * @param[in] field
 * @return field interpolation of the get operation
 * @remarks calling this routine for a coupling field that
 *          has a different exchange type than TARGET yields
 *          undefined results
 */
struct yac_interpolation *
yac_get_coupling_field_get_op_interpolation(struct coupling_field * field);

/**
 * sets the put operation for the coupling field
 * @param[in,out] field
 * @param[in]     event         event of the operation
 * @param[in]     interpolation interpolation to be executed
 * @remarks the coupling field takes control of the input data, it will be
 *          free/deleted by the coupling field
 */
void yac_set_coupling_field_put_op(
  struct coupling_field * field, struct event * event,
  struct yac_interpolation * interpolation);

/**
 * set the get operation for the coupling field
 * @param[in,out] field
 * @param[in]     event         event of the operation
 * @param[in]     interpolation interpolation to be executed
 * @remarks the coupling field takes control of the input data, it will be
 *          free/deleted by the coupling field
 */
void yac_set_coupling_field_get_op(
  struct coupling_field * field, struct event * event,
  struct yac_interpolation * interpolation);


/**
 * get the current datetime of the coupling field
 * @param[in] cpl_field
 * @return the datetime string
 *
 */
char* yac_coupling_field_get_datetime(
  struct coupling_field * cpl_field);

// destruction of fields
void yac_coupling_field_delete(struct coupling_field * cpl_field);

#endif // FIELDS_H

