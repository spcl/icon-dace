/**
 * @file instance.c
 * @brief Structs and interfaces to defined YAC instances
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

#ifndef INSTANCE_H
#define INSTANCE_H

#include "grid.h"
#include "component.h"
#include "yac_mpi.h"
#include "couple_config.h"

/** \example test_instance_parallel1.c
 * This example show how to set up a YAC instance. It uses three
 * processes.
 */

/** \example test_instance_parallel2.c
 * This example show how to set up a YAC instance. It uses nine
 * processes.
 */

/** \example test_instance_parallel3.c
 * This example show how to set up a YAC instance. It uses nine
 * processes.
 */

/** \example test_instance_parallel4.c
 * This example show how to set up a YAC instance. It uses nine
 * processes.
 */

struct yac_instance;
/**
 * Constructor for a yac instance.
 * @param[in] comm          MPI communicator that contains the processes of
 *                          all components that will be registered with this
 *                          yac instance
 * @return yac instance
 * @remark This is collectiv for all processes in comm.
 */
struct yac_instance * yac_instance_new(MPI_Comm comm);

/**
 * Dummy constructor, which can be called instead of \ref yac_instance_new.
 * @param[in] comm
 */
void yac_instance_dummy_new(MPI_Comm comm);


/**
 * Destructor for a yac instance
 * @param[in] instance yac instance to be deleted
 */
void yac_instance_delete(struct yac_instance * instance);

/**
 * Definition of job start and end datetime
 * @param[in] instance       yac instance
 * @param[in] start_datetime calendar job start datetime
 * @param[in] end_datetime   calendar job end datetime
 */
void yac_instance_def_datetime(
  struct yac_instance * instance, const char * start_datetime,
  const char * end_datetime );

/**
 * query routine for the start datetime of the job
 * @param[in] instance yac instance
 * @return start datetime
 */
char * yac_instance_get_start_datetime(struct yac_instance * instance);

/**
 * query routine for the end datetime of the job
 * @param[in] instance yac instance
 * @return end datetime
 */
char * yac_instance_get_end_datetime(struct yac_instance * instance);

/**
 * Get the coupling configuration data from a yac instance
 * @param[in] instance yac instance
 * @return coupling configuration data
 */
struct yac_couple_config * yac_instance_get_couple_config(
  struct yac_instance * instance);

/**
 * Sets the coupling configuration data for a yac instance
 * @param[in] instance      yac instance
 * @param[in] couple_config coupling configuration data
 */
void yac_instance_set_couple_config(
  struct yac_instance * instance,
    struct yac_couple_config * couple_config);

/**
 * Defines the components for a yac instance
 * @param[in] instance   yac instance
 * @param[in] comp_names names of components
 * @param[in] num_comps  number of entries in comp_names
 */
void yac_instance_def_components(
  struct yac_instance * instance, char const ** comp_names, size_t num_comps);

/**
 * Returns true if the components for this instance have already been defined
 * @param[in] instance yac instance
 */
int yac_instance_components_are_defined(
  struct yac_instance * instance);

/**
 * Adds a coupling field to a yac instance
 * @param[in] instance          yac instance
 * @param[in] field_name        name of the coupling field
 * @param[in] component_index   component index
 * @param[in] grid              grid
 * @param[in] interp_fields     interpolation fields
 * @param[in] num_interp_fields number of entries in interp_fields
 * @param[in] collection_size   collection size of field
 * @param[in] timestep          timestep at which put/get is called
 *                              for this field in ISO 8601 format
 * @return pointer to coupling field
 */
struct coupling_field * yac_instance_add_field(
  struct yac_instance * instance, char const * field_name,
  size_t component_index, struct yac_basic_grid * grid,
  struct interp_field * interp_fields, size_t num_interp_fields,
  int collection_size, char const * timestep);

/**
 * Defines a couple for a yac instance
 * @param[in] instance            yac instance
 * @param[in] src_comp_name       component name of the source component
 * @param[in] src_grid_name       grid name of the source grid
 * @param[in] src_field_name      field name of the source field
 * @param[in] tgt_comp_name       component name of the target component
 * @param[in] tgt_grid_name       grid name of the target grid
 * @param[in] tgt_field_name      field name of the target field
 * @param[in] coupling_period     time step for the coupling
 * @param[in] time_reduction      type for reducing multiple timesteps
 *                                (@see YAC_REDUCTION_TIME_NONE etc.)
 * @param[in] interp_stack_config interpolation stack config to be used
 * @param[in] src_lag             lag for this couple on the source component
 * @param[in] tgt_lag             lag for this couple on the target component
 * @param[in] weight_file_name    file name for the weights file.
 *                                `NULL` to disable this feature
 * @param[in] mapping_on_source   side where the mapping is computed.
 *                                Currently only source = 1 and target = 0 are allowed
 * @param[in] scale_factor        scale factor
 * @param[in] scale_summand       scale summand
 * @param[in] num_src_mask_names  number of source field mask names
 *                                ("0" if no source field mask names are
 *                                 provided)
 * @param[in] src_mask_names      array of source field mask names
 *                                ("NULL" if num_src_mask_names == 0)
 * @param[in] tgt_mask_name       target field mask name
 *                                ("NULL" if no target field mask name is
 *                                 provided)
 */
void yac_instance_def_couple(
  struct yac_instance * instance,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_period, int time_reduction,
  struct yac_interp_stack_config * interp_stack_config, int src_lag, int tgt_lag,
  const char* weight_file_name, int mapping_on_source,
  double scale_factor, double scale_summand, size_t num_src_mask_names,
  char const * const * src_mask_names, char const * tgt_mask_name);

/**
 * Gets a component defined in a yac instance
 * @param[in] instance          yac instance
 * @param[in] component_index   component index
 * @return component
 */
struct component * yac_instance_get_component(
  struct yac_instance * instance, size_t component_index);

/**
 * synchronizes the grid and field definitions
 @param[in] instance yac instance
*/
void yac_instance_sync_def(struct yac_instance * instance);

/**
 * initiates the generation of all data structures required for
 * interpolation
 * @param[in] instance yac instance
 */
void yac_instance_setup(struct yac_instance * instance);

/**
 * initiates the generation of all data structures required for
 * interpolation
 * @param[in] instance   yac instance
 * @param[in] emit_flags flags for configuring the generated coupling
 *                       configuration output
 *                       (\ref YAC_YAML_EMITTER_DEFAULT or
 *                        \ref YAC_YAML_EMITTER_JSON)
 * @return string containing coupling configuration
 */
char * yac_instance_setup_and_emit_config(
  struct yac_instance * instance, int emit_flags);

/**
 * returns a communicator containing all processes of the provided components
 *
 * @param[in]  instance   yac instance
 * @param[in]  comp_names component names
 * @param[in]  num_comp_names
 */
MPI_Comm yac_instance_get_comps_comm(
  struct yac_instance * instance,
  char const ** comp_names, size_t num_comp_names);

struct coupling_field* yac_instance_get_field(struct yac_instance * instance,
  const char * comp_name, const char* grid_name, const char * field_name);

#endif // INSTANCE_H
