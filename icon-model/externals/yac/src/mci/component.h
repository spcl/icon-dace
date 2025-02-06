// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMPONENT_H
#define COMPONENT_H

#include <mpi.h>

#include "couple_config.h"

/** \example test_component_config.c
* This contains a test of the component_config interface routines.
*/

struct yac_component_config;

struct yac_component_config * yac_component_config_new(
  struct yac_couple_config * couple_config, char const ** names,
  size_t num_names, MPI_Comm comm);
void yac_component_config_delete(struct yac_component_config * comp_config);

/**
 * Checks whether a component is locally available
 * @param[in] comp_config component configuration
 * @param[in] comp_name   component name
 * @return "1" if comp_name was provided to \ref yac_component_config_new \n
 *         "0" otherwise
 */
int yac_component_config_contains_component(
  struct yac_component_config * comp_config, char const * comp_name);

/**
 * creates a communicator that contains the processes of all
 * listed components
 * @param[in] comp_config component configuration
 * @param[in] names       list of component names
 * @param[in] num_names   number of component names
 * @return components communicator
 */
MPI_Comm yac_component_config_get_comps_comm(
  struct yac_component_config * comp_config,
  const char ** names, size_t num_names);

int yac_component_config_comp_size(
  struct yac_component_config * comp_config, char const * comp_name);

int yac_component_config_comp_rank(
  struct yac_component_config * comp_config, char const * comp_name);

#endif
