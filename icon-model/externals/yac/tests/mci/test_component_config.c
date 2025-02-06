// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "component.h"
#include "couple_config.h"

static struct yac_couple_config * generate_couple_config(
  char ** comp_names, size_t count);

int main (void) {

  // initialise MPI
  MPI_Init(NULL, NULL);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size != 8) {

    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }

  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------

  // names of all components
  char const * all_comp_names[] = {"comp_a", "comp_b", "comp_c"};
  enum {NUM_COMPS = sizeof(all_comp_names) / sizeof(all_comp_names[0])};
  // ranks for each component
  int comp_ranks[NUM_COMPS][4] = {{0, 1, 2, 4}, {0, 1, 3, 5}, {0, 2, 3, 6}};
  // number of ranks per component
  int num_comp_ranks[NUM_COMPS] = {4, 4, 4};

  // get the names of all components that are to be defined on the local process
  char const * local_comp_names[NUM_COMPS];
  size_t num_local_comps = 0;
  for (int i = 0; i < NUM_COMPS; ++i)
    for (int j = 0; j < num_comp_ranks[i]; ++j)
      if (rank == comp_ranks[i][j])
        local_comp_names[num_local_comps++] = all_comp_names[i];

  // generate dummy couple_config
  struct yac_couple_config * couple_config =
    generate_couple_config((char**)all_comp_names, NUM_COMPS);

  // generate component configuration
  struct yac_component_config * comp_config =
    yac_component_config_new(
      couple_config, (char const **)local_comp_names, num_local_comps,
      MPI_COMM_WORLD);

  // free memory in dummy couple_config
  yac_couple_config_delete(couple_config);

  //----------------------------------------------------------------------------
  // generate reference data
  //----------------------------------------------------------------------------

  // generate reference communicators for each component and for each
  // component pair
  MPI_Comm ref_comp_comm[NUM_COMPS], ref_comp_pair_comm[NUM_COMPS], ref_all_comps_comm;
  {
    MPI_Group world_group, comp_group[NUM_COMPS];
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    for (int i = 0; i < NUM_COMPS; ++i) {
      MPI_Group_incl(
        world_group, num_comp_ranks[i], comp_ranks[i], &comp_group[i]);
      MPI_Comm_create(MPI_COMM_WORLD, comp_group[i], &ref_comp_comm[i]);
    }
    for (int i = 0, k = 0; i < NUM_COMPS; ++i) {
      for (int j = i + 1; j < NUM_COMPS; ++j, ++k) {
        MPI_Group comp_pair_group;
        MPI_Group_union(comp_group[i], comp_group[j], &comp_pair_group);
        MPI_Comm_create(
          MPI_COMM_WORLD, comp_pair_group, &ref_comp_pair_comm[k]);
        MPI_Group_free(&comp_pair_group);
      }
    }
    {
      MPI_Group all_comps_group = MPI_GROUP_EMPTY;
      for (int i = 0; i < NUM_COMPS; ++i) {
        MPI_Group merge_group;
        MPI_Group_union(all_comps_group, comp_group[i], &merge_group);
        if (all_comps_group != MPI_GROUP_EMPTY)
          MPI_Group_free(&all_comps_group);
        all_comps_group = merge_group;
      }
      MPI_Comm_create(
        MPI_COMM_WORLD, all_comps_group, &ref_all_comps_comm);
      MPI_Group_free(&all_comps_group);
    }
    for (int i = 0; i < NUM_COMPS; ++i)
      MPI_Group_free(&comp_group[i]);
    MPI_Group_free(&world_group);
  }

  //----------------------------------------------------------------------------
  // testing
  //----------------------------------------------------------------------------

  for (int i = 0; i < NUM_COMPS; ++i) {
    int is_local = 0;
    int rank_idx;
    for (rank_idx = 0; (rank_idx < num_comp_ranks[i]) && !is_local; ++rank_idx)
      if (comp_ranks[i][rank_idx] == rank) is_local = 1;
    if (yac_component_config_contains_component(
          comp_config, all_comp_names[i]) != is_local)
      PUT_ERR("error in yac_component_config_contains_component");
    if (is_local) {
      if (yac_component_config_comp_size(comp_config, all_comp_names[i]) !=
          num_comp_ranks[i])
        PUT_ERR("error in yac_component_config_comp_size");
      if (yac_component_config_comp_rank(comp_config, all_comp_names[i]) !=
          rank_idx - 1)
        PUT_ERR("error in yac_component_config_comp_rank");
    }
  }
  if (yac_component_config_contains_component(comp_config, "dummy"))
    PUT_ERR("error in yac_component_config_contains_component");

  // encode into a single int which components are defined locally
  int comp_flags = 0;
  for (size_t i = 0; i < num_local_comps; ++i)
    comp_flags |= (1 << (local_comp_names[i][5] - 'a'));

  // check component pair communicators
  for (int i = 0, k = 0; i < NUM_COMPS; ++i) { // for all component
    for (int j = i + 1; j < NUM_COMPS; ++j, ++k) { // for all remaining components
      if ((comp_flags & (1 << i)) || (comp_flags & (1 << j))) {
        MPI_Comm comp_pair_comm =
          yac_component_config_get_comps_comm(
            comp_config,
            (char const *[]){all_comp_names[i], all_comp_names[j]}, 2);
        int compare_result;
        MPI_Comm_compare(
          ref_comp_pair_comm[k], comp_pair_comm, &compare_result);
        if (compare_result != MPI_CONGRUENT)
          PUT_ERR("error in yac_component_config_get_comps_comm");
        MPI_Comm_free(&comp_pair_comm);
      }
    }
  }

  // check communicators containing all processes
  if (ref_all_comps_comm != MPI_COMM_NULL) {
    MPI_Comm all_comps_comm =
      yac_component_config_get_comps_comm(comp_config, all_comp_names, NUM_COMPS);
    int compare_result;
    MPI_Comm_compare(ref_all_comps_comm, all_comps_comm, &compare_result);
    if (compare_result != MPI_CONGRUENT)
      PUT_ERR("error in yac_component_config_get_comps_comm");
    MPI_Comm_free(&all_comps_comm);
  }

  // check communicator for empty component list
  if (yac_component_config_get_comps_comm(
        comp_config, NULL, 0) != MPI_COMM_NULL)
    PUT_ERR("error in yac_component_config_get_comps_comm");

  //----------------------------------------------------------------------------
  // clean-up
  //----------------------------------------------------------------------------

  yac_component_config_delete(comp_config);

  if (ref_all_comps_comm != MPI_COMM_NULL)
    MPI_Comm_free(&ref_all_comps_comm);
  for (int i = 0; i < NUM_COMPS; ++i) {
    if (ref_comp_comm[i] != MPI_COMM_NULL)
      MPI_Comm_free(&ref_comp_comm[i]);
    if (ref_comp_pair_comm[i] != MPI_COMM_NULL)
      MPI_Comm_free(&ref_comp_pair_comm[i]);
  }

  // finalize MPI
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static struct yac_couple_config * generate_couple_config(
  char ** comp_names, size_t count) {

  struct yac_couple_config * couple_config = yac_couple_config_new();
  for (size_t i = 0; i < count; ++i)
    yac_couple_config_add_component(couple_config, comp_names[i]);
  return couple_config;
}
