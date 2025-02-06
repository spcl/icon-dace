// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils_mci.h"
#include "utils_common.h"
#include "yac_mpi_common.h"
#include "component.h"
#include "yac.h"

struct component_data {

  char const * name;
  MPI_Group group;
  int is_local;
};

struct yac_component_config {

  MPI_Comm comm;
  size_t num_comps;
  struct component_data comps[];
};

static int compare_component_data(void const * a, void const * b) {
  return strcmp(((struct component_data const *)a)->name,
                ((struct component_data const *)b)->name);
}

struct yac_component_config * yac_component_config_new(
  struct yac_couple_config * couple_config, char const ** names,
  size_t num_names, MPI_Comm comm_) {

  // get number of globally defined components
  size_t num_global_components =
    yac_couple_config_get_num_components(couple_config);

  // make a copy of the provided communicator to avoid interference
  // with other communication
  MPI_Comm comm;
  yac_mpi_call(MPI_Comm_dup(comm_, &comm), comm_);

  // allocate and initialise basic component configuration
  struct yac_component_config * comp_config =
    xmalloc(
      1 * sizeof(*comp_config) +
      num_global_components * sizeof(struct component_data));
  comp_config->comm = comm;
  comp_config->num_comps = num_global_components;
  for (size_t i = 0; i < num_global_components; ++i)
    comp_config->comps[i].name =
      strdup(yac_couple_config_get_component_name(couple_config, i));

  // sort global components by name --> ensure identical processing
  // order on all processes
  qsort(
    comp_config->comps, num_global_components,
    sizeof(*(comp_config->comps)), compare_component_data);

  // get the MPI group of the provided communicator
  MPI_Group world_group;
  yac_mpi_call(MPI_Comm_group(comm, &world_group), comm);

  // allocate temporary data
  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  int * rank_mask = xmalloc(2 * (size_t)comm_size * sizeof(*rank_mask));
  int * ranks = rank_mask + (size_t)comm_size;

  // for all sorted global component
  for (size_t i = 0; i < num_global_components; ++i) {

    struct component_data * comp = &(comp_config->comps[i]);
    char const * comp_name = comp->name;

    // search for current component in list of locally defined components
    int is_local = 0;
    for (size_t j = 0; (j < num_names) && !is_local; ++j)
      is_local = !strcmp(comp_name, names[j]);

    // determine which ranks are part of the current component
    yac_mpi_call(
      MPI_Allgather(
        &is_local, 1, MPI_INT, rank_mask, 1, MPI_INT, comm), comm);

    // generate list of all ranks included in the current component
    int rank_count = 0;
    for (int rank = 0; rank < comm_size; ++rank) {
      if (rank_mask[rank]) {
        ranks[rank_count] = rank;
        ++rank_count;
      }
    }

    // generate group containing all ranks from the current component
    MPI_Group comp_group;
    yac_mpi_call(
      MPI_Group_incl(
        world_group, rank_count, ranks, &comp_group), comm);

    comp->group = comp_group;
    comp->is_local = is_local;
  }

  // cleanup
  yac_mpi_call(MPI_Group_free(&world_group), comm);
  free(rank_mask);

  return comp_config;
}

static inline int compare_int(const void * a, const void * b) {

  int const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
}

static size_t get_comp_idx(
  char const * caller, struct yac_component_config * comp_config,
  char const * comp_name) {

  struct component_data * comps = comp_config->comps;
  size_t num_comps = comp_config->num_comps;
  size_t comp_idx = SIZE_MAX;
  for (size_t i = 0; (i < num_comps) && (comp_idx == SIZE_MAX); ++i)
    if (!strcmp(comps[i].name, comp_name)) comp_idx = i;
  YAC_ASSERT_F(
    comp_idx != SIZE_MAX,
    "ERROR(%s): invalid component name: \"%s\"", caller, comp_name)
  return comp_idx;
}

MPI_Comm yac_component_config_get_comps_comm(
  struct yac_component_config * comp_config,
  const char ** names, size_t num_names) {

  // if no component name was provided
  if (num_names == 0) return MPI_COMM_NULL;

  YAC_ASSERT_F(
    num_names < INT_MAX,
    "ERROR(yac_component_config_get_comps_comm): too many components (%zu)",
    num_names);

  // get the index of each component
  int * comp_idxs = xmalloc(num_names * sizeof(*comp_idxs));
  for (size_t i = 0; i < num_names; ++i)
    comp_idxs[i] =
      (int)get_comp_idx(
        "yac_component_config_get_comps_comm", comp_config, names[i]);

  // sort and remove duplicated component indices
  qsort(comp_idxs, num_names, sizeof(*comp_idxs), compare_int);
  yac_remove_duplicates_int(comp_idxs, &num_names);

  MPI_Group comps_group = MPI_GROUP_EMPTY;
  for (size_t i = 0; i < num_names; ++i) {
    int comp_idx = comp_idxs[i];
    MPI_Group comp_group = comp_config->comps[comp_idx].group;
    MPI_Group union_group;
    yac_mpi_call(
      MPI_Group_union(comps_group, comp_group, &union_group),
      comp_config->comm);
    if (comps_group != MPI_GROUP_EMPTY)
      yac_mpi_call(MPI_Group_free(&comps_group), comp_config->comm);
    comps_group = union_group;
  }

  int group_rank, group_size;
  yac_mpi_call(MPI_Group_rank(comps_group, &group_rank), comp_config->comm);
  yac_mpi_call(MPI_Group_size(comps_group, &group_size), comp_config->comm);
  YAC_ASSERT(
    group_rank != MPI_UNDEFINED,
    "ERROR(yac_component_config_get_comps_comm): "
    "local process not included in any component provided to this routine");

  // get rank (from comp_config->comm) of neighbouring ranks
  int group_neigh_ranks[3], neigh_ranks[3];
  group_neigh_ranks[0] = (group_rank + 1) % group_size;
  group_neigh_ranks[1] = (group_rank + group_size - 1) % group_size;
  group_neigh_ranks[2] = group_rank;
  MPI_Group comp_config_group;
  yac_mpi_call(
    MPI_Comm_group(comp_config->comm, &comp_config_group), comp_config->comm);
  yac_mpi_call(
    MPI_Group_translate_ranks(
      comps_group, 3, group_neigh_ranks, comp_config_group, neigh_ranks),
    comp_config->comm);
  MPI_Group_free(&comp_config_group);

  // exchange number of names with neighbouring processes
  int num_names_buffer = (int)num_names;
  int const tag = 0;
  yac_mpi_call(
    MPI_Sendrecv_replace(
      &num_names_buffer, 1, MPI_INT, neigh_ranks[0], tag,
      neigh_ranks[1], tag, comp_config->comm, MPI_STATUS_IGNORE),
    comp_config->comm);
  YAC_ASSERT_F(
    num_names_buffer == (int)num_names,
    "ERROR(yac_component_config_get_comps_comm): "
    "processes do not agree on number of component names "
    "(rank %d num_names %d != rank %d num_names %zu)",
    neigh_ranks[1], num_names_buffer, neigh_ranks[2], num_names);

  // exchange component indices with neighbouring processes
  int * comp_idxs_recv_buffer = xmalloc(num_names * sizeof(*comp_idxs_recv_buffer));
  yac_mpi_call(
    MPI_Sendrecv(
      comp_idxs, (int)num_names, MPI_INT, neigh_ranks[0], tag,
      comp_idxs_recv_buffer, (int)num_names, MPI_INT, neigh_ranks[1], tag,
      comp_config->comm, MPI_STATUS_IGNORE), comp_config->comm);
  for (size_t i = 0; i < num_names; ++i) {
    YAC_ASSERT_F(
      comp_idxs[i] == comp_idxs_recv_buffer[i],
      "ERROR(yac_component_config_get_comps_comm): "
      "processes do not agree on component indices "
      "(rank %d comp_idx[%zu] %d != rank %d comp_idx[%zu] %d)",
      neigh_ranks[1], i, comp_idxs[i],
      neigh_ranks[2], i, comp_idxs_recv_buffer[i]);
  }
  free(comp_idxs_recv_buffer);
  free(comp_idxs);

  MPI_Comm comps_comm;
  yac_mpi_call(
    MPI_Comm_create_group(
      comp_config->comm, comps_group, tag, &comps_comm), comp_config->comm);

  yac_mpi_call(MPI_Group_free(&comps_group), comp_config->comm);

  return comps_comm;
}

int yac_component_config_contains_component(
  struct yac_component_config * comp_config, char const * comp_name) {

  struct component_data * global_comps = comp_config->comps;
  size_t num_global_comps = comp_config->num_comps;

  for (size_t i = 0; i < num_global_comps; ++i)
    if (!strcmp(comp_name, global_comps[i].name))
      return global_comps[i].is_local;
  return 0;
}

static MPI_Group get_local_comp_group(
  char const * caller, struct yac_component_config * comp_config,
  char const * comp_name) {

  size_t comp_idx =
    get_comp_idx(caller, comp_config, comp_name);

  YAC_ASSERT_F(
    comp_config->comps[comp_idx].is_local,
    "ERROR(%s): component \"%s\" is defined but not on local process",
    caller, comp_name);

  return comp_config->comps[comp_idx].group;
}

int yac_component_config_comp_size(
  struct yac_component_config * comp_config, char const * comp_name){

  int size;
  MPI_Group_size(
    get_local_comp_group(
      "yac_component_config_comp_size", comp_config, comp_name), &size);
  return size;
}

int yac_component_config_comp_rank(
  struct yac_component_config * comp_config, char const * comp_name){

  int rank;
  MPI_Group_rank(
    get_local_comp_group(
      "yac_component_config_comp_rank", comp_config, comp_name), &rank);
  return rank;
}

void yac_component_config_delete(struct yac_component_config * comp_config) {

  if (comp_config == NULL) return;

  for (size_t i = 0; i < comp_config->num_comps; ++i) {
    free((void*)comp_config->comps[i].name);
    yac_mpi_call(
      MPI_Group_free(&(comp_config->comps[i].group)), comp_config->comm);
  }

  yac_mpi_call(MPI_Comm_free(&(comp_config->comm)), MPI_COMM_WORLD);
  free(comp_config);
}
