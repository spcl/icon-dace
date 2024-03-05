/**
 * @file component.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "yac_mpi.h"
#include "component.h"
#include "yac_interface.h"

// general component struct
struct component {

  // general component data

  char * name;
  MPI_Comm comm;
};

struct component_groups {
  char const * name;
  MPI_Group group;
};

struct yac_component_config {

  MPI_Comm comm;

  struct component ** components;
  size_t num_components;

  struct component_groups * comp_groups;
  size_t num_comp_groups;
};

static int compare_component_groups(void const * a, void const * b) {
  return strcmp(((struct component_groups const *)a)->name,
                ((struct component_groups const *)b)->name);
}

struct yac_component_config * yac_component_config_new(
  struct yac_couple_config * couple_config, char const ** names,
  size_t num_names, MPI_Comm comm) {

  size_t num_global_components =
    yac_couple_config_get_num_components(couple_config);

  struct yac_component_config * comp_config = xmalloc(1 * sizeof(*comp_config));

  yac_mpi_call(MPI_Comm_dup(comm, &(comp_config->comm)), comm);
  comm = comp_config->comm;

  // initialise components
  comp_config->components =
    xmalloc(num_names * sizeof(*(comp_config->components)));
  comp_config->num_components = num_names;
  for (size_t i = 0; i < num_names; ++i) {
    comp_config->components[i] =
      xmalloc(1 * sizeof(*(comp_config->components[i])));
    comp_config->components[i]->name = strdup(names[i]);
  }

  MPI_Group world_group;
  yac_mpi_call(MPI_Comm_group(comm, &world_group), comm);

  struct component_groups * comp_groups =
    ((comp_config->comp_groups =
        xmalloc(num_global_components * sizeof(*(comp_config->comp_groups)))));
  comp_config->num_comp_groups = num_global_components;

  // get the names of all components in couple_config
  for (size_t i = 0; i < num_global_components; ++i)
    comp_groups[i].name =
      strdup(yac_couple_config_get_component_name(couple_config, i));

  // sort comp_groups by name --> ensure identical processing order of
  // comp_groups on all processes
  qsort(comp_groups, num_global_components, sizeof(*comp_groups),
        compare_component_groups);

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  int * rank_mask = xmalloc((size_t)comm_size * sizeof(*rank_mask));
  int * ranks = xmalloc((size_t)comm_size * sizeof(*ranks));

  // for all sorted component groups
  for (size_t i = 0; i < num_global_components; ++i) {

    char const * comp_name = comp_groups[i].name;

    // search for current component in list of locally defined components
    struct component * comp = NULL;
    for (size_t j = 0; (j < num_names) && !comp; ++j)
      if (!strcmp(comp_name, names[j]))
        comp = comp_config->components[j];

    int is_in_comp = comp != NULL;

    // determine which rank is part of the current component
    yac_mpi_call(
      MPI_Allgather(
        &is_in_comp, 1, MPI_INT, rank_mask, 1, MPI_INT, comm), comm);

    // generate list of all ranks included in the current component
    int comp_size = 0;
    for (int rank = 0; rank < comm_size; ++rank) {
      if (rank_mask[rank]) {
        ranks[comp_size] = rank;
        ++comp_size;
      }
    }

    // generate group containing all ranks from the current component
    yac_mpi_call(
      MPI_Group_incl(
        world_group, comp_size, ranks, &(comp_groups[i].group)), comm);

    // generate communicator for all process that defined the current
    // component
    MPI_Comm comp_comm;
    yac_mpi_call(
      MPI_Comm_split(
        comm, is_in_comp, is_in_comp?0:MPI_UNDEFINED, &comp_comm), comm);
    if (is_in_comp) comp->comm = comp_comm;
  }

  free(ranks);
  free(rank_mask);
  yac_mpi_call(MPI_Group_free(&world_group), comm);

  return comp_config;
}

const char * yac_component_get_name(struct component * comp) {

  YAC_ASSERT(comp != NULL, "ERROR: invalid component pointer")

  return comp->name;
}

MPI_Comm yac_component_get_communicator(struct component * comp) {

  YAC_ASSERT(comp != NULL, "ERROR: invalid component pointer")

   return comp->comm;
}

static inline int compare_int(const void * a, const void * b) {

  int const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
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
  for (size_t i = 0, comp_idx; i < num_names; ++i) {
    comp_idx = SIZE_MAX;
    for (size_t j = 0;
         (j < comp_config->num_comp_groups) && (comp_idx == SIZE_MAX); ++j)
      if (!strcmp(comp_config->comp_groups[j].name, names[i])) comp_idx = j;
    YAC_ASSERT_F(
      comp_idx != SIZE_MAX,
      "ERROR(yac_component_config_get_comps_comm): "
      "invalid component name: \"%s\"", names[i]);
    comp_idxs[i] = (int)comp_idx;
  }

  // sort and remove duplicated component indices
  qsort(comp_idxs, num_names, sizeof(*comp_idxs), compare_int);
  yac_remove_duplicates_int(comp_idxs, &num_names);

  MPI_Group comps_group = MPI_GROUP_EMPTY;
  for (size_t i = 0; i < num_names; ++i) {
    int comp_idx = comp_idxs[i];
    MPI_Group comp_group = comp_config->comp_groups[comp_idx].group;
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

static void yac_component_delete(struct component * comp) {

  YAC_ASSERT(comp != NULL, "ERROR: invalid component pointer")

  free(comp->name);
  yac_mpi_call(MPI_Comm_free(&(comp->comm)), MPI_COMM_WORLD);
  free(comp);
}

struct component ** yac_component_config_get_components(
  struct yac_component_config * comp_config) {

  return comp_config->components;
}

size_t yac_component_config_get_num_components(
  struct yac_component_config * comp_config) {

  return comp_config->num_components;
}

void yac_component_config_delete(struct yac_component_config * comp_config) {

  if (comp_config == NULL) return;

  for (size_t i = 0; i < comp_config->num_components; ++i)
    yac_component_delete(comp_config->components[i]);
  free(comp_config->components);

  for (size_t i = 0; i < comp_config->num_comp_groups; ++i) {
    free((void*)comp_config->comp_groups[i].name);
    yac_mpi_call(
      MPI_Group_free(&(comp_config->comp_groups[i].group)), comp_config->comm);
  }
  free(comp_config->comp_groups);

  yac_mpi_call(MPI_Comm_free(&(comp_config->comm)), MPI_COMM_WORLD);
  free(comp_config);
}
