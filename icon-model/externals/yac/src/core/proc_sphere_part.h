// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PROC_SPHERE_PART_H
#define PROC_SPHERE_PART_H

#include "yac_types.h"
#include "geometry.h"

/** \example test_point_sphere_part.c
* This contains a test of the proc_sphere_part grid search algorithm.
*/

/** \example test_proc_sphere_part_parallel.c
* This contains a test of the proc_sphere_part grid search algorithm.
*/

enum NODE_TYPE {
  U_NODE = 1,
  T_NODE = 2,
};

// WARNING: before changing this datatype ensure that the MPI datatype created
// for this still matches its data layout
struct dist_vertex {
  double coord[3];
  enum NODE_TYPE node_type;
};

struct proc_sphere_part_node;

struct proc_sphere_part_node * yac_redistribute_vertices(
  struct dist_vertex ** vertices, size_t * num_vertices, MPI_Comm comm);
void yac_proc_sphere_part_node_delete(struct proc_sphere_part_node * node);
void yac_proc_sphere_part_do_point_search(
  struct proc_sphere_part_node * node, yac_coordinate_pointer search_coords,
  size_t count, int * ranks);
void yac_proc_sphere_part_do_bnd_circle_search(
  struct proc_sphere_part_node * node, struct bounding_circle bnd_circle,
  int * ranks, int * rank_count);
void yac_proc_sphere_part_get_neigh_ranks(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes,
  uint64_t min_size, int * send_flags, int * recv_flags,
  int comm_rank, int comm_size);

#endif // PROC_SPHERE_PART_H
