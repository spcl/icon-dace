/**
 * @file proc_sphere_part.h
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#ifndef PROC_SPHERE_PART_H
#define PROC_SPHERE_PART_H

#include "dist_grid.h"

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
struct dist_cell {
  double coord[3];
  enum NODE_TYPE node_type;
};

struct proc_sphere_part_node;

struct proc_sphere_part_node * yac_redistribute_cells(
  struct dist_cell ** cells, size_t * num_cells, MPI_Comm comm);
void yac_proc_sphere_part_node_delete(struct proc_sphere_part_node * node);
void yac_proc_sphere_part_do_point_search(
  struct proc_sphere_part_node * node, coordinate_pointer search_coords,
  size_t count, int * ranks);
void yac_proc_sphere_part_do_bnd_circle_search(
  struct proc_sphere_part_node * node, struct bounding_circle bnd_circle,
  int * ranks, int * rank_count);
void yac_proc_sphere_part_get_neigh_ranks(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes,
  uint64_t min_size, int * send_flags, int * recv_flags,
  int comm_rank, int comm_size);

#endif // PROC_SPHERE_PART_H
