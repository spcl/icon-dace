// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "proc_sphere_part.h"
#include "yac_mpi.h"

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  YAC_ASSERT(comm_size == 5, "ERROR wrong number of processes (has to be 5)")

  {
    // one process has no data
    double x_vertices[18] = {0,20,40,60,80,
                             100,120,140,160,180,
                             200,220,240,260,280,
                             300,320,340};
    double y_vertices[9] = {-80,-60,-40,-20,0,20,40,60,80};
    size_t local_start[5][2] = {{0,0},{7,0},{0,0},{0,3},{7,3}};
    size_t local_count[5][2] = {{10,6},{10,6},{0,0},{10,6},{10,6}};
    size_t num_vertices =
      local_count[comm_rank][0] * local_count[comm_rank][1];
    struct dist_vertex * vertices = xmalloc(num_vertices * sizeof(*vertices));

    for (size_t i = 0, k = 0; i < local_count[comm_rank][1]; ++i)
      for (size_t j = 0; j < local_count[comm_rank][0]; ++j, ++k)
        LLtoXYZ_deg(
          x_vertices[local_start[comm_rank][0]+j],
          y_vertices[local_start[comm_rank][1]+i],
          vertices[k].coord);

    struct proc_sphere_part_node * proc_sphere_part =
      yac_redistribute_vertices(&vertices, &num_vertices, MPI_COMM_WORLD);

    free(vertices);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    double x_vertices[] = {  0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
                           100,110,120,130,140,150,160,170,180,190,
                           200,210,220,230,240,250,260,270,280,290,
                           300,310,320,330,340,350};
    double y_vertices[] =
      {-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90};
    size_t local_start[5][2] = {{0,0},{0,5},{12,5},{24,5},{0,14}};
    size_t local_count[5][2] = {{36,5},{12,9},{12,9},{12,9},{36,5}};
    size_t num_vertices =
      local_count[comm_rank][0] * local_count[comm_rank][1];
    struct dist_vertex * vertices = xmalloc(num_vertices * sizeof(*vertices));

    for (size_t i = 0, k = 0; i < local_count[comm_rank][1]; ++i)
      for (size_t j = 0; j < local_count[comm_rank][0]; ++j, ++k)
        LLtoXYZ_deg(
          x_vertices[local_start[comm_rank][0]+j],
          y_vertices[local_start[comm_rank][1]+i],
          vertices[k].coord);

    struct proc_sphere_part_node * proc_sphere_part =
      yac_redistribute_vertices(&vertices, &num_vertices, MPI_COMM_WORLD);

    free(vertices);

    double x_search = 45.0, y_search = 80.0;
    struct bounding_circle bnd_circle;
    LLtoXYZ_deg(x_search, y_search, bnd_circle.base_vector);
    bnd_circle.inc_angle.sin = sin(3.13);
    bnd_circle.inc_angle.cos = cos(3.13);

    int search_ranks[5], rank_count;
    yac_proc_sphere_part_do_bnd_circle_search(
      proc_sphere_part, bnd_circle, search_ranks, &rank_count);

    if (rank_count != 5)
      PUT_ERR("error in yac_proc_sphere_part_do_bnd_circle_search");

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    double coords[5][2][3] =
      {{{1,0,0},{-1,0,0}},
       {{0,1,0}},
       {{0,-1,0}},
       {{0,0,1}},
       {{0,0,-1}}};
    size_t num_vertices[5] = {2,1,1,1,1};
    struct dist_vertex * vertices = xmalloc(2 * sizeof(*vertices));

    for (size_t i = 0; i < num_vertices[comm_rank]; ++i)
      for (int j = 0; j < 3; ++j)
        vertices[i].coord[j] = coords[comm_rank][i][j];

    struct proc_sphere_part_node * proc_sphere_part =
      yac_redistribute_vertices(
        &vertices, &num_vertices[comm_rank], MPI_COMM_WORLD);

    free(vertices);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    size_t num_vertices = 0;

    struct proc_sphere_part_node * proc_sphere_part =
      yac_redistribute_vertices(NULL, &num_vertices, MPI_COMM_WORLD);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
