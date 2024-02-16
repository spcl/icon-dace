/**
 * @file test_proc_sphere_part_parallel.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
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
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

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
    size_t num_cells = local_count[comm_rank][0] * local_count[comm_rank][1];
    struct dist_cell * cells = xmalloc(num_cells * sizeof(*cells));

    for (size_t i = 0, k = 0; i < local_count[comm_rank][1]; ++i)
      for (size_t j = 0; j < local_count[comm_rank][0]; ++j, ++k)
        LLtoXYZ_deg(
          x_vertices[local_start[comm_rank][0]+j],
          y_vertices[local_start[comm_rank][1]+i], cells[k].coord);

    struct proc_sphere_part_node * proc_sphere_part =
      yac_redistribute_cells(&cells, &num_cells, MPI_COMM_WORLD);

    free(cells);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    double coords[5][2][3] =
      {{{1,0,0},{-1,0,0}},
       {{0,1,0}},
       {{0,-1,0}},
       {{0,0,1}},
       {{0,0,-1}}};
    size_t num_cells[5] = {2,1,1,1,1};
    struct dist_cell * cells = xmalloc(2 * sizeof(*cells));

    for (int i = 0; i < num_cells[comm_rank]; ++i)
      for (int j = 0; j < 3; ++j)
        cells[i].coord[j] = coords[comm_rank][i][j];

    struct proc_sphere_part_node * proc_sphere_part =
      yac_redistribute_cells(&cells, &num_cells[comm_rank], MPI_COMM_WORLD);

    free(cells);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    size_t num_cells = 0;

    struct proc_sphere_part_node * proc_sphere_part =
      yac_redistribute_cells(NULL, &num_cells, MPI_COMM_WORLD);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
