/**
 * @file mpich_neighbor_alltoallw.c
 * @brief test for a defect in MPICH releases 4.0, 4.0.1 and 4.0.2
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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

/* run-time settings/expectations for configure-time checks
 * acx_mpirun_num_tasks=3
 * acx_mpi_expected_version >= 3.0
 */

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int my_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_size != 3) {
    fputs("ERROR: wrong number of processes (has to be 3)\n", stderr);
    exit(EXIT_FAILURE);
  }

  // build dist graph where rank 0 sends to rank 1 and 2

  static const int indegree[3] = {0,1,1};
  static const int sources[3][1] = {{-1},{0},{0}};
  static const int outdegree[3] = {2,0,0};
  static const int destinations[3][2] = {{1,2},{-1,-1},{-1,-1}};
  enum { reorder = 0 };
  MPI_Comm neigh_comm;

  MPI_Dist_graph_create_adjacent(
    MPI_COMM_WORLD, indegree[my_rank], sources[my_rank],  MPI_UNWEIGHTED,
    outdegree[my_rank], destinations[my_rank],  MPI_UNWEIGHTED,
    MPI_INFO_NULL, reorder, &neigh_comm);

  // send single integer from rank 0 to rank 1 and 2

  static const int send_buf[3][1] = {0, -1, -2};
  static const int sendcounts[3][2] = {{1,1},{-1,-1},{-1,-1}};
  static const MPI_Aint sdispls[3][2] = {{0,0},{-1,-1},{-1,-1}};
  static const MPI_Datatype sendtypes[3][2] = {{MPI_INT, MPI_INT},
                                  {MPI_DATATYPE_NULL, MPI_DATATYPE_NULL},
                                  {MPI_DATATYPE_NULL, MPI_DATATYPE_NULL}};
  int recv_buf[3][1] = {{-3}, {-4}, {-5}};
  static const int recvcounts[3][1] = {{-1},{1},{1}};
  static const MPI_Aint rdispls[3][1] = {{-1},{0},{0}};
  static const MPI_Datatype recvtypes[3][2] = {{MPI_DATATYPE_NULL},{MPI_INT},{MPI_INT}};

  MPI_Neighbor_alltoallw(
    (const void*)send_buf[my_rank], sendcounts[my_rank],
    sdispls[my_rank], sendtypes[my_rank],
    (void*)(recv_buf[my_rank]), recvcounts[my_rank],
    rdispls[my_rank], recvtypes[my_rank], neigh_comm);

  // check receive buffer

  static const int ref_recv_buf[3][1] = {{-3},{0},{0}};

  if (recv_buf[my_rank][0] != ref_recv_buf[my_rank][0]) {
    fputs("ERROR: wrong recv_buf\n", stderr);
    fprintf(stderr, "rank: %d, expected: %d, actual value: %d\n", my_rank,
            ref_recv_buf[my_rank][0], recv_buf[my_rank][0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * license-markup: "doxygen"
 * End:
 */
