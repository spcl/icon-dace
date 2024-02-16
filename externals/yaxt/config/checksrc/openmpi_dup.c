/**
 * @file openmpi_dup.c
 * @brief detects problem some versions of OpenMPI have with MPI_Type_dup
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>
#include <stdio.h>

/**
 * usage:

mpicc -std=gnu99 mpi_test.c && a.out

 * expected output: none

 * actual output on defective OpenMPI installations

recv_data[0]     = 131072
ref_recv_data[0] = 0
recv_data[1]     = -65536
ref_recv_data[1] = 2

*/
enum { recv_data_count = 2 };

void do_test(MPI_Datatype * sends, int ref_recv_data[recv_data_count]);

int main(void) {

  MPI_Init(NULL, NULL);

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size == 1) {

    {
      MPI_Datatype sends[2];
      sends[0] = MPI_INT;
      {
        int count = 1;
        int blocklength = 1;
        int array_of_displacements[] = {2};
        MPI_Type_create_indexed_block(count, blocklength, array_of_displacements,
                                      MPI_INT, &sends[1]);
        MPI_Type_commit(&sends[1]);
      }

      int ref_recv_data[recv_data_count] = {0, 2};

      do_test(sends, ref_recv_data);

      MPI_Type_free(&sends[1]);
    }

    {
      MPI_Datatype sends[2];

      // by definition of MPI_Type_dup, no MPI_Type_commit is required
      MPI_Type_dup(MPI_INT, &sends[0]);

      {
        int count = 1;
        int blocklength = 1;
        int array_of_displacements[] = {2};
        MPI_Type_create_indexed_block(count, blocklength, array_of_displacements,
                                      MPI_INT, &sends[1]);
        MPI_Type_commit(&sends[1]);
      }

      int ref_recv_data[recv_data_count] = {0, 2};

      do_test(sends, ref_recv_data);

      MPI_Type_free(&sends[1]);
      MPI_Type_free(&sends[0]);
    }

    {
      MPI_Datatype sends[2];

      {
        MPI_Type_dup(MPI_INT, &sends[0]);
        // by definition of MPI_Type_dup no MPI_Type_commit is required
      }
      {
        int count = 1;
        int blocklength = 1;
        int array_of_displacements[] = {1};
        MPI_Type_create_indexed_block(count, blocklength, array_of_displacements,
                                      MPI_INT, &sends[1]);
        MPI_Type_commit(&sends[1]);
      }

      int ref_recv_data[recv_data_count] = {0, 1};

      do_test(sends, ref_recv_data);

      MPI_Type_free(&sends[1]);
      MPI_Type_free(&sends[0]);
    }
  }

  MPI_Finalize();

  return 0;
}

void do_test(MPI_Datatype * sends, int ref_recv_data[recv_data_count]) {

  MPI_Datatype send_dt;

  {
    enum { count = 2 };
    int array_of_blocklengths[] = {1, 1};
    MPI_Aint array_of_displacements[] = {0, 0};
    MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
                           sends, &send_dt);
    MPI_Type_commit(&send_dt);
  }

  int send_data[3] = {0,1,2};
  int recv_data[2] = {-1,-1};

  MPI_Request request;

  MPI_Irecv(recv_data, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
  MPI_Send(send_data, 1, send_dt, 0, 0, MPI_COMM_WORLD);

  MPI_Waitall(1, &request, MPI_STATUSES_IGNORE);

  MPI_Type_free(&send_dt);

  for (int i = 0; i < recv_data_count; ++i) {

    if (recv_data[i] != ref_recv_data[i])
      printf("recv_data[%d]     = %d\n"
             "ref_recv_data[%d] = %d\n",
             i, recv_data[i], i, ref_recv_data[i]);
  }
}
