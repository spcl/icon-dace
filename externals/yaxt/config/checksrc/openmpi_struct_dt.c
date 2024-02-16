/**
 * @file openmpi_struct_dt.c
 * @brief demonstrates a problem some OpenMPI versions have with
 * transferring some data layouts
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
 * expected output for defective OpenMPI versions:
results_2[6]     = 8
ref_results_2[6] = 12
results_2[7]     = 9
ref_results_2[7] = 13

*/

static int
do_test(MPI_Datatype * recvs, MPI_Datatype * sends, int * inputs[2]);

int main(void) {

  MPI_Init(NULL, NULL);

  int rank, size;
  int ierror = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Datatype sends[2], recvs[2];

  {
    int count = 2, blocklen = 2, stride = 4;
    MPI_Type_vector(count, blocklen, stride, MPI_INT, &recvs[0]);
    MPI_Type_commit(&recvs[0]);
  }

  {
    int count = 1;
    int blocklength = 4;
    int array_of_displacements[] = {4};
    MPI_Type_create_indexed_block(count, blocklength, array_of_displacements,
                                  MPI_INT, &sends[0]);
    MPI_Type_commit(&sends[0]);
  }

  {
    int count = 1;
    int blocklength = 4;
    int array_of_displacements[] = {4};
    MPI_Type_create_indexed_block(count, blocklength, array_of_displacements,
                                  MPI_INT, &recvs[1]);
    MPI_Type_commit(&recvs[1]);
  }

  {
    int count = 2, blocklen = 2, stride = 4;
    MPI_Type_vector(count, blocklen, stride, MPI_INT, &sends[1]);
    MPI_Type_commit(&sends[1]);
  }

  {
    int raw_input[16] = {0,1,2,3,4,5,6,7,
                         8,9,10,11,12,13,14,15};
    int * input_1 = &raw_input[0], * input_2 = &raw_input[8];
    int * inputs[2] = {input_1, input_2};

    ierror = do_test(recvs, sends, inputs);
  }

  MPI_Type_free(&sends[1]);
  MPI_Type_free(&recvs[1]);
  MPI_Type_free(&sends[0]);
  MPI_Type_free(&recvs[0]);

  MPI_Finalize();

  return ierror;
}

static int
do_test(MPI_Datatype * recvs, MPI_Datatype * sends, int * inputs[2])
{

  int results_[3][8] = { {-1,-1,-1,-1,-1,-1,-1,-1},
                        {-1,-1,-1,-1,-1,-1,-1,-1},
                        {-1,-1,-1,-1,-1,-1,-1,-1} } ;
  int *results[2] = {results_[2], results_[0]};

  MPI_Datatype send_dt, recv_dt;

  {
    int count = 2;
    int array_of_blocklengths[] = {1, 1};
    MPI_Aint array_of_displacements[] = {0, ((unsigned char *)results[1] - (unsigned char *)results[0])};
    MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
                           recvs, &recv_dt);
    MPI_Type_commit(&recv_dt);
  }

  {
    int count = 2;
    int array_of_blocklengths[] = {1, 1};
    MPI_Aint array_of_displacements[] = {0, (inputs[1] - inputs[0]) * sizeof(int)};
    MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
                           sends, &send_dt);
    MPI_Type_commit(&send_dt);
  }

  MPI_Status status;

  MPI_Sendrecv(inputs[0], 1, send_dt, 0, 0,
               results[0], 1, recv_dt, 0, 0, MPI_COMM_SELF, MPI_STATUS_IGNORE);
  MPI_Type_free(&send_dt);
  MPI_Type_free(&recv_dt);

  static const int ref_results[2][8] = { {4,5,-1,-1,6,7,-1,-1},
                                         {-1,-1,-1,-1,8,9,12,13} };
  int retcode = 0;

  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 8; ++i) {

      retcode |= (results[j][i] != ref_results[j][i]);
      if (results[j][i] != ref_results[j][i])
        fprintf(stderr, "   results[%d][%d]     = %d\n"
                "!= ref_results[%d][%d] = %d\n",
                j, i, results[j][i], j, i, ref_results[j][i]);
    }
  return retcode;
}
