/**
 * @file intelmpi_datatype.c
 * @brief demonstrates a problem some IntelMPI and MVAPICH2 versions have with
 * transferring some data layouts
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
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
/**
 *
 * This program demonstrates a bug that occurs when using MPI datatypes.
 * The following MPI implementations exhibit this bug:
 * intelmpi/4.1.3.049
 * intelmpi/5.0.3.048
 * intelmpi/5.1.0.038
 * intelmpi/5.1.1.109
 * intelmpi/5.1.2.150
 * intelmpi/5.1.3.223
 * intelmpi/2017.0.098
 * intelmpi/2017.1.132
 * mvapich2/1.9b
 * mvapich2/2.1
 *
 * This problem does not occur with the following MPI implementations:
 * intelmpi/2017.2.174
 * intelmpi/2017.3.196
 * intelmpi/2018.0.128
 * intelmpi/2018.1.163
 * intelmpi/2018.4.274
 * intelmpi/2018.5.288
 * openmpi/1.8.4
 * openmpi/1.10.1
 * openmpi/1.10.2
 * openmpi/2.0.2
 *
 * Call the program by running:
 *   mpirun -n 2 test_mpich
 *
 * In case of a bug the program delivers an output similar to the following:
 *   i        0 ref_data        0 recv_data        1
 *   i        1 ref_data        0 recv_data        2
 *   i        2 ref_data        1 recv_data        3
 *   i        3 ref_data        2 recv_data        4
 *   number of mismatches 4
 *
 * Description of the underlying problem:
 * On rank 0 we have a 1D data array of size [NBLK * NIDX * VECTOR_SIZE * NLEV].
 * If stored in a 3D array, the layout would be [NBLK][NLEV][NIDX][VECTOR_SIZE].
 * The variables count, displs, lenghts describe points accessing the vectors of
 * single level. Values for each of these points for each level is sent from
 * rank 1 to rank 0.
 * Each point consists of VECTOR_SIZE elements. The data sent by rank 1 is
 * level-wise; data from different levels does not overlap in the send buffer.
 *
 * use 2 mpi processes for the test
 * acx_mpirun_num_tasks=2
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define VECTOR_SIZE (2)
#define NLEV (90)
#define NIDX (8)
#define NBLK (7)
#define TAG (0)

int count = 9;
int displs[9] = {   1, // idx 1 blk 0
                  720, // idx 0 blk 1
                 1440, // idx 0 blk 2
                 1447, // idx 7 blk 2
                 2160, // idx 0 blk 3
                 2166, // idx 6 blk 3
                 2880, // idx 0 blk 4
                 3600, // idx 0 blk 5
                 4320  // idx 0 blk 6
                };
int lenghts[9] = {7,
                  8,
                  1,
                  1,
                  5,
                  2,
                  8,
                  6,
                  8};

int main(void) {

  int mpi_size, mpi_rank, num_mismatches = 0;

  // initialize MPI
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_size != 2) {

    if (mpi_rank == 0) fputs("Wrong number of mpi procs\n", stderr);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // compute the number of points that will be exchanged
  // (not including VECTOR_SIZE and NLEV)
  unsigned num_exchange_indices = 0;
  for (int i = 0; i < count; ++i) num_exchange_indices += lenghts[i];

  if (mpi_rank == 0) {

    // create contiguous datatype for the vector
    MPI_Datatype datatype_vector;
    MPI_Type_contiguous(VECTOR_SIZE, MPI_DOUBLE, &datatype_vector);

    // generate datatype for all points of a single level
    MPI_Datatype datatype_1lev;
    MPI_Type_indexed(count, lenghts, displs, datatype_vector, &datatype_1lev);

    // create a datatype for receiving all levels at once (each level is
    // VECTOR_SIZE * NIDX elements appart)
    MPI_Aint lb, extent, dummy;
    MPI_Type_get_extent(datatype_1lev, &lb, &dummy);
    MPI_Type_get_extent(MPI_DOUBLE, &dummy, &extent);
    extent *= VECTOR_SIZE * NIDX;
    MPI_Datatype datatype_with_extent;
    MPI_Type_create_resized(datatype_1lev, lb, extent, &datatype_with_extent);
    MPI_Datatype datatype_nlev;
    MPI_Type_contiguous(NLEV, datatype_with_extent, &datatype_nlev);
    MPI_Type_commit(&datatype_nlev);

    // generate reference data
    double * ref_data =
      malloc(NIDX * NBLK * NLEV * VECTOR_SIZE * sizeof(*ref_data));
    for (int i = 0; i < VECTOR_SIZE * NLEV * NIDX * NBLK; ++i) ref_data[i] = -1.0;
    for (int lev = 0, v = 1; lev < NLEV; ++lev)
      for (int i = 0; i < count; ++i)
        for (int j = 0; j < lenghts[i]; ++j)
          for (int vector = 0; vector < VECTOR_SIZE; ++vector, ++v)
            ref_data[
              vector + VECTOR_SIZE * (displs[i] + j) + lev * VECTOR_SIZE * NIDX] =
                (double)(v);

    // allocate and initialize receive data
    double * recv_data =
      malloc(NIDX * NBLK * NLEV * VECTOR_SIZE * sizeof(*recv_data));
    for (int i = 0; i < VECTOR_SIZE * NLEV * NIDX * NBLK; ++i)
        recv_data[i] = -1.0;


    // receive data from rank 1
    MPI_Recv(recv_data, 1, datatype_nlev, 1, TAG, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    // // unpack data
    // // int pos = 0;
    // double * send_data =
      // malloc(VECTOR_SIZE * NLEV * num_exchange_indices * sizeof(*send_data));
    // for (int i = 0; i < VECTOR_SIZE * NLEV * num_exchange_indices; ++i)
      // send_data[i] = (double)(i+1);
    // MPI_Unpack(send_data, VECTOR_SIZE * NLEV * num_exchange_indices * sizeof(double),
               // &pos, recv_data, 1, datatype_nlev, MPI_COMM_WORLD);

    // check data
    for (int i = 0; i < NIDX * NBLK * NLEV * VECTOR_SIZE; ++i)
      if (ref_data[i] != recv_data[i]) {
        printf("i %8d ref_data %8.0lf recv_data %8.0lf\n", i, ref_data[i],
               recv_data[i]);
        ++num_mismatches;
      }
    // if (num_mismatches)
      printf("number of mismatches %d\n", num_mismatches);

    // free generated datatypes
    // free(send_data);
    free(ref_data);
    free(recv_data);
    MPI_Type_free(&datatype_nlev);
    MPI_Type_free(&datatype_with_extent);
    MPI_Type_free(&datatype_1lev);
    MPI_Type_free(&datatype_vector);

  } else {

    // send requested data
    double * send_data =
      malloc(VECTOR_SIZE * NLEV * num_exchange_indices * sizeof(*send_data));
    for (int i = 0; i < VECTOR_SIZE * NLEV * num_exchange_indices; ++i)
      send_data[i] = (double)(i+1);
    MPI_Send(send_data, VECTOR_SIZE * NLEV * num_exchange_indices, MPI_DOUBLE,
             0, TAG, MPI_COMM_WORLD);
    free(send_data);
  }

  // finalize down MPI
  MPI_Finalize();

  return (num_mismatches)?EXIT_FAILURE:EXIT_SUCCESS;
}
