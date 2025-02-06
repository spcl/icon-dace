/**
 * @file openmpi_datatype.c
 * @brief demonstrates a problem some OpenMPI versions have with
 * transferring some data layouts
 *
 * @copyright Copyright  (C)  2018 Moritz Hanke <hanke@dkrz.de>
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * usage:
mpicc openmpi_datatype.c && a.out

 * expected output for defective OpenMPI versions:
pack_buffer[1] = 1 != src_data[0] = 0
pack_buffer[2] = 2 != src_data[0] = 0
dst_data[1] = 1 != ref_dst_data[1] = 0
dst_data[2] = 2 != ref_dst_data[2] = 0

* it seems like there is a problem with send_dt
*/

int main(void) {

  size_t i;

  MPI_Init(NULL, NULL);

  MPI_Datatype send_dt;
  {
    MPI_Datatype base_types[2];
    MPI_Type_vector(3, 1, 0, MPI_INT, &base_types[0]);
    MPI_Type_indexed(1, (int[]){1}, (int[]){0}, MPI_INT, &base_types[1]);

    MPI_Type_create_struct(
      2, (int[]){1,1}, (MPI_Aint[]){0,12}, base_types, &send_dt);
    for (i = 0; i < sizeof(base_types) / sizeof(base_types[0]); ++i)
      MPI_Type_free(&base_types[i]);
    MPI_Type_commit(&send_dt);
  }

  enum {data_size = 4};

  const int send_transfer_pos[] = {0,0,0,3};

  // test exchange
  int src_data[data_size];
  int dst_data[data_size];
  int pack_buffer[data_size];
  int ref_dst_data[data_size];
  for (i = 0; i < data_size; ++i) src_data[i] = (int)i;
  for (i = 0; i < data_size; ++i) dst_data[i] = -1;
  for (i = 0; i < data_size; ++i) pack_buffer[i] = -1;
  for (i = 0; i < data_size; ++i) ref_dst_data[i] = -1;
  for (i = 0; i < sizeof(send_transfer_pos)/sizeof(send_transfer_pos[0]); ++i)
    ref_dst_data[i] = src_data[send_transfer_pos[i]];

   int pos = 0;
   MPI_Pack(src_data, 1, send_dt, pack_buffer, sizeof(pack_buffer), &pos,
            MPI_COMM_WORLD);

  MPI_Sendrecv(src_data, 1, send_dt, 0, 0,
               dst_data, 4, MPI_INT, 0, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  int err_count = 0;

  for (i = 0; i < sizeof(send_transfer_pos)/sizeof(send_transfer_pos[0]); ++i) {
    // interpretation of the pack_buffer is not covered by the MPI standard
    // but should work nevertheless
    if (pack_buffer[i] != src_data[send_transfer_pos[i]]) {
      fprintf(stderr, "pack_buffer[%zu] = %d != src_data[%d] = %d\n",
              i, pack_buffer[i], send_transfer_pos[i], ref_dst_data[i]);
      err_count++;
    }
  }

  for (i = 0; i < data_size; ++i) {
    if (dst_data[i] != ref_dst_data[i]) {
      fprintf(stderr, "dst_data[%zu] = %d != ref_dst_data[%zu] = %d\n",
              i, dst_data[i], i, ref_dst_data[i]);
      err_count++;
    }
  }

  // clean up
  MPI_Type_free(&send_dt);

  MPI_Finalize();

  return ((!err_count)?EXIT_SUCCESS:EXIT_FAILURE);
}