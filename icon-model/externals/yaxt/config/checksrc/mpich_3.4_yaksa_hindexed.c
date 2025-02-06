/**
 *
 * This program demonstrates a bug that occurs when using MPI datatypes.
 * The following MPI implementations exhibit this bug:
 * mpich/3.4.1 3.4 and 3.4.2
 *
 * This problem does not occur with the following MPI implementations:
 * intelmpi/2018.4.274 OpenMPI 4.0.5
 *
 * usage:
 *
 *   mpicc mpich_datatype.c && a.out
 *
 * expected output for defective MPICH versions:
 *
 *   pack_buffer[2]:  3 !=  2
 *   pack_buffer[3]: -2 !=  3
 *   pack_buffer[4]:  5 !=  3
 *   pack_buffer[5]:  6 !=  5
 *   pack_buffer[6]:  7 !=  6
 *   pack_buffer[7]: -3 !=  7
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

int main(int argc, char **argv) {

  // init mpi
  MPI_Init(&argc, &argv);

  // generate mpi datatype for packing
  MPI_Datatype pack_dt;
  {
    enum { count = 5 };
    static const int array_of_blocklengths[count] = {2, 2, 1, 3, 3};
    static const int array_of_displacements[count] = {1, 2, 3, 5, 9};
    MPI_Datatype oldtype = MPI_DOUBLE;
    MPI_Type_indexed(
      count, (int*)array_of_blocklengths, (int*)array_of_displacements,
      oldtype, &pack_dt);
    MPI_Type_commit(&pack_dt);
  }
  static const double inbuf[] = {-1,  1,  2,  3,
                                 -2,  5,  6,  7,
                                 -3,  9, 10, 11};

  static const double ref_pack_buffer[]
    = {1, 2,  2, 3,  3,  5, 6, 7,  9, 10, 11};
  enum {pack_count = sizeof(ref_pack_buffer) / sizeof(ref_pack_buffer[0])};
  int pack_size;
  MPI_Pack_size(1, pack_dt, MPI_COMM_WORLD, &pack_size);
  void *pack_buffer = malloc((size_t)pack_size);
  assert(pack_buffer);

  // use generated mpi datatype to pack data
  int pos = 0;
  MPI_Pack((void *)inbuf, 1, pack_dt, pack_buffer, pack_size,
           &pos, MPI_COMM_WORLD);

  // check contents of pack_buffer
  int err_flag = 0;
  for (int i = 0, pos = 0; i < pack_count; ++i) {
    double curr_pack_buffer_content;
    MPI_Unpack(pack_buffer, pack_size, &pos,
               &curr_pack_buffer_content, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    if (curr_pack_buffer_content != ref_pack_buffer[i]) {
      err_flag = 1;
      fprintf(stderr, "pack_buffer[%d]: %2.0lf != %2.0lf\n",
              i, curr_pack_buffer_content, ref_pack_buffer[i]);
    }
  }

  // clean up
  free(pack_buffer);
  MPI_Type_free(&pack_dt);

  // finalise mpi
  MPI_Finalize();

  return (err_flag)?(EXIT_FAILURE):(EXIT_SUCCESS);
}
