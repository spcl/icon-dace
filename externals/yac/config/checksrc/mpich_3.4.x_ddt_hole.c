#include <mpi.h>

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if (__GNUC__ == 11 || __GNUC__ == 12) && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

// acx_mpi_job_count=1

struct datatype {
  int i;
  uint64_t ui64;
};

static MPI_Datatype gen_ddt(MPI_Comm comm) {

  struct datatype dummy;
  MPI_Datatype ddt;
  MPI_Type_create_struct(
    2, (int[]){1,1},
    (MPI_Aint[])
      {(MPI_Aint)(intptr_t)(const void *)&(dummy.i) -
         (MPI_Aint)(intptr_t)(const void *)&dummy,
       (MPI_Aint)(intptr_t)(const void *)&(dummy.ui64) -
         (MPI_Aint)(intptr_t)(const void *)&dummy},
    (MPI_Datatype[]){MPI_INT, MPI_UINT64_T}, &ddt);
  MPI_Type_commit(&ddt);
  return ddt;
}

int main (void) {

  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Datatype ddt = gen_ddt(MPI_COMM_WORLD);

  int count = 20705;
  // in my tests for counts < 17403 the bug does not occure
  // int count = 17403;
  size_t buffer_size = count * sizeof(struct datatype);
  struct datatype * send_buffer = malloc(buffer_size);
  struct datatype * recv_buffer = malloc(buffer_size);

  memset(send_buffer, 0, buffer_size);
  memset(recv_buffer, 0, buffer_size);

  int tag = 0;
  MPI_Sendrecv(send_buffer, count, ddt, rank, tag,
               recv_buffer, count, ddt, rank, tag, 
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  free(recv_buffer);
  free(send_buffer);

  MPI_Type_free(&ddt);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
