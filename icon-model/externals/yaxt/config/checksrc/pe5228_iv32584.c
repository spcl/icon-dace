#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(void) {

  MPI_Init(NULL, NULL);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Datatype type;
  int displacement = 1;

  MPI_Type_create_indexed_block(1, 1, &displacement, MPI_INT, &type);
  MPI_Type_commit(&type);

  int in_data[3] = {0,1,2};
  int out_data = -1;
  MPI_Status status;

  MPI_Sendrecv(in_data, 1, type, my_rank, 0,
		       &out_data, 1, MPI_INT, my_rank, 0,
               MPI_COMM_WORLD, &status);

  if (out_data != 1) {
      printf("error in MPI_Sendrecv "
    		 "(out_data is %d but should be 1)\n", out_data);
      return EXIT_FAILURE;
  }

  MPI_Request request;

  out_data = -1;

  MPI_Irecv(&out_data, 1, MPI_INT, my_rank, 0, MPI_COMM_WORLD, &request);
  MPI_Send(in_data, 1, type, my_rank, 0, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);

  if (out_data != 1) {
      printf("error in MPI_Irecv/MPI_Send "
    		 "(out_data is %d but should be 1)\n", out_data);
      return EXIT_FAILURE;
  }
  MPI_Type_free(&type);

  MPI_Finalize();

  return 0;
}
