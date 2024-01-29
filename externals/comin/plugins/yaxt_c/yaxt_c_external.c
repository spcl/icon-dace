#include <stdio.h>

#include <yaxt.h>
#include <mpi.h>
#include <mpi_handshake.h>

double mean(double* data, int size){
  double sum = 0;
  for(int i=0; i<size; ++i)
    sum += data[i];
  return sum/size;
}

int main(int argc, char** argv){

  MPI_Init(&argc, &argv);

  const char* group_names[] = {"comin", "yaxt_c_external"};
  MPI_Comm group_comms[2];
  mpi_handshake(group_names, group_comms, 2, MPI_COMM_WORLD);

  MPI_Comm plugin_comm;
  const char* plugin_comm_name = "yaxt_c";
  mpi_handshake(&plugin_comm_name, &plugin_comm, 1, group_comms[0]);

  xt_initialize(plugin_comm);

  int ncells = atoi(argv[1]);
  struct Xt_stripe stripe = {.start=1, .stride=1, .nstrides=ncells};
  Xt_idxlist tgt_idxlist = xt_idxstripes_new( &stripe, 1 );
  Xt_idxlist empty_idxlist = xt_idxempty_new();

  Xt_xmap xmap = xt_xmap_all2all_new(
    empty_idxlist, tgt_idxlist, plugin_comm);

  Xt_redist yaxt_redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
  xt_xmap_delete(xmap);
  xt_idxlist_delete(empty_idxlist);
  xt_idxlist_delete(tgt_idxlist);

  double* buffer = malloc(ncells*sizeof(*buffer));
  int ntimesteps = atoi(argv[2]);
  for(int i = 0; i<ntimesteps; ++i){
    xt_redist_s_exchange1(yaxt_redist, NULL, buffer);
    // now we have the whole field available and can compute things
    printf("mean: %.2lf\n", mean(buffer, ncells));
  }
  MPI_Finalize();
  return 0;
}
