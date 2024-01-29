/*
  Example plugin to demonstrate the usage of yaxt in a ComIn plugin.
 */

#include <stdio.h>
#include <string.h>

#include <comin.h>
#include <mpi.h>
#include <yaxt.h>

#define ASSERT(COND) if(!(COND)){ fprintf(stderr, "assertion failed: %s\n", #COND); MPI_Abort(MPI_COMM_WORLD, 1); }

const int jg = 1;
int rank = -1;
MPI_Comm plugin_comm = MPI_COMM_NULL;
void* temp_var = NULL;
Xt_int* idxlist = NULL;
double* buffer = NULL;
int prognostic_cells = 0;
Xt_redist yaxt_redist;

void secondary_constructor();
void offload();
void destructor();

void comin_main(){
  int ierr = 0;

  rank = comin_parallel_get_host_mpi_rank();
  comin_callback_register(EP_SECONDARY_CONSTRUCTOR, secondary_constructor, &ierr);
  if(ierr != 0)
    comin_plugin_finish("yaxt_c_plugin", "callback registration (EP_SECONDARY_CONSTRUCTOR) failed");

  comin_callback_register(EP_ATM_WRITE_OUTPUT_BEFORE, offload, &ierr);
  if(ierr != 0)
    comin_plugin_finish("yaxt_c_plugin", "callback registration (EP_ATM_WRTITE_OUTPUT_BEFORE) failed");

  comin_callback_register(EP_DESTRUCTOR, destructor, &ierr);
  if(ierr != 0)
    comin_plugin_finish("yaxt_c_plugin", "callback registration (EP_DESTRUCTOR) failed");

  plugin_comm = MPI_Comm_f2c(comin_parallel_get_plugin_mpi_comm());

  if(!xt_initialized())
    xt_initialize(plugin_comm);

  int* decomp_domain = NULL;
  int decomp_domain_size[2] = {0,0};
  comin_descrdata_get_domain_cells_decomp_domain(jg, &decomp_domain, decomp_domain_size);

  int* glb_index = NULL;
  int glb_index_size = 0;
  comin_descrdata_get_domain_cells_glb_index(jg, &glb_index, &glb_index_size);

  ASSERT(decomp_domain_size[0]*decomp_domain_size[1] > glb_index_size);

  idxlist = malloc(glb_index_size*sizeof(*idxlist));
  prognostic_cells = 0;
  for(int i=0; i<glb_index_size; ++i){
    if(decomp_domain[i]==0){
      idxlist[prognostic_cells] = glb_index[i];
      prognostic_cells++;
    }
  }
  Xt_idxlist idxvec = xt_idxvec_new(idxlist, prognostic_cells);

  Xt_idxlist empty_idxlist = xt_idxempty_new();

  Xt_xmap xmap = xt_xmap_all2all_new(
    idxvec, empty_idxlist, plugin_comm);

  yaxt_redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
  xt_xmap_delete(xmap);
  xt_idxlist_delete(empty_idxlist);
  xt_idxlist_delete(idxvec);

  buffer = malloc(prognostic_cells*sizeof(*buffer));
}

void secondary_constructor(){
  int ep = EP_ATM_WRITE_OUTPUT_BEFORE;
  struct t_comin_var_descriptor var_desc;
  strcpy(var_desc.name, "temp");
  var_desc.id = 1;
  temp_var = comin_var_get(1, &ep, var_desc, COMIN_FLAG_READ);
}

void offload(){
  int ierr = 0;
  // this callback is only called in entrypoints outside of domain loops
  ASSERT(comin_current_get_domain_id() == COMIN_DOMAIN_OUTSIDE_LOOP);
  double* temp_ptr = comin_var_to_3d(temp_var); // ensures [jc, jb, jk]
  int shape[5];
  comin_var_get_shape(temp_var, shape, &ierr);
  int nlev = comin_descrdata_get_domain_nlev(jg);
  for(int i=0; i<prognostic_cells; ++i){
    int idx = idxlist[i];
    int idx0 = idx%shape[0];
    int idx2 = idx/shape[0];
    buffer[i] = temp_ptr[idx0 + nlev*shape[0] + idx2*shape[0]*shape[1]];
  }

  xt_redist_s_exchange1(yaxt_redist, buffer, NULL);
}

void destructor(){
  free(idxlist);
  free(buffer);
}
