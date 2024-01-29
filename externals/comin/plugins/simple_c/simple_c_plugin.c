/*
  Example plugin for the ICON Community Interface (ComIn)
  with basic (not MPI-parallel) callbacks and accessing variables and descriptive data structures.

  Note that in order to demonstrate ComIn's language interoperability,
  a similary plugin has been implemented in FORTRAN, see the subdirectory
  "simple".

  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

  SPDX-License-Identifier: BSD-3-Clause

  Please see the file LICENSE in the root of the source tree for this code.
  Where software is supplied by third parties, it is indicated in the
  headers of the routines.
*/


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <comin.h>

// points to the pointer to the buffer (buffer swapping)
void* pres, *simple_c_var, *simple_c_tracer;


void simple_c_constructor(){
  unsigned int major, minor, patch;
  comin_setup_get_version(&major, &minor, &patch);
  fprintf(stderr, "ComIn v%u.%u.%u simple_c_constructor called!\n", major, minor, patch);

  int ierr = 0;

  struct t_comin_var_descriptor desc;
  for(void* it=comin_var_get_descr_list_head();
      it != NULL;
      it = comin_var_get_descr_list_next(it)){
    comin_var_get_descr_list_var_desc(it, &desc, &ierr);
    if (ierr != 0) comin_plugin_finish("simple_c: simple_c_constructor", "comin_var_get_descr_list_var_desc failed for list_var_desc");
    fprintf(stderr, "found variable: %s on domain id %d\n", desc.name, desc.id);
  }

  struct t_comin_var_descriptor pres_d = {.name="pres", .id=1};
  int before_output = EP_ATM_WRITE_OUTPUT_BEFORE;
  pres = comin_var_get(1, &before_output, pres_d, COMIN_FLAG_READ);
  if(pres == NULL){
    fprintf(stderr, "Internal error!\n");
    abort();
  }

  struct t_comin_var_descriptor var_d = {.name="simple_c_var", .id=1};
  simple_c_var = comin_var_get(1, &before_output, var_d, COMIN_FLAG_WRITE);
  if(simple_c_var == NULL){
    fprintf(stderr, "Internal error!\n");
    abort();
  }

  struct t_comin_var_descriptor tracer_d = {.name="simple_c_tracer", .id=1};
  simple_c_tracer = comin_var_get(1, &before_output, tracer_d, COMIN_FLAG_WRITE);
  if(simple_c_tracer == NULL){
    fprintf(stderr, "Internal error!\n");
    abort();
  }

}

void simple_c_diagfct(){
  int jg, ierr = 0;
  fprintf(stderr, "simple_c_diagfct called!\n");
  jg = comin_current_get_domain_id();
  fprintf(stderr, "currently on domain %i\n", jg);
  int pres_shape[5], tracer_shape[5];
  comin_var_get_shape(pres, pres_shape, &ierr);
  if (ierr != 0) comin_plugin_finish("simple_c: comin_main", "comin_var_get_shape for pres");
  double* simple_c_var_data = comin_var_get_ptr(simple_c_var);
  double* pres_data = comin_var_get_ptr(pres);
  for(int i = 0; i<pres_shape[0]*pres_shape[1]*pres_shape[2]*pres_shape[3]*pres_shape[4]; ++i){
    simple_c_var_data[i] = pres_data[i] +42.;
  }
  comin_var_get_shape(simple_c_tracer, tracer_shape, &ierr);
  if (ierr != 0) comin_plugin_finish("simple_c: comin_main", "comin_var_get_shape failed for simple_c_tracer");
  double* simple_c_tracer_data = comin_var_to_3d(simple_c_tracer);
  for(int i = 0; i<tracer_shape[0]*tracer_shape[1]*tracer_shape[2]; ++i){
    simple_c_tracer_data[i] /= 1337.;
  }

}

void simple_c_destructor(){
  fprintf(stderr, "simple_c_destructor called!\n");
}


void comin_main(){
  int ierr = 0;
  int ilen = -1;
  const char* plugin_name = NULL;
  comin_current_get_plugin_name(&plugin_name, &ilen, &ierr);
  if (ierr != 0) comin_plugin_finish("simple_c: comin_main", "comin_current_get_plugin_info failed");

  int plugin_id = comin_current_get_plugin_id();
  fprintf(stderr, "plugin %s has id %d\n", plugin_name, plugin_id);
  struct t_comin_var_descriptor simple_var_d = {.name="simple_c_var", .id=1};
  comin_var_request_add(simple_var_d, true, &ierr);

  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "comin_var_request_add failed for simple_c_var.");
  comin_metadata_set_integer(simple_var_d, "zaxis_id", COMIN_ZAXIS_3D, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'zaxis_id' failed.");
  comin_metadata_set_logical(simple_var_d, "restart", false, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'restart' failed.");
  comin_metadata_set_logical(simple_var_d, "tracer", false, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'tracer' failed.");
  comin_metadata_set_integer(simple_var_d, "tracer_vlimit", 0, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'tracer_vlimit' failed.");
  comin_metadata_set_integer(simple_var_d, "tracer_hlimit", 0, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'tracer_hlimit' failed.");

  struct t_comin_var_descriptor simple_tracer_d = {.name="simple_c_tracer", .id=-1};
  comin_var_request_add(simple_tracer_d, false, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "comin_var_request_add failed for simple_c_tracer.");
  comin_metadata_set_integer(simple_tracer_d, "zaxis_id", COMIN_ZAXIS_3D, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'zaxis_id' failed.");
  comin_metadata_set_logical(simple_tracer_d, "restart", false, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'restart' failed.");
  comin_metadata_set_logical(simple_tracer_d, "tracer", true, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'tracer' failed.");
  comin_metadata_set_integer(simple_tracer_d, "tracer_vlimit", 0, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'tracer_vlimit' failed.");
  comin_metadata_set_integer(simple_tracer_d, "tracer_hlimit", 0, &ierr);
  if (ierr != 0)  comin_plugin_finish("simple_c: comin_main", "setting metadata 'tracer_hlimit' failed.");

  comin_callback_register(EP_SECONDARY_CONSTRUCTOR, &simple_c_constructor, &ierr);
  if (ierr != 0) comin_plugin_finish("simple_c: comin_main", "comin_callback_register failed for simple_c_constructor");
  comin_callback_register(EP_ATM_WRITE_OUTPUT_BEFORE, &simple_c_diagfct, &ierr);
  if (ierr != 0) comin_plugin_finish("simple_c: comin_main", "comin_callback_register failed for simple_c_diagfct");
  comin_callback_register(EP_DESTRUCTOR, &simple_c_destructor, &ierr);
  if (ierr != 0) comin_plugin_finish("simple_c: comin_main", "comin_callback_register failed for simple_c_destructor");
 
  char ep_name[MAX_LEN_EP_NAME+1];
  comin_callback_get_ep_name(EP_DESTRUCTOR, ep_name, &ierr);
  if (ierr != 0) comin_plugin_finish("simple_c: comin_main", "comin_callback_get_ep_name failed for ep_name");
  if (strncmp(ep_name, "EP_DESTRUCTOR", (size_t) 13)) {
    char output_text[255];
    sprintf(output_text,"Expected EP_DESTRUCTOR; got |%s|\n", ep_name);
    comin_plugin_finish("simple_c: comin_main", output_text);
  }
  // TODO: Access descriptive data structures

  /* access to comin descriptive global data (exemplary)*/
  int n_dom = comin_descrdata_get_global_n_dom();
  fprintf(stderr, "n_dom: %d\n", n_dom);
  int max_dom = comin_descrdata_get_global_max_dom();
  fprintf(stderr, "max_dom: %d\n", max_dom);
  int nproma = comin_descrdata_get_global_nproma();
  fprintf(stderr, "nproma: %d\n", nproma);
  int min_rlcell_int = comin_descrdata_get_global_min_rlcell_int();
  fprintf(stderr, "min_rlcell_int: %d\n", min_rlcell_int);
}
