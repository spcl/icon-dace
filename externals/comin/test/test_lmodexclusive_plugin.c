#include <string.h>
#include <stdio.h>
#include "comin.h"

void comin_main(){

  int ilen = -1;
  const char* plugin_name = NULL;
  comin_current_get_plugin_name(&plugin_name, &ilen);

  struct t_comin_var_descriptor foo_descr = {.id = 1};
  strncpy(foo_descr.name, "foo", MAX_LEN_VAR_NAME);

  struct t_comin_var_descriptor foo_exc_descr = {.id = 1};
  strncpy(foo_exc_descr.name, "foo_exc", MAX_LEN_VAR_NAME);

  int my_id = comin_current_get_plugin_id();
  if(my_id == 1){
    comin_var_request_add(foo_descr, false);

    comin_var_request_add(foo_exc_descr, true);
  }else{
    comin_error_set_errors_return(true);
    // 1. Variable was requested exclusive and is requested again
    comin_var_request_add(foo_exc_descr, false);
    if(comin_error_get() != COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE){
      comin_plugin_finish("test_lmodexclusive: comin_main",
                   "Request lmodexclusive check (variable exists lmodexclusive) broken");
    }else{
      fprintf(stderr, "Check successful: variables cannot be requested if exist with lmodexclusive=true\n");
    }

    // 2. Variable exists and is now requested exclusive
    comin_var_request_add(foo_descr, true);
    if(comin_error_get() != COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE){
      comin_plugin_finish("simple_fortran_plugin: comin_main",
                   "Request lmodexclusive check (variable exists) broken");
    }else{
      fprintf(stderr, "Check successful: variables cannot be requested lmodexclusive=true if exists\n");
    }
  }
}
