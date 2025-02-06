/**
   @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#include "comin.h"
#include <stdio.h>
#include <string.h>

/*
  This plugin checks several situation in which an error should occur and checks the error code.
 */

void secondary_constructor();
void check_var_get();

void dummy_callback(){}

void comin_main(){

  // comin_var_get is not allowed to be called from the primary
  // contructor
  struct t_comin_var_descriptor foo_descr = {.id = 1};
  strncpy(foo_descr.name, "foo", MAX_LEN_VAR_NAME);
  int context[1] = {EP_DESTRUCTOR};
  comin_error_set_errors_return(true);
  void* var_handle = comin_var_get(sizeof(context), context,
                                   foo_descr, COMIN_FLAG_READ);
  if(var_handle == NULL && comin_error_get() == COMIN_ERROR_VAR_GET_OUTSIDE_SECONDARY_CONSTRUCTOR)
    fprintf(stderr, "Check successful: comin_var_get returns NULL in primary constructor\n");
  else{
    comin_plugin_finish("test_errorcodes_plugin: comin_main", "comin_var_get did not return a null ptr in the primary constructor");
  }
    comin_error_set_errors_return(false);

  comin_callback_register(EP_SECONDARY_CONSTRUCTOR, secondary_constructor);
  comin_callback_register(EP_ATM_TIMELOOP_BEFORE, check_var_get);

}


void secondary_constructor(){
  // Check that no further callbacks can be registered
  comin_error_set_errors_return(true);
  comin_callback_register(EP_DESTRUCTOR, dummy_callback);
  if(comin_error_get() == COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR){
    fprintf(stderr, "Check successful: callbacks cannot be registered after primary constructor\n");
  }else{
    comin_plugin_finish("test_errorcodes_plugin: secondary_constructor",
                 "Callback register restriction broken");
  }

  // check that no further variables can be requested
  struct t_comin_var_descriptor bar_descr = {.id = 1};
  strncpy(bar_descr.name, "bar", MAX_LEN_VAR_NAME);
  comin_var_request_add(bar_descr, false);
  if(comin_error_get() == COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR){
    fprintf(stderr, "Check successful: variables cannot be requested after primary constructor\n");
  }else{
    comin_plugin_finish("test_errorcodes_plugin: secondary_constructor",
                 "Variable request restriction broken");
  }
  comin_error_set_errors_return(false);
}

void check_var_get(){
  struct t_comin_var_descriptor foo_descr = {.id = 1};
  strncpy(foo_descr.name, "foo", MAX_LEN_VAR_NAME);
  int context[1] = {EP_DESTRUCTOR};
  comin_error_set_errors_return(true);
  void* var_handle = comin_var_get(sizeof(context), context,
                                   foo_descr, COMIN_FLAG_READ);
  if(var_handle == NULL && comin_error_get() == COMIN_ERROR_VAR_GET_OUTSIDE_SECONDARY_CONSTRUCTOR){
    fprintf(stderr, "Check successful: variable pointers cannot be retrieved after secondary constructor\n");
  }else{
    comin_plugin_finish("test_errorcodes_plugin: test_vars", "Variable pointer retrieval restriction broken");
  }
  comin_error_set_errors_return(true);
}
