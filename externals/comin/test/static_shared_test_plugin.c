/**
   @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#include <stdio.h>

#include <comin.h>

// test some basic stuff to ensure interaction with host model works

int ep = EP_ATM_WRITE_OUTPUT_BEFORE;

void static_linking_test_foo_printer(){
#ifdef COMIN_STATIC_PLUGIN_IS_STATIC
  fprintf(stderr, "foo STATIC\n");
#else
  fprintf(stderr, "foo SHARED\n");
#endif
}

void static_linking_test_comin_main(){
  int ierr = 0;
  comin_callback_register(ep, static_linking_test_foo_printer, &ierr);
  if(ierr != 0)
    comin_plugin_finish("static_linking_test_comin_main", "cannot register callback");
}
