/**
   @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "comin.h"

void check_ep_name(int ep, char* reference){
  printf("Checking %s (%d) ... ", reference, ep);
  char ep_name[MAX_LEN_EP_NAME+1];
  comin_callback_get_ep_name(ep, ep_name);
  if (strcmp(reference, ep_name)){
    printf("Wrong ep_name returned (%s != %s)\n", reference, ep_name);
    exit(2);
  }
  printf("done\n");
}

#define CHECK_EP_NAME(NAME) check_ep_name(NAME, #NAME)

int main(){
  CHECK_EP_NAME(EP_SECONDARY_CONSTRUCTOR);
  CHECK_EP_NAME(EP_FINISH);
  CHECK_EP_NAME(EP_ATM_CONVECTION_AFTER);
  CHECK_EP_NAME(EP_DESTRUCTOR);
  return 0;
}
