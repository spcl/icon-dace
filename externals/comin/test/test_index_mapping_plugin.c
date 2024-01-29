/**
   @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#include <stdio.h>
#include "comin.h"

void comin_main(){
  /* Test local / global index mappings:

     We get the global cell index of this MPI rank's cell with
     local index 42... and then re-translate this index to the local
     index.
  */
  int local_idx = 42;
  int jg = 1;
  int* glb;
  int arr_size[1];
  comin_descrdata_get_domain_cells_glb_index(jg, &glb, arr_size);
  if( arr_size[0] > local_idx){
    int glb_index = glb[local_idx];
    int lcl = comin_descrdata_index_lookup_glb2loc_cell(1,glb_index) - 1; // subtract one due to fortran indexing
    fprintf(stderr, "global index: %d\n", glb_index);
    fprintf(stderr, "local index: %d (reference)    %d (looked up)\n", local_idx, lcl);
    if (lcl != local_idx)
      comin_plugin_finish("simple_fortran_constructor", "global index lookup check failed!");
  }
}
