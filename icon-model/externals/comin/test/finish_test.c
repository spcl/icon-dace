/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#include <stdbool.h>
#include <comin.h>
#include <stdio.h>

void setup(){
  // expected to fail
  comin_plugin_finish("setup", "This is a dummy error to test comin_plugin_finish!");
}
