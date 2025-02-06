/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#ifndef PYCOMIN_CALLBACKS_H
#define PYCOMIN_CALLBACKS_H

#define PY_SSIZE_T_CLEAN
#include <Python.h>

std::vector<PyMethodDef> py_comin_callbacks_methods();
void py_comin_generic_callback();

#endif
