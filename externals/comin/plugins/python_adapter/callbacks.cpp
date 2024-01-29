/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <vector>
#include <string>
#include <map>
#include <stdexcept>

#include "comin.h"

#include "callbacks.h"

// stores a vector of callable python objects for each plugin id and entry point.
// In contrast to the fortran and c interface we explicitly allow to
// register multiple callbacks for every entrypoint. This makes it
// easy to combine multiple python 3rd party modules. However, we do not make
// any guarantees about the order in which the callbacks are called.
static std::map<int,
  std::map<int, std::vector<PyObject*>>> callbacks;

void py_comin_generic_callback(){
  int ep = 0;
  ep = comin_current_get_ep ();
  int plugin_id = comin_current_get_plugin_id();
  for(PyObject* fun : callbacks[plugin_id][ep]){
    // call it (no arguments supplied)
    PyObject *rv = PyObject_CallObject(fun, 0);
    Py_CLEAR(rv);
    if(PyErr_Occurred()){
      PyErr_Print();
      std::abort();
    }
  }
}

static PyObject *
py_comin_callback_register(PyObject */*self*/, PyObject *args)
{
  int ierr = 0;
  int ientry_point;
  PyObject* callback;
  if (!PyArg_ParseTuple(args, "iO", &ientry_point, &callback)) {
    return NULL;
  }
  int entry_point = ientry_point;
  if (!PyCallable_Check(callback)) {
    return PyErr_Format(PyExc_TypeError, "A callable is required");
  }

  int plugin_id = comin_current_get_plugin_id();
  if(callbacks[plugin_id].find(entry_point) == callbacks[plugin_id].end() &&
     entry_point != EP_DESTRUCTOR){
    comin_callback_register(ientry_point, &py_comin_generic_callback, &ierr);
  }
  Py_INCREF(callback);
  callbacks[plugin_id][entry_point].push_back(callback);
  Py_RETURN_NONE;
}

static PyObject* py_comin_EP_DESCTRUCTOR(){
  return PyLong_FromLong((long)EP_DESTRUCTOR);
}

static PyObject* py_comin_callback_get_ep_name(PyObject* /*self*/,
                                               PyObject* args, PyObject* kwargs){
  int iep;
  static char const * kwlist[] = {(char*)"iep", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", (char**)&kwlist, &iep)) {
    return NULL;
  }
  char out_ep_name[MAX_LEN_EP_NAME+1];
  int ierr = 0;
  comin_callback_get_ep_name(iep, out_ep_name, &ierr);
  if(ierr != 0)
    return PyErr_Format(PyExc_TypeError, "comin_callback_get_ep_name failed (iep = %d)", iep);
  return Py_BuildValue("s", out_ep_name);
}

std::vector<PyMethodDef> py_comin_callbacks_methods() {
  return {
    {"_callback_register",  py_comin_callback_register, METH_VARARGS,
     "Registers callback to ICON"},
    {"_EP_DESTRUCTOR", (PyCFunction)py_comin_EP_DESCTRUCTOR, METH_NOARGS, ""},
    {"callback_get_ep_name", (PyCFunction)py_comin_callback_get_ep_name, METH_VARARGS | METH_KEYWORDS,
     "C function signature: void comin_callback_get_ep_name(int iep, char out_ep_name[MAX_LEN_EP_NAME+1], int* ierr)"},
  };
}
