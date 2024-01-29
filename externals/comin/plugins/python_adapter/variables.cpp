/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <limits>

#include "comin.h"

#include "util.h"
#include "variables.h"


static PyObject* py_comin_var_request_add(PyObject* /*self*/, PyObject* args, PyObject* kwargs)
{
  t_comin_var_descriptor var_descr;
  static char const* kwlist[] = {"var_descriptor", "lmodexclusive", NULL};
  char* name;
  Py_ssize_t len;
  bool lmodexclusive;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)p", (char**)&kwlist,
                                  &name, &len, &(var_descr.id), &lmodexclusive))
    return NULL;

  if(len > MAX_LEN_VAR_NAME)
    return PyErr_Format(PyExc_ValueError, "Variable name to long!");

  strncpy(var_descr.name, name, MAX_LEN_VAR_NAME);

  int ierr = 0;
  comin_var_request_add(var_descr, lmodexclusive, &ierr);
  if(ierr != 0)
    return PyErr_Format(PyExc_ValueError, "request_add_var failed (ierr=%d, name=%s, id=%d)",
                        ierr, var_descr.name, var_descr.id);
  Py_RETURN_NONE;
}

static PyObject* py_comin_metadata_set(PyObject* /*self*/, PyObject* args, PyObject* kwargs){
  t_comin_var_descriptor var_descr;
  char* key, *name;
  Py_ssize_t len;
  PyObject* val;
  static char const* kwlist[] = {"var_descriptor", "key", "value", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)sO", (char**)&kwlist,
                                  &name, &len, &(var_descr.id), &key, &val))
    return NULL;

  if(len > MAX_LEN_VAR_NAME)
    return PyErr_Format(PyExc_ValueError, "Variable name to long!");

  strncpy(var_descr.name, name, MAX_LEN_VAR_NAME);

  int ierr = 0;
  if(PyBool_Check(val)){
    comin_metadata_set_logical(var_descr, key, val==Py_True, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_set_logical (ierr=%d, name=%s, id=%d, key=%s)",
                          ierr, var_descr.name, var_descr.id, key);
    Py_RETURN_NONE;
  }
  if(PyLong_Check(val)){
    int overflow = 0;
    long val_l = PyLong_AsLongAndOverflow(val, &overflow);
    if (val_l > std::numeric_limits<int>::max() || val_l < std::numeric_limits<int>::min() || overflow!=0){
      return PyErr_Format(PyExc_ValueError, "value is out of the range of int");
    }
    comin_metadata_set_integer(var_descr, key, val_l, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_set_integer (ierr=%d, name=%s, id=%d, key=%s)",
                          ierr, var_descr.name, var_descr.id, key);
    Py_RETURN_NONE;
  }
  if(PyFloat_Check(val)){
    double val_d = PyFloat_AsDouble(val);
    if(val_d == -1. && PyErr_Occurred())
      return PyErr_Format(PyExc_ValueError, "value is not a float");
    comin_metadata_set_real(var_descr, key, val_d, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_set_real (ierr=%d, name=%s, id=%d, key=%s)",
                          ierr, var_descr.name, var_descr.id, key);
    Py_RETURN_NONE;
  }
  if(PyUnicode_Check(val)){
    const char* val_str = PyUnicode_AsUTF8(val);
    if(val_str == NULL)
      return PyErr_Format(PyExc_ValueError, "value cannot be converted to char* (UTF-8)");
    comin_metadata_set_character(var_descr, key, val_str, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_set_character (ierr=%d, name=%s, id=%d, key=%s)",
                          ierr, var_descr.name, var_descr.id, key);
    Py_RETURN_NONE;
  }
  PyTypeObject* type = val->ob_type;
  return PyErr_Format(PyExc_ValueError, "comin_metadata_set is not implemented for the provided object type (type=%s, name=%s, id=%d, key=%s)",
                      type->tp_name, var_descr.name, var_descr.id, key);
}

static PyObject* py_comin_metadata_get(PyObject* /*self*/, PyObject* args, PyObject* kwargs)
{
  t_comin_var_descriptor var_desc;
  static char const *kwlist[] = {"var_descriptor", "key", NULL};
  char* name;
  Py_ssize_t len;
  char* key;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)s", (char**)&kwlist,
                                  &name, &len, &(var_desc.id), &key))
    return NULL;
  strncpy(var_desc.name, name, len+1);
  int type = comin_metadata_get_typeid(key);
  int ierr = 0;
  if(type == COMIN_TYPEID_LOGICAL){
    _Bool val;
    comin_metadata_get_logical(var_desc, key, &val, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_get_logical failed (key=%s, ierr = %d)", key, ierr);
    if(val)
      Py_RETURN_TRUE;
    else
      Py_RETURN_FALSE;
  }
  if(type == COMIN_TYPEID_INTEGER){
    int val = 0;
    comin_metadata_get_integer(var_desc, key, &val, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_get_logical failed (key=%s, ierr = %d)", key, ierr);
    return PyLong_FromLong(val);
  }
  if(type == COMIN_TYPEID_REAL){
    double val = 0.;
    comin_metadata_get_real(var_desc, key, &val, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_get_real failed (key=%s, ierr = %d)", key, ierr);
    return PyFloat_FromDouble(val);
  }
  if(type == COMIN_TYPEID_CHARACTER){
    const char* val = NULL;
    int len = -1;
    comin_metadata_get_character(var_desc, key, &val, &len, &ierr);
    if(ierr != 0)
      return PyErr_Format(PyExc_ValueError, "comin_metadata_get_character failed (key=%s, ierr = %d)", key, ierr);
    return PyUnicode_FromStringAndSize(val, len);
  }
  return PyErr_Format(PyExc_ValueError, "Undefined datatype for key %s",
                      key);
}

static PyObject* py_comin_var_get(PyObject* /*self*/, PyObject* args, PyObject* kwargs){
  PyObject* context;
  t_comin_var_descriptor var_desc;
  int flag = -1;

  static char const *kwlist[] = {"context", "var_descriptor", "flag", NULL};
  char* name;
  Py_ssize_t len;
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O(s#i)|i", (char**)&kwlist,
                                  &context, &name, &len, &(var_desc.id), &flag))
    return NULL;
  strncpy(var_desc.name, name, len+1);

  if(!PyList_Check(context))
    return PyErr_Format(PyExc_ValueError, "PyCominVar requires list of integers as contexts");
  std::vector<int> icontext(PyList_Size(context));
  for(size_t i=0; i<icontext.size();++i)
    icontext[i] = (int)PyLong_AsLong(PyList_GetItem(context, i));

  void* var = comin_var_get(icontext.size(), icontext.data(), var_desc, flag);
  if(var == NULL)
    return PyErr_Format(PyExc_ValueError, "Variable does not exist (name=%s, id=%d)",
                        var_desc.name, var_desc.id);
  return PyCapsule_New(var, "var", NULL);
}

static PyObject* py_comin_var_get_buffer(PyObject* self, PyObject* args, PyObject* kwargs){
  PyObject* handle_cap;
  static char const *kwlist[] = {"handle", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                  &handle_cap))
    return NULL;
  int ierr = 0;
  void* handle = PyCapsule_GetPointer(handle_cap, "var");

  double* ptr = comin_var_get_ptr(handle);
  if(ptr == NULL)
    return PyErr_Format(PyExc_ValueError, "comin_var_get_ptr failed");
  int shape[5];
  comin_var_get_shape(handle, shape, &ierr);
  if(ierr != 0)
    return PyErr_Format(PyExc_ValueError, "comin_var_get_shape failed (ierr=%d)", ierr);
  Py_buffer buffer;
  fill_buffer<double>(&buffer, ptr, shape, 5, 0);
  buffer.obj = self;
  return PyMemoryView_FromBuffer(&buffer);
}

static PyObject* py_comin_var_get_pos(PyObject* /*self*/, PyObject* args, PyObject* kwargs){
  PyObject* handle_cap;
  static char const *kwlist[] = {"handle", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                  &handle_cap))
    return NULL;
  void* handle = PyCapsule_GetPointer(handle_cap, "var");
  int pos_jc, pos_jk, pos_jb, pos_jn, ierr = 0;
  comin_var_get_pos(handle, &pos_jc, &pos_jk, &pos_jb, &pos_jn, &ierr);
  if(ierr != 0)
    return PyErr_Format(PyExc_ValueError, "comin_var_get_pos failed (ierr = %d)", ierr);
  return Py_BuildValue("iiii", pos_jc, pos_jk, pos_jb, pos_jn);
}

static PyObject* py_comin_var_get_ncontained(PyObject* /*self*/, PyObject* args, PyObject* kwargs){
  PyObject* handle_cap;
  static char const *kwlist[] = {"handle", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                  &handle_cap))
    return NULL;
  void* handle = PyCapsule_GetPointer(handle_cap, "var");
  int ncontained, ierr = 0;
  comin_var_get_ncontained(handle, &ncontained, &ierr);
  if(ierr != 0)
    return PyErr_Format(PyExc_ValueError, "comin_var_get_ncontained failed (ierr = %d)", ierr);
  return PyLong_FromLong(ncontained);
}

static PyObject* py_comin_var_get_descr_list_head()
{
  void* head = comin_var_get_descr_list_head();
  return PyCapsule_New(head, "var_descr_list", NULL);
}

static PyObject* py_comin_var_get_descr_list_next(PyObject* /*self*/, PyObject* args, PyObject* kwargs)
{
  PyObject* py_current;
  static char const *kwlist[] = {"current", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                  &py_current))
    return NULL;
  void* current = PyCapsule_GetPointer(py_current, "var_descr_list");
  void* next = comin_var_get_descr_list_next(current);
  if(next == NULL)
    Py_RETURN_NONE;
  else
    return PyCapsule_New(next, "var_descr_list", NULL);
}

static PyObject* py_comin_var_get_descr_list_var_desc(PyObject* /*self*/, PyObject* args, PyObject* kwargs)
{
  PyObject* py_current;
  static char const *kwlist[] = {"current", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                  &py_current))
    return NULL;
  void* current = PyCapsule_GetPointer(py_current, "var_descr_list");
  int ierr = 0;
  t_comin_var_descriptor var_desc;
  comin_var_get_descr_list_var_desc(current, &var_desc, &ierr);
  if(ierr != 0)
    return PyErr_Format(PyExc_ValueError, "comin_var_get_descr_list_var_desc failed (ierr = %d)", ierr);
  return Py_BuildValue("si", var_desc.name, var_desc.id);
}

std::vector<PyMethodDef> py_comin_variables_methods() {
  return {
    {"var_request_add", (PyCFunction)py_comin_var_request_add,
     METH_VARARGS | METH_KEYWORDS,
     "Request the host model to add a variable, arguments: name string, domain id, [metadata key], [metadata value], ..."},
    {"metadata_get", (PyCFunction)py_comin_metadata_get,
     METH_VARARGS | METH_KEYWORDS, "retrieve metadata, arguments: name string, domain id, metadata key"},
    {"_metadata_set", (PyCFunction)py_comin_metadata_set,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"_var_get", (PyCFunction)py_comin_var_get,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"_var_get_buffer", (PyCFunction)py_comin_var_get_buffer,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"_var_get_pos", (PyCFunction)py_comin_var_get_pos,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"_var_get_ncontained", (PyCFunction)py_comin_var_get_ncontained,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"_var_get_descr_list_head", (PyCFunction)py_comin_var_get_descr_list_head,
     METH_NOARGS, ""},
    {"_var_get_descr_list_next", (PyCFunction)py_comin_var_get_descr_list_next,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"_var_get_descr_list_var_desc", (PyCFunction)py_comin_var_get_descr_list_var_desc,
     METH_VARARGS | METH_KEYWORDS, ""},
  };
}
