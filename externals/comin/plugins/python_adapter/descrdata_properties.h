/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#ifndef PY_COMIN_DESCRDATA_PROPERTIES
#define PY_COMIN_DESCRDATA_PROPERTIES

#include "comin.h"

static PyObject* BuildValue_FromPtr(std::string ctype, void* data){
  if(ctype == "double") return Py_BuildValue("f", *(double*)data);
  if(ctype == "int") return Py_BuildValue("i", *(int*)data);
  return PyErr_Format(PyExc_ValueError, "Unknown data type: %s", ctype);
}

// helper function to generate dict from properties array
static PyObject* py_comin_descrdata_properties_to_dict(const comin_descrdata_property_t* properties){
  PyObject* dict = PyDict_New();
  const comin_descrdata_property_t* prop = properties;
  while(prop->name != NULL){
    PyDict_SetItemString(dict, prop->name, PyCapsule_New((void*)prop, "descrdata_property", NULL));
    prop++;
  }
  return dict;
}

static PyObject* py_comin_descrdata_eval_property(PyObject* self, PyObject* args, PyObject* kwargs){
  PyObject *prop_cap;
  int jg = 0;
  static char const *kwlist[] = {"property", "jg", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|i", (char**)&kwlist, &prop_cap, &jg))
    return NULL;
  comin_descrdata_property_t* prop = (comin_descrdata_property_t*)PyCapsule_GetPointer(prop_cap, "descrdata_property");
  if (prop == NULL)
    return NULL;
  int ierr = 0;
  void* dataptr;
  if (std::string(prop->datatype) == "void"){ // return the dict of subproperties
    return py_comin_descrdata_properties_to_dict(prop->subtypes);
  }else if(prop->ndims == 0){ // return scalar
    if(std::string(prop->datatype) == "int"){
      int val;
      if(prop->has_jg)
        val = ((int(*)(int))prop->get_function)(jg);
      else
        val = ((int(*)())prop->get_function)();
      return PyLong_FromLong(val);
    }
    if(std::string(prop->datatype) == "double"){
      double val;
      if(prop->has_jg)
        val = ((double(*)(int))prop->get_function)(jg);
      else
        val = ((double(*)())prop->get_function)();
      return PyFloat_FromDouble(val);
    }
    if(std::string(prop->datatype) == "bool"){
      bool val;
      if(prop->has_jg)
        val = ((bool(*)(int))prop->get_function)(jg);
      else
        val = ((bool(*)())prop->get_function)();
      if(val) Py_RETURN_TRUE;
      else Py_RETURN_FALSE;
    }
    return PyErr_Format(PyExc_ValueError, "datatype %s for scalar value %s not implemented", prop->datatype, prop->name);
  }else{
    std::vector<int> arr_size(prop->ndims);
    if(prop->has_jg)
      ((void(*)(int, void**, int*, int*))prop->get_function)(jg, &dataptr, arr_size.data(), &ierr);
    else
      ((void(*)(void**, int*, int*))prop->get_function)(&dataptr, arr_size.data(), &ierr);
    if (dataptr == NULL)
      return PyErr_Format(PyExc_ValueError, "%s get function returned NULL", prop->name);
    if(prop->ndims == 0){ // return the value directly
      return BuildValue_FromPtr(prop->datatype, dataptr);
    }
    if(std::string(prop->datatype) == "char"  && prop->ndims == 1)
      return Py_BuildValue("s#", dataptr, arr_size[0]);
    // return MemoryView
    Py_buffer buffer;
    if(std::string(prop->datatype) == "double"){
      fill_buffer<double>(&buffer, dataptr, arr_size.data(), prop->ndims, 1);
    }else if(std::string(prop->datatype) == "int"){
      fill_buffer<int>(&buffer, dataptr, arr_size.data(), prop->ndims, 1);
    }else if(std::string(prop->datatype) == "int8_t"){
      fill_buffer<int8_t>(&buffer, dataptr, arr_size.data(), prop->ndims, 1);
    }else
      return PyErr_Format(PyExc_ValueError, "datatype %s for property %s not implemented", prop->datatype, prop->name);
    buffer.obj = self;
    return PyMemoryView_FromBuffer(&buffer);
  }
}

#endif
