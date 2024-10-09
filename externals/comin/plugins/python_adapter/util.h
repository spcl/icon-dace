/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#ifndef UTIL_H
#define UTIL_H

#include <comin_global.inc>
#include <comin.h>
#include <Python.h>

// format characters - see: https://docs.python.org/3/library/struct.html#format-characters
template<class T>
const char* get_format_char();

template<>
inline const char* get_format_char<double>(){return "d";}
template<>
inline const char* get_format_char<float>(){return "f";}
template<>
inline const char* get_format_char<int>(){return "i";}
template<>
inline const char* get_format_char<int8_t>(){return "i1";}

template<class T>
void fill_buffer(Py_buffer* buffer, void* mem, int* shape, int ndims, int readonly = 1){
  buffer->obj = NULL;
  buffer->buf = (void*)mem;
  buffer->len = sizeof(T);
  for(int i = 0; i<ndims; buffer->len *= shape[i++]);
  buffer->readonly = readonly;
  buffer->itemsize = sizeof(T);
  buffer->format = (char*)get_format_char<T>(); // T
  buffer->ndim = ndims;
  buffer->shape = new Py_ssize_t[ndims];
  for(int i=0; i<ndims;++i)
    buffer->shape[i] = shape[i];

  buffer->strides = new Py_ssize_t[ndims];
  buffer->strides[0] = sizeof(T);
  for (int i=1; i<buffer->ndim; i++)
    buffer->strides[i] = buffer->strides[i-1] * buffer->shape[i-1];
  buffer->suboffsets = NULL;
  buffer->internal = NULL;
}

#endif
