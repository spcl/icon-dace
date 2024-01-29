/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <filesystem>

#include <dlfcn.h>

#include "comin.h"

#include "comin.py.h"

#include "callbacks.h"
#include "variables.h"
#include "descrdata.h"

static struct PyModuleDef py_comin_module = {
    PyModuleDef_HEAD_INIT,
    "_comin",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    NULL /*m_methods - will be set in main before calling PyModule_Create*/
};

static std::vector<PyMethodDef> pyCominMethods;

PyMODINIT_FUNC
PyInit_comin(void)
{
  // collect methods from other cpp files:
  for(const auto& m : py_comin_callbacks_methods())
    pyCominMethods.push_back(m);
  for(const auto& m : py_comin_variables_methods())
    pyCominMethods.push_back(m);
  for(const auto& m : py_comin_descrdata_methods())
    pyCominMethods.push_back(m);
  pyCominMethods.push_back({NULL, NULL, 0, NULL}); /* Sentinel */

  py_comin_module.m_methods = pyCominMethods.data();
  PyObject* pSelf = PyModule_Create(&py_comin_module);
  if (pSelf == NULL)
    return NULL;

  return pSelf;
}

static int py_comin_instance_counter = 0;

extern "C" {
  void comin_main() {

    using namespace std::string_literals;

    int ierr = 0;
    int mpi_rank = comin_parallel_get_host_mpi_rank();
    if (mpi_rank == 0)
      std::cerr << "setup_python_adapter" << std::endl;

    if(!Py_IsInitialized()) {
      // this is a workaround for a python problem that occurs if python
      // is embedded in a library that is loaded with dlopen.
      // See for example:
      // - https://bugs.python.org/issue4434
      // - https://stackoverflow.com/questions/8302810/undefined-symbol-in-c-when-loading-a-python-shared-library
      // - https://stackoverflow.com/questions/64295279/so-type-plugin-with-embedded-python-interpreter-how-to-link-load-libpython
      char libname[17];
      sprintf(libname, "libpython3.%d.so", PY_MINOR_VERSION);
      void* handle = dlopen(libname, RTLD_LAZY | RTLD_GLOBAL);
      if(handle == nullptr){
        if (mpi_rank == 0) {
          std::cerr << "Cannot load " << libname << ": " << dlerror() << std::endl;
          comin_plugin_finish(__func__, "Plugin cannot be loaded.");
        }
        else{
          if (mpi_rank == 0)  std::cerr << "Python adapter: " << libname << " loaded!" << std::endl;
        }
      }
      dlclose(handle);

      PyImport_AppendInittab("_comin", &PyInit_comin);
      Py_Initialize();
      std::string comin_py_c((char*)comin_py, comin_py_len);
      PyObject* compiled_comin =
        Py_CompileString(comin_py_c.c_str(), "comin.py", Py_file_input);
      PyImport_ExecCodeModule("comin", compiled_comin);
    }

    py_comin_instance_counter++;
    comin_callback_register(EP_DESTRUCTOR,
                            [](){
                              py_comin_generic_callback();
                              py_comin_instance_counter--;
                              if(py_comin_instance_counter == 0)
                                Py_Finalize();
                            }, &ierr);
    if (ierr != 0)
      return;

    // Dictionary to store the plugins global variables
    // independently from other python plugins
    PyObject* globals = PyDict_New();

    int ilen = -1;
    const char* plugin_options_c = NULL;
    comin_current_get_plugin_options(&plugin_options_c, &ilen, &ierr);
    std::string plugin_options(plugin_options_c, ilen);

    if (ierr != 0)
      return;
    try {
      if (mpi_rank == 0)  std::cerr << "Running python script " << plugin_options << std::endl;
      FILE* file = fopen(plugin_options.c_str(), "r");
      if (file == NULL) throw std::runtime_error("Cannot read "s + plugin_options);
      else {
        PyObject* pyResult = PyRun_File(file, plugin_options.c_str(), Py_file_input, globals, globals);
        Py_XDECREF(pyResult);
        fclose(file);
        if (PyErr_Occurred()){
          PyErr_Print();
          throw std::runtime_error("Error while executing "s + plugin_options);
        }
      }
    } catch ( std::exception& err ) {
      std::cerr << err.what() << std::endl;
      comin_plugin_finish(__func__, "Error while executing script");
      return;
    }
  }
}
