if (COMIN_ENABLE_EXAMPLES)
  add_subdirectory(simple_fortran)
  add_subdirectory(yaxt_fortran)
  add_subdirectory(calc_water_column)
  add_subdirectory(yaxt_c)
  add_subdirectory(simple_c)
  add_subdirectory(yac_input)

  if (BUILD_TESTING)

  comin_add_test(NAME multi_plugin)
  comin_test_add_plugin(TEST multi_plugin
    NAME "simple_c_plugin"
    PLUGIN_LIBRARY "$<TARGET_FILE:simple_c_plugin>")
  comin_test_add_plugin(TEST multi_plugin
    NAME "simple_fortran_plugin"
    PLUGIN_LIBRARY "$<TARGET_FILE:simple_fortran_plugin>")

  endif(BUILD_TESTING)
endif()

if (COMIN_ENABLE_PYTHON_ADAPTER)
  add_subdirectory(python_adapter)
endif()
