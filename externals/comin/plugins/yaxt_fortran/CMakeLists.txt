find_package(YAXT)

if(YAXT_FOUND)
  add_library(yaxt_fortran_plugin MODULE yaxt_fortran_plugin.F90)
  target_link_libraries(yaxt_fortran_plugin ComIn YAXT::YAXT_Fortran)

  comin_add_test(NAME yaxt_fortran_plugin
    REFERENCE_OUTPUT yaxt_fortran_plugin_reference)
  comin_test_add_plugin(TEST yaxt_fortran_plugin
    PLUGIN_LIBRARY "$<TARGET_FILE:yaxt_fortran_plugin>"
    COMM "yaxt_fortran_plugin_comm")
else()
  message(STATUS "disabling yaxt_fortran_plugin")
endif()
