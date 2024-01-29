add_library(simple_fortran_plugin MODULE simple_fortran_plugin.F90)
target_link_libraries(simple_fortran_plugin ComIn)

comin_add_test(NAME simple_fortran
  REFERENCE_OUTPUT simple_fortran_reference)
comin_test_add_plugin(TEST simple_fortran
  NAME "simple_fortran_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:simple_fortran_plugin>")
