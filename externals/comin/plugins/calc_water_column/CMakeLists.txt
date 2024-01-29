add_library(calc_water_column_plugin MODULE calc_water_column_plugin.F90)
target_link_libraries(calc_water_column_plugin ComIn)

comin_add_test(NAME calc_water_column
  REFERENCE_OUTPUT calc_water_column_reference)
comin_test_add_plugin(TEST calc_water_column
  NAME "calc_water_column_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:calc_water_column_plugin>")
