add_library(simple_c_plugin MODULE simple_c_plugin.c)
target_link_libraries(simple_c_plugin ComIn)

comin_add_test(NAME simple_c
  REFERENCE_OUTPUT simple_c_reference)
comin_test_add_plugin(TEST simple_c
  PLUGIN_LIBRARY "$<TARGET_FILE:simple_c_plugin>")

# add a static variant of this plugin to demonstrate static linking of plugins
# The minimal example is linked to this plugin
add_library(simple_c_plugin_static STATIC simple_c_plugin.c)

# HACK: rename `comin_main` to `simple_c_comin_main` to be unique when statically linked
target_compile_definitions(simple_c_plugin_static PRIVATE comin_main=simple_c_comin_main)
target_link_libraries(simple_c_plugin_static ComIn)
