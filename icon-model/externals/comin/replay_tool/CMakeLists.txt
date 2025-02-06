find_package(YAC)
find_package(NetCDF REQUIRED COMPONENTS C Fortran)

add_subdirectory(src)

add_library(comin_run_recorder_plugin MODULE comin_run_recorder_plugin.F90 src/comin_descrdata_save.F90)
target_link_libraries(comin_run_recorder_plugin MPI::MPI_Fortran ComIn NetCDF::NetCDF_Fortran utils)

add_executable(comin_replay comin_replay.F90 src/comin_descrdata_load.F90)
add_executable(${PROJECT_NAME}::comin_replay ALIAS comin_replay)

set_target_properties(comin_replay PROPERTIES
  EXPORT_NAME "${PROJECT_NAME}::comin_replay"
)

if(YAC_FOUND)
  # We need to link YAC before NetCDF because CMake does not resolve
  # the linking order correctly as both YAC and NetCDF_Fortran depend
  # on NetCDF_C
  target_link_libraries(comin_replay YAC::YAC)
  target_compile_definitions(comin_replay PRIVATE -DENABLE_YAC)
endif()
target_link_libraries(comin_replay MPI::MPI_Fortran ComIn NetCDF::NetCDF_Fortran utils)

install(TARGETS comin_replay EXPORT "${PROJECT_NAME}Targets")

add_library(comin_var_recorder_plugin MODULE comin_var_recorder_plugin.F90)
target_link_libraries(comin_var_recorder_plugin ComIn NetCDF::NetCDF_Fortran utils)

add_library(comin_var_replay_plugin MODULE comin_var_replay_plugin.F90)
add_library(${PROJECT_NAME}::comin_var_replay_plugin ALIAS comin_var_replay_plugin)

set_target_properties(comin_var_replay_plugin PROPERTIES
  EXPORT_NAME "${PROJECT_NAME}::comin_var_replay_plugin"
)

target_link_libraries(comin_var_replay_plugin
  PUBLIC ComIn NetCDF::NetCDF_Fortran
  PRIVATE utils
)

install(
  TARGETS comin_var_replay_plugin
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

### testing ###
# record a run of the minimal example and replay it

if(BUILD_TESTING)

comin_add_test(NAME record
  NUM_PROCS 2)
comin_test_add_plugin(TEST record
  NAME "comin_run_recorder_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:comin_run_recorder_plugin>"
)
comin_test_add_plugin(TEST record
  NAME comin_var_recorder_plugin
  PLUGIN_LIBRARY "$<TARGET_FILE:comin_var_recorder_plugin>"
  OPTIONS "pres,u,v,qr,qi,qv,qs,qc")

set_tests_properties(record PROPERTIES
  FIXTURES_SETUP minimal_example_record
)
get_test_property(record WORKING_DIRECTORY MINIMAL_EXAMPLE_RECORD_WD)

comin_add_replay_test(NAME replay
  REPLAY_DATA_PATH "${MINIMAL_EXAMPLE_RECORD_WD}/"
  NUM_PROCS 2
)
comin_test_add_plugin(TEST replay
  NAME "comin_run_recorder_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:comin_run_recorder_plugin>"
)
comin_test_add_plugin(TEST replay
  NAME comin_var_recorder_plugin
  PLUGIN_LIBRARY "$<TARGET_FILE:comin_var_recorder_plugin>"
  OPTIONS "pres,u,v,qr,qi,qv,qs,qc")

set_tests_properties(replay PROPERTIES
  FIXTURES_REQUIRED minimal_example_record
)
set_tests_properties(replay PROPERTIES
  FIXTURES_SETUP replay_record
)
get_test_property(replay WORKING_DIRECTORY REPLAY_RECORD_WD)

add_test(NAME record_replay_check
  COMMAND bash -c
  "diff <(${NetCDF_NCDUMP} ${MINIMAL_EXAMPLE_RECORD_WD}/0.nc) <(${NetCDF_NCDUMP} ${REPLAY_RECORD_WD}/0.nc) && \
     diff <(${NetCDF_NCDUMP} ${MINIMAL_EXAMPLE_RECORD_WD}/1.nc) <(${NetCDF_NCDUMP} ${REPLAY_RECORD_WD}/1.nc)"
)
set_tests_properties(record_replay_check PROPERTIES
  FIXTURES_REQUIRED "minimal_example_record;replay_record"
)

# compare the variable records.
# ignore pos_jn. replay can not ensure that it is the same.
add_test(NAME var_record_replay_check
  COMMAND bash -c
  "diff <(${NetCDF_NCDUMP} ${MINIMAL_EXAMPLE_RECORD_WD}/vars_0.nc | grep -v pos_jn) <(${NetCDF_NCDUMP} ${REPLAY_RECORD_WD}/vars_0.nc | grep -v pos_jn) && \
     diff <(${NetCDF_NCDUMP} ${MINIMAL_EXAMPLE_RECORD_WD}/vars_1.nc | grep -v pos_jn) <(${NetCDF_NCDUMP} ${REPLAY_RECORD_WD}/vars_1.nc | grep -v pos_jn)"
)
set_tests_properties(var_record_replay_check PROPERTIES
  FIXTURES_REQUIRED "minimal_example_record;replay_record"
)

endif(BUILD_TESTING)
