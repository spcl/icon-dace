find_package(YAXT)

if(YAXT_FOUND)
  add_library(yaxt_c_plugin MODULE yaxt_c_plugin.c)
  target_link_libraries(yaxt_c_plugin ComIn::ComIn YAXT::YAXT_C MPI::MPI_C)

  add_executable(yaxt_c_external yaxt_c_external.c ${CMAKE_SOURCE_DIR}/src/mpi_handshake.c)
  target_link_libraries(yaxt_c_external YAXT::YAXT_C MPI::MPI_C)
  target_include_directories(yaxt_c_external PRIVATE  ${CMAKE_SOURCE_DIR}/include)

  comin_add_test(NAME yaxt_c
    NUM_PROCS 2)
  comin_test_add_plugin(TEST yaxt_c
    NAME yaxt_c_plugin
    PLUGIN_LIBRARY "$<TARGET_FILE:yaxt_c_plugin>"
    COMM "yaxt_c")
  comin_test_add_external_process(TEST yaxt_c
    COMMAND "$<TARGET_FILE:yaxt_c_external> 5126 3")
endif()
