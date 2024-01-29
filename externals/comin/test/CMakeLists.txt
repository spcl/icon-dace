add_executable(mpi_communicator_test_external
  mpi_communicator_test_external.c
  ${CMAKE_SOURCE_DIR}/src/mpi_handshake.c)
target_link_libraries(mpi_communicator_test_external MPI::MPI_C)

add_library(mpi_communicator_test_plugin MODULE mpi_communicator_test.c)
target_link_libraries(mpi_communicator_test_plugin ComIn::ComIn MPI::MPI_C)

comin_add_test(NAME mpi_communicator
  NUM_PROCS 2
  REFERENCE_OUTPUT "mpi_communicator_test_reference"
)
comin_test_add_plugin(TEST mpi_communicator
  NAME "mpi_communicator_test_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:mpi_communicator_test_plugin>"
  PRIMARY_CONSTRUCTOR "mpi_communicator_test_setup"
  COMM "mpi_communicator_test_comm")
comin_test_add_external_process(TEST mpi_communicator
  NUM_PROCS 1
  COMMAND "$<TARGET_FILE:mpi_communicator_test_external> mpi_communicator_test_comm"
)

comin_add_test(NAME mpi_communicator2
  NUM_PROCS 2
  EXTERNAL_PROCESSES
  REFERENCE_OUTPUT "mpi_communicator_test2_reference"
)
comin_test_add_plugin(TEST mpi_communicator2
  NAME "mpi_communicator_test_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:mpi_communicator_test_plugin>"
  PRIMARY_CONSTRUCTOR "mpi_communicator_test_setup"
  COMM "mpi_communicator_test_comm")
comin_test_add_plugin(TEST mpi_communicator2
  NAME "mpi_communicator_test_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:mpi_communicator_test_plugin>"
  PRIMARY_CONSTRUCTOR "mpi_communicator_test_setup"
  COMM "mpi_communicator_test_comm2")
comin_test_add_external_process(TEST mpi_communicator2
  NUM_PROCS 1
  COMMAND "$<TARGET_FILE:mpi_communicator_test_external> mpi_communicator_test_comm"
)
comin_test_add_external_process(TEST mpi_communicator2
  NUM_PROCS 2
  COMMAND "$<TARGET_FILE:mpi_communicator_test_external> mpi_communicator_test_comm2"
)


add_library(test_parallel_plugin MODULE test_parallel_plugin.F90)
target_link_libraries(test_parallel_plugin ComIn MPI::MPI_Fortran)

comin_add_test(NAME parallel
  NUM_PROCS 2
  REFERENCE_OUTPUT parallel_reference)
comin_test_add_plugin(TEST parallel
  NAME "test_parallel_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:test_parallel_plugin>"
  PRIMARY_CONSTRUCTOR "test_parallel_plugin_setup")

add_library(finish_test MODULE finish_test.c)
target_link_libraries(finish_test ComIn::ComIn)

comin_add_test(NAME finish_test)
  comin_test_add_plugin(TEST finish_test
    NAME "finish_test"
    PLUGIN_LIBRARY "$<TARGET_FILE:finish_test>"
    PRIMARY_CONSTRUCTOR "setup")
  set_tests_properties(finish_test PROPERTIES
    PASS_REGULAR_EXPRESSION "This is a dummy error to test comin_plugin_finish!"
  )

add_executable(test_ep_names test_ep_names.c)
target_link_libraries(test_ep_names ComIn)
add_test(test_ep_names test_ep_names)

add_library(test_lmodexclusive_plugin MODULE test_lmodexclusive_plugin.c)
target_link_libraries(test_lmodexclusive_plugin ComIn)
comin_add_test(NAME test_lmodexclusive
  REFERENCE_OUTPUT test_lmodexclusive_reference)
comin_test_add_plugin(TEST test_lmodexclusive
  NAME "test_lmodexclusive_plugin1"
  PLUGIN_LIBRARY "$<TARGET_FILE:test_lmodexclusive_plugin>")
comin_test_add_plugin(TEST test_lmodexclusive
  NAME "test_lmodexclusive_plugin2"
  PLUGIN_LIBRARY "$<TARGET_FILE:test_lmodexclusive_plugin>")

add_library(test_errorcodes_plugin MODULE test_errorcodes_plugin.c)
target_link_libraries(test_errorcodes_plugin ComIn)
comin_add_test(NAME test_errorcodes
  REFERENCE_OUTPUT test_errorcodes_reference)
comin_test_add_plugin(TEST test_errorcodes
  NAME "test_errorcodes_plugin"
  PLUGIN_LIBRARY "$<TARGET_FILE:test_errorcodes_plugin>")

add_library(test_index_mapping_plugin MODULE test_index_mapping_plugin.c)
target_link_libraries(test_index_mapping_plugin ComIn)
comin_add_test(NAME test_index_mapping
  NUM_PROCS 2
  REFERENCE_OUTPUT test_index_mapping_reference)
comin_test_add_plugin(TEST test_index_mapping
  NAME "test_index_mapping"
  PLUGIN_LIBRARY "$<TARGET_FILE:test_index_mapping_plugin>")

if(COMIN_ENABLE_STATIC_LINKING_TEST)
  add_library(static_plugin STATIC static_shared_test_plugin.c)
  target_link_libraries(static_plugin ComIn)
  target_compile_definitions(static_plugin PRIVATE -DCOMIN_STATIC_PLUGIN_IS_STATIC)

  add_library(shared_plugin MODULE static_shared_test_plugin.c)
  target_link_libraries(shared_plugin ComIn)

  # link static_plugin into minimal example and ensure that the needed symbols
  # are exported to the dynamic symbol table
  target_link_libraries(minimal_example
    -Wl,--dynamic-list=${CMAKE_CURRENT_SOURCE_DIR}/dynamic_list.sym
    -Wl,-ustatic_linking_test_comin_main
    static_plugin)

  comin_add_test(NAME static_linking_test
    REFERENCE_OUTPUT static_linking_test_reference)
  comin_test_add_plugin(TEST static_linking_test
    NAME "static_plugin"
    PRIMARY_CONSTRUCTOR "static_linking_test_comin_main")
  comin_test_add_plugin(TEST static_linking_test
    NAME "shared_plugin"
    PLUGIN_LIBRARY "$<TARGET_FILE:shared_plugin>"
    PRIMARY_CONSTRUCTOR "static_linking_test_comin_main")
endif()
