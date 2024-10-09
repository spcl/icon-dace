find_package(YAC)

# link example directory to make relative references to the python files
file(CREATE_LINK ${PROJECT_SOURCE_DIR}/examples ${CMAKE_CURRENT_BINARY_DIR}/examples SYMBOLIC)

comin_add_test(NAME simple_python_plugin
  REFERENCE_OUTPUT simple_python_plugin_reference)
comin_test_add_plugin(TEST simple_python_plugin
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../examples/simple_python_plugin.py --hello=world 'argument with spaces'")

if( YAC_FOUND )
  comin_add_replay_test(NAME yac_python_example
    REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
    NUM_PROCS 3
    REFERENCE_OUTPUT yac_python_example_reference
  )
  comin_test_add_plugin(TEST yac_python_example
    NAME "python_adapter"
    PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
    OPTIONS "../examples/yac.py")
endif()

if( NOT MPI4PY_FOUND )
  # check if a working mpi4py module exists
  execute_process(COMMAND
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${Python3_EXECUTABLE} -c
    "from mpi4py import MPI; print(MPI.COMM_WORLD.size if MPI.COMM_WORLD.rank == 0 else '')"
    OUTPUT_VARIABLE mpi4py_config_test_output
    RESULT_VARIABLE mpi4py_config_test_result
    ERROR_VARIABLE mpi4py_config_test_error
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(mpi4py_config_test_result)
    message(WARNING "mpi4py is not working properly. Tests with mpi4py are disabled")
    set(MPI4PY_FOUND FALSE)
  else()
    if(mpi4py_config_test_output MATCHES "2")
      set(MPI4PY_FOUND TRUE)
    else()
      message(WARNING "mpi4py is not working properly. Tests with mpi4py are disabled")
      set(MPI4PY_FOUND FALSE)
    endif()
  endif()
endif()

if(MPI4PY_FOUND)
  file(CREATE_LINK ${PROJECT_SOURCE_DIR}/test/test_mpi.py ${CMAKE_CURRENT_BINARY_DIR}/test_mpi.py SYMBOLIC)
  comin_add_test(NAME python_test_mpi
    NUM_PROCS 2
    REFERENCE_OUTPUT test_mpi_reference)
  comin_test_add_plugin(TEST python_test_mpi
    NAME "python_adapter"
    PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
    OPTIONS "../test_mpi.py"
    COMM "python_mpi_test_comm")
  comin_test_add_external_process(TEST python_test_mpi
    COMMAND "$<TARGET_FILE:mpi_communicator_test_external> python_mpi_test_comm")
endif()

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/variables.py" "${CMAKE_CURRENT_BINARY_DIR}/variables.py" SYMBOLIC)
comin_add_replay_test(NAME python_variables
  REPLAY_DATA "test_nwp_R02B04_R02B05_nest"
  NUM_PROCS 3
  REFERENCE_OUTPUT python_variables_reference)
comin_test_add_plugin(TEST python_variables
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../variables.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/descrdata.py" "${CMAKE_CURRENT_BINARY_DIR}/descrdata.py" SYMBOLIC)
comin_add_test(NAME python_descrdata
  REFERENCE_OUTPUT python_descrdata_reference)
comin_test_add_plugin(TEST python_descrdata
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../descrdata.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/generate_python_api_doc.py" "${CMAKE_CURRENT_BINARY_DIR}/generate_python_api_doc.py" SYMBOLIC)
comin_add_test(NAME python_generate_python_api_doc
  REFERENCE_OUTPUT python_generate_python_api_doc_reference)
comin_test_add_plugin(TEST python_generate_python_api_doc
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../generate_python_api_doc.py")

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/finish.py" "${CMAKE_CURRENT_BINARY_DIR}/finish.py" SYMBOLIC)
comin_add_test(NAME python_finish)
comin_test_add_plugin(TEST python_finish
  NAME "python_adapter"
  PLUGIN_LIBRARY $<TARGET_FILE:python_adapter>
  OPTIONS "../finish.py")
set_tests_properties(python_finish PROPERTIES WILL_FAIL ON)
