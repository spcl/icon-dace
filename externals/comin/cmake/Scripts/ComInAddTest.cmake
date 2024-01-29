if(COMMAND comin_add_test)
    return()
endif()

if(NOT MPI_FOUND)
  find_package(MPI)
endif()

if(MPI_FOUND)
  # Check if we are using OpenMPI
  set(COMIN_TEST_ENABLE_REFERENCE_CHECKS TRUE CACHE BOOL "Enable reference checking for comin tests")
  if(COMIN_TEST_ENABLE_REFERENCE_CHECKS)
    if(MPI_C_LIBRARY_VERSION_STRING MATCHES "Open MPI v[2-9]\\.[0-9]\\.[0-9]")
      set(COMIN_IS_OPENMPI TRUE CACHE INTERNAL "Do we use OpenMPI?")
      message(STATUS "OpenMPI found! Test output is compared against references!")
    else()
      message(WARNING "You are not using OpenMPI. Test output cannot be compared against reference! Disabling COMIN_TEST_ENABLE_REFERENCE_CHECKS")
      set(COMIN_TEST_ENABLE_REFERENCE_CHECKS FALSE CACHE BOOL "Enable reference checking for comin tests" FORCE)
    endif()
  endif()
else()
  message(WARNING "No MPI found. For using the ComIn testing facilities MPI is needed.")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND COMIN_TEST_ENABLE_REFERENCE_CHECKS)
  message(WARNING "You are not using GFortran. Test output cannot be compared against reference! Disabling COMIN_TEST_ENABLE_REFERENCE_CHECKS")
  set(COMIN_TEST_ENABLE_REFERENCE_CHECKS FALSE CACHE BOOL "Enable reference checking for comin tests" FORCE)
endif()

add_custom_target(update_test_references
  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red --bold "CAUTION: Please check if the updated references contain only changes that are intended!")

# DOCUMENTATION CAN BE FOUND IN doc/cmake.md
# PLEASE DON'T FORGET TO KEEP THEM UP TO DATE!
function(comin_add_test)
  if(NOT MPI_FOUND)
    return()
  endif()
  set(options)
  set(oneValueArgs NAME NUM_PROCS REFERENCE_OUTPUT)
  set(multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT ARGS_NUM_PROCS)
    set(ARGS_NUM_PROCS 1)
  endif()

  set(test_dir ${CMAKE_CURRENT_BINARY_DIR}/${ARGS_NAME})
  file(MAKE_DIRECTORY ${test_dir})

  add_test(NAME ${ARGS_NAME}
    COMMAND sh ${test_dir}/run.sh
    WORKING_DIRECTORY ${test_dir})

  # we have created the test target
  if(NOT BUILD_TESTING)
    return()
  endif()

  # add a dummy target to collect properties
  add_custom_target(_comin_test_${ARGS_NAME})
  set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY PLUGIN_COUNT 0)
  set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY NAMELIST "&comin_nml\n")
  set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY NUM_PROCS ${ARGS_NUM_PROCS})
  set_property(TARGET _comin_test_${ARGS_NAME} PROPERTY EXTERNAL_PROCESSES "")

  # generate run.sh
  set(run_sh "#!/bin/sh\n")
  string(APPEND run_sh
    "cd ${test_dir}\n"
    "set -e\n"
    "ln -f -s $<TARGET_FILE_DIR:ComIn::minimal_example>/grids\n"
    "${MPIEXEC}")
  if(${COMIN_IS_OPENMPI})
    string(APPEND run_sh " --output-filename ${test_dir}/output ")
  endif()

  string(APPEND run_sh " ${MPIEXEC_NUMPROC_FLAG} $<TARGET_PROPERTY:_comin_test_${ARGS_NAME},NUM_PROCS> ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ComIn::minimal_example> ${MPIEXEC_POSTFLAGS} ./grids")

  # add external processes
  string(APPEND run_sh "$<TARGET_GENEX_EVAL:_comin_test_${ARGS_NAME},$<TARGET_PROPERTY:_comin_test_${ARGS_NAME},EXTERNAL_PROCESSES>> \n")

  # check the references
  if(ARGS_REFERENCE_OUTPUT AND COMIN_TEST_ENABLE_REFERENCE_CHECKS)
    get_filename_component(ref_path ${ARGS_REFERENCE_OUTPUT} REALPATH BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
    string(APPEND run_sh "diff -r ${test_dir}/output ${ref_path} \n")

    add_custom_target(update_test_references_${ARGS_NAME})
    add_custom_command(TARGET update_test_references_${ARGS_NAME}
      COMMAND ${CMAKE_COMMAND} -E copy_directory ${test_dir}/output ${ref_path})
    add_dependencies(update_test_references update_test_references_${ARGS_NAME})
  endif()

  file(GENERATE OUTPUT ${test_dir}/run.sh
    CONTENT ${run_sh})

  # generate master.nml
  file(GENERATE OUTPUT ${test_dir}/master.nml
    CONTENT "$<TARGET_GENEX_EVAL:_comin_test_${ARGS_NAME},$<TARGET_PROPERTY:_comin_test_${ARGS_NAME},NAMELIST>>/\n")

endfunction()

# DOCUMENTATION CAN BE FOUND IN doc/cmake.md
# PLEASE DON'T FORGET TO KEEP THEM UP TO DATE!
function(comin_test_add_plugin)
  if(NOT BUILD_TESTING)
    return()
  endif()
  set(options)
  set(oneValueArgs TEST NAME PLUGIN_LIBRARY PRIMARY_CONSTRUCTOR OPTIONS COMM)
  set(multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(ARGS_UNPARSED_ARGUMENTS)
    message(WARNING "Unknown arguments: ${ARGS_UNPARSED_ARGUMENTS}")
  endif()
  if(NOT ARGS_TEST)
    message(FATAL_ERROR "No test provided")
  endif()

  if( NOT ARGS_PLUGIN_LIBRARY AND NOT ARGS_PRIMARY_CONSTRUCTOR)
    message(FATAL_ERROR "Either PLUGIN_LIBRARY or PRIMARY_CONSTRUCTOR must be provided")
  endif()

  get_property(PLUGIN_COUNT TARGET _comin_test_${ARGS_TEST} PROPERTY PLUGIN_COUNT)
  MATH(EXPR PLUGIN_COUNT "${PLUGIN_COUNT}+1")
  set_property(TARGET _comin_test_${ARGS_TEST} PROPERTY PLUGIN_COUNT ${PLUGIN_COUNT})
  set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
    "  plugin_list(${PLUGIN_COUNT})%name = \"${ARGS_NAME}\"\n")
  if(ARGS_PLUGIN_LIBRARY) # optional
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%plugin_library = \"${ARGS_PLUGIN_LIBRARY}\"\n")
  endif()
  if(ARGS_PRIMARY_CONSTRUCTOR)
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%primary_constructor = \"${ARGS_PRIMARY_CONSTRUCTOR}\"\n")
  endif()
  if(ARGS_OPTIONS) # optional
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%options = \"${ARGS_OPTIONS}\"\n")
  endif()
  if(ARGS_COMM) # optional
    set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY NAMELIST
      "  plugin_list(${PLUGIN_COUNT})%comm = \"${ARGS_COMM}\"\n")
  endif()
endfunction()

# DOCUMENTATION CAN BE FOUND IN doc/cmake.md
# PLEASE DON'T FORGET TO KEEP THEM UP TO DATE!
function(comin_test_add_external_process)
  if(NOT BUILD_TESTING)
    return()
  endif()
  set(options)
  set(oneValueArgs TEST NUM_PROCS COMMAND)
  set(multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(ARGS_UNPARSED_ARGUMENTS)
    message(WARNING "Unknown arguments: ${ARGS_UNPARSED_ARGUMENTS}")
  endif()
  if(NOT ARGS_TEST)
    message(FATAL_ERROR "No test provided")
  endif()
  if(NOT ARGS_COMMAND)
    message(FATAL_ERROR "No command provided")
  endif()
  if(NOT ARGS_NUM_PROCS)
    set(ARGS_NUM_PROCS 1)
  endif()
  set_property(TARGET _comin_test_${ARGS_TEST} APPEND_STRING PROPERTY EXTERNAL_PROCESSES
    "  : ${MPIEXEC_NUMPROC_FLAG} ${ARGS_NUM_PROCS} ${ARGS_COMMAND}")
endfunction()
