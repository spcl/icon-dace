cmake_minimum_required(VERSION 3.17)

project(ComIn
  VERSION 0.1.0
  LANGUAGES Fortran C
  DESCRIPTION "Community Interface")

configure_file(
  ${CMAKE_SOURCE_DIR}/include/comin_version.inc.in
  ${CMAKE_SOURCE_DIR}/include/comin_version.inc @ONLY)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
endif(NOT CMAKE_BUILD_TYPE)

include(GNUInstallDirs)
list( APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR} )
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Scripts/")
add_subdirectory(cmake)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build Testing" OFF) # turn testing off by default as it introduces a lot of dependencies
option(COMIN_ENABLE_EXAMPLES "Enable the example plugins" ON)
option(COMIN_ENABLE_PYTHON_ADAPTER "Enable the python adapter" OFF) # introduces python as a dependency
option(COMIN_ENABLE_STATIC_LINKING_TEST "Enable the static linking test" OFF) # may yield problems in non-default environments
include(cmake/Scripts/SetDLLibs.cmake)

include(CTest)
if(NOT BUILD_TESTING)
  add_custom_target(test)
endif()

# Set the Fortran compiler
enable_language(Fortran)
set(Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
find_package(MPI REQUIRED COMPONENTS C Fortran)

set(SOURCE_FILES
  src/mpi_handshake.c
  src/mo_mpi_handshake.F90
  src/comin_c_utils.F90
  src/comin_setup_utils.F90
  src/comin_setup.F90
  src/comin_plugin_interface.F90
  src/comin_parallel_types.F90
  src/comin_parallel.F90
  src/comin_host_interface.F90
  src/comin_errhandler_constants.F90
  src/comin_errhandler_types.F90
  src/comin_errhandler.F90
  src/comin_callback_types.F90
  src/comin_callback.F90
  src/comin_variable_types.F90
  src/comin_metadata.F90
  src/comin_variable.F90
  src/comin_descrdata_types.F90
  src/comin_descrdata.F90
  src/comin_plugin_types.F90
  src/comin_state.F90
  src/comin_setup_constants.F90)


# Create the library
add_library(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran ${COMIN_DL_LIBS})

set_target_properties(${PROJECT_NAME} PROPERTIES
  OUTPUT_NAME "comin"
  SOVERSION "${PROJECT_VERSION_MAJOR}")

target_include_directories(${PROJECT_NAME} PRIVATE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)

set_target_properties(${PROJECT_NAME}
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
  POSITION_INDEPENDENT_CODE 1
)

target_include_directories(${PROJECT_NAME} PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran)
target_include_directories(${PROJECT_NAME} PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILER_INCLUDE_DIRS}>)

include(ComInAddTest)

IF(BUILD_TESTING)
  add_subdirectory(minimal_example)
  add_subdirectory(test)
endif()

install(
  TARGETS ${PROJECT_NAME}
  EXPORT "${PROJECT_NAME}Targets"
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  DIRECTORY
  "${Fortran_MODULE_DIRECTORY}/"
  "include/"
  TYPE INCLUDE
)

export(
  EXPORT "${PROJECT_NAME}Targets"
  FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
  NAMESPACE ComIn::
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

install(
  EXPORT "${PROJECT_NAME}Targets"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake"
  NAMESPACE ComIn::
)

install(
  FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake"
)

# the following lines are necessary to make the `find_package`s in the subdirectories happy
add_library(ComIn::ComIn ALIAS ComIn)
list(APPEND CMAKE_PREFIX_PATH "${PROJECT_BINARY_DIR}")

add_subdirectory(plugins)

add_subdirectory(doc)
