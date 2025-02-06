# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

cmake_minimum_required(VERSION 3.18)

project(
  fortran-support
  VERSION 0.1.0
  LANGUAGES Fortran C)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build tests" ON)

option(FS_ENABLE_BACKTRACE_TEST "Test backtrace function" ON)
option(FS_ENABLE_OMP "Build with OpenMP support" OFF)
option(FS_ENABLE_OPENACC "Build with OpenACC support" OFF)
option(FS_ENABLE_MIXED_PRECISION "Use mixed precision" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE
        STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo"
        FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
endif(NOT CMAKE_BUILD_TYPE)

include(GNUInstallDirs)

add_subdirectory(src)

# Allow for 'make test' even if the tests are disabled:
enable_testing()
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

export(EXPORT "${PROJECT_NAME}-targets"
       FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake")

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)

install(EXPORT "${PROJECT_NAME}-targets"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")

include("${PROJECT_SOURCE_DIR}/cmake/list_sources.cmake")

# The following targets are only added if this is a parent project
if("${CMAKE_PROJECT_NAME}" STREQUAL "${PROJECT_NAME}")
  # Collect source files of the project:
  list_sources(format_sources)
  # Exclude source files of the external projects:
  list(FILTER format_sources EXCLUDE REGEX "/_deps/")
  # Exclude the configuration header file, which does not have the GENERATED
  # property ON because it is automatically generated at the configure time and
  # not at the build time:
  list(REMOVE_ITEM format_sources "${CMAKE_CURRENT_BINARY_DIR}/src/config.h")

  # Exclude source files (not yet automatically) generated with Ragel:
  file(GLOB ragel_sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.rl")
  list(TRANSFORM ragel_sources REPLACE "rl$" "c")
  list(REMOVE_ITEM format_sources ${ragel_sources})

  # Collect C/C++ source files:
  set(c_format_sources "${format_sources}")
  list(FILTER c_format_sources INCLUDE REGEX "\.(c|h|cpp|hpp)$")

  add_custom_target(
    format-c
    COMMAND clang-format --style=file -i ${c_format_sources}
    COMMENT "Format C codes")

  # Collect Fortran source files:
  set(fortran_format_sources "${format_sources}")
  list(FILTER fortran_format_sources INCLUDE REGEX "\.(f|F|f90|F90)$")

  add_custom_target(
    format-fortran
    COMMAND # We do not use a config file but specify the formatting arguments
            # on the command line because the feature requires an extra Python
            # package and does not work for certain arguments even in the latest
            # release (e.g., see
            # https://github.com/pseewald/fprettify/issues/94):
            fprettify --case 2 2 2 1 --indent 2 --whitespace 2 --strip-comments
            ${fortran_format_sources}
    COMMENT "Format Fortran codes")

  add_custom_target(
    format
    DEPENDS format-c format-fortran
    COMMENT "Format codes")
endif()
