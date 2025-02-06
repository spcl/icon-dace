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
  math-support
  VERSION 0.1.0
  LANGUAGES Fortran)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build tests" ON)

option(MS_ENABLE_MIXED_PRECISION "Enable mixed precision" OFF)

# NITPICKING: Do we really need this option? Somebody told that it is needed for
# CPU/GPU but does not have a negative impact on GPU/CPU.
option(MS_ENABLE_LOOP_EXCHANGE "Enable loop exchange" OFF)

# NITPICKING: same here:
option(MS_ENABLE_DIM_SWAP "Enable dimension swap" OFF)

option(MS_ENABLE_OPENACC "Enable OpenACC support" OFF)

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

# The contents of the INTERFACE_INCLUDE_DIRECTORIES of IMPORTED targets are
# treated as SYSTEM includes by default (see
# https://cmake.org/cmake/help/v3.0/prop_tgt/NO_SYSTEM_FROM_IMPORTED.html). This
# happens if CMake knows how to provide system includes to the compiler in use.
# This knowledge for the NVHPC compilers was add in CMake 3.22.0 (see
# https://gitlab.kitware.com/cmake/cmake/-/merge_requests/6691). Unfortunately,
# the knowledge was incomplete because NVHPC Fortran 21.7 and older do not
# support the '-isystem' flag. This was fixed only in CMake 3.27.0 (see
# https://gitlab.kitware.com/cmake/cmake/-/merge_requests/8588). The following
# circumvents the issue for older NVHPC and older CMake:
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" AND CMAKE_Fortran_COMPILER_VERSION
                                                  VERSION_LESS "21.7")
  unset(CMAKE_INCLUDE_SYSTEM_FLAG_Fortran)
endif()

# Dependencies

# avoid duplicate definition of target fortran-support
if(NOT TARGET fortran-support)
  find_package(fortran-support CONFIG)
  if(NOT fortran-support_FOUND)
    include(FetchContent)

    message(CHECK_START "${PROJECT_NAME}: Fetching the fortran-support library")

    FetchContent_Declare(
      libfortran-support
      GIT_REPOSITORY
        git@gitlab.dkrz.de:icon-libraries/libfortran-support.git
      GIT_TAG 1651516a214fe0fa4e5afba1abf58cc777066572)

    set(save_BUILD_TESTING "${BUILD_TESTING}")
    set(BUILD_TESTING
        OFF
        CACHE BOOL "Do not need to build testing for fortran-support" FORCE)

    FetchContent_MakeAvailable(libfortran-support)
    message(CHECK_PASS "done")

    set(BUILD_TESTING
        "${save_BUILD_TESTING}"
        CACHE BOOL "Copy back the original value" FORCE)

  endif()
endif()

include("${PROJECT_SOURCE_DIR}/cmake/check_macro.cmake")
check_macro_defined(HAS_OPENACC_MACRO _OPENACC QUIET)

if(MS_ENABLE_OPENACC)
  # If _OPENACC is defined, assume that the required compiler flags are already
  # provided, e.g. in CMAKE_Fortran_FLAGS:
  if(NOT HAS_OPENACC_MACRO)
    # On LUMI, we only have OpenACC_Fortran, but no OpenACC_C
    find_package(OpenACC QUIET)
    if(NOT OpenACC_Fortran_FOUND)
      message(FATAL_ERROR "Could NOT find OpenACC_Fortran.")
    endif()
  endif()
else()
  if(HAS_OPENACC_MACRO)
    message(
      FATAL_ERROR
        "The OpenACC support is disabled but the _OPENACC macro is defined")
  endif()
endif()

# math-support

add_subdirectory(src)

# Allow for 'make test' even if the tests are disabled:
enable_testing()
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# Install

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO)

export(EXPORT "${PROJECT_NAME}-targets"
       FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake")

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)

install(EXPORT "${PROJECT_NAME}-targets"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")

# The following targets are only added if this is the parent project
if("${CMAKE_PROJECT_NAME}" STREQUAL "${PROJECT_NAME}")

  include("${PROJECT_SOURCE_DIR}/cmake/list_sources.cmake")

  # Collect source files of the project:
  list_sources(format_sources)
  # Exclude source files of the external projects:
  list(FILTER format_sources EXCLUDE REGEX "/_deps/")
  # Collect Fortran source files:
  set(fortran_format_sources "${format_sources}")
  list(FILTER fortran_format_sources INCLUDE REGEX "\.(f|F|f90|F90)$")
  # fprettify fails when formatting long mantissa numbers
  list(REMOVE_ITEM fortran_format_sources
       "${CMAKE_CURRENT_SOURCE_DIR}/src/mo_math_constants.f90")

  add_custom_target(
    format
    COMMAND # We do not use a config file but specify the formatting arguments
            # on the command line because the feature requires an extra Python
            # package and does not work for certain arguments even in the latest
            # release (e.g., see
            # https://github.com/pseewald/fprettify/issues/94):
            fprettify --case 2 2 2 1 --indent 2 --whitespace 2 --strip-comments
            ${fortran_format_sources}
    COMMENT "Format Fortran codes")

endif()
