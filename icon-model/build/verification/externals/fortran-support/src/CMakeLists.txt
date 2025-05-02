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

include(CheckIncludeFiles)
check_include_files("execinfo.h" HAVE_EXECINFO_H)
check_include_files("link.h" HAVE_LINK_H)
check_include_files("unwind.h" HAVE_UNWIND_H)
check_include_files("sys/resource.h" HAVE_GETRUSAGE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/config.h)

if(FS_ENABLE_OMP)
  find_package(OpenMP QUIET)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

add_library(
  fortran-support
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
  mo_exception.F90
  mo_expression.F90
  mo_fortran_tools.F90
  mo_hash_table.F90
  mo_io_units.F90
  mo_namelist.F90
  mo_octree.F90
  mo_random_number_generators.F90
  mo_simple_dump.F90
  mo_util_backtrace.F90
  mo_util_file.F90
  mo_util_libc.F90
  mo_util_nml.F90
  mo_util_rusage.F90
  mo_util_sort.F90
  mo_util_stride.F90
  mo_util_string.F90
  mo_util_string_parse.F90
  mo_util_system.F90
  mo_util_table.F90
  mo_util_texthash.F90
  mo_util_timer.F90
  nml_annotate.c
  util_arithmetic_expr.c
  util_arithmetic_expr.h
  util_backtrace.c
  util_backtrace.h
  util_file.c
  util_file.h
  util_hash.c
  util_hash.h
  util_stride.c
  util_string_parse.c
  util_sysinfo.c
  util_system.c
  util_timer.c)

add_library(${PROJECT_NAME}::fortran-support ALIAS fortran-support)

set(Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(
  fortran-support
  PROPERTIES Fortran_PREPROCESS ON
             Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
             C_STANDARD 99
             C_STANDARD_REQUIRED ON
             EXPORT_NAME ${PROJECT_NAME}::fortran-support)

if(FS_ENABLE_MIXED_PRECISION)
  target_compile_definitions(fortran-support PRIVATE __MIXED_PRECISION)
endif()

include("${PROJECT_SOURCE_DIR}/cmake/check_macro.cmake")
check_macro_defined(HAS_OPENACC_MACRO _OPENACC QUIET)
if(FS_ENABLE_OPENACC)
  # If _OPENACC is defined, assume that the required compiler flags are already
  # provided, e.g. in CMAKE_Fortran_FLAGS:
  if(NOT HAS_OPENACC_MACRO)
    # On LUMI, we only have OpenACC_Fortran, but no OpenACC_C
    find_package(OpenACC QUIET)
    if(NOT OpenACC_Fortran_FOUND)
      message(FATAL_ERROR "Could NOT find OpenACC_Fortran.")
    endif()
    target_compile_options(fortran-support PRIVATE ${OpenACC_Fortran_OPTIONS})
    # This make sures that unit tests (FortUTF) compiles without the need of
    # passing OpenACC compile option.
    target_link_libraries(fortran-support PRIVATE OpenACC::OpenACC_Fortran)
  endif()
else()
  if(HAS_OPENACC_MACRO)
    message(
      FATAL_ERROR
        "The OpenACC support is disabled but the _OPENACC macro is defined")
  endif()
endif()

target_include_directories(
  fortran-support
  PUBLIC
    # Path to the Fortran modules:
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>
    $<INSTALL_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_INSTALL_INCLUDEDIR}>>
  INTERFACE
    # Path to the internal C/C++ headers (for testing): Requires CMake 3.15+ for
    # multiple compile languages
    # https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:C,CXX>:${CMAKE_CURRENT_SOURCE_DIR}>>
  PRIVATE
    # Path to config.h (for C and C++ only): Requires CMake 3.15+ for multiple
    # compile languages
    # https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:C,CXX>:${CMAKE_CURRENT_BINARY_DIR}>>)

install(TARGETS fortran-support EXPORT "${PROJECT_NAME}-targets")

install(DIRECTORY "${Fortran_MODULE_DIRECTORY}/" TYPE INCLUDE)
