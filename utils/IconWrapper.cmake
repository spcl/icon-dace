# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# This is a wrapper file that simplifies the integration of CMake-based bundled
# libraries into ICON build system.

# Do not impose any additional version constraints:
cmake_minimum_required(VERSION
  ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION})

# Do not impose any additional language constraints:
project(IconWrapper LANGUAGES NONE)

# Set the path to the root source directory of the bundled library:
set(ICON_DEP_SOURCE_DIR "" CACHE PATH
  "Path to the source directory of the ICON dependency")
if(NOT ICON_DEP_SOURCE_DIR)
  # By default, the root source directory of the bundled library is expected to
  # be the parent directory of the current build directory:
  get_filename_component(ICON_DEP_SOURCE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}" DIRECTORY)
endif()

# Set the path to the root build directory of the bundled library:
set(dep_build_dir "_")

add_subdirectory("${ICON_DEP_SOURCE_DIR}" "${dep_build_dir}")

# Allow for 'make test' even if the bundled library does not provide any tests
# or the tests are disabled:
enable_testing()

# Set the list of libraries (targets) that ICON depends on:
set(ICON_DEP_LIBRARIES "" CACHE STRING
  "List of libraries (targets) that ICON needs to be linked to")

# Optionally work around the fact the CMake does not recognize the NEC compiler
# family and requires some help in figuring out the implicit linker flags of
# the compilers:
set(ICON_TRY_NEC_WORKAROUNDS ON CACHE BOOL
  "Help CMake recognize NEC compilers and their implicit linker flags")

if(ICON_DEP_LIBRARIES)
  # If the list of dependencies is not empty, we activate the linker flag
  # extraction mechanism. First, all languages that are enabled for the bundled
  # library must be enabled for the current project:
  get_property(dep_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  enable_language(${dep_languages})

  # Make sure that Fortran is enabled:
  if(NOT "Fortran" IN_LIST dep_languages)
    enable_language(Fortran)
  endif()

  if(ICON_TRY_NEC_WORKAROUNDS)
    # Apply the workarounds to the currently enabled languages:
    get_property(all_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
    set(nec_languages)

    # The users might be confused with the NEC-related messages when they do not
    # use NEC compilers:
    set(CMAKE_REQUIRED_QUIET 1)

    # CMake is expected to mistake the NEC compiler for the GNU one:
    if("C" IN_LIST all_languages AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
      include(CheckSymbolExists)
      check_symbol_exists(__NEC__ "" ICON_C_is_NEC)
      if(ICON_C_is_NEC)
        set(abi_C_file "${CMAKE_ROOT}/Modules/CMakeCCompilerABI.c")
        list(APPEND nec_languages "C")
      endif()
    endif()

    # CMake is expected to mistake the NEC compiler for the GNU one:
    if("CXX" IN_LIST all_languages AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      include(CheckCXXSymbolExists)
      check_cxx_symbol_exists(__NEC__ "" ICON_CXX_is_NEC)
      if(ICON_CXX_is_NEC)
        set(abi_CXX_file "${CMAKE_ROOT}/Modules/CMakeCXXCompilerABI.cpp")
        list(APPEND nec_languages "CXX")
      endif()
    endif()

    # CMake is expected to mistake the NEC C compiler for the GNU one:
    if("Fortran" IN_LIST all_languages AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      file(WRITE "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90"
        "      program main
      implicit none
#ifdef __NEC__
      integer a
#else
      choke me
#endif
#ifndef __NEC__
      choke me
#else
      integer b
#endif
      a = 4
      b = 2
      end
")
      try_compile(ICON_Fortran_is_NEC "${PROJECT_BINARY_DIR}"
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90")
      if(ICON_Fortran_is_NEC)
        set(abi_Fortran_file "${CMAKE_ROOT}/Modules/CMakeFortranCompilerABI.F")
        list(APPEND nec_languages "Fortran")
      endif()
    endif()

    # Re-run the checks for internal compiler flags:
    if(nec_languages)
      message(STATUS "Applying the NEC workarounds")
      include(CMakeDetermineCompilerABI)
      # NEC compilers use 'nld' instead of 'ld' for linking:
      set(CMAKE_LINK_STARTFILE "nld")
      foreach(language IN LISTS nec_languages)
        message(STATUS "The real ${language} compiler ID is NEC")
        unset(CMAKE_${language}_ABI_COMPILED)
        CMAKE_DETERMINE_COMPILER_ABI(${language} "${abi_${language}_file}")
      endforeach()
    endif()
  endif()

  # Create a dummy executable (we have to set a source file to avoid complains
  # from CMake, so we give it this file):
  add_executable(IconFlags EXCLUDE_FROM_ALL "${CMAKE_CURRENT_LIST_FILE}")

  # The executable depends on libraries that ICON needs to be linked to:
  target_link_libraries(IconFlags PRIVATE "${ICON_DEP_LIBRARIES}")

  # To be able to extract all linker flags, we create a custom language:
  set_target_properties(IconFlags PROPERTIES LINKER_LANGUAGE IconFlagLanguage)

  # Define a list of variables that the custom language should inherit from
  # Fortran language to make the <LINK_LIBRARIES> below complete:
  set(Fortran_variables
    CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES
    CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)

  # Append the RPATH flags if requested:
  set(ICON_ENABLE_RPATHS ON CACHE BOOL "Enable runtime library search paths")
  if(ICON_ENABLE_RPATHS)
    list(APPEND Fortran_variables CMAKE_EXECUTABLE_RUNTIME_Fortran_FLAG)
  endif()

  # Copy the variables from Fortran to the custom language:
  foreach(Fortran_variable IN LISTS Fortran_variables)
    string(REPLACE Fortran IconFlagLanguage
      IconFlagLanguage_variable ${Fortran_variable})
    set(${IconFlagLanguage_variable} ${${Fortran_variable}})
  endforeach()

  # Set the custom language linker command so that when the resulting link.txt
  # file is executed with a shell interpreter, it reports the libraries that
  # ICON needs to be linked to:
  set(CMAKE_IconFlagLanguage_LINK_EXECUTABLE
    "<CMAKE_COMMAND> -E echo '<LINK_LIBRARIES>'")
endif()
