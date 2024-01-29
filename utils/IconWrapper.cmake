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
