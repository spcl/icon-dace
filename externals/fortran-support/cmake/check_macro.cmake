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

# cmake-format: off
# check_macro_defined(<output>
#                     <macro>
#                     [LANG] <lang>
#                     [QUIET])
# cmake-format: on
# -----------------------------------------------------------------------------
# Sets <output> to ON or OFF depending on whether <macro> is already defined via
# CMAKE_${lang}_FLAGS or not. However, this only works when <lang> is Fortran
# and fails otherwise.
#
# <lang> can be set using LANG. <lang> is set as Fortran by default.
#
function(check_macro_defined output macro)

  cmake_parse_arguments(PARSE_ARGV 1 ARG "QUIET" "LANG" "")

  if(NOT ARG_LANG)
    set(ARG_LANG "Fortran")
  endif()

  if(NOT ARG_LANG STREQUAL "Fortran")
    message(FATAL_ERROR "check_macro_defined supports only LANG Fortran")
  endif()

  if(NOT ARG_QUIET)
    message(CHECK_START "Checking whether ${macro} is defined")
  endif()

  # Write a simple Fortran program that checks for ${macro}
  set(check_source_code
      "
      program main
        implicit none
#ifdef ${macro}
        integer a
#else
#endif
#ifndef ${macro}
#else
        integer b
#endif
        a = 4
        b = 2
      end
")

  # Write the Fortran code to a file
  file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90"
       "${check_source_code}")

  # Try to compile the program
  try_compile(result "${CMAKE_BINARY_DIR}"
              "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90")

  if(NOT ARG_QUIET)
    if(${result})
      message(CHECK_PASS "yes")
    else()
      message(CHECK_PASS "no")
    endif()
  endif()

  set(${output}
      ${result}
      PARENT_SCOPE)

endfunction()
