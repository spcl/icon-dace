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

add_library(
  math-support
  mo_gridman_constants.f90
  mo_lib_grid_geometry_info.f90
  mo_lib_loopindices.f90
  mo_math_constants.f90
  mo_math_types.f90
  mo_math_utilities.F90)

add_library(${PROJECT_NAME}::math-support ALIAS math-support)

set(Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(
  math-support PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
                          EXPORT_NAME ${PROJECT_NAME}::math-support)

if(MS_ENABLE_LOOP_EXCHANGE)
  # NITPICKING: can this be set for the Fortran compiler only?
  target_compile_definitions(math-support PRIVATE __LOOP_EXCHANGE)
endif()

if(MS_ENABLE_MIXED_PRECISION)
  target_compile_definitions(math-support PRIVATE __MIXED_PRECISION)
endif()

if(MS_ENABLE_DIM_SWAP)
  target_compile_definitions(math-support PRIVATE __SWAPDIM)
endif()

if(MS_ENABLE_OPENACC)
  # If _OPENACC is defined, assume that the required compiler flags are already
  # provided, e.g. in CMAKE_Fortran_FLAGS:
  if(NOT HAS_OPENACC_MACRO)
    target_compile_options(math-support PRIVATE ${OpenACC_Fortran_OPTIONS})
    # This make sures that unit tests (FortUTF) compiles without the need of
    # passing OpenACC compile option.
    target_link_libraries(math-support PRIVATE OpenACC::OpenACC_Fortran)
  endif()
endif()

target_include_directories(
  math-support
  PUBLIC
    # Path to the Fortran modules:
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>
    $<INSTALL_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_INSTALL_INCLUDEDIR}>>
)

target_link_libraries(math-support PRIVATE fortran-support::fortran-support)

install(TARGETS math-support EXPORT "${PROJECT_NAME}-targets")

install(DIRECTORY "${Fortran_MODULE_DIRECTORY}/" TYPE INCLUDE)
