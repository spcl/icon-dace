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

project(microphysics_1mom_schemes VERSION 0.1.0 LANGUAGES Fortran)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_OPENACC "Build with OpenACC support" OFF)
option(BUILD_LOOP_EXCHANGE "Build with the loop exchange feature" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
endif()

set(BUILD_TESTING OFF CACHE BOOL "Disable test suites of the dependencies")

set(ICON_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../.." CACHE PATH
  "Path to the root source directory of ICON")

#### Dependencies
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

find_package(NetCDF COMPONENTS Fortran REQUIRED)

if(BUILD_OPENACC)
  find_package(OpenACC REQUIRED)
endif()

find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_subdirectory("${ICON_SOURCE_DIR}/externals/fortran-support" fortran-support)

add_subdirectory("${ICON_SOURCE_DIR}/externals/math-support" math-support)

#### Targets
set(MICROPHYSICS_PATH "${ICON_SOURCE_DIR}/src/granules/microphysics_1mom_schemes")
set(SATAD_PATH "${ICON_SOURCE_DIR}/src/atm_phy_schemes")

# Add the library with the specified files
# (keep in sync with ${ICON_SOURCE_DIR}/.gitlab-ci.yml):
add_library(microphysics_1mom_schemes
  ${SATAD_PATH}/mo_satad.f90
  ${SATAD_PATH}/mo_lookup_tables_constants.f90
  ${MICROPHYSICS_PATH}/gscp_graupel.f90
  ${MICROPHYSICS_PATH}/gscp_data.f90
  ${MICROPHYSICS_PATH}/gscp_kessler.f90
  ${MICROPHYSICS_PATH}/gscp_cloudice.f90
  ${MICROPHYSICS_PATH}/gscp_ice.f90
  ${MICROPHYSICS_PATH}/microphysics_1mom_schemes.f90
  ${ICON_SOURCE_DIR}/src/shared/mo_kind.f90
  ${ICON_SOURCE_DIR}/src/shared/mo_physical_constants.f90
)

target_link_libraries(microphysics_1mom_schemes
  PRIVATE
    fortran-support::fortran-support
    math-support::math-support
)

add_executable(microphysics_1mom_driver
  microphysics_1mom_driver.f90
)

target_link_libraries(microphysics_1mom_driver
  PRIVATE
    microphysics_1mom_schemes
    fortran-support::fortran-support
    NetCDF::NetCDF_Fortran
)

if (BUILD_OPENACC)
  target_compile_options(microphysics_1mom_driver PRIVATE ${OpenACC_Fortran_OPTIONS})
  target_link_libraries(microphysics_1mom_driver PRIVATE OpenACC::OpenACC_Fortran)
  target_compile_options(microphysics_1mom_schemes PRIVATE ${OpenACC_Fortran_OPTIONS})
  target_link_libraries(microphysics_1mom_schemes PRIVATE OpenACC::OpenACC_Fortran)
endif()

if (BUILD_LOOP_EXCHANGE)
  target_compile_definitions(microphysics_1mom_schemes PRIVATE _LOOP_EXCHANGE=1)
endif()

set_target_properties(
  microphysics_1mom_schemes
  microphysics_1mom_driver
  PROPERTIES
    Fortran_PREPROCESS ON
)

#### Custom targets and tests
set(probtest_CMD ${Python3_EXECUTABLE} ${ICON_SOURCE_DIR}/externals/probtest/probtest.py)

enable_testing()

include(ExternalData)
set(remote_file_name "aes-new-gr_moderate-dt30s_atm_3d_ml_20080801T000000Z.nc")
set(ExternalData_URL_TEMPLATES
  "https://swift.dkrz.de/v1/dkrz_f23c4ba9-4a6c-4d70-8739-7bca74b37234/muphys_data/${remote_file_name}")
ExternalData_expand_arguments(fetch_input local_file_name "DATA{input-data.nc}")
ExternalData_Add_Target(fetch_input)
set_property(TARGET fetch_input PROPERTY EXCLUDE_FROM_ALL TRUE)

add_test(NAME fetch_input
  COMMAND "${CMAKE_COMMAND}"
    --build "${CMAKE_CURRENT_BINARY_DIR}"
    --config "$<CONFIG>"
    --target fetch_input
)
set_tests_properties(fetch_input PROPERTIES FIXTURES_SETUP fetch_input)

add_custom_target(tolerance_probtest_init
  COMMAND ${probtest_CMD} init
    --codebase-install ${CMAKE_CURRENT_BINARY_DIR}
    --template-name ${PROJECT_SOURCE_DIR}/granule.jinja
    --experiment-name 1mom
    --member-num 9
  COMMENT "Running probtest initialization for the tolerance data update"
)

add_test(NAME init_probtest
  COMMAND ${probtest_CMD} init
    --codebase-install ${CMAKE_CURRENT_BINARY_DIR}
    --template-name ${PROJECT_SOURCE_DIR}/granule.jinja
    --experiment-name 1mom
)
set_tests_properties(init_probtest PROPERTIES FIXTURES_SETUP init_probtest)

add_custom_target(tolerance_probtest_perturb
  COMMAND ${probtest_CMD} perturb
  COMMENT "Running probtest perturbation for the tolerance data update"
  DEPENDS fetch_input tolerance_probtest_init
)

add_custom_target(tolerance_probtest_ensemble
  COMMAND ${probtest_CMD} run-ensemble
  COMMENT "Running probtest ensemble for the tolerance data update"
  DEPENDS microphysics_1mom_driver tolerance_probtest_perturb
)

set(all_schemes cloudice cloudice2mom graupel kessler)
# define number of tests for each scheme
# (the numbers are hard-coded in the driver)
set(cloudice2mom_num_tests 3)
set(cloudice_num_tests 3)
set(graupel_num_tests 3)
set(kessler_num_tests 1)

set(run_all_schemes_body)
set(tolerance_dependencies)
set(tolerance_ref_files)
set(tolerance_tol_files)
foreach(scheme IN LISTS all_schemes)
  set(run_all_schemes_body "${run_all_schemes_body}
${CMAKE_CURRENT_BINARY_DIR}/microphysics_1mom_driver ${scheme}")

  add_test(NAME run_${scheme}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/microphysics_1mom_driver ${scheme}
  )
  set_tests_properties(run_${scheme}
    PROPERTIES
      FIXTURES_REQUIRED "fetch_input;init_probtest"
      FIXTURES_SETUP run_${scheme}
  )

  foreach(case RANGE 1 ${${scheme}_num_tests})
    add_custom_target(tolerance_probtest_stats_${scheme}_${case}
      COMMAND ${probtest_CMD} stats
        --ensemble
        --file-id "NetCDF" "${scheme}_${case}.nc"
        --stats-file-name "${scheme}_${case}_{member_id}"
      COMMENT "Running probtest stats for ${scheme} (${case}) for the tolerance data update"
      DEPENDS tolerance_probtest_init tolerance_probtest_ensemble
    )

    list(APPEND tolerance_ref_files "${scheme}_${case}_ref")

    add_custom_target(tolerance_probtest_tolerance_${scheme}_${case}
      COMMAND "${CMAKE_COMMAND}" -E make_directory tolerance
      COMMAND ${probtest_CMD} tolerance
        --stats-file-name "${scheme}_${case}_{member_id}"
        --tolerance-file-name "tolerance/${scheme}_${case}"
      COMMENT "Running probtest tolerance for ${scheme} (${case}) for the tolerance data update"
      DEPENDS tolerance_probtest_stats_${scheme}_${case}
    )

    list(APPEND tolerance_tol_files "tolerance/${scheme}_${case}")

    add_custom_target(tolerance_probtest_check_${scheme}_${case}
      COMMAND ${probtest_CMD} check
        --input-file-ref "${scheme}_${case}_ref"
        --input-file-cur "${scheme}_${case}_1"
        --tolerance-file-name "tolerance/${scheme}_${case}"
      COMMENT "Running probtest check for ${scheme} (${case}) for the tolerance data update"
      DEPENDS tolerance_probtest_tolerance_${scheme}_${case}
    )

    list(APPEND tolerance_dependencies tolerance_probtest_tolerance_${scheme}_${case})

    add_test(NAME stats_${scheme}_${case}
      COMMAND ${probtest_CMD} stats
        --no-ensemble
        --file-id "NetCDF" "${scheme}_${case}.nc"
        --stats-file-name "${scheme}_${case}"
    )
    set_tests_properties(stats_${scheme}_${case}
      PROPERTIES
        FIXTURES_REQUIRED "init_probtest;run_${scheme}"
        FIXTURES_SETUP stats_${scheme}_${case}
    )

    add_test(NAME check_${scheme}_${case}
      COMMAND ${probtest_CMD} check
        --input-file-ref "${PROJECT_SOURCE_DIR}/reference/${scheme}_${case}_ref"
        --input-file-cur "${scheme}_${case}"
        --tolerance-file-name "${PROJECT_SOURCE_DIR}/tolerance/${scheme}_${case}"
    )
    set_tests_properties(check_${scheme}_${case}
      PROPERTIES
        FIXTURES_REQUIRED stats_${scheme}_${case}
    )
  endforeach()
endforeach()

add_custom_target(tolerance
  COMMAND "${CMAKE_COMMAND}" -E copy
    ${tolerance_ref_files} ${PROJECT_SOURCE_DIR}/reference
  COMMAND "${CMAKE_COMMAND}" -E copy
    ${tolerance_tol_files} ${PROJECT_SOURCE_DIR}/tolerance
  DEPENDS ${tolerance_dependencies}
)

# We have to support CMake 3.18. Otherwise, we would create the file and assign
# the required permission in-place:
file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/run_all_schemes.sh"
"#!/bin/bash
set -e
${run_all_schemes_body}
")
file(
  COPY "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/run_all_schemes.sh"
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
  FILE_PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)
