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

enable_language(CXX)

find_package(GTest 1.13)
if(NOT GTest_FOUND)
  message(CHECK_START "Fetching external GTest")
  if("${CMAKE_VERSION}" VERSION_GREATER_EQUAL "3.24")
    cmake_policy(SET CMP0135 NEW)
  endif()
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.tar.gz
    URL_HASH MD5=95b29f0038ec84a611df951d74d99897)
  set(INSTALL_GTEST
      OFF
      CACHE BOOL "Disable installation of googletest")
  set(BUILD_GMOCK
      OFF
      CACHE BOOL "Do not build the googlemock subproject")
  FetchContent_MakeAvailable(googletest)
  # GTest puts its build artifacts to ${CMAKE_BINARY_DIR} instead of
  # ${CMAKE_CURRENT_BINARY_DIR}. Therefore, we have to make it behave:
  set_target_properties(
    gtest_main gtest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${googletest_BINARY_DIR}/bin"
               LIBRARY_OUTPUT_DIRECTORY "${googletest_BINARY_DIR}/lib"
               ARCHIVE_OUTPUT_DIRECTORY "${googletest_BINARY_DIR}/lib"
               PDB_OUTPUT_DIRECTORY "${googletest_BINARY_DIR}/bin"
               COMPILE_PDB_OUTPUT_DIRECTORY "${googletest_BINARY_DIR}/lib")
  message(CHECK_PASS "done")
endif()

include("${PROJECT_SOURCE_DIR}/cmake/gtest_helper.cmake")

if(FS_ENABLE_BACKTRACE_TEST)
  add_icon_c_test(UtilBacktraceTest ctest_util_backtrace.cpp)
endif()

file(
  WRITE "${CMAKE_CURRENT_BINARY_DIR}/test.namelist"
  "&info_nml\n"
  " file_name = 'test.namelist'\n"
  " keywords  = 'unit test','C language','nml_annotate'\n"
  "/\n"
  "&test_nml\n"
  " numbers         = 0,1,2,3\n"
  " floating_points = -7.1,1.0,3.7,8.9\n"
  " booleans        = .TRUE., .FALSE., .true.\n"
  " chars           = 'unit', 'test'\n"
  "/\n")
file(
  WRITE "${CMAKE_CURRENT_BINARY_DIR}/result.namelist"
  "\n"
  "NAMELIST info_nml\n"
  "    file_name                                   'test.namelist'\n"
  "    keywords                                    'unit test', 'C language', 'nml_annotate'\n"
  "\n"
  "NAMELIST test_nml\n"
  "    numbers                                     0, 1, 2, 3\n"
  "    floating_points                             -7.1, 1.0, 3.7, 8.9\n"
  "    booleans                                    .TRUE., .FALSE., .true.\n"
  "    chars                                       'unit', 'test'\n")
add_icon_c_test(UtilNmlAnnotateTest ctest_nml_annotate.cpp)
add_icon_c_test(UtilArithmeticExprTest ctest_util_arithmetic_expr.cpp)
add_icon_c_test(UtilHashTest ctest_util_hash.cpp)
add_icon_c_test(UtilStrideTest ctest_util_stride.cpp)
add_icon_c_test(UtilStringParseTest ctest_util_string_parse.cpp)
add_icon_c_test(UtilSysinfoTest ctest_util_sysinfo.cpp)
add_icon_c_test(UtilSystemTest ctest_util_system.cpp)
add_icon_c_test(UtilTimerTest ctest_util_timer.cpp)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/util_file_test.txt"
     "This is a test file for unit tests for util_file.c\n")
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/util_file_test.txt"
     "${CMAKE_CURRENT_BINARY_DIR}/util_file_link.txt" SYMBOLIC)
add_icon_c_test(UtilFileTest ctest_util_file.cpp)
