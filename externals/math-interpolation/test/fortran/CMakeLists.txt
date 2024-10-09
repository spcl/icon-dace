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

message(CHECK_START "Fetching external FortUTF")
if("${CMAKE_VERSION}" VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW)
endif()
include(FetchContent)
FetchContent_Declare(
  fortutf
  URL https://github.com/artemis-beta/FortUTF/archive/694afc80a5362d0711a99f248a25b1b5edd668f5.tar.gz
  URL_HASH MD5=1cc90b741c89436b4711b6dd1ea05980)
FetchContent_MakeAvailable(fortutf)
message(CHECK_PASS "done")

set(FORTUTF_PROJECT_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(FORTUTF_PROJECT_SRC_LIBRARY math-interpolation::math-interpolation)
# set(FORTUTF_PROJECT_SRC_FILES helpers.f90)

include("${fortutf_SOURCE_DIR}/cmake/fortutf.cmake")

# We have to prevent FortUTF from generating invalid Fortran code:
set(name_holder ${PROJECT_NAME})
string(REPLACE "-" "_" PROJECT_NAME ${PROJECT_NAME})
# FortUTF_Find_Tests()
set(PROJECT_NAME ${name_holder})
# add_test(NAME FortUTF_UnitTest COMMAND math_interpolation_Tests)
