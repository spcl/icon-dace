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
# add_icon_c_test(<test>
#                 <source>)
# cmake-format: on
# -----------------------------------------------------------------------------
# Compiles a test executable with the name <test> using the source code
# <source>. The googletest and fortran-support libraries will be linked
# automatically.
#
# The C++ standard is set to C++17.
#
macro(add_icon_c_test test_name file_name)
  add_executable("CTest_${test_name}" ${file_name})
  target_link_libraries(
    "CTest_${test_name}" PRIVATE fortran-support::fortran-support
                                 GTest::gtest_main stdc++fs)
  add_test(NAME "CTest_${test_name}" COMMAND "CTest_${test_name}")
  set_property(TEST "CTest_${test_name}" PROPERTY LABELS C)
  set_target_properties("CTest_${test_name}"
                        PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
endmacro()
