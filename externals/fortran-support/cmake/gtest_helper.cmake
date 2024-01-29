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

macro(add_icon_c_test test_name file_names)
    add_executable("CTest_${test_name}" ${file_names})
    target_link_libraries("CTest_${test_name}" PRIVATE fortran-support::fortran-support GTest::gtest_main)
    add_test(NAME "CTest_${test_name}" COMMAND "CTest_${test_name}")
    set_property(TEST "CTest_${test_name}" PROPERTY LABELS C)
    set_target_properties("CTest_${test_name}" PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON)
endmacro()
