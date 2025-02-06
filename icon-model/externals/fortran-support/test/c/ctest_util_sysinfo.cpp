// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

#include <gtest/gtest.h>

#include <util_sysinfo.h>

class UtilSysinfoTest : public ::testing::Test {};

TEST_F(UtilSysinfoTest, CanGetUserName) {
    char name[100];
    int actual_len;

    EXPECT_NO_THROW({ util_user_name(name, &actual_len); });
    EXPECT_GT(actual_len, 0);
}

TEST_F(UtilSysinfoTest, CanGetOsSystem) {
    char name[100];
    int actual_len;

    EXPECT_NO_THROW({ util_os_system(name, &actual_len); });
    EXPECT_GT(actual_len, 0);
}

TEST_F(UtilSysinfoTest, CanGetNodeName) {
    char name[100];
    int actual_len;

    EXPECT_NO_THROW({ util_node_name(name, &actual_len); });
    EXPECT_GT(actual_len, 0);
}

TEST_F(UtilSysinfoTest, CanGetMaxResidentSetSize) {
    int maxrss;

    EXPECT_NO_THROW({ util_get_maxrss(&maxrss); });
}

TEST_F(UtilSysinfoTest, CanGetCompilerRelease) {
    char release_str[100];
    int rstr_len;

    EXPECT_NO_THROW({ util_compiler_release(release_str, &rstr_len); });
    EXPECT_GT(rstr_len, 0);
}

TEST_F(UtilSysinfoTest, CanGetProcessId) {
    long int pid;

    EXPECT_NO_THROW({ util_c_getpid(&pid); });
}
