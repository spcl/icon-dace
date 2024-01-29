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

#include <string>
#include <iostream>

#include <util_file.h>

static std::string working_dir = ".";

class UtilFileTest : public ::testing::Test {};


TEST_F(UtilFileTest, FileIsLink) {
    std::string file_notlink = working_dir + "/util_file_test.txt";
    std::string file_islink = working_dir + "/util_file_link.txt";
    std::string file_error = working_dir + "/util_file_noexist.txt";
    char *file_cstr;

    file_cstr = &file_notlink[0];
    EXPECT_EQ(util_islink(file_cstr), 0);
    file_cstr = &file_islink[0];
    EXPECT_EQ(util_islink(file_cstr), 1);
    file_cstr = &file_error[0];
    EXPECT_EQ(util_islink(file_cstr), -1);
}


TEST_F(UtilFileTest, CanGetFileSize) {
    std::string file = working_dir + "/util_file_test.txt";

    char *file_cstr;
    file_cstr = &file[0];

    EXPECT_EQ(util_filesize(file_cstr), 51);
}


TEST_F(UtilFileTest, CheckFileWritable) {
    std::string file = working_dir + "/util_file_test.txt";

    char *file_cstr;
    file_cstr = &file[0];

    EXPECT_EQ(util_file_is_writable(file_cstr), 1);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    if (argc > 1) {
        working_dir = argv[1];
    }

    std::cout << "Working directory: " << working_dir << std::endl;

    return RUN_ALL_TESTS();
}
