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
#include <filesystem>

#include <util_file.h>

static std::string working_dir;

class UtilFileTest : public ::testing::Test {
  protected:
    static bool ContainsSubstring(const std::string &str,
                                  const std::string &substr) {
        return str.find(substr) != std::string::npos;
    }
};


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

TEST_F(UtilFileTest, CanCreateTmpFile) {
    int filename_len;
    char filename_cstr[30];

    filename_len = util_create_tmpfile(filename_cstr, 30);

    std::cout << "Temporary file " << filename_cstr << " created." << std::endl;

    std::string filename(filename_cstr);

    EXPECT_PRED2(ContainsSubstring, filename, "/tmp/icon");
    EXPECT_EQ(filename_len, 24);

    EXPECT_EQ(util_create_tmpfile(filename_cstr, 10), -1);

    std::filesystem::remove(filename);
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

TEST_F(UtilFileTest, CanCreateSymlink) {
    std::string target_file = working_dir + "/util_file_test.txt";
    std::string new_link = working_dir + "/util_file_new_link.txt";

    EXPECT_EQ(createSymlink(target_file.c_str(), new_link.c_str()), 0);

    char *new_link_cstr = &new_link[0];
    EXPECT_EQ(util_islink(new_link_cstr), 1);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::filesystem::path path = argv[0];

    working_dir = path.parent_path();

    std::cout << "Working directory: " << working_dir << std::endl;

    return RUN_ALL_TESTS();
}
