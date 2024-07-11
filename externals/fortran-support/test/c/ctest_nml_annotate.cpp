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

#include <fstream>
#include <filesystem>

#include <nml_annotate.h>

static std::string working_dir;

class NmlAnnotateTest : public ::testing::Test {
  protected:
    void EXPECT_FILE_EQ(const std::string &file1, const std::string &file2) {
        std::fstream fs1, fs2;
        fs1.open(file1);
        fs2.open(file2);

        std::string line_in_file1, line_in_file2;
        bool error = false;
        int line_num = 0;
        while (fs1) {
            std::getline(fs1, line_in_file1);
            std::getline(fs2, line_in_file2);
            EXPECT_EQ(line_in_file1, line_in_file2)
                << "line " << line_num << " is different";
            line_num++;
        }

        fs1.close();
        fs2.close();
    }
};

TEST_F(NmlAnnotateTest, CanAnnotateNamelist) {
    std::string namelist_file = working_dir + "/test.namelist";
    std::string annotate_file = working_dir + "/annotated.namelist";
    std::string result_file = working_dir + "/result.namelist";

    char *namelist_file_cstr, *annotate_file_cstr;
    namelist_file_cstr = &namelist_file[0];
    annotate_file_cstr = &annotate_file[0];
    EXPECT_NO_THROW(
        { util_annotate_nml(namelist_file_cstr, annotate_file_cstr); });

    EXPECT_FILE_EQ(annotate_file, result_file);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::filesystem::path path = argv[0];

    working_dir = path.parent_path();

    std::cout << "Working directory: " << working_dir << std::endl;

    return RUN_ALL_TESTS();
}
