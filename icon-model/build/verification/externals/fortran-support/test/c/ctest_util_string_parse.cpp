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

#include <vector>
#include <string>

#include <util_string_parse.h>

class UtilStringParseTest : public ::testing::Test {
  public:
    void EXPECT_VECTOR_EQ(std::vector<int> &v1, std::vector<int> &v2) {
        EXPECT_EQ(v1.size(), v2.size());
        for (int i = 0; i < v1.size(); ++i) {
            EXPECT_EQ(v1[i], v2[i])
                << "Vector element mismatch at index i = " << i;
        }
    }
};

TEST_F(UtilStringParseTest, CanParseIntList) {
    std::string parse_line = "1,3,5";
    std::vector<int> output(11);
    int ierr;

    const char *char_array = parse_line.c_str();

    EXPECT_NO_THROW({
        do_parse_intlist(char_array, output.size(), 10, output.data(), &ierr);
    });
}

TEST_F(UtilStringParseTest, ParseIntListIsCorrect) {
    int nlev = 10;

    // comma and semicolon both separates the numbers
    std::string parse_line = "1,2,3;nlev";
    std::vector<int> result = { 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1 };

    // One extra index [0] unused because Fortran index starts from 1
    std::vector<int> output(nlev + 1);
    int ierr;

    const char *char_array = parse_line.c_str();

    do_parse_intlist(char_array, output.size(), nlev, output.data(), &ierr);

    EXPECT_VECTOR_EQ(output, result);
}

TEST_F(UtilStringParseTest, ParseIntListIsCorrect2) {
    int nlev = 10;

    // comma and semicolon both separates the numbers
    std::string parse_line = "1;3,4...7";
    std::vector<int> result = { 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0 };

    // One extra index [0] unused because Fortran index starts from 1
    std::vector<int> output(nlev + 1);
    int ierr;

    const char *char_array = parse_line.c_str();

    do_parse_intlist(char_array, output.size(), nlev, output.data(), &ierr);

    EXPECT_VECTOR_EQ(output, result);
}

TEST_F(UtilStringParseTest, ParseIntListIsCorrect3) {
    int nlev = 30;

    // comma and semicolon both separates the numbers
    std::string parse_line = "1,3,5...10,20...nlev";
    // clang-format off
    std::vector<int> result = { 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    // clang-format on

    // One extra index [0] unused because Fortran index starts from 1
    std::vector<int> output(nlev + 1);
    int ierr;

    const char *char_array = parse_line.c_str();

    do_parse_intlist(char_array, output.size(), nlev, output.data(), &ierr);

    EXPECT_VECTOR_EQ(output, result);
}

TEST_F(UtilStringParseTest, ParseIntListIsCorrect4) {
    int nlev = 30;

    // comma and semicolon both separates the numbers
    std::string parse_line = "1,2, 10 ...22;2;16-(3+11), N-2,16-(2+10);5";
    // clang-format off
    std::vector<int> result = { 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 
                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 0, 0, 0, 0, 0, 1, 0, 0 };
    // clang-format on

    // One extra index [0] unused because Fortran index starts from 1
    std::vector<int> output(nlev + 1);
    int ierr;

    const char *char_array = parse_line.c_str();

    do_parse_intlist(char_array, output.size(), nlev, output.data(), &ierr);

    EXPECT_VECTOR_EQ(output, result);
}
