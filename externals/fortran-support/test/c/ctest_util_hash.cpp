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

#include <stdexcept>
#include <string>

#include <util_hash.h>

class UtilHashCTest : public ::testing::Test {};

TEST_F(UtilHashCTest, CanCallHashword) {
    std::string s = "Unittest";

    EXPECT_NO_THROW({ util_hashword(&s[0], 0, 0); });
}

TEST_F(UtilHashCTest, HashwordIsCorrect1) {
    std::string s = "Unittest";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 345061529);
}

TEST_F(UtilHashCTest, HashwordIsCorrect2) {
    std::string s = "UnittestFrameworkC";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 1976263765);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l0) {
    std::string s = "";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 3735928559);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l1) {
    std::string s = "U";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 2658498343);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l2) {
    std::string s = "UT";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 1367754127);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l3) {
    std::string s = "UiT";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 1170274080);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l4) {
    std::string s = "UniT";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 4137270609);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l5) {
    std::string s = "UniTt";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 2664783129);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l6) {
    std::string s = "UniTet";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 604414756);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l7) {
    std::string s = "UniTest";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 3603701657);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l8) {
    std::string s = "UnitTest";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 2518181024);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l9) {
    std::string s = "UnitTests";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 1752525959);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l10) {
    std::string s = "CUnitTests";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 693311643);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l11) {
    std::string s = "UnitTestTMP";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 1109590718);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l12) {
    std::string s = "UnitTestForC";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 3182545988);
}

TEST_F(UtilHashCTest, HashwordIsCorrect_l13) {
    std::string s = "UnitTestFrCPP";

    uint32_t hashword = util_hashword(&s[0], s.length(), 0);

    EXPECT_EQ(hashword, 3082204807);
}
