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
