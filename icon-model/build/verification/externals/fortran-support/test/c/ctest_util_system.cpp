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

#include <util_system.h>

class UtilSystemTest : public ::testing::Test {};

TEST_F(UtilSystemTest, CanAbort) {
    EXPECT_EXIT({ util_abort(); }, testing::ExitedWithCode(1), "");
}

TEST_F(UtilSystemTest, CanExit) {
    EXPECT_EXIT({ util_exit(1); }, testing::ExitedWithCode(1), "");
}

TEST_F(UtilSystemTest, CanExit2) {
    EXPECT_EXIT({ util_exit(2); }, testing::ExitedWithCode(2), "");
}

TEST_F(UtilSystemTest, CanCallSystemDate) {
    char sys[] = "date";

    EXPECT_NO_THROW({ util_system(sys); });
}
