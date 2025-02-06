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

#include <util_backtrace.h>

class UtilBacktraceCTest : public ::testing::Test {
  protected:
    static bool ContainsSubstring(const std::string &str,
                                  const std::string &substr) {
        return str.find(substr) != std::string::npos;
    }
};

TEST_F(UtilBacktraceCTest, CanCallBacktrace) {
    EXPECT_NO_THROW({ util_backtrace(); });
}

TEST_F(UtilBacktraceCTest, PrintsCorrectBacktrace) {
    testing::internal::CaptureStdout();

    util_backtrace();

    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_PRED2(ContainsSubstring, output, "util_backtrace");
    EXPECT_PRED2(ContainsSubstring, output, "libfortran-support");
}
