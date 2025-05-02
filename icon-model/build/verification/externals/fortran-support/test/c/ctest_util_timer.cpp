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

#include <util_timer.h>

class UtilTimerTest : public ::testing::Test {};

TEST_F(UtilTimerTest, CanGetUserTime) {
    double user_time, system_time;

    EXPECT_NO_THROW({ util_cputime(&user_time, &system_time); });
}

TEST_F(UtilTimerTest, CanGetWallTime) {
    EXPECT_NO_THROW({ double walltime = util_walltime(); });
}

TEST_F(UtilTimerTest, CanGetTimeOfDay) {
    EXPECT_NO_THROW({ double walltime = util_gettimeofday(); });
}

TEST_F(UtilTimerTest, CanInitializeFineTimer) {
    EXPECT_NO_THROW({ util_init_real_time(); });
}

TEST_F(UtilTimerTest, RealTimeSizeIsCorrect) {
    int rt_size;

    util_get_real_time_size(&rt_size);

    EXPECT_EQ(rt_size, 8);
}

TEST_F(UtilTimerTest, CanReadRealTime) {
    void *time_ptr = ::operator new(1);

    EXPECT_NO_THROW({ util_read_real_time(time_ptr); });
}

TEST_F(UtilTimerTest, CanGetDiffRealTime) {
    double t1 = 999, t2 = 1000, t_diff;
    void *ptr1 = &t1, *ptr2 = &t2;

    util_diff_real_time(ptr1, ptr2, &t_diff);

    EXPECT_EQ(t_diff, 1);
}
