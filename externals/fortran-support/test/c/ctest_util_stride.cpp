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
#include <experimental/random>

#include <util_stride.h>

class UtilStrideTest : public ::testing::Test {};

TEST_F(UtilStrideTest, CanGetPointerDiff) {
    double *a, *b;
    EXPECT_NO_THROW({ util_get_ptrdiff(a, b); });
}

TEST_F(UtilStrideTest, PointerDiffIsCorrect) {
    float f_array[10];
    double d_array[10];

    EXPECT_EQ(util_get_ptrdiff(&f_array[0], &f_array[1]), 4);
    EXPECT_EQ(util_get_ptrdiff(&f_array[2], &f_array[7]), 20);
    EXPECT_EQ(util_get_ptrdiff(&f_array[6], &f_array[3]), 12);
    EXPECT_EQ(util_get_ptrdiff(&d_array[0], &d_array[1]), 8);
    EXPECT_EQ(util_get_ptrdiff(&d_array[0], &d_array[9]), 72);
    EXPECT_EQ(util_get_ptrdiff(&d_array[9], &d_array[0]), 72);
    EXPECT_EQ(util_get_ptrdiff(&d_array[5], &d_array[2]), 24);
}

TEST_F(UtilStrideTest, CanGet1DStride) {
    int stride;

    float f_array[100];
    double d_array[100];

    util_stride_1d(&stride, sizeof(float), &f_array[15], &f_array[25]);
    EXPECT_EQ(stride, 10);
    util_stride_1d(&stride, sizeof(double), &d_array[37], &d_array[81]);
    EXPECT_EQ(stride, 44);
}

TEST_F(UtilStrideTest, CanGet1DStride2) {
    int stride;

    float f_array[1000];
    double d_array[1000];

    int p1 = std::experimental::randint(0, 499);
    int p2 = std::experimental::randint(500, 999);

    util_stride_1d(&stride, sizeof(float), &f_array[p1], &f_array[p2]);
    EXPECT_EQ(stride, p2 - p1);
    util_stride_1d(&stride, sizeof(double), &d_array[p1], &d_array[p2]);
    EXPECT_EQ(stride, p2 - p1);
}

TEST_F(UtilStrideTest, CanGet2DStride) {
    int stride[2];

    float f_array[100];
    double d_array[100];

    util_stride_2d(&stride[0], sizeof(float), &f_array[15], &f_array[25],
                   &f_array[45]);
    EXPECT_EQ(stride[0], 10);
    EXPECT_EQ(stride[1], 30);
    util_stride_2d(&stride[0], sizeof(double), &d_array[37], &d_array[48],
                   &d_array[81]);
    EXPECT_EQ(stride[0], 11);
    EXPECT_EQ(stride[1], 44);
}

TEST_F(UtilStrideTest, CanGet2DStride2) {
    int stride[2];

    float f_array[1000];
    double d_array[1000];

    int p1 = std::experimental::randint(0, 499);
    int p2 = std::experimental::randint(500, 999);
    int p3 = std::experimental::randint(500, 999);

    util_stride_2d(&stride[0], sizeof(float), &f_array[p1], &f_array[p2],
                   &f_array[p3]);
    EXPECT_EQ(stride[0], p2 - p1);
    EXPECT_EQ(stride[1], p3 - p1);
    util_stride_2d(&stride[0], sizeof(double), &d_array[p1], &d_array[p2],
                   &d_array[p3]);
    EXPECT_EQ(stride[0], p2 - p1);
    EXPECT_EQ(stride[1], p3 - p1);
}
