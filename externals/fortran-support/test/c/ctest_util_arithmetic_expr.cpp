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

#include <util_arithmetic_expr.h>

class UtilArithmeticExprTest : public ::testing::Test {
  protected:
    const char *op_parentheses = "(";
    const char *op_greater = ">";
    const char *op_lesser = "<";
    const char *op_addition = "+";
    const char *op_subtraction = "-";
    const char *op_division = "/";
    const char *op_multiplication = "*";
    const char *op_power = "^";

    void EXPECT_FCT_EQ(char *string, const char *fct) {
        EXPECT_PRED2([](auto string, auto fct) { return strcpy(string, fct); },
                     string, fct);
    }
};

TEST_F(UtilArithmeticExprTest, CanRunPriority) {
    EXPECT_NO_THROW(int x = priority(*op_addition));
}

TEST_F(UtilArithmeticExprTest, PriorityIsCorrect) {
    EXPECT_EQ(priority(*op_parentheses), 0);
    EXPECT_EQ(priority(*op_greater), 1);
    EXPECT_EQ(priority(*op_lesser), 1);
    EXPECT_EQ(priority(*op_addition), 2);
    EXPECT_EQ(priority(*op_subtraction), 2);
    EXPECT_EQ(priority(*op_division), 3);
    EXPECT_EQ(priority(*op_multiplication), 3);
    EXPECT_EQ(priority(*op_power), 4);
}

TEST_F(UtilArithmeticExprTest, CanRunLeftAssociative) {
    EXPECT_NO_THROW(int l = left_associative(*op_addition));
}

TEST_F(UtilArithmeticExprTest, LeftAssociativeIsCorrect) {
    EXPECT_EQ(left_associative(*op_parentheses), 0);
    EXPECT_EQ(left_associative(*op_greater), 1);
    EXPECT_EQ(left_associative(*op_lesser), 1);
    EXPECT_EQ(left_associative(*op_addition), 1);
    EXPECT_EQ(left_associative(*op_subtraction), 1);
    EXPECT_EQ(left_associative(*op_division), 1);
    EXPECT_EQ(left_associative(*op_multiplication), 1);
    EXPECT_EQ(left_associative(*op_power), 1);
}

TEST_F(UtilArithmeticExprTest, CanDoParseInfix) {
    std::string expression = "1+2*3/4-5*6+7/8";
    const char *char_array = expression.c_str();

    t_list queue;

    do_parse_infix(char_array, &queue);
}

TEST_F(UtilArithmeticExprTest, ParseInfixIsCorrect) {
    std::string expression = "20*10/20";
    const char *char_array = expression.c_str();

    t_list queue;

    do_parse_infix(char_array, &queue);

    // TODO: add more concrete examples for do_parse_infix
    EXPECT_EQ(queue.list[2].op, '*');
    EXPECT_EQ(queue.list[4].op, '/');

    EXPECT_EQ(queue.list[0].val, 20);
    EXPECT_EQ(queue.list[1].val, 10);
    EXPECT_EQ(queue.list[3].val, 20);
}

TEST_F(UtilArithmeticExprTest, CanGetFCTName) {
    char string[4];
    int ierr;

    ierr = get_fctname(0, string);
    ierr = get_fctname(1, string);
    ierr = get_fctname(8, string);
}

TEST_F(UtilArithmeticExprTest, GetFCTNameIsCorrect) {
    char string[5];
    int ierr;

    for (int i = 0; i < NUM_FCT; ++i) {
        ierr = get_fctname(i, string);
        EXPECT_EQ(ierr, 0);
        EXPECT_FCT_EQ(string, fct_name[i]);
    }
}

TEST_F(UtilArithmeticExprTest, GetFCTNameThrowsError) {
    char string[5];
    int ierr;

    ierr = get_fctname(-1, string);
    EXPECT_EQ(ierr, 1);

    ierr = get_fctname(9, string);
    EXPECT_EQ(ierr, 1);
}
