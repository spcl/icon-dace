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
#include <cmath>

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

    void EXPECT_T_ITEM_VAL_EQ(t_item item, double val) {
        EXPECT_PRED2([](auto item,
                        auto val) { return item.type == 1 && item.val == val; },
                     item, val);
    }

    void EXPECT_T_ITEM_OP_EQ(t_item item, char op) {
        EXPECT_PRED2(
            [](auto item, auto op) { return item.type == 2 && item.op == op; },
            item, op);
    }

    void EXPECT_T_ITEM_FCT_EQ(t_item item, fct_type fct) {
        EXPECT_PRED2([](auto item,
                        auto fct) { return item.type == 3 && item.fct == fct; },
                     item, fct);
    }

    void EXPECT_T_ITEM_VAR_EQ(t_item item, char *field) {
        EXPECT_PRED2(
            [](auto item, auto field) {
                return item.type == 0 && strcpy(item.field, field);
            },
            item, field);
    }

    void EXPECT_FCT_NAME_EQ(char *string, const char *fct_name) {
        EXPECT_PRED2(
            [](auto string, auto fct_name) { return strcpy(string, fct_name); },
            string, fct_name);
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
    std::string expression = "20 * 10 / 20";
    const char *char_array = expression.c_str();

    t_list queue;

    do_parse_infix(char_array, &queue);

    // Reverse Polish notation: 20 10 * 20 /
    EXPECT_T_ITEM_VAL_EQ(queue.list[0], 20);
    EXPECT_T_ITEM_VAL_EQ(queue.list[1], 10);
    EXPECT_T_ITEM_OP_EQ(queue.list[2], '*');
    EXPECT_T_ITEM_VAL_EQ(queue.list[3], 20);
    EXPECT_T_ITEM_OP_EQ(queue.list[4], '/');
}

TEST_F(UtilArithmeticExprTest, ParseInfixIsCorrect2) {
    std::string expression = "3 + 4 * 2 / (1 - 5) ^ 2 ^ 3";
    const char *char_array = expression.c_str();

    t_list queue;

    do_parse_infix(char_array, &queue);

    // Reverse Polish notation: 3 4 2 * 1 5 - 2 ^ 3 ^ / +
    EXPECT_T_ITEM_VAL_EQ(queue.list[0], 3);
    EXPECT_T_ITEM_VAL_EQ(queue.list[1], 4);
    EXPECT_T_ITEM_VAL_EQ(queue.list[2], 2);
    EXPECT_T_ITEM_OP_EQ(queue.list[3], '*');
    EXPECT_T_ITEM_VAL_EQ(queue.list[4], 1);
    EXPECT_T_ITEM_VAL_EQ(queue.list[5], 5);
    EXPECT_T_ITEM_OP_EQ(queue.list[6], '-');
    EXPECT_T_ITEM_VAL_EQ(queue.list[7], 2);
    EXPECT_T_ITEM_OP_EQ(queue.list[8], '^');
    EXPECT_T_ITEM_VAL_EQ(queue.list[9], 3);
    EXPECT_T_ITEM_OP_EQ(queue.list[10], '^');
    EXPECT_T_ITEM_OP_EQ(queue.list[11], '/');
    EXPECT_T_ITEM_OP_EQ(queue.list[12], '+');
}

TEST_F(UtilArithmeticExprTest, ParseInfixIsCorrect3) {
    std::string expression = "sin ( max ( 2, 3 ) / 3 * pi )";
    const char *char_array = expression.c_str();

    t_list queue;

    do_parse_infix(char_array, &queue);

    // Reverse Polish notation: 2 3 max 3 / pi * sin
    EXPECT_T_ITEM_VAL_EQ(queue.list[0], 2);
    EXPECT_T_ITEM_VAL_EQ(queue.list[1], 3);
    EXPECT_T_ITEM_FCT_EQ(queue.list[2], MAX);
    EXPECT_T_ITEM_VAL_EQ(queue.list[3], 3);
    EXPECT_T_ITEM_OP_EQ(queue.list[4], '/');
    EXPECT_T_ITEM_VAL_EQ(queue.list[5], M_PI);
    EXPECT_T_ITEM_OP_EQ(queue.list[6], '*');
    EXPECT_T_ITEM_FCT_EQ(queue.list[7], SIN);
}

TEST_F(UtilArithmeticExprTest, ParseInfixIsCorrect4) {
    std::string expression = "if([z_sfc] > 2, [z_sfc], 0)";
    const char *char_array = expression.c_str();

    t_list queue;

    do_parse_infix(char_array, &queue);

    // Reverse Polish notation: z_sfc 2 > z_sfc 0 if
    EXPECT_T_ITEM_VAR_EQ(queue.list[0], (char *) "z_sfc");
    EXPECT_T_ITEM_VAL_EQ(queue.list[1], 2);
    EXPECT_T_ITEM_OP_EQ(queue.list[2], '>');
    EXPECT_T_ITEM_VAR_EQ(queue.list[3], (char *) "z_sfc");
    EXPECT_T_ITEM_VAL_EQ(queue.list[4], 0);
    EXPECT_T_ITEM_FCT_EQ(queue.list[5], IF);
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
        EXPECT_FCT_NAME_EQ(string, fct_name[i]);
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
