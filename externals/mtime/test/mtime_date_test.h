// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#ifndef _MTIME_DATE_TEST_H
#define _MTIME_DATE_TEST_H

#include <check.h>
#include "mtime_date.h"

void add_mtime_date_test_to_suite(Suite *suite);

/*** SPECIAL ASSERT FUNCTIONS ***/
void assertDate(const char *input_string, int64_t year, int month, int day, const char *expected_output_string,
                const char *expected_basic_string);
void assertInvalidDate(const char *input_string);

#endif
