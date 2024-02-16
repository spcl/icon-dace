// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#ifndef _MTIME_TIME_TEST_H
#define _MTIME_TIME_TEST_H

#include <check.h>
#include "mtime_time.h"

void add_mtime_time_test_to_suite(Suite* suite);

/*** SPECIAL ASSERT FUNCTIONS ***/
void assertTime(const char* input_string, int hour, int minute, int second, int ms, const char* expected_output_string, const char* expected_basic_string);
void assertInvalidTime(const char* input_string);

#endif
