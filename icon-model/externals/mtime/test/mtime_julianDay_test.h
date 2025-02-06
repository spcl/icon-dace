// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#ifndef _MTIME_JULIANDAY_TEST_H
#define _MTIME_JULIANDAY_TEST_H

#include <check.h>
#include <stdint.h>
#include "mtime_julianDay.h"

void add_mtime_julianDay_test_to_suite(Suite *suite);

/*** SPECIAL ASSERT FUNCTIONS ***/
void assertJulianDay(int64_t day, int64_t ms, const char *expected_string);
void assertJulianDelta(char sign, int64_t day, int64_t ms, char expected_sign, int64_t expected_day, int64_t expected_ms);
void assertAddJulianDelta(int64_t original_day, int64_t original_ms, int64_t delta_day, int64_t delta_ms, int64_t expected_day,
                          int64_t expected_ms);
void assertSubtractJulianDelta(int64_t original_day, int64_t original_ms, int64_t delta_day, int64_t delta_ms, int64_t expected_day,
                               int64_t expected_ms);
void assertSubtractJulianDay(int64_t day1, int64_t ms1, int64_t day2, int64_t ms2, char expected_delta_sign,
                             int64_t expected_delta_day, int64_t expected_delta_ms);

#endif
