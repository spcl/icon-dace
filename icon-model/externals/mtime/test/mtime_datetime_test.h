// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#ifndef _MTIME_DATETIME_TEST_H
#define _MTIME_DATETIME_TEST_H

#include <check.h>
#include <stdint.h>
#include "mtime_datetime.h"
#include "mtime_date.h"

void add_mtime_datetime_test_to_suite(Suite *suite);

/*** SPECIAL ASSERT FUNCTIONS ***/
void assertDateTime(const char *input_string, int64_t year, int month, int day, int hour, int minute, int second, int ms,
                    const char *expected_output_string, const char *expected_basic_string);
void assertInvalidDateTime(const char *input_string);
void assertNoDaysInMonth(const char *date_string, int expected_days);
void assertNoDaysInYear(const char *date_string, int expected_days);
void assertDayOfYear(const char *date_string, int expected_day);
void assertNoOfSecondsInMonth(const char *date_string, int64_t expected_seconds);
void assertNoOfSecondsInDay(const char *date_string, int expected_days);
void assertJulianDayFromDatetime(const char *date_string, int expected_day, int expected_ms);

#endif
