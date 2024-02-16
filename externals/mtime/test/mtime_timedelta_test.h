// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#ifndef _MTIME_TIMEDELTA_TEST_H
#define _MTIME_TIMEDELTA_TEST_H

#include <check.h>
#include <stdint.h>
#include "mtime_timedelta.h"

void add_mtime_timedelta_test_to_suite(Suite* suite);

/*** SPECIAL ASSERT FUNCTIONS ***/
void assertTimeDelta(const char* input_string, bool std_form, char sign, int64_t year, int month, int day, int hour, int minute, int second, int ms, const char* expected_output_string);
void assertInvalidTimeDelta(const char* input_string);
void assertTimeDeltaToJulianDelta(const char* timedelta_string, const char* base_datetime_string, char expected_jd_sign, int64_t expected_jd_day, int64_t expected_jd_ms);
void assertTimeDeltaToJulianDeltaWithoutExpectations(const char* timedelta_string, const char* base_datetime_string, char expected_jd_sign);

void assertGetTimeDeltaFromDateTime (const char* dt1_string, const char* dt2_string, const char* expected_td_string);
void assertTimeDeltaToJulianDeltaToTimeDelta (const char* base_dt_string, const char* td_string);
#endif
