// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include "mtime_time_test.h"

#include <stdio.h>
#include <stdlib.h>
#include "mtime_calendar.h"
#include "mtime_datetime.h"

typedef struct _time* Time;

START_TEST(test_create_time_from_valid_strings)
{
	// Normal times
	assertTime("05:42:23.109", 5, 42, 23, 109, "05:42:23.109", "054223.109");
	assertTime("05:42:23.109", 5, 42, 23, 109, "05:42:23.109", "054223.109");
	assertTime("05:42:23", 5, 42, 23, 0, "05:42:23.000", "054223.000");
	assertTime("05:42", 5, 42, 0, 0, "05:42:00.000", "054200.000");
	assertTime("05", 5, 0, 0, 0, "05:00:00.000", "050000.000");
	assertTime("054223.109", 5, 42, 23, 109, "05:42:23.109", "054223.109");
	assertTime("054223", 5, 42, 23, 0, "05:42:23.000", "054223.000");
	assertTime("0542", 5, 42, 0, 0, "05:42:00.000", "054200.000");
	assertTime("17:00:07.9", 17, 0, 7, 900, "17:00:07.900", "170007.900");

	// Boundary times
	assertTime("00:00:00.0", 0, 0, 0, 0, "00:00:00.000", "000000.000");
	assertTime("00:00:00.1", 0, 0, 0, 100, "00:00:00.100", "000000.100");
	assertTime("01:00:00.0", 1, 0, 0, 0, "01:00:00.000", "010000.000");
	assertTime("23:59:59.999", 23, 59, 59, 999, "23:59:59.999", "235959.999");
}
END_TEST

START_TEST(test_invalid_strings)
{
	assertInvalidTime(NULL);
	assertInvalidTime("");
	assertInvalidTime("abc");
	assertInvalidTime("05:42:23.1000");
	assertInvalidTime("05:42:60.109");
	assertInvalidTime("05:60:23.109");
	assertInvalidTime("24:42:23.109");
	assertInvalidTime("24:00:00.000");
	assertInvalidTime("5:42");
	// And many more...
}
END_TEST

START_TEST(test_constructAndCopyTime)
{
	int hour = 13;
	int minute = 4;
	int second = 29;
	int ms = 678;
	
	Time time_src = newRawTime(hour, minute, second, ms);
	Time time_dest = constructAndCopyTime(time_src);
	ck_assert(time_dest != NULL);
	ck_assert_int_eq(hour, time_dest->hour);
	ck_assert_int_eq(minute, time_dest->minute);
	ck_assert_int_eq(second, time_dest->second);
	ck_assert_int_eq(ms, time_dest->ms);
	
	deallocateTime(time_src);
	deallocateTime(time_dest);
}
END_TEST

START_TEST(test_replaceTime)
{
	int hour = 13;
	int minute = 4;
	int second = 29;
	int ms = 678;
	
	Time time_src = newRawTime(hour, minute, second, ms);
	Time time_dest = newTime("00:00:00.0");
	replaceTime(time_src, time_dest);
	ck_assert(time_dest != NULL);
	ck_assert_int_eq(hour, time_dest->hour);
	ck_assert_int_eq(minute, time_dest->minute);
	ck_assert_int_eq(second, time_dest->second);
	ck_assert_int_eq(ms, time_dest->ms);

	deallocateTime(time_src);
	deallocateTime(time_dest);
}
END_TEST

START_TEST(test_timeToPosixString)
{
	char tmp[MAX_DATE_STR_LEN];
	char* format1 = "%X";
	char* format2 = "%l:%M:%S";
	
	Time time = newTime("14:04:59.666");
	ck_assert_str_eq("14:04:59", timeToPosixString(time, tmp, format1));
	ck_assert_str_eq(" 2:04:59", timeToPosixString(time, tmp, format2));
	
	deallocateTime(time);
}
END_TEST

static void setup(void)
{
	initCalendar(PROLEPTIC_GREGORIAN);
}

static void teardown(void)
{
	freeCalendar();
}

void add_mtime_time_test_to_suite(Suite* suite)
{
    TCase *tcase = tcase_create("mtime_time_test");
    suite_add_tcase(suite, tcase);
    tcase_add_checked_fixture(tcase, setup, teardown);
    tcase_add_test(tcase, test_create_time_from_valid_strings);
    tcase_add_test(tcase, test_invalid_strings);
    tcase_add_test(tcase, test_constructAndCopyTime);
    tcase_add_test(tcase, test_replaceTime);
    tcase_add_test(tcase, test_timeToPosixString);
}

/*** SPECIAL ASSERT FUNCTIONS ***/

void assertTime(const char* input_string, int hour, int minute, int second, int ms, const char* expected_output_string, const char* expected_basic_string)
{
	const char* format = "Parsing of time string \"%s\" failed.";
	size_t length = snprintf(NULL, 0, format, input_string) + 1;
	char* msg = malloc(length);
	snprintf(msg, length, format, input_string);

	char tmp[MAX_DATE_STR_LEN];

	Time time = newTime(input_string);
	ck_assert_msg(time != NULL, msg);
	ck_assert_str_eq(expected_output_string, timeToString(time, tmp));
	ck_assert_str_eq(expected_output_string, tmp);
	ck_assert_int_eq(hour, time->hour);
	ck_assert_int_eq(minute, time->minute);
	ck_assert_int_eq(second, time->second);
	ck_assert_int_eq(ms, time->ms);
	ck_assert_str_eq(expected_basic_string, timeToBasicString(time, tmp));
	ck_assert_str_eq(expected_basic_string, tmp);

	Time raw_time = newRawTime(hour, minute, second, ms);
	ck_assert(raw_time != NULL);
	ck_assert_str_eq(expected_output_string, timeToString(raw_time, tmp));
	ck_assert_int_eq(hour, raw_time->hour);
	ck_assert_int_eq(minute, raw_time->minute);
	ck_assert_int_eq(second, raw_time->second);
	ck_assert_int_eq(ms, raw_time->ms);
	ck_assert_str_eq(expected_basic_string, timeToBasicString(raw_time, tmp));

	Time time_from_string = newTime(timeToString(raw_time, tmp));
	ck_assert(time_from_string != NULL);
	ck_assert_str_eq(expected_output_string, timeToString(time_from_string, tmp));

	Time time_from_basic_string = newTime(timeToBasicString(raw_time, tmp));
	ck_assert(time_from_basic_string != NULL);
	ck_assert_str_eq(expected_basic_string, timeToBasicString(time_from_basic_string, tmp));

	deallocateTime(time);
	deallocateTime(raw_time);
	deallocateTime(time_from_string);
	deallocateTime(time_from_basic_string);
	free(msg);
}

void assertInvalidTime(const char* input_string)
{
	Time time = newTime(input_string);

	if (time != NULL)
	{
		const char* format = "Invalid time string \"%s\" didn't result in null pointer but in time \"%s\".";
		char output_string[MAX_DATE_STR_LEN];
		timeToString(time, output_string);

		size_t length = snprintf(NULL, 0, format, input_string, output_string) + 1;
		char* msg = malloc(length);
		snprintf(msg, length, format, input_string, output_string);
		deallocateTime(time);

		ck_abort_msg(msg);
		free(msg); //Never executed?
	}
}
