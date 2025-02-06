// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include "mtime_julianDay_test.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "mtime_calendar.h"

typedef struct _julianday *JulianDay;
typedef struct _juliandelta *JulianDelta;

START_TEST(test_create_julianDay_ProlepticGregorian)
{
  //	assertJulianDay(-784350575246, 0, "-784350575246.0"); //FAILS
  //	assertJulianDay(-784350575246, 43200000 - 1, "-784350575246.43199999"); //FAILS
  assertJulianDay(-784350575246, 43200000, "-784350575246.43200000");
  assertJulianDay(-784350575246, 43200000 + 1, "-784350575246.43200001");
  assertJulianDay(-784350575246, 86400000, "-784350575246.86400000");
  assertJulianDay(-42, 42, "-42.42");
  assertJulianDay(0, 0, "0.0");
  assertJulianDay(0, 86400000, "0.86400000");
  assertJulianDay(89653545, 786, "89653545.786");
  assertJulianDay(784354017364, 0, "784354017364.0");
  assertJulianDay(784354017364, 43200000 - 1, "784354017364.43199999");
  //	assertJulianDay(784354017364, 43200000, "784354017364.43200000"); //FAILS
  //	assertJulianDay(784354017364, 43200000 + 1, "784354017364.43200001"); //FAILS
  //	assertJulianDay(784354017364, 86400000, "784354017364.86400000"); //FAILS
}
END_TEST

START_TEST(test_create_julianDay_YearOf365Days)
{
  //	assertJulianDay(-783831531521, 0, "-783831531521.0"); //FAILS
  //	assertJulianDay(-783831531521, 43200000 - 1, "-783831531521.43199999"); //FAILS
  assertJulianDay(-783831531521, 43200000, "-783831531521.43200000");
  assertJulianDay(-783831531521, 43200000 + 1, "-783831531521.43200001");
  assertJulianDay(-783831531521, 86400000, "-783831531521.86400000");
  assertJulianDay(-42, 42, "-42.42");
  assertJulianDay(0, 0, "0.0");
  assertJulianDay(0, 86400000, "0.86400000");
  assertJulianDay(89653545, 786, "89653545.786");
  assertJulianDay(783831531519, 0, "783831531519.0");
  assertJulianDay(783831531519, 43200000 - 1, "783831531519.43199999");
  //	assertJulianDay(783831531519, 43200000, "783831531519.43200000"); //FAILS
  //	assertJulianDay(783831531519, 43200000 + 1, "783831531519.43200001"); //FAILS
  //	assertJulianDay(783831531519, 86400000, "783831531519.86400000"); //FAILS
}
END_TEST

START_TEST(test_create_julianDay_YearOf360Days)
{
  //	assertJulianDay(-773094113281, 0, "-773094113281.0"); //FAILS
  //	assertJulianDay(-773094113281, 43200000 - 1, "-773094113281.43199999"); //FAILS
  assertJulianDay(-773094113281, 43200000, "-773094113281.43200000");
  assertJulianDay(-773094113281, 43200000 + 1, "-773094113281.43200001");
  assertJulianDay(-773094113281, 86400000, "-773094113281.86400000");
  assertJulianDay(-42, 42, "-42.42");
  assertJulianDay(0, 0, "0.0");
  assertJulianDay(0, 86400000, "0.86400000");
  assertJulianDay(89653545, 786, "89653545.786");
  assertJulianDay(773094113279, 0, "773094113279.0");
  assertJulianDay(773094113279, 43200000 - 1, "773094113279.43199999");
  //	assertJulianDay(773094113279, 43200000, "773094113279.43200000"); //FAILS
  //	assertJulianDay(773094113279, 43200000 + 1, "773094113279.43200001"); //FAILS
  //	assertJulianDay(773094113279, 86400000, "773094113279.86400000"); //FAILS
}
END_TEST

START_TEST(test_create_julianDelta_ProlepticGregorian)
{
  assertJulianDelta('-', -1568704592610, -86400000, '-', -1568704592611, 0);
  assertJulianDelta('-', -1568704592610, -43200000 + 1, '-', -1568704592610, -43200000 + 1);
  assertJulianDelta('-', -1568704592610, -43200000, '-', -1568704592610, -43200000);
  assertJulianDelta('-', -1568704592610, -43200000 - 1, '-', -1568704592610, -43200000 - 1);
  assertJulianDelta('-', -1568704592610, 0, '-', -1568704592610, 0);
  assertJulianDelta('-', -42, -42, '-', -42, -42);
  assertJulianDelta('-', 0, -86400000, '-', -1, 0);
  assertJulianDelta('-', 0, 0, '-', 0, 0);
  assertJulianDelta('+', 0, 0, '+', 0, 0);
  assertJulianDelta('+', 0, 86400000, '+', 1, 0);
  assertJulianDelta('+', 89653545, 786, '+', 89653545, 786);
  assertJulianDelta('+', 1568704592610, 0, '+', 1568704592610, 0);
  assertJulianDelta('+', 1568704592610, 43200000 - 1, '+', 1568704592610, 43200000 - 1);
  assertJulianDelta('+', 1568704592610, 43200000, '+', 1568704592610, 43200000);
  assertJulianDelta('+', 1568704592610, 43200000 + 1, '+', 1568704592610, 43200000 + 1);
  assertJulianDelta('+', 1568704592610, 86400000, '+', 1568704592611, 0);
}
END_TEST

START_TEST(test_create_julianDelta_YearOf365Days)
{
  assertJulianDelta('-', -1567663063040, -86400000, '-', -1567663063041, 0);
  assertJulianDelta('-', -1567663063040, -43200000 + 1, '-', -1567663063040, -43200000 + 1);
  assertJulianDelta('-', -1567663063040, -43200000, '-', -1567663063040, -43200000);
  assertJulianDelta('-', -1567663063040, -43200000 - 1, '-', -1567663063040, -43200000 - 1);
  assertJulianDelta('-', -1567663063040, 0, '-', -1567663063040, 0);
  assertJulianDelta('-', -42, -42, '-', -42, -42);
  assertJulianDelta('-', 0, -86400000, '-', -1, 0);
  assertJulianDelta('-', 0, 0, '-', 0, 0);
  assertJulianDelta('+', 0, 0, '+', 0, 0);
  assertJulianDelta('+', 0, 86400000, '+', 1, 0);
  assertJulianDelta('+', 89653545, 786, '+', 89653545, 786);
  assertJulianDelta('+', 1567663063040, 0, '+', 1567663063040, 0);
  assertJulianDelta('+', 1567663063040, 43200000 - 1, '+', 1567663063040, 43200000 - 1);
  assertJulianDelta('+', 1567663063040, 43200000, '+', 1567663063040, 43200000);
  assertJulianDelta('+', 1567663063040, 43200000 + 1, '+', 1567663063040, 43200000 + 1);
  assertJulianDelta('+', 1567663063040, 86400000, '+', 1567663063041, 0);
}
END_TEST

START_TEST(test_create_julianDelta_YearOf360Days)
{
  assertJulianDelta('-', -1546188226560, -86400000, '-', -1546188226561, 0);
  assertJulianDelta('-', -1546188226560, -43200000 + 1, '-', -1546188226560, -43200000 + 1);
  assertJulianDelta('-', -1546188226560, -43200000, '-', -1546188226560, -43200000);
  assertJulianDelta('-', -1546188226560, -43200000 - 1, '-', -1546188226560, -43200000 - 1);
  assertJulianDelta('-', -1546188226560, 0, '-', -1546188226560, 0);
  assertJulianDelta('-', -42, -42, '-', -42, -42);
  assertJulianDelta('-', 0, -86400000, '-', -1, 0);
  assertJulianDelta('-', 0, 0, '-', 0, 0);
  assertJulianDelta('+', 0, 0, '+', 0, 0);
  assertJulianDelta('+', 0, 86400000, '+', 1, 0);
  assertJulianDelta('+', 89653545, 786, '+', 89653545, 786);
  assertJulianDelta('+', 1546188226560, 0, '+', 1546188226560, 0);
  assertJulianDelta('+', 1546188226560, 43200000 - 1, '+', 1546188226560, 43200000 - 1);
  assertJulianDelta('+', 1546188226560, 43200000, '+', 1546188226560, 43200000);
  assertJulianDelta('+', 1546188226560, 43200000 + 1, '+', 1546188226560, 43200000 + 1);
  assertJulianDelta('+', 1546188226560, 86400000, '+', 1546188226561, 0);
}
END_TEST

START_TEST(test_addJulianDelta_ProlepticGregorian)
{
  assertAddJulianDelta(0, 40000000, 0, 40000000, 0, 80000000);
  assertAddJulianDelta(0, 0, 0, 86400000 - 1, 0, 86400000 - 1);
  assertAddJulianDelta(0, 0, 0, 86400000, 1, 0);
  assertAddJulianDelta(0, 60000000, 0, 60000000, 1, 33600000);
  assertAddJulianDelta(0, 60000000, 42, 60000000, 43, 33600000);
  assertAddJulianDelta(-789643, 60000000, 42, 60000000, -789600, 33600000);
  assertAddJulianDelta(789557, 60000000, 42, 60000000, 789600, 33600000);
  assertAddJulianDelta(0, 0, 773094113279, 43200000 - 1, 773094113279, 43200000 - 1);
  assertAddJulianDelta(386547056640, 43200000 - 1, 386547056638, 86400000, 773094113279, 43200000 - 1);
  assertAddJulianDelta(-773094113281, 43200000, 773094113281, 43200000 - 1, 0, 86400000 - 1);
  assertAddJulianDelta(-773094113281, 43200000, 773094113280, 43200000, 0, 0);
  assertAddJulianDelta(-773094113281, 43200000, 1557448130644, 86400000 - 1, 784354017364, 43200000 - 1);
}
END_TEST

START_TEST(test_subtractJulianDelta_ProlepticGregorian)
{
  assertSubtractJulianDelta(0, 40000000, 0, -40000000, 0, 0);
  assertSubtractJulianDelta(0, 0, 0, -86400000 + 1, -1, 1);
  assertSubtractJulianDelta(0, 0, 0, -86400000, -1, 0);
  assertSubtractJulianDelta(0, 20000000, 0, -60000000, -1, 46400000);
  assertSubtractJulianDelta(0, 20000000, -42, -60000000, -43, 46400000);
  assertSubtractJulianDelta(789643, 20000000, -42, -60000000, 789600, 46400000);
  assertSubtractJulianDelta(-789557, 20000000, -42, -60000000, -789600, 46400000);
  assertSubtractJulianDelta(0, 0, -773094113281, 43200000, -773094113281, 43200000);
  assertSubtractJulianDelta(-386547056640, 43200000, -386547056640, -86400000, -773094113281, 43200000);
  assertSubtractJulianDelta(784354017364, 43200000 - 1, -784354017364, -43200000, -1, 86400000 - 1);
  assertSubtractJulianDelta(784354017364, 43200000 - 1, -784354017364, -43200000 + 1, 0, 0);
  assertSubtractJulianDelta(784354017364, 43200000 - 1, -1557448130644, -86400000 + 1, -773094113281, 43200000);
}
END_TEST

START_TEST(test_subtractJulianDay_ProlepticGregorian)
{
  assertSubtractJulianDay(0, 40000000, 0, 0, '+', 0, 40000000);
  assertSubtractJulianDay(0, 0, -1, 1, '+', 0, 86400000 - 1);
  assertSubtractJulianDay(0, 0, -1, 0, '+', 1, 0);
  assertSubtractJulianDay(0, 20000000, -1, 46400000, '+', 0, 60000000);
  assertSubtractJulianDay(0, 20000000, -43, 46400000, '+', 42, 60000000);
  assertSubtractJulianDay(789643, 20000000, 789600, 46400000, '+', 42, 60000000);
  assertSubtractJulianDay(-789557, 20000000, -789600, 46400000, '+', 42, 60000000);
  assertSubtractJulianDay(0, 0, -773094113281, 43200000, '+', 773094113280, 43200000);
  assertSubtractJulianDay(-386547056640, 43200000, -773094113281, 43200000, '+', 386547056641, 0);
  assertSubtractJulianDay(784354017364, 43200000 - 1, -1, 86400000 - 1, '+', 784354017364, 43200000);
  assertSubtractJulianDay(784354017364, 43200000 - 1, 0, 0, '+', 784354017364, 43200000 - 1);
  assertSubtractJulianDay(784354017364, 43200000 - 1, -773094113281, 43200000, '+', 1557448130644, 86400000 - 1);

  assertSubtractJulianDay(0, 40000000, 0, 80000000, '-', 0, -40000000);
  assertSubtractJulianDay(0, 0, 0, 86400000 - 1, '-', 0, -86400000 + 1);
  assertSubtractJulianDay(0, 0, 1, 0, '-', -1, 0);
  assertSubtractJulianDay(0, 60000000, 1, 33600000, '-', 0, -60000000);
  assertSubtractJulianDay(0, 60000000, 43, 33600000, '-', -42, -60000000);
  assertSubtractJulianDay(-789643, 60000000, -789600, 33600000, '-', -42, -60000000);
  assertSubtractJulianDay(789557, 60000000, 789600, 33600000, '-', -42, -60000000);
  assertSubtractJulianDay(0, 0, 773094113279, 43200000 - 1, '-', -773094113279, -43200000 + 1);
  assertSubtractJulianDay(386547056640, 43200000 - 1, 773094113279, 43200000 - 1, '-', -386547056639, 0);
  assertSubtractJulianDay(-773094113281, 43200000, 0, 86400000 - 1, '-', -773094113281, -43200000 + 1);
  assertSubtractJulianDay(-773094113281, 43200000, 0, 0, '-', -773094113280, -43200000);
  assertSubtractJulianDay(-773094113281, 43200000, 784354017364, 43200000 - 1, '-', -1557448130644, -86400000 + 1);
}
END_TEST

static void
setup_ProlepticGregorian(void)
{
  initCalendar(PROLEPTIC_GREGORIAN);
}

static void
setup_YearOf365Days(void)
{
  initCalendar(YEAR_OF_365_DAYS);
}

static void
setup_YearOf360Days(void)
{
  initCalendar(YEAR_OF_360_DAYS);
}

static void
teardown(void)
{
  freeCalendar();
}

void
add_mtime_julianDay_test_to_suite(Suite *suite)
{
  TCase *tcase_ProlepticGregorian = tcase_create("mtime_julianDay_test_ProlepticGregorian");
  suite_add_tcase(suite, tcase_ProlepticGregorian);
  tcase_add_checked_fixture(tcase_ProlepticGregorian, setup_ProlepticGregorian, teardown);
  tcase_add_test(tcase_ProlepticGregorian, test_create_julianDay_ProlepticGregorian);
  tcase_add_test(tcase_ProlepticGregorian, test_create_julianDelta_ProlepticGregorian);
  tcase_add_test(tcase_ProlepticGregorian, test_addJulianDelta_ProlepticGregorian);
  tcase_add_test(tcase_ProlepticGregorian, test_subtractJulianDelta_ProlepticGregorian);
  tcase_add_test(tcase_ProlepticGregorian, test_subtractJulianDay_ProlepticGregorian);

  TCase *tcase_YearOf365Days = tcase_create("mtime_julianDay_test_YearOf365Days");
  suite_add_tcase(suite, tcase_YearOf365Days);
  tcase_add_checked_fixture(tcase_YearOf365Days, setup_YearOf365Days, teardown);
  tcase_add_test(tcase_YearOf365Days, test_create_julianDay_YearOf365Days);
  tcase_add_test(tcase_YearOf365Days, test_create_julianDelta_YearOf365Days);

  TCase *tcase_YearOf360Days = tcase_create("mtime_julianDay_test_YearOf360Days");
  suite_add_tcase(suite, tcase_YearOf360Days);
  tcase_add_checked_fixture(tcase_YearOf360Days, setup_YearOf360Days, teardown);
  tcase_add_test(tcase_YearOf360Days, test_create_julianDay_YearOf360Days);
  tcase_add_test(tcase_YearOf360Days, test_create_julianDelta_YearOf360Days);
}

/*** SPECIAL ASSERT FUNCTIONS ***/

void
assertJulianDay(int64_t day, int64_t ms, const char *expected_string)
{
  const char *format1 = "Creation of julianday failed for (%lld, %lld).";
  size_t length1 = snprintf(NULL, 0, format1, day, ms) + 1;
  char *msg1 = malloc(length1);
  snprintf(msg1, length1, format1, day, ms);

  JulianDay julianday = newJulianDay(day, ms);
  ck_assert_msg(julianday != NULL, msg1);

  const char *format2 = "newJulianDay failed: expected:(%lld, %lld) actual:(%lld, %lld).";
  size_t length2 = snprintf(NULL, 0, format2, day, ms, julianday->day, julianday->ms) + 1;
  char *msg2 = malloc(length2);
  snprintf(msg2, length2, format2, day, ms, julianday->day, julianday->ms);

  ck_assert_msg(day == julianday->day, msg2);
  ck_assert_msg(ms == julianday->ms, msg2);

  char output_string[MAX_JULIANDAY_STR_LEN];
  char *output_string_ptr = juliandayToString(julianday, output_string);

  ck_assert_str_eq(expected_string, output_string);
  ck_assert_str_eq(expected_string, output_string_ptr);

  deallocateJulianDay(julianday);
  free(msg1);
  free(msg2);
}

void
assertJulianDelta(char sign, int64_t day, int64_t ms, char expected_sign, int64_t expected_day, int64_t expected_ms)
{
  const char *format1 = "Creation of juliandelta failed for (%c, %lld, %lld).";
  size_t length1 = snprintf(NULL, 0, format1, sign, day, ms) + 1;
  char *msg1 = malloc(length1);
  snprintf(msg1, length1, format1, sign, day, ms);

  JulianDelta juliandelta = newJulianDelta(sign, day, ms);
  ck_assert_msg(juliandelta != NULL, msg1);

  const char *format2 = "newJulianDelta failed: expected:(%c, %lld, %lld) actual:(%c, %lld, %lld).";
  size_t length2
      = snprintf(NULL, 0, format2, expected_sign, expected_day, expected_ms, juliandelta->sign, juliandelta->day, juliandelta->ms)
        + 1;
  char *msg2 = malloc(length2);
  snprintf(msg2, length2, format2, expected_sign, expected_day, expected_ms, juliandelta->sign, juliandelta->day, juliandelta->ms);

  ck_assert_msg(expected_sign == juliandelta->sign, msg2);
  ck_assert_msg(expected_day == juliandelta->day, msg2);
  ck_assert_msg(expected_ms == juliandelta->ms, msg2);

  deallocateJulianDelta(juliandelta);
  free(msg1);
  free(msg2);
}

void
assertAddJulianDelta(int64_t original_day, int64_t original_ms, int64_t delta_day, int64_t delta_ms, int64_t expected_day,
                     int64_t expected_ms)
{
  JulianDay originalday = newJulianDay(original_day, original_ms);
  ck_assert(originalday != NULL);

  JulianDelta delta = newJulianDelta('+', delta_day, delta_ms);
  ck_assert(delta != NULL);

  JulianDay resultday = addJulianDeltaToJulianDay(originalday, delta, newJulianDay(0, 0));
  ck_assert(delta != NULL);

  const char *format = "addJulianDeltaToJulianDay failed: expected:(%lld, %lld) actual:(%lld, %lld).";
  size_t length = snprintf(NULL, 0, format, expected_day, expected_ms, resultday->day, resultday->ms) + 1;
  char *msg = malloc(length);
  snprintf(msg, length, format, expected_day, expected_ms, resultday->day, resultday->ms);

  ck_assert_msg(expected_day == resultday->day, msg);
  ck_assert_msg(expected_ms == resultday->ms, msg);

  deallocateJulianDay(originalday);
  deallocateJulianDelta(delta);
  deallocateJulianDay(resultday);
  free(msg);
}

void
assertSubtractJulianDelta(int64_t original_day, int64_t original_ms, int64_t delta_day, int64_t delta_ms, int64_t expected_day,
                          int64_t expected_ms)
{
  JulianDay originalday = newJulianDay(original_day, original_ms);
  ck_assert(originalday != NULL);

  JulianDelta delta = newJulianDelta('-', delta_day, delta_ms);
  ck_assert(delta != NULL);

  JulianDay resultday = subtractJulianDeltaFromJulianDay(originalday, delta, newJulianDay(0, 0));
  ck_assert(delta != NULL);

  const char *format = "subtractJulianDelta failed: expected:(%lld, %lld) actual:(%lld, %lld).";
  size_t length = snprintf(NULL, 0, format, expected_day, expected_ms, resultday->day, resultday->ms) + 1;
  char *msg = malloc(length);
  snprintf(msg, length, format, expected_day, expected_ms, resultday->day, resultday->ms);

  ck_assert_msg(expected_day == resultday->day, msg);
  ck_assert_msg(expected_ms == resultday->ms, msg);

  deallocateJulianDay(originalday);
  deallocateJulianDelta(delta);
  deallocateJulianDay(resultday);
  free(msg);
}

void
assertSubtractJulianDay(int64_t day1, int64_t ms1, int64_t day2, int64_t ms2, char expected_delta_sign, int64_t expected_delta_day,
                        int64_t expected_delta_ms)
{
  JulianDay julianday1 = newJulianDay(day1, ms1);
  ck_assert(julianday1 != NULL);

  JulianDay julianday2 = newJulianDay(day2, ms2);
  ck_assert(julianday2 != NULL);

  JulianDelta delta = subtractJulianDay(julianday1, julianday2, newJulianDelta('+', 0, 0));
  ck_assert(delta != NULL);

  const char *format = "subtractJulianDay failed: expected:(%c, %lld, %lld) actual:(%c, %lld, %lld).";
  size_t length
      = snprintf(NULL, 0, format, expected_delta_sign, expected_delta_day, expected_delta_ms, delta->sign, delta->day, delta->ms)
        + 1;
  char *msg = malloc(length);
  snprintf(msg, length, format, expected_delta_sign, expected_delta_day, expected_delta_ms, delta->sign, delta->day, delta->ms);

  ck_assert_msg(expected_delta_sign == delta->sign, msg);
  ck_assert_msg(expected_delta_day == delta->day, msg);
  ck_assert_msg(expected_delta_ms == delta->ms, msg);

  deallocateJulianDay(julianday1);
  deallocateJulianDay(julianday2);
  deallocateJulianDelta(delta);
  free(msg);
}
