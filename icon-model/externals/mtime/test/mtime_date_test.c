// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include "mtime_date_test.h"

#include <stdio.h>
#include <stdlib.h>
#include "mtime_calendar.h"

typedef struct _date *Date;

START_TEST(test_create_date_from_valid_strings)
{
  // Normal dates of today and negative years
  assertDate("2014-04-07", 2014, 4, 7, "2014-04-07", "20140407");
  assertDate("000002014-04-07", 2014, 4, 7, "2014-04-07", "20140407");
  assertDate("+2014-04-07", 2014, 4, 7, "2014-04-07", "20140407");
  assertDate("-2014-04-07", -2014, 4, 7, "-2014-04-07", "-20140407");
  assertDate("-002014-04-07", -2014, 4, 7, "-2014-04-07", "-20140407");
  assertDate("20140407", 2014, 4, 7, "2014-04-07", "20140407");
  assertDate("+20140407", 2014, 4, 7, "2014-04-07", "20140407");
  assertDate("-20140407", -2014, 4, 7, "-2014-04-07", "-20140407");
  assertDate("2014-04", 2014, 4, 1, "2014-04-01", "20140401");

  // Leap years and negative leap years
  assertDate("2012-02-29", 2012, 2, 29, "2012-02-29", "20120229");
  assertDate("-2012-02-29", -2012, 2, 29, "-2012-02-29", "-20120229");

  // Dates from early dates and negative years
  assertDate("0014-04-07", 14, 4, 7, "14-04-07", "00140407");
  assertDate("-0014-04-07", -14, 4, 7, "-14-04-07", "-00140407");
  assertDate("14-04-07", 14, 4, 7, "14-04-07", "00140407");
  assertDate("-14-04-07", -14, 4, 7, "-14-04-07", "-00140407");
  assertDate("00140407", 14, 4, 7, "14-04-07", "00140407");
  assertDate("-00140407", -14, 4, 7, "-14-04-07", "-00140407");
  assertDate("14-04", 14, 4, 1, "14-04-01", "00140401");

  // Dates from future and negative years
  assertDate("30165-04-07", 30165, 4, 7, "30165-04-07", "301650407");
  assertDate("-30165-04-07", -30165, 4, 7, "-30165-04-07", "-301650407");

  // Boundary dates
  assertDate("-2147483648-01-01", -2147483648, 1, 1, "-2147483648-01-01", "-21474836480101");
  assertDate("2147483647-12-31", 2147483647, 12, 31, "2147483647-12-31", "21474836471231");

  // Around year zero
  assertDate("-1-12-31", -1, 12, 31, "-1-12-31", "-00011231");
  assertDate("0-01-01", 0, 1, 1, "0-01-01", "00000101");
  assertDate("0-02-29", 0, 2, 29, "0-02-29", "00000229");
  assertDate("1-01-01", 1, 1, 1, "1-01-01", "00010101");
  assertDate("-00011231", -1, 12, 31, "-1-12-31", "-00011231");
  assertDate("00000101", 0, 1, 1, "0-01-01", "00000101");
  assertDate("+00000101", 0, 1, 1, "0-01-01", "00000101");
  assertDate("-00000101", 0, 1, 1, "0-01-01", "00000101");
  assertDate("00000229", 0, 2, 29, "0-02-29", "00000229");
  assertDate("00010101", 1, 1, 1, "1-01-01", "00010101");

  // Posix boundaries
  assertDate("1582-10-15", 1582, 10, 15, "1582-10-15", "15821015");
  assertDate("9999-12-31", 9999, 12, 31, "9999-12-31", "99991231");
  assertDate("1582-10-14", 1582, 10, 14, "1582-10-14", "15821014");
  assertDate("10000-01-01", 10000, 01, 01, "10000-01-01", "100000101");
}
END_TEST

START_TEST(test_invalid_strings)
{
  assertInvalidDate(NULL);
  assertInvalidDate("");
  assertInvalidDate("abc");
  assertInvalidDate("2014-05-32");
  assertInvalidDate("2014-13-15");
  assertInvalidDate("-2147483649-05-15");
  assertInvalidDate("2147483648-05-15");
  assertInvalidDate("140407");   // Not allowed by ISO
  assertInvalidDate("-140407");  // Not allowed by ISO
                                 // And many more...
}
END_TEST

START_TEST(test_constructAndCopyDate)
{
  int year = 2014;
  int month = 04;
  int day = 29;

  Date date_src = newRawDate(year, month, day);
  Date date_dest = constructAndCopyDate(date_src);
  ck_assert(date_dest != NULL);
  ck_assert_int_eq(year, date_dest->year);
  ck_assert_int_eq(month, date_dest->month);
  ck_assert_int_eq(day, date_dest->day);

  int year_neg = -2014;

  Date date_src_neg = newRawDate(year_neg, month, day);
  Date date_dest_neg = constructAndCopyDate(date_src_neg);
  ck_assert(date_dest_neg != NULL);
  ck_assert_int_eq(year_neg, date_dest_neg->year);
  ck_assert_int_eq(month, date_dest_neg->month);
  ck_assert_int_eq(day, date_dest_neg->day);

  deallocateDate(date_src);
  deallocateDate(date_dest);
  deallocateDate(date_src_neg);
  deallocateDate(date_dest_neg);
}
END_TEST

START_TEST(test_replaceDate)
{
  int year = 2014;
  int month = 04;
  int day = 29;

  Date date_src = newRawDate(year, month, day);
  Date date_dest = newDate("0-01-01");
  replaceDate(date_src, date_dest);
  ck_assert(date_dest != NULL);
  ck_assert_int_eq(year, date_dest->year);
  ck_assert_int_eq(month, date_dest->month);
  ck_assert_int_eq(day, date_dest->day);

  int year_neg = -2014;

  Date date_src_neg = newRawDate(year_neg, month, day);
  Date date_dest_neg = newDate("0-01-01");
  replaceDate(date_src_neg, date_dest_neg);
  ck_assert(date_dest_neg != NULL);
  ck_assert_int_eq(year_neg, date_dest_neg->year);
  ck_assert_int_eq(month, date_dest_neg->month);
  ck_assert_int_eq(day, date_dest_neg->day);

  deallocateDate(date_src);
  deallocateDate(date_dest);
  deallocateDate(date_src_neg);
  deallocateDate(date_dest_neg);
}
END_TEST

START_TEST(test_dateToPosixString)
{
  char *tmp = malloc(MAX_DATE_STR_LEN);
  char *format1 = "%x";
  char *format2 = "%d.%m.%Y";

  Date dateA = newDate("2012-02-29");
  ck_assert_str_eq("02/29/12", dateToPosixString(dateA, tmp, format1));
  ck_assert_str_eq("29.02.2012", dateToPosixString(dateA, tmp, format2));

  Date dateB = newDate("1582-10-15");
  ck_assert_str_eq("10/15/82", dateToPosixString(dateB, tmp, format1));
  ck_assert_str_eq("15.10.1582", dateToPosixString(dateB, tmp, format2));

  Date dateC = newDate("9999-12-31");
  ck_assert_str_eq("12/31/99", dateToPosixString(dateC, tmp, format1));
  ck_assert_str_eq("31.12.9999", dateToPosixString(dateC, tmp, format2));

  deallocateDate(dateA);
  deallocateDate(dateB);
  deallocateDate(dateC);
  free(tmp);
}
END_TEST

START_TEST(test_compareDate_1)
{
  Date date1 = newDate("-2147483648-01-01");
  Date date2 = newDate("-2012-02-29");
  Date date3 = newDate("-1-12-31");
  Date date4 = newDate("0-01-01");
  Date date5 = newDate("2012-02-29");
  Date date6 = newDate("2147483647-01-01");

  ck_assert_int_eq(equal_to, compareDate(date1, date1));
  ck_assert_int_eq(less_than, compareDate(date1, date2));
  ck_assert_int_eq(less_than, compareDate(date1, date3));
  ck_assert_int_eq(less_than, compareDate(date1, date4));
  ck_assert_int_eq(less_than, compareDate(date1, date5));
  ck_assert_int_eq(less_than, compareDate(date1, date6));

  ck_assert_int_eq(greater_than, compareDate(date2, date1));
  ck_assert_int_eq(equal_to, compareDate(date2, date2));
  ck_assert_int_eq(less_than, compareDate(date2, date3));
  ck_assert_int_eq(less_than, compareDate(date2, date4));
  ck_assert_int_eq(less_than, compareDate(date2, date5));
  ck_assert_int_eq(less_than, compareDate(date2, date6));

  ck_assert_int_eq(greater_than, compareDate(date3, date1));
  ck_assert_int_eq(greater_than, compareDate(date3, date2));
  ck_assert_int_eq(equal_to, compareDate(date3, date3));
  ck_assert_int_eq(less_than, compareDate(date3, date4));
  ck_assert_int_eq(less_than, compareDate(date3, date5));
  ck_assert_int_eq(less_than, compareDate(date3, date6));

  ck_assert_int_eq(greater_than, compareDate(date4, date1));
  ck_assert_int_eq(greater_than, compareDate(date4, date2));
  ck_assert_int_eq(greater_than, compareDate(date4, date3));
  ck_assert_int_eq(equal_to, compareDate(date4, date4));
  ck_assert_int_eq(less_than, compareDate(date4, date5));
  ck_assert_int_eq(less_than, compareDate(date4, date6));

  ck_assert_int_eq(greater_than, compareDate(date5, date1));
  ck_assert_int_eq(greater_than, compareDate(date5, date2));
  ck_assert_int_eq(greater_than, compareDate(date5, date3));
  ck_assert_int_eq(greater_than, compareDate(date5, date4));
  ck_assert_int_eq(equal_to, compareDate(date5, date5));
  ck_assert_int_eq(less_than, compareDate(date5, date6));

  ck_assert_int_eq(greater_than, compareDate(date6, date1));
  ck_assert_int_eq(greater_than, compareDate(date6, date2));
  ck_assert_int_eq(greater_than, compareDate(date6, date3));
  ck_assert_int_eq(greater_than, compareDate(date6, date4));
  ck_assert_int_eq(greater_than, compareDate(date6, date5));
  ck_assert_int_eq(equal_to, compareDate(date6, date6));

  deallocateDate(date1);
  deallocateDate(date2);
  deallocateDate(date3);
  deallocateDate(date4);
  deallocateDate(date5);
  deallocateDate(date6);
}
END_TEST

START_TEST(test_compareDate_2)
{
  Date dateA = newDate("2012-02-29");
  Date dateB = newDate("2012-02-28");
  Date dateC = newDate("2012-01-29");
  Date dateD = newDate("2011-02-29");

  ck_assert_int_eq(greater_than, compareDate(dateA, dateB));
  ck_assert_int_eq(greater_than, compareDate(dateA, dateC));
  ck_assert_int_eq(greater_than, compareDate(dateA, dateD));

  ck_assert_int_eq(less_than, compareDate(dateB, dateA));
  ck_assert_int_eq(greater_than, compareDate(dateB, dateC));
  ck_assert_int_eq(greater_than, compareDate(dateB, dateD));

  ck_assert_int_eq(less_than, compareDate(dateC, dateA));
  ck_assert_int_eq(less_than, compareDate(dateC, dateB));
  ck_assert_int_eq(greater_than, compareDate(dateC, dateD));

  ck_assert_int_eq(less_than, compareDate(dateD, dateA));
  ck_assert_int_eq(less_than, compareDate(dateD, dateB));
  ck_assert_int_eq(less_than, compareDate(dateD, dateC));

  deallocateDate(dateA);
  deallocateDate(dateB);
  deallocateDate(dateC);
  deallocateDate(dateD);
}
END_TEST

static void
setup(void)
{
  initCalendar(PROLEPTIC_GREGORIAN);
}

static void
teardown(void)
{
  freeCalendar();
}

void
add_mtime_date_test_to_suite(Suite *suite)
{
  TCase *tcase = tcase_create("mtime_date_test");
  suite_add_tcase(suite, tcase);
  tcase_add_checked_fixture(tcase, setup, teardown);
  tcase_add_test(tcase, test_create_date_from_valid_strings);
  tcase_add_test(tcase, test_invalid_strings);
  tcase_add_test(tcase, test_constructAndCopyDate);
  tcase_add_test(tcase, test_replaceDate);
  tcase_add_test(tcase, test_dateToPosixString);
  tcase_add_test(tcase, test_compareDate_1);
  tcase_add_test(tcase, test_compareDate_2);
}

/*** SPECIAL ASSERT FUNCTIONS ***/

void
assertDate(const char *input_string, int64_t year, int month, int day, const char *expected_output_string,
           const char *expected_basic_string)
{
  const char *format = "Parsing of date string \"%s\" failed.";
  size_t length = snprintf(NULL, 0, format, input_string) + 1;
  char *msg = malloc(length);
  snprintf(msg, length, format, input_string);

  char tmp[MAX_DATE_STR_LEN];

  Date date = newDate(input_string);
  ck_assert_msg(date != NULL, msg);
  ck_assert_str_eq(expected_output_string, dateToString(date, tmp));
  ck_assert_str_eq(expected_output_string, tmp);
  ck_assert_int_eq(year, date->year);
  ck_assert_int_eq(month, date->month);
  ck_assert_int_eq(day, date->day);
  ck_assert_str_eq(expected_basic_string, dateToBasicString(date, tmp));
  ck_assert_str_eq(expected_basic_string, tmp);

  Date raw_date = newRawDate(year, month, day);
  ck_assert(raw_date != NULL);
  ck_assert_str_eq(expected_output_string, dateToString(raw_date, tmp));
  ck_assert_int_eq(year, raw_date->year);
  ck_assert_int_eq(month, raw_date->month);
  ck_assert_int_eq(day, raw_date->day);
  ck_assert_str_eq(expected_basic_string, dateToBasicString(raw_date, tmp));

  Date date_from_string = newDate(dateToString(raw_date, tmp));
  ck_assert(date_from_string != NULL);
  ck_assert_str_eq(expected_output_string, dateToString(date_from_string, tmp));

  if (year >= -9999 && year <= 9999)
    {
      Date date_from_basic_string = newDate(dateToBasicString(raw_date, tmp));
      ck_assert(date_from_basic_string != NULL);
      ck_assert_str_eq(expected_basic_string, dateToBasicString(date_from_basic_string, tmp));
      deallocateDate(date_from_basic_string);
    }

  deallocateDate(date);
  deallocateDate(raw_date);
  deallocateDate(date_from_string);
  free(msg);
}

void
assertInvalidDate(const char *input_string)
{
  Date date = newDate(input_string);

  if (date != NULL)
    {
      const char *format = "Invalid date string \"%s\" didn't result in null pointer but in date \"%s\".";
      char output_string[MAX_DATE_STR_LEN];
      dateToString(date, output_string);

      size_t length = snprintf(NULL, 0, format, input_string, output_string) + 1;
      char *msg = malloc(length);
      snprintf(msg, length, format, input_string, output_string);
      deallocateDate(date);

      ck_abort_msg(msg);
      free(msg);  // Never executed?
    }
}
