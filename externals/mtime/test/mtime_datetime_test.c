// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include "mtime_datetime_test.h"

#include <stdio.h>
#include <stdlib.h>
#include "mtime_calendar.h"
#include "mtime_date.h"
#include "mtime_julianDay.h"

typedef struct _datetime* DateTime;
typedef struct _date* Date;

START_TEST(test_create_datetime_from_valid_strings)
{
	// Normal datetimes of today and negative years
	assertDateTime("1980-04-07T05:42:23.109", 1980, 4, 7, 5, 42, 23, 109, "1980-04-07T05:42:23.109", "19800407T054223.109");
	assertDateTime("1982-04-07 05:42:23.109", 1982, 4, 7, 5, 42, 23, 109, "1982-04-07T05:42:23.109", "19820407T054223.109");
	assertDateTime("1983-04-07T05:42:23.109", 1983, 4, 7, 5, 42, 23, 109, "1983-04-07T05:42:23.109", "19830407T054223.109");
	assertDateTime("1984-04-07 05:42:23", 1984, 4, 7, 5, 42, 23, 0, "1984-04-07T05:42:23.000", "19840407T054223.000");
	assertDateTime("1985-04-07 05:42", 1985, 4, 7, 5, 42, 0, 0, "1985-04-07T05:42:00.000", "19850407T054200.000");
	assertDateTime("1986-04-07 05", 1986, 4, 7, 5, 0, 0, 0, "1986-04-07T05:00:00.000", "19860407T050000.000");
	assertDateTime("1987-04-07", 1987, 4, 7, 0, 0, 0, 0, "1987-04-07T00:00:00.000", "19870407T000000.000");
	assertDateTime("19890407T054200.456", 1989, 4, 7, 5, 42, 0, 456, "1989-04-07T05:42:00.456", "19890407T054200.456");
	assertDateTime("19900407 054200.456", 1990, 4, 7, 5, 42, 0, 456, "1990-04-07T05:42:00.456", "19900407T054200.456");
	assertDateTime("19910407", 1991, 4, 7, 0, 0, 0, 0, "1991-04-07T00:00:00.000", "19910407T000000.000");
	assertDateTime("00001992-04-07T05:42:23.109", 1992, 4, 7, 5, 42, 23, 109, "1992-04-07T05:42:23.109", "19920407T054223.109");
	assertDateTime("-1995-04-07T05:42:23.109", -1995, 4, 7, 5, 42, 23, 109, "-1995-04-07T05:42:23.109", "-19950407T054223.109");
	assertDateTime("-1997-04-07 05:42:23.109", -1997, 4, 7, 5, 42, 23, 109, "-1997-04-07T05:42:23.109", "-19970407T054223.109");
	assertDateTime("-1998-04-07T05:42:23.109", -1998, 4, 7, 5, 42, 23, 109, "-1998-04-07T05:42:23.109", "-19980407T054223.109");
	assertDateTime("-1999-04-07T05:42:23", -1999, 4, 7, 5, 42, 23, 0, "-1999-04-07T05:42:23.000", "-19990407T054223.000");
	assertDateTime("-2000-04-07T05:42", -2000, 4, 7, 5, 42, 0, 0, "-2000-04-07T05:42:00.000", "-20000407T054200.000");
	assertDateTime("-2001-04-07T05", -2001, 4, 7, 5, 0, 0, 0, "-2001-04-07T05:00:00.000", "-20010407T050000.000");
	assertDateTime("-2002-04-07", -2002, 4, 7, 0, 0, 0, 0, "-2002-04-07T00:00:00.000", "-20020407T000000.000");
	assertDateTime("-20040407T054200.456", -2004, 4, 7, 5, 42, 0, 456, "-2004-04-07T05:42:00.456", "-20040407T054200.456");
	assertDateTime("-20040407T054200", -2004, 4, 7, 5, 42, 0, 0, "-2004-04-07T05:42:00.000", "-20040407T054200.000");
	assertDateTime("-20040407T0542", -2004, 4, 7, 5, 42, 0, 0, "-2004-04-07T05:42:00.000", "-20040407T054200.000");
	assertDateTime("-20040407T05", -2004, 4, 7, 5, 0, 0, 0, "-2004-04-07T05:00:00.000", "-20040407T050000.000");
	assertDateTime("-20050407 054200.456", -2005, 4, 7, 5, 42, 0, 456, "-2005-04-07T05:42:00.456", "-20050407T054200.456");
	assertDateTime("-20060423", -2006, 4, 23, 0, 0, 0, 0, "-2006-04-23T00:00:00.000", "-20060423T000000.000");
	assertDateTime("-00002007-04-07T05:42:23.109", -2007, 4, 7, 5, 42, 23, 109, "-2007-04-07T05:42:23.109", "-20070407T054223.109");

	// Leap years and negative leap years
	assertDateTime("2012-02-29", 2012, 2, 29, 0, 0, 0, 0, "2012-02-29T00:00:00.000", "20120229T000000.000");
	assertDateTime("-2012-02-29", -2012, 2, 29, 0, 0, 0, 0, "-2012-02-29T00:00:00.000", "-20120229T000000.000");

	// Datetimes from early dates and negative years
	assertDateTime("352-04-07T05:42:23.101", 352, 4, 7, 5, 42, 23, 101, "352-04-07T05:42:23.101", "03520407T054223.101");
	assertDateTime("352-04-07 05:42:23.103", 352, 4, 7, 5, 42, 23, 103, "352-04-07T05:42:23.103", "03520407T054223.103");
	assertDateTime("352-04-07T05:42:23.104", 352, 4, 7, 5, 42, 23, 104, "352-04-07T05:42:23.104", "03520407T054223.104");
	assertDateTime("352-04-07 05:42:05", 352, 4, 7, 5, 42, 5, 0, "352-04-07T05:42:05.000", "03520407T054205.000");
	assertDateTime("352-04-07 05:06", 352, 4, 7, 5, 6, 0, 0, "352-04-07T05:06:00.000", "03520407T050600.000");
	assertDateTime("352-04-07 07", 352, 4, 7, 7, 0, 0, 0, "352-04-07T07:00:00.000", "03520407T070000.000");
	assertDateTime("352-04-08", 352, 4, 8, 0, 0, 0, 0, "352-04-08T00:00:00.000", "03520408T000000.000");
	assertDateTime("03520407T054200.110", 352, 4, 7, 5, 42, 0, 110, "352-04-07T05:42:00.110", "03520407T054200.110");
	assertDateTime("03520407 054200.111", 352, 4, 7, 5, 42, 0, 111, "352-04-07T05:42:00.111", "03520407T054200.111");
	assertDateTime("03520412", 352, 4, 12, 0, 0, 0, 0, "352-04-12T00:00:00.000", "03520412T000000.000");
	assertDateTime("0000352-04-07T05:42:23.113", 352, 4, 7, 5, 42, 23, 113, "352-04-07T05:42:23.113", "03520407T054223.113");
	assertDateTime("-352-04-07T05:42:23.116", -352, 4, 7, 5, 42, 23, 116, "-352-04-07T05:42:23.116", "-03520407T054223.116");
	assertDateTime("-03520407T054200.117", -352, 4, 7, 5, 42, 0, 117, "-352-04-07T05:42:00.117", "-03520407T054200.117");

	// Datetimes from future and _negative years
	assertDateTime("568028-04-07T05:42:23.101", 568028, 4, 7, 5, 42, 23, 101, "568028-04-07T05:42:23.101", "5680280407T054223.101");
	assertDateTime("568028-04-07 05:42:23.103", 568028, 4, 7, 5, 42, 23, 103, "568028-04-07T05:42:23.103", "5680280407T054223.103");
	assertDateTime("568028-04-07T05:42:23.104", 568028, 4, 7, 5, 42, 23, 104, "568028-04-07T05:42:23.104", "5680280407T054223.104");
	assertDateTime("568028-04-07 05:42:05", 568028, 4, 7, 5, 42, 5, 0, "568028-04-07T05:42:05.000", "5680280407T054205.000");
	assertDateTime("568028-04-07 05:06", 568028, 4, 7, 5, 6, 0, 0, "568028-04-07T05:06:00.000", "5680280407T050600.000");
	assertDateTime("568028-04-07 07", 568028, 4, 7, 7, 0, 0, 0, "568028-04-07T07:00:00.000", "5680280407T070000.000");
	assertDateTime("568028-04-08", 568028, 4, 8, 0, 0, 0, 0, "568028-04-08T00:00:00.000", "5680280408T000000.000");
	assertDateTime("0000568028-04-07T05:42:23.113", 568028, 4, 7, 5, 42, 23, 113, "568028-04-07T05:42:23.113", "5680280407T054223.113");
	assertDateTime("-568028-04-07T05:42:23.114", -568028, 4, 7, 5, 42, 23, 114, "-568028-04-07T05:42:23.114", "-5680280407T054223.114");
	
	//Boundary datetimes
	assertDateTime("-2147483648-01-01 00:00:00.000", -2147483648, 1, 1, 0, 0, 0, 0, "-2147483648-01-01T00:00:00.000", "-21474836480101T000000.000");
	assertDateTime("2147483647-12-31 23:59:59.999", 2147483647, 12, 31, 23, 59, 59, 999, "2147483647-12-31T23:59:59.999", "21474836471231T235959.999");
	
	// Around year zero
	assertDateTime("-1-12-31 23:59:59.999", -1, 12, 31, 23, 59, 59, 999, "-1-12-31T23:59:59.999", "-00011231T235959.999");
	assertDateTime("0-01-01 00:00:00.000", 0, 1, 1, 0, 0 ,0 ,0, "0-01-01T00:00:00.000", "00000101T000000.000");
	assertDateTime("0-02-29 13:28:46.832", 0, 2, 29, 13, 28, 46, 832, "0-02-29T13:28:46.832", "00000229T132846.832");
	assertDateTime("1-01-01 00:00:00.000", 1, 1, 1, 0, 0 ,0 ,0, "1-01-01T00:00:00.000", "00010101T000000.000");
	assertDateTime("-00011231T235959.999", -1, 12, 31, 23, 59, 59, 999, "-1-12-31T23:59:59.999", "-00011231T235959.999");
	assertDateTime("00000101 000000.000", 0, 1, 1, 0, 0 ,0 ,0, "0-01-01T00:00:00.000", "00000101T000000.000");
	assertDateTime("00000229 132846.832", 0, 2, 29, 13, 28, 46, 832, "0-02-29T13:28:46.832", "00000229T132846.832");
	assertDateTime("00010101 000000.000", 1, 1, 1, 0, 0 ,0 ,0, "1-01-01T00:00:00.000", "00010101T000000.000");

	// Posix boundaries
	assertDateTime("1582-10-15 00:00:00.000", 1582, 10, 15, 0, 0, 0, 0, "1582-10-15T00:00:00.000", "15821015T000000.000");
	assertDateTime("9999-12-31 23:59:59.999", 9999, 12, 31, 23, 59, 59, 999, "9999-12-31T23:59:59.999", "99991231T235959.999");
	assertDateTime("1582-10-14 23:59:59.999", 1582, 10, 14, 23, 59, 59, 999, "1582-10-14T23:59:59.999", "15821014T235959.999");
	assertDateTime("10000-01-01 00:00:00.000", 10000, 01, 01, 0, 0, 0, 0, "10000-01-01T00:00:00.000", "100000101T000000.000");
}
END_TEST

START_TEST(test_invalid_strings)
{
	assertInvalidDateTime(NULL);
	assertInvalidDateTime("");
	assertInvalidDateTime("abc");
	assertInvalidDateTime("2014-05-15T13:28:56.1000");
	assertInvalidDateTime("2014-05-15T13:28:60.123");
	assertInvalidDateTime("2014-05-15T13:60:56.123");
	assertInvalidDateTime("2014-05-15T24:28:56.123");
	assertInvalidDateTime("2014-05-32T13:28:56.132");
	assertInvalidDateTime("2014-13-15T13:28:56.132");
	assertInvalidDateTime("-2147483649-05-15T13:28:56.132");
	assertInvalidDateTime("2147483648-05-15T13:28:56.132");
	assertInvalidDateTime("2014-05-15T24:00:00.000");
	// And many more...
}
END_TEST

START_TEST(test_constructAndCopyDateTime)
{
	int year = 2014;
	int month = 04;
	int day = 29;
	int hour = 13;
	int minute = 4;
	int second = 29;
	int ms = 678;
	
	DateTime datetime_src = newRawDateTime(year, month, day, hour, minute, second, ms);
	DateTime datetime_dest = constructAndCopyDateTime(datetime_src);
	ck_assert(datetime_dest != NULL);
	ck_assert_int_eq(year, datetime_dest->date.year);
	ck_assert_int_eq(month, datetime_dest->date.month);
	ck_assert_int_eq(day, datetime_dest->date.day);
	ck_assert_int_eq(hour, datetime_dest->time.hour);
	ck_assert_int_eq(minute, datetime_dest->time.minute);
	ck_assert_int_eq(second, datetime_dest->time.second);
	ck_assert_int_eq(ms, datetime_dest->time.ms);	
	
	int year_neg = -2014;
	
	DateTime datetime_src_neg = newRawDateTime(year_neg, month, day, hour, minute, second, ms);
	DateTime datetime_dest_neg = constructAndCopyDateTime(datetime_src_neg);
	ck_assert(datetime_dest_neg != NULL);
	ck_assert_int_eq(year_neg, datetime_dest_neg->date.year);
	ck_assert_int_eq(month, datetime_dest_neg->date.month);
	ck_assert_int_eq(day, datetime_dest_neg->date.day);
	ck_assert_int_eq(hour, datetime_dest_neg->time.hour);
	ck_assert_int_eq(minute, datetime_dest_neg->time.minute);
	ck_assert_int_eq(second, datetime_dest_neg->time.second);
	ck_assert_int_eq(ms, datetime_dest_neg->time.ms);	
	
	deallocateDateTime(datetime_src);
	deallocateDateTime(datetime_dest);
	deallocateDateTime(datetime_src_neg);
	deallocateDateTime(datetime_dest_neg);
}
END_TEST

START_TEST(test_replaceDatetime)
{
	int year = 2014;
	int month = 04;
	int day = 29;
	int hour = 13;
	int minute = 4;
	int second = 29;
	int ms = 678;
	
	DateTime datetime_src = newRawDateTime(year, month, day, hour, minute, second, ms);
	DateTime datetime_dest = newDateTime("0-01-01T00:00:00.000");
	replaceDatetime(datetime_src, datetime_dest);
	ck_assert(datetime_dest != NULL);
	ck_assert_int_eq(year, datetime_dest->date.year);
	ck_assert_int_eq(month, datetime_dest->date.month);
	ck_assert_int_eq(day, datetime_dest->date.day);
	ck_assert_int_eq(hour, datetime_dest->time.hour);
	ck_assert_int_eq(minute, datetime_dest->time.minute);
	ck_assert_int_eq(second, datetime_dest->time.second);
	ck_assert_int_eq(ms, datetime_dest->time.ms);	
	
	int year_neg = -2014;
	
	DateTime datetime_src_neg = newRawDateTime(year_neg, month, day, hour, minute, second, ms);
	DateTime datetime_dest_neg = newDateTime("0-01-01T00:00:00.000");
	replaceDatetime(datetime_src_neg, datetime_dest_neg);
	ck_assert(datetime_dest_neg != NULL);
	ck_assert_int_eq(year_neg, datetime_dest_neg->date.year);
	ck_assert_int_eq(month, datetime_dest_neg->date.month);
	ck_assert_int_eq(day, datetime_dest_neg->date.day);
	ck_assert_int_eq(hour, datetime_dest_neg->time.hour);
	ck_assert_int_eq(minute, datetime_dest_neg->time.minute);
	ck_assert_int_eq(second, datetime_dest_neg->time.second);
	ck_assert_int_eq(ms, datetime_dest_neg->time.ms);	
	
	deallocateDateTime(datetime_src);
	deallocateDateTime(datetime_dest);
	deallocateDateTime(datetime_src_neg);
	deallocateDateTime(datetime_dest_neg);
}
END_TEST

START_TEST(test_datetimeToPosixString)
{
	char tmp[MAX_DATE_STR_LEN];
	char* format1 = "%x";
	char* format2 = "%X";
	char* format3 = "%d.%m.%Y %H:%M:%S";
	
	DateTime datetimeA = newDateTime("2012-02-29T14:04:59.666");
	ck_assert_str_eq("02/29/12", datetimeToPosixString(datetimeA, tmp, format1));
	ck_assert_str_eq("14:04:59", datetimeToPosixString(datetimeA, tmp, format2));
	ck_assert_str_eq("29.02.2012 14:04:59", datetimeToPosixString(datetimeA, tmp, format3));	
	
	DateTime datetimeB = newDateTime("1582-10-15T00:00:00.000");
	ck_assert_str_eq("10/15/82", datetimeToPosixString(datetimeB, tmp, format1));
	ck_assert_str_eq("00:00:00", datetimeToPosixString(datetimeB, tmp, format2));
	ck_assert_str_eq("15.10.1582 00:00:00", datetimeToPosixString(datetimeB, tmp, format3));	
	
	DateTime datetimeC = newDateTime("9999-12-31T23:59:59.999");
	ck_assert_str_eq("12/31/99", datetimeToPosixString(datetimeC, tmp, format1));
	ck_assert_str_eq("23:59:59", datetimeToPosixString(datetimeC, tmp, format2));
	ck_assert_str_eq("31.12.9999 23:59:59", datetimeToPosixString(datetimeC, tmp, format3));	
	
	deallocateDateTime(datetimeA);
	deallocateDateTime(datetimeB);
	deallocateDateTime(datetimeC);
}
END_TEST

START_TEST(test_testYearIsLeapYear)
{
	ck_assert(testYearIsLeapYear(-2147483648));
	ck_assert(!testYearIsLeapYear(-2147483647));
	ck_assert(testYearIsLeapYear(-2147483600));
	ck_assert(!testYearIsLeapYear(-2147483500));
	ck_assert(!testYearIsLeapYear(-2014));
	ck_assert(testYearIsLeapYear(-2012));
	ck_assert(testYearIsLeapYear(-2000));
	ck_assert(!testYearIsLeapYear(-1800));
	ck_assert(testYearIsLeapYear(-400));
	ck_assert(!testYearIsLeapYear(-100));
	ck_assert(testYearIsLeapYear(-4));
	ck_assert(testYearIsLeapYear(0));
	ck_assert(testYearIsLeapYear(4));
	ck_assert(!testYearIsLeapYear(100));
	ck_assert(testYearIsLeapYear(400));
	ck_assert(!testYearIsLeapYear(1800));
	ck_assert(testYearIsLeapYear(2000));
	ck_assert(testYearIsLeapYear(2012));
	ck_assert(!testYearIsLeapYear(2014));
	ck_assert(!testYearIsLeapYear(2147483500));
	ck_assert(testYearIsLeapYear(2147483600));
	ck_assert(!testYearIsLeapYear(2147483647));
}
END_TEST

START_TEST(test_compareDatetime_1)
{
	DateTime datetime1 = newDateTime("-2147483648-01-01 00:00:00.000");
	DateTime datetime2 = newDateTime("-2012-02-29 11:12:13.214");
	DateTime datetime3 = newDateTime("-1-12-31 23:59:59.999");
	DateTime datetime4 = newDateTime("0-01-01 00:00:00.000");
	DateTime datetime5 = newDateTime("2012-02-29 11:12:13.214");
	DateTime datetime6 = newDateTime("2147483647-01-01 23:59:59.999");
	
	ck_assert_int_eq(equal_to,     compareDatetime(datetime1, datetime1));
	ck_assert_int_eq(less_than,    compareDatetime(datetime1, datetime2));
	ck_assert_int_eq(less_than,    compareDatetime(datetime1, datetime3));
	ck_assert_int_eq(less_than,    compareDatetime(datetime1, datetime4));
	ck_assert_int_eq(less_than,    compareDatetime(datetime1, datetime5));
	ck_assert_int_eq(less_than,    compareDatetime(datetime1, datetime6));
	
	ck_assert_int_eq(greater_than, compareDatetime(datetime2, datetime1));
	ck_assert_int_eq(equal_to,     compareDatetime(datetime2, datetime2));
	ck_assert_int_eq(less_than,    compareDatetime(datetime2, datetime3));
	ck_assert_int_eq(less_than,    compareDatetime(datetime2, datetime4));
	ck_assert_int_eq(less_than,    compareDatetime(datetime2, datetime5));
	ck_assert_int_eq(less_than,    compareDatetime(datetime2, datetime6));
	
	ck_assert_int_eq(greater_than, compareDatetime(datetime3, datetime1));
	ck_assert_int_eq(greater_than, compareDatetime(datetime3, datetime2));
	ck_assert_int_eq(equal_to,     compareDatetime(datetime3, datetime3));
	ck_assert_int_eq(less_than,    compareDatetime(datetime3, datetime4));
	ck_assert_int_eq(less_than,    compareDatetime(datetime3, datetime5));
	ck_assert_int_eq(less_than,    compareDatetime(datetime3, datetime6));
	
	ck_assert_int_eq(greater_than, compareDatetime(datetime4, datetime1));
	ck_assert_int_eq(greater_than, compareDatetime(datetime4, datetime2));
	ck_assert_int_eq(greater_than, compareDatetime(datetime4, datetime3));
	ck_assert_int_eq(equal_to,     compareDatetime(datetime4, datetime4));
	ck_assert_int_eq(less_than,   compareDatetime(datetime4, datetime5));
	ck_assert_int_eq(less_than,    compareDatetime(datetime4, datetime6));
	
	ck_assert_int_eq(greater_than, compareDatetime(datetime5, datetime1));
	ck_assert_int_eq(greater_than, compareDatetime(datetime5, datetime2));
	ck_assert_int_eq(greater_than, compareDatetime(datetime5, datetime3));
	ck_assert_int_eq(greater_than, compareDatetime(datetime5, datetime4));
	ck_assert_int_eq(equal_to,     compareDatetime(datetime5, datetime5));
	ck_assert_int_eq(less_than,    compareDatetime(datetime5, datetime6));
	
	ck_assert_int_eq(greater_than, compareDatetime(datetime6, datetime1));
	ck_assert_int_eq(greater_than, compareDatetime(datetime6, datetime2));
	ck_assert_int_eq(greater_than, compareDatetime(datetime6, datetime3));
	ck_assert_int_eq(greater_than, compareDatetime(datetime6, datetime4));
	ck_assert_int_eq(greater_than, compareDatetime(datetime6, datetime5));
	ck_assert_int_eq(equal_to,     compareDatetime(datetime6, datetime6));
	
	deallocateDateTime(datetime1);
	deallocateDateTime(datetime2);
	deallocateDateTime(datetime3);
	deallocateDateTime(datetime4);
	deallocateDateTime(datetime5);
	deallocateDateTime(datetime6);
}
END_TEST

START_TEST(test_compareDatetime_2)
{
	DateTime datetimeA = newDateTime("2012-02-29 20:15:46.567");
	DateTime datetimeB = newDateTime("2012-02-29 20:15:46.566");
	DateTime datetimeC = newDateTime("2012-02-29 20:15:45.567");
	DateTime datetimeD = newDateTime("2012-02-29 20:14:46.567");
	DateTime datetimeE = newDateTime("2012-02-29 19:15:46.567");
	DateTime datetimeF = newDateTime("2012-02-28 20:15:46.567");
	DateTime datetimeG = newDateTime("2012-01-29 20:15:46.567");
	DateTime datetimeH = newDateTime("2011-02-29 20:15:46.567");
	
	ck_assert_int_eq(greater_than, compareDatetime(datetimeA, datetimeB));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeA, datetimeC));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeA, datetimeD));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeA, datetimeE));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeA, datetimeF));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeA, datetimeG));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeA, datetimeH));
	
	ck_assert_int_eq(less_than,    compareDatetime(datetimeB, datetimeA));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeB, datetimeC));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeB, datetimeD));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeB, datetimeE));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeB, datetimeF));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeB, datetimeG));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeB, datetimeH));
	
	ck_assert_int_eq(less_than,    compareDatetime(datetimeC, datetimeA));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeC, datetimeB));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeC, datetimeD));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeC, datetimeE));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeC, datetimeF));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeC, datetimeG));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeC, datetimeH));

	ck_assert_int_eq(less_than,    compareDatetime(datetimeD, datetimeA));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeD, datetimeB));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeD, datetimeC));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeD, datetimeE));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeD, datetimeF));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeD, datetimeG));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeD, datetimeH));
	
	ck_assert_int_eq(less_than,    compareDatetime(datetimeE, datetimeA));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeE, datetimeB));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeE, datetimeC));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeE, datetimeD));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeE, datetimeF));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeE, datetimeG));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeE, datetimeH));
	
	ck_assert_int_eq(less_than,    compareDatetime(datetimeF, datetimeA));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeF, datetimeB));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeF, datetimeC));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeF, datetimeD));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeF, datetimeE));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeF, datetimeG));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeF, datetimeH));
	
	ck_assert_int_eq(less_than,    compareDatetime(datetimeG, datetimeA));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeG, datetimeB));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeG, datetimeC));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeG, datetimeD));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeG, datetimeE));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeG, datetimeF));
	ck_assert_int_eq(greater_than, compareDatetime(datetimeG, datetimeH));
	
	ck_assert_int_eq(less_than,    compareDatetime(datetimeH, datetimeA));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeH, datetimeB));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeH, datetimeC));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeH, datetimeD));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeH, datetimeE));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeH, datetimeF));
	ck_assert_int_eq(less_than,    compareDatetime(datetimeH, datetimeG));
	
	deallocateDateTime(datetimeA);
	deallocateDateTime(datetimeB);
	deallocateDateTime(datetimeC);
	deallocateDateTime(datetimeD);
	deallocateDateTime(datetimeE);
	deallocateDateTime(datetimeF);
	deallocateDateTime(datetimeG);
	deallocateDateTime(datetimeH);
}
END_TEST

START_TEST(test_getNoOfDaysInMonthDateTime_ProlepticGregorian)
{
	assertNoDaysInMonth("2014-01-01", 31);
	assertNoDaysInMonth("2014-02-01", 28);
	assertNoDaysInMonth("2014-03-01", 31);
	assertNoDaysInMonth("2014-04-01", 30);
	assertNoDaysInMonth("2014-05-01", 31);
	assertNoDaysInMonth("2014-06-01", 30);
	assertNoDaysInMonth("2014-07-01", 31);
	assertNoDaysInMonth("2014-08-01", 31);
	assertNoDaysInMonth("2014-09-01", 30);
	assertNoDaysInMonth("2014-10-01", 31);
	assertNoDaysInMonth("2014-11-01", 30);
	assertNoDaysInMonth("2014-12-01", 31);
	
	assertNoDaysInMonth("-2012-02-01", 29); 
	assertNoDaysInMonth("-2000-02-01", 29); 
	assertNoDaysInMonth("-1800-02-01", 28); 
	assertNoDaysInMonth("0-02-01",     29); 
	assertNoDaysInMonth("1800-02-01",  28); 
	assertNoDaysInMonth("2000-02-01",  29); 
	assertNoDaysInMonth("2012-02-01",  29); 
}
END_TEST

START_TEST(test_getNoOfDaysInMonthDateTime_YearOf365Days)
{
	assertNoDaysInMonth("2014-01-01", 31);
	assertNoDaysInMonth("2014-02-01", 28);
	assertNoDaysInMonth("2014-03-01", 31);
	assertNoDaysInMonth("2014-04-01", 30);
	assertNoDaysInMonth("2014-05-01", 31);
	assertNoDaysInMonth("2014-06-01", 30);
	assertNoDaysInMonth("2014-07-01", 31);
	assertNoDaysInMonth("2014-08-01", 31);
	assertNoDaysInMonth("2014-09-01", 30);
	assertNoDaysInMonth("2014-10-01", 31);
	assertNoDaysInMonth("2014-11-01", 30);
	assertNoDaysInMonth("2014-12-01", 31);
	
	assertNoDaysInMonth("-2012-02-01", 28); 
	assertNoDaysInMonth("-2000-02-01", 28); 
	assertNoDaysInMonth("-1800-02-01", 28); 
	assertNoDaysInMonth("0-02-01",     28); 
	assertNoDaysInMonth("1800-02-01",  28); 
	assertNoDaysInMonth("2000-02-01",  28); 
	assertNoDaysInMonth("2012-02-01",  28); 
}
END_TEST

START_TEST(test_getNoOfDaysInMonthDateTime_YearOf360Days)
{
	assertNoDaysInMonth("2014-01-01", 30);
	assertNoDaysInMonth("2014-02-01", 30);
	assertNoDaysInMonth("2014-03-01", 30);
	assertNoDaysInMonth("2014-04-01", 30);
	assertNoDaysInMonth("2014-05-01", 30);
	assertNoDaysInMonth("2014-06-01", 30);
	assertNoDaysInMonth("2014-07-01", 30);
	assertNoDaysInMonth("2014-08-01", 30);
	assertNoDaysInMonth("2014-09-01", 30);
	assertNoDaysInMonth("2014-10-01", 30);
	assertNoDaysInMonth("2014-11-01", 30);
	assertNoDaysInMonth("2014-12-01", 30);
	
	assertNoDaysInMonth("-2012-02-01", 30); 
	assertNoDaysInMonth("-2000-02-01", 30); 
	assertNoDaysInMonth("-1800-02-01", 30); 
	assertNoDaysInMonth("0-02-01",     30); 
	assertNoDaysInMonth("1800-02-01",  30); 
	assertNoDaysInMonth("2000-02-01",  30); 
	assertNoDaysInMonth("2012-02-01",  30); 
}
END_TEST

START_TEST(test_getNoOfDaysInYearDateTime_ProlepticGregorian)
{
	assertNoDaysInYear("-2147483648-01-01", 366);
	assertNoDaysInYear("-2147483647-01-01", 365); 
	assertNoDaysInYear("-2147483600-01-01", 366); 
	assertNoDaysInYear("-2147483500-01-01", 365); 
	assertNoDaysInYear("-2014-01-01", 365); 
	assertNoDaysInYear("-2012-01-01", 366); 
	assertNoDaysInYear("-2000-01-01", 366); 
	assertNoDaysInYear("-1800-01-01", 365); 
	assertNoDaysInYear("-400-01-01", 366); 
	assertNoDaysInYear("-100-01-01", 365); 
	assertNoDaysInYear("-4-01-01", 366); 
	assertNoDaysInYear("0-01-01", 366); 
	assertNoDaysInYear("4-01-01", 366); 
	assertNoDaysInYear("100-01-01", 365); 
	assertNoDaysInYear("400-01-01", 366); 
	assertNoDaysInYear("1800-01-01", 365); 
	assertNoDaysInYear("2000-01-01", 366); 
	assertNoDaysInYear("2012-01-01", 366); 
	assertNoDaysInYear("2014-01-01", 365); 
	assertNoDaysInYear("2147483500-01-01", 365); 
	assertNoDaysInYear("2147483600-01-01", 366); 
	assertNoDaysInYear("2147483644-01-01", 366); 
	assertNoDaysInYear("2147483647-01-01", 365); 
}
END_TEST

START_TEST(test_getNoOfDaysInYearDateTime_YearOf365Days)
{
	assertNoDaysInYear("-2147483648-01-01", 365);
	assertNoDaysInYear("-2147483647-01-01", 365); 
	assertNoDaysInYear("-2147483600-01-01", 365); 
	assertNoDaysInYear("-2147483500-01-01", 365); 
	assertNoDaysInYear("-2014-01-01", 365); 
	assertNoDaysInYear("-2012-01-01", 365); 
	assertNoDaysInYear("-2000-01-01", 365); 
	assertNoDaysInYear("-1800-01-01", 365); 
	assertNoDaysInYear("-400-01-01", 365); 
	assertNoDaysInYear("-100-01-01", 365); 
	assertNoDaysInYear("-4-01-01", 365); 
	assertNoDaysInYear("0-01-01", 365); 
	assertNoDaysInYear("4-01-01", 365); 
	assertNoDaysInYear("100-01-01", 365); 
	assertNoDaysInYear("400-01-01", 365); 
	assertNoDaysInYear("1800-01-01", 365); 
	assertNoDaysInYear("2000-01-01", 365); 
	assertNoDaysInYear("2012-01-01", 365); 
	assertNoDaysInYear("2014-01-01", 365); 
	assertNoDaysInYear("2147483500-01-01", 365); 
	assertNoDaysInYear("2147483600-01-01", 365); 
	assertNoDaysInYear("2147483644-01-01", 365); 
	assertNoDaysInYear("2147483647-01-01", 365); 
}
END_TEST

START_TEST(test_getNoOfDaysInYearDateTime_YearOf360Days)
{
	assertNoDaysInYear("-2147483648-01-01", 360);
	assertNoDaysInYear("-2147483647-01-01", 360); 
	assertNoDaysInYear("-2147483600-01-01", 360); 
	assertNoDaysInYear("-2147483500-01-01", 360); 
	assertNoDaysInYear("-2014-01-01", 360); 
	assertNoDaysInYear("-2012-01-01", 360); 
	assertNoDaysInYear("-2000-01-01", 360); 
	assertNoDaysInYear("-1800-01-01", 360); 
	assertNoDaysInYear("-400-01-01", 360); 
	assertNoDaysInYear("-100-01-01", 360); 
	assertNoDaysInYear("-4-01-01", 360); 
	assertNoDaysInYear("0-01-01", 360); 
	assertNoDaysInYear("4-01-01", 360); 
	assertNoDaysInYear("100-01-01", 360); 
	assertNoDaysInYear("400-01-01", 360); 
	assertNoDaysInYear("1800-01-01", 360); 
	assertNoDaysInYear("2000-01-01", 360); 
	assertNoDaysInYear("2012-01-01", 360); 
	assertNoDaysInYear("2014-01-01", 360); 
	assertNoDaysInYear("2147483500-01-01", 360); 
	assertNoDaysInYear("2147483600-01-01", 360); 
	assertNoDaysInYear("2147483644-01-01", 360); 
	assertNoDaysInYear("2147483647-01-01", 360); 
}
END_TEST

START_TEST(test_getDayOfYearFromDateTime_ProlepticGregorian)
{
	assertDayOfYear("2014-01-01", 1);
	assertDayOfYear("2014-02-28", 59);
	assertDayOfYear("2014-03-01", 60);
	assertDayOfYear("2014-05-07", 127);
	assertDayOfYear("2014-12-30", 364);
	assertDayOfYear("2014-12-31", 365);
	
	assertDayOfYear("-2014-01-01", 1);
	assertDayOfYear("-2014-02-28", 59);
	assertDayOfYear("-2014-03-01", 60);
	assertDayOfYear("-2014-05-07", 127);
	assertDayOfYear("-2014-12-30", 364);
	assertDayOfYear("-2014-12-31", 365);
	
	assertDayOfYear("2012-01-01", 1);
	assertDayOfYear("2012-02-28", 59);
	assertDayOfYear("2012-02-29", 60);
	assertDayOfYear("2012-03-01", 61);
	assertDayOfYear("2012-05-07", 128);
	assertDayOfYear("2012-12-30", 365);
	assertDayOfYear("2012-12-31", 366);
}
END_TEST

START_TEST(test_getDayOfYearFromDateTime_YearOf365Days)
{
	assertDayOfYear("2014-01-01", 1);
	assertDayOfYear("2014-02-28", 59);
	assertDayOfYear("2014-03-01", 60);
	assertDayOfYear("2014-05-07", 127);
	assertDayOfYear("2014-12-30", 364);
	assertDayOfYear("2014-12-31", 365);
	
	assertDayOfYear("-2014-01-01", 1);
	assertDayOfYear("-2014-02-28", 59);
	assertDayOfYear("-2014-03-01", 60);
	assertDayOfYear("-2014-05-07", 127);
	assertDayOfYear("-2014-12-30", 364);
	assertDayOfYear("-2014-12-31", 365);
	
	assertDayOfYear("2012-01-01", 1);
	assertDayOfYear("2012-02-28", 59);
	assertDayOfYear("2012-03-01", 60);
	assertDayOfYear("2012-05-07", 127);
	assertDayOfYear("2012-12-30", 364);
	assertDayOfYear("2012-12-31", 365);
}
END_TEST

START_TEST(test_getDayOfYearFromDateTime_YearOf360Days)
{
	assertDayOfYear("2014-01-01", 1);
	assertDayOfYear("2014-02-28", 58);
	assertDayOfYear("2014-03-01", 61);
	assertDayOfYear("2014-05-07", 127);
	assertDayOfYear("2014-12-30", 360);
	
	assertDayOfYear("-2014-01-01", 1);
	assertDayOfYear("-2014-02-28", 58);
	assertDayOfYear("-2014-03-01", 61);
	assertDayOfYear("-2014-05-07", 127);
	assertDayOfYear("-2014-12-30", 360);
	
	assertDayOfYear("2012-01-01", 1);
	assertDayOfYear("2012-02-28", 58);
	assertDayOfYear("2012-02-29", 59);
	assertDayOfYear("2012-03-01", 61);
	assertDayOfYear("2012-05-07", 127);
	assertDayOfYear("2012-12-30", 360);
}
END_TEST

START_TEST(test_getNoOfSecondsElapsedInMonthDateTime_ProlepticGregorian)
{
	assertNoOfSecondsInMonth("2014-01-01 00:00:00.000", 0);
	assertNoOfSecondsInMonth("2014-01-01 00:00:00.999", 0);
	assertNoOfSecondsInMonth("2014-01-01 00:00:01.000", 1);
	assertNoOfSecondsInMonth("2014-01-01 01:00:00.444", 3600);
	assertNoOfSecondsInMonth("2014-01-01 23:59:59.999", 86399);
	assertNoOfSecondsInMonth("2014-01-02 00:00:00.000", 86400);
	assertNoOfSecondsInMonth("2014-01-02 00:00:00.000", 86400);
	assertNoOfSecondsInMonth("2014-01-31 23:59:59.999", 2678399);
	assertNoOfSecondsInMonth("2014-09-17 16:45:23.798", 1442723);
	assertNoOfSecondsInMonth("-2014-09-17 16:45:23.798", 1442723);
	assertNoOfSecondsInMonth("2014-12-31 23:59:59.999", 2678399); 
	assertNoOfSecondsInMonth("2012-02-29 23:59:59.999", 2505599);
}
END_TEST

START_TEST(test_getNoOfSecondsElapsedInMonthDateTime_YearOf365Days)
{
	assertNoOfSecondsInMonth("2014-01-01 00:00:00.000", 0);
	assertNoOfSecondsInMonth("2014-01-01 00:00:00.999", 0);
	assertNoOfSecondsInMonth("2014-01-01 00:00:01.000", 1);
	assertNoOfSecondsInMonth("2014-01-01 01:00:00.444", 3600);
	assertNoOfSecondsInMonth("2014-01-01 23:59:59.999", 86399);
	assertNoOfSecondsInMonth("2014-01-02 00:00:00.000", 86400);
	assertNoOfSecondsInMonth("2014-01-02 00:00:00.000", 86400);
	assertNoOfSecondsInMonth("2014-01-31 23:59:59.999", 2678399);
	assertNoOfSecondsInMonth("2014-09-17 16:45:23.798", 1442723);
	assertNoOfSecondsInMonth("-2014-09-17 16:45:23.798", 1442723);
	assertNoOfSecondsInMonth("2014-12-31 23:59:59.999", 2678399); 
}
END_TEST

START_TEST(test_getNoOfSecondsElapsedInMonthDateTime_YearOf360Days)
{
	assertNoOfSecondsInMonth("2014-01-01 00:00:00.000", 0);
	assertNoOfSecondsInMonth("2014-01-01 00:00:00.999", 0);
	assertNoOfSecondsInMonth("2014-01-01 00:00:01.000", 1);
	assertNoOfSecondsInMonth("2014-01-01 01:00:00.444", 3600);
	assertNoOfSecondsInMonth("2014-01-01 23:59:59.999", 86399);
	assertNoOfSecondsInMonth("2014-01-02 00:00:00.000", 86400);
	assertNoOfSecondsInMonth("2014-01-30 23:59:59.999", 2591999);
	assertNoOfSecondsInMonth("2014-09-17 16:45:23.798", 1442723);
	assertNoOfSecondsInMonth("-2014-09-17 16:45:23.798", 1442723);
	assertNoOfSecondsInMonth("2014-12-30 23:59:59.999", 2591999); 
	assertNoOfSecondsInMonth("2014-02-30 23:59:59.999", 2591999);
}
END_TEST

START_TEST(test_getNoOfSecondsElapsedInDayDateTime_ProlepticGregorian)
{
	assertNoOfSecondsInDay("2014-01-01 00:00:00.000", 0);
	assertNoOfSecondsInDay("2014-01-01 00:00:00.999", 0);
	assertNoOfSecondsInDay("2014-01-01 00:00:01.000", 1);
	assertNoOfSecondsInDay("2014-01-01 01:00:00.444", 3600);
	assertNoOfSecondsInDay("2014-01-01 23:59:59.999", 86399);
	assertNoOfSecondsInDay("2014-09-17 16:45:23.798", 60323);
	assertNoOfSecondsInDay("-2014-09-17 16:45:23.798", 60323);
	assertNoOfSecondsInDay("2014-12-31 23:59:59.999", 86399); 
	assertNoOfSecondsInDay("2012-02-29 23:59:59.999", 86399);
}
END_TEST

START_TEST(test_getNoOfSecondsElapsedInDayDateTime_YearOf365Days)
{
	assertNoOfSecondsInDay("2014-01-01 00:00:00.000", 0);
	assertNoOfSecondsInDay("2014-01-01 00:00:00.999", 0);
	assertNoOfSecondsInDay("2014-01-01 00:00:01.000", 1);
	assertNoOfSecondsInDay("2014-01-01 01:00:00.444", 3600);
	assertNoOfSecondsInDay("2014-01-01 23:59:59.999", 86399);
	assertNoOfSecondsInDay("2014-09-17 16:45:23.798", 60323);
	assertNoOfSecondsInDay("-2014-09-17 16:45:23.798", 60323);
	assertNoOfSecondsInDay("2014-12-31 23:59:59.999", 86399); 
}
END_TEST

START_TEST(test_getNoOfSecondsElapsedInDayDateTime_YearOf360Days)
{
	assertNoOfSecondsInDay("2014-01-01 00:00:00.000", 0);
	assertNoOfSecondsInDay("2014-01-01 00:00:00.999", 0);
	assertNoOfSecondsInDay("2014-01-01 00:00:01.000", 1);
	assertNoOfSecondsInDay("2014-01-01 01:00:00.444", 3600);
	assertNoOfSecondsInDay("2014-01-01 23:59:59.999", 86399);
	assertNoOfSecondsInDay("2014-09-17 16:45:23.798", 60323);
	assertNoOfSecondsInDay("-2014-09-17 16:45:23.798", 60323);
	assertNoOfSecondsInDay("2014-12-31 23:59:59.999", 86399); 
	assertNoOfSecondsInDay("2014-02-30 23:59:59.999", 86399);
}
END_TEST

START_TEST(test_getJulianDayFromDateTime_ProlepticGregorian)
{
	// No idea, how this should work
//	assertJulianDay("-2147483648-01-01 00:00:00.000", 0, 0);
//	assertJulianDay("-2014-09-17 16:45:23.798", 60323);
//	assertJulianDay("-01-12-31 23:59:59.999", 0, 0);
//	assertJulianDay("0-01-01 00:00:00.000", 0, 0);
//	assertJulianDay("2014-09-17 16:45:23.798", 60323);
//	assertJulianDay("2147483647-12-31 23:59:59.999", 0, 0);
}
END_TEST

START_TEST(test_getJulianDayFromDateTime_YearOf365Days)
{
	// No idea, how this should work
//	assertJulianDay("-2147483648-01-01 00:00:00.000", 0, 0);
//	assertJulianDay("-2014-09-17 16:45:23.798", 60323);
//	assertJulianDay("-01-12-31 23:59:59.999", 0, 0);
//	assertJulianDay("0-01-01 00:00:00.000", 0, 0);
//	assertJulianDay("2014-09-17 16:45:23.798", 60323);
//	assertJulianDay("2147483647-12-31 23:59:59.999", 0, 0);
}
END_TEST

START_TEST(test_getJulianDayFromDateTime_YearOf360Days)
{
	// No idea, how this should work
//	assertJulianDay("-2147483648-01-01 00:00:00.000", 0, 0);
//	assertJulianDay("-2014-09-17 16:45:23.798", 60323);
//	assertJulianDay("-01-12-31 23:59:59.999", 0, 0);
//	assertJulianDay("0-01-01 00:00:00.000", 0, 0);
//	assertJulianDay("2014-09-17 16:45:23.798", 60323);
//	assertJulianDay("2147483647-12-31 23:59:59.999", 0, 0);
}
END_TEST

START_TEST(test_getDateTimeIsInRange)
{
	DateTime datetime1 = newDateTime("-2147483648-01-01 00:00:00.000");
	DateTime datetime2 = newDateTime("-2012-02-29 11:12:13.214");
	DateTime datetime3 = newDateTime("-1-12-31 23:59:59.999");
	DateTime datetime4 = newDateTime("0-01-01 00:00:00.000");
	DateTime datetime5 = newDateTime("2012-02-29 11:12:13.214");
	DateTime datetime6 = newDateTime("2147483647-01-01 23:59:59.999");
	
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime1, datetime1, datetime6));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime2, datetime1, datetime6));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime3, datetime1, datetime6));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime4, datetime1, datetime6));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime5, datetime1, datetime6));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime6, datetime1, datetime6));
	
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime1, datetime2, datetime3));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime2, datetime2, datetime3));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime3, datetime2, datetime3));
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime4, datetime2, datetime3));
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime5, datetime2, datetime3));
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime6, datetime2, datetime3));
	
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime1, datetime2, datetime4));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime2, datetime2, datetime4));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime3, datetime2, datetime4));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime4, datetime2, datetime4));
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime5, datetime2, datetime4));
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime6, datetime2, datetime4));
	
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime1, datetime4, datetime2));
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime2, datetime4, datetime2));
//	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime3, datetime4, datetime2)); //FAILS (according to documentation)
//	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime4, datetime4, datetime2)); //FAILS (according to documentation)
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime5, datetime4, datetime2));
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime6, datetime4, datetime2));
	
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime1, datetime5, datetime5));
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime2, datetime5, datetime5));
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime3, datetime5, datetime5));
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime4, datetime5, datetime5));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime5, datetime5, datetime5));
	ck_assert_int_eq(greater_than, getDateTimeIsInRange(datetime6, datetime5, datetime5));
	
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime1, datetime5, datetime6));
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime2, datetime5, datetime6));
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime3, datetime5, datetime6));
	ck_assert_int_eq(less_than,    getDateTimeIsInRange(datetime4, datetime5, datetime6));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime5, datetime5, datetime6));
	ck_assert_int_eq(equal_to,     getDateTimeIsInRange(datetime6, datetime5, datetime6));
	
	deallocateDateTime(datetime1);
	deallocateDateTime(datetime2);
	deallocateDateTime(datetime3);
	deallocateDateTime(datetime4);
	deallocateDateTime(datetime5);
	deallocateDateTime(datetime6);
}
END_TEST

START_TEST(test_convertDateToDateTime)
{
	int year = 2014;
	int month = 04;
	int day = 29;

	Date date = newRawDate(year, month, day);
	DateTime datetime = newDateTime("0-01-01T12:34:50.123");
	datetime = convertDateToDateTime(date, datetime);
	ck_assert(datetime != NULL);
	ck_assert_int_eq(year, datetime->date.year);
	ck_assert_int_eq(month, datetime->date.month);
	ck_assert_int_eq(day, datetime->date.day);
	ck_assert_int_eq(0, datetime->time.hour);
	ck_assert_int_eq(0, datetime->time.minute);
	ck_assert_int_eq(0, datetime->time.second);
	ck_assert_int_eq(0, datetime->time.ms);

	int year_neg = -2014;

	Date date_neg = newRawDate(year_neg, month, day);
	DateTime datetime_neg = newDateTime("0-01-01T12:34:50.123");
	datetime_neg = convertDateToDateTime(date_neg, datetime_neg);
	ck_assert(datetime_neg != NULL);
	ck_assert_int_eq(year_neg, datetime_neg->date.year);
	ck_assert_int_eq(month, datetime_neg->date.month);
	ck_assert_int_eq(day, datetime_neg->date.day);
	ck_assert_int_eq(0, datetime_neg->time.hour);
	ck_assert_int_eq(0, datetime_neg->time.minute);
	ck_assert_int_eq(0, datetime_neg->time.second);
	ck_assert_int_eq(0, datetime_neg->time.ms);

	deallocateDate(date);
	deallocateDateTime(datetime);
	deallocateDate(date_neg);
	deallocateDateTime(datetime_neg);
}
END_TEST

START_TEST(test_convertDateTimeToDate)
{
	int year = 2014;
	int month = 04;
	int day = 29;

	DateTime datetime = newRawDateTime(year, month, day, 23, 42, 7, 109);
	Date date = newDate("0-01-01");
	date = convertDateTimeToDate(datetime, date);
	ck_assert(date != NULL);
	ck_assert_int_eq(year, date->year);
	ck_assert_int_eq(month, date->month);
	ck_assert_int_eq(day, date->day);

	int year_neg = -2014;

	DateTime datetime_neg = newRawDateTime(year_neg, month, day, 23, 42, 7, 109);
	Date date_neg = newDate("0-01-01");
	date_neg = convertDateTimeToDate(datetime_neg, date_neg);
	ck_assert(date_neg != NULL);
	ck_assert_int_eq(year_neg, date_neg->year);
	ck_assert_int_eq(month, date_neg->month);
	ck_assert_int_eq(day, date_neg->day);

	deallocateDate(date);
	deallocateDateTime(datetime);
	deallocateDate(date_neg);
	deallocateDateTime(datetime_neg);
}
END_TEST

static void setup_ProlepticGregorian(void)
{
	initCalendar(PROLEPTIC_GREGORIAN);
}

static void setup_YearOf365Days(void)
{
	initCalendar(YEAR_OF_365_DAYS);
}

static void setup_YearOf360Days(void)
{
	initCalendar(YEAR_OF_360_DAYS);
}

static void teardown(void)
{
	freeCalendar();
}

void add_mtime_datetime_test_to_suite(Suite* suite)
{
    TCase *tcase_ProlepticGregorian = tcase_create("mtime_datetime_test_ProlepticGregorian");
    suite_add_tcase(suite, tcase_ProlepticGregorian);
    tcase_add_checked_fixture(tcase_ProlepticGregorian, setup_ProlepticGregorian, teardown);
    tcase_add_test(tcase_ProlepticGregorian, test_create_datetime_from_valid_strings);
    tcase_add_test(tcase_ProlepticGregorian, test_invalid_strings);
    tcase_add_test(tcase_ProlepticGregorian, test_constructAndCopyDateTime);
    tcase_add_test(tcase_ProlepticGregorian, test_replaceDatetime);
    tcase_add_test(tcase_ProlepticGregorian, test_datetimeToPosixString);
    tcase_add_test(tcase_ProlepticGregorian, test_testYearIsLeapYear);
    tcase_add_test(tcase_ProlepticGregorian, test_compareDatetime_1);
    tcase_add_test(tcase_ProlepticGregorian, test_compareDatetime_2);
    tcase_add_test(tcase_ProlepticGregorian, test_getNoOfDaysInMonthDateTime_ProlepticGregorian);
    tcase_add_test(tcase_ProlepticGregorian, test_getNoOfDaysInYearDateTime_ProlepticGregorian);
    tcase_add_test(tcase_ProlepticGregorian, test_getDayOfYearFromDateTime_ProlepticGregorian);
    tcase_add_test(tcase_ProlepticGregorian, test_getNoOfSecondsElapsedInMonthDateTime_ProlepticGregorian);
    tcase_add_test(tcase_ProlepticGregorian, test_getNoOfSecondsElapsedInDayDateTime_ProlepticGregorian);
    tcase_add_test(tcase_ProlepticGregorian, test_getJulianDayFromDateTime_ProlepticGregorian);
    tcase_add_test(tcase_ProlepticGregorian, test_getDateTimeIsInRange);
    tcase_add_test(tcase_ProlepticGregorian, test_convertDateToDateTime);
    tcase_add_test(tcase_ProlepticGregorian, test_convertDateTimeToDate);

    TCase *tcase_YearOf365Days = tcase_create("mtime_datetime_test_YearOf365Days");
    suite_add_tcase(suite, tcase_YearOf365Days);
    tcase_add_checked_fixture(tcase_YearOf365Days, setup_YearOf365Days, teardown);
    tcase_add_test(tcase_YearOf365Days, test_getNoOfDaysInMonthDateTime_YearOf365Days);
    tcase_add_test(tcase_YearOf365Days, test_getNoOfDaysInYearDateTime_YearOf365Days);
    tcase_add_test(tcase_YearOf365Days, test_getDayOfYearFromDateTime_YearOf365Days);
    tcase_add_test(tcase_YearOf365Days, test_getNoOfSecondsElapsedInMonthDateTime_YearOf365Days);
    tcase_add_test(tcase_YearOf365Days, test_getNoOfSecondsElapsedInDayDateTime_YearOf365Days);
    tcase_add_test(tcase_YearOf365Days, test_getJulianDayFromDateTime_YearOf365Days);

    TCase *tcase_YearOf360Days = tcase_create("mtime_datetime_test_YearOf360Days");
    suite_add_tcase(suite, tcase_YearOf360Days);
    tcase_add_checked_fixture(tcase_YearOf360Days, setup_YearOf360Days, teardown);
    tcase_add_test(tcase_YearOf360Days, test_getNoOfDaysInMonthDateTime_YearOf360Days);
    tcase_add_test(tcase_YearOf360Days, test_getNoOfDaysInYearDateTime_YearOf360Days);
    tcase_add_test(tcase_YearOf360Days, test_getDayOfYearFromDateTime_YearOf360Days);
    tcase_add_test(tcase_YearOf360Days, test_getNoOfSecondsElapsedInMonthDateTime_YearOf360Days);
    tcase_add_test(tcase_YearOf360Days, test_getNoOfSecondsElapsedInDayDateTime_YearOf360Days);
    tcase_add_test(tcase_YearOf360Days, test_getJulianDayFromDateTime_YearOf360Days);
}

/*** SPECIAL ASSERT FUNCTIONS ***/

void assertDateTime(const char* input_string, int64_t year, int month, int day, int hour, int minute, int second, int ms, const char* expected_output_string, const char* expected_basic_string)
{
	const char* format = "Parsing of datetime string \"%s\" failed.";
	size_t length = snprintf(NULL, 0, format, input_string) + 1;
	char* msg = malloc(length);
	snprintf(msg, length, format, input_string);

	char tmp[MAX_DATE_STR_LEN];

	DateTime datetime = newDateTime(input_string);
	ck_assert_msg(datetime != NULL, msg);
	ck_assert_str_eq(expected_output_string, datetimeToString(datetime, tmp));
	ck_assert_str_eq(expected_output_string, tmp);
	ck_assert_int_eq(year, datetime->date.year);
	ck_assert_int_eq(month, datetime->date.month);
	ck_assert_int_eq(day, datetime->date.day);
	ck_assert_int_eq(hour, datetime->time.hour);
	ck_assert_int_eq(minute, datetime->time.minute);
	ck_assert_int_eq(second, datetime->time.second);
	ck_assert_int_eq(ms, datetime->time.ms);
	ck_assert_str_eq(expected_basic_string, datetimeToBasicString(datetime, tmp));
	ck_assert_str_eq(expected_basic_string, tmp);

	DateTime raw_datetime = newRawDateTime(year, month, day, hour, minute, second, ms);
	ck_assert(raw_datetime != NULL);
	ck_assert_str_eq(expected_output_string, datetimeToString(raw_datetime, tmp));
	ck_assert_int_eq(year, raw_datetime->date.year);
	ck_assert_int_eq(month, raw_datetime->date.month);
	ck_assert_int_eq(day, raw_datetime->date.day);
	ck_assert_int_eq(hour, raw_datetime->time.hour);
	ck_assert_int_eq(minute, raw_datetime->time.minute);
	ck_assert_int_eq(second, raw_datetime->time.second);
	ck_assert_int_eq(ms, raw_datetime->time.ms);
	ck_assert_str_eq(expected_basic_string, datetimeToBasicString(raw_datetime, tmp));

	DateTime datetime_from_string = newDateTime(datetimeToString(raw_datetime, tmp));
	ck_assert(datetime_from_string != NULL);
	ck_assert_str_eq(expected_output_string, datetimeToString(datetime_from_string, tmp));

	if (year >= -9999 && year <= 9999)
	{
		DateTime datetime_from_basic_string = newDateTime(datetimeToBasicString(raw_datetime, tmp));
		ck_assert(datetime_from_basic_string != NULL);
		ck_assert_str_eq(expected_basic_string, datetimeToBasicString(datetime_from_basic_string, tmp));
		deallocateDateTime(datetime_from_basic_string);
	}

	deallocateDateTime(datetime);
	deallocateDateTime(raw_datetime);
	deallocateDateTime(datetime_from_string);
	free(msg);
}


void assertInvalidDateTime(const char* input_string)
{
	DateTime datetime = newDateTime(input_string);

	if (datetime != NULL)
	{
		const char* format = "Invalid datetime string \"%s\" didn't result in null pointer but in datetime \"%s\".";
		char output_string[MAX_DATE_STR_LEN];
		datetimeToString(datetime, output_string);

		size_t length = snprintf(NULL, 0, format, input_string, output_string) + 1;
		char* msg = malloc(length);
		snprintf(msg, length, format, input_string, output_string);
		deallocateDateTime(datetime);

		ck_abort_msg(msg);
		free(msg); //Never executed?
	}
}

void genericAssertDatetimeToInteger(const char* datetime_string, int expected_value, int getNoOf(DateTime dt), const char* function_name)
{
	DateTime datetime = newDateTime(datetime_string);
	int actual_value = getNoOf(datetime);

	char cal[MAX_CALENDAR_STR_LEN];
	calendarToString(cal);

	const char* format = "%s for \"%s\" failed: expected:%d, actual:%d [%s].";
	size_t length = snprintf(NULL, 0, format, function_name, datetime_string, expected_value, actual_value, cal) + 1;
	char* msg = malloc(length);
	snprintf(msg, length, format, function_name, datetime_string, expected_value, actual_value, cal);

	ck_assert_msg(expected_value == actual_value, msg);

	deallocateDateTime(datetime);
	free(msg);
}

void assertNoDaysInMonth(const char* datetime_string, int expected_days)
{
	genericAssertDatetimeToInteger(datetime_string, expected_days, getNoOfDaysInMonthDateTime, "getNoOfDaysInMonthDateTime");
}

void assertNoDaysInYear(const char* datetime_string, int expected_days)
{
	genericAssertDatetimeToInteger(datetime_string, expected_days, getNoOfDaysInYearDateTime, "getNoOfDaysInYearDateTime");
}

void assertDayOfYear(const char* datetime_string, int expected_day)
{
	genericAssertDatetimeToInteger(datetime_string, expected_day, getDayOfYearFromDateTime, "getDayOfYearFromDateTime");
}

void assertNoOfSecondsInDay(const char* datetime_string, int expected_days)
{
	genericAssertDatetimeToInteger(datetime_string, expected_days, getNoOfSecondsElapsedInDayDateTime, "getNoOfSecondsElapsedInDayDateTime");
}

void assertNoOfSecondsInMonth(const char* datetime_string, int64_t expected_seconds)
{
	DateTime datetime = newDateTime(datetime_string);
	int64_t actual_seconds = getNoOfSecondsElapsedInMonthDateTime(datetime);

	char cal[MAX_CALENDAR_STR_LEN];
	calendarToString(cal);

	const char* format = "getNoOfSecondsElapsedInMonthDateTime for \"%s\" failed: expected:%lld, actual:%lld [%s].";
	size_t length = snprintf(NULL, 0, format, datetime_string, expected_seconds, actual_seconds, cal) + 1;
	char* msg = malloc(length);
	snprintf(msg, length, format, datetime_string, expected_seconds, actual_seconds, cal);

	ck_assert_msg(expected_seconds == actual_seconds, msg);

	deallocateDateTime(datetime);
	free(msg);
}

void assertJulianDayFromDatetime(const char* datetime_string, int expected_day, int expected_ms)
{
	DateTime datetime = newDateTime(datetime_string);
	struct _julianday* julianday = newJulianDay(0, 0);
	julianday = getJulianDayFromDateTime(datetime, julianday);

	const char* format1 = "Converting of datetime \"%s\" to Julian Day failed.";
	size_t length1 = snprintf(NULL, 0, format1, datetime_string) + 1;
	char* msg1 = malloc(length1);
	snprintf(msg1, length1, format1, datetime_string);

	ck_assert_msg(julianday != NULL, msg1);

	char cal[MAX_CALENDAR_STR_LEN];
	calendarToString(cal);

	const char* format2 = "getJulianDayFromDateTime for \"%s\" failed: expected:(%lld, %lld) actual:(%lld, %lld) [%s].";
	size_t length2 = snprintf(NULL, 0, format2, datetime_string, expected_day, expected_ms, julianday->day, julianday->ms, cal) + 1;
	char* msg2 = malloc(length2);
	snprintf(msg2, length2, format2, datetime_string, expected_day, expected_ms, julianday->day, julianday->ms, cal);

	ck_assert_msg(expected_day == julianday->day, msg2);
	ck_assert_msg(expected_ms == julianday->ms, msg2);

	deallocateDateTime(datetime);
	deallocateJulianDay(julianday);
	free(msg1);
	free(msg2);
}

