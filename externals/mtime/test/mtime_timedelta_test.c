// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include "mtime_timedelta_test.h"

#include <stdio.h>
#include <stdlib.h>
#include "mtime_calendar.h"
#include "mtime_datetime.h"
#include "mtime_julianDay.h"
#include "mtime_eventHandling.h"

typedef struct _timedelta* TimeDelta;
typedef struct _datetime* DateTime;
typedef struct _juliandelta* JulianDelta;

START_TEST(test_create_timedelta_from_valid_strings)
{
        assertTimeDelta("P0Y", true, '+', 0, 0, 0, 0, 0, 0, 0, "PT00.000S");
        assertTimeDelta("P0M", true, '+', 0, 0, 0, 0, 0, 0, 0, "PT00.000S");
        assertTimeDelta("P0D", true, '+', 0, 0, 0, 0, 0, 0, 0, "PT00.000S");
        assertTimeDelta("PT0S", true, '+', 0, 0, 0, 0, 0, 0, 0, "PT00.000S");                
        assertTimeDelta("P42Y", true, '+', 42, 0, 0, 0, 0, 0, 0, "P42Y");
	assertTimeDelta("P6M", true, '+', 0, 6, 0, 0, 0, 0, 0, "P06M");
	assertTimeDelta("P23D", true, '+', 0, 0, 23, 0, 0, 0, 0, "P23D");
	assertTimeDelta("P109D", false, '+', 0, 0, 109, 0, 0, 0, 0, "P109D");
	assertTimeDelta("+P109D", false, '+', 0, 0, 109, 0, 0, 0, 0, "P109D");
	assertTimeDelta("PT13H", true, '+', 0, 0, 0, 13, 0, 0, 0, "PT13H");
	assertTimeDelta("PT36H", false, '+', 0, 0, 0, 36, 0, 0, 0, "PT36H");
	assertTimeDelta("PT0036H", false, '+', 0, 0, 0, 36, 0, 0, 0, "PT36H");
	assertTimeDelta("PT15M", true, '+', 0, 0, 0, 0, 15, 0, 0, "PT15M");
	assertTimeDelta("-PT15M", true, '-', 0, 0, 0, 0, 15, 0, 0, "-PT15M");
	assertTimeDelta("PT75M", false, '+', 0, 0, 0, 0, 75, 0, 0, "PT75M");
	assertTimeDelta("PT8S", true, '+', 0, 0, 0, 0, 0, 8, 0, "PT08.000S");
	assertTimeDelta("PT08S", true, '+', 0, 0, 0, 0, 0, 8, 0, "PT08.000S");
	assertTimeDelta("+PT8S", true, '+', 0, 0, 0, 0, 0, 8, 0, "PT08.000S");
	assertTimeDelta("PT57S", true, '+', 0, 0, 0, 0, 0, 57, 0, "PT57.000S");
	//TODO:assertTimeDelta("-PT1024S", false, '-', 0, 0, 0, 0, 0, 1024, 0, "-PT1024.000S");
	assertTimeDelta("PT0.033S", true, '+', 0, 0, 0, 0, 0, 0, 33, "PT00.033S");

	assertTimeDelta("P504Y4M12D", true, '+', 504, 4, 12, 0, 0, 0, 0, "P504Y04M12D");
	assertTimeDelta("+P504Y4M12D", true, '+', 504, 4, 12, 0, 0, 0, 0, "P504Y04M12D");
	assertTimeDelta("-P504Y4M12D", true, '-', 504, 4, 12, 0, 0, 0, 0, "-P504Y04M12D");

	assertTimeDelta("PT18H1M13.123S", true, '+', 0, 0, 0, 18, 1, 13, 123, "PT18H01M13.123S");
	assertTimeDelta("+PT18H1M13.123S", true, '+', 0, 0, 0, 18, 1, 13, 123, "PT18H01M13.123S");
	assertTimeDelta("-PT18H1M13.123S", true, '-', 0, 0, 0, 18, 1, 13, 123, "-PT18H01M13.123S");

	assertTimeDelta("PT1M13.123S", true, '+', 0, 0, 0, 0, 1, 13, 123, "PT01M13.123S");
	assertTimeDelta("PT18H1M13.123S", true, '+', 0, 0, 0, 18, 1, 13, 123, "PT18H01M13.123S");
	assertTimeDelta("P12DT18H1M13.123S", true, '+', 0, 0, 12, 18, 1, 13, 123, "P12DT18H01M13.123S");
	assertTimeDelta("P4M12DT18H1M13.123S", true, '+', 0, 4, 12, 18, 1, 13, 123, "P04M12DT18H01M13.123S");
	assertTimeDelta("P504Y4M12DT18H1M13.123S", true, '+', 504, 4, 12, 18, 1, 13, 123, "P504Y04M12DT18H01M13.123S");
	assertTimeDelta("+P504Y4M12DT18H1M13.123S", true, '+', 504, 4, 12, 18, 1, 13, 123, "P504Y04M12DT18H01M13.123S");
	assertTimeDelta("-P504Y4M12DT18H1M13.123S", true, '-', 504, 4, 12, 18, 1, 13, 123, "-P504Y04M12DT18H01M13.123S");
	assertTimeDelta("-P504Y4M12DT18H1M13.000S", true, '-', 504, 4, 12, 18, 1, 13, 000, "-P504Y04M12DT18H01M13.000S");
	assertTimeDelta("-P504Y4M12DT18H1M", true, '-', 504, 4, 12, 18, 1, 0, 000, "-P504Y04M12DT18H01M");
	assertTimeDelta("-P504Y4M12DT18H", true, '-', 504, 4, 12, 18, 0, 0, 000, "-P504Y04M12DT18H");
	assertTimeDelta("-P504Y4M12D", true, '-', 504, 4, 12, 0, 0, 0, 000, "-P504Y04M12D");
	assertTimeDelta("-P504Y4M", true, '-', 504, 4, 0, 0, 0, 0, 000, "-P504Y04M");
	assertTimeDelta("P4M12DT18H1M", true, '+', 0, 4, 12, 18, 1, 0, 0, "P04M12DT18H01M");
	assertTimeDelta("P12DT18H", true, '+', 0, 0, 12, 18, 0, 0, 0, "P12DT18H");

	assertTimeDelta("P0Y0M0DT0H0M0.000S", true, '+', 0, 0, 0, 0, 0, 0, 0, "PT00.000S");
	assertTimeDelta("P0YT0H", true, '+', 0, 0, 0, 0, 0, 0, 0, "PT00.000S");
	assertTimeDelta("P1Y1M1DT1H1M1.001S", true, '+', 1, 1, 1, 1, 1, 1, 1, "P1Y01M01DT01H01M01.001S");
	assertTimeDelta("P2147483647Y", true, '+', 2147483647, 0, 0, 0, 0, 0, 0, "P2147483647Y");
	assertTimeDelta("-P2147483648Y", true, '-', 2147483648, 0, 0, 0, 0, 0, 0, "-P2147483648Y");
	assertTimeDelta("P2147483647D", false, '+', 0, 0, 2147483647, 0, 0, 0, 0, "P2147483647D");
	assertTimeDelta("-P2147483647D", false, '-', 0, 0, 2147483647, 0, 0, 0, 0, "-P2147483647D");
	//TODO:assertTimeDelta("PT2147483647H", false, '+', 0, 0, 0, 2147483647, 0, 0, 0, "PT2147483647H");
	//TODO:assertTimeDelta("-PT2147483647H", false, '-', 0, 0, 0, 2147483647, 0, 0, 0, "-PT2147483647H");
	//TODO:assertTimeDelta("PT2147483647M", false, '+', 0, 0, 0, 0, 2147483647, 0, 0, "PT2147483647M");
	//TODO:assertTimeDelta("-PT2147483647M", false, '-', 0, 0, 0, 0, 2147483647, 0, 0, "-PT2147483647M");
        //TODO:assertTimeDelta("PT2147483647S", false, '+', 0, 0, 0, 0, 0, 2147483647, 0, "PT2147483647.000S");
	//TODO:assertTimeDelta("-PT2147483647S", false, '-', 0, 0, 0, 0, 0, 2147483647, 0, "-PT2147483647.000S");
	assertTimeDelta("PT0.999S", true, '+', 0, 0, 0, 0, 0, 0, 999, "PT00.999S");
	assertTimeDelta("-PT0.999S", true, '-', 0, 0, 0, 0, 0, 0, 999, "-PT00.999S");
	assertTimeDelta("PT86760S", false, '+', 0, 0, 0, 0, 0, 86760, 0, "PT86760.000S");        
}
END_TEST

START_TEST(test_invalid_strings)
{
	assertInvalidTimeDelta(NULL);
	assertInvalidTimeDelta("");
	assertInvalidTimeDelta("abc");
	assertInvalidTimeDelta("P3Z");
	assertInvalidTimeDelta("P42M");
	// And many more...
}
END_TEST

START_TEST(test_constructAndCopyTimeDelta)
{
	int64_t year = 504;
	int month = 4;
	int day = 12;
	int hour = 18;
	int minute = 1;
	int second = 13;
	int ms = 123;

	TimeDelta timedelta_src = newRawTimeDelta('+', year, month, day, hour, minute, second, ms);
	TimeDelta timedelta_dest = constructAndCopyTimeDelta(timedelta_src);
	ck_assert(timedelta_dest != NULL);
	ck_assert(timedelta_dest->flag_std_form);
	ck_assert_int_eq('+', timedelta_dest->sign);
	ck_assert_int_eq(year, timedelta_dest->year);
	ck_assert_int_eq(month, timedelta_dest->month);
	ck_assert_int_eq(day, timedelta_dest->day);
	ck_assert_int_eq(hour, timedelta_dest->hour);
	ck_assert_int_eq(minute, timedelta_dest->minute);
	ck_assert_int_eq(second, timedelta_dest->second);
	ck_assert_int_eq(ms, timedelta_dest->ms);
        
	TimeDelta timedelta_src_neg = newRawTimeDelta('-', year, month, day, hour, minute, second, ms);
        TimeDelta timedelta_dest_neg = constructAndCopyTimeDelta(timedelta_src_neg);
	ck_assert(timedelta_dest_neg != NULL);
	ck_assert(timedelta_dest_neg->flag_std_form);
	ck_assert_int_eq('-', timedelta_dest_neg->sign);
	ck_assert_int_eq(year, timedelta_dest_neg->year);
	ck_assert_int_eq(month, timedelta_dest_neg->month);
	ck_assert_int_eq(day, timedelta_dest_neg->day);
	ck_assert_int_eq(hour, timedelta_dest_neg->hour);
	ck_assert_int_eq(minute, timedelta_dest_neg->minute);
	ck_assert_int_eq(second, timedelta_dest_neg->second);
	ck_assert_int_eq(ms, timedelta_dest_neg->ms);

        TimeDelta timedelta_src_nonstd;
        TimeDelta timedelta_dest_nonstd;

        timedelta_src_nonstd = newRawTimeDelta('-', 0, 0, 2000, 0, 0, 0, 0);
        timedelta_dest_nonstd = constructAndCopyTimeDelta(timedelta_src_nonstd);
	ck_assert(timedelta_dest_nonstd != NULL);
	ck_assert(!timedelta_dest_nonstd->flag_std_form);
	ck_assert_int_eq('-', timedelta_dest_nonstd->sign);
	ck_assert_int_eq(0, timedelta_dest_nonstd->year);
	ck_assert_int_eq(0, timedelta_dest_nonstd->month);
	ck_assert_int_eq(2000, timedelta_dest_nonstd->day);
	ck_assert_int_eq(0, timedelta_dest_nonstd->hour);
	ck_assert_int_eq(0, timedelta_dest_nonstd->minute);
	ck_assert_int_eq(0, timedelta_dest_nonstd->second);
	ck_assert_int_eq(0, timedelta_dest_nonstd->ms);

	deallocateTimeDelta(timedelta_src_nonstd);
	deallocateTimeDelta(timedelta_dest_nonstd);

        timedelta_src_nonstd = newRawTimeDelta('-', 0, 0, 0, 0, 0, 5760000, 0);
        timedelta_dest_nonstd = constructAndCopyTimeDelta(timedelta_src_nonstd);
	ck_assert(timedelta_dest_nonstd != NULL);
	ck_assert(!timedelta_dest_nonstd->flag_std_form);
	ck_assert_int_eq('-', timedelta_dest_nonstd->sign);
	ck_assert_int_eq(0, timedelta_dest_nonstd->year);
	ck_assert_int_eq(0, timedelta_dest_nonstd->month);
	ck_assert_int_eq(0, timedelta_dest_nonstd->day);
	ck_assert_int_eq(0, timedelta_dest_nonstd->hour);
	ck_assert_int_eq(0, timedelta_dest_nonstd->minute);
	ck_assert_int_eq(5760000, timedelta_dest_nonstd->second);
	ck_assert_int_eq(0, timedelta_dest_nonstd->ms);

	deallocateTimeDelta(timedelta_src);
	deallocateTimeDelta(timedelta_dest);
	deallocateTimeDelta(timedelta_src_neg);
	deallocateTimeDelta(timedelta_dest_neg);
	deallocateTimeDelta(timedelta_src_nonstd);
	deallocateTimeDelta(timedelta_dest_nonstd);
}
END_TEST

START_TEST(test_replaceTimeDelta)
{
	int64_t year = 504;
	int month = 4;
	int day = 12;
	int hour = 18;
	int minute = 1;
	int second = 13;
	int ms = 123;

	TimeDelta timedelta_src = newRawTimeDelta('+', year, month, day, hour, minute, second, ms);
	TimeDelta timedelta_dest = newTimeDelta("PT0S");
	replaceTimeDelta(timedelta_src, timedelta_dest);
	ck_assert(timedelta_dest != NULL);
	ck_assert_int_eq('+', timedelta_dest->sign);
	ck_assert_int_eq(year, timedelta_dest->year);
	ck_assert_int_eq(month, timedelta_dest->month);
	ck_assert_int_eq(day, timedelta_dest->day);
	ck_assert_int_eq(hour, timedelta_dest->hour);
	ck_assert_int_eq(minute, timedelta_dest->minute);
	ck_assert_int_eq(second, timedelta_dest->second);
	ck_assert_int_eq(ms, timedelta_dest->ms);

	TimeDelta timedelta_src_neg = newRawTimeDelta('-', year, month, day, hour, minute, second, ms);
	TimeDelta timedelta_dest_neg = newTimeDelta("P0Y");
	replaceTimeDelta(timedelta_src_neg, timedelta_dest_neg);
	ck_assert(timedelta_dest_neg != NULL);
	ck_assert_int_eq('-', timedelta_dest_neg->sign);
	ck_assert_int_eq(year, timedelta_dest_neg->year);
	ck_assert_int_eq(month, timedelta_dest_neg->month);
	ck_assert_int_eq(day, timedelta_dest_neg->day);
	ck_assert_int_eq(hour, timedelta_dest_neg->hour);
	ck_assert_int_eq(minute, timedelta_dest_neg->minute);
	ck_assert_int_eq(second, timedelta_dest_neg->second);
	ck_assert_int_eq(ms, timedelta_dest_neg->ms);

	TimeDelta timedelta_src_nonstd = newRawTimeDelta('-', 0, 0, 2000, 0, 0, 0, 0);
	TimeDelta timedelta_dest_nonstd = newTimeDelta("P0D");
	replaceTimeDelta(timedelta_src_nonstd, timedelta_dest_nonstd);
	ck_assert(timedelta_dest_nonstd != NULL);
	ck_assert(!timedelta_dest_nonstd->flag_std_form);
	ck_assert_int_eq('-', timedelta_dest_nonstd->sign);
	ck_assert_int_eq(0, timedelta_dest_nonstd->year);
	ck_assert_int_eq(0, timedelta_dest_nonstd->month);
	ck_assert_int_eq(2000, timedelta_dest_nonstd->day);
	ck_assert_int_eq(0, timedelta_dest_nonstd->hour);
	ck_assert_int_eq(0, timedelta_dest_nonstd->minute);
	ck_assert_int_eq(0, timedelta_dest_nonstd->second);
	ck_assert_int_eq(0, timedelta_dest_nonstd->ms);

	deallocateTimeDelta(timedelta_src);
	deallocateTimeDelta(timedelta_dest);
	deallocateTimeDelta(timedelta_src_neg);
	deallocateTimeDelta(timedelta_dest_neg);
	deallocateTimeDelta(timedelta_src_nonstd);
	deallocateTimeDelta(timedelta_dest_nonstd);
}
END_TEST

START_TEST(test_timeDeltaToJulianDelta_ProlepticGregorian_1)
{
	assertTimeDeltaToJulianDelta("P0Y", "0-01-01T00:00:00.000", '+', 0, 0);

	assertTimeDeltaToJulianDelta("P1Y", "-1-01-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y", "-1-01-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y", "-1-02-28T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y", "-1-03-01T00:00:00.000", '+', 366, 0);
	assertTimeDeltaToJulianDelta("P1Y", "-1-12-31T00:00:00.000", '+', 366, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "0-01-01T00:00:00.000", '+', 366, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "0-02-28T00:00:00.000", '+', 366, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "0-02-29T00:00:00.000", '+', 366, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "0-03-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "0-12-31T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "1-01-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "1-01-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "1-02-28T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "1-03-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "1-12-31T00:00:00.000", '+', 365, 0);

	assertTimeDeltaToJulianDelta("-P1Y", "-1-01-01T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "-1-01-01T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "-1-02-28T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "-1-03-01T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "-1-12-31T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "0-01-01T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "0-02-28T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "0-02-29T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "0-03-01T00:00:00.000", '-', -366, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "0-12-31T00:00:00.000", '-', -366, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "1-01-01T00:00:00.000", '-', -366, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "1-01-01T00:00:00.000", '-', -366, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "1-02-28T00:00:00.000", '-', -366, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "1-03-01T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y",  "1-12-31T00:00:00.000", '-', -365, 0);

	assertTimeDeltaToJulianDelta("P2Y",  "0-01-01T00:00:00.000", '+', 731, 0);
	assertTimeDeltaToJulianDelta("P6Y",  "0-01-01T00:00:00.000", '+', 2192, 0);
	assertTimeDeltaToJulianDelta("P26Y",  "0-01-01T00:00:00.000", '+', 9497, 0);
	assertTimeDeltaToJulianDelta("P126Y",  "0-01-01T00:00:00.000", '+', 46021, 0);
	assertTimeDeltaToJulianDelta("P406Y",  "0-01-01T00:00:00.000", '+', 148289, 0);
	assertTimeDeltaToJulianDelta("P10406Y",  "0-01-01T00:00:00.000", '+', 3800714, 0);
	assertTimeDeltaToJulianDelta("P10407Y",  "0-01-01T00:00:00.000", '+', 3801079, 0);
	assertTimeDeltaToJulianDelta("P10408Y",  "0-01-01T00:00:00.000", '+', 3801444, 0);
	assertTimeDeltaToJulianDelta("P10409Y",  "0-01-01T00:00:00.000", '+', 3801810, 0);
	assertTimeDeltaToJulianDelta("P10409Y",  "0-03-01T00:00:00.000", '+', 3801809, 0);
	assertTimeDeltaToJulianDelta("P10409Y",  "1-01-01T00:00:00.000", '+', 3801809, 0);
	assertTimeDeltaToJulianDelta("-P10409Y",  "1-01-01T00:00:00.000", '-', -3801810, 0);
	assertTimeDeltaToJulianDelta("-P10409Y",  "1201-01-01T00:00:00.000", '-', -3801810, 0);
	assertTimeDeltaToJulianDelta("-P10409Y",  "-1199-01-01T00:00:00.000", '-', -3801810, 0);
	assertTimeDeltaToJulianDelta("-P10409Y",  "0-01-01T00:00:00.000", '-', -3801809, 0);
	assertTimeDeltaToJulianDeltaWithoutExpectations("-P2147483648Y",  "0-01-01T00:00:00.000", '-');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P2147483647Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("-P2147483648Y",  "2147483647-01-01T00:00:00.000", '-');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P2147483647Y",  "-2147483648-01-01T00:00:00.000", '+');

	assertTimeDeltaToJulianDelta("-P1Y",  "2147483647-01-01T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("P1Y",  "-2147483648-01-01T00:00:00.000", '+', 366, 0);
}
END_TEST

START_TEST(test_timeDeltaToJulianDelta_ProlepticGregorian_2)
{
	assertTimeDeltaToJulianDelta("P01M", "1-01-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-02-01T00:00:00.000", '+', 28, 0);
	assertTimeDeltaToJulianDelta("P01M", "0-02-01T00:00:00.000", '+', 29, 0);
	assertTimeDeltaToJulianDelta("P01M", "0-03-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-04-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-05-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-06-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-07-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-08-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-09-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-10-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-11-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-12-01T00:00:00.000", '+', 31, 0);

	assertTimeDeltaToJulianDelta("-P01M", "1-01-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-02-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-03-01T00:00:00.000", '-', -28, 0);
	assertTimeDeltaToJulianDelta("-P01M", "0-03-01T00:00:00.000", '-', -29, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-04-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-05-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-06-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-07-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-08-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-09-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-10-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-11-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-12-01T00:00:00.000", '-', -30, 0);
}
END_TEST

START_TEST(test_timeDeltaToJulianDelta_ProlepticGregorian_3)
{
	assertTimeDeltaToJulianDelta("P01D", "1-01-01T00:00:00.000", '+', 1, 0);
	assertTimeDeltaToJulianDelta("PT01H", "1-01-01T00:00:00.000", '+', 0, 3600000);
	assertTimeDeltaToJulianDelta("PT01M", "1-01-01T00:00:00.000", '+', 0, 60000);
	assertTimeDeltaToJulianDelta("PT01S", "1-01-01T00:00:00.000", '+', 0, 1000);
	assertTimeDeltaToJulianDelta("PT00.001S", "1-01-01T00:00:00.000", '+', 0, 1);

	assertTimeDeltaToJulianDelta("-P01D", "1-01-01T00:00:00.000", '-', -1, 0);
	assertTimeDeltaToJulianDelta("-PT01H", "1-01-01T00:00:00.000", '-', 0, -3600000);
	assertTimeDeltaToJulianDelta("-PT01M", "1-01-01T00:00:00.000", '-', 0, -60000);
	assertTimeDeltaToJulianDelta("-PT01S", "1-01-01T00:00:00.000", '-', 0, -1000);
	assertTimeDeltaToJulianDelta("-PT00.001S", "1-01-01T00:00:00.000", '-', 0, -1);

	assertTimeDeltaToJulianDelta("P25Y", "1-01-01T00:00:00.000", '+', 9131, 0);
	assertTimeDeltaToJulianDelta("P25Y4M", "1-01-01T00:00:00.000", '+', 9251, 0);
	assertTimeDeltaToJulianDelta("P25Y4M12D", "1-01-01T00:00:00.000", '+', 9263, 0);
	assertTimeDeltaToJulianDelta("P25Y4M12DT06H23M", "1-01-01T00:00:00.000", '+', 9263, 22980000);
	assertTimeDeltaToJulianDelta("P25Y4M12DT06H23M56S", "1-01-01T00:00:00.000", '+', 9263, 23036000);
	assertTimeDeltaToJulianDelta("P25Y4M12DT06H23M56.132S", "0-01-01T00:00:00.000", '+', 9264, 23036132);
	assertTimeDeltaToJulianDelta("P25Y4M12DT06H23M56.132S", "0-03-01T00:00:00.000", '+', 9265, 23036132);
	assertTimeDeltaToJulianDelta("-P25Y4M12DT06H23M56.132S", "0-03-01T00:00:00.000", '-', -9264, -23036132);

	assertTimeDeltaToJulianDelta("P366D", "0-01-01T00:00:00.000", '+', 366, 0);
	assertTimeDeltaToJulianDelta("P366D", "0-03-01T00:00:00.000", '+', 366, 0);
	assertTimeDeltaToJulianDelta("P366D", "1-01-01T00:00:00.000", '+', 366, 0);

	assertTimeDeltaToJulianDelta("-P366D", "0-01-01T00:00:00.000", '-', -366, 0);
	assertTimeDeltaToJulianDelta("-P366D", "0-03-01T00:00:00.000", '-', -366, 0);
	assertTimeDeltaToJulianDelta("-P366D", "1-01-01T00:00:00.000", '-', -366, 0);

	assertTimeDeltaToJulianDelta("PT24H", "1-01-01T00:00:00.000", '+', 1, 0);
	assertTimeDeltaToJulianDelta("PT30H", "1-01-01T00:00:00.000", '+', 1, 21600000);
	assertTimeDeltaToJulianDelta("PT1440M", "1-01-01T00:00:00.000", '+', 1, 0);
	assertTimeDeltaToJulianDelta("PT4335M", "1-01-01T00:00:00.000", '+', 3, 900000);
	assertTimeDeltaToJulianDelta("PT86400S", "1-01-01T00:00:00.000", '+', 1, 0);
	assertTimeDeltaToJulianDelta("PT3653389.498S", "1-01-01T00:00:00.000", '+', 42, 24589498);

	assertTimeDeltaToJulianDelta("-PT24H", "1-01-01T00:00:00.000", '-', -1, 0);
	assertTimeDeltaToJulianDelta("-PT30H", "1-01-01T00:00:00.000", '-', -1, -21600000);
	assertTimeDeltaToJulianDelta("-PT1440M", "1-01-01T00:00:00.000", '-', -1, 0);
	assertTimeDeltaToJulianDelta("-PT4335M", "1-01-01T00:00:00.000", '-', -3, -900000);
	assertTimeDeltaToJulianDelta("-PT86400S", "1-01-01T00:00:00.000", '-', -1, 0);
	assertTimeDeltaToJulianDelta("-PT3653389.498S", "1-01-01T00:00:00.000", '-', -42, -24589498);
}
END_TEST

START_TEST(test_timeDeltaToJulianDelta_YearOf365Days)
{
	assertTimeDeltaToJulianDelta("P0Y", "0-01-01T00:00:00.000", '+', 0, 0);
	assertTimeDeltaToJulianDelta("P1Y", "0-01-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("P1Y", "1-01-01T00:00:00.000", '+', 365, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "0-01-01T00:00:00.000", '-', -365, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "1-01-01T00:00:00.000", '-', -365, 0);

	assertTimeDeltaToJulianDelta("P01M", "1-01-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-02-01T00:00:00.000", '+', 28, 0);
	assertTimeDeltaToJulianDelta("P01M", "0-02-01T00:00:00.000", '+', 28, 0);
	assertTimeDeltaToJulianDelta("P01M", "0-03-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-04-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-05-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-06-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-07-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-08-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-09-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-10-01T00:00:00.000", '+', 31, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-11-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-12-01T00:00:00.000", '+', 31, 0);

	assertTimeDeltaToJulianDelta("-P01M", "1-01-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-02-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-03-01T00:00:00.000", '-', -28, 0);
	assertTimeDeltaToJulianDelta("-P01M", "0-03-01T00:00:00.000", '-', -28, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-04-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-05-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-06-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-07-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-08-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-09-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-10-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-11-01T00:00:00.000", '-', -31, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-12-01T00:00:00.000", '-', -30, 0);

}
END_TEST

START_TEST(test_timeDeltaToJulianDelta_YearOf360Days)
{
	assertTimeDeltaToJulianDelta("P0Y", "0-01-01T00:00:00.000", '+', 0, 0);
	assertTimeDeltaToJulianDelta("P1Y", "0-01-01T00:00:00.000", '+', 360, 0);
	assertTimeDeltaToJulianDelta("P1Y", "1-01-01T00:00:00.000", '+', 360, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "0-01-01T00:00:00.000", '-', -360, 0);
	assertTimeDeltaToJulianDelta("-P1Y", "1-01-01T00:00:00.000", '-', -360, 0);

	assertTimeDeltaToJulianDelta("P01M", "1-01-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-02-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "0-02-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "0-03-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-04-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-05-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-06-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-07-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-08-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-09-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-10-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-11-01T00:00:00.000", '+', 30, 0);
	assertTimeDeltaToJulianDelta("P01M", "1-12-01T00:00:00.000", '+', 30, 0);

	assertTimeDeltaToJulianDelta("-P01M", "1-01-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-02-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-03-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "0-03-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-04-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-05-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-06-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-07-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-08-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-09-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-10-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-11-01T00:00:00.000", '-', -30, 0);
	assertTimeDeltaToJulianDelta("-P01M", "1-12-01T00:00:00.000", '-', -30, 0);
}
END_TEST

START_TEST(test_timeDeltaToJulianDelta_highDeltas)
{
	assertTimeDeltaToJulianDeltaWithoutExpectations("P0Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P1Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P9Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P99Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P999Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P9999Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P99999Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P999999Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P9999999Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P99999999Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P999999999Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("P2147483647Y",  "0-01-01T00:00:00.000", '+');
	assertTimeDeltaToJulianDeltaWithoutExpectations("-P2147483647Y",  "0-01-01T00:00:00.000", '-');
	assertTimeDeltaToJulianDeltaWithoutExpectations("-P2147483648Y",  "0-01-01T00:00:00.000", '-');
}
END_TEST

/**
 * check the string representation of mtime variables of type
 * datetime and timedelta
 */
START_TEST(test_timeDeltaAddLeadingZero)
{
  char tds[MAX_TIMEDELTA_STR_LEN];

  /* string representation od time delta adds '0' */
  TimeDelta oneDay= newTimeDelta("P1D");
  ck_assert_str_eq("P01D",timedeltaToString(oneDay,tds));

  TimeDelta oneMonth = newTimeDelta("PT1M");
  ck_assert_str_eq("PT01M",timedeltaToString(oneMonth,tds));

  TimeDelta what = newTimeDelta("P1Y6M");
  ck_assert_str_eq("P1Y06M",timedeltaToString(what,tds));

}
END_TEST

/**
 * check the handling of milliseconds
 */
START_TEST(test_timeDeltaMilliseconds)
{
  char tds[MAX_TIMEDELTA_STR_LEN];

  /* check milliseconds made with Raw-Constructor */
  DateTime rawDateTime = newRawDateTime(1999, 12, 12, 12, 12,12 , 1);
  ck_assert(12  == rawDateTime->time.minute);
  ck_assert(12  == rawDateTime->time.hour);
  ck_assert(1.0 == rawDateTime->time.ms);

  /* check milliseconds make with Sting-Constructor */
  DateTime tenth = newDateTime("1888-12-12T22:58:44.1");
  datetimeToString(tenth,tds);
  ck_assert_str_eq("1888-12-12T22:58:44.100",tds);
  ck_assert(100.0 == tenth->time.ms);

  /*  check with ms value lower than 100 */
  DateTime hundred = newDateTime("1888-12-12T22:58:44.01");
  ck_assert(10.0 == hundred->time.ms);

  /*  check with ms value lower than 10 */
  DateTime milli;
  milli = newDateTime("1888-12-12T22:58:44.001");
  ck_assert(1.0 == milli->time.ms);

  milli = newDateTime("1888-12-12T22:58:44.004");
  ck_assert(4.0 == milli->time.ms);
}
END_TEST

/**
 * test the addition of a month uppon the and of january
 *
 * the goal is to mimic a list of dates, which are at the end of a month
 */
START_TEST(test_mtime_add_months_atEnd)
{
  int i;
  struct _timedelta *td, *offset;
  struct _datetime *start, *ref;
#ifdef DEBUG
  char curs[MAX_DATETIME_STR_LEN];
#endif

  start  = newDateTime("2001-01-01T00:00:00");
  ref    = newDateTime("1999-12-31T23:50:00");
  td     = newTimeDelta("P1M");
  offset = newTimeDelta("-PT10M");

  ref    = start;

  struct _datetime* cur = constructAndCopyDateTime(ref);
#ifdef DEBUG
  printf("\n#==================================================================\n");
  printf("Results from test: test_mtime_add_months_atEnd\n");
  printf("Start date is: %s\n",datetimeToString(cur, curs));
#endif

  for (i = 0; i < 14; i++)
  {
    /* compute month shift */
    addTimeDeltaToDateTime(cur,td,cur);
    /* subtract 10 minutes */
    struct _datetime* cur_backshift = constructAndCopyDateTime(cur);
    addTimeDeltaToDateTime(cur_backshift,offset,cur_backshift);

#ifdef EBUG
    printf("             Current date is: %s (i=%d)\n",datetimeToString(cur_backshift, curs),i);
#endif

    switch (i%12) {
      case 0: case 2: case 4: case 6: case 7: case 9: case 11:
        ck_assert_int_eq(31, cur_backshift->date.day);
        break;
      case 1:
        ck_assert_int_eq(28, cur_backshift->date.day);
        break;
      case 3: case 5: case 8: case 10:
        ck_assert_int_eq(30, cur_backshift->date.day);
        break;
      default:
        printf(" error - this should never happen\n");
        break;
    }
  }
}
END_TEST

/**
 * test the addition of a month uppon the and of january: the defaut behaviour
 * i.e. february is skipped
 */
START_TEST(test_mtime_add_months)
{

  int i;
  struct _timedelta *td;
  struct _datetime *ref;
  char curs[MAX_DATETIME_STR_LEN];

  ref    = newDateTime("1999-12-31T23:50:00");
  td     = newTimeDelta("P1M");

  struct _datetime* cur = constructAndCopyDateTime(ref);

#ifdef EBUG
  printf("\n#==================================================================\n");
  printf("Results from test: test_mtime_add_months\n");
  printf("Start date is: %s\n",datetimeToString(cur, curs));
#endif
  for (i = 0; i < 14; i++)
  {
    /* compute month shift */
    addTimeDeltaToDateTime(cur,td,cur);

#ifdef EBUG
    printf("             Current date is: %s (i=%d)\n",datetimeToString(cur, curs),i);
#endif
    if (i == 0 ) ck_assert_str_eq("2000-01-31T23:50:00.000",datetimeToString(cur,curs));
    if (i == 1 ) ck_assert_str_eq("2000-03-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 2 ) ck_assert_str_eq("2000-04-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 3 ) ck_assert_str_eq("2000-05-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 4 ) ck_assert_str_eq("2000-06-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 5 ) ck_assert_str_eq("2000-07-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 6 ) ck_assert_str_eq("2000-08-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 7 ) ck_assert_str_eq("2000-09-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 8 ) ck_assert_str_eq("2000-10-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 9 ) ck_assert_str_eq("2000-11-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 10) ck_assert_str_eq("2000-12-02T23:50:00.000",datetimeToString(cur,curs));
    if (i == 11) ck_assert_str_eq("2001-01-02T23:50:00.000",datetimeToString(cur,curs));
  }
}
END_TEST

/**
 * check adding negative time periods
 */
START_TEST(test_add_negative_delta)
{
  struct _datetime *dummy;
  struct _datetime *init = newDateTime("1111-11-11T00:00:00.000");
  struct _timedelta *delta = newTimeDelta("PT10M");
  struct _timedelta *minusdelta = newTimeDelta("-PT10M");
  char init_s[MAX_DATETIME_STR_LEN];

  ck_assert_str_eq("1111-11-11T00:10:00.000",datetimeToString(addTimeDeltaToDateTime(init,delta,init),init_s));
  ck_assert_str_eq("1111-11-11T00:20:00.000",datetimeToString(addTimeDeltaToDateTime(init,delta,init),init_s));
  ck_assert_str_eq("1111-11-11T00:10:00.000",datetimeToString(addTimeDeltaToDateTime(init,minusdelta,init),init_s));
  ck_assert_str_eq("1111-11-11T00:00:00.000",datetimeToString(addTimeDeltaToDateTime(init,minusdelta,init),init_s));
}
END_TEST

/**
 * test events created with offset to trigger 10 minutes before end of months
 */
START_TEST(test_event_at_end_of_month)
{
  /* some lines stolen from model_integrate.c */
  char current_time[MAX_DATETIME_STR_LEN];  // Again, we use *STR_LEN.
  char current_step[MAX_TIMEDELTA_STR_LEN]; // and here too.
  struct _timedelta *timestep  = newTimeDelta("PT10M");
  struct _datetime *dummy;
  struct _datetime *start_date = newDateTime("2013-01-01T00:00:00.000");
  struct _datetime *stop_date  = newDateTime("2014-03-03T14:00:00.000");
  struct _datetime *model_time = newDateTime("2013-01-01T00:00:00.000");

  struct _event* out = newEvent("outputEvent",
                                "2013-01-01T00:00:00.000",
                                "2013-01-01T00:00:00.000",
                                "2019-01-01T00:00:00.000",
                                "P1M",
                                 NULL);
  struct _event* outAtEnd = newEvent("outputEventAtEnd",
                                     "2012-12-31T23:50:00.000",
                                     "2013-01-01T00:00:00.000",
                                     "2019-01-01T00:00:00.000",
                                     "P1M",
                                     NULL);

#ifdef EBUG
  printf ( "\n#==================================================================\n" );
  printf ( "Output of test: test_event_at_end_of_month\n" );
  printf("Model time step: %s\n", timedeltaToString(timestep, current_step));
  printf("Model start time: %s\n\n", datetimeToString(model_time, current_time));
#endif

  do
    {
      /* check if outputEvent is active */
      dummy = addTimeDeltaToDateTime(model_time,newTimeDelta("PT10"),dummy);
      bool activeOut      = isCurrentEventActive(out,model_time, NULL, NULL);
      bool activeOutShift = isCurrentEventActive(out,dummy , NULL, NULL);
      bool activeOutAtEnd = isCurrentEventActive(outAtEnd,model_time, NULL, NULL);
#ifdef EBUG
      if (activeOut || activeOutAtEnd || activeOutShift) printf("Model time: %s  ", datetimeToString(model_time, current_time));
      if (activeOut)                   printf("outputEvent active! ");
      if (activeOutAtEnd)              printf("outputEventAtEnd active!");
      if (activeOutShift)              printf("outputShifted    active!");
      if (activeOut || activeOutAtEnd ) printf("\n");
#endif

      addTimeDeltaToDateTime(model_time, timestep, model_time); /* Increment time by timestep. */


    } while (compareDatetime(model_time, stop_date) <= 0);
}
END_TEST

START_TEST(test_getTimeDeltaFromDateTime)
{
        // These tests are taken from ac66aeaef2dde828aa31b63f88076070a8282b49
        assertGetTimeDeltaFromDateTime("1979-03-01T00:00:00.000", "1979-01-01T01:00:00.000", "P01M27DT23H");
	assertGetTimeDeltaFromDateTime("1979-07-01T00:00:00.000", "1979-01-01T01:00:00.000", "P05M29DT23H");
	assertGetTimeDeltaFromDateTime("1979-12-01T00:00:00.000", "1979-01-01T01:00:00.000", "P10M29DT23H");
	assertGetTimeDeltaFromDateTime("1980-01-01T00:00:00.000", "1979-01-01T01:00:00.000", "P11M30DT23H");

	// These tests are taken from 5bcba6591c38fbb9f7a3a5f63c86963d8a066595
	assertGetTimeDeltaFromDateTime("2017-07-31T00:00:00.000", "2017-07-01T00:00:00.000", "P30D");
	assertGetTimeDeltaFromDateTime("2017-08-01T00:00:00.000", "2017-07-01T00:00:00.000", "P1M");

	// FIXME: This needs to be fulfilled for amip restart to work fine.
	assertGetTimeDeltaFromDateTime("1981-01-01T00:00:00.000", "1980-01-01T00:10:00.000", "P11M30DT23H50M");
}
END_TEST

START_TEST(test_timeDeltaToJulianDeltaToTimeDelta)
{
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-01-01T01:00:00.000", "P01M27DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-02-01T01:00:00.000", "P01M27DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-03-01T00:00:00.000", "P01M27DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-04-01T01:00:00.000", "P01M27DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-05-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1979-06-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1979-07-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1979-08-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1979-09-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1979-10-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1979-11-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1979-12-01T01:00:00.000", "P01M27DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1980-01-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1980-02-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1980-03-01T01:00:00.000", "P01M27DT23H");
        assertTimeDeltaToJulianDeltaToTimeDelta("1980-04-01T01:00:00.000", "P01M27DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-01-01T01:00:00.000", "P05M29DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-07-01T00:00:00.000", "P05M29DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-01-01T01:00:00.000", "P10M29DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-12-01T00:00:00.000", "P10M29DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1979-01-01T01:00:00.000", "P11M30DT23H");
	assertTimeDeltaToJulianDeltaToTimeDelta("1980-01-01T00:00:00.000", "P11M30DT23H");

	assertTimeDeltaToJulianDeltaToTimeDelta("2017-07-01T00:00:00.000", "P30D");
	assertTimeDeltaToJulianDeltaToTimeDelta("2017-07-31T00:00:00.000", "P30D");
	assertTimeDeltaToJulianDeltaToTimeDelta("2017-07-01T00:00:00.000", "P1M");
	assertTimeDeltaToJulianDeltaToTimeDelta("2017-08-01T00:00:00.000", "P1M");

	assertTimeDeltaToJulianDeltaToTimeDelta("1980-01-01T00:10:00.000", "P11M30DT23H50M");
	assertTimeDeltaToJulianDeltaToTimeDelta("1981-01-01T00:00:00.000", "P11M30DT23H50M");
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

void add_mtime_timedelta_test_to_suite(Suite* suite)
{
	TCase *tcase_ProlepticGregorian = tcase_create("mtime_time_test_ProlepticGregorian");
	suite_add_tcase(suite, tcase_ProlepticGregorian);
	tcase_add_checked_fixture(tcase_ProlepticGregorian, setup_ProlepticGregorian, teardown);
	tcase_add_test(tcase_ProlepticGregorian, test_create_timedelta_from_valid_strings);
	tcase_add_test(tcase_ProlepticGregorian, test_invalid_strings);
	tcase_add_test(tcase_ProlepticGregorian, test_constructAndCopyTimeDelta);
	tcase_add_test(tcase_ProlepticGregorian, test_replaceTimeDelta);
	tcase_add_test(tcase_ProlepticGregorian, test_timeDeltaToJulianDelta_ProlepticGregorian_1);
	tcase_add_test(tcase_ProlepticGregorian, test_timeDeltaToJulianDelta_ProlepticGregorian_2);
	tcase_add_test(tcase_ProlepticGregorian, test_timeDeltaToJulianDelta_ProlepticGregorian_3);
	tcase_add_test(tcase_ProlepticGregorian, test_timeDeltaToJulianDelta_highDeltas);
	tcase_add_test(tcase_ProlepticGregorian, test_getTimeDeltaFromDateTime);
        tcase_add_test(tcase_ProlepticGregorian, test_timeDeltaToJulianDeltaToTimeDelta);

	TCase *tcase_YearOf365Days = tcase_create("mtime_time_test_YearOf365Days");
	suite_add_tcase(suite, tcase_YearOf365Days);
	tcase_add_checked_fixture(tcase_YearOf365Days, setup_YearOf365Days, teardown);
	tcase_add_test(tcase_YearOf365Days, test_timeDeltaToJulianDelta_YearOf365Days);

	TCase *tcase_YearOf360Days = tcase_create("mtime_time_test_YearOf360Days");
	suite_add_tcase(suite, tcase_YearOf360Days);
	tcase_add_checked_fixture(tcase_YearOf360Days, setup_YearOf360Days, teardown);
	tcase_add_test(tcase_YearOf360Days, test_timeDeltaToJulianDelta_YearOf360Days);

        TCase *tcase_timedeltaCompute = tcase_create("mtime_timedelta_computations");
	tcase_add_checked_fixture(tcase_timedeltaCompute, setup_ProlepticGregorian, teardown);
        suite_add_tcase(suite,tcase_timedeltaCompute);
        tcase_add_test(tcase_timedeltaCompute,test_timeDeltaAddLeadingZero);
        tcase_add_test(tcase_timedeltaCompute,test_mtime_add_months);
        tcase_add_test(tcase_timedeltaCompute,test_mtime_add_months_atEnd);
        tcase_add_test(tcase_timedeltaCompute,test_add_negative_delta);
        tcase_add_test(tcase_timedeltaCompute,test_event_at_end_of_month);
        tcase_add_test(tcase_timedeltaCompute,test_timeDeltaMilliseconds);
}

/*** SPECIAL ASSERT FUNCTIONS ***/

void assertTimeDelta(const char* input_string, bool std_form, char sign, int64_t year, int month, int day, int hour, int minute, int second, int ms, const char* expected_output_string)
{
	const char* format = "Parsing of timedelta string \"%s\" failed.";
	size_t length = snprintf(NULL, 0, format, input_string) + 1;
	char* msg = malloc(length);
	snprintf(msg, length, format, input_string);

	char tmp[MAX_DATE_STR_LEN];

	TimeDelta timedelta = newTimeDelta(input_string);
	ck_assert_msg(timedelta != NULL, msg);
	ck_assert_str_eq(expected_output_string, timedeltaToString(timedelta, tmp));
	ck_assert_str_eq(expected_output_string, tmp);
	ck_assert_int_eq(std_form, timedelta->flag_std_form);
	ck_assert_int_eq(sign, timedelta->sign);
	ck_assert_int_eq(year, timedelta->year);
	ck_assert_int_eq(month, timedelta->month);
	ck_assert_int_eq(day, timedelta->day);
	ck_assert_int_eq(hour, timedelta->hour);
	ck_assert_int_eq(minute, timedelta->minute);
	ck_assert_int_eq(second, timedelta->second);
	ck_assert_int_eq(ms, timedelta->ms);

	TimeDelta raw_timedelta = newRawTimeDelta(sign, year, month, day, hour, minute, second, ms);
	ck_assert(raw_timedelta != NULL);
	ck_assert_str_eq(expected_output_string, timedeltaToString(raw_timedelta, tmp));
	ck_assert_int_eq(std_form, timedelta->flag_std_form);
	ck_assert_int_eq(sign, raw_timedelta->sign);
	ck_assert_int_eq(year, raw_timedelta->year);
	ck_assert_int_eq(month, raw_timedelta->month);
	ck_assert_int_eq(day, raw_timedelta->day);
	ck_assert_int_eq(hour, raw_timedelta->hour);
	ck_assert_int_eq(minute, raw_timedelta->minute);
	ck_assert_int_eq(second, raw_timedelta->second);
	ck_assert_int_eq(ms, raw_timedelta->ms);

	TimeDelta timedelta_from_string = newTimeDelta(timedeltaToString(raw_timedelta, tmp));
	ck_assert(timedelta_from_string != NULL);
	ck_assert_str_eq(expected_output_string, timedeltaToString(timedelta_from_string, tmp));

	deallocateTimeDelta(timedelta);
	free(msg);
}

void assertInvalidTimeDelta(const char* input_string)
{
	TimeDelta timedelta = newTimeDelta(input_string);

	if (timedelta != NULL)
	{
		const char* format = "Invalid timedelta string \"%s\" didn't result in null pointer but in timedelta \"%s\".";
		char output_string[MAX_DATE_STR_LEN];
		timedeltaToString(timedelta, output_string);

		size_t length = snprintf(NULL, 0, format, input_string, output_string) + 1;
		char* msg = malloc(length);
		snprintf(msg, length, format, input_string, output_string);
		deallocateTimeDelta(timedelta);

		ck_abort_msg(msg);
		free(msg); //Never executed?
	}
}

void assertTimeDeltaToJulianDelta(const char* timedelta_string, const char* base_datetime_string, char expected_jd_sign, int64_t expected_jd_day, int64_t expected_jd_ms)
{
	TimeDelta timedelta = newTimeDelta(timedelta_string);
	ck_assert(timedelta != NULL);

	DateTime datetime = newDateTime(base_datetime_string);
	ck_assert(datetime != NULL);

	JulianDelta juliandeltaTemp = newJulianDelta('+', 0, 0);
	JulianDelta juliandelta = timeDeltaToJulianDelta(timedelta, datetime, juliandeltaTemp);
	ck_assert(datetime != NULL);
	ck_assert(juliandeltaTemp == juliandelta);

	const char* format = "timeDeltaToJulianDelta failed for \"%s\" + \"%s\": expected:(%c, %lld, %lld) actual:(%c, %lld, %lld).";
	size_t length = snprintf(NULL, 0, format, base_datetime_string, timedelta_string, expected_jd_sign, expected_jd_day, expected_jd_ms, juliandelta->sign, juliandelta->day, juliandelta->ms) + 1;
	char* msg = malloc(length);
	snprintf(msg, length, format, base_datetime_string, timedelta_string, expected_jd_sign, expected_jd_day, expected_jd_ms, juliandelta->sign, juliandelta->day, juliandelta->ms);

	ck_assert_msg(expected_jd_sign == juliandelta->sign, msg);
	ck_assert_msg(expected_jd_day == juliandelta->day, msg);
	ck_assert_msg(expected_jd_ms == juliandelta->ms, msg);

	deallocateTimeDelta(timedelta);
	deallocateDateTime(datetime);
	deallocateJulianDelta(juliandelta);
}

void assertTimeDeltaToJulianDeltaWithoutExpectations(const char* timedelta_string, const char* base_datetime_string, char expected_jd_sign)
{
	TimeDelta timedelta = newTimeDelta(timedelta_string);
	ck_assert(timedelta != NULL);

	DateTime datetime = newDateTime(base_datetime_string);
	ck_assert(datetime != NULL);

	JulianDelta juliandeltaTemp = newJulianDelta('+', 0, 0);
	JulianDelta juliandelta = timeDeltaToJulianDelta(timedelta, datetime, juliandeltaTemp);
	ck_assert(datetime != NULL);
	ck_assert(juliandeltaTemp == juliandelta);

	ck_assert(expected_jd_sign == juliandelta->sign);

	deallocateTimeDelta(timedelta);
	deallocateDateTime(datetime);
	deallocateJulianDelta(juliandelta);
}


// It seems we need this tested too.
void assertGetTimeDeltaFromDateTime (const char* dt1_string, const char* dt2_string, const char* expected_td_string)
{
	struct _datetime* dt1 = newDateTime(dt1_string);
	ck_assert(dt1 != NULL);

	struct _datetime* dt2 = newDateTime(dt2_string);
	ck_assert(dt2 != NULL);

	struct _timedelta td;
        void* err = getTimeDeltaFromDateTime(dt1, dt2, &td);
	ck_assert(err != NULL);

	struct _timedelta* expected_td = newTimeDelta(expected_td_string);
	ck_assert(expected_td != NULL);

	char str[MAX_TIMEDELTA_STR_LEN];
	timedeltaToString(&td, str);
	ck_assert_msg(compareTimeDelta(&td, expected_td) == equal_to, "base datetime %s expected %s != %s returned", dt2_string, expected_td_string, str);

	deallocateDateTime(dt1);
	deallocateDateTime(dt2);
	deallocateTimeDelta(expected_td);
}

void assertTimeDeltaToJulianDeltaToTimeDelta (const char* base_dt_string, const char* td_string)
{
	void* err;

	struct _datetime* base_dt = newDateTime(base_dt_string);
	ck_assert(base_dt != NULL);

	struct _timedelta* td = newTimeDelta(td_string);
	ck_assert(td != NULL);

	struct _juliandelta jd;
	err = timeDeltaToJulianDelta(td, base_dt, &jd);
	ck_assert(err != NULL);

	struct _timedelta res_td;
	err = julianDeltaToTimeDelta(&jd, base_dt, &res_td);
	ck_assert(err != NULL);

	char str[MAX_TIMEDELTA_STR_LEN];
	timedeltaToString(&res_td, str);

	ck_assert_msg(compareTimeDelta(td, &res_td) == equal_to, "base datetime %s original %s != %s returned",
			base_dt_string, td_string, str);

	//FIXME: Deallocate.
}

// FIXME: Fill with content
void assertJulianDeltaToTimeDelta ()
{
	
}

