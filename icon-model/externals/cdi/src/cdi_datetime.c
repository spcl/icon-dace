/* DO NOT REMOVE the config.h include file under any circumstances,
 * it's very much needed on some platforms */
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
/* DO NOT REMOVE the above config.h include file under any
 * circumstances as long as it's the autoconf configuration header
 * used to build this package. When it's missing on some platforms,
 * some poor person has to do long, tedious debugging sessions, where
 * struct offsets almost imperceptibly change from one file to the
 * next to find out what happened */

#include "cdi_datetime.h"
#include <stdio.h>
#include <stdlib.h>

// ==================================================================
#ifdef __cplusplus
extern "C"
{
#endif

  // clang-format off
void
cdiDecodeDate(int date, int *year, int *month, int *day)
{
  int iyear = date / 10000;
  *year = iyear;
  int idate = date - iyear * 10000;
  if (idate < 0) idate = -idate;
  int imonth = idate / 100;
  *month = imonth;
  *day = idate - imonth * 100;
}

int
cdiEncodeDate(int year, int month, int day)
{
  int iyear = abs(year);
  int date = iyear * 10000 + month * 100 + day;
  if (year < 0) date = -date;

  return date;
}

void
cdiDecodeTime(int time, int *hour, int *minute, int *second)
{
  int ihour = time / 10000, itime = time - ihour * 10000, iminute = itime / 100;
  *hour = ihour;
  *minute = iminute;
  *second = itime - iminute * 100;
}

int
cdiEncodeTime(int hour, int minute, int second)
{
  return hour * 10000 + minute * 100 + second;
}
  // clang-format on

#ifdef __cplusplus
}
#endif
// ==================================================================

CdiDate
cdiDate_set(int64_t date)
{
  int64_t iyear = date / 10000;
  int year = (int) iyear;
  int64_t idate = date - iyear * 10000;
  if (idate < 0) idate = -idate;
  int64_t imonth = idate / 100;
  int month = (int) imonth;
  int day = (int) (idate - imonth * 100);

  CdiDate cdiDate;
  cdiDate.year = year;
  cdiDate.month = (short) month;
  cdiDate.day = (short) day;

  return cdiDate;
}

CdiTime
cdiTime_set(int time)
{
  int hour, minute, second, ms = 0;
  cdiDecodeTime(time, &hour, &minute, &second);

  CdiTime cdiTime;
  cdiTime.hour = (short) hour;
  cdiTime.minute = (short) minute;
  cdiTime.second = (short) second;
  cdiTime.ms = (short) ms;

  return cdiTime;
}

CdiDateTime
cdiDateTime_set(int64_t date, int time)
{
  CdiDateTime cdiDateTime;
  cdiDateTime.date = cdiDate_set(date);
  cdiDateTime.time = cdiTime_set(time);

  return cdiDateTime;
}

int64_t
cdiDate_get(CdiDate cdiDate)
{
  int64_t iyear = abs(cdiDate.year);
  int64_t date = iyear * 10000 + cdiDate.month * 100 + cdiDate.day;
  if (cdiDate.year < 0) date = -date;

  return date;
}

int
cdiTime_get(CdiTime cdiTime)
{
  return cdiEncodeTime(cdiTime.hour, cdiTime.minute, cdiTime.second);
}

CdiDate
cdiDate_encode(int year, int month, int day)
{
  CdiDate cdiDate;
  cdiDate.year = year;
  cdiDate.month = (short) month;
  cdiDate.day = (short) day;

  return cdiDate;
}

void
cdiDate_decode(CdiDate cdiDate, int *year, int *month, int *day)
{
  *year = cdiDate.year;
  *month = cdiDate.month;
  *day = cdiDate.day;
}

CdiTime
cdiTime_encode(int hour, int minute, int second, int ms)
{
  CdiTime cdiTime;
  cdiTime.hour = (short) hour;
  cdiTime.minute = (short) minute;
  cdiTime.second = (short) second;
  cdiTime.ms = (short) ms;

  return cdiTime;
}

void
cdiTime_decode(CdiTime cdiTime, int *hour, int *minute, int *second, int *ms)
{
  *hour = cdiTime.hour;
  *minute = cdiTime.minute;
  *second = cdiTime.second;
  *ms = cdiTime.ms;
}

void
cdiDate_init(CdiDate *cdiDate)
{
  cdiDate->year = 0;
  cdiDate->month = 0;
  cdiDate->day = 0;
}

void
cdiTime_init(CdiTime *cdiTime)
{
  cdiTime->hour = 0;
  cdiTime->minute = 0;
  cdiTime->second = 0;
  cdiTime->ms = 0;
}

void
cdiDateTime_init(CdiDateTime *cdiDateTime)
{
  cdiDate_init(&cdiDateTime->date);
  cdiTime_init(&cdiDateTime->time);
}

bool
cdiDate_isEQ(CdiDate cdiDate1, CdiDate cdiDate2)
{
  // clang-format off
  return (cdiDate1.year  == cdiDate2.year
       && cdiDate1.month == cdiDate2.month
       && cdiDate1.day   == cdiDate2.day);
  // clang-format on
}

bool
cdiTime_isEQ(CdiTime cdiTime1, CdiTime cdiTime2)
{
  // clang-format off
  return (cdiTime1.hour   == cdiTime2.hour
       && cdiTime1.minute == cdiTime2.minute
       && cdiTime1.second == cdiTime2.second
       && cdiTime1.ms     == cdiTime2.ms);
  // clang-format on
}

bool
cdiDateTime_isEQ(CdiDateTime cdiDateTime1, CdiDateTime cdiDateTime2)
{
  // clang-format off
  return (cdiDateTime1.date.year   == cdiDateTime2.date.year
       && cdiDateTime1.date.month  == cdiDateTime2.date.month
       && cdiDateTime1.date.day    == cdiDateTime2.date.day
       && cdiDateTime1.time.hour   == cdiDateTime2.time.hour
       && cdiDateTime1.time.minute == cdiDateTime2.time.minute
       && cdiDateTime1.time.second == cdiDateTime2.time.second
       && cdiDateTime1.time.ms     == cdiDateTime2.time.ms);
  // clang-format on
}

bool
cdiDateTime_isNE(CdiDateTime cdiDateTime1, CdiDateTime cdiDateTime2)
{
  return !cdiDateTime_isEQ(cdiDateTime1, cdiDateTime2);
}

bool
cdiDateTime_isLT(CdiDateTime cdiDateTime1, CdiDateTime cdiDateTime2)
{
  int64_t date1 = cdiDate_get(cdiDateTime1.date);
  int64_t date2 = cdiDate_get(cdiDateTime2.date);
  int time1 = cdiTime_get(cdiDateTime1.time);
  int time2 = cdiTime_get(cdiDateTime2.time);
  return (date1 < date2 || (date1 == date2 && time1 < time2));
}

bool
cdiDateTime_isNull(CdiDateTime cdiDateTime)
{
  // clang-format off
  return (cdiDateTime.date.year == 0
       && cdiDateTime.date.month == 0
       && cdiDateTime.date.day == 0
       && cdiDateTime.time.hour == 0
       && cdiDateTime.time.minute == 0
       && cdiDateTime.time.second == 0
       && cdiDateTime.time.ms == 0);
  // clang-format on
}

#define DATE_FORMAT "%5.4d-%2.2d-%2.2d"
#define TIME_FORMAT "%2.2d:%2.2d:%2.2d"

const char *
CdiDateTime_string(CdiDateTime cdiDateTime)
{
  int year, month, day;
  cdiDate_decode(cdiDateTime.date, &year, &month, &day);
  int hour, minute, second, ms;
  cdiTime_decode(cdiDateTime.time, &hour, &minute, &second, &ms);

  static char datetimeString[64];
  snprintf(datetimeString, sizeof(datetimeString), DATE_FORMAT "T" TIME_FORMAT, year, month, day, hour, minute, second);

  return datetimeString;
}
