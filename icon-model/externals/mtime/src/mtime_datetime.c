// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @file mtime_datetime.c
 *
 * @brief DateTime and some operations supported on DateTime.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

#include "mtime_datetime.h"

#include "mtime_calendar.h"
#include "mtime_julianDay.h"
#include "mtime_iso8601.h"

/* MIN allowed year : -2147483648
   MAX allowed year :  2147483647
   Upto MilliSecond resolution supported. */

/**
 * @brief Construct new DateTime using an ISO 8601 conforming string.
 *
 *        MAX allowed year = 2147483647 and MIN allowed year = -2147483648.
 *
 *         Strings of the following form are valid:
 *
 *         <table>
 *         <tr><th> Extended Format: <th>
 *         <tr><td> +/-([Y]+)-MM-DDTHH:MM:SS.SSS  <td> [Y]+ stands for year can be of any char length (>= 1)
 *                                                     but must be in range: [MIN allowed year, MAX allowed year].
 *                                                     Sign can be +/- or empty with default +.
 *                                                     Also, .SSS can be of length 1-3
 *         <tr><td> +/-([Y]+)-MM-DDTHH:MM:SS      <td>
 *         <tr><td> +/-([Y]+)-MM-DDTHH:MM         <td>
 *         <tr><td> +/-([Y]+)-MM-DDTHH            <td>
 *         <tr><td> +/-([Y]+)-MM-DD               <td>
 *         <tr><td> +/-([Y]+)-MM                  <td>
 *         </table>
 *
 *         <table>
 *         <tr><th> Basic Format: <th>
 *         <tr><td> +/-YYYYMMDDTHHMMSS.SSS        <td> Note that Year is 4 char long.
 *                                                     Sign can be +/- or empty with default +.
 *                                                     Also, .SSS can be of length 1-3
 *         <tr><td> +/-YYYYMMDDTHHMMSS            <td>
 *         <tr><td> +/-YYYYMMDDTHHMM              <td>
 *         <tr><td> +/-YYYYMMDDTHH                <td>
 *         <tr><td> +/-YYYYMMDD                   <td>
 *         </table>
 *
 *         In essence, strings follow the guidelines in ISO8601 standard (but do not support all formats).
 *
 * @param  dts
 *         A pointer to char. The string should contain parameters with which DateTime is created.
 *         A string literal can be accepted.
 *
 * @return dt
 *         A pointer to a filled datetime.
 */

struct _datetime *
newDateTime(const char *dts)
{
  struct _datetime *dt = NULL;
  if ((dts != NULL) && (getCalendarType()))
    {
      /* Intialize a dummy ISO8601 object. */
      struct iso8601_datetime isoDt = { .sign_of_year = '+' };

      /* Check ISO8601 compliance. */
      if (verify_string_datetime(dts, &isoDt) == DATETIME_MATCH
          && (dt = (struct _datetime *) calloc(1, sizeof(struct _datetime))) != NULL)
        {
          dt->date.year = isoDt.year;
          dt->date.month = isoDt.month;
          dt->date.day = isoDt.day;
          dt->time.hour = isoDt.hour;
          dt->time.minute = isoDt.minute;
          dt->time.second = isoDt.second;
          dt->time.ms = isoDt.ms;
        }
    }
  return dt;
}

/**
 * @brief Construct new DateTime using 'raw' numerical values.
 *
 * @param  _year
 *         An "int64_t" value denoting the year part of DateTime.
 * @param  _month
 *         An "int" value denoting the month part of DateTime.
 * @param  _day
 *         An "int" value denoting the day part of DateTime.
 * @param  _hour
 *         An "int" value denoting the hour part of DateTime.
 * @param  _minute
 *         An "int" value denoting the minute part of DateTime.
 * @param  _second
 *         An "int" value denoting the second part of DateTime.
 * @param  _ms
 *         An "int" value denoting the milli-second part of DateTime.
 *
 * @return dt
 *         A pointer to a filled DateTime.
 */

struct _datetime *
newRawDateTime(int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms)
{
  char dts[MAX_DATETIME_STR_LEN];
  memset(dts, 0, sizeof(dts));

  snprintf(dts, MAX_DATETIME_STR_LEN, "%" PRIi64 "-%02d-%02dT%02d:%02d:%02d.%03d", _year, _month, _day, _hour, _minute, _second,
           _ms);

  /* Resuse string interface and create object. */
  struct _datetime *dt = newDateTime(dts);

  return dt;
}

/**
 * @brief Copy the values and construct a new datetime.
 *
 * @param  dt
 *         A pointer to struct _datetime. Values of dt are used to initialize the new datetime being created.
 *
 * @return _dt
 *         A pointer to an initialized DateTime.
 */

struct _datetime *
constructAndCopyDateTime(struct _datetime *dt)
{
  if (dt != NULL)
    return newRawDateTime(dt->date.year, dt->date.month, dt->date.day, dt->time.hour, dt->time.minute, dt->time.second,
                          dt->time.ms);
  else
    return NULL;
}

/**
 * @brief Destructor of DateTime.
 *
 * @param  dt
 *         A pointer to struct _datetime. dt is deallocated.
 */

void
deallocateDateTime(struct _datetime *dt)
{
  if (dt != NULL)
    {
      free(dt);
      dt = NULL;
    }
}

/**
 * @brief Compare two datetimes and return (dt1 > dt2) OR (dt1 = dt2) OR (dt1 < dt2).
 *
 * @param  dt1
 *         A pointer to struct _datetime.
 *
 * @param  dt2
 *         A pointer to struct _datetime.
 *
 * @return boolean
 *         if (dt1 > dt2), return greater_than. If (dt1 == dt2), return equal_to. If (dt1 < dt2), return less_than.
 *         Return compare_error indicating error.
 */

compare_return_val
compareDatetime(struct _datetime *dt1, struct _datetime *dt2)
{
  if ((dt1 != NULL) && (dt2 != NULL))
    {

      /* Initialize. */
      compare_return_val boolean = compare_error;

      if (dt1->date.year == dt2->date.year)
        {
          if (dt1->date.month == dt2->date.month)
            {
              if (dt1->date.day == dt2->date.day)
                {
                  if (dt1->time.hour == dt2->time.hour)
                    {
                      if (dt1->time.minute == dt2->time.minute)
                        {
                          if (dt1->time.second == dt2->time.second)
                            {
                              if (dt1->time.ms == dt2->time.ms)
                                {
                                  boolean = equal_to;
                                }
                              else if (dt1->time.ms > dt2->time.ms)
                                {
                                  boolean = greater_than;
                                }
                              else
                                {
                                  boolean = less_than;
                                }
                            }
                          else if (dt1->time.second > dt2->time.second)
                            {
                              boolean = greater_than;
                            }
                          else
                            {
                              boolean = less_than;
                            }
                        }
                      else if (dt1->time.minute > dt2->time.minute)
                        {
                          boolean = greater_than;
                        }
                      else
                        {
                          boolean = less_than;
                        }
                    }
                  else if (dt1->time.hour > dt2->time.hour)
                    {
                      boolean = greater_than;
                    }
                  else
                    {
                      boolean = less_than;
                    }
                }
              else if (dt1->date.day > dt2->date.day)
                {
                  boolean = greater_than;
                }
              else
                {
                  boolean = less_than;
                }
            }
          else if (dt1->date.month > dt2->date.month)
            {
              boolean = greater_than;
            }
          else
            {
              boolean = less_than;
            }
        }
      else if (dt1->date.year > dt2->date.year)
        {
          boolean = greater_than;
        }
      else
        {
          boolean = less_than;
        }

      return boolean;
    }
  else
    return compare_error;
}

/**
 * @brief COPY a DateTime object.
 *
 * Routine replaceDateTime copies the contents of source DateTime into a Destination DateTime object.
 *
 * @param  dtsrc
 *         A pointer to struct _datetime. Copy "FROM" DateTime object.
 *
 * @param  dtdest
 *         A pointer to struct _datetime. Copy "TO" DateTime object.
 *
 * @return dtdest
 *         A pointer to 'copied' DateTime Object.
 */

struct _datetime *
replaceDatetime(struct _datetime *dtsrc, struct _datetime *dtdest)
{
  if ((dtdest != NULL) && (dtsrc != NULL))
    {

      dtdest->date.year = dtsrc->date.year;
      dtdest->date.month = dtsrc->date.month;
      dtdest->date.day = dtsrc->date.day;
      dtdest->time.hour = dtsrc->time.hour;
      dtdest->time.minute = dtsrc->time.minute;
      dtdest->time.second = dtsrc->time.second;
      dtdest->time.ms = dtsrc->time.ms;

      return dtdest;
    }
  else
    return NULL;
}

/*
 * NOTE: Internal and not doxyfied.
 *
 * @brief Copy the contents of Date into a DateTime object.
 *
 * Routine convertDateToDateTime copies the contents of date object into a datetime object
 * and sets the time parameters to 0. For e.g A date 1999-12-01 can be copied into a DateTime
 * object resulting in 1999-12-01T00:00:00.000. This function can be used by Date objects
 * to get an equivalent DateTime object and then call routines supporting DateTime.
 *
 * @param  d
 *         A pointer to struct _date. Copy "FROM" date object.
 *
 * @param  dt_return
 *         A pointer to struct _datetime. Copy "TO" datetime object.
 *
 * @return dt_return
 *         A pointer to the copied datetime object.
 */

struct _datetime *
convertDateToDateTime(struct _date *d, struct _datetime *dt_return)
{
  if ((d != NULL) && (dt_return != NULL))
    {

      dt_return->date.year = d->year;
      dt_return->date.month = d->month;
      dt_return->date.day = d->day;

      /* Add dummy time. */
      dt_return->time.hour = 0;
      dt_return->time.minute = 0;
      dt_return->time.second = 0;
      dt_return->time.ms = 0;

      return dt_return;
    }
  else
    return NULL;
}

/*
 * Note: Internal and not doxyfied.
 *
 * @brief Copy the contents of DateTime into a Date object.
 *
 * Routine convertDateTimeToDate copies the date elements of DateTime object(Year,Month,Day)
 * into a Date object and discards the time elements (Hour, minute, second, ms). For e.g
 * A DateTime 1999-12-01T00:00:00.000 can be converted into a Date Object resulting in
 * 1999-12-01. (This function is used to 'go-back-to' Date format! Date objects, which convert
 * to DateTime using the routine convertDateToDateTime can return to the original Date
 * format after relevant processing on the meta-DateTime object.)
 *
 * @param  dt
 *         A pointer to struct _datetime. Copy "FROM" datetime object.
 *
 * @param  d_return
 *         A pointer to struct _date. Copy "TO" date object.
 *
 * @return d_return
 *         A pointer to the copied date object.
 */

struct _date *
convertDateTimeToDate(struct _datetime *dt, struct _date *d_return)
{
  if ((dt != NULL) && (d_return != NULL))
    {

      /* Copy Date elements. Discard Time elements. */
      d_return->year = dt->date.year;
      d_return->month = dt->date.month;
      d_return->day = dt->date.day;

      return d_return;
    }
  else
    return NULL;
}

/*! \endcond */

/**
 * @brief Get the 'day-of-year' value of a DateTime.
 *
 * Routine getDayOfYearFromDateTime returns Day of Year for the DateTime. This routine supports
 * all Calendar types.
 *
 * For eg. the day of year value for 2001-10-15T00:00:00.000 will be 288 for Gregorian Calendar.
 * Similarly, this value will be 285 for Calendar of type 360 day-Calendar.
 *
 * @param  dt
 *         A pointer to struct _datetime. Retrieve the 'day-of-year' from this DT object.
 *
 * @return doy
 *         Integer value of doy. The value depends on the calendar type. Zero indicates error.
 */

int
getDayOfYearFromDateTime(struct _datetime *dt)
{
  if (dt != NULL)
    {
      int doy = 0;

      /* Get current Calendar type. */
      calendarType ct = getCalendarType();

      if (ct == PROLEPTIC_GREGORIAN)
        {
          if (testYearIsLeapYear(dt->date.year))
            {
              /* Leap year. */
              doy = nofDaysAfterARGMonthsInLeapYear[dt->date.month - 1] + dt->date.day;
            }
          else
            {
              /* Non-leap year. */
              doy = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1] + dt->date.day;
            }
        }
      else if (ct == YEAR_OF_365_DAYS)
        {
          /* Non-leap year characteristics. */
          doy = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1] + dt->date.day;
        }
      else if (ct == YEAR_OF_360_DAYS)
        {
          /* Each month has 30 days. */
          doy = (dt->date.month - 1) * NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 + dt->date.day;
        }

      return doy;
    }
  else
    return 0;
}

/**
 * @brief Get nod (number of Days) in the month of DateTime.
 *
 * Routine getNoOfDaysInMonthDateTime returns number of days in the month of DateTime. This routine
 * supports all calendar types.
 *
 * For eg. the number of days for 2001-10-15T00:00:00.000 will be 31 for Gregorian Calendar.
 * Similarly, this value will be 30 for Calendar of type 360 day-Calendar.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return nod
 *         Integer value of nod. The value depends on the month and the calendar type. Zero indicates error.
 */
// TODO on Luis: Is this function doing the right thing?

int
getNoOfDaysInMonthDateTime(struct _datetime *dt)
{
  if (dt != NULL)
    {
      int nod = 0;

      /* Get current Calendar type. */
      calendarType ct = getCalendarType();

      if (ct == PROLEPTIC_GREGORIAN)
        {
          if (testYearIsLeapYear(dt->date.year))
            {
              /* Leap year. */
              nod = nofDaysAfterARGMonthsInLeapYear[dt->date.month] - nofDaysAfterARGMonthsInLeapYear[dt->date.month - 1];
            }
          else
            {
              /* Non leap year. */
              nod = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month] - nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1];
            }
        }
      else if (ct == YEAR_OF_365_DAYS)
        {
          /* Non-leap year. */
          nod = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month] - nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1];
        }
      else if (ct == YEAR_OF_360_DAYS)
        {
          /* Each month has 30 days. */
          nod = NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360;
        }

      return nod;
    }
  else
    return 0;
}

/**
 * @brief Get number of days in the Year of DateTime.
 *
 * Routine getNoOfDaysInYearDateTime returns number of days in the Year of DateTime. This routine
 * supports all calendar types.
 *
 * Number of days returned will depend on the calendar type and if applicable, leap v/s non leap year.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return nod
 *         Integer value of nod. The value depends on the year and the calendar type. Zero indicates error.
 */
// TODO on Luis: Is this function doing the right thing?
int
getNoOfDaysInYearDateTime(struct _datetime *dt)
{
  if (dt != NULL)
    {
      int nod = 0;

      /* Get current Calendar type. */
      calendarType ct = getCalendarType();

      if (ct == PROLEPTIC_GREGORIAN)
        {
          if (testYearIsLeapYear(dt->date.year))
            {
              /* Leap year. */
              nod = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Non leap year. */
              nod = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }
        }
      else if (ct == YEAR_OF_365_DAYS)
        {
          /* Non Leap year. */
          nod = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
        }
      else if (ct == YEAR_OF_360_DAYS)
        {
          /* Year has 360 days. */
          nod = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360;
        }

      return nod;
    }
  else
    return 0;
}

/**
 * @brief Get number of seconds elapsed in the month of DateTime.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return no_of_seconds
 *         int64_t value of no_of_seconds. -1 indicates error.
 */

int64_t
getNoOfSecondsElapsedInMonthDateTime(struct _datetime *dt)
{
  if (dt != NULL)
    return ((dt->date.day - 1) * (NO_OF_MS_IN_A_DAY / NO_OF_MS_IN_A_SECOND)
            + dt->time.hour * (NO_OF_MS_IN_A_HOUR / NO_OF_MS_IN_A_SECOND)
            + dt->time.minute * (NO_OF_MS_IN_A_MINUTE / NO_OF_MS_IN_A_SECOND) + dt->time.second);
  else
    return -1;
}

/**
 * @brief Get number of seconds elapsed in the day of DateTime.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return no_of_seconds
 *         int value of no_of_seconds. -1 indicates error.
 */

int
getNoOfSecondsElapsedInDayDateTime(struct _datetime *dt)
{
  if (dt != NULL)
    return (dt->time.hour * (NO_OF_MS_IN_A_HOUR / NO_OF_MS_IN_A_SECOND)
            + dt->time.minute * (NO_OF_MS_IN_A_MINUTE / NO_OF_MS_IN_A_SECOND) + dt->time.second);
  else
    return -1;
}

/**
 * @brief Get the Julian Day from DateTime.
 *
 * The routine getJulianDayFromDateTime returns the equivalent Julian date to DateTime. Internally
 * it calls translation routines based on Calendar type.
 *
 * @param  dt
 *         A pointer to struct _datetime. The DT's value is converted to julian day value.
 *
 ** @param  jd
 *         A pointer to struct _julianday. JD where the converted value is stored.
 *
 * @return jd
 *         A pointer to struct _julianday containing a copy of the JD value corresponding to the DT.
 */

struct _julianday *
getJulianDayFromDateTime(struct _datetime *dt, struct _julianday *jd)
{
  /* Function pointer date2julian points to the correct translation routine. */
  return date2julian(dt, jd);
}

/**
 * @brief Get the DateTime from Julian Day.
 *
 * The routine getDateTimeFromJulianDay returns the equivalent DateTime to Julian date. Internally
 * it calls translation routines based on Calendar type.
 *
 * @param  jd
 *         A pointer to struct _julianday. The JD's value is converted to julian day value.
 *
 ** @param  dt
 *         A pointer to struct _datetime. The DT where the converted value is stored.
 *
 * @return jd
 *         A pointer to struct _datetime containing a copy of the DT value corresponding to the JD.
 */

struct _datetime *
getDateTimeFromJulianDay(struct _julianday *jd, struct _datetime *dt)
{
  /* Function pointer julian2date points to the correct translation routine. */
  return julian2date(jd, dt);
}

/**
 * @brief Get DateTime as a string.
 *
 * datetimeToString returns a string in IS08601 compliant (and extended) format.
 *
 * @param  dt
 *         A pointer to struct _datetime. The datetime to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where datetime is to be written.
 *
 * @return toStr
 *         A pointer to the string containing datetime.
 */

char *
datetimeToString(struct _datetime *dt, char *toStr)
{
  if ((dt != NULL) && (toStr != NULL))
    {

      memset(toStr, '\0', MAX_DATETIME_STR_LEN);

      snprintf(toStr, MAX_DATETIME_STR_LEN, "%" PRIi64 "-%02d-%02dT%02d:%02d:%02d.%03d", dt->date.year, dt->date.month,
               dt->date.day, dt->time.hour, dt->time.minute, dt->time.second, dt->time.ms);

      return toStr;
    }
  else
    return NULL;
}

/**
 * @brief Get DateTime as a basic string.
 *
 * datetimeToBasicString returns a string in IS08601 compliant format.
 *
 * @param  dt
 *         A pointer to struct _datetime. The datetime to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where datetime is to be written.
 *
 * @return toStr
 *         A pointer to the string containing datetime.
 */

char *
datetimeToBasicString(struct _datetime *dt, char *toStr)
{
  if ((dt != NULL) && (toStr != NULL))
    {

      memset(toStr, '\0', MAX_DATETIME_STR_LEN);

      if (dt->date.year < 0)
        snprintf(toStr, MAX_DATETIME_STR_LEN, "% 05d%02d%02dT%02d%02d%02d.%03d", (int) dt->date.year, dt->date.month, dt->date.day,
                 dt->time.hour, dt->time.minute, dt->time.second,
                 dt->time.ms);  // TODO on LUIS: Should we block when year is > 9999 format?
      else
        snprintf(toStr, MAX_DATETIME_STR_LEN, "%04d%02d%02dT%02d%02d%02d.%03d", (int) dt->date.year, dt->date.month, dt->date.day,
                 dt->time.hour, dt->time.minute, dt->time.second, dt->time.ms);

      return toStr;
    }
  else
    return NULL;
}

/**
 * @brief Get DateTime in 'struct tm' format and return as a string.
 *
 * Only dates between and including 1582-10-15 TO 9999-12-31 supported.
 *
 * @param  dt
 *         A pointer to struct _datetime. The datetime to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where datetime is to be written.
 *
 * @param  fmtString
 *         A pointer to char. Desired Format string. CRITICAL: Inappropriate fmt string will cause dump.
 *
 * @return toStr
 *         A pointer to the string containing datetime.
 */

char *
datetimeToPosixString(struct _datetime *dt, char *toStr, char *fmtString)
{
  if ((dt != NULL) && (toStr != NULL) && (fmtString != NULL))
    {
      struct tm tm_info;

      /*Range check. Return with NULL indicating Error */
      if (dt->date.year < POSIXSTRING_YEAR_LOWER_BOUND)
        {
          return NULL;
        }
      else if (((dt->date.year == POSIXSTRING_YEAR_LOWER_BOUND) && (dt->date.month < POSIXSTRING_MONTH_LOWER_BOUND))
               || ((dt->date.year == POSIXSTRING_YEAR_LOWER_BOUND) && (dt->date.month == POSIXSTRING_MONTH_LOWER_BOUND)
                   && (dt->date.day < POSIXSTRING_DAY_LOWER_BOUND)))
        {
          return NULL;
        }
      else if (dt->date.year > POSIXSTRING_YEAR_UPPER_BOUND)
        {
          return NULL;
        }

      tm_info.tm_sec = dt->time.second;
      tm_info.tm_min = dt->time.minute;
      tm_info.tm_hour = dt->time.hour;

      tm_info.tm_mday = dt->date.day;
      /* Range of month is from 0 to 11. */
      tm_info.tm_mon = dt->date.month - 1;
      /* tm's year is w.r.t the year 1900. */
      tm_info.tm_year = (int) (dt->date.year - 1900);

      memset(toStr, '\0', MAX_DATETIME_STR_LEN);
      strftime(toStr, MAX_DATETIME_STR_LEN, fmtString,
               &tm_info); /*CRITICAL: Careful what you wish for. Incompatiable fmtString will cause core dump*/

      return toStr;
    }
  else
    return NULL;
}

/* Let [A,B] denote a range. Let X be a DateTime. If X < A , return less_than. If X > B, return greater_than.
   Else X lies within [A,B] => return equal_to.
   If B is NULL, B can be thought of as infinity: If X < A , return less_than else return equal_to. */
compare_return_val
getDateTimeIsInRange(struct _datetime *dtRef, struct _datetime *dtStart, struct _datetime *dtEnd)
{
  if ((dtRef != NULL) && (dtStart != NULL))  // dtEnd can be NULL.
    {
      /* Init flag. */
      compare_return_val flag = compare_error;

      if (compareDatetime(dtRef, dtStart) == less_than)
        {
          /* X < A*/
          flag = less_than;
        }
      else if (compareDatetime(dtRef, dtEnd) == greater_than)
        {
          /* X > B */
          flag = greater_than;
        }
      else
        {
          /* A <= X <= B. Notice that if dtEnd is NULL, then compareDatetime()
             will return with error != greater_than and hence we return equal_to as desired. */
          flag = equal_to;
        }
      return flag;
    }
  else
    return compare_error;
}
