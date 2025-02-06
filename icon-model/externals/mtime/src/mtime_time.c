// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @addtogroup CBindings libmtime C language bindings
 *
 * @file mtime_time.c
 *
 * @brief Time and some operations supported on Time.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mtime_time.h"

#include "mtime_calendar.h"
#include "mtime_iso8601.h"

/* Upto MilliSecond resolution supported. */

/**
 * @brief Construct new Time using an ISO 8601 conforming string.
 *
 *         Strings of the following form are valid:
 *
 *         <table>
 *         <tr><th> Extended Format: <th>
 *         <tr><td> HH:MM:SS.SSS <td> .SSS can be of length 1-3
 *         <tr><td> HH:MM:SS     <td>
 *         <tr><td> HH:MM        <td>
 *         <tr><td> HH           <td>
 *         </table>
 *
 *         <table>
 *         <tr><th> Basic Format: </th>
 *         <tr><td> HHMMSS.SSS  <td> .SSS can be of length 1-3
 *         <tr><td> HHMMSS      <td>
 *         <tr><td> HHMM        <td>
 *         <tr><td> HH          <td>
 *         </table>
 *
 *         In essence, strings follow the guidelines in ISO8601 standard (but do not support all formats).
 *
 * @param  ts
 *         A pointer to char. The string contains parameters with which Time is created.
 *
 * @return t
 *         A pointer to a filled Time.
 */

struct _time *
newTime(const char *ts)
{
  if ((ts != NULL) && (getCalendarType() && (strlen(ts) != 0)))
    {
      /* Convert ts to dts by appending dummy Date 0-01-01 for testing with verify_string_datetime(). */
      char *dts = (char *) calloc(MAX_DATETIME_STR_LEN, sizeof(char));
      if (dts == NULL) return NULL;

      if (strstr(ts, ":")) /* Extended Format.  */
        strcpy(dts, "0000-01-01T");
      else /* Basic Format.     */
        strcpy(dts, "00000101T");

      strncat(dts, ts, MAX_DATETIME_STR_LEN - 8);

      struct iso8601_datetime *isoDt = new_iso8601_datetime('+', 0, 0, 0, 0, 0, 0, 0);
      if (isoDt == NULL)
        {
          free(dts);
          dts = NULL;
          return NULL;
        }

      /* Verify ISO 8601 compliance. */
      if (verify_string_datetime(dts, isoDt) != DATETIME_MATCH)
        {
          free(dts);
          dts = NULL;
          deallocate_iso8601_datetime(isoDt);
          isoDt = NULL;
          return NULL;
        }

      struct _time *t = (struct _time *) calloc(1, sizeof(struct _time));
      if (t == NULL)
        {
          free(dts);
          dts = NULL;
          deallocate_iso8601_datetime(isoDt);
          isoDt = NULL;
          return NULL;
        }

      t->hour = isoDt->hour;
      t->minute = isoDt->minute;
      t->second = isoDt->second;
      t->ms = isoDt->ms;

      // Cleanup.
      free(dts);
      dts = NULL;
      deallocate_iso8601_datetime(isoDt);

      return t;
    }
  else
    return NULL;
}

/**
 * @brief Construct new Time using 'raw' numerical values.
 *
 * @param  _hour
 *        An "int" value denoting the hour part of time.
 * @param  _minute
 *        An "int" value denoting the minute part of time.
 * @param  _second
 *        An "int" value denoting the second part of time.
 * @param  _ms
 *        An "int" value denoting the milli-second part of time.
 *
 * @return t
 *         A pointer to a filled Time.
 */

struct _time *
newRawTime(int _hour, int _minute, int _second, int _ms)
{
  char *ts = (char *) calloc(MAX_TIME_STR_LEN, sizeof(char));
  if (ts == NULL) return NULL;

  snprintf(ts, MAX_TIME_STR_LEN, "%02d:%02d:%02d.%03d", _hour, _minute, _second, _ms);

  /* Reuse the string interface. */
  struct _time *t = newTime(ts);

  free(ts);
  ts = NULL;

  return t;
}

/**
 * @brief Copy the values and construct a new Time.
 *
 * @param  t
 *         A pointer to struct _time. Values of t are used to initialize the new time being created.
 *
 * @return _t
 *         A pointer to an initialized Time object.
 */

struct _time *
constructAndCopyTime(struct _time *t)
{
  if (t != NULL)
    return newRawTime(t->hour, t->minute, t->second, t->ms);
  else
    return NULL;
}

/**
 * @brief Destructor of Time.
 *
 * @param  t
 *         A pointer to struct _time. t is deallocated.
 */

void
deallocateTime(struct _time *t)
{
  if (t != NULL)
    {
      free(t);
      t = NULL;
    }
}

/**
 * @brief COPY a time object.
 *
 * Routine replaceTime copies the contents of source Time into a Destination Time object.
 *
 * @param  tsrc
 *         A pointer to struct _time. Copy "FROM" time object.
 *
 * @param  tdest
 *        A pointer to struct _time. Copy "TO" time object.
 *
 * @return tdest
 *         A pointer to 'copied' time Object.
 */

struct _time *
replaceTime(struct _time *tsrc, struct _time *tdest)
{
  if ((tdest != NULL) && (tsrc != NULL))
    {

      tdest->hour = tsrc->hour;
      tdest->minute = tsrc->minute;
      tdest->second = tsrc->second;
      tdest->ms = tsrc->ms;

      return tdest;
    }
  else
    return NULL;
}

/**
 * @brief Get time as an extended string.
 *
 * timetoString returns a string in IS08601 compliant (and extended) format.
 *
 * @param  t
 *         A pointer to struct _time. The time to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where time is to be written.
 *
 * @return toStr
 *         A pointer to the string containing time.
 */

char *
timeToString(struct _time *t, char *toStr)
{
  if ((t != NULL) && (toStr != NULL))
    {
      memset(toStr, '\0', MAX_TIME_STR_LEN);

      snprintf(toStr, MAX_TIME_STR_LEN, "%02d:%02d:%02d.%03d", t->hour, t->minute, t->second, t->ms);

      return toStr;
    }
  else
    return NULL;
}

/**
 * @brief Get time as a basic string.
 *
 * timetoString returns a string in IS08601 compliant format.
 *
 * @param  t
 *         A pointer to struct _time. The time to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where time is to be written.
 *
 * @return toStr
 *         A pointer to the string containing time.
 */

char *
timeToBasicString(struct _time *t, char *toStr)
{
  if ((t != NULL) && (toStr != NULL))
    {
      memset(toStr, '\0', MAX_TIME_STR_LEN);

      snprintf(toStr, MAX_TIME_STR_LEN, "%02d%02d%02d.%03d", t->hour, t->minute, t->second, t->ms);

      return toStr;
    }
  else
    return NULL;
}

/**
 * @brief Get Time in 'struct tm' format and return as a string.
 *
 * @param  t
 *         A pointer to struct _time. The time to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where time is to be written.
 *
 * @param  fmtString
 *         A pointer to char. Desired Format string. CRITICAL: Inappropriate fmt string will cause dump.
 *
 * @return toStr
 *         A pointer to the string containing time.
 */

char *
timeToPosixString(struct _time *t, char *toStr, char *fmtString)
{
  if ((t != NULL) && (toStr != NULL))
    {
      struct tm tm_info;

      tm_info.tm_sec = t->second;
      tm_info.tm_min = t->minute;
      tm_info.tm_hour = t->hour;

      memset(toStr, '\0', MAX_TIME_STR_LEN);
      strftime(toStr, MAX_TIME_STR_LEN, fmtString,
               &tm_info); /*CRITICAL: Careful what you wish for. Incompatiable fmtString will cause core dump*/

      return toStr;
    }
  else
    return NULL;
}
