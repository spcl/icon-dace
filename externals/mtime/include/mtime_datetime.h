// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_datetime.h
 *
 * @brief DateTime and some operations supported on DateTime.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 */

#ifndef _MTIME_DATETIME_H
#define _MTIME_DATETIME_H

#include <stdint.h>
#include <stdbool.h>

#include "mtime_date.h"
#include "mtime_time.h"

/**
 * @struct _datetime
 *
 * @brief struct _datetime contains a struct _date and a struct _time element.
 */

struct _datetime
{
  struct _date date;  ///< Date elements.
  struct _time time;  ///< Time elements.
};

struct _datetime *newDateTime(const char *datetime_string);

struct _datetime *newRawDateTime(int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms);

struct _datetime *constructAndCopyDateTime(struct _datetime *dt);

void deallocateDateTime(struct _datetime *dt);

compare_return_val compareDatetime(struct _datetime *dt1, struct _datetime *dt2);

struct _datetime *replaceDatetime(struct _datetime *dtsrc, struct _datetime *dtdest);

char *datetimeToString(struct _datetime *dt, char *toStr);

char *datetimeToBasicString(struct _datetime *dt, char *toStr);

char *datetimeToPosixString(struct _datetime *dt, char *toStr, char *fmtString);

int getNoOfDaysInMonthDateTime(struct _datetime *dt);

int getNoOfDaysInYearDateTime(struct _datetime *dt);

/*! \cond PRIVATE */
static inline bool
testYearIsLeapYear(int64_t year)
{
  bool isLeapYear = !(year % 400) || ((year % 100) && !(year % 4));
  return isLeapYear;
}

struct _datetime *convertDateToDateTime(struct _date *d, struct _datetime *dt_return);

struct _date *convertDateTimeToDate(struct _datetime *dt, struct _date *d_return);
/*! \endcond */

int getDayOfYearFromDateTime(struct _datetime *currentdt);

int64_t getNoOfSecondsElapsedInMonthDateTime(struct _datetime *dt);

int getNoOfSecondsElapsedInDayDateTime(struct _datetime *dt);

struct _julianday *getJulianDayFromDateTime(struct _datetime *dt, struct _julianday *jd);

struct _datetime *getDateTimeFromJulianDay(struct _julianday *jd, struct _datetime *dt);

compare_return_val getDateTimeIsInRange(struct _datetime *dtRef, struct _datetime *dtStart, struct _datetime *dtEnd);

/**
 * @}
 */

#endif
