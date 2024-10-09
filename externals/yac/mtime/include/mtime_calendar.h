// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_calendar.h
 *
 * @brief Singleton Calendar connecting all supported calendar types.
 *
 * @details
 *
 * Three calendar types are provided:
 *
 *   - a proleptic Gregorian calendar
 *   - a calendar with 365 days per year without leap years
 *   - a calendar with 360 days per year and each month having 30 days
 *
 * To use this library, a call to initCalendar() with the respective selector (enum calendarType) must be
 * done first. The implementation is based on a singleton concept meaning that only one calendar
 * can be active at a time. To release a calendar a call to freeCalendar() has to be done.
 *
 *
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology
 * @author  Rahul Sinha, Max Planck Institute for Meteorology
 *
 * @date March 2013
 *
 * @note Calendar type, once initialized, should not be changed.
 */

#ifndef _MTIME_CALENDAR_H
#define _MTIME_CALENDAR_H

extern const int NO_OF_SEC_IN_A_DAY;
extern const int NO_OF_SEC_IN_A_HOUR;
extern const int NO_OF_SEC_IN_A_MINUTE;

extern const int NO_OF_MS_IN_A_DAY;
extern const int NO_OF_MS_IN_HALF_DAY;
extern const int NO_OF_MS_IN_A_HOUR;
extern const int NO_OF_MS_IN_A_MINUTE;
extern const int NO_OF_MS_IN_A_SECOND;

/* #define NO_OF_SEC_IN_A_DAY               86400 */
/* #define NO_OF_SEC_IN_A_HOUR              3600 */
/* #define NO_OF_SEC_IN_A_MINUTE            60 */

/* #define NO_OF_MS_IN_A_DAY 		86400000 */
/* #define NO_OF_MS_IN_HALF_DAY 		43200000 */
/* #define NO_OF_MS_IN_A_HOUR 		3600000 */
/* #define NO_OF_MS_IN_A_MINUTE 		60000 */
/* #define NO_OF_MS_IN_A_SECOND 		1000 */

#define NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 30

#define NO_OF_MONTHS_IN_A_YEAR 12
#define NO_OF_HOURS_IN_A_DAY 24

#define NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360 360
#define NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365 365
#define NO_OF_DAYS_IN_A_LEAP_YEAR 366

#define NO_OF_DAYS_IN_400_YEARS 146097

#define YEAR_UPPER_BOUND 2147483647
#define YEAR_LOWER_BOUND -2147483648

#define POSIXSTRING_DAY_LOWER_BOUND 15
#define POSIXSTRING_MONTH_LOWER_BOUND 10
#define POSIXSTRING_YEAR_LOWER_BOUND 1582
#define POSIXSTRING_YEAR_UPPER_BOUND 9999

/// provides a string length for toString.
#define MAX_CALENDAR_STR_LEN 32
#define MAX_DATETIME_STR_LEN 32
#define MAX_DATE_STR_LEN 32
#define MAX_TIME_STR_LEN 32
#define MAX_TIMEDELTA_STR_LEN 32
#define MAX_JULIANDAY_STR_LEN 32

struct _datetime;
struct _julianday;

/* Used as dummy string to init TimeDeltas in the lib. */
extern const char initDummyTDString[];
extern const char initDummyDString[];
extern const char initDummyDTString[];

extern const int nofDaysAfterARGMonthsInLeapYear[13];
extern const int nofDaysAfterARGMonthsInNonLeapYear[13];

extern const int nofDaysInARGMonthInLeapYear[];
extern const int nofDaysInARGMonthIn365DayYear[];
extern const int nofDaysInARGMonthIn360DayYear[];

extern const int monthSpecificDeltaInMonthsLeapyear[12][13];
extern const int monthSpecificDeltaInMonths365[12][13];
extern const int monthSpecificDeltaInMonths360[12][13];

typedef enum
{
  equal_to = 0,
  greater_than = 1,
  less_than = -1,
  compare_error = -128
} compare_return_val;

/**
 * @enum calendarType
 *
 * @brief enum calendarType lists the calendarTypes supported. The values are used for selecting calendars.
 *
 */

typedef enum
{
  CALENDAR_NOT_SET = 0,     ///< calendar is not defined yet.
  PROLEPTIC_GREGORIAN = 1,  ///< proleptic Gregorian calendar.
  YEAR_OF_365_DAYS = 2,     ///< 365 day year without leap years.
  YEAR_OF_360_DAYS = 3      ///< 360 day year with 30 day months.
} calendarType;

///@brief Function pointer connecting the Calendar to Julian routine. The pointed-to Function depends on the selected Calendar type.
extern struct _julianday *(*date2julian)(struct _datetime *date, struct _julianday *julian);

///@brief Function pointer connecting Julian to Calendar routine. The pointed to Function depends on the selected Calendar type.
extern struct _datetime *(*julian2date)(struct _julianday *julian, struct _datetime *date);

void initCalendar(calendarType ct);

void freeCalendar(void);

calendarType getCalendarType(void);

char *calendarToString(char *calendar);

/**
 * @}
 */
#endif
