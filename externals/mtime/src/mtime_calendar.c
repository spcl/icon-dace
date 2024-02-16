// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @addtogroup CBindings libmtime C language bindings
 *
 * @file mtime_calendar.c
 * 
 * @brief Singleton Calendar connecting all supported calendar types.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note Calendar type, once initialized, should not be changed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "mtime_calendar.h"

#include "mtime_calendarGregorian.h"
#include "mtime_calendar365day.h"
#include "mtime_calendar360day.h"

const int __attribute__ ((aligned (8))) NO_OF_SEC_IN_A_DAY    =    86400;
const int __attribute__ ((aligned (8))) NO_OF_SEC_IN_A_HOUR   =     3600;
const int __attribute__ ((aligned (8))) NO_OF_SEC_IN_A_MINUTE =       60;

const int __attribute__ ((aligned (8))) NO_OF_MS_IN_A_DAY 	= 86400000;
const int __attribute__ ((aligned (8))) NO_OF_MS_IN_HALF_DAY 	= 43200000;
const int __attribute__ ((aligned (8))) NO_OF_MS_IN_A_HOUR 	=  3600000;
const int __attribute__ ((aligned (8))) NO_OF_MS_IN_A_MINUTE 	=    60000;
const int __attribute__ ((aligned (8))) NO_OF_MS_IN_A_SECOND 	=     1000;

/* Used as dummy string to init data-structures in the lib. */
//WARNING: DO NOT EDIT THESE STRINGS (unless you know what you are doing ofcourse). This will break the lib.
const char initDummyDTString[] = "0-01-01T00:00:00.000";
const char initDummyDString[]  = "0-01-01";
const char initDummyTDString[] = "PT00.000S";


/* Number of days after 0,1,2...12 months in a Leap year.*/
const int nofDaysAfterARGMonthsInLeapYear[13] =
  { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };

/* Number of days after 0,1,2...12 months in a non-Leap year*/
const int nofDaysAfterARGMonthsInNonLeapYear[13] =
  { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };


/* #days in each month: Jan = 31, Feb = 28 ... */
const int nofDaysInARGMonthInLeapYear[] =
  { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const int nofDaysInARGMonthIn365DayYear[] =
  { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const int nofDaysInARGMonthIn360DayYear[] =
  { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 };


/* Each row represents months 1,2,3..12; Each column entry is the number
 of days after 0,1,2,3...12 months from the corresponding month. */
const int monthSpecificDeltaInMonthsLeapyear[12][13] =
  {
    { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 },
    { 0, 29, 60, 90, 121, 151, 182, 213, 243, 274, 304, 335, 366 },
    { 0, 31, 61, 92, 122, 153, 184, 214, 245, 275, 306, 337, 366 },
    { 0, 30, 61, 91, 122, 153, 183, 214, 244, 275, 306, 335, 366 },
    { 0, 31, 61, 92, 123, 153, 184, 214, 245, 276, 305, 336, 366 },
    { 0, 30, 61, 92, 122, 153, 183, 214, 245, 274, 305, 335, 366 },
    { 0, 31, 62, 92, 123, 153, 184, 215, 244, 275, 305, 336, 366 },
    { 0, 31, 61, 92, 122, 153, 184, 213, 244, 274, 305, 335, 366 },
    { 0, 30, 61, 91, 122, 153, 182, 213, 243, 274, 304, 335, 366 },
    { 0, 31, 61, 92, 123, 152, 183, 213, 244, 274, 305, 336, 366 },
    { 0, 30, 61, 92, 121, 152, 183, 213, 243, 274, 305, 335, 366 },
    { 0, 31, 62, 91, 122, 152, 183, 213, 244, 275, 305, 336, 366 } };

const int monthSpecificDeltaInMonths365[12][13] =
  {
    { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 },
    { 0, 28, 59, 89, 120, 150, 181, 212, 242, 273, 303, 334, 365 },
    { 0, 31, 61, 92, 122, 153, 184, 214, 245, 275, 306, 337, 365 },
    { 0, 30, 61, 91, 122, 153, 183, 214, 244, 275, 306, 334, 365 },
    { 0, 31, 61, 92, 123, 153, 184, 214, 245, 276, 304, 335, 365 },
    { 0, 30, 61, 92, 122, 153, 183, 214, 245, 273, 304, 334, 365 },
    { 0, 31, 62, 92, 123, 153, 184, 215, 243, 274, 304, 335, 365 },
    { 0, 31, 61, 92, 122, 153, 184, 212, 243, 273, 304, 334, 365 },
    { 0, 30, 61, 91, 122, 153, 181, 212, 242, 273, 303, 334, 365 },
    { 0, 31, 61, 92, 123, 151, 182, 212, 243, 273, 304, 335, 365 },
    { 0, 30, 61, 92, 120, 151, 182, 212, 242, 273, 304, 334, 365 },
    { 0, 31, 62, 90, 121, 151, 182, 212, 243, 274, 304, 335, 365 } };

const int monthSpecificDeltaInMonths360[12][13] =
  {
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 },
    { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360 } };


/* Function pointers connecting translation routines according to calendarType. */
struct _julianday *
(*date2julian)(struct _datetime *date, struct _julianday *julian);
struct _datetime *
(*julian2date)(struct _julianday *julian, struct _datetime *date);

/* Static variable describing the calendar type. */
static calendarType _calendar_type = CALENDAR_NOT_SET;
/* Flag to disallow change in calendar type.  */
static bool calendar_initialized = false;

/**
 * @brief To initialize a new calendar
 *
 * Calendar init is done at the very begining to use this calendar lib. Init intializes the 
 * calendar to:
 * - PROLEPTIC_GREGORIAN
 * - YEAR_OF_365_DAYS
 * - OR YEAR_OF_360_DAYS Calendar.
 * 
 * The calendar type and hence it's behaviour (Calendar to Julian conversion and vice versa) is fixed for the lifetime of the calendar object.
 * Attempts to change the calendar type at runtime is discouraged. The lib has built-in checks to reject change attempts.
 * Calendar can be "re-initialized" after calling freeCalendar() but is not advised. 
 * MANTRA: Know what you are doing before you do it and do it right the first time.
 *
 * @param  ct
 *         An object of struct calendarType. 
 *	   ct 'fixes' the calendar type by setting function pointers to point to functions designed
 *	   to work with that calendarType.
 */
void
initCalendar(calendarType ct)
{
  //DO not allow calendar type to be modified. If calendar was previously set, override the user specified ct with sytem stored type.
  if (calendar_initialized) 
    {
      ct = _calendar_type;
    }

  switch (ct)
    {
  case PROLEPTIC_GREGORIAN:
    date2julian = getJulianFromDateGregorian;
    julian2date = getDateGregorianFromJulian;
    _calendar_type = PROLEPTIC_GREGORIAN;
    calendar_initialized = true;
    break;
  case YEAR_OF_365_DAYS:
    date2julian = getJulianFromDate365;
    julian2date = getDate365FromJulian;
    _calendar_type = YEAR_OF_365_DAYS;
    calendar_initialized = true;
    break;
  case YEAR_OF_360_DAYS:
    date2julian = getJulianFromDate360;
    julian2date = getDate360FromJulian;
    _calendar_type = YEAR_OF_360_DAYS;
    calendar_initialized = true;
    break;
  case CALENDAR_NOT_SET:
  default:
    break;
    }

  return;
}


/** 
 * @brief called to discard the selected calendar type.
 *
 * RESETs the calendar. Should be performed ONLY at exit! 
 * WARNING: Freeing the calendar and re-assigning a new calendar type not supported.
 */
void
freeCalendar(void)
{
  _calendar_type = CALENDAR_NOT_SET;
  calendar_initialized = false;
  return;
}


/**
 * @brief To query the current calendar type.
 * 
 * Gets the calendarType. 
 *
 * @return _calendar_type
 *         A member of enum calendarType. Denotes the currently set calendar type.
 */
calendarType
getCalendarType(void)
{
  return _calendar_type;
}


/** 
 * @brief convert the calendar identifier into a human readable string. 
 *
 * Gets a string describing the calendar in use. 
 *
 * @param  calendar
 *         A pointer to char. A string where the currently set calendar type will be written.
 */

char*
calendarToString(char *calendar)
{
  if (calendar != NULL)
    {
      memset(calendar,'\0',MAX_CALENDAR_STR_LEN);
      switch (_calendar_type)
        {
          case PROLEPTIC_GREGORIAN:
            strcpy(calendar, "Proleptic Gregorian");
            break;
          case YEAR_OF_365_DAYS:
            strcpy(calendar, "Year of 365 days");
            break;
          case YEAR_OF_360_DAYS:
            strcpy(calendar, "Year of 360 days");
            break;
          case CALENDAR_NOT_SET:
          default:
            strcpy(calendar, "Not defined");
            break;
        } 
      return calendar;
    }
  else
    return NULL;
}
