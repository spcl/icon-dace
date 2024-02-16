/**
 * @file mtime_calendar365day.c
 * 
 * @brief Convert julian dates to Calendar-with-365-days dates and vice-versa.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note Calendar-with-365-days has a non-leap year characteristic and has 365 days in a calendar year. 
 * Number of days in 12 months is assumed to be: {  31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }.
 * Julian Day (0,0) corresponds to Calendar-with-365-day's day 0-01-01T12:00:00.000.
 */

#include <stdint.h>
#include <stddef.h>

#include "mtime_calendar365day.h"
#include "mtime_julianDay.h"
#include "mtime_datetime.h"


/**
 * @brief Convert Julian Date to Calendar-with-365-days date.
 *
 * The getDate365FromJulian routine accepts a date on Julian axis and translates it to a date
 * on Calendar-of-365-days. For eg. If Julian Day values are day = 1 and ms = 0, the corresponding 
 * Calendar-of-365-days values will be Year = 0, Month = 1, Day = 2, Hour = 12, Minute = 0, 
 * Second = 0 and MS = 0.
 *
 * The routine expects non NULL parameters. In case Julian-Day or DateTime
 * is NULL, the routine returns with a NULL indicating failure.
 *
 * @param  jd 
 *         A pointer to struct _julianday. JD contains the date to be 'translated' to type d365.
 * @param  d365
 *         A pointer to struct _datetime. The translated date is copied into d365.
 * @return d365
 *         A pointer to the filled d365. 
 */

struct _datetime*
getDate365FromJulian(struct _julianday *jd, struct _datetime* d365)
{
  /* Check for malformed parameters */
  if ((jd != NULL) && (d365 != NULL)){

  int i, days;
  int64_t jday, jms;

  jday = jd->day;
  jms = jd->ms;

  /* Handle the 12 hour skew between two time axis. */
  jms = jms + NO_OF_MS_IN_HALF_DAY;
  if (jms >= NO_OF_MS_IN_A_DAY)
    {
      jday = jday + 1;
      jms = jms - NO_OF_MS_IN_A_DAY;
    }

  /* Julian is symmetric about 0.0 JD while d360 is not. */
  if (jday >= 0)
    {
      d365->date.year = (jday) / NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;

      days = (int)(jday % NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365);
      for (i = (NO_OF_MONTHS_IN_A_YEAR - 1); i >= 0; i--)
        {
          if (days >= nofDaysAfterARGMonthsInNonLeapYear[i])
            {
              d365->date.month = i + 1;
              d365->date.day = days - nofDaysAfterARGMonthsInNonLeapYear[i] + 1;

              break;
            }
        }
    }
  else
    {
      d365->date.year = (jday / NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365) - 1;
      days = (int)(-jday % NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365);
      for (i = (NO_OF_MONTHS_IN_A_YEAR -1); i >= 0; i--)
        {
          if ( (NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365 - days + 1) > nofDaysAfterARGMonthsInNonLeapYear[i])
            {
              d365->date.month = i + 1;
              // Avoid divide by zero.
              if (nofDaysAfterARGMonthsInNonLeapYear[i] == 0)
                d365->date.day = (NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365 - days + 1);
              else
                d365->date.day = (NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365 - days + 1) % nofDaysAfterARGMonthsInNonLeapYear[i];

              break;
            }
        }

      /* This will be true whenever days is 0 */
      if (d365->date.day == (nofDaysInARGMonthIn365DayYear[i] + 1))
        {
          d365->date.month = d365->date.month + 1;
          if (d365->date.month == NO_OF_MONTHS_IN_A_YEAR + 1)
            {
              d365->date.year = d365->date.year + 1;
              d365->date.month = 1;
            }
          d365->date.day = 1;
        }
    }

  /* Santize for range. */
  if((d365->date.year > YEAR_UPPER_BOUND) || (d365->date.year < YEAR_LOWER_BOUND))
    {
      /* ERROR: Exceeds allowed year range. */
      return NULL;
    }


  d365->time.hour = (int) jms / NO_OF_MS_IN_A_HOUR;
  d365->time.minute = ((int) jms - d365->time.hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
  d365->time.second = ((int) jms - d365->time.hour * NO_OF_MS_IN_A_HOUR - d365->time.minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
  d365->time.ms = (int) jms - d365->time.hour * NO_OF_MS_IN_A_HOUR - d365->time.minute * NO_OF_MS_IN_A_MINUTE - d365->time.second * NO_OF_MS_IN_A_SECOND;

  return d365;
}
else
  return NULL;
}


/**
 * @brief Convert Calendar-with-365-days date to Julian Date.
 *
 * The getJulianFromDate365 routine accepts a date on Calendar-of-365-days axis and translates it to a date
 * on Julian axis. For eg. If Calendar-of-365-days values are Year = 0, Month = 1, Day = 2, Hour = 12, Minute = 0, 
 * Second = 0 and MS = 0 , corresponding Julian Day values will be day = 1 and ms = 0.
 *
 * The routine expects non NULL parameters. In case Julian-Day or DateTime
 * is NULL, the routine returns with a NULL indicating failure.
 *
 * @param  d365
 *         A pointer to struct _datetime. d365 contains the date to be translated.
 * @param  jd 
 *         A pointer to struct _julianday. The translated date is copied to JD.
 * @return jd
 *         A pointer to the filled jd. 
 */

struct _julianday*
getJulianFromDate365(struct _datetime *d365, struct _julianday* jd)
{
  /* Check for malformed parameters */
  if ((d365 != NULL) && (jd != NULL)){

  /* Julian is symmetric about 0.0 JD while d365 is not. 
   The code looks exactly the same in both if and else but it can be instructive to see how the two execution
   paths calculate the respective values. */
  if (d365->date.year < 0)
    {
      jd->day = d365->date.year * NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365 + nofDaysAfterARGMonthsInNonLeapYear[d365->date.month - 1] + d365->date.day - 1;
    }
  else
    {
      jd->day = d365->date.year * NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365 + nofDaysAfterARGMonthsInNonLeapYear[d365->date.month - 1] + d365->date.day - 1;
    }

  /* Handle the 12 hour skew between two time axis. */
  if (d365->time.hour < 12)
    {
      jd->ms = NO_OF_MS_IN_HALF_DAY + (NO_OF_MS_IN_A_HOUR * d365->time.hour);
      jd->day = jd->day - 1;
    }
  else
    {
      jd->ms = NO_OF_MS_IN_A_HOUR * (d365->time.hour - 12);
    }

  jd->ms = jd->ms + NO_OF_MS_IN_A_MINUTE * d365->time.minute + NO_OF_MS_IN_A_SECOND * d365->time.second + d365->time.ms;

  return jd;
}
else
  return NULL;
}
