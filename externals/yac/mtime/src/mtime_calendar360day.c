/**
 * @file mtime_calendar360day.c
 * 
 * @brief Convert julian dates to Calendar-with-360-days dates and vice-versa.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note Calendar-with-360-days has 30 days in 'every' month and hence 360 days in every calendar year. 
 *       Also, Julian Day (0,0) corresponds to Calendar-with-360-day's day 0-01-01T12:00:00.000.
 */

#include <stdint.h>
#include <stddef.h>

#include "mtime_calendar360day.h"
#include "mtime_julianDay.h"
#include "mtime_datetime.h"


/**
 * @brief Convert Julian Date to Calendar-with-360-days date.
 *
 * The getDate360FromJulian routine accepts a date on Julian axis and translates it to a date
 * on Calendar-of-360-days. For eg. If Julian Day values are day = 1 and ms = 0, the corresponding 
 * Calendar-of-360-days values will be Year = 0, Month = 1, Day = 2, Hour = 12, Minute = 0, 
 * Second = 0 and MS = 0.
 *
 * The routine expects non NULL parameters. In case Julian-Day or DateTime
 * is NULL, the routine returns with a NULL indicating failure.
 *
 * @param  jd 
 *         A pointer to struct _julianday. JD contains the date to be 'translated' to type d360.
 * @param  d360
 *	  A pointer to struct _datetime. The translated date is copied into d360.
 * @return d360
 * 	  A pointer to the filled d360.
 */

struct _datetime*
getDate360FromJulian(struct _julianday *jd, struct _datetime* d360)
{
  /* Check for malformed parameters */
  if ((jd != NULL) && (d360 != NULL)){

  int days;
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
      d360->date.year = (jday) / NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360;

      days = (int)(jday % NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360);
      d360->date.month = (int)(days / NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360) + 1;
      d360->date.day = (int)(days % NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360) + 1;
    }
  else
    {
      d360->date.year = (jday / NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360) - 1;

      days = (int)(-jday % NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360);
      d360->date.month = NO_OF_MONTHS_IN_A_YEAR - (days / NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360);
      d360->date.day = NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 - (days % NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 - 1);

      /* This will be true whenever days is divisible by NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 */
      if (d360->date.day == (NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 + 1))
        {
          d360->date.month = d360->date.month + 1;
          if (d360->date.month == (NO_OF_MONTHS_IN_A_YEAR + 1))
            {
              d360->date.year = d360->date.year + 1;
              d360->date.month = 1;
            }
          d360->date.day = 1;
        }
    }

  /* Santize for range. */
  if((d360->date.year > YEAR_UPPER_BOUND) || (d360->date.year < YEAR_LOWER_BOUND))
    {
      /* ERROR: Exceeds allowed year range. */
      return NULL;
    }

  d360->time.hour = (int) jms / NO_OF_MS_IN_A_HOUR;
  d360->time.minute = ((int) jms - d360->time.hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
  d360->time.second = ((int) jms - d360->time.hour * NO_OF_MS_IN_A_HOUR - d360->time.minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
  d360->time.ms = (int) jms - d360->time.hour * NO_OF_MS_IN_A_HOUR - d360->time.minute * NO_OF_MS_IN_A_MINUTE - d360->time.second * NO_OF_MS_IN_A_SECOND;

  return d360;
}
else
  return NULL;
}

  /**
   * @brief Convert Calendar-with-360-days date to Julian Date.
   *
   * The getJulianFromDate360 routine accepts a date on Calendar-of-360-days axis and translates it to a date
   * on Julian axis. For eg. If Calendar-of-360-days values are Year = 0, Month = 1, Day = 2, Hour = 12, Minute = 0, 
   * Second = 0 and MS = 0 , corresponding Julian Day values will be day = 1 and ms = 0.
   *
   * The routine expects non NULL parameters. In case Julian-Day or DateTime
   * is NULL, the routine returns with a NULL indicating failure.
   *
   * @param  d360
   *         A pointer to struct _datetime. d360 contains the date to be translated.
   * @param  jd 
   *         A pointer to struct _julianday. The translated date is copied to JD.
   * @return jd
   *         A pointer to the filled jd. 
   */

struct _julianday*
getJulianFromDate360(struct _datetime *d360, struct _julianday *jd)
{
  /* Check for malformed parameters */
  if ((d360 != NULL) && (jd != NULL)){

  /* Julian is symmetric about 0.0 JD while d360 is not. */
  if (d360->date.year < 0)
    {
      jd->day =    d360->date.year * NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360 
		+ (d360->date.month) * NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 
		+ (d360->date.day - NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360) - 1;
    }
  else
    {
      jd->day =    d360->date.year * NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360 
		+ (d360->date.month - 1) * NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 
		+ (d360->date.day - 1);
    }

  /* Handle the 12 hour skew between two time axis. */
  if (d360->time.hour < 12)
    {
      jd->ms 	= NO_OF_MS_IN_HALF_DAY + (NO_OF_MS_IN_A_HOUR * d360->time.hour);
      jd->day 	= jd->day - 1;
    }
  else
    {
      jd->ms = NO_OF_MS_IN_A_HOUR * (d360->time.hour - 12);
    }

  jd->ms = jd->ms + NO_OF_MS_IN_A_MINUTE * d360->time.minute 
		  + NO_OF_MS_IN_A_SECOND * d360->time.second 
		  + d360->time.ms;

  return jd;
}
else
  return NULL;
}
