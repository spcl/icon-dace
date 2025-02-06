// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @file mtime_calendarGregorian.c
 *
 * @brief Convert julian dates to Gregorian-Calendar dates and vice-versa.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note TODO: LUIS, you might want to add the details here. And may be add some comments in the code?
 */

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

#include "mtime_calendarGregorian.h"
#include "mtime_julianDay.h"
#include "mtime_datetime.h"

/* Number of days after 0,1,2...12 months. */
static const int msn[12] = { 0, 31, 61, 92, 122, 153, 184, 214, 245, 275, 306, 337 };

/* Compute floor division */
static int64_t
ifloor(int64_t ir1, int64_t ir2)
{
  lldiv_t r = lldiv(ir1, ir2);
  return r.quot - (r.rem < 0);
}

/**
 * @brief Convert Julian Date to Gregorian date.
 *
 * The getDateGregorianFromJulian routine accepts a date on Julian axis and translates it to a date
 * on Gregorian. For eg. If Julian Day values are day = 1721060 and ms = 0, the corresponding
 * Gregorian values will be Year = 0, Month = 1, Day = 1, Hour = 12, Minute = 0,
 * Second = 0 and MS = 0.
 *
 * The routine expects non NULL parameters. In case Julian-Day or DateTime
 * is NULL, the routine returns with a NULL indicating failure.
 *
 * @param  jd
 *         A pointer to struct _julianday. JD contains the date to be 'translated' to type gregorian.
 * @param  gd
 *         A pointer to struct _datetime. The translated date is copied into gd.
 * @return gd
 *         A pointer to the filled gd.
 */

struct _datetime *
getDateGregorianFromJulian(struct _julianday *jd, struct _datetime *gd)
{
  if ((jd != NULL) && (gd != NULL))
    {
      int64_t tmp2, tmp3;

      int64_t tjd, tms;
      int64_t a, b, c, d, j, m, y;

      tjd = jd->day;
      tms = jd->ms;
      if (tms >= 43200000)
        {
          tjd = tjd + 1;
          tms = tms - 86400000;
        }
      j = tjd - 1721119;
      d = 100 * j - 25;
      tmp2 = 3652425;
      a = ifloor(d, tmp2);
      tmp2 = 4;
      b = a - ifloor(a, tmp2);
      tmp2 = 100 * b + d;
      tmp3 = 36525;
      y = ifloor(tmp2, tmp3);
      tmp2 = 4;
      c = b + j - 365 * y - ifloor(y, tmp2);
      m = (5 * c + 456) / 153;
      d = c - msn[m - 3];
      if (m > 12)
        {
          y = y + 1;
          m = m - 12;
        }

      gd->date.year = y;
      // TODO on Luis: Is this ok?
      if ((gd->date.year > YEAR_UPPER_BOUND) || (gd->date.year < YEAR_LOWER_BOUND))
        {
          /* ERROR: Exceeds allowed year range. */
          return NULL;
        }

      gd->date.month = (int) m;
      gd->date.day = (int) d;
      gd->time.hour = (int) ifloor(tms, 3600000);
      gd->time.minute = ((int) tms - gd->time.hour * 3600000) / 60000;
      gd->time.second = ((int) tms - gd->time.hour * 3600000 - gd->time.minute * 60000) / 1000;
      gd->time.ms = (int) tms - gd->time.hour * 3600000 - gd->time.minute * 60000 - gd->time.second * 1000;
      gd->time.hour = gd->time.hour + 12;

      return gd;
    }
  else
    return NULL;
}

/**
 * @brief Convert Gregorian date to Julian Date.
 *
 * The getJulianFromDateGregorian routine accepts a date on Gregorian axis and translates it to a date
 * on Julian axis. For eg. If Gregorian values are Year = -4713, Month = 11, Day = 24, Hour = 12, Minute = 0,
 * Second = 0 and MS = 0 , corresponding Julian Day values will be day = 0 and ms = 0.
 *
 * The routine expects non NULL parameters. In case Julian-Day or DateTime
 * is NULL, the routine returns with a NULL indicating failure.
 *
 * @param  gd
 *         A pointer to struct _datetime. gd contains the date to be translated.
 * @param  jd
 *         A pointer to struct _julianday. The translated date is copied to JD.
 * @return jd
 *         A pointer to the filled jd.
 */

struct _julianday *
getJulianFromDateGregorian(struct _datetime *gd, struct _julianday *jd)
{
  if ((gd != NULL) && (jd != NULL))
    {
      int64_t tmp2, tmp3, tmp4;

      int64_t d, m, y;

      y = gd->date.year;
      m = gd->date.month;
      d = gd->date.day;

      if (m < 3)
        {
          m = m + 12;
          y = y - 1;
        }
      tmp2 = 4;
      tmp3 = 100;
      tmp4 = 400;
      jd->day = d + msn[m - 3] + 365 * y + ifloor(y, tmp2) - ifloor(y, tmp3) + ifloor(y, tmp4) + 1721119;
      jd->ms = 0;
      jd->ms = jd->ms + 3600000 * (gd->time.hour - 12);
      jd->ms = jd->ms + 60000 * gd->time.minute;
      jd->ms = jd->ms + 1000 * gd->time.second;
      jd->ms = jd->ms + gd->time.ms;
      if (jd->ms < 0)
        {
          jd->ms = 86400000 + jd->ms;
          jd->day = jd->day - 1;
        }

      return jd;
    }
  else
    return NULL;
}
