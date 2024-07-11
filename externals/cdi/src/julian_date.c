#include "julian_date.h"
#include <math.h>

// convert Julian calendar day into year, months, day
static void
decode_julday(int calendar, int64_t julianDay,  // Julian day number to convert
              int *year,                        // Gregorian year (out)
              int *mon,                         // Gregorian month (1-12) (out)
              int *day)                         // Gregorian day (1-31) (out)
{
  int64_t a = julianDay;

  double b = floor((a - 1867216.25) / 36524.25);
  double c = a + b - floor(b / 4) + 1525;

  if (calendar == CALENDAR_STANDARD || calendar == CALENDAR_GREGORIAN)
    if (a < 2299161) c = a + 1524;

  double d = floor((c - 122.1) / 365.25);
  double e = floor(365.25 * d);
  double f = floor((c - e) / 30.6001);

  *day = (int) (c - e - floor(30.6001 * f));
  *mon = (int) (f - 1 - 12 * floor(f / 14));
  *year = (int) (d - 4715 - floor((7 + *mon) / 10));
}

// convert year, month, day into Julian calendar day
static int64_t
encode_julday(int calendar, int year, int month, int day)
{
  int iy = (month <= 2) ? year - 1 : year;
  int im = (month <= 2) ? month + 12 : month;
  int ib = (iy < 0) ? ((iy + 1) / 400 - (iy + 1) / 100) : (iy / 400 - iy / 100);

  if (calendar == CALENDAR_STANDARD || calendar == CALENDAR_GREGORIAN)
    {
      if (year > 1582 || (year == 1582 && (month > 10 || (month == 10 && day >= 15))))
        {
          // 15th October 1582 AD or later
        }
      else
        {
          // 4th October 1582 AD or earlier
          ib = -2;
        }
    }

  int64_t julianDay = (int64_t) (floor(365.25 * iy) + (int64_t) (30.6001 * (im + 1)) + ib + 1720996.5 + day + 0.5);

  return julianDay;
}

int64_t
date_to_julday(int calendar, int64_t date)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);

  return encode_julday(calendar, year, month, day);
}

int64_t
julday_to_date(int calendar, int64_t julianDay)
{
  int year, month, day;
  decode_julday(calendar, julianDay, &year, &month, &day);

  return cdiEncodeDate(year, month, day);
}

int
time_to_sec(int time)
{
  int hour, minute, second;
  cdiDecodeTime(time, &hour, &minute, &second);

  int seconds = hour * 3600 + minute * 60 + second;

  return seconds;
}

int
sec_to_time(int secofday)
{
  int hour = secofday / 3600;
  int minute = secofday / 60 - hour * 60;
  int second = secofday - hour * 3600 - minute * 60;

  return cdiEncodeTime(hour, minute, second);
}

double
secofday_encode(CdiTime time)
{
  int hour = time.hour;
  int minute = time.minute;
  int second = time.second;
  return hour * 3600 + minute * 60 + second + time.ms / 1000.0;
}

CdiTime
secofday_decode(double secondOfDay)
{
  CdiTime time;

  double secondOfDayIntegral;
  time.ms = lround(modf(secondOfDay, &secondOfDayIntegral) * 1000);

  int fullSeconds = lrint(secondOfDayIntegral);

  int hour = fullSeconds / 3600;
  int minute = fullSeconds / 60 - hour * 60;
  int second = fullSeconds - hour * 3600 - minute * 60;

  time.hour = hour;
  time.minute = minute;
  time.second = second;

  return time;
}

static int64_t
calendarDay_encode(int calendar, CdiDate date)
{
  int dpy = calendar_dpy(calendar);

  if (dpy == 360 || dpy == 365 || dpy == 366)
    return encode_calday(dpy, date.year, date.month, date.day);
  else
    return encode_julday(calendar, date.year, date.month, date.day);
}

static CdiDate
calendarDay_decode(int calendar, int64_t julday)
{
  int year, month, day;
  int dpy = calendar_dpy(calendar);

  if (dpy == 360 || dpy == 365 || dpy == 366)
    decode_calday(dpy, julday, &year, &month, &day);
  else
    decode_julday(calendar, julday, &year, &month, &day);

  CdiDate date;
  date.year = year;
  date.month = month;
  date.day = day;

  return date;
}

JulianDate
julianDate_encode(int calendar, CdiDateTime dt)
{
  JulianDate julianDate;

  julianDate.julianDay = calendarDay_encode(calendar, dt.date);
  julianDate.secondOfDay = secofday_encode(dt.time);

  return julianDate;
}

CdiDateTime
julianDate_decode(int calendar, JulianDate julianDate)
{
  CdiDateTime dt;

  dt.date = calendarDay_decode(calendar, julianDate.julianDay);
  dt.time = secofday_decode(julianDate.secondOfDay);

  return dt;
}

static void
adjust_seconds(JulianDate *julianDate)
{
  double SecondsPerDay = 86400.0;

  while (julianDate->secondOfDay >= SecondsPerDay)
    {
      julianDate->secondOfDay -= SecondsPerDay;
      julianDate->julianDay++;
    }

  while (julianDate->secondOfDay < 0.0)
    {
      julianDate->secondOfDay += SecondsPerDay;
      julianDate->julianDay--;
    }
}

// add seconds to julianDate
JulianDate
julianDate_add_seconds(JulianDate julianDate, int64_t seconds)
{
  julianDate.secondOfDay += seconds;

  adjust_seconds(&julianDate);

  return julianDate;
}

// add julianDate1 and julianDate2
JulianDate
julianDate_add(JulianDate julianDate1, JulianDate julianDate2)
{
  JulianDate julianDate;
  julianDate.julianDay = julianDate1.julianDay + julianDate2.julianDay;
  julianDate.secondOfDay = julianDate1.secondOfDay + julianDate2.secondOfDay;

  adjust_seconds(&julianDate);

  return julianDate;
}

// subtract julianDate2 from julianDate1
JulianDate
julianDate_sub(JulianDate julianDate1, JulianDate julianDate2)
{
  JulianDate julianDate;
  julianDate.julianDay = julianDate1.julianDay - julianDate2.julianDay;
  julianDate.secondOfDay = julianDate1.secondOfDay - julianDate2.secondOfDay;

  adjust_seconds(&julianDate);

  return julianDate;
}

double
julianDate_to_seconds(JulianDate julianDate)
{
  return julianDate.julianDay * 86400.0 + julianDate.secondOfDay;
}

#ifdef TEST2
int
main(void)
{
  int calendar = CALENDAR_STANDARD;
  int factor = 86400;
  int value = 30;

  int year = 1979;
  int month = 1;
  int day = 15;
  int hour = 12;
  int minute = 30;
  int second = 17;
  int ms = 0;

  CdiDateTime dt;
  dt.date = cdiDate_encode(year, month, day);
  dt.time = cdiTime_encode(hour, minute, second, ms);
  printf("%d/%02d/%02d %02d:%02d:%02d.%03d\n", dt.date.year, dt.date.month, dt.date.day, dt.time.hour, dt.time.minute,
         dt.time.second, dt.time.ms);

  JulianDate julianDate = julianDate_encode(calendar, dt);

  dt = julianDate_decode(calendar, julianDate);
  printf("%d/%02d/%02d %02d:%02d:%02d.%03d   %d %g\n", dt.date.year, dt.date.month, dt.date.day, dt.time.hour, dt.time.minute,
         dt.time.second, dt.time.ms, (int) julianDate.julianDay, julianDate.secondOfDay);

  for (int i = 0; i < 420; i++)
    {
      dt = julianDate_decode(calendar, julianDate);
      printf("%2d %d/%02d/%02d %02d:%02d:%02d.%03d\n", i, dt.date.year, dt.date.month, dt.date.day, dt.time.hour, dt.time.minute,
             dt.time.second, dt.time.ms);
      julianDate = julianDate_add_seconds(julianDate, value * factor);
    }

  return 0;
}
#endif
