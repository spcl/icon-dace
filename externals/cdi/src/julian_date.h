#ifndef JULIAN_DATE_H
#define JULIAN_DATE_H

#include "cdi_datetime.h"
#include "calendar.h"

// clang-format off

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  int64_t julianDay;
  double secondOfDay;
} JulianDate;

JulianDate julianDate_encode(int calendar, CdiDateTime dt);
CdiDateTime julianDate_decode(int calendar, JulianDate julianDate);
JulianDate julianDate_add_seconds(JulianDate julianDate, int64_t seconds);
JulianDate julianDate_add(JulianDate julianDate1, JulianDate julianDate2);
JulianDate julianDate_sub(JulianDate julianDate1, JulianDate julianDate2);
double julianDate_to_seconds(JulianDate julianDate);

double secofday_encode(CdiTime time);
CdiTime secofday_decode(double secondOfDay);

#ifdef __cplusplus
}
#endif

// clang-format on

#endif /* JULIAN_DATE_H */
