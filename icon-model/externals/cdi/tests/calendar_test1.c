#include <stdio.h>
#include <stdlib.h>

#include "cdi.h"
#include "julian_date.h"

int
main(void)
{
  int calendar = CALENDAR_STANDARD;
  int j = 0;

  // 1 - Check valid range of years
  {
    int nmin = 11000;
    int vdate0 = -80001201;
    int vtime0 = 120500;

    // printf("start time: %8d %4d\n", vdate0, vtime0);

    CdiDateTime dt = cdiDateTime_set(vdate0, vtime0);

    for (int i = 0; i < nmin; i++)
      {
        JulianDate julianDate = julianDate_encode(calendar, dt);

        dt = julianDate_decode(calendar, julianDate);
        int vdate = (int) cdiDate_get(dt.date);
        int vtime = cdiTime_get(dt.time);

        if (vdate0 != vdate || vtime0 != vtime)
          fprintf(stderr, "%4d %8d %4d %8d %4d %9d %g\n", ++j, vdate0, vtime0, vdate, vtime, (int) julianDate.julianDay,
                  julianDate.secondOfDay);

        dt.date.year++;
        julianDate = julianDate_encode(calendar, dt);

        dt = julianDate_decode(calendar, julianDate);
        vdate0 = (int) cdiDate_get(dt.date);
        vtime0 = cdiTime_get(dt.time);
      }

    // printf("stop time: %8d %4d\n", vdate0, vtime0);
  }
  // 2 - Check time increment of one minute
  {
    int nmin = 120000;
    int ijulinc = 60;
    int vdate0 = 20001201;
    int vtime0 = 0;

    // printf("start time: %8d %4d\n", vdate0, vtime0);

    CdiDateTime dt = cdiDateTime_set(vdate0, vtime0);
    JulianDate julianDate = julianDate_encode(calendar, dt);

    for (int i = 0; i < nmin; i++)
      {
        int year, mon, day, hour, minute, second;
        cdiDecodeDate(vdate0, &year, &mon, &day);
        cdiDecodeTime(vtime0, &hour, &minute, &second);

        if (++minute >= 60)
          {
            minute = 0;
            if (++hour >= 24)
              {
                hour = 0;
                if (++day >= 32)
                  {
                    day = 1;
                    if (++mon >= 13)
                      {
                        mon = 1;
                        year++;
                      }
                  }
              }
          }

        vdate0 = cdiEncodeDate(year, mon, day);
        vtime0 = cdiEncodeTime(hour, minute, second);

        julianDate = julianDate_add_seconds(julianDate, ijulinc);
        dt = julianDate_decode(calendar, julianDate);
        int vdate = (int) cdiDate_get(dt.date);
        int vtime = cdiTime_get(dt.time);

        if (vdate0 != vdate || vtime0 != vtime)
          fprintf(stderr, "%4d %8d %4d %8d %4d %9d %g\n", ++j, vdate0, vtime0, vdate, vtime, (int) julianDate.julianDay,
                  julianDate.secondOfDay);
      }

    // printf("stop time: %8d %4d\n", vdate0, vtime0);
  }
  return j == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
