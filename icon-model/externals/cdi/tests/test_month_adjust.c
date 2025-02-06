#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "normalize_month.h"

int
main(void)
{
  static const int vals[] = { -12, -17, 17, 11, -769652651 };
  enum
  {
    numVals = sizeof(vals) / sizeof(vals[0])
  };
  for (size_t i = 0; i < numVals; ++i)
    {
      enum
      {
        startYear = 1900
      };
      int month = vals[i];
      struct YearMonth ym = normalize_month(startYear, month);
      if ((long) ym.year * 12 + ym.month != (long) startYear * 12 + month)
        {
          fprintf(stderr,
                  "problem: month=%d, ym.month=%d, ym.year=%d\n"
                  "(long)ym.year * 12 + ym.month = %ld\n"
                  "(long)startYear * 12 + month = %ld\n",
                  month, ym.month, ym.year, (long) ym.year * 12 + ym.month, (long) startYear * 12 + month);
          abort();
        }
    }
  {
    struct timeval tv;
    if (gettimeofday(&tv, NULL))
      {
        perror("failed to get time for random number generator initialization");
        exit(EXIT_FAILURE);
      }
    srandom((unsigned) (tv.tv_sec ^ tv.tv_usec));
  }
  for (size_t j = 0; j < 1000000; ++j)
    {
      int year = (int) (random() - RAND_MAX / 2), month = (int) (random() - RAND_MAX / 2);
      struct YearMonth ym = normalize_month(year, month);
      if ((long) ym.year * 12 + ym.month != (long) year * 12 + month)
        {
          fprintf(stderr,
                  "problem: month=%d, ym.month=%d, ym.year=%d\n"
                  "(long)ym.year * 12 + ym.month = %ld\n"
                  "(long)year * 12 + month = %ld\n",
                  month, ym.month, ym.year, (long) ym.year * 12 + ym.month, (long) year * 12 + month);
          abort();
        }
    }
  return EXIT_SUCCESS;
}
