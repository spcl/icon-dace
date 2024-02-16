#ifndef NORMALIZE_MONTH_H
#define NORMALIZE_MONTH_H

#include <stdlib.h>

struct YearMonth
{
  int year, month;
};

/* normalizes month to range [1,12] and adjusts year accordingly */
static inline struct YearMonth
normalize_month(int year, int month)
{
  div_t modres = div(month - 1, 12);
  year += modres.quot - ((month < 1) & (modres.rem != 0));
  return (struct YearMonth){ .year = year, .month = (modres.rem + 12) % 12 + 1 };
}

#endif
