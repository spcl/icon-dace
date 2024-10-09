// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include <stdio.h>
#include <stdlib.h>

#include <mtime_calendar.h>
#include <mtime_timedelta.h>

int
main(void)
{
  initCalendar(PROLEPTIC_GREGORIAN);

  {
    struct _timedelta *td1 = newTimeDelta("PT10M");
    struct _timedelta *td2 = newTimeDelta("PT2H");

    int64_t remainder, quot;

    remainder = moduloTimedelta(td2, td1, &quot);

    fprintf(stderr, "quot = %ld, rem = %ld\n", quot, remainder);
  }

  return 0;
}
