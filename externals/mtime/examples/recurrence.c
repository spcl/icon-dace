// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/*
 * strings to parse:
 *
 * time intervals <interval>
 *
 * <interval> := <start>/<end>
 * <interval> := <start>/<duration>
 * <interval> := <duration>/<end>
 * <interval> := <duration>
 *
 * repeating time intervals
 *
 * R<n>/<interval>
 *
 * ommiting n means indefinit repetitions
 *
 *  1. split string by /
 *  2. assign components to sub parts
 *  3. need to scan the repititor with strtol with repetitor+1 as start
 *  4. get the others properly scaned by the ragel parser
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>

#include <mtime_calendar.h>
#include <mtime_datetime.h>
#include <mtime_timedelta.h>
#include <mtime_utilities.h>

void testCase(const char *caseX)
{
  char repetitor[MAX_REPETITION_STR_LEN];
  char start[MAX_REPETITION_STR_LEN];
  char end[MAX_REPETITION_STR_LEN];
  char duration[MAX_REPETITION_STR_LEN];
  
  struct _timedelta *d;
  struct _datetime *s, *e;
  
  char dstring[MAX_DATETIME_STR_LEN];  // Again, we use *STR_LEN.
  char tdstring[MAX_TIMEDELTA_STR_LEN]; // and here too.
    
  printf("Input %s\n", caseX);
  splitRepetitionString(caseX, repetitor, start, end, duration);
  if (repetitor[0] != '\0') 
    {
      printf("Repetitor: %s\n", repetitor);
      printf("Repetitions %d\n", getRepetitions(repetitor));
    }
  if (start[0] != '\0') 
    {
      printf("Start: %s\n", start);
      s = newDateTime(start);
      printf("Start rev: %s\n", datetimeToString(s, dstring));
      deallocateDateTime(s);
    }
  if (end[0] != '\0') 
    {
      printf("End: %s\n", end);
      e = newDateTime(end);
      printf("End rev: %s\n", datetimeToString(e, dstring));
      deallocateDateTime(e);
    }
  if (duration[0] != '\0') 
    {
      printf("Duration: %s\n", duration);
      d = newTimeDelta(duration);
      printf("Duration rev: %s\n", timedeltaToString(d, tdstring));
      deallocateTimeDelta(d);
    }
  printf("\n\n");
  
  return;
}

const char *cases[] = {
  "R2/20130304T000000/20130504T030405",
  "R3/8000101T230516/P10D",
  "R/PT2M/-34560101T040404",
  "R4/P1Y",
  "R012/2013-03-04T00:00:00.546/2013-05-04T03:04:05",
  "R3/800-01-01T23:05:16/P10D",
  "R/PT2M/-347856-01-01T04:04:04",
  "R4/P01Y",
  NULL
};

int main(int argc, char *argv[])
{
  initCalendar(PROLEPTIC_GREGORIAN);

  int i = 0;
  while (cases[i] != NULL)
    {
      testCase(cases[i]);
      i++;
    }

  return 0;
}
