//#define YAC_VERBOSE_EVENT
/**
 * @file event.c
 *
 * @copyright Copyright  (C)  2013 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "event.h"
#include "utils.h"

#define SECONDS_PER_DAY 86400

int yac_number_of_events;

struct event {
  int is_set;
  struct _timedelta *model_time_step;
  struct _timedelta *coupling_time_step;
  struct _timedelta *dummy_time_step;
  int lag;
  enum yac_reduction_type time_operation;
  enum yac_action_type action;
  enum yac_time_unit_type time_unit;
  struct _datetime *startdate;
  struct _datetime *stopdate;
  struct _datetime *currentdate;
  struct _datetime *offset;
  struct _event* coupling;
};

/* ---------------------------------------------------------------------- */

struct event * yac_event_new () {

  /* allocate memory for one event */
  struct event * event = xmalloc(1 * sizeof (*event));

  /* initialise the event struct */

  event->is_set         = 0;
  event->lag            = 0;

  event->time_operation = TIME_NONE;

  ++yac_number_of_events;

  return event;
}

/* ---------------------------------------------------------------------- */

void yac_event_add ( struct event * event,
                     char const * delta_model_time,
                     char const * delta_coupling_time,
                     int lag,
                     enum yac_reduction_type time_operation,
                     const char * startdate,
                     const char * stopdate) {

  YAC_ASSERT(
    getCalendarType() != CALENDAR_NOT_SET,
    "ERROR(yac_event_add): no calendar has been defined")
  YAC_ASSERT(!event->is_set,
    "ERROR(yac_event_add): event is already initialised.");

  /* set the event triggering parameters */

  char PTstring[MAX_TIMEDELTA_STR_LEN];

  event->is_set         = 1;
  event->lag            = lag;
  event->time_operation = time_operation;

  /* create dummy_event */

  event->dummy_time_step = newTimeDelta(getPTStringFromMS(0, PTstring));

  /* model time step and coupling time step need to be converted */
  event->model_time_step = newTimeDelta(delta_model_time);
  event->coupling_time_step = newTimeDelta(delta_coupling_time);

  YAC_ASSERT_F(
    event->model_time_step,
    "ERROR(yac_event_add): failed to generate model timestep from \"%s\" "
    "(has to be in ISO 8601 format)", delta_model_time)
  YAC_ASSERT_F(
    event->coupling_time_step,
    "ERROR(yac_event_add): failed to generate model timestep from \"%s\""
    "(has to be in ISO 8601 format)", delta_coupling_time)

  // if the model time step and the coupling time step are identical, a time
  // operation does not make sense
  if (compareTimeDelta(event->model_time_step, event->coupling_time_step)
      == equal_to) {
    event->time_operation = TIME_NONE;
  }

  calendarType set_calender = getCalendarType ();

  YAC_ASSERT(
    set_calender != CALENDAR_NOT_SET, "ERROR: Calender is not set.")

  event->startdate = newDateTime(startdate);
  YAC_ASSERT(event->startdate != NULL, "ERROR in event definition")

  event->currentdate = newDateTime(startdate);
  YAC_ASSERT(event->currentdate != NULL, "ERROR in event definition")

  event->stopdate = newDateTime(stopdate);
  YAC_ASSERT(event->stopdate != NULL, "ERROR in event definition")

  char timedelta_str[MAX_TIMEDELTA_STR_LEN];

  event->coupling =
    newEvent(
      "Coupling", (char *) startdate, (char *) startdate, (char *) stopdate,
      timedeltaToString(event->coupling_time_step, timedelta_str), NULL );

  YAC_ASSERT(event->coupling != NULL, "ERROR in event definition")

#ifdef YAC_VERBOSE_EVENT
  char datetime_str[MAX_DATETIME_STR_LEN];
  printf ("event add for event ID %d \n", event_id );
  printf ("- start   date:    %s\n",
          datetimeToString(event->startdate,datetime_str) );
  printf ("- current date:    %s\n",
          datetimeToString(event->currentdate,datetime_str) );
  printf ("- stop    date:    %s\n",
          datetimeToString(event->stopdate,datetime_str) );
  printf ("- model    time step: %s\n",
          timedeltaToString(event->model_time_step,timedelta_str) );
  printf ("- coupling time step: %s\n",
          timedeltaToString(event->coupling_time_step,timedelta_str) );
#endif

  YAC_ASSERT(event->lag >= 0, "ERROR: lag has to be greater than zero")

  /* add lag times delta(t) to the current date such that we can detect
     a restart event which has to be performed lag times delta(t) before
     the final date is reached. */
  if ( event->lag > 0 ) {
    for ( int i = 0; i < event->lag; ++i )
      event->currentdate =
        addTimeDeltaToDateTime (
          event->currentdate, event->model_time_step, event->currentdate );
#ifdef YAC_VERBOSE_EVENT
    printf ("- model   time lag:  %d\n", event->lag );
    printf ("- new current date:  %s\n",
            datetimeToString(event->currentdate,datetime_str) );
#endif
  }
}

/* ---------------------------------------------------------------------- */

void yac_event_update ( struct event * event ) {

  YAC_ASSERT(event->is_set,
    "ERROR(yac_event_update): event has not yet been initialised");

  event->currentdate =
    addTimeDeltaToDateTime(
      event->currentdate, event->model_time_step, event->currentdate);
}

/* ---------------------------------------------------------------------- */

enum yac_action_type yac_event_check ( struct event * event ) {

  /* Return value action:
     --------------------

     action = 0  : no action to be performed
     action = 1  : time operation required (only put)
     action = 2  : coupling action required
     action = 3  : restart file needs to be written
     action = 4-7: last put/get for restart
     action = 8  : date is out of bound, beyond restart

     -------------------------------------------------------------------- */

#ifdef YAC_VERBOSE_EVENT
  char datetime_str[MAX_DATETIME_STR_LEN];
  printf ("- current date: %s\n",
          datetimeToString(event->currentdate,datetime_str) );
  printf ("- stop    date: %s\n",
          datetimeToString(event->stopdate,datetime_str) );
#endif

  compare_return_val check_stopdate =
    compareDatetime(event->currentdate,event->stopdate);

  YAC_ASSERT(
    (check_stopdate == less_than) ||
    (check_stopdate == equal_to) ||
    (check_stopdate == greater_than),
    "ERROR(yac_event_check): compareDatetime(currentdate,stopdate) failed")

  enum yac_action_type action;

  /* Detect out of bound */
  if (check_stopdate == greater_than) {

#ifdef YAC_VERBOSE_EVENT
     printf ("event check for event ID %d delivers OUT_OF_BOUND\n", event_id );
#endif
     action = OUT_OF_BOUND;

  } else {

    compare_return_val  check_startdate =
      compareDatetime(event->currentdate,event->startdate);
    YAC_ASSERT(
      check_startdate != less_than,
      "ERROR(yac_event_check): compareDatetime(currentdate,startdate) failed; "
      "current date is before start date")
    YAC_ASSERT(
      (check_startdate == equal_to) ||
      (check_startdate == greater_than),
      "ERROR(yac_event_check): compareDatetime(currentdate,startdate) failed")

    /* Detect active event */
    if (isCurrentEventActive(event->coupling,
                              event->currentdate,
                              event->dummy_time_step,
                              event->dummy_time_step)) {

      /* Detect restart event */
      if (check_stopdate == equal_to) {

#ifdef YAC_VERBOSE_EVENT
        printf ("event check for event ID %d delivers RESTART\n", event_id );
#endif
        action = RESTART;

      } else {

#ifdef YAC_VERBOSE_EVENT
        printf ("event check for event ID %d delivers COUPLING\n", event_id );
#endif
        action = COUPLING;
      }
    } else {

      action = (event->time_operation == TIME_NONE)?NONE:REDUCTION;
    }
  }

  event->action = action;
  return action;
}

/* ---------------------------------------------------------------------- */

void yac_event_delete ( struct event * event ) {

  deallocateDateTime(event->startdate);
  deallocateDateTime(event->currentdate);
  deallocateDateTime(event->stopdate);

  deallocateTimeDelta(event->dummy_time_step);
  deallocateTimeDelta(event->model_time_step);
  deallocateTimeDelta(event->coupling_time_step);

  deallocateEvent(event->coupling);

  free(event);

  --yac_number_of_events;
}

/* ---------------------------------------------------------------------- */

int yac_get_event_lag ( struct event * event ) {

  return event->lag;
}

/* ---------------------------------------------------------------------- */

int yac_get_event_time_operation ( struct event * event ) {

  return event->time_operation;
}

/* ---------------------------------------------------------------------- */

char * yac_get_event_coupling_timestep (
  struct event * event, char * timedelta_str ) {

  memset(timedelta_str,'\0',MAX_TIMEDELTA_STR_LEN);

  timedeltaToString(event->coupling_time_step,timedelta_str);

  return timedelta_str;
}


/* ---------------------------------------------------------------------- */

char * yac_get_event_model_timestep (
  struct event * event, char * timedelta_str ) {

  memset(timedelta_str,'\0',MAX_TIMEDELTA_STR_LEN);

  timedeltaToString(event->model_time_step,timedelta_str);

  return timedelta_str;
}

char * yac_get_event_current_datetime (
  struct event * event, char* datetime_str) {
  datetimeToString(event->currentdate, datetime_str);
  return datetime_str;
}

static int64_t str2int64(char const * time) {

  char * endptr;
  int64_t time_int64_t = (int64_t)strtol(time, &endptr, 10);
  YAC_ASSERT_F(
    *endptr == '\0', "ERROR(str2int64): invalid time string '%s'", time)
  return time_int64_t;
}

char const * yac_time_to_ISO(
  char const * time, enum yac_time_unit_type time_unit) {

  YAC_ASSERT(
    (time_unit == C_MILLISECOND) ||
    (time_unit == C_SECOND) ||
    (time_unit == C_MINUTE) ||
    (time_unit == C_HOUR) ||
    (time_unit == C_ISO_FORMAT),
    "ERROR(yac_time_to_ISO): unsupported time unit")

  YAC_ASSERT(
    getCalendarType() != CALENDAR_NOT_SET,
    "ERROR(yac_time_to_ISO): calendar has not yet been set")

  static char PT_String_buffer[MAX_TIMEDELTA_STR_LEN];

  switch (time_unit) {

    default:
    case (C_MILLISECOND):
      getPTStringFromMS(str2int64(time), PT_String_buffer);
      break;
    case (C_SECOND):
      getPTStringFromSeconds(str2int64(time), PT_String_buffer);
      break;
    case (C_MINUTE):
      getPTStringFromMinutes(str2int64(time), PT_String_buffer);
      break;
    case (C_HOUR):
      getPTStringFromHours(str2int64(time), PT_String_buffer);
      break;
    case (C_ISO_FORMAT):
      strncpy(PT_String_buffer, time, MAX_TIMEDELTA_STR_LEN-1);
      break;
  };

  return PT_String_buffer;
}
