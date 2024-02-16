//#define VERBOSE
/**
 * @file test_events.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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
#include "tests.h"
#include "event.h"
#include "utils.h"

int main (void) {

  // basic setup
  initCalendar(PROLEPTIC_GREGORIAN);

  {
    struct {
      enum yac_time_unit_type time_unit;
      char const * value;
    } tests[] =
      {{.time_unit = C_MILLISECOND,
        .value = "3600000"},
      {.time_unit = C_SECOND,
        .value = "3600"},
        {.time_unit = C_MINUTE,
        .value = "60"},
        {.time_unit = C_HOUR,
        .value = "1"},
        {.time_unit = C_ISO_FORMAT,
        .value = "PT01H"}};
    enum {NUM_TESTS = sizeof(tests) / sizeof(tests[0])};
    for (size_t i = 0; i < NUM_TESTS; ++i)
      if (strcmp(
            yac_time_to_ISO(
              tests[i].value, tests[i].time_unit), "PT01H"))
        PUT_ERR("ERROR in yac_time_to_ISO");
  }

  {
    struct event * event = yac_event_new ();

    char const * model_time_step    = "5000";  // 5s per time step
    char const * coupling_time_step = "15000"; // 15s per coupling time step

    enum yac_time_unit_type time_unit = C_MILLISECOND;
    int lag = 0;
    enum yac_reduction_type time_operation = TIME_AVERAGE;

    // 60s runtime
    char const * start_date = "1850-01-01T00:00:00";
    char const * stop_date =  "1850-01-01T00:01:00";

    char const * model_time_step_iso =
      strdup(yac_time_to_ISO(model_time_step, time_unit));
    char const * coupling_time_step_iso =
      strdup(yac_time_to_ISO(coupling_time_step, time_unit));
    yac_event_add(
      event, model_time_step_iso, coupling_time_step_iso,
      lag, time_operation, start_date, stop_date);

    char time_step_buffer[MAX_TIMEDELTA_STR_LEN];
    if (strcmp(
          yac_get_event_coupling_timestep(event, time_step_buffer),
          coupling_time_step_iso))
      PUT_ERR("ERROR in yac_get_event_coupling_timestep");
    if (strcmp(
          yac_get_event_model_timestep(event, time_step_buffer),
          model_time_step_iso))
      PUT_ERR("ERROR in yac_get_event_model_timestep");
    if (yac_get_event_lag(event) != lag)
      PUT_ERR("ERROR in yac_get_event_lag");

    free((void*)coupling_time_step_iso);
    free((void*)model_time_step_iso);

    enum {NUM_TIME_STEPS = 15};

    enum yac_action_type ref_action[NUM_TIME_STEPS] =
      {COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       RESTART, OUT_OF_BOUND, OUT_OF_BOUND};

    for (int t = 0; t < NUM_TIME_STEPS; ++t) {

      if (yac_event_check(event) != ref_action[t]) PUT_ERR("wrong action\n");
      yac_event_update(event);
    }

    yac_event_delete(event);
  }

  {
    struct event * event = yac_event_new ();

    char const * model_time_step    = "5000";  // 5s per time step
    char const * coupling_time_step = "15000"; // 15s per coupling time step

    enum yac_time_unit_type time_unit = C_MILLISECOND;
    int lag = 0;
    enum yac_reduction_type time_operation = TIME_NONE;

    // 60s runtime
    char const * start_date = "1850-01-01T00:00:00";
    char const * stop_date =  "1850-01-01T00:01:00";

    char const * model_time_step_iso =
      strdup(yac_time_to_ISO(model_time_step, time_unit));
    char const * coupling_time_step_iso =
      strdup(yac_time_to_ISO(coupling_time_step, time_unit));
    yac_event_add(
      event, model_time_step_iso, coupling_time_step_iso,
      lag, time_operation, start_date, stop_date);

    char time_step_buffer[MAX_TIMEDELTA_STR_LEN];
    if (strcmp(
          yac_get_event_coupling_timestep(event, time_step_buffer),
          coupling_time_step_iso))
      PUT_ERR("ERROR in yac_get_event_coupling_timestep");
    if (strcmp(
          yac_get_event_model_timestep(event, time_step_buffer),
          model_time_step_iso))
      PUT_ERR("ERROR in yac_get_event_model_timestep");
    if (yac_get_event_lag(event) != lag)
      PUT_ERR("ERROR in yac_get_event_lag");

    free((void*)coupling_time_step_iso);
    free((void*)model_time_step_iso);

    enum {NUM_TIME_STEPS = 15};

    enum yac_action_type ref_action[NUM_TIME_STEPS] =
      {COUPLING, NONE, NONE,
       COUPLING, NONE, NONE,
       COUPLING, NONE, NONE,
       COUPLING, NONE, NONE,
       RESTART, OUT_OF_BOUND, OUT_OF_BOUND};

    for (int t = 0; t < NUM_TIME_STEPS; ++t) {

      if (yac_event_check(event) != ref_action[t]) PUT_ERR("wrong action\n");
      yac_event_update(event);
    }

    yac_event_delete(event);
  }

  {
    struct event * event = yac_event_new ();

    char * model_time_step = "PT30M";
    char * coupling_time_step = "PT60M";

    enum yac_time_unit_type time_unit = C_ISO_FORMAT;
    int lag = 0;
    enum yac_reduction_type time_operation = TIME_AVERAGE;

    // 60s runtime
    char const * start_date = "1850-01-01T00:00:00";
    char const * stop_date =  "1850-01-01T06:00:00";

    char const * model_time_step_iso =
      strdup(yac_time_to_ISO(model_time_step, time_unit));
    char const * coupling_time_step_iso =
      strdup(yac_time_to_ISO(coupling_time_step, time_unit));
    yac_event_add(
      event, model_time_step_iso, coupling_time_step_iso,
      lag, time_operation, start_date, stop_date);

    char time_step_buffer[MAX_TIMEDELTA_STR_LEN];
    if (strcmp(
          yac_get_event_coupling_timestep(event, time_step_buffer),
          coupling_time_step_iso))
      PUT_ERR("ERROR in yac_get_event_coupling_timestep");
    if (strcmp(
          yac_get_event_model_timestep(event, time_step_buffer),
          model_time_step_iso))
      PUT_ERR("ERROR in yac_get_event_model_timestep");
    if (yac_get_event_lag(event) != lag)
      PUT_ERR("ERROR in yac_get_event_lag");

    free((void*)coupling_time_step_iso);
    free((void*)model_time_step_iso);

    enum {NUM_TIME_STEPS = 16};

    enum yac_action_type ref_action[NUM_TIME_STEPS] =
      {COUPLING, REDUCTION, COUPLING, REDUCTION,
       COUPLING, REDUCTION, COUPLING, REDUCTION,
       COUPLING, REDUCTION, COUPLING, REDUCTION,
       RESTART, OUT_OF_BOUND, OUT_OF_BOUND, OUT_OF_BOUND};

    for (int t = 0; t < NUM_TIME_STEPS; ++t) {

      if (yac_event_check(event) != ref_action[t]) PUT_ERR("wrong action\n");
      yac_event_update(event);
    }

    yac_event_delete(event);
  }

  freeCalendar();

  return TEST_EXIT_CODE;
}


