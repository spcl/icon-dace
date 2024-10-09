// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tests.h"
#include "event.h"

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

  {
    char const * model_time_step    = "5000";  // 5s per time step
    char const * coupling_time_step = "15000"; // 15s per coupling time step

    enum yac_time_unit_type time_unit = C_MILLISECOND;
    enum yac_reduction_type time_operation = TIME_AVERAGE;

    // 60s runtime
    char const * start_date = "1850-01-01T00:00:00";
    char const * stop_date =  "1850-01-01T00:01:00";

    enum {MIN_LAG = -6, MAX_LAG = 3};
    enum {NUM_LAGS = MAX_LAG - MIN_LAG + 1};
    enum {NUM_TIME_STEPS = 15};

    enum yac_action_type ref_action[NUM_TIME_STEPS - MIN_LAG] =
      {COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       COUPLING, REDUCTION, REDUCTION,
       RESTART, OUT_OF_BOUND, OUT_OF_BOUND};

    for (int lag = MIN_LAG; lag <= MAX_LAG; ++lag) {

      struct event * event = yac_event_new ();

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

      for (int t = lag; t < NUM_TIME_STEPS; ++t) {

        enum yac_action_type action = yac_event_check(event);
        if (action != ref_action[t - MIN_LAG])
          PUT_ERR("wrong action\n");
        yac_event_update(event);
      }

      yac_event_delete(event);
    }
  }

  freeCalendar();

  return TEST_EXIT_CODE;
}


