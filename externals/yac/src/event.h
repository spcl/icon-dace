/**
 * @file event.h
 * @brief Structs and interfaces to handle events
 *
 * Event are defined based on the user input to trigger actions
 * for sending and receiving data (put and get action) during the
 * data exchange of (interpolated) fields.
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
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

#ifndef EVENT_H
#define EVENT_H

#include "mtime_datetime.h"
#include "mtime_timedelta.h"
#include "mtime_eventHandling.h"
#include "couple_config.h"

/** \example test_events.c
 * This show how to work with events.
 */

enum yac_action_type {
  NONE               = 0,
  REDUCTION          = 1,
  COUPLING           = 2,
  RESTART            = 3,
  GET_FOR_RESTART    = 4,
  PUT_FOR_RESTART    = 5,
  GET_FOR_CHECKPOINT = 6,
  PUT_FOR_CHECKPOINT = 7,
  OUT_OF_BOUND       = 8,
};

extern int yac_number_of_events;
struct event;

/**
 * Event constructor
 * @return event event
 */
struct event * yac_event_new();

/**
  * Fill in the information for an event
  * @param[in]  event              event
  * @param[in]  model_time_step    model time step as ISO 8061 string
  * @param[in]  coupling_time_step coupling period as ISO 8061 string
  * @param[in]  lag                time lag (in model time steps)
  * @param[in]  time_operation     online processing of coupling fields
  *                                /average/minimum/maximum/accumulate/instant
  * @param[in]  start_date         start date of the job
  * @param[in]  stop_date          end date of the job
  */
void yac_event_add ( struct event * event,
                     char const *model_time_step,
                     char const *coupling_time_step,
                     int lag,
                     enum yac_reduction_type time_operation,
                     const char * start_date,
                     const char * stop_date);

/**
 * Delete an event
 * @param[in] event event
 */
void yac_event_delete ( struct event * event );

/**
 * Check for an event
 * @param[in] event event
 */
enum yac_action_type yac_event_check ( struct event * event );

/**
 * Add delta time to an event
 * @param[in] event event
 */
void yac_event_update ( struct event * event );

/**
 * Get the time lag of an event
 * @param[in] event event
 */
int yac_get_event_lag ( struct event * event );

/**
 * Get the time operation of an event
 * @param[in] event event
 */
int yac_get_event_time_operation ( struct event * event );

/**
 * Get the coupling time step for an event
 * @param[in] event         event
 * @param[in] timedelta_str string memory
 */
char * yac_get_event_coupling_timestep (
  struct event * event, char * timedelta_str );

/**
 * Get the model time step for an event
 * @param[in] event         event
 * @param[in] timedelta_str string memory
 */
char * yac_get_event_model_timestep (
  struct event * event, char * timedelta_str );

/**
 * Get the current datetime for an event
 * @param[in] event         event
 * @param[in] datetime_str string memory
 * @return pointer to datetime_str on success, NULL otherwise
 * @remark the buffer associated to datetime_str has to have at least a
 *         size of MAX_DATETIME_STR_LEN
 */
char * yac_get_event_current_datetime(
  struct event * event, char * datetime_str);

/**
 * Converts time string to ISO 8601 time string
 * @param[in] time      time string
 * @param[in] time_unit unit of time string
 * @return ISO 8601 time string
 */
char const * yac_time_to_ISO(
  char const * time, enum yac_time_unit_type time_unit);

#endif
