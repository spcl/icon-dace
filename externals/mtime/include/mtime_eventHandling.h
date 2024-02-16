// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_eventHandling.h
 *
 * @brief  Definition of the basic event type and its methods. 
 *         Definition of Even-Groups which contains a list of events.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 * 
 * @date March 2013
 */


#ifndef _MTIME_EVENTHANDLING_H
#define _MTIME_EVENTHANDLING_H

#include <stdint.h>
#include <stdbool.h>


///provides a Maximum string length for Group Names.
#define MAX_EVENT_STR_LEN 132 // TODO- On Luis. Need to agree on how a event_string should look.
///provides a Maximum string length for Group Names.
#define MAX_GROUPNAME_STR_LEN 132
///provides a Maximum string length for Event names.
#define MAX_EVENTNAME_STR_LEN 132

struct _datetime;
struct _timedelta;
struct _event;


/**
 * @struct _eventGroup
 *
 * @brief struct _eventGroup defines an Event-group. Each event group has an associated list of events.
 *        Event group is a place holder to 'group' events based on some user defined charteristics.
 *              
 */

struct _eventGroup
{
  int64_t               eventGroupId;           ///< Event Group's ID.
  char*                 eventGroupName;         ///< Event Group's name.

  struct _event*        rootEvent;              ///< Pointer to a List of events in this group.
};

struct _eventGroup*
newEventGroup(const char* _eventGroupName);

void
deallocateEventGroup(struct _eventGroup* eg);

bool
addNewEventToEventGroup(struct _event* e, struct _eventGroup* eg);

bool
removeEventFromEventGroup(char* eventName, struct _eventGroup* eg);

int64_t
getEventGroupId(struct _eventGroup* eg);

char*
getEventGroupName(struct _eventGroup* eg, char* gname);

struct _event*
getEventGroupRootEvent(struct _eventGroup* eg);


/**
 * @struct _event
 *
 * @brief struct _event defines events. Events are set to be triggered at pre-specified intervals. 
 *              
 */
/* Note: all updates here must be mirrored in src/mtime_c_bindings.f90! */
struct _event
{
  struct _event*                nextEventInGroup;            ///< Pointer to the next event in a given Event Group.

  int64_t                       eventId;                     ///< Auto generated. (For future use. Not in use yet).
  char*                         eventName;                   ///< Event's name.

  struct _datetime*             eventsLastEvaluationDateTime;///< Last evaluation datetime to allow for multiple queries. 
  struct _datetime*             eventReferenceDateTime;      ///< Anchor datetime. Can be NULL; Initialized to eventFirstDateTime. .

  struct _datetime*             eventFirstDateTime;          ///< Start datetime of the event. Can be NULL; Initialized to 0-01-01T00:00:00.000.
  struct _datetime*             eventLastDateTime;           ///< Last datetime of the event. Can be NULL; stays a NULL pointer: Logically NULL.

  struct _timedelta*            eventInterval;               ///< TimeDelta between succesive triggers. MUST be specified.
  struct _timedelta*            eventOffset;                 ///< Shift eventReferenceDateTime by 'eventOffset' . Can be NULL; stays a NULL pointer: Logically 0.

  bool                          neverTriggerEvent;           ///< for cases an event is set by an interval of 0
  
  bool                          triggerCurrentEvent;         ///< is this event active?

  bool                          nextEventIsFirst;            ///< Is the next scheduled event the first to be triggered?
  bool                          lastEventWasFinal;           ///< Was the previously triggered event final?
 
  //Properties of the event being triggered.
  /* The following flags make sense IFF the events are CURRENTLY in a triggered state i.e triggerCurrentEvent==true.Else the behavior is undefined.*/
  bool                          eventisFirstInDay;           ///< Is the event being triggered (i.e event when triggerCurrentEvent is true) First-In-Day?
  bool                          eventisFirstInMonth;         ///< Is the event being triggered (i.e event when triggerCurrentEvent is true) First-In-Month?
  bool                          eventisFirstInYear;          ///< Is the event being triggered (i.e event when triggerCurrentEvent is true) First-In-Year?
  bool                          eventisLastInDay;            ///< Is the event being triggered (i.e event when triggerCurrentEvent is true) Last-In-Day?
  bool                          eventisLastInMonth;          ///< Is the event being triggered (i.e event when triggerCurrentEvent is true) Last-In-Month?
  bool                          eventisLastInYear;           ///< Is the event being triggered (i.e event when triggerCurrentEvent is true) Last-In-Year?

  struct _datetime*             triggerNextEventDateTime;       ///< Trigger the event next at this DT.
  struct _datetime*             triggeredPreviousEventDateTime; ///< Last Trigger event happened at this DT. Initialized with "0-01-01T00:00:00.000".
};

struct _event*
newEvent(const char* _eventName, const char* _eventReferenceDate, const char* _eventFirstDate, const char* _eventLastDate, const char* _eventInterval,  const char* _eventOffset);

struct _event*
newEventWithDataType(const char* _eventName, struct _datetime* _eventReferenceDate, struct _datetime* _eventFirstDate, struct _datetime* _eventLastDate, struct _timedelta* _eventInterval, struct _timedelta* _eventOffset);

void
deallocateEvent(struct _event* e);

struct _event*
constructAndCopyEvent(struct _event* e);

bool
isCurrentEventActive(struct _event* e, struct _datetime* current_dt, struct _timedelta* plus_slack, struct _timedelta* minus_slack);

bool
iseventNextInNextDay(struct _event* e);

bool
iseventNextInNextMonth(struct _event* e);

bool
iseventNextInNextYear(struct _event* e);

char*
eventToString(struct _event* e, char* string);

struct _datetime*
getTriggerNextEventAtDateTime(struct _event* e, struct _datetime* dt_current, struct _datetime* dt_return);

struct _datetime*
getTriggeredPreviousEventAtDateTime(struct _event* e, struct _datetime* dt_return);

struct _event*
getNextEventInGroup(struct _event* e);

int64_t
getEventId(struct _event* e);

char*
getEventName(struct _event* e, char* ename);

struct _datetime*
getEventReferenceDateTime(struct _event* e);

struct _datetime*
getEventFirstDateTime(struct _event* e);

struct _datetime*
getEventLastDateTime(struct _event* e);

struct _timedelta*
getEventInterval(struct _event* e);

bool
getNextEventIsFirst(struct _event* e);

bool
getEventisFirstInDay(struct _event* e);

bool
getEventisFirstInMonth(struct _event* e);

bool
getEventisFirstInYear(struct _event* e);

bool
getEventisLastInDay(struct _event* e);

bool
getEventisLastInMonth(struct _event* e);

bool
getEventisLastInYear(struct _event* e);

/**
 * @}
 */

#endif
