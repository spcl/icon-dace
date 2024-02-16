/**
 * @file mtime_eventHandling.c
 *
 * @brief Event-groups which contains a list of events.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "mtime_eventHandling.h"

#include "mtime_datetime.h"
#include "mtime_timedelta.h"
#include "mtime_julianDay.h"
#include "mtime_eventList.h"
#include "mtime_iso8601.h"

/* Local static functions. */
static struct _datetime*
getEventFirstTriggerDateTime(struct _datetime*, struct _timedelta*, struct _timedelta*, struct _datetime*, struct _datetime*);

// The IDs are unique only when they are live. Once a group/event has been deleted, the IDs will be reused.
// Currently, the IDs are generated but not used.
static int64_t EVENTGROUPID = 0;
static int64_t EVENTID = 0;


/**
 * @brief Construct new event-Group using a string.
 *
 * @param  egn
 *         A pointer to char. This string contains the name of event group.
 *
 * @return eg
 *         A pointer to a initialized event-Group.
 *
 */

struct _eventGroup*
newEventGroup(const char* egn)
{
if ((egn != NULL) && (getCalendarType()))
{
  struct _eventGroup* eg = (struct _eventGroup*)calloc(1,sizeof(struct _eventGroup));
  if (eg == NULL)
    return NULL ;

  /* Copy Group name. */
  eg->eventGroupName = (char*)calloc(MAX_GROUPNAME_STR_LEN,sizeof(char));
  if (eg->eventGroupName == NULL )
    {
      free(eg);
      eg = NULL;
      return NULL ;
    }
  strncpy(eg->eventGroupName, egn, MAX_GROUPNAME_STR_LEN-1);
  eg->eventGroupName[MAX_GROUPNAME_STR_LEN - 1] = '\0';

  /* Generate eventGroupId. */
  EVENTGROUPID = EVENTGROUPID + 1;
  eg->eventGroupId = EVENTGROUPID;

  /* Initialize a NULL pointer. Stores a List of events associated with 'this' Event Group. */
  eg->rootEvent = NULL;

  return eg;
}
else
  return NULL;
}


/**
 * @brief Destructor of EventGroup.
 *
 * @param  eg
 *         A pointer to struct _eventGroup. eg is deallocated.
 */

void
deallocateEventGroup(struct _eventGroup* eg)
{
  if (eg != NULL )
    {
      EVENTGROUPID = EVENTGROUPID - 1;

      if(eg->eventGroupName)
        {
          free(eg->eventGroupName);
          eg->eventGroupName = NULL;
        }

      /* Deallocate all events in the list. */
      deallocateEventsInGroup(eg);

      free(eg);
      eg = NULL;
    }
}


/**
 * @brief Add new event to an eventgroup.
 *
 * @param  e
 *         A pointer to struct _event. The event to be added.
 *
 * @param  eg
 *         A pointer to struct _eventGroup. The eventGroup where the event is added.
 *
 * @return bool
 *         true/false indicating success or failure of addition.
 */

bool
addNewEventToEventGroup(struct _event* e, struct _eventGroup* eg)
{
  if ((e != NULL) && (eg != NULL))
    return addNewEventToGroup(e,eg);
  else
    return false;
}


/**
 * @brief Remove event from eventgroup. CRITICAL: Also, deallocate the event.
 *
 * @param  en
 *         A pointer to char. The name of event to be removed.
 *
 * @param  eg
 *         A pointer to struct _eventGroup. The eventGroup to which this event belongs.
 *
 * @return bool
 *         true/false indicating success or failure of removal.
 */

bool
removeEventFromEventGroup(char* en, struct _eventGroup* eg)
{
  if ((en != NULL )&& (eg != NULL) )
    return removeEventWithNameFromGroup(en,eg);
  else
    return false;
}

int64_t
getEventGroupId(struct _eventGroup* eg)
{
  if (eg != NULL)
    return eg->eventGroupId;
  else
    return 0;
}

char*
getEventGroupName(struct _eventGroup* eg, char* gname)
{
  if ((eg != NULL) && (gname != NULL))
    {
      strncpy(gname, eg->eventGroupName, MAX_GROUPNAME_STR_LEN-1);
      gname[MAX_GROUPNAME_STR_LEN-1] = '\0';
      return gname;
    }
  else
    return NULL;
}

struct _event*
getEventGroupRootEvent(struct _eventGroup* eg)
{
  if (eg != NULL)
    return eg->rootEvent;
  else
    return NULL;
}


  /**
   * @brief Construct new event using strings.
   *
   * To initialize an event, the event name and event interval must be specified. The reference date
   * (also known as the anchor-date), start date and event last date are optional. If first date is not
   * specified, first date is initialized as 0-01-01T00:00:00.000. If reference date is not specified, first date
   * is copied to reference date. If event last date is not defined, it is left as NULL and is logically equivalent
   * to Infinity.
   *
   * Notice that events trigger every _eventReferenceDT + Integral multiple of _eventInterval strictly in [_eventFirstDT,_eventLastDT].
   *
   * @param  _en
   *         A pointer to char. This string contains the name of event. String name must be defined.
   *
   * @param  _eventReferenceDT
   *         A pointer to char. This string contains the Anchor date. Events are triggered at Anchor date + (N * eventInterval) where N is in Z.
   *         This can be NULL in which case Anchor Date = First Date.
   *
   * @param  _eventFirstDT
   *         A pointer to char. This string contains the Starting datetime (A) of the DateTime range [A-B] in which trigger is allowed.
   *         Can be NULL.
   *
   * @param  _eventLastDT
   *         A pointer to char. This string contains the Ending datetime (B) of the DateTime range [A-B] in which trigger is allowed.
   *         Can be NULL. If defined, events will not trigger beyond this point. Else, equivalent to infinity.
   *
   * @param  _eventOffset
   *         A pointer to char. Adds a logical shift to _eventReferenceDT and the sifted value is used as the anchor date.
   *         Can be NULL; Logically equivalent to 0 shift.
   *
   * @param  _eventInterval
   *         A pointer to char. This string contains the timestep. This must be defined for every event.
   *
   * @return e
   *         A pointer to an initialized event.
   *
   */

struct _event*
newEvent(const char* _en,
	 const char* _eventReferenceDT,
	 const char* _eventFirstDT,
	 const char* _eventLastDT,
	 const char* _eventInterval,
	 const char* _eventOffset)
{
  if ((_en != NULL) &&  (_eventInterval != NULL) && getCalendarType())
  {
  struct _event* e = (struct _event*)calloc(1,sizeof(struct _event));
  if (e == NULL )
    return NULL;

  /* Initialize a null pointer. Connectes to the 'next' event in the list of Events in an Event Group.  */
  e->nextEventInGroup = NULL;

  /* Copy event name. */
  e->eventName = (char*)calloc(MAX_EVENTNAME_STR_LEN,sizeof(char));
  if (e->eventName == NULL )
    goto cleanup_and_exit;
  strncpy(e->eventName, _en, MAX_EVENTNAME_STR_LEN-1);
  e->eventName[MAX_EVENTNAME_STR_LEN - 1] = '\0';

  /* Generate eventId. */
  EVENTID = EVENTID + 1;
  e->eventId = EVENTID;

  /* Last evaluation date */
  e->eventsLastEvaluationDateTime = NULL;

  /* Start Date. */
  if (_eventFirstDT)
    e->eventFirstDateTime = newDateTime(_eventFirstDT);
  else
    e->eventFirstDateTime = newDateTime(initDummyDTString); // Default start date = "0-01-01T00:00:00.000".

  if (e->eventFirstDateTime == NULL )
    goto cleanup_and_exit;


  /* Anchor date. */
  if (_eventReferenceDT)
      e->eventReferenceDateTime = newDateTime(_eventReferenceDT);
  else
      e->eventReferenceDateTime = constructAndCopyDateTime(e->eventFirstDateTime); // If Anchor date is not defined, Anchor Date = First Date.

   if (e->eventReferenceDateTime == NULL )
     goto cleanup_and_exit;


  /* Last date. */
  if (_eventLastDT)
    {
      e->eventLastDateTime = newDateTime(_eventLastDT);
      if (e->eventLastDateTime == NULL )
        goto cleanup_and_exit;
    }
  else
    {
      e->eventLastDateTime = NULL; // Logically equivalent to Inf.
    }

  /* _eventInterval must be defined. */
  e->eventInterval = newTimeDelta(_eventInterval);
  if (e->eventInterval == NULL)
    goto cleanup_and_exit;

  if (e->eventInterval->year == 0 && e->eventInterval->month == 0 && e->eventInterval->day == 0
      && e->eventInterval->hour == 0 && e->eventInterval->minute == 0 && e->eventInterval->second == 0 && e->eventInterval->ms == 0)
    {
      e->neverTriggerEvent   = true;

      e->triggerCurrentEvent = false;

      e->nextEventIsFirst    = false;
      e->lastEventWasFinal   = false;

      e->eventisFirstInDay   = false;
      e->eventisFirstInMonth = false;
      e->eventisFirstInYear  = false;
      e->eventisLastInDay    = false;
      e->eventisLastInMonth  = false;
      e->eventisLastInYear   = false;

      /* Removed to prevent segfaults ... deallocateTimeDelta(e->eventInterval); */

      e->eventOffset = NULL;
      e->triggeredPreviousEventDateTime = NULL;
      e->triggerNextEventDateTime = NULL;

      return e;
    }

  if (_eventOffset)
    {
      e->eventOffset = newTimeDelta(_eventOffset);
      if (e->eventOffset == NULL)
	goto cleanup_and_exit;
    }
  else
    e->eventOffset = NULL; //Logically equivalent to 0.

  /* Initialize with 'some' (arbitrary) value. */
  e->triggeredPreviousEventDateTime = newDateTime(initDummyDTString);
  if (e->triggeredPreviousEventDateTime == NULL)
    goto cleanup_and_exit;

  e->triggerNextEventDateTime = newDateTime(initDummyDTString);
  /* triggerNextEventDateTime holds the DateTime when an event should be triggered subjectively (i.e assuming real time starts at -Inf and moves fwd.)
     At init (right now), triggerNextEventDateTime stores the DateTime of First-Ever trigger.
     Events trigger at e->eventReferenceDateTime (+ e-> eventOffset ) + N * e->eventInterval where N is a positive, negative integer.
     Hence, first trigger happens at the nearest possible trigger >= e->eventFirstDateTime. */
  if (!getEventFirstTriggerDateTime(e->eventFirstDateTime, e->eventInterval, e->eventOffset, e->eventReferenceDateTime, e->triggerNextEventDateTime))
    goto cleanup_and_exit;

  /* Init the Flags. */
  e->neverTriggerEvent          = false;

  e->triggerCurrentEvent 	= false;

  e->nextEventIsFirst 		= true;
  e->lastEventWasFinal          = false;

  e->eventisFirstInDay 		= false;
  e->eventisFirstInMonth 	= false;
  e->eventisFirstInYear 	= false;
  e->eventisLastInDay 		= false;
  e->eventisLastInMonth  	= false;
  e->eventisLastInYear 		= false;

  return e;

  cleanup_and_exit:
		if (e)
		  {
		    if (e->eventName)
		      {
			free(e->eventName);
			e->eventName = NULL;
		      }
		    deallocateDateTime(e->eventReferenceDateTime);
		    e->eventReferenceDateTime = NULL;
		    deallocateDateTime(e->eventFirstDateTime);
		    e->eventFirstDateTime = NULL;
		    deallocateDateTime(e->eventLastDateTime);
		    e->eventLastDateTime = NULL;
		    deallocateTimeDelta(e->eventInterval);
		    e->eventInterval = NULL;
                    deallocateTimeDelta(e->eventOffset);
                    e->eventOffset = NULL;
		    deallocateDateTime(e->triggeredPreviousEventDateTime);
		    e->triggeredPreviousEventDateTime = NULL;
                    deallocateDateTime(e->triggerNextEventDateTime);
                    e->triggerNextEventDateTime = NULL;
		    free(e);
		    e = NULL;
		  }
		return NULL ;
  }
  else
    return NULL;
}


  /**
   * @brief Construct new event using data-types.
   *
   * To initialize an event, the event name and event interval must be specified. The reference date
   * (also known as the anchor-date), start date and event last date are optional. If first date is not
   * specified, first date is initialized as 0-01-01T00:00:00.000. If reference date is not specified, first date
   * is copied to reference date. If event last date is not defined, it is left as NULL and is logically equivalent
   * to Infinity.
   *
   * Notice that events trigger every _eventReferenceDT + Integral multiple of _eventInterval strictly in [_eventFirstDT,_eventLastDT].
   *
   * @param  _en
   *         A pointer to char. This string contains the name of event.
   * @param  _eventReferenceDT
   *         A pointer to struct _datetime. This pointer contains the First Reference date also called anchor date.
   *         This can be NULL, but if defined acts as the true starting datetime overriding e->eventFirstDateTime.
   * @param  _eventFirstDT
   *         A pointer to struct _datetime. This pointer contains the Starting datetime. This must be defined for every Event.
   * @param  _eventLastDT
   *         A pointer to struct _datetime. This pointer contains the Ending datetime.
   *         This can be NULL. If defined, events will not trigger beyond this point.
   * @param  _eventInterval
   *         A pointer to struct _timedelta. This pointer contains the timestep. This must be defined for every event.
   * @return e
   *         A pointer to an initialized event.
   *
   */

struct _event*
newEventWithDataType(const char* _en,
			struct _datetime*  _eventReferenceDT,
			struct _datetime*  _eventFirstDT,
			struct _datetime*  _eventLastDT,
			struct _timedelta* _eventInterval,
			struct _timedelta* _eventOffset)
{
  if ((_en != NULL) &&  (_eventInterval != NULL) && getCalendarType())
  {
  struct _event* e = (struct _event*)calloc(1,sizeof(struct _event));
  if (e == NULL )
    return NULL;

  /* Initialize a null pointer. Connectes to the 'next' event in the list of Events in an Event Group.  */
  e->nextEventInGroup = NULL;

  /* Copy event name. */
  e->eventName = (char*)calloc(MAX_EVENTNAME_STR_LEN,sizeof(char));
  if (e->eventName == NULL )
    goto cleanup_and_exit;
  strncpy(e->eventName, _en, MAX_EVENTNAME_STR_LEN-1);
  e->eventName[MAX_EVENTNAME_STR_LEN - 1] = '\0';

  /* Generate eventId. */
  EVENTID = EVENTID + 1;
  e->eventId = EVENTID;

  /* Start Date. */
  if (_eventFirstDT)
    e->eventFirstDateTime = constructAndCopyDateTime(_eventFirstDT);
  else
    e->eventFirstDateTime = newDateTime(initDummyDTString); // Default start date = "0-01-01T00:00:00.000".

  if (e->eventFirstDateTime == NULL )
    goto cleanup_and_exit;


  /* Anchor date. */
  if (_eventReferenceDT)
      e->eventReferenceDateTime = constructAndCopyDateTime(_eventReferenceDT);
  else
      e->eventReferenceDateTime = constructAndCopyDateTime(e->eventFirstDateTime); // If Anchor date is not defined, Anchor Date = First Date.

   if (e->eventReferenceDateTime == NULL )
     goto cleanup_and_exit;


  /* Last date. */
  if (_eventLastDT)
      e->eventLastDateTime = constructAndCopyDateTime(_eventLastDT);
  else
      e->eventLastDateTime = NULL; // Logically equivalent to Inf.


  if (_eventInterval->year == 0 && _eventInterval->month == 0 && _eventInterval->day == 0
      && _eventInterval->hour == 0 && _eventInterval->minute == 0 && _eventInterval->second == 0 && _eventInterval->ms == 0)
    {
      e->neverTriggerEvent   = true;

      e->triggerCurrentEvent = false;

      e->nextEventIsFirst    = false;
      e->lastEventWasFinal   = false;

      e->eventisFirstInDay   = false;
      e->eventisFirstInMonth = false;
      e->eventisFirstInYear  = false;
      e->eventisLastInDay    = false;
      e->eventisLastInMonth  = false;
      e->eventisLastInYear   = false;

      /* Added to prevent segfaults ... */
      e->eventInterval = constructAndCopyTimeDelta(_eventInterval);
      e->eventOffset = NULL;
      e->triggeredPreviousEventDateTime = NULL;
      e->triggerNextEventDateTime = NULL;

      return e;
    }

  /* _eventInterval must be defined. */
  e->eventInterval = constructAndCopyTimeDelta(_eventInterval);
  if (e->eventInterval == NULL)
    goto cleanup_and_exit;


  /* Event trigger offset. */
  if (_eventOffset)
    e->eventOffset = constructAndCopyTimeDelta(_eventOffset);
  else
    e->eventOffset = NULL;

  /* Initialize with 'some'(arbitrary) value. */
  e->triggeredPreviousEventDateTime = newDateTime(initDummyDTString);
  if (e->triggeredPreviousEventDateTime == NULL)
    goto cleanup_and_exit;

  e->triggerNextEventDateTime = newDateTime(initDummyDTString);
  /* triggerNextEventDateTime holds the DateTime when an event should be triggered subjectively (i.e assuming real time starts at -Inf and moves fwd.)
     At init (right now), triggerNextEventDateTime stores the DateTime of First-Ever trigger.
     Events trigger at e->eventReferenceDateTime (+ e->eventOffset) + N*e->eventInterval where N is a positive, negative or zero integer.
     Hence, first trigger happens at the nearest possible trigger >= e->eventFirstDateTime. */
  if (!getEventFirstTriggerDateTime(e->eventFirstDateTime, e->eventInterval, e->eventOffset, e->eventReferenceDateTime, e->triggerNextEventDateTime))
    goto cleanup_and_exit;

  /* Init the Flags. */
  e->neverTriggerEvent          = false;

  e->triggerCurrentEvent 	= false;

  e->nextEventIsFirst 		= true;
  e->lastEventWasFinal          = false;

  e->eventisFirstInDay 		= false;
  e->eventisFirstInMonth 	= false;
  e->eventisFirstInYear 	= false;
  e->eventisLastInDay 		= false;
  e->eventisLastInMonth  	= false;
  e->eventisLastInYear 		= false;

  return e;

  cleanup_and_exit:
		if (e)
		  {
		    if (e->eventName)
		      {
			free(e->eventName);
			e->eventName = NULL;
		      }
		    deallocateDateTime(e->eventReferenceDateTime);
		    e->eventReferenceDateTime = NULL;
		    deallocateDateTime(e->eventFirstDateTime);
		    e->eventFirstDateTime = NULL;
		    deallocateDateTime(e->eventLastDateTime);
		    e->eventLastDateTime = NULL;
		    deallocateTimeDelta(e->eventInterval);
		    e->eventInterval = NULL;
                    deallocateTimeDelta(e->eventOffset);
                    e->eventOffset = NULL;
		    deallocateDateTime(e->triggeredPreviousEventDateTime);
		    e->triggeredPreviousEventDateTime = NULL;
                    deallocateDateTime(e->triggerNextEventDateTime);
                    e->triggerNextEventDateTime = NULL;
		    free(e);
		    e = NULL;
		  }
		return NULL ;
  }
  else
    return NULL;
}


/**
 * @brief Destructor of Event.
 *
 * @param  e
 *         A pointer to struct _event. e is deallocated.
 */

void
deallocateEvent(struct _event* e)
{
  if (e != NULL )
    {
      EVENTID = EVENTID - 1;

      if(e->eventName)
        {
          free(e->eventName);
          e->eventName = NULL;
       }

      if (e->eventsLastEvaluationDateTime != NULL)
	deallocateDateTime(e->eventsLastEvaluationDateTime);
      deallocateDateTime(e->eventReferenceDateTime);
      deallocateDateTime(e->eventFirstDateTime);
      deallocateDateTime(e->eventLastDateTime);
      deallocateDateTime(e->triggerNextEventDateTime);
      deallocateDateTime(e->triggeredPreviousEventDateTime);

      deallocateTimeDelta(e->eventInterval);
      deallocateTimeDelta(e->eventOffset);

      /* Warning: Do not free the NextEventInGroup. That's not your job. */
      e->nextEventInGroup = NULL;

      free(e);
      e = NULL;
    }
}

struct _event*
constructAndCopyEvent(struct _event* ev)
{
  struct _event* ev_copy;

  if ( ev != NULL )
    {
      ev_copy = newEventWithDataType(ev->eventName,
                                     ev->eventReferenceDateTime,
                                     ev->eventFirstDateTime,
                                     ev->eventLastDateTime,
                                     ev->eventInterval,
                                     ev->eventOffset);

      ev_copy->eventId = ev->eventId;
      ev_copy->nextEventInGroup = ev->nextEventInGroup;

      ev_copy->neverTriggerEvent =   ev->neverTriggerEvent;
      ev_copy->triggerCurrentEvent = ev->triggerCurrentEvent;

      ev_copy->nextEventIsFirst =    ev->nextEventIsFirst;
      ev_copy->lastEventWasFinal =   ev->lastEventWasFinal;


      ev_copy->eventisFirstInDay =   ev->eventisFirstInDay;
      ev_copy->eventisFirstInMonth = ev->eventisFirstInMonth;
      ev_copy->eventisFirstInYear =  ev->eventisFirstInYear;
      ev_copy->eventisLastInDay =    ev->eventisLastInDay;
      ev_copy->eventisLastInMonth =  ev->eventisLastInMonth;
      ev_copy->eventisLastInYear =   ev->eventisLastInYear;

      ev_copy->triggeredPreviousEventDateTime = constructAndCopyDateTime(ev->triggeredPreviousEventDateTime);
      ev_copy->triggerNextEventDateTime = constructAndCopyDateTime(ev->triggerNextEventDateTime);

      return ev_copy;
    }
    else
      return NULL;
}




/* INTERNAL FUNCTION. */
//Get if trigger is true. Trigger true in [T-minus_slack,T+plus_slack].
static
compare_return_val
isTriggerTimeInRange(struct _datetime* current_dt, struct _datetime* triggerNextEventDateTime, struct _timedelta* plus_slack, struct _timedelta* minus_slack)
{
  int cmp_val_flag = -128;

  /* Make a local copy of slack to avoid updating the user supplied timedeltas. */
  struct _timedelta* plus_slack_local  = NULL;
  struct _timedelta* minus_slack_local = NULL;

  if (plus_slack)
          plus_slack_local  = constructAndCopyTimeDelta(plus_slack);
  else
          plus_slack_local  = newTimeDelta(initDummyTDString); //Translates to 0.

  if (minus_slack)
          minus_slack_local = constructAndCopyTimeDelta(minus_slack);
  else
          minus_slack_local = newTimeDelta(initDummyTDString); //Translates to 0.



  /* If plus_slack_local is defined, return the status of current_dt vis-a-vis trigger time +/- allowed delta.  */
  if((plus_slack_local->sign == '+') && ((minus_slack_local->sign == '+')))
    {
      int upper_val_flag = -128;
      int lower_val_flag = -128;

      struct _datetime *dt_upperbound              = newDateTime(initDummyDTString);
      struct _datetime *dt_lowerbound              = newDateTime(initDummyDTString);

      minus_slack_local->sign = '-'; /* Change sign to obtain substraction. */


      upper_val_flag = compareDatetime(current_dt,addTimeDeltaToDateTime(triggerNextEventDateTime,plus_slack_local,dt_upperbound));
      lower_val_flag = compareDatetime(current_dt,addTimeDeltaToDateTime(triggerNextEventDateTime,minus_slack_local,dt_lowerbound));

      if (
           (upper_val_flag == less_than || upper_val_flag == equal_to)
           &&
           (lower_val_flag == greater_than || lower_val_flag == equal_to)
         )
        {
          cmp_val_flag = equal_to;
        }
      else if (upper_val_flag == greater_than)
        {
          cmp_val_flag = greater_than;
        }
      else if (lower_val_flag == less_than)
        {
          cmp_val_flag = less_than;
        }

      //Cleaup.
      deallocateDateTime(dt_upperbound);
      deallocateDateTime(dt_lowerbound);
    }
  else /* If slack is malformed (negative sign is not permitted), return as normal (follow exact match for equal).  */
    {
      cmp_val_flag = compareDatetime(current_dt,triggerNextEventDateTime);
    }

  deallocateTimeDelta(plus_slack_local);
  deallocateTimeDelta(minus_slack_local);

  return cmp_val_flag;
}


/**
 * @brief Check if this event is active by comparing event's trigger time with current_dt.
 *
 *        The current_dt must lie in event's trigger time
 *        (subject to optional specified slack: [Trigger_time - minus_slack, Trigger_time + plus_slack]. Slacks can be NULL. Always inclusive.)
 *        The lib has no built-in clock but relies on isCurrentEventActive(.) being called (polled) from the application at fixed intervals.
 *
 * @param  event
 *         A pointer to struct _event. This is the event being tested.
 *
 * @param  current_dt
 *         A pointer to struct _datetime. This is the 'current' datetime of the system.
 *
 * @param  plus_slack
 *         A pointer to struct _timedelta. Events are triggered between [actual_trigger_time, actual_trigger_time + plus_slack].
 * 	   Can be NULL; logically 0. Sign MUST be '+'
 *
 * @param  minus_slack
 *         A pointer to struct _timedelta. Events are triggered between [actual_trigger_time - minus_slack, actual_trigger_time].
 *         Can be NULL; logically 0. Sign MUST be '+'
 *
 * @return bool
 *         true/false indicating if the event is active.
 */

bool
isCurrentEventActive(struct _event* event, struct _datetime* current_dt, struct _timedelta* plus_slack, struct _timedelta* minus_slack)
{

  if ((event != NULL ) && (current_dt != NULL)) { // allowed slacks can be NULL.

    if (event->neverTriggerEvent)
      {
        return false;
      }

    /* If Event has been triggered for a datetime and checked once again return that state */

    if (event->eventsLastEvaluationDateTime != NULL)
      {
	if (compareDatetime(current_dt, event->eventsLastEvaluationDateTime) == equal_to)
	  {
	    return true;
	  }
      }

    /* Reset flags. */
    event->triggerCurrentEvent        = false;

    event->eventisFirstInDay          = false;
    event->eventisFirstInMonth        = false;
    event->eventisFirstInYear         = false;
    event->eventisLastInDay           = false;
    event->eventisLastInMonth         = false;
    event->eventisLastInYear          = false;

    //++++++++ Events trigger [start-end]. Slack is ignored outside this range. ++++++++++//
    /* If current_dt is yet not in trigger range, Return on false. */
    if (compareDatetime(current_dt,event->eventFirstDateTime) == less_than)
      return false;
    /* If last date is defined, check if current_dt is ahead of event->eventLastDateTime. Return on false. */
    if (event->eventLastDateTime && (compareDatetime(current_dt,event->eventLastDateTime) == greater_than))
      {
	/* If current_dt > event->eventLastDateTime and the event has never triggered, event->lastEventWasFinal should stay false. */
	if(!event->nextEventIsFirst)
	  event->lastEventWasFinal = true;  //No further trigger possible.

	return false;
      }


    /* In case the current_dt is ahead of event->triggerNextEventDateTime, we need to update the event->triggerNextEventDateTime
       to the next possible trigger time or else the events will never trigger. */
    if (isTriggerTimeInRange(current_dt, event->triggerNextEventDateTime, plus_slack, NULL) == greater_than)
      {
	/* Get the first trigger time and update. */
	struct _timedelta* modulo_td = newTimeDelta(initDummyTDString);
	moduloTimeDeltaFromDateTime(event->triggerNextEventDateTime, event->eventInterval, current_dt, modulo_td);
	addTimeDeltaToDateTime(current_dt,modulo_td,event->triggerNextEventDateTime);
	/* Cleanup. */
	deallocateTimeDelta(modulo_td);
      }

    /* Check if trigger time is now. Trigger allowed with a slack provided [start-end] condition is met. */
    if( (!event->lastEventWasFinal) && (isTriggerTimeInRange(current_dt, event->triggerNextEventDateTime, plus_slack, minus_slack) == equal_to) )
      {
	/* If current Datetime is equal to next trigger datetime, Event is active. */

	/* If the event being triggred is the first event to be triggered, it is all of FirstIn*  */
	/* If previous triggered event had a different day/month/year, this event must be FirstIn* */
	if( (event->nextEventIsFirst == true) || (iseventNextInNextDay(event) == true))
	  {
	    event->eventisFirstInDay = true;
	  }
	if( (event->nextEventIsFirst == true)  || (iseventNextInNextMonth(event) == true))
	  {
	    event->eventisFirstInMonth = true;
	  }
	if( (event->nextEventIsFirst == true)  || (iseventNextInNextYear(event) == true))
	  {
	    event->eventisFirstInYear = true;
	  }

	/* Set previous-triggered-datetime to now. */
	replaceDatetime(current_dt, event->triggeredPreviousEventDateTime);

	int cmp_val_flag = -128;
	struct _datetime* tmp_dt = newDateTime(initDummyDTString);

	/* Set the new next-trigger datetime. If eventLastDateTime is defined, then update should not happen
	   if current_dt + eventInterval > eventLastDateTime */
	if (
            !( event->eventLastDateTime
               &&
               (
		( cmp_val_flag = isTriggerTimeInRange ( addTimeDeltaToDateTime(current_dt, event->eventInterval,tmp_dt),
							event->eventLastDateTime,
							NULL,
							minus_slack
							)
		  ) == greater_than
		 )
	       )
	     )
	  {
	    addTimeDeltaToDateTime(event->triggerNextEventDateTime,event->eventInterval,event->triggerNextEventDateTime);
	  }
	else
	  event->lastEventWasFinal = true;

	deallocateDateTime(tmp_dt); // Cleanup.

	/* Set event.*/
	event->triggerCurrentEvent = true;


	/* If the future event (not the current event being triggered) has a different day/month/year, current event must be LastIn*.
	   Notice that triggerNextEventDateTime is now the FUTURE event after the update above. */
	if (
	    (	event->eventLastDateTime
		&&
		(cmp_val_flag == greater_than)
		)
	    ||
	    (iseventNextInNextDay(event) == true)
	    )
	  {
	    event->eventisLastInDay = true;
	  }

	if(
	   ( event->eventLastDateTime
	     &&
	     (cmp_val_flag == greater_than)
	     )
	   ||
	   (iseventNextInNextMonth(event) == true))
	  {
	    event->eventisLastInMonth = true;
	  }

	if (
	    ( event->eventLastDateTime
	      &&
	      (cmp_val_flag == greater_than)
              )
	    ||
	    (iseventNextInNextYear(event) == true))
	  {
	    event->eventisLastInYear = true;
	  }

	/* Reset indicating this is no longer true. */
	event->nextEventIsFirst = false;

	if (event->eventsLastEvaluationDateTime == NULL)
	  {
	    event->eventsLastEvaluationDateTime = newDateTime(initDummyDTString);
	  }

	event->eventsLastEvaluationDateTime = replaceDatetime(current_dt, event->eventsLastEvaluationDateTime);

	return true;
      }

    return false;
  }
  else

    return false;
}


/* Is next event in next day(s) of previous event. */
bool
iseventNextInNextDay(struct _event* e)
{
  if( e != NULL )
    {
      bool ret = false;

      struct _date* d_next = newDate(initDummyDString);
      struct _date* d_prev = newDate(initDummyDString);
      d_next = convertDateTimeToDate(e->triggerNextEventDateTime,d_next);
      d_prev = convertDateTimeToDate(e->triggeredPreviousEventDateTime,d_prev);

      if(compareDate(d_next,d_prev) == greater_than)
	{
	  ret = true;
	}
      else
	{
	  ret = false;
	}

      deallocateDate(d_next);
      deallocateDate(d_prev);

      return ret;
    }
  else
    return false;
}


/* Is next event in next month(s) of previous event. */
bool
iseventNextInNextMonth(struct _event* e)
{
  if( e != NULL )
    {
      bool ret = false;

      struct _date* d_next = newDate(initDummyDString);
      struct _date* d_prev = newDate(initDummyDString);
      d_next = convertDateTimeToDate(e->triggerNextEventDateTime,d_next);
      d_prev = convertDateTimeToDate(e->triggeredPreviousEventDateTime,d_prev);

      d_next->day = 1;
      d_prev->day = 1;

      if(compareDate(d_next,d_prev) == greater_than)
	{
	  ret = true;
	}
      else
	{
	  ret = false;
	}

      deallocateDate(d_next);
      deallocateDate(d_prev);

      return ret;
    }
  else
    return false;
}

/* Is next event in next year(s) of previous event. */
bool
iseventNextInNextYear(struct _event* e)
{
  if( e != NULL )
    {
      bool ret = false;

      struct _date* d_next = newDate(initDummyDString);
      struct _date* d_prev = newDate(initDummyDString);
      d_next = convertDateTimeToDate(e->triggerNextEventDateTime,d_next);
      d_prev = convertDateTimeToDate(e->triggeredPreviousEventDateTime,d_prev);

      d_next->day = 1;
      d_next->month = 1;

      d_prev->day = 1;
      d_prev->month = 1;


      if(compareDate(d_next,d_prev) == greater_than)
	{
	  ret = true;
	}
      else
      {
	ret = false;
      }

      deallocateDate(d_next);
      deallocateDate(d_prev);

      return ret;
    }
  else
    return false;
}


/**
 * @brief Get Event as a string.
 *
 * @param  e
 *         A pointer to struct _event. The event to be converted to string.
 *
 * @param  string
 *         A pointer to char. String where event is to be written.
 *
 * @return string
 *         A pointer to the string containing event description.
 */
//TODO on Luis. Exact format.
char*
eventToString(struct _event* e, char* string)
{
  if ((e != NULL) && (string != NULL))
    {
      memset(string,'\0',MAX_EVENT_STR_LEN);
      if (e->eventName != NULL)
	{
      	  sprintf(string, "%s", e->eventName);
          return string;
	}
      else
	{
	   //TODO.
	   ;
	}
      return NULL;
    }
  else
    return NULL;
}


/* Calculate the first trigger time. First trigger is the nearest allowed trigger time >= start_dt.
   Triggers are allowd only ref_dt + N * timestep where N can be positive, negative or 0. */
static
struct _datetime*
getEventFirstTriggerDateTime(struct _datetime* start_dt, struct _timedelta* timestep, struct _timedelta* offset, struct _datetime* ref_dt, struct _datetime* first_trigger_dt)
{
  if ((start_dt != NULL) && (timestep != NULL) && (ref_dt != NULL) && (first_trigger_dt != NULL) ) //offset can be NULL.
    {

      if (timestep->year != 0 && timestep->month == 0 && timestep->day == 0
          && timestep->hour == 0 && timestep->minute == 0 && timestep->second == 0 && timestep->ms == 0)
        {
          /* NOTE: years only is a trivial case. */

          /* Get anchor. ref + offset is the real anchor. */
          struct _datetime *anchor = constructAndCopyDateTime(ref_dt);
          addTimeDeltaToDateTime(ref_dt, offset, anchor); // Note: If offset is null, no addition takes place.

	  if (compareDatetime(anchor, start_dt) == greater_than)
	    {
	      replaceDatetime(anchor, first_trigger_dt);
	    }
	  else
	    {
	      /* Determine difference between anchor and start year */
	      int64_t differenceOfYears = start_dt->date.year - anchor->date.year; // result always >= 0.
	      int64_t yearsToAdd = differenceOfYears % timestep->year;

	      /* We only need to update the year */
	      replaceDatetime(anchor, first_trigger_dt);
	      first_trigger_dt->date.year += yearsToAdd;
	    }

          /* Cleanup. */
          deallocateDateTime(anchor);

          /* And return. */
          return first_trigger_dt;
        }

      else if (timestep->year == 0 && timestep->month != 0 && timestep->day == 0
               && timestep->hour == 0 && timestep->minute == 0 && timestep->second == 0 && timestep->ms == 0)
        {
          /* NOTE: months only is again a trivial case. */

          /* Get anchor. ref + offset is the real anchor. */
          struct _datetime *anchor = constructAndCopyDateTime(ref_dt);
          addTimeDeltaToDateTime(ref_dt, offset, anchor); // Note: If offset is null, no addition takes place.

	  if (compareDatetime(anchor, start_dt) == greater_than)
	    {
	      replaceDatetime(anchor, first_trigger_dt);
	    }
	  else
	    {
	      /* Determine difference between anchor and start in months */
	      int differenceOfMonths = 12*(start_dt->date.year - anchor->date.year) + start_dt->date.month - anchor->date.month;

	      int yearsToAdd = differenceOfMonths / 12;
	      int monthsToAdd = differenceOfMonths % timestep->month;

	      /* We only need to update the year and month */
	      replaceDatetime(anchor, first_trigger_dt);
	      first_trigger_dt->date.year += yearsToAdd;
	      first_trigger_dt->date.month += monthsToAdd;
	    }

          /* Cleanup. */
          deallocateDateTime(anchor);

          /* And return. */
          return first_trigger_dt;
        }

      else if (timestep->year == 0 && timestep->month == 0)

        {
          /* NOTE: This code is higly 'rigged', to a point where it might feel micromanaged.
             This is to speed-up the event-init process. Without the hack, the code
             will be too slow and hence the ends justify the means.
          */

          struct _timedelta* timestep_bkw;

          /* Get start date in julian. */
          struct _julianday* start_jd = newJulianDay(0, 0);
          start_jd = date2julian(start_dt, start_jd);


          /* Get timedelta in juliandelta. */
          struct _juliandelta* timestep_jd = newJulianDelta('+', 0, 0);
          timestep_jd = timeDeltaToJulianDelta(timestep, ref_dt, timestep_jd);


          /* For speed-up */
          struct _juliandelta* timestep_fastfwd_jd = newJulianDelta('+', 0, 0);
          struct _juliandelta* timestep_fastbkw_jd = newJulianDelta('+', 0, 0);


          /* Get anchor. ref + offset is the real anchor. */
          struct _datetime* anchor = constructAndCopyDateTime(ref_dt);
          addTimeDeltaToDateTime(ref_dt, offset, anchor); // Note: If offset is null, no addition takes place.

          /* Get anchor in julian. */
          struct _julianday* anchor_jd = newJulianDay(0, 0);
          anchor_jd = date2julian(anchor, anchor_jd);

          /* Optimization hack: Calculate an approx metric are_dates_too_far and speed up the jumps.  */
          int64_t are_dates_too_far = 0;
          if (timestep_jd->day)
            are_dates_too_far = (start_jd->day - anchor_jd->day)/(timestep_jd->day);
          else if (timestep_jd->ms)
            are_dates_too_far = (start_jd->day - anchor_jd->day)/((float)timestep_jd->ms/NO_OF_MS_IN_A_DAY);
          //else ... well, should never happen. If it does, the initialized value of zero persists.


          int lambda = 1; //speed-up const. Default is 1 or no speedup.
          if (are_dates_too_far>10 || are_dates_too_far<-10)
            lambda = 100000; //speed up if start-date and anchor are 'too-far' away.


          /* Fast-Fwd */
          timestep_fastfwd_jd->day = lambda*timestep_jd->day;
          timestep_fastfwd_jd->ms  = lambda*timestep_jd->ms;

          while( timestep_fastfwd_jd->ms >= NO_OF_MS_IN_A_DAY )
            {
              timestep_fastfwd_jd->day++;
              timestep_fastfwd_jd->ms = timestep_fastfwd_jd->ms - NO_OF_MS_IN_A_DAY;
            }


          /* We need to Loop backwards: Create a timestep replica and change the sign to negative to travel-back-in-time. */
          timestep_bkw = constructAndCopyTimeDelta(timestep);
          timestep_bkw->sign = '-';

          struct _juliandelta* timestep_bkw_jd = newJulianDelta('+', 0, 0);
          timeDeltaToJulianDelta(timestep_bkw, ref_dt, timestep_bkw_jd);

          /* Fast-Bkwd */
          timestep_fastbkw_jd->day = lambda*timestep_bkw_jd->day;
          timestep_fastbkw_jd->ms  = lambda*timestep_bkw_jd->ms;

          while( timestep_fastbkw_jd->ms <=  (-1)*(NO_OF_MS_IN_A_DAY) )
            {
              timestep_fastbkw_jd->day--;
              timestep_fastbkw_jd->ms = timestep_fastbkw_jd->ms + NO_OF_MS_IN_A_DAY;
            }


          switch (compareDatetime(start_dt,ref_dt))
            {
            case greater_than: /* start_dt > ref_dt */

              if (lambda != 1)
                {

                  /* Jump very fast and reach the start date quickly. */
                  do
                    {
                      anchor_jd->day = anchor_jd->day + timestep_fastfwd_jd->day;
                      anchor_jd->ms  = anchor_jd->ms  + timestep_fastfwd_jd->ms;

                      if (anchor_jd->ms >= NO_OF_MS_IN_A_DAY)
                        {
                          anchor_jd->day = anchor_jd->day + 1;
                          anchor_jd->ms = anchor_jd->ms - NO_OF_MS_IN_A_DAY;
                        }

                    } while ( !((anchor_jd->day > start_jd->day)  ||  (anchor_jd->day == start_jd->day && anchor_jd->ms > start_jd->ms)) );


                  /* But I jumped one time too much. Move back */
                  anchor_jd->day = anchor_jd->day - timestep_fastfwd_jd->day;
                  anchor_jd->ms  = anchor_jd->ms  - timestep_fastfwd_jd->ms;

                  if (anchor_jd->ms < 0)
                    {
                      anchor_jd->day = anchor_jd->day - 1;
                      anchor_jd->ms = anchor_jd->ms + NO_OF_MS_IN_A_DAY;
                    }
                }

              /* I am close. Now determine the actual time. */
              do
                {
                  anchor_jd->day = anchor_jd->day + timestep_jd->day;
                  anchor_jd->ms = anchor_jd->ms + timestep_jd->ms;

                  if (anchor_jd->ms >= NO_OF_MS_IN_A_DAY)
                    {
                      anchor_jd->day = anchor_jd->day + 1;
                      anchor_jd->ms = anchor_jd->ms - NO_OF_MS_IN_A_DAY;
                    }

                } while ( !((anchor_jd->day > start_jd->day)  ||  (anchor_jd->day == start_jd->day && anchor_jd->ms >= start_jd->ms)) );

              /* anchor_jd is now the true event-trigger-time. Set it to anchor. */
              julian2date(anchor_jd,anchor);

              break;

            case equal_to: /* start_dt == ref_dt */

              break;

            case less_than: /* start_dt < ref_dt */


              /* Jump very fast bkwd and reach the start date quickly. */
              do
                {
                  anchor_jd->day = anchor_jd->day + timestep_fastbkw_jd->day;
                  anchor_jd->ms = anchor_jd->ms + timestep_fastbkw_jd->ms;

                  if ( anchor_jd->ms < 0 )
                    {
                      anchor_jd->day = anchor_jd->day - 1;
                      anchor_jd->ms = anchor_jd->ms + NO_OF_MS_IN_A_DAY;
                    }

                } while(! (anchor_jd->day < start_jd->day  ||  (anchor_jd->day == start_jd->day && anchor_jd->ms < start_jd->ms)) );


              /* I jumped one time too much. Move forward. */
              anchor_jd->day = anchor_jd->day + timestep_fastfwd_jd->day;
              anchor_jd->ms = anchor_jd->ms + timestep_fastfwd_jd->ms;

              if (anchor_jd->ms >= NO_OF_MS_IN_A_DAY)
                {
                  anchor_jd->day = anchor_jd->day + 1;
                  anchor_jd->ms = anchor_jd->ms - NO_OF_MS_IN_A_DAY;
                }


              /* I am close. Get the real time. */
              do
                {
                  anchor_jd->day = anchor_jd->day + timestep_bkw_jd->day;
                  anchor_jd->ms = anchor_jd->ms + timestep_bkw_jd->ms;

                  if ( anchor_jd->ms < 0 )
                    {
                      anchor_jd->day = anchor_jd->day - 1;
                      anchor_jd->ms = anchor_jd->ms + NO_OF_MS_IN_A_DAY;
                    }

                } while(! (anchor_jd->day < start_jd->day  ||  (anchor_jd->day == start_jd->day && anchor_jd->ms < start_jd->ms)) );


              /* I jumped one time too much. Move forward.  */
              anchor_jd->day = anchor_jd->day + timestep_jd->day;
              anchor_jd->ms = anchor_jd->ms + timestep_jd->ms;

              if (anchor_jd->ms >= NO_OF_MS_IN_A_DAY)
                {
                  anchor_jd->day = anchor_jd->day + 1;
                  anchor_jd->ms = anchor_jd->ms - NO_OF_MS_IN_A_DAY;
                }

              /* anchor_jd now has the true trigger time. Copy it to anchor. */
              julian2date(anchor_jd,anchor);

            default:

              ; //Should never happen. Don't return NULL just yet as deallocateDateTime(anchor) still needs to be called.

            }

          /* Copy the contents to target. */
          replaceDatetime(anchor,first_trigger_dt);


          /* Cleanup. */

          deallocateDateTime(anchor);
          deallocateJulianDay(start_jd);
          deallocateJulianDay(anchor_jd);
          deallocateJulianDelta(timestep_jd);
          deallocateTimeDelta(timestep_bkw);
          deallocateJulianDelta(timestep_bkw_jd);

          deallocateJulianDelta(timestep_fastfwd_jd);
          deallocateJulianDelta(timestep_fastbkw_jd);

          /* And return. */
          return first_trigger_dt;
        }
      else
        {
          /* NOTE: general case updates by explicit adding. */

          /* Get anchor. ref + offset is the real anchor. */
          struct _datetime *anchor = constructAndCopyDateTime(ref_dt);
          addTimeDeltaToDateTime(ref_dt, offset, anchor); // Note: If offset is null, no addition takes place.

          while (compareDatetime(anchor, start_dt) == less_than)
            {
              addTimeDeltaToDateTime(anchor, timestep, anchor);
            }
          /* Copy the contents to target. */
          replaceDatetime(anchor,first_trigger_dt);

          /* Cleanup. */
          deallocateDateTime(anchor);

          /* And return. */
          return first_trigger_dt;
        }
    }
  else
    return NULL;
}


/**
  * @brief Get the Datetime when 'this' event will be triggered next.
  *
  * WARNING: The value returned is with-respect-to current_dt and not a true copy of triggerNextEventDateTime in the event data structure.
  *
  * @param  e
  *         A pointer to struct _event. This is the event being queried.
  *
  * @param  dt_return
  *         A pointer to struct _datetime. The next trigger datetime is copied here.
  *
  * @return dt_return
  *         A pointer to DateTime with next-trigger datetime.
  */

struct _datetime*
getTriggerNextEventAtDateTime(struct _event* e, struct _datetime* current_dt, struct _datetime* dt_return)
{
  if ((e != NULL) && (current_dt != NULL) && (dt_return != NULL))
    {
       /* If last date is defined, check if current_dt is ahead of event->eventLastDateTime.
          Return on NULL indicating no future trigger event scheduled.*/
      if (e->eventLastDateTime && (compareDatetime(current_dt,e->eventLastDateTime) == greater_than))
        return NULL;

      /* This if should evaluate to false normally. The only occasion when it does not is if the
         current_dt > triggerNextEventDateTime and isCurrentEventActive(.) has never been called  */
      if (
           (e->nextEventIsFirst)
           &&
           (compareDatetime(current_dt,e->triggerNextEventDateTime) == greater_than)
         )
        {
              struct _timedelta* modulo_td = newTimeDelta("PT00.000S");
              // Get the first trigger time and return (WARNING: Do not update e->triggerNextEventDateTime here!).
              moduloTimeDeltaFromDateTime(e->triggerNextEventDateTime, e->eventInterval, current_dt, modulo_td);
              addTimeDeltaToDateTime(current_dt,modulo_td,dt_return);
              // Cleanup.
              deallocateTimeDelta(modulo_td);
              return dt_return;
        }
      // else ...
      dt_return = replaceDatetime(e->triggerNextEventDateTime, dt_return);
      return dt_return;
    }
  else
    return NULL;
}

/**
  * @brief Get the Datetime when 'this' event was triggered last.
  *
  * @param  e
  *         A pointer to struct _event. This is the event being queried.
  *
  * @param  dt_return
  *         A pointer to struct _datetime. The last trigger datetime is copied here.
  *
  * @return dt_return
  *         A pointer to DT with previous-trigger datetime.
  */
/* NOTE: If the event was never tiggered, default value of 0-01-01T00:00:00.000 is returned. */
struct _datetime*
getTriggeredPreviousEventAtDateTime(struct _event* e, struct _datetime* dt_return)
{
  if ((e != NULL )&& (dt_return != NULL) )
    {
      /* No trigger ever happened? */
      if (e->nextEventIsFirst)
        return NULL;

      replaceDatetime(e->triggeredPreviousEventDateTime, dt_return);
      return dt_return;
    }
  else
    return NULL;
}


struct _event*
getNextEventInGroup(struct _event* e)
{
  if (e != NULL)
    return e->nextEventInGroup;
  else
    return NULL;
}

int64_t
getEventId(struct _event* e)
{
  if (e != NULL)
    return e->eventId;
  else
    return 0;
}

char*
getEventName(struct _event* e, char* ename)
{
  if ((e != NULL) && (ename != NULL))
    {
      strncpy(ename, e->eventName, MAX_EVENTNAME_STR_LEN-1);
      ename[MAX_EVENTNAME_STR_LEN-1] = '\0';
      return ename;
    }
  else
    return NULL;
}

struct _datetime*
getEventReferenceDateTime(struct _event* e)
{
  if (e != NULL)
    return e->eventReferenceDateTime;
  else
    return NULL;
}

struct _datetime*
getEventFirstDateTime(struct _event* e)
{
  if (e != NULL)
    return e->eventFirstDateTime;
  else
    return NULL;
}

struct _datetime*
getEventLastDateTime(struct _event* e)
{
  if (e != NULL)
    return e->eventLastDateTime;
  else
    return NULL;
}

struct _timedelta*
getEventInterval(struct _event* e)
{
  if (e != NULL)
    return e->eventInterval;
  else
    return NULL;
}

bool
getNextEventIsFirst(struct _event* e)
{
  if (e != NULL)
    return e->nextEventIsFirst;
  else
    return false;
}

bool
getEventisFirstInDay(struct _event* e)
{
  if (e != NULL)
    return e->eventisFirstInDay;
  else
    return false;
}

bool
getEventisFirstInMonth(struct _event* e)
{
  if (e != NULL)
    return e->eventisFirstInMonth;
  else
    return false;
}

bool
getEventisFirstInYear(struct _event* e)
{
  if (e != NULL)
    return e->eventisFirstInYear;
  else
    return false;
}

bool
getEventisLastInDay(struct _event* e)
{
  if (e != NULL)
    return e->eventisLastInDay;
  else
    return false;
}

bool
getEventisLastInMonth(struct _event* e)
{
  if (e != NULL)
    return e->eventisLastInMonth;
  else
    return false;
}

bool
getEventisLastInYear(struct _event* e)
{
  if (e != NULL)
    return e->eventisLastInYear;
  else
    return false;
}
