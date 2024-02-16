// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "mtime_calendar.h"
#include "mtime_datetime.h"
#include "mtime_date.h"
#include "mtime_time.h"
#include "mtime_timedelta.h"
#include "mtime_eventHandling.h"
#include "mtime_eventList.h"

#include "mtime_julianDay.h"

int
main(int argc, char *argv[])
{
  /* CRITICAL: Must be done first to use this lib.
     Init calanedar. 
     Arguments {PROLEPTIC_GREGORIAN,YEAR_OF_365_DAYS,YEAR_OF_360_DAYS} specify calendar type. */
  initCalendar(PROLEPTIC_GREGORIAN); // OR initCalendar(YEAR_OF_360_DAYS); OR initCalendar(YEAR_OF_365_DAYS);

  /* Check calendar type currently set. */
  char calendar_in_use[MAX_CALENDAR_STR_LEN]; // ALWAYS use *STR_LEN for size of strings.
  if(calendarToString(calendar_in_use))       // Good idea to check for return on NULL to capture errors being returned. 
    printf("Calendar in use: %s\n", calendar_in_use);
  else
    printf("ERROR: Calendar initialization failed.");



  /* Simulate a clock. Shows how addTimeDeltaToDateTime() works. */
  char current_time[MAX_DATETIME_STR_LEN];  // Again, we use *STR_LEN.
  char current_step[MAX_TIMEDELTA_STR_LEN]; // and here too.


  //PT strings can be directly used for initialization OR functions like getPTStringFrom* called.
  //struct _timedelta *timestep = newTimeDelta("PT12M"); 				 // time step : 720 s (12 minutes).        Use this OR
  struct _timedelta *timestep  = newTimeDelta(getPTStringFromSeconds(720,current_step)); // time step : 720 s (12 minutes).               this.
  if(!timestep)	// Strogly recommended to check if new* succeded!
    printf("ERROR: Could not initialize Model time step.");

  struct _datetime *start_date = newDateTime("2013-09-01T02:10:00.000");                 // start date.
  if (!start_date) // Strogly recommended to check if new* succeded!
    printf("ERROR: Could not initialize Model start time.");

  struct _datetime *stop_date  = newDateTime("2013-09-03T14:00:00.000");                 // stop  date.
  if (!stop_date) // Strogly recommended to check if new* succeded!
    printf("ERROR: Could not initialize Model stop time.");

  struct _datetime *model_time = newDateTime("2013-09-01T02:10:00.000");	         // current date.
  if (!model_time) // Strogly recommended to check if new* succeded!
    printf("ERROR: Could not initialize Model (current) time.");

  printf("Model time step: %s\n", timedeltaToString(timestep, current_step));
  printf("Model start time: %s\n\n", datetimeToString(model_time, current_time));

  do
    {
      printf("Model time: %s\n", datetimeToString(model_time, current_time)); /* Display. */
      addTimeDeltaToDateTime(model_time, timestep, model_time); /* Increment time by timestep. */
    } while (compareDatetime(model_time, stop_date) <= 0);

  printf("\nModel Finished. \n Current time after Exit: %s\n", datetimeToString(model_time, current_time));


  printf("\n\n***************************************************************************************************************\n\n");


  /* Create event and test events addition to event groups. */
  printf("\nCreate event and test events addition to event groups.\n\n");

  
  struct _eventGroup* testAllowedGroup = newEventGroup("EventGroupToTestAllowedEvents");
  if (testAllowedGroup)
    printf("\nGroup Name - EventGroupToTestAllowedEvents:\n");

  /* Minimum requirement: Name and Timestep.
     If start data is supplied as NULL, start initialized to 0-01-01T00:00:00.000 
     If Reference date is not supplied, Reference date is initialized to start date (0-01-01T00:00:00.000 if start date itself is NULL).
     If end date is NULL, it is logically equivalent to infinity. The data-structure will not be allocated any memory and will be a NULL pointer.
  */

  struct _event* Ocean = newEvent("Ocean",NULL,NULL,NULL, "PT01M",NULL);
  if(Ocean && addNewEventToEventGroup(Ocean,testAllowedGroup))
    printf("\nCreated and Added event Ocean.");

  struct _event* InEarth = newEvent("InEarth", NULL, NULL, "2050-01-01T00:10:00.000", "P01Y",NULL);
  if(InEarth && addNewEventToEventGroup(InEarth,testAllowedGroup))
    printf("\nCreated and Added event InEarth.");

  struct _event* Systems = newEvent("Systems", NULL, "2050-01-01T00:05:00.000", NULL, "PT30.000S", NULL);
  if(Systems && addNewEventToEventGroup(Systems,testAllowedGroup))
    printf("\nCreated and Added event Systems.");

  struct _event* iWillFail1 = newEvent(NULL, "2050-01-01T00:06:00.000", "2050-01-01T00:05:00.000", "2050-01-01T00:10:00.000", "PT30.000S",NULL);
  if(!iWillFail1)
    printf("\niWillFail1 event can not be created. Name can not be NULL. Name must be unique. ");

  struct _event* iWillFail2 = newEvent("iWillFail2", "2050-01-01T00:06:00.000", "2050-01-01T00:05:00.000", "2050-01-01T00:10:00.000", NULL,NULL);
  if(!iWillFail2)
    printf("\niWillFail2 event can not be created. Time step can not be null. ");


  /* Create events and add events to an event group. */
  struct _eventGroup* someGroup = newEventGroup("SomeEvent");
  if(someGroup)
    printf("\n\nGroup Name - SomeEvent:\n"); 

  struct _event* Fisher = newEvent("Fisher",NULL,"2050-01-01T00:00:00.000",NULL, "PT01M",NULL);
  if(Fisher && addNewEventToEventGroup(Fisher,someGroup))
    printf("\nCreated and Added event Fisher.");

  struct _event* Pearson = newEvent("Pearson", "2050-01-01T00:00:45.000", "2050-01-01T00:05:00.000", "2050-01-01T00:10:00.000", "PT30.000S",NULL);
  if(Pearson && addNewEventToEventGroup(Pearson,someGroup))
    printf("\nCreated and Added event Pearson.");

  /* Event handle and name can be different.  */
  struct _event* Cox = newEvent("Cox", "2050-01-01T00:16:00.000", "2050-01-01T00:05:00.000", "2050-01-01T00:10:00.000", "PT30.000S",NULL);
  if(Cox && addNewEventToEventGroup(Cox,someGroup))
    printf("\nCreated and Added event Cox.");
  

  /* List all events in a group. */
  printf("\n\n\nList all events in the group ");
  char event_Gname[MAX_GROUPNAME_STR_LEN];
  printf("%s:\n",  getEventGroupName(someGroup, event_Gname));
  char event_name[MAX_EVENTNAME_STR_LEN];
  //Start at root of Group and loop over.
  struct _event* e = getEventGroupRootEvent(someGroup);
  while ( e != NULL )
  {
    printf("\nEvent: %s",eventToString(e,event_name));
    e = e->nextEventInGroup;
  }



  /* Create a dummy clock and poll all events in a group to see if some event is active. */
  printf("\n\n\nGet Next Event trigger times:\n");
  //struct _timedelta* dummy_current_time_step 	  = newTimeDelta("PT01.000S"); 					// Use this OR
  struct _timedelta*   dummy_current_time_step    = newTimeDelta(getPTStringFromSeconds(1,current_step));	//     this.
  struct _datetime*    dummy_current_date         = newDateTime("2050-01-01T00:00:00.000");
  struct _datetime*    dummy_current_stop_date 		  = newDateTime("2050-01-02T00:13:30.000");
  if (!(dummy_current_time_step && dummy_current_date && dummy_current_stop_date)) // We said this before, always check when calling new*
    printf("\nDummy initializations failed.\n");

  /*  Test cases explaining the relationship between eventReferenceDateTime (Anchor), eventFirstDateTime (Start) and eventLastDateTime (Stop). */
  char trigger_time[MAX_DATETIME_STR_LEN];
  struct _datetime* next_trigger_dt             = newDateTime(initDummyDTString);
  /* Fisher: Ref date = NULL. => First trigger at start date. */ 
  getTriggerNextEventAtDateTime(Fisher,dummy_current_date,next_trigger_dt);
  datetimeToString(next_trigger_dt, trigger_time); 
  printf("\nFirst trigger of Event Fisher: %s  (Expect - 2050-01-01T00:00:00.000)",trigger_time);
  /* Pearson: Ref date != NULL. First trigger at nearest trigger DT >= start date such that trigger only allowed at Anchor + N * interval. */
  getTriggerNextEventAtDateTime(Pearson,dummy_current_date,next_trigger_dt);
  datetimeToString(next_trigger_dt, trigger_time);
  printf("\nFirst trigger of Event Pearson: %s (Expect - 2050-01-01T00:05:15.000)",trigger_time);
  /* Cox: Ref date != NULL. Even if anchor date is into the future, the above rule applies (i.e N can be negative (zero as well))*/
  getTriggerNextEventAtDateTime(Cox,dummy_current_date,next_trigger_dt);
  datetimeToString(next_trigger_dt, trigger_time);
  printf("\nFirst trigger of Event Cox: %s (Expect - 2050-01-01T00:05:00.000)",trigger_time);


  /* An example of an acesible API: Check if events scheduled trigger is first. */
  if (getNextEventIsFirst(Fisher))
    printf("\n\nConfirming scheduled trigger of Cox at %s is first! \n\n",trigger_time);

 
  printf("\n\n*******************************************************************************************************\n\n\n\n");

  printf("Run a Clock and poll all events in Group. Report if event is active:\n");
  //Timer.
  while  (compareDatetime(dummy_current_date,dummy_current_stop_date) <= 0)
    {
      printf("\nSystem time: %s", datetimeToString(dummy_current_date, current_time));
      //Start at root of Group.
      e = someGroup->rootEvent;
      //Loop over all events in an event group.
      while ( e != NULL )
        {
	  printf("\nCheck event %s .... ",eventToString(e,event_name));
          if (isCurrentEventActive(e,dummy_current_date, NULL, NULL))
            {
              printf("\n\n %s is active. ",eventToString(e,event_name));
              /* The following flags make sense IFF the events are CURRENTLY in a triggered state 
                 i.e isCurrentEventActive(.) is true. Else, the behavior is undefined.*/
              if (getEventisFirstInDay(e))
                printf("This is the FIRST TRIGGER of %s THIS DAY", event_name); 
              if (getEventisFirstInMonth(e))               
                printf(", THIS MONTH ");
              if (getEventisFirstInYear(e))               
                printf("and THIS YEAR.");
              if (getEventisLastInDay(e))
                printf("This is the LAST TRIGGER of %s THIS DAY", event_name);
              if (getEventisLastInMonth(e))
                printf(", THIS MONTH ");
              if (getEventisLastInYear(e))
                printf("and THIS YEAR.");
              printf("\n");
            }
	  /* reevaluate to check if that works */
	  if (isCurrentEventActive(e,dummy_current_date, NULL, NULL))
	    {
	      printf("\n %s is expected to be still active. \n\n",eventToString(e,event_name));
	    }
          e = e->nextEventInGroup;
	}
      addTimeDeltaToDateTime(dummy_current_date,dummy_current_time_step,dummy_current_date);
    }
   

  printf("\n\n\n*****************************************************\n\n");

  /* Now that the clock has crossed the End times (if defined), we expect to recieve NULL getTriggerNextEventAtDateTime(.). If not defined,  
     the next trigger should be returned. */
  printf("\n\n\n");
  printf("Now that the clock has crossed the End times (if defined),\n"); 
  printf("  we expect to recieve NULL on getTriggerNextEventAtDateTime(.),\n");
  printf("  ff not defined, the next trigger should be returned:\n");
  
  if(getTriggerNextEventAtDateTime(Fisher,dummy_current_date,next_trigger_dt))
    {
      datetimeToString(next_trigger_dt, trigger_time);
      printf("\nNext trigger of Event Fisher: %s  (Expect - 2050-01-01T00:14:00.000)",trigger_time);   
    }
  else
     printf("\nNext trigger of Event Fisher: NULL");

  if(getTriggerNextEventAtDateTime(Pearson,dummy_current_date,next_trigger_dt))
    {
      datetimeToString(next_trigger_dt, trigger_time);
      printf("\nNext trigger of Event Pearson: %s ",trigger_time);
    }
  else
    printf("\nNext trigger of Event Pearson: NULL  (Expect - NULL)");

  if(getTriggerNextEventAtDateTime(Cox,dummy_current_date,next_trigger_dt))
    {
      datetimeToString(next_trigger_dt, trigger_time);
      printf("\nNext trigger of Event Cox: %s  (Expect - <empty>)\n",trigger_time);
    }
  else
    printf("\nNext trigger of Event Cox: NULL  (Expect - NULL)\n");

  /* Cleanup before exit. */
  deallocateEventGroup(someGroup); // Will clean up all events in the group.
  deallocateEventGroup(testAllowedGroup); // Will clean up all events in the group.
  deallocateTimeDelta(timestep);
  deallocateTimeDelta(dummy_current_time_step);
  deallocateDateTime(dummy_current_stop_date);
  deallocateDateTime(dummy_current_date);
  deallocateDateTime(model_time);
  deallocateDateTime(stop_date);
  deallocateDateTime(start_date);
  deallocateDateTime(next_trigger_dt);

  /* CRITICAL: Must be called only at the end or the lib will enter an inconsitent state.*/
  freeCalendar();
  return 0;
}
