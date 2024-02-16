/*! \cond PRIVATE */
/**
 * @file mtime_eventList.c
 *
 * @brief _eventList data structure supports 'Event-groups'; Each group stores multiple events as a link list of nodes.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note All functions in this file are internal.
 */


// FIXME: This FILE is CALLed mtime_eventList.c, but the STRUCT as well the fortran TYPE is CALLed eventGroup!

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mtime_eventList.h"

#include "mtime_eventHandling.h"


/* Deallocate all Event's in Group. */
void
deallocateEventsInGroup(struct _eventGroup* eg)
{
  if (eg != NULL )
    { 
      struct _event* conductor = eg->rootEvent;
      /* Loop over and deallocate all non NULL events in the list. */
      while ((eg->rootEvent != NULL ))
        {
          eg->rootEvent = eg->rootEvent->nextEventInGroup;
	  deallocateEvent(conductor);
          conductor = eg->rootEvent;
        }
    }
}


/* Add new Event to Group. Attach event at the end of list. */
bool
addNewEventToGroup(struct _event* ev, struct _eventGroup* eg)
{
  if ( (ev != NULL ) && (eg != NULL) ){
  if(eg->rootEvent != NULL)
    {
      /* Non empty Group. */

      struct _event* conductor = eg->rootEvent;
      /* Loop over and get to the end of lisr. */
      while(conductor->nextEventInGroup != NULL)
        {
          conductor = conductor->nextEventInGroup;
        }

      /* Add new node. */
      conductor->nextEventInGroup = ev;
    }
  else
    {
      /* This is an empty Group.*/
      eg->rootEvent = ev;
    }

  return true;
}
return false;
}


/* Remove node with a given name */
bool
removeEventWithNameFromGroup(char* nodeName, struct _eventGroup* eg)
{
  if ( (nodeName != NULL ) && (eg != NULL) ){
  struct _event* conductor_fwd = eg->rootEvent;
  struct _event* conductor_bkd = eg->rootEvent;

  if(conductor_fwd == NULL)
    {
      /*Empty group. Can't delete a ghost!*/
      return false;
    }

  /* Delete root. */
  if(!strcmp(conductor_fwd->eventName,nodeName))
    {
      /* Set root to root->next. */
      eg->rootEvent = eg->rootEvent->nextEventInGroup;
      /* Deallocate root. */
      deallocateEvent(conductor_bkd);
      return true;
    }

  if(eg->rootEvent->nextEventInGroup != NULL)
    {
      /* There are atleast 2 objects in the Group. Restart from root->next. */
      conductor_fwd = conductor_fwd->nextEventInGroup;

      /* Loop over and find a name match. */
      while(strcmp(conductor_fwd->eventName,nodeName))
        {
          conductor_fwd = conductor_fwd->nextEventInGroup;
          conductor_bkd = conductor_bkd->nextEventInGroup;

	  /* No name match. Can't delete. */
          if(conductor_fwd == NULL)
            {
              return false;
            }
        }
    }
  else
    {
      /* Root node, the only node, was not a match */
      return false;
    }

  /* The while exited with a match. Delete conductor_fwd. */
  conductor_bkd->nextEventInGroup = conductor_fwd->nextEventInGroup;
  deallocateEvent(conductor_fwd);

  return true;
}
else
return false;
}

/*! \endcond */
