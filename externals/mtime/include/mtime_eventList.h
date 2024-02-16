// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/*! \cond PRIVATE */
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_eventList.h
 *
 * @brief _eventList data structure supports 'Event-groups'; Each group stores multiple events as a link list of nodes.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 * @note All functions in this file are internal.
 */

#ifndef _MTIME_EVENTLIST_H
#define _MTIME_EVENTLIST_H

#include <stdbool.h> 

struct _eventGroup;
struct _event;

void
deallocateEventsInGroup(struct _eventGroup* eg);

bool
addNewEventToGroup(struct _event* ev, struct _eventGroup* eg);

bool
removeEventWithNameFromGroup(char* nodeName, struct _eventGroup* eg);

/**
 * @}
 */
#endif

/*! \endcond */
