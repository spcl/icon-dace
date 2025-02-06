// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_time.h
 *
 * @brief Time and some operations supported on Time.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology
 *
 * @date March 2013
 *
 */

#ifndef _MTIME_TIME_H
#define _MTIME_TIME_H

#include <stdint.h>

/**
 * @struct _time
 *
 * @brief struct _time containing usual time parameters.
 */

struct _time
{
  int hour;    ///< hour part of time.
  int minute;  ///< minute part of time.
  int second;  ///< second part of time.
  int ms;      ///< milli-second part of time. 0<=ms<=999.
};

struct _time *newTime(const char *time_string);

struct _time *newRawTime(int _hour, int _minute, int _second, int _ms);

struct _time *constructAndCopyTime(struct _time *t);

void deallocateTime(struct _time *t);

struct _time *replaceTime(struct _time *, struct _time *);

char *timeToString(struct _time *, char *);

char *timeToBasicString(struct _time *, char *);

char *timeToPosixString(struct _time *t, char *toStr, char *fmtString);

/**
 * @}
 */
#endif
