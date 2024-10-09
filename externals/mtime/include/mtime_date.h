// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_date.h
 *
 * @brief Date and some operations supported on Date.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 */

#ifndef _MTIME_DATE_H
#define _MTIME_DATE_H

#include "mtime_calendar.h"

/**
 * @struct _date
 *
 * @brief struct _date containing usual date parameters.
 */

struct _date
{
  int64_t year;  ///< Year of date. Can be both positive and negative.
  int month;     ///< Month of date.
  int day;       ///< day of date.
};

struct _date *newDate(const char *ds);

struct _date *newRawDate(int64_t _year, int _month, int _day);

struct _date *constructAndCopyDate(struct _date *d);

void deallocateDate(struct _date *d);

/*! \cond PRIVATE */
compare_return_val compareDate(struct _date *, struct _date *);
/*! \endcond  */

struct _date *replaceDate(struct _date *, struct _date *);

char *dateToString(struct _date *, char *ds);

char *dateToBasicString(struct _date *, char *ds);

char *dateToPosixString(struct _date *d, char *toStr, char *fmtString);

/**
 * @}
 */
#endif
