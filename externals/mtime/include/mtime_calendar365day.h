// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/*! \cond PRIVATE */
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_calendar365day.h
 *
 * @brief Convert julian dates to Calendar-with-365-days dates and vice-versa.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 * @note Calendar-with-365-days has a non-leap year characteristic and has 365 days in a calendar year.
 * Number of days in the 12 months is given by: {  31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }.
 * Julian Day (0,0) corresponds to Calendar-with-365-day's day 0-01-01T12:00:00.000Z.
 *
 */

#ifndef _MTIME_CALENDAR365DAY_H
#define _MTIME_CALENDAR365DAY_H

struct _datetime;
struct _julianday;

struct _datetime *getDate365FromJulian(struct _julianday *jd, struct _datetime *gd);
struct _julianday *getJulianFromDate365(struct _datetime *gd, struct _julianday *jd);

/**
 * @}
 */
#endif

/*! \endcond */
