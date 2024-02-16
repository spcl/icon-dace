// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/*! \cond PRIVATE */
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_calendar360day.h
 * 
 * @brief Convert julian dates to Calendar-with-360-days dates and vice-versa.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 * @note Calendar-with-360-days has 30 days in 'every' month and hence 360 days in every calendar year. 
 *       Also, Julian Day (0,0) corresponds to Calendar-with-360-day's day 0-01-01T12:00:00.000Z.
 */

#ifndef _MTIME_CALENDAR360DAY_H
#define _MTIME_CALENDAR360DAY_H

struct _datetime;
struct _julianday;

struct _datetime *
getDate360FromJulian(struct _julianday *jd, struct _datetime* gd);
struct _julianday *
getJulianFromDate360(struct _datetime* gd, struct _julianday *jd);

/**
 * @}
 */
#endif

/*! \endcond */
