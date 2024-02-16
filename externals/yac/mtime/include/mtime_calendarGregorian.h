/*! \cond PRIVATE */
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_calendarGregorian.h
 * 
 * @brief Convert julian dates to Gregorian-Calendar dates and vice-versa.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 * @note TODO- On Luis.
 */


#ifndef _MTIME_CALENDARGREGORIAN_H
#define _MTIME_CALENDARGREGORIAN_H

struct _datetime;
struct _julianday;

struct _datetime *
getDateGregorianFromJulian(struct _julianday *jd, struct _datetime* gd);
struct _julianday *
getJulianFromDateGregorian(struct _datetime* gd, struct _julianday *jd);

/**
 * @}
 */
#endif

/*! \endcond */
