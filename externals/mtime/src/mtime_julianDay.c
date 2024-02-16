// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/**
 * @file mtime_julianDay.c
 *
 * @brief Julian Day Calendar and some operations supported on julian dates.
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
#include <inttypes.h>


#include"mtime_calendar.h"
#include"mtime_julianDay.h"

#define JULIAN_DAY_UPPER_BOUND_360 		 773094113279
#define JULIAN_DAY_LOWER_BOUND_360              -773094113281
#define JULIAN_DAY_UPPER_BOUND_365 		 783831531519
#define JULIAN_DAY_LOWER_BOUND_365              -783831531521
#define JULIAN_DAY_UPPER_BOUND_GREGORIAN 	 784354017364
#define JULIAN_DAY_LOWER_BOUND_GREGORIAN 	-784350575246


/**
 * @brief Construct new Julian-Date object.
 *
 * @param  _day
 *         An "int64_t" value denoting the Day part of JD.
 * @param  _ms
 *         An "int64_t" value denoting the milli-second part of Date.
 * @return jd
 *         A pointer to a initialized Julian-Date object. 
 */

struct _julianday*
newJulianDay(int64_t _day, int64_t _ms)
{
  /* Get Calendar type. */
  calendarType ct = getCalendarType();
  if(ct)
  {
  struct _julianday* jd = (struct _julianday*)calloc(1,sizeof(struct _julianday));
  if (jd == NULL )
    return NULL ;

  /* Return on NULL if out of bound. Bounds are defined according to calendar type 
     and common bound on year; i.e, (year > 2147483647) and (year < -2147483648). 
  */
  switch (ct)
    {
      case PROLEPTIC_GREGORIAN:
        if(   
            (	
		(_day > JULIAN_DAY_UPPER_BOUND_GREGORIAN) 
		|| 
		(	
			(_day == JULIAN_DAY_UPPER_BOUND_GREGORIAN) 
			&& 
			(_ms >= NO_OF_MS_IN_HALF_DAY)
		)
	    )
            || 
            (
		(_day < JULIAN_DAY_LOWER_BOUND_GREGORIAN) 
		|| 
		(
			(_day == JULIAN_DAY_LOWER_BOUND_GREGORIAN) 
			&& 
			(_ms < NO_OF_MS_IN_HALF_DAY)
		)
	    )
          )
          return NULL;
        break;
      case YEAR_OF_365_DAYS:
        if( 
            (
		(_day > JULIAN_DAY_UPPER_BOUND_365) 
		|| 
		(
			(_day == JULIAN_DAY_UPPER_BOUND_365)
			&&
			(_ms >= NO_OF_MS_IN_HALF_DAY)
		)
	    )
            ||
            (
		(_day < JULIAN_DAY_LOWER_BOUND_365) 
		|| 
		(	
			(_day == JULIAN_DAY_LOWER_BOUND_365) 
			&& 
			(_ms < NO_OF_MS_IN_HALF_DAY)
		)
	    )
          ) 
	  return NULL;
        break;
      case YEAR_OF_360_DAYS:
        if( 
            (
		(_day > JULIAN_DAY_UPPER_BOUND_360) 
		|| 
		(
			(_day == JULIAN_DAY_UPPER_BOUND_360) 
			&& 
			(_ms >= NO_OF_MS_IN_HALF_DAY)
		)
	    )
            ||
	    (
		(_day < JULIAN_DAY_LOWER_BOUND_360) 
		|| 
		(
			(_day == JULIAN_DAY_LOWER_BOUND_360) 
			&& 
			(_ms < NO_OF_MS_IN_HALF_DAY)
		)
	    )
          )
          return NULL;
        break;
      case CALENDAR_NOT_SET:
      default:
        //Should never be here.
        return NULL;
    }

  jd->day = _day;
  jd->ms = _ms;

  return jd;
  }
  else
    return NULL;
}

/**
 * @brief COPY a JulianDay object.
 *
 * Routine replaceJulianDay copies the contents of source JulianDay into a destination JulianDay object.
 *
 * @param  jdsrc         
 *         A pointer to struct _julianday. Copy "FROM" JulianDay object.
 *
 * @param  jddest
 *         A pointer to struct _julianday. Copy "TO" JulianDay object.
 *
 * @return jddest
 *         A pointer to 'copied' JulianDay Object.
 */

struct _julianday*
replaceJulianday(struct _julianday* jdsrc, struct _julianday* jddest)
{
  if ((jddest != NULL) && (jdsrc != NULL)) {

    jddest->day = jdsrc->day;
    jddest->ms  = jdsrc->ms;

    return jddest;
  }
  else
    return NULL;
}

/*! \cond PRIVATE */
/* Internal data structure. */
/* Create a new Julian Delta. */
/* Note: A negative juliandelta is represented in the following way: -P01DT00.500S  = jd2->sign = '-', jd2->day = -30, jd2->ms = -500. */

struct _juliandelta*
newJulianDelta(char _sign, int64_t _day, int64_t _ms)
{
  if (_sign != '+' && _sign != '-')
  {
    return NULL;
  }
  if (getCalendarType())
    {
      struct _juliandelta* jd = (struct _juliandelta *)calloc(1,sizeof(struct _juliandelta));
      if (jd == NULL )
        return NULL ;

      imaxdiv_t refactor_ms = imaxdiv((intmax_t) _ms, (intmax_t) NO_OF_MS_IN_A_DAY);

      jd->sign = _sign;
      jd->day = _day + refactor_ms.quot;
      jd->ms = refactor_ms.rem;
      
      return jd;
    }
  else
    return NULL;
}
/*! \endcond */


/**
 * @brief Destructor of Julian-Date.
 *
 * @param  jd
 *         A pointer to struct _julianday. jd is deallocated.
 */

void
deallocateJulianDay(struct _julianday* jd)
{
  if (jd != NULL )
    {
      free(jd);
      jd = NULL;
    }
}


/* Destructor of JulianDelta.*/

void
deallocateJulianDelta(struct _juliandelta* jd)
{
  if (jd != NULL )
    {
      free(jd);
      jd = NULL;
    }
}

struct _julianday*
addJulianDeltaToJulianDay(struct _julianday* jd1, struct _juliandelta* jd2, struct _julianday* jd)
{
  if ((jd1 != NULL ) && (jd2 != NULL) && (jd != NULL) ){

  jd->day = jd1->day + jd2->day;
  jd->ms = jd1->ms + jd2->ms;

  if (jd->ms >= NO_OF_MS_IN_A_DAY)
    {
      jd->day = jd->day + 1;
      jd->ms = jd->ms - NO_OF_MS_IN_A_DAY;
    }
  return jd;
}
else
return NULL;
}


struct _julianday*
substractJulianDeltaFromJulianDay(struct _julianday* jd1, struct _juliandelta* jd2, struct _julianday* jd)
{
  if ((jd1 != NULL )&& (jd2 != NULL) && (jd != NULL) ){

  /* substractJulianDelta 'looks like' addJulianDelta but it is not. 
     Usually, when substractJulianDelta is called, jd2 will have every element with a negative sign. 
     This will ensure substraction in the proper sense. */

  jd->day = jd1->day + jd2->day;
  jd->ms = jd1->ms + jd2->ms;

  if ( jd->ms < 0 )
    {
      jd->ms = NO_OF_MS_IN_A_DAY + jd->ms;
      jd->day = jd->day - 1;
    }

  return jd;
}
else
  return NULL;
}

int
compareJulianDay(struct _julianday* jd1, struct _julianday* jd2)
{
  int iret = compare_error;
    
  if (jd1->day == jd2->day)
    {
      if (jd1->ms == jd2->ms)
        {
          iret = equal_to;
        }
      else if (jd1->ms > jd2->ms)
        {
          iret = greater_than;
          
        }
      else
        {
          iret = less_than;
        }
    }
  else if (jd1->day < jd2->day)
    {
      iret = less_than;
    }
  else
    {
      iret = greater_than;
    }
    
  return iret;
}

/*! \cond PRIVATE */
/*Internal function.*/
/* Return jd1 >= jd2 */

bool
internal_compareJulianDay(struct _julianday* jd1, struct _julianday* jd2)
{
  if (jd1->day > jd2->day)
    return true;
  else if (jd2->day > jd1->day)
    return false;
  else if (jd1->ms > jd2->ms)
    return true;
  else if (jd2->ms > jd1->ms)
    return false;
  else
    return true;
}
/*! \endcond */

/**
 * @brief Get the JulianDelta between two JulianDays as (jd1-jd2).
 *
 * @param  jd1
 *         A pointer to struct _julianday. 
 *
 * @param  jd2
 *         A pointer to struct _julianday.
 *
 * @param  jd
 *         A pointer to struct _juliandelta. Copy the result of (jd1 - jd2) in jd.
 *
 * @return jd
 *         A pointer to Juliandelta containing the result of substraction.
 */

struct _juliandelta*
substractJulianDay(struct _julianday* jd1, struct _julianday* jd2, struct _juliandelta* jd)
{
  if ((jd1 != NULL )&& (jd2 != NULL) && (jd != NULL) ){

  if ( internal_compareJulianDay(jd1,jd2) == true )
    {
      /*JD1 >= JD2*/

      jd->sign = '+';
      jd->day = jd1->day - jd2->day;
      jd->ms = jd1->ms - jd2->ms;

      if(jd->ms < 0)
        {
          jd->ms = NO_OF_MS_IN_A_DAY + jd->ms;
          jd->day = jd->day - 1;
        }
    }
  else
    {
      /*  JD2 > JD1 */

      /*Recursive call with jd1 and jd2 switched and then negate the values.*/
      jd = substractJulianDay(jd2,jd1,jd);

      jd->sign = '-';
      jd->day = (-1)*jd->day;
      jd->ms = (-1)*jd->ms;
    }

  return jd;
}
else
  return NULL;
}


/**
 * @brief Get Julian as a string.
 *
 * juliandayToString returns a string in Day.MS format.
 *
 * @param  jd
 *         A pointer to struct _julianday. The julian value to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where julian value is to be written.
 *
 * @return toStr
 *         A pointer to the string containing julian value.
 */

char*
juliandayToString(struct _julianday* jd, char* toStr)
  {
    if ((jd != NULL )&& ( toStr != NULL ) ){

      memset(toStr,'\0',MAX_JULIANDAY_STR_LEN);
 
      snprintf(toStr,MAX_JULIANDAY_STR_LEN,"%" PRIi64 ".%" PRIi64, jd->day, jd->ms );

      return toStr;
    }
    else
      return NULL;
  }

