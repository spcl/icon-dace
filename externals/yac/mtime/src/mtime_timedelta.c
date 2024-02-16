/**
 * @file mtime_timedelta.c
 * 
 * @brief TimeDelta and some operations supported on TimeDelta.
 * 
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>

#include "mtime_timedelta.h"

#include "mtime_julianDay.h"
#include "mtime_datetime.h"
#include "mtime_iso8601.h"



/* MIN allowed year : -2147483648
   MAX allowed year :  2147483647
   Upto MilliSecond resolution supported. */



/**
 * @brief Construct new TimeDelta using an ISO 8601 conforming string.
 *
 * @param  tds
 *         A pointer to char. The string should contain parameters with which TimeDelta is to be created.
 *
 * @return td
 *         A pointer to a filled TimeDelta. 
 */

struct _timedelta*
newTimeDelta(const char* tds)
{
  if ((tds != NULL) && (getCalendarType()))
    {
      char duration_type_flag = 0;
      /* Verify string is ISO8601 compliant. */
      struct iso8601_duration* isoDuration = new_iso8601_duration('+', 0, 0, 0, 0, 0, 0, 0);
      if (isoDuration == NULL )
        return NULL ;

      if (
           ((duration_type_flag = verify_string_duration(tds, isoDuration)) != DURATION_MATCH_STD)
           && 
           (duration_type_flag != DURATION_MATCH_LONG)
         )
        {
          deallocate_iso8601_duration(isoDuration);
          isoDuration = NULL;
          return NULL ;
        }

      /* Create TimeDelta object. */
      struct _timedelta* td = (struct _timedelta *)calloc(1,sizeof(struct _timedelta));
      if (td == NULL )
        {
          deallocate_iso8601_duration(isoDuration);
          isoDuration = NULL;
          return NULL ;
        }

      duration_type_flag == 2?(td->flag_std_form = 1):(td->flag_std_form = 0);
      /* IMPORTANT: Negative/Positive time delta is indicated using td->sign (-/+). year,month,day..etc are always positive integers or 0. */
      td->sign          = isoDuration->sign;
      td->year          = isoDuration->year;
      td->month         = isoDuration->month;
      td->day           = isoDuration->day;
      td->hour          = isoDuration->hour;
      td->minute        = isoDuration->minute;
      td->second        = isoDuration->second;
      td->ms            = isoDuration->ms;

      //Cleanup.
      deallocate_iso8601_duration(isoDuration);
      isoDuration = NULL;       

      return td;
    }
  else
    return NULL ;
}


/**
 * @brief Construct new TimeDelta using 'raw' numerical values.
 *
 * @param  _sign
 *         An char value denoting positive('+') or negative('-') TimeDelta. 
 * @param  _year
 *         An "int64_t" value denoting the year part of TimeDelta.
 * @param  _month
 *         An "int" value denoting the month part of TimeDelta.
 * @param  _day
 *         An "int" value denoting the day part of TimeDelta.
 * @param  _hour
 *         An "int" value denoting the hour part of TimeDelta.
 * @param  _minute
 *         An "int" value denoting the minute part of TimeDelta.
 * @param  _second
 *         An "int" value denoting the second part of TimeDelta.
 * @param  _ms
 *         An "int" value denoting the milli-second part of TimeDelta.
 *
 * @return td
 *         A pointer to a filled TimeDelta. 
 */

struct _timedelta *
newRawTimeDelta(char _sign, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms)
{
  struct _timedelta* td;
  
  td = (struct _timedelta *) calloc(1, sizeof(struct _timedelta));
  if (td == NULL )
    {
      return NULL ;
    }

  if ( (_day < 32) && (_second < 60) && (_minute < 60) && (_hour < 24) )
    {
      td->flag_std_form = true;
    }
  else
    {
      td->flag_std_form = false;
    }
  
  td->sign   = _sign;
  td->year   = _year;
  td->month  = _month;
  td->day    = _day;
  td->hour   = _hour;
  td->minute = _minute;
  td->second = _second;
  td->ms     = _ms;
  
  return td;
}

/**
 * @brief Copy the values and construct a new TimeDelta.
 *
 * @param  td
 *         A pointer to struct _timedelta. Values of td are used to initialize the new timedelta being created. 
 *
 * @return _td
 *         A pointer to an initialized TimeDelta object. 
 */

struct _timedelta*
constructAndCopyTimeDelta(struct _timedelta* td)
{
  if ( td != NULL )
    {
      struct _timedelta* td_new;
      td_new = (struct _timedelta *) calloc(1, sizeof(struct _timedelta));
      if (td_new == NULL )
        {
          return NULL ;
        }

      td_new->flag_std_form = td->flag_std_form;
      td_new->sign          = td->sign;
      td_new->year          = td->year;
      td_new->month         = td->month;
      td_new->day           = td->day;
      td_new->hour          = td->hour;
      td_new->minute        = td->minute;
      td_new->second        = td->second;
      td_new->ms            = td->ms;

      return td_new;
    }
  else
    return NULL;
}


/**
 * @brief Destructor of TimeDelta.
 *
 * @param  td
 *         A pointer to struct _timedelta. td is deallocated.
 */

void
deallocateTimeDelta(struct _timedelta* td)
{
  if (td != NULL )
    {
      free(td);
      td = NULL;
    }
}


/**
 * @brief Compare two timedelta and return (td1 > td2) OR (td1 = td2) OR (td1 < td2).
 *
 * @param  td1
 *         A pointer to struct _timedelta.
 *
 * @param  td2
 *         A pointer to struct _timedelta.
 *
 * @return boolean
 *         if (td1 > td2), return greater_than. If (td1 == td2), return equal_to. If (td1 < td2), return less_than. 
 *         Return compare_error indicating error.
 */

compare_return_val
compareTimeDelta(struct _timedelta* td1, struct _timedelta* td2)
{

  if ((td1 != NULL) && (td2 != NULL))
    {

      /* If signs are different, we already know the order. */
      if (td1->sign == '+' && td2->sign == '-' )
	return greater_than;
      else if (td1->sign == '-' && td2->sign == '+' )
	return less_than;

      compare_return_val boolean = compare_error;

      /* 
       * add special case for seconds larger than 59 by creating
       * temporary/-ies, otherwise proceed with normal processing 
       */
      if ( (td1->second < 59) && td2->second < 59)  
	{	 
	  if(td1->year == td2->year)
	    {
	      if(td1->month == td2->month)
		{
		  if(td1->day == td2->day)
		    {
		      if(td1->hour == td2->hour)
			{
			  if(td1->minute == td2->minute)
			    {
			      if(td1->second == td2->second)
				{
				  if(td1->ms == td2->ms)
				    {
				      boolean = equal_to;
				    }
				  else if(td1->ms > td2->ms)
				    {
				      boolean = greater_than;
				    }
				  else if(td1->ms < td2->ms)
				    {
				      boolean = less_than;
				    }
				}
			      else if(td1->second > td2->second)
				{
				  boolean = greater_than;
				}
			      else if(td1->second < td2->second)
				{
				  boolean = less_than;
				}
			    }
			  else if(td1->minute > td2->minute)
			    {
			      boolean = greater_than;
			    }
			  else if(td1->minute < td2->minute)
			    {
			      boolean = less_than;
			    }       
			}
		      else if(td1->hour > td2->hour)
			{
			  boolean = greater_than;
			}
		      else if(td1->hour < td2->hour)
			{
			  boolean = less_than;
			}
		    }
		  else  if(td1->day > td2->day)
		    {
		      boolean = greater_than;
		    }
		  else  if(td1->day < td2->day)
		    {
		      boolean = less_than;                    
		    }
		}
	      else if (td1->month > td2->month)
		{
		  boolean = greater_than;
		}
	      else if (td1->month < td2->month)
		{
		  boolean = less_than;
		}
	    }
	  else if (td1->year > td2->year)
	    {
	      boolean = greater_than;
	    }
	  else if (td1->year < td2->year)
	    {
	      boolean = less_than;
	    }
	  }
      else
	{
	  struct _timedelta td1_tmp, td2_tmp;
	  if ( (td1->year == 0) && (td2->year == 0))
	    {
	      if ( (td1->month == 0) && (td2->month == 0))
		{
		  int second;
		  
		  td1_tmp.day    = td1->day;
		  td1_tmp.hour   = td1->hour;
		  td1_tmp.minute = td1->minute;
		  td1_tmp.second = td1->second;
		  td1_tmp.ms     = td1->ms;

		  second = td1_tmp.second;
		  td1_tmp.day    += second/86400;
		  td1_tmp.hour   += (second%86400)/3600;
		  td1_tmp.minute += ((second%86400)%3600)/60;
		  td1_tmp.second  = (((second%86400)%3600)%60);

		  td2_tmp.day    = td2->day;
		  td2_tmp.hour   = td2->hour;
		  td2_tmp.minute = td2->minute;
		  td2_tmp.second = td2->second;
		  td2_tmp.ms     = td2->ms;

		  second = td2_tmp.second;
		  td2_tmp.day    += second/86400;
		  td2_tmp.hour   += (second%86400)/3600;
		  td2_tmp.minute += ((second%86400)%3600)/60;
		  td2_tmp.second  = (((second%86400)%3600)%60);

		  if(td1_tmp.day == td2_tmp.day)
		    {
		      if(td1_tmp.hour == td2_tmp.hour)
			{
			  if(td1_tmp.minute == td2_tmp.minute)
			    {
			      if(td1_tmp.second == td2_tmp.second)
				{
				  if(td1_tmp.ms == td2_tmp.ms)
				    {
				      boolean = equal_to;
				    }
				  else if(td1_tmp.ms > td2_tmp.ms)
				    {
				      boolean = greater_than;
				    }
				  else if(td1_tmp.ms < td2_tmp.ms)
				    {
				      boolean = less_than;
				    }
				}
			      else if(td1_tmp.second > td2_tmp.second)
				{
				  boolean = greater_than;
				}
			      else if(td1_tmp.second < td2_tmp.second)
				{
				  boolean = less_than;
				}
			    }
			  else if(td1_tmp.minute > td2_tmp.minute)
			    {
			      boolean = greater_than;
			    }
			  else if(td1_tmp.minute < td2_tmp.minute)
			    {
			      boolean = less_than;
			    }       
			}
		      else if(td1_tmp.hour > td2_tmp.hour)
			{
			  boolean = greater_than;
			}
		      else if(td1_tmp.hour < td2_tmp.hour)
			{
			  boolean = less_than;
			}
		    }
		  else  if(td1_tmp.day > td2_tmp.day)
		    {
		      boolean = greater_than;
		    }
		  else  if(td1_tmp.day < td2_tmp.day)
		    {
		      boolean = less_than;                    
		    }
		}
	      else
		{
		  return compare_error;
		}
	    }
	  else
	    {
	      return compare_error;	      
	    }
	}
      
      if (td1->sign == '+' && td2->sign == '+' )
	return boolean;
      else if (td1->sign == '-' && td2->sign == '-' )
	return boolean*(-1);
    }
  
  return compare_error;
}


/**
 * @brief COPY a TimeDelta object.
 *
 * Routine replaceTimeDelta copies the contents of source TimeDelta into a Destination TimeDelta object.
 *
 * @param  tdsrc         
 *         A pointer to struct _timedelta. Copy "FROM" TimeDelta object.
 *
 * @param  tddest
 *         A pointer to struct _timedelta. Copy "TO" TimeDelta object.
 *
 * @return tddest
 *         A pointer to 'copied' TimeDelta Object.
 */

struct _timedelta*
replaceTimeDelta(struct _timedelta* tdsrc, struct _timedelta* tddest)
{
  if ((tddest != NULL) && (tdsrc != NULL)){

  tddest->flag_std_form = tdsrc->flag_std_form;

  tddest->sign          = tdsrc->sign;

  tddest->year          = tdsrc->year;
  tddest->month         = tdsrc->month;
  tddest->day           = tdsrc->day;
  tddest->hour          = tdsrc->hour;
  tddest->minute        = tdsrc->minute;
  tddest->second        = tdsrc->second;
  tddest->ms            = tdsrc->ms;

  return tddest;
}
else
  return NULL;
}


static 
struct _juliandelta*
localTimeDeltaToJulianDelta_NonStandardTimeDelta(struct _timedelta* td, struct _datetime* base_dt, struct _juliandelta* jd_return)
  {

    if ( td->sign == '-' )
      {
        /* Negative TimeDelta.  A negative juliandelta is represented in the following 
           way: -P01DT00.500S  = jd2->sign = '-', jd2->day = -30, jd2->ms = -500. */

        jd_return->sign = '-';

        /* Non standard Months not supported. */

        /* Day. */
        jd_return->day = (-1) * (int64_t)td->day;
        /* Rest. */
        jd_return->ms   = (-1) * (int64_t)td->hour   * NO_OF_MS_IN_A_HOUR
                        + (-1) * (int64_t)td->minute * NO_OF_MS_IN_A_MINUTE
                        + (-1) * (int64_t)td->second * NO_OF_MS_IN_A_SECOND
                        + (-1) * td->ms; /* No need to promote ms to int64_t; ms greater than 999 not supported even in non std long form. */

        while (jd_return->ms <= (-1) *NO_OF_MS_IN_A_DAY)
          {
            jd_return->ms = jd_return->ms + NO_OF_MS_IN_A_DAY;
            jd_return->day--;
          }

      }
    else if ( td->sign == '+' )
      {
        /* Positive TimeDelta. */

        jd_return->sign = '+';

        /* Non standard Months not supported. */

        /* Day. */
        jd_return->day = (int64_t)td->day;
        /* Rest. */
        jd_return->ms    = (int64_t)td->hour   * NO_OF_MS_IN_A_HOUR
                         + (int64_t)td->minute * NO_OF_MS_IN_A_MINUTE
                         + (int64_t)td->second * NO_OF_MS_IN_A_SECOND
                         + td->ms; /* No need to promote ms to int64_t; ms greater than 999 not supported even in non std long form. */

        while (jd_return->ms >= NO_OF_MS_IN_A_DAY)
          {
            jd_return->ms = jd_return->ms - NO_OF_MS_IN_A_DAY;
            jd_return->day++;
          }
      }
    else
      return NULL; /* ERROR: TD sign not defined. Should never happen. */

    return jd_return;
  }


static
struct _juliandelta*
localTimeDeltaToJulianDelta_StandardTimeDelta_CalType360OR365(struct _timedelta* td, struct _datetime* base_dt, struct _juliandelta* jd_return)
{
  /* No of days in a year. */
  int ndiny;

  /* Pointer to an array of month_specific_delta_in_months. */
  const int (*msdinm)[NO_OF_MONTHS_IN_A_YEAR+1];


  /* Initialize delta to 0.*/
  jd_return->day = 0;
  jd_return->ms = 0;

  switch (getCalendarType())
    {
      case YEAR_OF_365_DAYS:

      msdinm = monthSpecificDeltaInMonths365;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
      break;


      case YEAR_OF_360_DAYS:

      msdinm = monthSpecificDeltaInMonths360;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360;
      break;


      default:
      break;
    }

  if ( td->sign == '-' )
    {
      /* Negative TimeDelta.  A negative juliandelta is represented in the following 
         way: -P01DT00.500S  = jd2->sign = '-', jd2->day = -30, jd2->ms = -500. */

      jd_return->sign = '-';
      /* Year. */
      jd_return->day = (-1) * td->year * ndiny;
      /* Month. */
      jd_return->day = jd_return->day + (-1) * ( ndiny - (msdinm[base_dt->date.month - 1][NO_OF_MONTHS_IN_A_YEAR - td->month]) );
      /* Day. */
      jd_return->day = jd_return->day + (-1) * td->day;
      /* Rest. */
      jd_return->ms     = (-1) * (int64_t)td->hour   * NO_OF_MS_IN_A_HOUR
                        + (-1) * (int64_t)td->minute * NO_OF_MS_IN_A_MINUTE
                        + (-1) * (int64_t)td->second * NO_OF_MS_IN_A_SECOND
                        + (-1) * td->ms;
    }
  else if ( td->sign == '+' )
    {
      /* Positive TimeDelta. */

      jd_return->sign = '+';
      /* Year.  */
      jd_return->day = td->year * ndiny;
      /* Month. No of days in a TimeDelta depends on Base date.*/
      jd_return->day = jd_return->day + msdinm[base_dt->date.month - 1][ td->month ] ;
      /* Day. */
      jd_return->day = jd_return->day + td->day;
      /* Rest. */
      jd_return->ms     = (int64_t)td->hour   * NO_OF_MS_IN_A_HOUR
                        + (int64_t)td->minute * NO_OF_MS_IN_A_MINUTE
                        + (int64_t)td->second * NO_OF_MS_IN_A_SECOND
                        + td->ms;
    }
  else
    return NULL; /* ERROR: TD sign not defined. Should never happen. */

  return jd_return;
}


static
struct _juliandelta*
localTimeDeltaToJulianDelta_StandardTimeDelta_CalTypeGREGORIAN(struct _timedelta* td, struct _datetime* base_dt, struct _juliandelta* jd_return)
{
  /* Gregorian will have 366 days or 365 days depending on Leap year or Non-Leap year. */


  /* No of days in a year. */
  int ndiny;

  /* Pointer to an array of month_specific_delta_in_months. */
  const int (*msdinm)[NO_OF_MONTHS_IN_A_YEAR+1];


  /* Initialize delta to 0.*/
  jd_return->day = 0;
  jd_return->ms = 0;

      if ( td->sign == '+' )
        {
          jd_return->sign = '+';

          /* In final year, the crucial point is the month of february. */
          if (
                ( 
                        (testYearIsLeapYear(base_dt->date.year + td->year + 1) ) 
                        && 
                        (base_dt->date.month >= 3) 
                )
                ||
                ( 
                        (testYearIsLeapYear(base_dt->date.year + td->year)) 
                        && 
                        (base_dt->date.month < 3) 
                )
            )
            {
              /* If the (base year + delta year) is a year before a leap year and base month is >= 3 
               OR 
               base year + delta year is a leap year and month is < 3
               => An addition of leap-year specific delta for each month.
              */
              msdinm = monthSpecificDeltaInMonthsLeapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Otherwise addition of non-leap year specific month. */
              msdinm = monthSpecificDeltaInMonths365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          /* Fast-Fwd >= 400 */
          int64_t numberOf400YearPeriods = td->year / 400;
          int64_t i =  base_dt->date.year + numberOf400YearPeriods * 400;
          jd_return->day = jd_return->day + numberOf400YearPeriods * NO_OF_DAYS_IN_400_YEARS;


          /* The year from (target date - 399) to base_date + delta years - 1 */
          for (; i < base_dt->date.year + td->year; i++)
            {
              if (
                        ((testYearIsLeapYear(i + 1)) && (base_dt->date.month >= 3))
                        ||
                        ((testYearIsLeapYear(i)) && (base_dt->date.month < 3))
                 )
                {
                  /* If the next year is a leap year and month is >= 3 OR 
                   this year is a leap year and month is less than 3 
                   => delta of 1 year corresponds to 366 day julian delta. 
                  */
                  jd_return->day = jd_return->day + NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  jd_return->day = jd_return->day + NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }
            }

          jd_return->day = jd_return->day + msdinm[base_dt->date.month - 1][td->month] 
                                          + td->day;

          jd_return->ms =  (int64_t)td->hour   * NO_OF_MS_IN_A_HOUR 
                        +  (int64_t)td->minute * NO_OF_MS_IN_A_MINUTE 
                        +  (int64_t)td->second * NO_OF_MS_IN_A_SECOND 
                        + td->ms;
        }
      else if ( td->sign == '-' )
        {
          /* A negative juliandelta is represented in the following way: 
             -P01DT00.500S  = jd2->sign = '-', jd2->day = -30, jd2->ms = -500.  */

          jd_return->sign = '-';

          /* In final year, the crucial point is the month of february. */
          if (
                ((testYearIsLeapYear(base_dt->date.year - td->year - 1)) && (base_dt->date.month < 3))
                ||
                ((testYearIsLeapYear(base_dt->date.year - td->year)) && (base_dt->date.month >= 3)))
            {
              /* If the (base year - delta year) is a year after leap year and base month is < 3 
               OR 
               base year - delta year is a leap year and month is >= 3
               => A substraction of leap-year specific delta for each month.
              */
              msdinm = monthSpecificDeltaInMonthsLeapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Otherwise. */
              msdinm = monthSpecificDeltaInMonths365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          /* Fast-Fwd >= 400 */
          int64_t numberOf400YearPeriods = td->year / 400;
          int64_t i =  base_dt->date.year - numberOf400YearPeriods * 400;
          jd_return->day = jd_return->day - numberOf400YearPeriods * NO_OF_DAYS_IN_400_YEARS;


          /* The year from (target date + 399) to base_date - delta years + 1 */
          for (; i > base_dt->date.year - td->year; i--)
            {
              if (
                        ((testYearIsLeapYear(i - 1)) && (base_dt->date.month < 3))
                        ||
                        ((testYearIsLeapYear(i)) && (base_dt->date.month >= 3))
                 )
                {
                  /* If the previous year is a leap year and month is < 3 OR 
                   this year is a leap year and month is >= 3 
                   => delta of 1 year corresponds to 366 day julian delta. 
                   */
                  jd_return->day = jd_return->day - NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  jd_return->day = jd_return->day - NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }
            }
          jd_return->day = jd_return->day + (-1) * (ndiny - msdinm[base_dt->date.month - 1][NO_OF_MONTHS_IN_A_YEAR - td->month])
                                          + (-1) * td->day;

          jd_return->ms = (-1) *  (int64_t)td->hour * NO_OF_MS_IN_A_HOUR 
                        + (-1) *  (int64_t)td->minute * NO_OF_MS_IN_A_MINUTE 
                        + (-1) *  (int64_t)td->second * NO_OF_MS_IN_A_SECOND 
                        + (-1) * td->ms;

        }
      else
        return NULL; /* ERROR: Sign not set. Should never happen. */

   return jd_return;
}




/* Internal function. */
/* Converts TimeDelta value to a lib defined juliandelta value. Juliadelta depends on the calendar type. 
 Notice that TimeDelta is not uniquely defined but depends on the definition of corresponding DateTime 
 (Referred to as Anchor date in the following).

 The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
 positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
 Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
 a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
 so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
 anchor date month; a delta of 2 month corresponds to as many days as in the month before the
 anchor date month and the month before that and so on. 
 
 For eg, TimeDelta of P01M, when the 'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a
 Julian Delta of positive 28 days ( #days in February ) and 0 ms. Also, TimeDelta of P02M, when the 
 'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of positive 28+31 days
 ( #days in February PLUS #days in March) and 0 ms. Similarly, a TimeDelta of -P01M with 'Anchor' 
 DateTime of 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31 days 
 ( #days in January ) and 0 ms. Likewise, TimeDelta of -P02M with 'Anchor' DateTime of 
 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31+31 days ( #days in January
 PLUS #days in December) and 0 ms. The same logic is used for all calendar types (with corresponding 
 days in months according to current calendar type.)
 */

struct _juliandelta*
timeDeltaToJulianDelta(struct _timedelta* td, struct _datetime* base_dt, struct _juliandelta* jd_return)
{
  if ((td != NULL) && (base_dt != NULL) && (jd_return != NULL)){

    if ( !td->flag_std_form ) /* Non standard long string: PT2000S */
      {
        localTimeDeltaToJulianDelta_NonStandardTimeDelta(td, base_dt, jd_return);
      }
    else /* Standard Form eg. PT20M30S  */
      {
        switch (getCalendarType())
          {
            case YEAR_OF_365_DAYS:
            case YEAR_OF_360_DAYS:
        
            localTimeDeltaToJulianDelta_StandardTimeDelta_CalType360OR365(td, base_dt, jd_return);
            break;


            case PROLEPTIC_GREGORIAN:

            localTimeDeltaToJulianDelta_StandardTimeDelta_CalTypeGREGORIAN(td, base_dt, jd_return);
            break;


            default:
            return NULL; /* ERROR: Calendar type not set. */

          }

      }
    return jd_return;  
  }
  else
    return NULL;
}


/* Internal function. */
/* Converts a lib defined Julian delta value to a TimeDelta value. TimeDelta depends on the calendar type. 
 Notice that TimeDelta is not uniquely defined but depends on the definition of corresponding DateTime 
 (Referred to as Anchor date in the following).

 The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
 positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
 Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
 a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
 so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
 anchor date month; a delta of 2 month corresponds to as many days as in the month before the
 anchor date month and the month before that and so on. 

 For eg, TimeDelta of P01M, when the 'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a
 Julian Delta of positive 28 days ( #days in February ) and 0 ms. Also, TimeDelta of P02M, when the 
 'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of positive 28+31 days
 ( #days in February PLUS #days in March) and 0 ms. Similarly, a TimeDelta of -P01M with 'Anchor' 
 DateTime of 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31 days 
 ( #days in January ) and 0 ms. Likewise, TimeDelta of -P02M with 'Anchor' DateTime of 
 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31+31 days ( #days in January
 PLUS #days in December) and 0 ms. The same logic is used for all calendar types (with corresponding 
 days in months according to current calendar type).
*/

struct _timedelta*
julianDeltaToTimeDelta(struct _juliandelta* jd, struct _datetime* base_dt, struct _timedelta* td_return)
{
  if ((jd != NULL) && (base_dt != NULL) && (td_return != NULL)){

  int i = 0;
  /* No of days in a year. */
  int ndiny;

  /* Pointer to an array of no of days in month. */
  /* Pointer to an array of month_specific_delta_in_months. */
  const int (*msdinm)[NO_OF_MONTHS_IN_A_YEAR+1];

  /* We simply assume this routines returns a "normalized" timedelta. */
  td_return->flag_std_form = 1;

  /* Set parameter according to calender. */
  switch (getCalendarType())
    {
    case YEAR_OF_365_DAYS:

      msdinm = monthSpecificDeltaInMonths365;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
      break;

    case YEAR_OF_360_DAYS:

      msdinm = monthSpecificDeltaInMonths360;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360;
      break;

    case PROLEPTIC_GREGORIAN:
      /* Handle all Gregorian related code here. */

      /* Gregorian will have 366 days and 365 days depending on Leap year. */

      if ( jd->sign == '-' )
        {
          /* Negative TimeDelta.  A negative juliandelta is represented in the following 
             way: -P01DT00.500S  = jd2->sign = '-', jd2->day = -30, jd2->ms = -500. */  

          td_return->sign = '-';

          /* No of days in the final year */
          int64_t delta_final_year;
          int64_t days = (-1)*jd->day;
          /* Set counter to base year and then iterate backward to get to the final year.
             For each loop forward, increment year by 1.
          */
          int64_t j = base_dt->date.year;
          /* Initialize to 0. */
          td_return->year = 0;

          /* Fast-Fwd >= 400 */
          if (days >= NO_OF_DAYS_IN_400_YEARS)
            {
              int64_t numberOf400YearPeriods = days / NO_OF_DAYS_IN_400_YEARS;
              td_return->year = td_return->year + numberOf400YearPeriods * 400;
              j = j - numberOf400YearPeriods * 400;
              days = days - numberOf400YearPeriods * NO_OF_DAYS_IN_400_YEARS;
            }

          do
            {

              /* Loop over and get to the final year by substracting 366/365 days depending
                 on leap/non-leap year. For each substraction, increment year by 1.
              */

              /* The crucial point is month of february. */
              delta_final_year = days;
              if (
                  ( (testYearIsLeapYear(j)) && (base_dt->date.month >= 3) )
                  ||
                  ( (testYearIsLeapYear(j-1)) && (base_dt->date.month < 3) )
                  )
                {
                  /* If year is leap year and base month is >= 3
                     OR
                     next year is a leap year and month is < 3
                     => delta of 1 year corresponds to 366 day julian delta.
                  */
                  days = days - NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  days = days - NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }
              td_return->year++;
              j--;
            } while (days >= 0);

          /* The loop ran one time too much. */

          if (days < 0)
            {
              td_return->year--;
              j++;
            }

          /* In final year, the crucial point is the month of february. */

          if (
              ((testYearIsLeapYear(j)) && (base_dt->date.month >= 3))
              ||
              ((testYearIsLeapYear(j-1)) && (base_dt->date.month < 3))
              )
            {
              /* If final year is a leap year and base month is >= 3
                 OR
                 year preceding final year is a leap year and month is < 3
                 => An addition of leap-year specific delta for each month.
              */
              msdinm = monthSpecificDeltaInMonthsLeapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Otherwise. */
              msdinm = monthSpecificDeltaInMonths365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          for (i = 1; i <= NO_OF_MONTHS_IN_A_YEAR; i++)
            {
              if (delta_final_year < msdinm[base_dt->date.month - 1][i])
                {
                  /* Month  */
                  td_return->month = i - 1;
                  /* Day. */
                  td_return->day = (int) delta_final_year - msdinm[base_dt->date.month - 1][i - 1];
                  break;
                }
            }

          /* Time. */
          td_return->hour       = ((-1)*(int) jd->ms) / NO_OF_MS_IN_A_HOUR;
          td_return->minute     = ((-1)*(int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
          td_return->second     = ((-1)*(int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
          td_return->ms         =  (-1)*(int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
        }
      else if (jd->sign == '+')
        {
          td_return->sign = '+';

          /* No days in the final year. */
          int64_t delta_final_year;
          int64_t days = jd->day;
          /* Set counter to base year and then loop back to get to the final year.
             For each loop back, increment year by 1.
          */
          int64_t j = base_dt->date.year;
          /* Initialize. */
          td_return->year = 0;

          /* Fast-Fwd >= 400 */
          if (days >= NO_OF_DAYS_IN_400_YEARS)
            {
              int64_t numberOf400YearPeriods = days / NO_OF_DAYS_IN_400_YEARS;
              td_return->year = td_return->year + numberOf400YearPeriods * 400;
              j = j - numberOf400YearPeriods * 400;
              days = days - numberOf400YearPeriods * NO_OF_DAYS_IN_400_YEARS;
            }

          do
            {
              /* Loop over and get the year by substracting 366/365 days depending
                 on leap/non-leap year. For each substraction, increment year by 1.
              */

              /* The crucial point is month of february. */
              delta_final_year = days;
              if (
                  ((testYearIsLeapYear(j+1)) && (base_dt->date.month >= 3))
                  ||
                  ((testYearIsLeapYear(j)) && (base_dt->date.month < 3))
                  )
                {
                  /* If next year is leap year and base month is < 3
                     OR
                     this year is a leap year and month is < 3
                     => delta of 1 year corresponds to 366 day julian delta.
                  */
                  days = days - NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  days = days - NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }

              td_return->year++;
              j++;
            } while (days >= 0);

          /* The loop ran one time too much. */
          if (days < 0)
            {
              td_return->year--;
              j--;
            }

          /* In final year, the crucial point is the month of february. */
          if (
              ((testYearIsLeapYear(j+1)) && (base_dt->date.month >= 3))
              ||
              ((testYearIsLeapYear(j)) && (base_dt->date.month < 3))
              )
            {
              /* If final year is a leap year and base month is >= 3
                 OR
                 final year's previous year is a leap year and month is < 3
                 => An addition of leap-year specific delta for each month.
              */
              msdinm = monthSpecificDeltaInMonthsLeapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Otherwise. */
              msdinm = monthSpecificDeltaInMonths365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          td_return->month = 0;
          td_return->day = (int) delta_final_year;
          for (i = NO_OF_MONTHS_IN_A_YEAR; i > 0; i--)
            {
              if (delta_final_year >= msdinm[base_dt->date.month - 1][i])
                {
                  td_return->month = i;
                  td_return->day = (int) delta_final_year - msdinm[base_dt->date.month - 1][i];
                  break;
                }
            }

          /* Time */
          td_return->hour       =  (int) jd->ms / NO_OF_MS_IN_A_HOUR;
          td_return->minute     = ((int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
          td_return->second     = ((int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
          td_return->ms         =  (int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
        }
        td_return->flag_std_form = 1;

      return td_return;

      default:
      return NULL; /* Calender type not defined. */
    }

  /* Handle 360 and 365 calendar type here as they are similar. */

  /* 360 and 365 day calendars */
  if (jd->sign == '+')
    {
      /* Positive delta. */
      td_return->sign = '+';

      /* Year. */
      td_return->year = jd->day / ndiny;

      int delta_final_year = jd->day % ndiny;

      for (i = NO_OF_MONTHS_IN_A_YEAR; i > 0; i--)
        {
          if ((int) delta_final_year < (ndiny - msdinm[base_dt->date.month - 1][i - 1]))
            {
              /* Month. */
              td_return->month = NO_OF_MONTHS_IN_A_YEAR - i;
              /* Day. */
              td_return->day = (int) delta_final_year - (ndiny - msdinm[base_dt->date.month - 1][i]);
              break;
            }
        }

      /* Time. */
      td_return->hour   =  (int) jd->ms / NO_OF_MS_IN_A_HOUR;
      td_return->minute = ((int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
      td_return->second = ((int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
      td_return->ms     =  (int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
    }
  else if (jd->sign == '-')
    {
      /* Negative TimeDelta.  A negative juliandelta is represented in the following 
         way: -P01DT00.500S  = jd2->sign = '-', jd2->day = -30, jd2->ms = -500. */

      /* Negative delta. */
      td_return->sign = '-';

      /* Year. */
      td_return->year = ((-1)*jd->day) / ndiny;

      int delta_final_year = ((-1)*jd->day) % ndiny;

      for (i = 1; i <= NO_OF_MONTHS_IN_A_YEAR; i++)
        {
          if ((int)delta_final_year < msdinm[base_dt->date.month - 1][i])
            {
              /* Month. */
              td_return->month = i - 1;
              /* Day.  */
              td_return->day = (int) delta_final_year - msdinm[base_dt->date.month - 1][i - 1];
              break;
            }
        }

      /* Time. */
      td_return->hour   =  (-1)*(int) jd->ms / NO_OF_MS_IN_A_HOUR;
      td_return->minute = ((-1)*(int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
      td_return->second = ((-1)*(int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
      td_return->ms     =  (-1)*(int) jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
    } 
  else
  return NULL; /* ERROR: Sign of julian delta not defined. */

  return td_return;
}
else
  return NULL;
}
/*! \endcond */


struct _divisionquotienttimespan*
divideTimeDeltaInSeconds(struct _timedelta* dividend, struct _timedelta* divisor, struct _divisionquotienttimespan* quo_ret)
{
  if ((dividend != NULL) && (divisor != NULL) && (quo_ret != NULL))
    {
      if ((dividend->year == 0) && (dividend->month == 0) &&  (divisor->year == 0) && (divisor->month == 0))
        {
          intmax_t numerator   = (intmax_t) (((int64_t) dividend->day * 86400 + 
                                                        dividend->hour * 3600 + 
                                                        dividend->minute * 60 + 
                                                        dividend->second) * 1000 + 
                                                        dividend->ms);
          intmax_t denominator = (intmax_t) (((int64_t) divisor->day * 86400 + 
                                                        divisor->hour * 3600 + 
                                                        divisor->minute * 60 + 
                                                        divisor->second ) * 1000 + 
                                                        divisor->ms);

          if (denominator == 0) /* Division by zero is illegal. */
            return NULL;

          imaxdiv_t div = imaxdiv(numerator, denominator);
           
          quo_ret->quotient        = (int64_t) div.quot;          
          quo_ret->remainder_in_ms = (int64_t) div.rem;
          return quo_ret;
        }
    }
  return NULL;
}

/**
 * @brief division by an interval given in of seconds.
 *
 * @param  refdt
 *         Reference date (a pointer to struct _datetime).
 * @param  dt
 *         A pointer to struct _datetime.
 * @param  intvlsec
 *         Interval given in seconds.
 *
 * @return result of division. NULL indicates error.
 */
struct _divisionquotienttimespan*
divideDatetimeDifferenceInSeconds(struct _datetime* dt1, struct _datetime* dt2, struct _timedelta* divisor, struct _divisionquotienttimespan* quo_ret)
{
  if ((dt1 != NULL) && (dt2 != NULL) && (divisor != NULL) && (quo_ret != NULL))
    {
      if ((divisor->year == 0) && (divisor->month == 0))
        {
          struct _julianday* jd1 = newJulianDay(0, 0);
          if (jd1 != NULL) jd1 = date2julian(dt1, jd1);

          struct _julianday* jd2 = newJulianDay(0, 0);
          if (jd2 != NULL) jd2 = date2julian(dt2, jd2);

          intmax_t denominator = (intmax_t) (((int64_t) divisor->day  * 86400 + divisor->hour * 3600 +divisor->minute  * 60 + divisor->second ) * 1000 + divisor->ms);

          if (denominator == 0) /* Division by zero is illegal. */
            return NULL;

          struct _juliandelta* jd = newJulianDelta('+', 0, 0);
          jd = substractJulianDay(jd1, jd2, jd);  
          intmax_t numerator = (intmax_t) (jd->day * 86400000 + jd->ms);
          
          imaxdiv_t div = imaxdiv(numerator, denominator);
          
          quo_ret->quotient        = (int64_t) div.quot;          
          quo_ret->remainder_in_ms = (int64_t) div.rem;
          
          if (jd1 != NULL)  deallocateJulianDay(jd1);
          if (jd2 != NULL)  deallocateJulianDay(jd2);
          if (jd  != NULL)  deallocateJulianDelta(jd);
          return quo_ret;
        }
    }
  return NULL;
}

/**
 * @brief division of two differences in datetimes.
 *
 * @param  dt1_dividend, dt2_dividend, dt1_divisor, dt2_divisor
 *         Reference date (a pointer to struct _datetime).
 *
 * @param  intvlsec
 *         Interval given in seconds.
 *
 * @return result of division. NULL indicates error.
 */
struct _divisionquotienttimespan*
divideTwoDatetimeDiffsInSeconds(struct _datetime* dt1_dividend, struct _datetime* dt2_dividend, struct _datetime* dt1_divisor, struct _datetime* dt2_divisor, int64_t* denominator_ret,  struct _divisionquotienttimespan* quo_ret)
{
  if ((dt1_dividend != NULL) && (dt2_dividend != NULL) && 
      (dt1_divisor != NULL) && (dt2_divisor != NULL) && (quo_ret != NULL))
    {
          // dividend
          struct _julianday* jd1_dividend = newJulianDay(0, 0);
          if (jd1_dividend != NULL) jd1_dividend = date2julian(dt1_dividend, jd1_dividend);

          struct _julianday* jd2_dividend = newJulianDay(0, 0);
          if (jd2_dividend != NULL) jd2_dividend = date2julian(dt2_dividend, jd2_dividend);

          struct _juliandelta* jd_dividend = newJulianDelta('+', 0, 0);
          jd_dividend = substractJulianDay(jd1_dividend, jd2_dividend, jd_dividend);  
          intmax_t numerator = (intmax_t) (jd_dividend->day * 86400000 + jd_dividend->ms);

          // divisor
          struct _julianday* jd1_divisor = newJulianDay(0, 0);
          if (jd1_divisor != NULL) jd1_divisor = date2julian(dt1_divisor, jd1_divisor);

          struct _julianday* jd2_divisor = newJulianDay(0, 0);
          if (jd2_divisor != NULL) jd2_divisor = date2julian(dt2_divisor, jd2_divisor);

          struct _juliandelta* jd_divisor = newJulianDelta('+', 0, 0);
          jd_divisor = substractJulianDay(jd1_divisor, jd2_divisor, jd_divisor);  
          intmax_t denominator = (intmax_t) (jd_divisor->day * 86400000 + jd_divisor->ms);
          
          imaxdiv_t div = imaxdiv(numerator, denominator);
          
          quo_ret->quotient        = (int64_t) div.quot;          
          quo_ret->remainder_in_ms = (int64_t) div.rem;
          *denominator_ret = (int64_t) (denominator / 1000);
          
          if (jd1_dividend != NULL)  deallocateJulianDay(jd1_dividend);
          if (jd2_dividend != NULL)  deallocateJulianDay(jd2_dividend);
          if (jd_dividend  != NULL)  deallocateJulianDelta(jd_dividend);
          
          if (jd1_divisor != NULL)  deallocateJulianDay(jd1_divisor);
          if (jd2_divisor != NULL)  deallocateJulianDay(jd2_divisor);
          if (jd_divisor  != NULL)  deallocateJulianDelta(jd_divisor);
          
          return quo_ret;          
    }
  else 
    return NULL;
}

/**
 * @brief Get the TimeDelta between two Dates d1 and d2 as (d1-d2).
 *
 * Routine getTimeDeltaFromDate 'substracts' two Dates and returns the TimeDelta between
 * them. Internally, Dates are converted to DateTimes and then delta is calculated using 
 * getTimeDeltaFromDateTime().
 * 
 * This routine  handles all supported Calendar types; i.e. the translation from Calendar date 
 * to Julian date and conversion from Julian Delta to normal TimeDetla is Calendar-type dependent. 
 * For eg. for Calendar type Gregorian, the TimeDelta between 2001-02-01 and 2001-01-01 will be 1 month. 
 * Similarly, for Calendar of type 360-Day-Calendar, the TimeDelta will be 1 month. It must be noted 
 * however, that the two dates differ by 31 and 30 days respectively.
 *
 * @param  d1
 *         A pointer to struct _date. 
 *
 * @param  d2
 *         A pointer to struct _date.
 *
 * @param  td_return
 *         A pointer to struct _timedelta. Copy the result of (d1 - d2) in td_return.
 *
 * @return td_return
 *         A pointer to TimeDelta containing the result of substraction.
 */

struct _timedelta*
getTimeDeltaFromDate(struct _date* d1, struct _date* d2, struct _timedelta* td_return)
{
  if ((d1 != NULL )&& (d2 != NULL) && (td_return != NULL) ){

  /* Convert Date to datetime and resuse the DateTime interface to calculate time delta. */
  struct _datetime *dt1 = newDateTime(initDummyDTString);
  if (dt1 == NULL)
    return NULL;
  dt1 = convertDateToDateTime(d1, dt1);

  struct _datetime *dt2 = newDateTime(initDummyDTString);
  if ( dt2 == NULL )
    {
      deallocateDateTime(dt1);
      return NULL;
    }
  dt2 = convertDateToDateTime(d2, dt2);

  /* Call the Datetime function to get TD. dt1 - dt2 */
  td_return = getTimeDeltaFromDateTime(dt1,dt2,td_return);

  /* Cleanup. */
  deallocateDateTime(dt1);
  deallocateDateTime(dt2);

  return td_return;
}
else
  return NULL;
}


/**
 * @brief Get the TimeDelta between two DateTimes dt1 and dt2 as (dt1-dt2).
 *
 * Routine getTimeDeltaFromDateTime 'substracts' two DateTime's and returns the TimeDelta between
 * them. Each datetime is converted to an equivalent Julian Date. Substraction is then performed
 * on Julian axis. The "Julian delta" is finally converted back to normal calendar delta. 
 * 
 * This routine handles all supported Calendar types; i.e. the translation from Calendar date 
 * to Julian date and conversion from Julian Delta to normal TimeDetla is Calendar-type dependent. 
 * For eg. for Calendar type Gregorian, the TimeDelta between 2001-02-01T00:00:00.000 and 
 * 2001-01-01T00:00:00.000 will be 1 month. Similarly, for Calendar of type 360-Day-Calendar, 
 * the TimeDelta will be 1 month. It must be noted however, that the two dates differ by 31 and 
 * 30 days respectively. 
 *
 * @param  dt1
 *         A pointer to struct _datetime. 
 *
 * @param  dt2
 *         A pointer to struct _datetime.
 *
 * @param  td_return
 *         A pointer to struct _timedelta. Copy the result of (dt1 - dt2) in td_return.
 *
 * @return td_return
 *        A pointer to TimeDelta containing the result of substraction.
 */

struct _timedelta*
getTimeDeltaFromDateTime(struct _datetime* dt1, struct _datetime* dt2, struct _timedelta* td_return)
{
  if ((dt1 != NULL )&& (dt2 != NULL) && (td_return != NULL) ){
 
  /* Convert dt1 to Julian. */
  struct _julianday* jd1 = newJulianDay(0, 0);
  if (jd1 == NULL)
    return NULL;
  jd1 = date2julian(dt1, jd1);

  /* Convert dt2 to Julian. */
  struct _julianday* jd2 = newJulianDay(0, 0);
  if (jd2 == NULL)
    {
      deallocateJulianDay(jd1);
      return NULL;
    }
  jd2 = date2julian(dt2, jd2);

  /* Calculat Delta on Julian axis. "RULE: A - B = Delta". If A > B, Delta is positive. */
  struct _juliandelta* jd = newJulianDelta('+', 0, 0);
  if (jd == NULL)
    {
      deallocateJulianDay(jd1);
      deallocateJulianDay(jd2);
      return NULL;
    }

  /* Substract the 2 dates on julian axis. */
  jd = substractJulianDay(jd1, jd2, jd);

  /* Convert Julian-delta to TimeDelta. */
  td_return = julianDeltaToTimeDelta(jd, dt2, td_return);

  /* Cleanup. */
  deallocateJulianDay(jd1);
  deallocateJulianDay(jd2);
  deallocateJulianDelta(jd);

  return td_return;
}
else
  return NULL;
}


/**
 * @brief Get total number of milliseconds in timedelta.
 *
 * Routine getTotalMilliSecondsTimeDelta returns the total number of milliseconds in TimeDelta. 
 * Notice that TimeDelta is not uniquely defined but depends on the definition of corresponding 
 * DateTime. TimeDelta is first converted to corresponding delta on the Julian axis. Julian delta 
 * is finally converted to the correct millisecond value.
 *
 * @param  td
 *         A pointer to struct _timedelta. Retrieve the number of milliseconds in this TD object.
 *
 * @param  base_dt
 *         A pointer to struct _datetime. Reference Datetime for the TD.
 *
 * @return totalmilliSeconds
 *         Integer value of totalmilliSeconds. 0 indicates error. TODO on Luis: Is this ok?
 */

int64_t
getTotalMilliSecondsTimeDelta(struct _timedelta* td, struct _datetime* base_dt)
{
  if ((td != NULL )&& (base_dt != NULL) ){

  int64_t totalmilliSeconds = 0;

  struct _juliandelta* jd = newJulianDelta('+', 0, 0);
  if ( jd == NULL )
    return 0;
  jd = timeDeltaToJulianDelta(td, base_dt, jd);
  if ( jd == NULL )
    {
      deallocateJulianDelta(jd);
      return 0;
    }

  totalmilliSeconds = jd->day * NO_OF_MS_IN_A_DAY + jd->ms;

  deallocateJulianDelta(jd);

  return totalmilliSeconds;
}
else
  return 0;
}


/**
 * @brief Get total number of seconds in timedelta.
 *
 * Routine getTotalSecondsTimeDelta returns the total number of seconds in TimeDelta. Notice that TimeDelta 
 * is not uniquely defined but depends on the definition of corresponding DateTime. Internally, number of seconds
 * is calculated by calling the routine getTotalMilliSecondsTimeDelta() and then converting the millisecond value 
 * to seconds by dividing it by 1000.
 *
 * @param  td
 *         A pointer to struct _timedelta. Retrieve the number of seconds in this TD object.
 *
 * @param  base_dt
 *         A pointer to struct _datetime. Reference Datetime for the TD.
 *
 * @return totalSeconds
 *         Integer value of totalSeconds. 0 indicates error.
 */

int64_t
getTotalSecondsTimeDelta(struct _timedelta* td, struct _datetime* base_dt)
{
  if ((td != NULL ) && (base_dt != NULL))
    return getTotalMilliSecondsTimeDelta(td, base_dt) / NO_OF_MS_IN_A_SECOND;
  else
    return 0;
}


/**
 * @brief Get TimeDelta as an extended string.
 *
 * timedeltaToString returns a string in IS08601 compliant (and extended) format.
 *
 * @param  td
 *         A pointer to struct _timedelta. The timedelta to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where timedelta is to be written.
 *
 * @return toStr
 *         A pointer to the string containing timedelta.
 */

char *
timedeltaToString(struct _timedelta* td, char* toStr)
{
  if ((td != NULL )&& (toStr != NULL) ){

  memset(toStr,'\0',MAX_TIMEDELTA_STR_LEN);

  if (td->sign == '-')
    sprintf(toStr,"%cP",td->sign);
  else if (td->sign == '+')
    strcpy (toStr,"P");
  else
    return NULL; /* ERROR: TD sign not set. Should never happen. */
  
  if (td->year != 0)
    sprintf(&(toStr[strlen(toStr)]),"%" PRIi64 "Y",td->year);

  if (td->month != 0)
    sprintf(&(toStr[strlen(toStr)]),"%02dM",td->month);

  if (td->day != 0)
    sprintf(&(toStr[strlen(toStr)]),"%02dD",td->day);

  sprintf(&(toStr[strlen(toStr)]),"T");

  if (td->hour != 0)
    sprintf(&(toStr[strlen(toStr)]),"%02dH",td->hour);

  if (td->minute != 0)
    sprintf(&(toStr[strlen(toStr)]),"%02dM",td->minute);

  if ((td->second != 0) || (td->ms != 0))
    {
      if (td->second > 59)
        {
          /* reset and generate string with seconds only */
          memset(toStr,'\0',MAX_TIMEDELTA_STR_LEN);

          if (td->sign == '-')
            sprintf(toStr,"%cPT",td->sign);
          else if (td->sign == '+')
            strcpy (toStr,"PT");
          else
            return NULL; /*ERROR: TD sign not set. Should never happen. */

          sprintf(&(toStr[strlen(toStr)]),"%d.%03dS",td->second, td->ms);
        }
      else
        {
          sprintf(&(toStr[strlen(toStr)]),"%02d.%03dS",td->second,td->ms);
        }
    }


  //Discard T if all time values are 0.
  if(toStr[strlen(toStr)-1] == 'T')
    toStr[strlen(toStr)-1] = '\0';

  //Return PT00.000S if all delta values are 0.
  if(toStr[strlen(toStr)-1] == 'P')
    strcat(toStr,"T00.000S");

  return toStr;
}
else
  return NULL;
}


/**
* @brief Add timedelta to Date.
*
* Routine addTimeDeltaToDate adds a timedelta to a Date and returns the new Date. Both Date 
* and TimeDetla are first converted to corresponding values on the Julian axis. Addition is performed on
* the julian axis and the resulting Julian Date is converted back to the corrsponding Date. 
*
* The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
* positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
* Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
* a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
* so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
* anchor date month; a delta of 2 month corresponds to as many days as in the month before the
* anchor date month and the month before that and so on.
*
* @param  d
*         A pointer to struct _date. The base date.
*
* @param  td
*         A pointer to struct _timedelta. The time delta to be added to d.
*
* @param  d_return
*         A pointer to struct _date. The result of addition is copied here.
*
* @return d_return
*         A pointer to the struct _date contianing the result of addition.
*/

struct _date *
addTimeDeltaToDate(struct _date* d, struct _timedelta* td, struct _date* d_return)
{
  if ((d != NULL )&& (td != NULL) && (d_return != NULL) ){

  /* Convert Date to Datetime and reuse the DateTime interface for Calculating the sum.*/
  struct _datetime *dt = newDateTime(initDummyDTString);
  if ( dt == NULL )
    return NULL;

  dt = convertDateToDateTime(d, dt);

  struct _datetime *dt_return = newDateTime(initDummyDTString);
  if ( dt_return == NULL )
  {
    deallocateDateTime(dt);
    return NULL;
  }

  /* Call the DateTime interface to calculate the new Datetime. */
  dt_return = addTimeDeltaToDateTime(dt,td, dt_return);

  /* Get Date from Datetime. */
  d_return = convertDateTimeToDate(dt_return, d_return);

  deallocateDateTime(dt);
  deallocateDateTime(dt_return);

  return d_return;
  }
  else
    return NULL;
}


/**
 * @brief Add timedelta to DateTime.
 *
 * Routine addTimeDeltaToDateTime adds a timedelta to a DateTime and returns the new DateTime. Both DateTime 
 * and TimeDetla are first converted to corresponding values on the Julian axis. Addition is performed on
 * the julian axis and the resulting Julian Date is converted back to the corrsponding DateTime.
 *
 * The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
 * positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
 * Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
 * a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
 * so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
 * anchor date month; a delta of 2 month corresponds to as many days as in the month before the
 * anchor date month and the month before that and so on.
 *
 * @param  dt
 *         A pointer to struct _datetime. The base datetime.
 *
 * @param  td
 *         A pointer to struct _timedelta. The time delta to be added to dt.
 *
 * @param  dt_return
 *         A pointer to struct _datetime. The result of addition is copied here.
 *
 * @return dt_return
 *         A pointer to the struct _datetime contianing the result of addition.
 */

struct _datetime*
addTimeDeltaToDateTime(struct _datetime* dt, struct _timedelta* td, struct _datetime* dt_return)
{
  if ((dt != NULL ) && (td != NULL) && (dt_return != NULL)) {

    /* Convert base datetime to Julian. */
    struct _julianday* jd1 = newJulianDay(0, 0);
    if ( jd1 == NULL)
      return NULL;
    jd1 = date2julian(dt, jd1);

    /* Get julian delta. */
    struct _juliandelta* jd2 = newJulianDelta('+', 0, 0);
    if ( jd2 == NULL )
      {
        deallocateJulianDay(jd1);
        return NULL;
      }
    jd2 = timeDeltaToJulianDelta(td, dt, jd2);

    struct _julianday* jd = newJulianDay(0, 0);
    if ( jd == NULL )
      {
        deallocateJulianDay(jd1);
        deallocateJulianDelta(jd2);
        return NULL;
      }

    if ( td->sign == '+' )
      {
        jd = addJulianDeltaToJulianDay(jd1, jd2, jd);
      }
    else if ( td->sign == '-' )
      {
        jd = substractJulianDeltaFromJulianDay(jd1, jd2, jd);
      }
    else
      return NULL; /* ERROR: Sign of timedelta is not defined. Should never happen. */

    /* Get the Datetime */
    dt_return = julian2date(jd, dt_return);

    deallocateJulianDay(jd1);
    deallocateJulianDelta(jd2);
    deallocateJulianDay(jd);

    return dt_return;
  }
  else
    return NULL;
}


/**
 * @brief Get the timedelta between current_dt and start_dt plus next integral-multiple-of-timestep (timedelta).
 *
 * Routine moduloTimeDeltaFromDateTime returns the timedelta between the current DateTime (current_dt) and the event's next-trigger time.
 * The next trigger time is defined as the the Anchor DateTime (start_dt) + N * TimeDelta(timestep) 
 * where N is the minimum positive integer for which this sum is >= Current DateTime. In case 
 * Anchor DateTime > Current DateTime, TimeDelta is calculated as start_dt - current_dt.
 * 
 * Notice that this TimeDelta will always be positive.
 *
 * @param  start_dt
 *         A pointer to struct _datetime. The base datetime.
 *
 * @param  timestep
 *         A pointer to struct _timedelta. delta between two consecutive triggers.
 *
 * @param  current_dt
 *         A pointer to struct _datetime. The Current Date time.
 *
 * @param  modulo_td
 *         A pointer to struct _timedelta. The timedelta between 'current datetime' and 'Start Datetime plus next 
 *         integral-multiple-of-timestep' is copied here.
 *
 * @return modulo_td
 *         A pointer to the struct _timedelta contianing the modulo timedelta. If Start time is in the future, returns on start_time - current_time.
 */

struct _timedelta*
moduloTimeDeltaFromDateTime(struct _datetime* start_dt, struct _timedelta* timestep, struct _datetime* current_dt, struct _timedelta* modulo_td)
{
  if ((start_dt != NULL )&& (timestep != NULL) && (current_dt != NULL) && (modulo_td != NULL) ){

  struct _datetime* dt_tmp = newDateTime(initDummyDTString);
  if ( dt_tmp == NULL )
    return NULL;

  if (compareDatetime(start_dt,current_dt)==(less_than))
    {
      /* Loop over */
      replaceDatetime(start_dt,dt_tmp);
      while(compareDatetime((dt_tmp = addTimeDeltaToDateTime(dt_tmp,timestep,dt_tmp)),current_dt) == less_than);

      /* Return n*dt_tmp - current_dt */
      modulo_td = getTimeDeltaFromDateTime(dt_tmp,current_dt,modulo_td);
    }
  else
    {
      /* Start time is in the future, return start_time - current_time. */
      modulo_td = getTimeDeltaFromDateTime(start_dt,current_dt,modulo_td);
    }

  deallocateDateTime(dt_tmp);

  return modulo_td;
  }
  else
    return NULL;
}


/**
 * @brief Return the element-wise product of a scalar and a timedelta.
 *
 * elementwiseScalarMultiplyTimeDelta multiplies scalar lambda with each element of timedelta and returns the result in scaled_td.
 * Scalar can be both positive and negative and so can the timedelta. The timedelta can not have days,months or years however: Only
 * Timedeltas upto hours should call this routine. Also scaled_td->hour >= 24 will lead to an error.
 *
 *
 * @param  base_td
 *         A pointer to struct _timedelta. The base timedelta.
 *
 * @param  lambda
 *         A scalar to be multiplied.
 *
 * @param  scaled_td
 *         A pointer to struct _timedelta. The element-wise product of lambda and base_td is stored here.
 *
 * @return scaled_td
 *        A pointer to struct _timedelta. The filled structure containing the scaled timedelta values.    
 */
struct _timedelta*
elementwiseScalarMultiplyTimeDelta(struct _timedelta* base_td, int64_t lambda, struct _timedelta* scaled_td)
{
  if ((base_td != NULL) && (scaled_td != NULL) )
    {
      /* Scalar multiplication not supported for TimeDeltas consisting of day/month/year. */
      if ( base_td->day > 0 || base_td->month > 0 || base_td->year > 0 )
        return NULL;

      /* Reset scaled_td to 0. */
      memset(scaled_td,0,sizeof(struct _timedelta));

      /* Scalar can be positive or negative. */
      if (((lambda < 0) && (base_td->sign == '+')) || ((lambda > 0) && (base_td->sign == '-')))
        scaled_td->sign = '-';
      else
        scaled_td->sign = '+';

      /* Sign already handled above. Make lambda positive. */
      if (lambda < 0)
        lambda *= -1;


      /* Multiply each element by scalar. */

      scaled_td->ms += (int) lambda*base_td->ms;
      while ( scaled_td->ms >= NO_OF_MS_IN_A_SECOND )
        {
          scaled_td->ms -= NO_OF_MS_IN_A_SECOND;
          scaled_td->second += 1;
        }
 
      scaled_td->second += (int) lambda*base_td->second;
      while ( scaled_td->second >= 60 )
        {       
          scaled_td->second -= 60;
          scaled_td->minute += 1;
        }
   
      scaled_td->minute += (int) lambda*base_td->minute;
      while ( scaled_td->minute >= 60 )
        {
          scaled_td->minute -= 60;
          scaled_td->hour += 1;
        }

      scaled_td->hour += (int) lambda*base_td->hour;

     /* Scalar multiplication can not give a value in excess of 24 hours. */
     if ( scaled_td->hour >= NO_OF_HOURS_IN_A_DAY )
       {
         /* ERROR: Return on NULL. */
         return NULL;
       }

      scaled_td->day    = 0;
      scaled_td->month  = 0;
      scaled_td->year   = 0;

      return scaled_td;
        
    }
  else
    return NULL;
}


/**
 * @brief Return the element-wise product of a REAL scalar and a timedelta.
 *
 * elementwiseScalarMultiplyTimeDeltaDP multiplies scalar lambda with
 * each element of timedelta and returns the result in scaled_td.
 * Scalar can be both positive and negative and so can the
 * timedelta. The timedelta can not have days,months or years however:
 * Only Timedeltas upto hours should call this routine. Also
 * scaled_td->hour >= 24 will lead to an error.
 *
 *
 * @param  base_td
 *         A pointer to struct _timedelta. The base timedelta.
 *
 * @param  lambda
 *         A scalar to be multiplied (double precision floating-point value).
 *
 * @param  scaled_td
 *         A pointer to struct _timedelta. The element-wise product of lambda and base_td is stored here.
 *
 * @return scaled_td
 *        A pointer to struct _timedelta. The filled structure containing the scaled timedelta values.    
 */
struct _timedelta*
elementwiseScalarMultiplyTimeDeltaDP(struct _timedelta* base_td, double lambda, struct _timedelta* scaled_td)
{
  if ((base_td != NULL) && (scaled_td != NULL) )
    {
      /* Scalar multiplication not supported for TimeDeltas consisting of day/month/year. */
      if ( base_td->day > 0 || base_td->month > 0 || base_td->year > 0 )
        return NULL;
      /* Reset scaled_td to 0. */
      memset(scaled_td,0,sizeof(struct _timedelta));
      /* Scalar can be positive or negative. */
      if (((lambda < 0.) && (base_td->sign == '+')) || ((lambda > 0.) && (base_td->sign == '-')))
        scaled_td->sign = '-';
      else
        scaled_td->sign = '+';
      /* Sign already handled above. Make lambda positive. */
      if (lambda < 0.)
        lambda *= -1.;

      /* Multiply each element by scalar. */

      scaled_td->hour = (int) lambda*base_td->hour;
      /* Scalar multiplication can not give a value in excess of 24 hours. */
      if ( scaled_td->hour >= NO_OF_HOURS_IN_A_DAY )
        {
          /* ERROR: Return on NULL. */
          return NULL;
        }
      double remainder_minutes = 60.*(lambda*base_td->hour - scaled_td->hour);
      scaled_td->minute = (int) (lambda*base_td->minute + remainder_minutes);
      double remainder_seconds = 60.*(lambda*base_td->minute + remainder_minutes - scaled_td->minute);
      scaled_td->second += (int) (lambda*base_td->second + remainder_seconds);
      double remainder_ms      = 1000.*((lambda*base_td->second + remainder_seconds) - scaled_td->second);
      scaled_td->ms = (int) (lambda*base_td->ms + remainder_ms);
      scaled_td->day    = 0;
      scaled_td->month  = 0;
      scaled_td->year   = 0;
      return scaled_td;
    }
  else
    return NULL;
}


/**
 * @brief Return the element-wise sum of two timedeltas.
 *
 * elementwiseAddTimeDeltatoTimeDelta adds two timedeltas elementwise and returns the result.
 * Timedeltas being added must be of the same sign; Substraction is not supported. 
 * The timedelta can not have days,months or years however: Only Timedeltas upto hours should call this routine. 
 * Also td_return->hour >= 24 will lead to an error.
 *
 *
 * @param  td1
 *         A pointer to struct _timedelta.
 *
 * @param  td2
 *         A pointer to struct _timedelta.
 *
 * @param  td_return
 *         A pointer to struct _timedelta. The element-wise sum of td1 and td2 is stored here.
 *
 * @return td_return
 *         A pointer to struct _timedelta. The filled structure containing the added timedelta values.   
 */
struct _timedelta*
elementwiseAddTimeDeltatoTimeDelta(struct _timedelta* td1, struct _timedelta* td2,  struct _timedelta* td_return)
{
  if ( (td1 != NULL) && (td2 != NULL) && (td_return != NULL) )
    {
      /* TD addition not supported for TimeDeltas consisting of day/month/year. */
      if ( td1->day > 0 || td1->month > 0 || td1->year > 0 || td2->day > 0 || td2->month > 0 || td2->year > 0)
        return NULL;

      /*Reset td_return to 0.*/
      memset(td_return,0,sizeof(struct _timedelta));

      if(td1->sign == td2->sign)
        {
          /* If signs match, do add. */
          td_return->sign = td1->sign;

          td_return->ms += (td1->ms + td2->ms);  
          if ( td_return->ms >= NO_OF_MS_IN_A_SECOND )
            {
              td_return->ms -= NO_OF_MS_IN_A_SECOND;
              td_return->second += 1;
            }

          td_return->second += (td1->second + td2->second);
          if ( td_return->second >= 60 )
            {
	      int remainder, amount;
	      remainder = td_return->second % 60;
	      amount    = td_return->second / 60; 
              td_return->second = remainder;
              td_return->minute = amount;
            }

          td_return->minute += (td1->minute + td2->minute);
          if ( td_return->minute >= 60 )
            {
	      int remainder, amount;
	      remainder = td_return->minute % 60;
	      amount    = td_return->minute / 60; 
              td_return->minute = remainder;
              td_return->hour = amount;
            }

          td_return->hour += (td1->hour + td2->hour);
          if ( td_return->hour >= NO_OF_HOURS_IN_A_DAY )
            {
              /* Hour-sum can not be allowed to exceed 24 hours. */
              return NULL;
            }

          td_return->day        = 0;
          td_return->month      = 0;
          td_return->year       = 0;
        } 
      else
        {
          /* Substraction not supported. */
          return NULL; 
        }

      return td_return;
    } 
  else
    return NULL;

}

/**
 * @brief Returns the remainder of timedelta a modulo timedelta p  
 *
 * moduloTimedelta(struct _timedelta *a, struct _timdelta *p) returns the remainder of
 * a modulo p. The function is restricted to input timedeltas less than 29 days.
 *
 * @param  a
 *         a valid timedelta objact <  P29D.
 *
 * @param  p
 *         a valid timedelta objact <  P29D.
 *
 * @param  p
 *         on return contains the quotient
 *
 * @return rem 
 *         returns remainder of the division of the numerator a by
 *         the denominator p. quotient * p + rem shall equal a.
 */

int64_t
moduloTimedelta(struct _timedelta *a, struct _timedelta *p, int64_t *quot)
{
  struct _datetime *dt = newDateTime("0001-01-01");
    
  struct _juliandelta *jdd1 = newJulianDelta('+', (int64_t) 0, (int64_t) 0);
  struct _juliandelta *jdd2 = newJulianDelta('+', (int64_t) 0, (int64_t) 0);
        
  jdd1 = timeDeltaToJulianDelta(a, dt, jdd1);
  jdd2 = timeDeltaToJulianDelta(p, dt, jdd2);

  int64_t d1, d2;
  
  d1 = (int64_t) 86400000 * jdd1->day + jdd1->ms;
  d2 = (int64_t) 86400000 * jdd2->day + jdd2->ms;

  ldiv_t d = ldiv(d1, d2);

  *quot = d.quot;

  return d.rem;
}

/**
 * @brief Return a PT String corresponding to arbitrary number of milliseconds.
 *
 * getPTStringFromMS() translates ms values to ISO 8601 compliant timedelta string. 
 * Conversion of ms >= 86400000 and  ms <= -86400000 not supported.
 *
 * @param  _ms
 *         An int64_t value to be translated.
 *
 * @param  PTstr
 *         A pointer to char. Translated string is written here.
 *
 * @return PTstr
 *         A pointer to char. The translated TimeDelta string.   
 */

char*
getPTStringFromMS(int64_t _ms, char* PTstr)
{
   /* Reuse the juliandelta to TimeDelta conversion routine. */
   
   /* Create a _juliandelta object and copy the _ms to the jd->ms. */
   struct _juliandelta* jd = NULL; 
   /* if ((_ms >= NO_OF_MS_IN_A_DAY) || (_ms <= ((-1)*NO_OF_MS_IN_A_DAY))) */
   /*   { */
   /*     /\* ERROR: Conversion greater than 23:59:59:999 not supported. *\/ */
   /*     return NULL; */
   /*   } */

   if (_ms >= 0)
     jd = newJulianDelta('+', 0, _ms); 
   else
     jd = newJulianDelta('-', 0, _ms);

   /* Create dummy variables for julianDeltaToTimeDelta() */
   struct _datetime* dumm_base_dt       = newDateTime(initDummyDTString);
   struct _timedelta* dummy_td_return   = newTimeDelta(initDummyTDString);
   
   /* Get the translated TimeDelta and return the corresponding string. */
   PTstr = timedeltaToString(julianDeltaToTimeDelta(jd, dumm_base_dt, dummy_td_return), PTstr);
     
   deallocateDateTime(dumm_base_dt);
   deallocateTimeDelta(dummy_td_return);
   if(jd)
     deallocateJulianDelta(jd);

   return PTstr;
}


/**
 * @brief Return a PT String corresponding to arbitrary number of seconds.
 *
 * getPTStringFromSeconds() translates second values to ISO 8601 compliant timedelta string. 
 * Conversion of s >= 86400 and  s <= -86400 not supported.
 *
 * @param  _s
 *         An int64_t value to be translated.
 *
 * @param  PTstr
 *         A pointer to char. Translated string is written here.
 *
 * @return PTstr
 *         A pointer to char. The translated TimeDelta string.   
 */

char*
getPTStringFromSeconds(int64_t _s, char* PTstr)
{
  return getPTStringFromMS( _s * NO_OF_MS_IN_A_SECOND, PTstr);
}


/**
 * @brief Return a PT String corresponding to arbitrary number of seconds.
 *
 * getPTStringFromSecondsFloat() translates second values to ISO 8601 compliant timedelta string. 
 * Conversion of s >= 86400 and  s <= -86400 not supported. Returned PT string has a precision of 3
 * digits.
 *
 * @param  _s
 *         Float value to be translated.
 *
 * @param  PTstr
 *         A pointer to char. Translated string is written here.
 *
 * @return PTstr
 *         A pointer to char. The translated TimeDelta string.   
 */

char*
getPTStringFromSecondsFloat(float _s, char* PTstr)
{
  return getPTStringFromMS( roundf(_s * NO_OF_MS_IN_A_SECOND), PTstr);
}


/**
 * @brief Return a PT String corresponding to arbitrary number of seconds.
 *
 * getPTStringFromSecondsDouble() translates second values to ISO 8601 compliant timedelta string. 
 * Conversion of s >= 86400 and  s <= -86400 not supported. Returned PT string has a precision of 3 digits.
 *
 * @param  _s
 *         Double value to be translated.
 *
 * @param  PTstr
 *         A pointer to char. Translated string is written here.
 *
 * @return PTstr
 *         A pointer to char. The translated TimeDelta string.   
 */

char*
getPTStringFromSecondsDouble(double _s, char* PTstr)
{
  return getPTStringFromMS( round(_s * NO_OF_MS_IN_A_SECOND), PTstr);
}


/**
 * @brief Return a PT String corresponding to arbitrary number of minutes.
 *
 * getPTStringFromMinutes() translates minutes values to ISO 8601 compliant timedelta string. 
 * Conversion of m >= 1440 and  m <= -1440 not supported.
 *
 * @param  _m
 *         An int64_t value to be translated.
 *
 * @param  PTstr
 *         A pointer to char. Translated string is written here.
 *
 * @return PTstr
 *         A pointer to char. The translated TimeDelta string.   
 */

char*
getPTStringFromMinutes(int64_t _m, char* PTstr)
{
  return getPTStringFromMS( _m * NO_OF_MS_IN_A_MINUTE, PTstr);
}


/**
 * @brief Return a PT String corresponding to arbitrary number of Hours.
 *
 * getPTStringFromHours() translates hour values to ISO 8601 compliant timedelta string. 
 * Conversion of h >= 24 and  ms <= -24 not supported.
 *
 * @param  _h
 *         An int64_t value to be translated.
 *
 * @param  PTstr
 *         A pointer to char. Translated string is written here.
 *
 * @return PTstr
 *         A pointer to char. The translated TimeDelta string.   
 */

char*
getPTStringFromHours(int64_t _h, char* PTstr)
{
  return getPTStringFromMS( _h * NO_OF_MS_IN_A_HOUR, PTstr);
}
