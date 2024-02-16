/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_timedelta.h
 * 
 * @brief TimeDelta and some operations supported on TimeDelta.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology. 
 * 
 * @date March 2013
 * 
 */

#ifndef _MTIME_TIMEDELTA_H
#define _MTIME_TIMEDELTA_H

#include <stdint.h>
#include <stdbool.h>

#include "mtime_calendar.h"


struct _datetime;
struct _date;
struct _julianday;
struct _juliandelta;


/**
 * @struct _timedelta
 *
 * @brief Struct _timedelta containing timedelta and sign of year parameter. 
 */

struct _timedelta
{
  int          flag_std_form; ///< Is timedelta formatted in standard form (eg. PT30M10S) or long form (eg. PT3000S)? 
  char 		sign;	///< sign of time delta. Sign can be '+' or '-'.

  int64_t 	year;	///< Year part of timedelta.
  int 		month;	///< Month part of timedelta.
  int 		day;	///< Day part of timedelta.
  int 		hour;	///< Hour part of timedelta.
  int 		minute;	///< Minute part of timedelta.
  int 		second;	///< Second part of timedelta.
  int 		ms;	///< Milli-Second part of timedelta.
};


/**
 * @struct _divisionquotienttimespan
 *
 * @brief Struct _divisionquotienttimespan is used for storing division of two time-delta results. 
 */

struct _divisionquotienttimespan
{
  int64_t quotient;          ///< Quotient of two timedeltas.
  int64_t remainder_in_ms;   ///< Remainder in milli seconds of two timedeltas division.
};

struct _timedelta*
newTimeDelta(const char* timedelta_string);

struct _timedelta*
newRawTimeDelta(char _sign, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms);

struct _timedelta*
constructAndCopyTimeDelta(struct _timedelta* td);

void
deallocateTimeDelta(struct _timedelta* td);

compare_return_val
compareTimeDelta(struct _timedelta* td1, struct _timedelta* td2);

struct _timedelta*
replaceTimeDelta(struct _timedelta* tdsrc, struct _timedelta* tddest);

struct _juliandelta*
timeDeltaToJulianDelta(struct _timedelta *td, struct _datetime *dt, struct _juliandelta *jd);

struct _timedelta*
julianDeltaToTimeDelta(struct _juliandelta* jd, struct _datetime *dt, struct _timedelta* td_return);

struct _divisionquotienttimespan*
divideTimeDeltaInSeconds(struct _timedelta* dividend, struct _timedelta* divisor, struct _divisionquotienttimespan* quo_ret);

struct _divisionquotienttimespan*
divideTwoDatetimeDiffsInSeconds(struct _datetime* dt1_dividend, struct _datetime* dt2_dividend,struct _datetime* dt1_divisor, struct _datetime* dt2_divisor, int64_t * denominator_ret,  struct _divisionquotienttimespan* quo_ret);

struct _divisionquotienttimespan*
divideDatetimeDifferenceInSeconds(struct _datetime* dt1, struct _datetime* dt2, struct _timedelta* divisor, struct _divisionquotienttimespan* quo_ret);

struct _timedelta*
getTimeDeltaFromDate(struct _date*, struct _date*, struct _timedelta*);

struct _timedelta*
getTimeDeltaFromDateTime(struct _datetime* dt1, struct _datetime* dt2, struct _timedelta* td_return);

int64_t
getTotalMilliSecondsTimeDelta(struct _timedelta *td, struct _datetime* dt);

int64_t
getTotalSecondsTimeDelta(struct _timedelta *td, struct _datetime* dt);

char*
timedeltaToString(struct _timedelta * td, char *toString);

struct _date*
addTimeDeltaToDate(struct _date* dt, struct _timedelta* td, struct _date* dt_return);

struct _datetime*
addTimeDeltaToDateTime(struct _datetime* dt, struct _timedelta* td, struct _datetime* dt_return);

struct _timedelta*
moduloTimeDeltaFromDateTime(struct _datetime* start_dt, struct _timedelta* timestep, struct _datetime* current_dt, struct _timedelta* modulo_td);

struct _timedelta*
elementwiseScalarMultiplyTimeDelta(struct _timedelta* base_td, int64_t lambda, struct _timedelta* scaled_td);

struct _timedelta*
elementwiseAddTimeDeltatoTimeDelta(struct _timedelta* td1, struct _timedelta* td2,  struct _timedelta* td_return);

int64_t
moduloTimedelta(struct _timedelta *a, struct _timedelta *p, int64_t *quot);

char*
getPTStringFromMS(int64_t _ms, char* PTstr);

char*
getPTStringFromSeconds(int64_t _s, char* PTstr);

char*
getPTStringFromSecondsFloat(float _s, char* PTstr);

char*
getPTStringFromSecondsDouble(double _s, char* PTstr);

char*
getPTStringFromMinutes(int64_t _m, char* PTstr);

char*
getPTStringFromHours(int64_t _h, char* PTstr);

/**
 * @}
 */
#endif
