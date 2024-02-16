# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
/*! \cond PRIVATE */
/**
 * @brief ISO 8601_2004 complaint Time.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note USAGE: Compile this rl file and generate iso8601.c file. 
	  Compile the iso8601.c file using a C compiler. iso8601.h file needs to be edited seperately. 
          match_found = 1 => DATE/DATETIME. match_found = 2 => Duration. Else non-compliant string and hence REJECT.
	  Due to application requirements, current implementation allows year in the range 2147483647 and -2147483648 only!
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>

#include "mtime_iso8601.h"

#define MAX_BUFFER_LENGTH 132

#define YEAR_UPPER_BOUND 2147483647L
#define YEAR_LOWER_BOUND -2147483648L

//#define SECOND_UPPER_BOUND 86399 

/* Allowed year range = 2147483647 TO -2147483648  */
bool RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = false;

/* Allowed second max = 86399 */
bool RAISE_SECOND_UPPER_LIMIT_EXCEPTION = false;

%%{
    machine date_machine;
    write data;
}%%


struct internal_datetime
  {
    char            sign_of_year;
    int64_t         year;
    int             month;
    int             day;
    int             hour;
    int             minute;
    int             second;
    int 	    ms;
  };


static 
void 
date_machine( char *str, ISO8601_STATUS* stat, struct internal_datetime* dtObj, struct iso8601_duration* duObj)
  {
    char *p = str, *pe = str + strlen( str );
    char *ts, *te = 0;
    int cs;

    RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = false;	
    RAISE_SECOND_UPPER_LIMIT_EXCEPTION = false;	

    %%{
	action _start 
	  {
	    ts = p;
	  }                
	
	action date_match
	  {
	    *stat = DATETIME_MATCH;
	  }

	action duration_match_iso_normal
	  {
	    *stat = DURATION_MATCH_STD;
	  }

         action duration_match_long_form
          {
            *stat = DURATION_MATCH_LONG;
          }

	action year_sign_match
	  {
	    dtObj->sign_of_year = fc;
	  }
			
	action duration_sign_match
	  {
	    duObj->sign = fc;
	  }

	action year_year_copy
	  {
	    te = p+1; 
	    /* Reset ts to point to begining of string. */
	    ts = str;                      

	    char _year[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _year, ts, (size_t)(te-ts));
	    _year[MAX_BUFFER_LENGTH-1] = '\0';
	    
	    /* To ensure strtol works. */
	    if ( _year[0] == '$')
	      _year[0] = '-';

            long _yearl;
	    char *end;
	    _yearl = strtol(_year, &end, 10);
	                
           if (end == _year)
              {
                // fprintf(stderr, "%s: not a decimal number\n", _year);                
              }
            /*
             * Ignore this case, as a - is trailing alwyas.
             *
             * else if ('\0' != *end)
             *   {
             *     fprintf(stderr, "%s: extra characters at end of input: %s\n", _year, end);
             *   }
             */
            else if ((LONG_MIN == _yearl || LONG_MAX == _yearl) && ERANGE == errno)
              {
                // fprintf(stderr, "%s out of range of type long\n", _year);
                RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = true;
              }
            else if ((_yearl > YEAR_UPPER_BOUND) || (_yearl < YEAR_LOWER_BOUND))
              {
                // fprintf(stderr, "%s out of range of user defined year range\n", _year);
                RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = true;
              }
            else
              {
                // fprintf(stderr, "Correct year %s \n", _year);
                dtObj->year = (int64_t) _yearl;
              }
	  }

	action year_month_copy               
	  {
	    te = p+1;
	    char _month[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _month, ts, (size_t)(te-ts));
	    _month[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->month = atoi(_month);
	  }

	action year_day_copy                 
	  {
	    te = p+1;
	    char _day[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _day, ts, (size_t)(te-ts));
	    _day[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->day = atoi(_day);
	  }

	action year_hour_copy                
	  {
	    te = p+1;
	    char _hour[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _hour, ts, (size_t)(te-ts));
	    _hour[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->hour = atoi(_hour);
	  }

	action year_minute_copy              
	  {
	    te = p+1;
	    char _minute[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _minute, ts, (size_t)(te-ts));
	    _minute[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->minute = atoi(_minute);
	  }

	action year_second_copy              
	  {
	    te = p+1;
	    char _second[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _second, ts, (size_t)(te-ts));
	    _second[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->second = atoi(_second);                
	  }

	action year_ms_copy          
	  {
	    te = p+1;
	    char _ms[8] = {'\0'};
	    strncpy( _ms, ts, (size_t)(te-ts));
	    _ms[8-1] = '\0';
            if(strlen(_ms) == 1)
              dtObj->ms = atoi(_ms)*100;
            else if(strlen(_ms) == 2)
              dtObj->ms = atoi(_ms)*10;
            else
              dtObj->ms = atoi(_ms);
          }


	action duration_year_copy            
	  {
	    te = p;                 
	    char _du_year[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_year, ts, (size_t)(te-ts));
	    _du_year[MAX_BUFFER_LENGTH-1] = '\0';

            long _yearl;
	    char *end;
            _yearl = strtol(_du_year,&end, 10);

            if (end == _du_year)
              {
                // fprintf(stderr, "%s: not a decimal number\n", _du_year);                
              }
            /*
             * Ignore this case, as a - is trailing alwyas.
             *
             * else if ('\0' != *end)
             *   {
             *     fprintf(stderr, "%s: extra characters at end of input: %s\n", _du_year, end);
             *   }
             */
            else if ((LONG_MIN == _yearl || LONG_MAX == _yearl) && ERANGE == errno)
              {
                // fprintf(stderr, "%s out of range of type long\n", _du_year);
                RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = true;
              }
            else if (_yearl > (YEAR_UPPER_BOUND + 1)) // abs(YEAR_LOWER_BOUND) ...
              {
                // fprintf(stderr, "%s out of range of user defined year range\n", _du_year);
                RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = true;
              }
            else
              {
                // fprintf(stderr, "Correct year %s \n", _du_year);
                duObj->year = (int64_t) _yearl;
              }
	  }

	action duration_month_copy           
	  {
	    te = p;
	    char _du_month[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_month, ts, (size_t)(te-ts));
	    _du_month[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->month = atoi(_du_month);
	  }

	action duration_day_copy             
	  {
	    te = p;
	    char _du_day[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_day, ts, (size_t)(te-ts));
	    _du_day[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->day = atoi(_du_day);     
	  }

	action duration_hour_copy            
	  {
	    te = p;
	    char _du_hour[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_hour, ts, (size_t)(te-ts));
	    _du_hour[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->hour = atoi(_du_hour);
	  }

	action duration_minute_copy          
	  {
	    te = p;
	    char _du_minute[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_minute, ts, (size_t)(te-ts));
	    _du_minute[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->minute = atoi(_du_minute);
	  }

        action duration_second_copy          
          {
            te = p;
            char _du_second[MAX_BUFFER_LENGTH] = {'\0'};
            char* _du_ms;
            int _ms;

            strncpy( _du_second, ts, (size_t)(te-ts));
            _du_second[MAX_BUFFER_LENGTH-1] = '\0';
            duObj->second = atoi(_du_second);
 
            if(strstr(_du_second,"."))
              {
                _du_ms = (strstr(_du_second,".")+1);
                if(_du_ms[0] == '-')
                  _du_ms = _du_ms + 1;
 
                _ms = atoi(_du_ms);
 
                if(strlen(_du_ms) == 1)
                  duObj->ms = _ms*100;
                else if(strlen(_du_ms) == 2)
                  duObj->ms = _ms*10;
                else
                  duObj->ms = _ms;        
              }               
          }

	#BIG: Year can be of any length but not of size 4 (which case is handeled later). String MUST be in extended format.

	Big_Year = 	
			(   
				(('+'|'$') @year_sign_match)?

				(([0-9]+'-') >_start @year_year_copy)
			);
 


	Big_rest_year =
			(	
				(( [0][1-9] | [1][0-2] ) >_start @year_month_copy)

				( '-' (	( [0][1-9] | [12][0-9] | [3][01] ) ) >_start @year_day_copy )?
			);
	
	Big_Time = ('T'|'t'|space)

			(  
				( ([01][0-9] | [2][0-3])  >_start @year_hour_copy )  
				(   
					 ( (':') ([0-5][0-9]) >_start @year_minute_copy  )
					 (  	
						( (':') ([0-5][0-9])  >_start @year_second_copy) 
						(
							( ('.'|',')([0-9]{1,3}) >_start @year_ms_copy)
						)?
					 )? 
				)? 
				(	
					('Z'|'z') >_start @duration_sign_match 
				)?
			);



	
	#Year which is exactly 4 charachters long. Non-extended (Basic) format.

	Year = 	
		(    
			(('+'|'$') @year_sign_match)?     
			( ([0-9]{4}) >_start @year_year_copy)    
		) ;

	rest_year = 
			(    
					(([0][1-9]|[1][0-2]) >_start @year_month_copy) 

					(([12][0-9]|[0][1-9]|[3][01]) >_start @year_day_copy)?  
			);

	Time = ('T'|'t'|space)

		(   
				(  	
					( ([01][0-9]|[2][0-3]) >_start @year_hour_copy )  
					( 
						( ([0-5][0-9]) >_start @year_minute_copy)  
						( 
							( ([0-5][0-9]) >_start @year_second_copy)
							(
								( ('.'|',')([0-9]{1,3}) >_start @year_ms_copy)
							)?
						)? 
					)? 
					(       
                                        	('Z'|'z') >_start @duration_sign_match
                                	)?
				) 
		);
		 

	

	#Date can not be: 1. Just Year and Month OR 2. Just Year in non-extended (basic) format.
	Date_can_not_be_YearMonth_exception = 
	  (
	    	('+'|'$')?
		( [0-9]{4} ) ( [0][1-9] | [1][0-2] )
	  ) ;
	Date_can_not_be_Year_exception = 
	  (
	   	('+'|'$')?
		( [0-9]{4} )
	  ) ;
	

	#If string contains T, t or ' ', enforce specifying month and day.
	If_T_or_space_then_full_date_only_exception = 
	  (
		('+'|'$')?  
		(
			 ( [0-9]{4} ('T'|'t'|space) any* ) 
			|
			 ( [0-9]{4}[0-9]{2}('T'|'t'|space) any* )
			|
			 ( [0-9]+'-'[0-9]{2}('T'|'t'|space) any* )
		)			  
	  ) ;

        #In Big Year, Reject YYYY...YY- (Only year ending in '-').
        If_Big_Year_ends_in_dash =
          (
                         ( any*'-' )
          ) ;



	#Negative duration is represented as -P01D instead of P-01D. Negative represented as '$'
	duration_standard_iso_type  =	
			(     
				(('$'|'+') >_start @duration_sign_match)?
				('P'|'p')
				( 
					(	
						([0-9]+ ('Y'|'y')) >_start @duration_year_copy 	
					)? 
					(
						(( [0-9] | [0][0-9] | [1][0-1] ) ('M'|'m')) >_start @duration_month_copy
					)? 
					(	
						(( [0-9] | [0][0-9] | [12][0-9] | [3][01] ) ('D'|'d')) >_start @duration_day_copy
					)?
				) 
				(
					('T'|'t') 
					(
							(	
								(( [0-9] | [01][0-9] | [2][0-3] ) ('H'|'h')) >_start @duration_hour_copy
							)? 
							(
								(( [0-9] | [0-5][0-9]) ('M'|'m')) >_start @duration_minute_copy
							)? 
							( 
								(( [0-9] | [0-5][0-9])('.'[0-9]{1,3})?('S'|'s')) >_start @duration_second_copy
							)? 
                                        ) - ('T'|'t')	
				)? 
			) - ('P'|'p');



	# Negative represented as '$'
        duration_long_day  =
                        (
                                (('$'|'+') >_start @duration_sign_match)?
                                ('P'|'p')
                                (
                                	(( [0-9]+ ) ('D'|'d')) >_start @duration_day_copy
                                )
                        );


        # Negative represented as '$'
        duration_long_hour  =
                        (
                                (('$'|'+') >_start @duration_sign_match)?
                                ('PT'|'pt')
                                (
                                	(( [0-9]+) ('H'|'h')) >_start @duration_hour_copy
                                )
                        );

        # Negative represented as '$'
        duration_long_minute  =
                        (
                                (('$'|'+') >_start @duration_sign_match)?
                                ('PT'|'pt')
                                (
                                        (( [0-9]+) ('M'|'m')) >_start @duration_minute_copy
                                )
                        );


        # Negative represented as '$'
        duration_long_second  =
                        (
                                (('$'|'+') >_start @duration_sign_match)?
                                ('PT'|'pt')
                                (
                                         (( [0-9]+)('.'[0-9]{1,3})?('S'|'s')) >_start @duration_second_copy
                                )
                        );




	
	main := 
		(
			(
		
				(       
					(	(Year (rest_year)? (Time)? '\n') 
						- 
						(Date_can_not_be_YearMonth_exception '\n') 
						- 	
						(Date_can_not_be_Year_exception '\n')
					)
					|
					((Big_Year (Big_rest_year)? (Big_Time)? '\n') - (If_Big_Year_ends_in_dash '\n') )
				) 
				- 
				(If_T_or_space_then_full_date_only_exception '\n')
		
			) @date_match
		)
		|
		(
                        (
			  (
                                (duration_long_day '\n')
                                |
                                (duration_long_hour '\n')
                                |
                                (duration_long_minute '\n')
                                |
                                (duration_long_second '\n')

			   )@duration_match_long_form
                        )
                        |
                        (
                                (duration_standard_iso_type '\n') @duration_match_iso_normal
                        )

		);              


	# Initialize and execute.
	write init;
	write exec;
    }%%
    
}

/* Internal function which calls the state machine. */
static 
ISO8601_STATUS 
get_date_time(const char* buffer, struct iso8601_datetime* datetimeObj, struct iso8601_duration* durationObj)
  {
    /* Create a local buffer and copy the string to be tested. */
    char buf[MAX_BUFFER_LENGTH] = {'\0'};
    strncpy(buf,buffer,MAX_BUFFER_LENGTH);
    buf[MAX_BUFFER_LENGTH-1] = '\0';

    /* Initialize Success or Failure flag. */
    ISO8601_STATUS stat = FAILURE;

    /* Placeholder for values of DateTime and Duration. */
    struct internal_datetime dtObj = {0};
    struct iso8601_duration duObj = {0};

    /* Initialize month and day to 1. In case these values are not specified in the buffer
       the default value should be 1. For eg. Date 1999-12 should return Year 1999, Month
       12 and DATE as "1".
    */
    dtObj.month = 1;
    dtObj.day = 1;

    /* Ragel expects \n at the end of the string. */
    char* replace = strchr(buf, '\0');
    *replace = '\n';

    /*The fact that '-' sign is used to denote negative years as well as a seperator in datetime 
      causes the regex to fail in certain scenarios. The fix ( hack? ) is to replace the - sign 
      with a '$' sign and copy back the value after processing. 
    */
    if(buf[0] == '-')
      buf[0] = '$';


    /* Execute Ragel Machine. */
    date_machine(buf,&stat,&dtObj,&duObj);

    /* stat contains the type of match. */
    if((stat == DATETIME_MATCH) && (RAISE_YEAR_OUT_OF_BOUND_EXCEPTION == false))
      {       
	/* Set sign of year. */
	if(dtObj.sign_of_year == '$')
	  datetimeObj->sign_of_year = '-';
	else
	  datetimeObj->sign_of_year = '+';

	/* Set data. */
	datetimeObj->year 	= dtObj.year;
	datetimeObj->month  	= dtObj.month;
	datetimeObj->day 	= dtObj.day;               
	datetimeObj->hour 	= dtObj.hour;
	datetimeObj->minute 	= dtObj.minute;
	datetimeObj->second 	= dtObj.second;
	datetimeObj->ms 	= dtObj.ms;
      }               
    else if((stat == DURATION_MATCH_STD || stat == DURATION_MATCH_LONG) && (RAISE_YEAR_OUT_OF_BOUND_EXCEPTION == false))
      {
       if (stat == DURATION_MATCH_STD)    /* STD: eg. P01Y05M */
         durationObj->flag_std_form = 1;
       else				  /*LONG: eg. P17M    */
         durationObj->flag_std_form = 0;
  
	/* Set sign of duration. */
	if (duObj.sign == '$')
	  durationObj->sign = '-';
	else
	  durationObj->sign = '+';

	/* Set rest. */
	durationObj->year   = duObj.year;
	durationObj->month  = duObj.month;
	durationObj->day    = duObj.day;
	durationObj->hour   = duObj.hour;
	durationObj->minute = duObj.minute;
	durationObj->second = duObj.second;
	durationObj->ms     = duObj.ms;
      }
    else
      {
        stat = FAILURE;
      }

    return stat;
  }


/*Check DateTime string compliance and get DateTime values. */
ISO8601_STATUS 
verify_string_datetime(const char* test_string,struct iso8601_datetime* dummy_isoDtObj)
  {
    ISO8601_STATUS stat = FAILURE;
    struct iso8601_duration* dummy_isoDObj = new_iso8601_duration('+',0,0,0,0,0,0,0);
    if (dummy_isoDObj == NULL)
      return FAILURE;

    stat = get_date_time(test_string, dummy_isoDtObj, dummy_isoDObj);

    deallocate_iso8601_duration(dummy_isoDObj);

    return stat;
  }

/*Check TimeDelta string compliance and get duration values.*/
ISO8601_STATUS
verify_string_duration(const char* test_string, struct iso8601_duration* dummy_isoDObj)
  {
    ISO8601_STATUS stat = FAILURE;
    struct iso8601_datetime* dummy_isoDtObj = new_iso8601_datetime('+',0,0,0,0,0,0,0);
    if ( dummy_isoDtObj == NULL)
      return FAILURE;

    stat = get_date_time(test_string, dummy_isoDtObj, dummy_isoDObj);

    deallocate_iso8601_datetime(dummy_isoDtObj);

    return stat;
  }


struct iso8601_datetime* 
new_iso8601_datetime( char _sign_of_year, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms)
  {
    struct iso8601_datetime* isoDtObj = (struct iso8601_datetime*)calloc(1,sizeof(struct iso8601_datetime));
    if (isoDtObj == NULL)
      return NULL;

    isoDtObj->sign_of_year 	= _sign_of_year;
    isoDtObj->year 		= _year;
    isoDtObj->month 		= _month;
    isoDtObj->day 		= _day;
    isoDtObj->hour 		= _hour;
    isoDtObj->minute 		= _minute;
    isoDtObj->second 		= _second;
    isoDtObj->ms 		= _ms;

    return isoDtObj;
  }


void 
deallocate_iso8601_datetime(struct iso8601_datetime* iso8601_datetimeObj)
  {
    if ( iso8601_datetimeObj != NULL)
      {
	free(iso8601_datetimeObj);
	iso8601_datetimeObj = NULL;
      }
  }

struct iso8601_duration* 
new_iso8601_duration( char _sign, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms )
  {
    struct iso8601_duration* isoDObj = (struct iso8601_duration*)calloc(1,sizeof(struct iso8601_duration));
    if (isoDObj == NULL)
      return NULL;

    isoDObj->sign 	= _sign;
    isoDObj->year 	= _year;
    isoDObj->month 	= _month;
    isoDObj->day 	= _day;
    isoDObj->hour 	= _hour;
    isoDObj->minute 	= _minute;
    isoDObj->second 	= _second;
    isoDObj->ms 	= _ms;

    return isoDObj;
  }

void 
deallocate_iso8601_duration(struct iso8601_duration* iso8601_durationObj)
  {
    if ( iso8601_durationObj != NULL )
      {
	free(iso8601_durationObj);
	iso8601_durationObj = NULL;
      }
  }

/*! \endcond */
