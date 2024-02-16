// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//

#line 1 "mtime_iso8601.rl"
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


#line 42 "mtime_iso8601.c"
static const char _date_machine_actions[] = {
	0, 1, 0, 1, 1, 1, 2, 1, 
	3, 1, 6, 1, 7, 1, 8, 1, 
	9, 1, 10, 1, 11, 1, 12, 1, 
	13, 1, 14, 1, 15, 1, 16, 1, 
	17, 1, 18, 2, 0, 5, 2, 0, 
	12, 2, 3, 2, 3, 4, 0, 5
	
};

static const short _date_machine_key_offsets[] = {
	0, 0, 6, 10, 13, 20, 23, 25, 
	29, 31, 33, 37, 39, 41, 46, 48, 
	53, 58, 61, 62, 64, 66, 68, 72, 
	74, 80, 82, 84, 86, 89, 92, 97, 
	101, 104, 111, 115, 124, 127, 129, 134, 
	136, 141, 143, 146, 151, 156, 165, 173, 
	181, 187, 188, 197, 203, 207, 210, 217, 
	226, 228, 232, 236, 238, 239, 246, 251, 
	258, 263, 268, 273, 276, 287, 296, 304, 
	308, 310, 314, 322, 328, 334, 337, 347, 
	353, 361, 371, 378, 387, 389, 393, 397, 
	399, 400, 409, 418, 420, 424, 428, 430, 
	435, 440, 451, 460, 469, 478, 486, 494, 
	494, 497
};

static const char _date_machine_trans_keys[] = {
	36, 43, 80, 112, 48, 57, 80, 112, 
	48, 57, 45, 48, 57, 32, 48, 49, 
	84, 116, 9, 13, 50, 48, 49, 48, 
	57, 10, 58, 90, 122, 48, 53, 48, 
	57, 10, 58, 90, 122, 48, 53, 48, 
	57, 10, 44, 46, 90, 122, 48, 57, 
	10, 90, 122, 48, 57, 10, 90, 122, 
	48, 57, 10, 90, 122, 10, 48, 51, 
	49, 57, 10, 45, 48, 51, 49, 50, 
	49, 57, 10, 32, 84, 116, 9, 13, 
	48, 57, 48, 49, 48, 50, 45, 48, 
	57, 45, 48, 57, 45, 48, 49, 50, 
	57, 45, 48, 49, 57, 45, 48, 57, 
	45, 48, 51, 49, 50, 52, 57, 45, 
	48, 49, 57, 10, 32, 45, 84, 116, 
	9, 13, 48, 57, 50, 48, 49, 48, 
	57, 10, 90, 122, 48, 53, 48, 57, 
	10, 90, 122, 48, 53, 48, 51, 45, 
	48, 57, 45, 48, 49, 50, 57, 45, 
	48, 50, 51, 57, 10, 48, 49, 50, 
	51, 84, 116, 52, 57, 68, 77, 89, 
	100, 109, 121, 48, 57, 68, 77, 89, 
	100, 109, 121, 48, 57, 68, 89, 100, 
	121, 48, 57, 10, 10, 48, 49, 50, 
	51, 84, 116, 52, 57, 68, 77, 100, 
	109, 48, 57, 68, 77, 100, 109, 10, 
	84, 116, 50, 48, 49, 51, 53, 54, 
	57, 46, 72, 77, 83, 104, 109, 115, 
	48, 57, 48, 57, 83, 115, 48, 57, 
	83, 115, 48, 57, 83, 115, 10, 46, 
	72, 77, 83, 104, 109, 115, 10, 48, 
	53, 54, 57, 46, 77, 83, 109, 115, 
	48, 57, 46, 77, 83, 109, 115, 10, 
	48, 53, 54, 57, 46, 83, 115, 48, 
	57, 46, 83, 115, 46, 72, 77, 83, 
	104, 109, 115, 48, 51, 52, 57, 46, 
	72, 77, 83, 104, 109, 115, 48, 57, 
	10, 51, 84, 116, 48, 50, 52, 57, 
	68, 100, 48, 57, 68, 100, 68, 100, 
	48, 49, 68, 77, 100, 109, 48, 49, 
	50, 57, 68, 77, 100, 109, 48, 57, 
	68, 77, 100, 109, 48, 49, 10, 84, 
	116, 68, 77, 89, 100, 109, 121, 48, 
	49, 50, 57, 68, 89, 100, 121, 48, 
	57, 68, 77, 89, 100, 109, 121, 48, 
	57, 68, 77, 89, 100, 109, 121, 48, 
	49, 50, 57, 50, 48, 49, 51, 53, 
	54, 57, 46, 72, 77, 83, 104, 109, 
	115, 48, 57, 48, 57, 83, 115, 48, 
	57, 83, 115, 48, 57, 83, 115, 10, 
	46, 72, 77, 83, 104, 109, 115, 48, 
	57, 46, 72, 77, 83, 104, 109, 115, 
	48, 57, 48, 57, 83, 115, 48, 57, 
	83, 115, 48, 57, 83, 115, 10, 48, 
	53, 54, 57, 10, 48, 53, 54, 57, 
	46, 72, 77, 83, 104, 109, 115, 48, 
	51, 52, 57, 46, 72, 77, 83, 104, 
	109, 115, 48, 57, 46, 72, 77, 83, 
	104, 109, 115, 48, 57, 10, 48, 49, 
	50, 51, 84, 116, 52, 57, 48, 49, 
	50, 51, 84, 116, 52, 57, 48, 49, 
	50, 51, 84, 116, 52, 57, 50, 48, 
	49, 50, 48, 49, 0
};

static const char _date_machine_single_lengths[] = {
	0, 4, 2, 1, 5, 1, 0, 4, 
	0, 0, 4, 0, 0, 5, 0, 3, 
	3, 3, 1, 0, 0, 2, 2, 0, 
	4, 0, 0, 0, 1, 1, 3, 2, 
	1, 3, 2, 5, 1, 0, 3, 0, 
	3, 0, 1, 1, 1, 7, 6, 6, 
	4, 1, 7, 4, 4, 3, 1, 7, 
	0, 2, 2, 2, 1, 7, 1, 5, 
	5, 1, 3, 3, 7, 7, 4, 2, 
	2, 2, 4, 4, 4, 3, 6, 4, 
	6, 6, 1, 7, 0, 2, 2, 2, 
	1, 7, 7, 0, 2, 2, 2, 1, 
	1, 7, 7, 7, 7, 6, 6, 0, 
	1, 1
};

static const char _date_machine_range_lengths[] = {
	0, 1, 1, 1, 1, 1, 1, 0, 
	1, 1, 0, 1, 1, 0, 1, 1, 
	1, 0, 0, 1, 1, 0, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 2, 1, 2, 1, 1, 1, 1, 
	1, 1, 1, 2, 2, 1, 1, 1, 
	1, 0, 1, 1, 0, 0, 3, 1, 
	1, 1, 1, 0, 0, 0, 2, 1, 
	0, 2, 1, 0, 2, 1, 2, 1, 
	0, 1, 2, 1, 1, 0, 2, 1, 
	1, 2, 3, 1, 1, 1, 1, 0, 
	0, 1, 1, 1, 1, 1, 0, 2, 
	2, 2, 1, 1, 1, 1, 1, 0, 
	1, 1
};

static const short _date_machine_index_offsets[] = {
	0, 0, 6, 10, 13, 20, 23, 25, 
	30, 32, 34, 39, 41, 43, 49, 51, 
	56, 61, 65, 67, 69, 71, 74, 78, 
	80, 86, 88, 90, 92, 95, 98, 103, 
	107, 110, 116, 120, 128, 131, 133, 138, 
	140, 145, 147, 150, 154, 158, 167, 175, 
	183, 189, 191, 200, 206, 211, 215, 220, 
	229, 231, 235, 239, 242, 244, 252, 256, 
	263, 269, 273, 278, 282, 292, 301, 308, 
	312, 315, 319, 326, 332, 338, 342, 351, 
	357, 365, 374, 379, 388, 390, 394, 398, 
	401, 403, 412, 421, 423, 427, 431, 434, 
	438, 442, 452, 461, 470, 479, 487, 495, 
	496, 499
};

static const char _date_machine_indicies[] = {
	0, 0, 3, 4, 2, 1, 5, 6, 
	2, 1, 7, 8, 1, 9, 10, 11, 
	9, 9, 9, 1, 13, 12, 1, 14, 
	1, 15, 16, 17, 17, 1, 18, 1, 
	19, 1, 15, 20, 17, 17, 1, 21, 
	1, 22, 1, 15, 23, 23, 17, 17, 
	1, 24, 1, 15, 17, 17, 25, 1, 
	15, 17, 17, 26, 1, 15, 17, 17, 
	1, 15, 1, 14, 1, 27, 1, 15, 
	28, 1, 29, 31, 30, 1, 32, 1, 
	33, 9, 9, 9, 9, 1, 32, 1, 
	32, 1, 27, 1, 7, 34, 1, 7, 
	35, 1, 7, 36, 37, 38, 1, 7, 
	38, 39, 1, 7, 38, 1, 7, 40, 
	42, 41, 38, 1, 7, 38, 43, 1, 
	45, 44, 7, 44, 44, 44, 38, 1, 
	47, 46, 1, 48, 1, 15, 17, 17, 
	49, 1, 50, 1, 15, 17, 17, 21, 
	1, 48, 1, 7, 43, 1, 7, 43, 
	38, 1, 7, 39, 38, 1, 51, 52, 
	53, 54, 55, 57, 58, 56, 1, 60, 
	61, 62, 60, 61, 62, 59, 1, 60, 
	61, 62, 60, 61, 62, 63, 1, 64, 
	62, 64, 62, 63, 1, 65, 1, 51, 
	66, 67, 68, 69, 58, 58, 70, 1, 
	72, 61, 72, 61, 71, 1, 72, 61, 
	72, 61, 1, 51, 58, 58, 1, 74, 
	73, 75, 76, 1, 77, 79, 80, 81, 
	79, 80, 81, 78, 1, 82, 1, 81, 
	81, 83, 1, 81, 81, 84, 1, 81, 
	81, 1, 51, 1, 77, 79, 80, 81, 
	79, 80, 81, 1, 51, 85, 86, 1, 
	77, 80, 81, 80, 81, 87, 1, 77, 
	80, 81, 80, 81, 1, 51, 88, 89, 
	1, 77, 81, 81, 90, 1, 77, 81, 
	81, 1, 77, 79, 80, 81, 79, 80, 
	81, 78, 87, 1, 77, 79, 80, 81, 
	79, 80, 81, 87, 1, 51, 92, 58, 
	58, 91, 93, 1, 72, 72, 94, 1, 
	72, 72, 1, 72, 72, 94, 1, 72, 
	61, 72, 61, 71, 94, 1, 72, 61, 
	72, 61, 94, 1, 72, 61, 72, 61, 
	94, 1, 95, 58, 58, 1, 60, 61, 
	62, 60, 61, 62, 59, 96, 1, 60, 
	62, 60, 62, 63, 1, 60, 61, 62, 
	60, 61, 62, 96, 1, 60, 61, 62, 
	60, 61, 62, 96, 63, 1, 98, 97, 
	99, 100, 1, 101, 103, 104, 105, 103, 
	104, 105, 102, 1, 106, 1, 105, 105, 
	107, 1, 105, 105, 108, 1, 105, 105, 
	1, 95, 1, 101, 103, 104, 105, 103, 
	104, 105, 109, 1, 110, 111, 112, 113, 
	111, 112, 113, 109, 1, 114, 1, 113, 
	113, 115, 1, 113, 113, 116, 1, 113, 
	113, 1, 95, 85, 86, 1, 95, 88, 
	89, 1, 101, 103, 104, 105, 103, 104, 
	105, 102, 117, 1, 101, 111, 104, 105, 
	111, 104, 105, 109, 1, 101, 103, 104, 
	105, 103, 104, 105, 117, 1, 51, 52, 
	53, 54, 55, 58, 57, 56, 1, 52, 
	53, 54, 55, 57, 58, 56, 1, 52, 
	53, 54, 55, 58, 57, 56, 1, 1, 
	13, 12, 1, 47, 46, 1, 0
};

static const char _date_machine_trans_targs[] = {
	2, 0, 3, 101, 102, 45, 100, 4, 
	28, 5, 20, 27, 6, 19, 7, 103, 
	8, 18, 9, 10, 11, 12, 13, 14, 
	15, 16, 17, 21, 22, 23, 25, 26, 
	24, 104, 29, 30, 31, 44, 32, 33, 
	34, 42, 43, 35, 36, 105, 37, 41, 
	38, 39, 40, 103, 46, 78, 80, 81, 
	47, 82, 54, 47, 77, 70, 50, 48, 
	49, 103, 51, 74, 75, 76, 52, 52, 
	53, 55, 68, 69, 61, 56, 61, 62, 
	65, 60, 57, 58, 59, 63, 64, 64, 
	66, 67, 67, 71, 73, 72, 72, 103, 
	79, 83, 97, 99, 89, 84, 89, 95, 
	96, 88, 85, 86, 87, 90, 91, 49, 
	49, 49, 92, 93, 94, 98
};

static const char _date_machine_trans_actions[] = {
	44, 0, 1, 0, 0, 0, 0, 9, 
	0, 0, 1, 1, 1, 1, 15, 3, 
	0, 35, 1, 17, 0, 1, 19, 0, 
	38, 21, 21, 11, 0, 1, 1, 1, 
	13, 3, 0, 9, 1, 1, 0, 11, 
	1, 1, 1, 13, 0, 3, 1, 1, 
	15, 1, 17, 5, 1, 1, 1, 1, 
	1, 0, 0, 0, 27, 25, 23, 0, 
	27, 7, 1, 1, 1, 1, 1, 0, 
	27, 1, 1, 1, 1, 0, 0, 29, 
	31, 33, 0, 0, 0, 1, 1, 0, 
	1, 1, 0, 1, 1, 1, 0, 41, 
	0, 1, 1, 1, 1, 0, 0, 29, 
	31, 33, 0, 0, 0, 0, 0, 29, 
	31, 33, 0, 0, 0, 0
};

static const int date_machine_start = 1;
static const int date_machine_first_final = 103;
static const int date_machine_error = 0;

static const int date_machine_en_main = 1;


#line 41 "mtime_iso8601.rl"



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

    
#line 325 "mtime_iso8601.c"
	{
	cs = date_machine_start;
	}

#line 330 "mtime_iso8601.c"
	{
	int _klen;
	unsigned int _trans;
	const char *_acts;
	unsigned int _nacts;
	const char *_keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_keys = _date_machine_trans_keys + _date_machine_key_offsets[cs];
	_trans = _date_machine_index_offsets[cs];

	_klen = _date_machine_single_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += (unsigned int)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _date_machine_range_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += (unsigned int)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _date_machine_indicies[_trans];
	cs = _date_machine_trans_targs[_trans];

	if ( _date_machine_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _date_machine_actions + _date_machine_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 70 "mtime_iso8601.rl"
	{
	    ts = p;
	  }
	break;
	case 1:
#line 75 "mtime_iso8601.rl"
	{
	    *stat = DATETIME_MATCH;
	  }
	break;
	case 2:
#line 80 "mtime_iso8601.rl"
	{
	    *stat = DURATION_MATCH_STD;
	  }
	break;
	case 3:
#line 85 "mtime_iso8601.rl"
	{
            *stat = DURATION_MATCH_LONG;
          }
	break;
	case 4:
#line 90 "mtime_iso8601.rl"
	{
	    dtObj->sign_of_year = (*p);
	  }
	break;
	case 5:
#line 95 "mtime_iso8601.rl"
	{
	    duObj->sign = (*p);
	  }
	break;
	case 6:
#line 100 "mtime_iso8601.rl"
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
	break;
	case 7:
#line 147 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _month[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _month, ts, (size_t)(te-ts));
	    _month[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->month = atoi(_month);
	  }
	break;
	case 8:
#line 156 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _day[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _day, ts, (size_t)(te-ts));
	    _day[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->day = atoi(_day);
	  }
	break;
	case 9:
#line 165 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _hour[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _hour, ts, (size_t)(te-ts));
	    _hour[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->hour = atoi(_hour);
	  }
	break;
	case 10:
#line 174 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _minute[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _minute, ts, (size_t)(te-ts));
	    _minute[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->minute = atoi(_minute);
	  }
	break;
	case 11:
#line 183 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _second[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _second, ts, (size_t)(te-ts));
	    _second[MAX_BUFFER_LENGTH-1] = '\0';
	    dtObj->second = atoi(_second);                
	  }
	break;
	case 12:
#line 192 "mtime_iso8601.rl"
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
	break;
	case 13:
#line 207 "mtime_iso8601.rl"
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
	break;
	case 14:
#line 247 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_month[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_month, ts, (size_t)(te-ts));
	    _du_month[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->month = atoi(_du_month);
	  }
	break;
	case 15:
#line 256 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_day[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_day, ts, (size_t)(te-ts));
	    _du_day[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->day = atoi(_du_day);     
	  }
	break;
	case 16:
#line 265 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_hour[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_hour, ts, (size_t)(te-ts));
	    _du_hour[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->hour = atoi(_du_hour);
	  }
	break;
	case 17:
#line 274 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_minute[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_minute, ts, (size_t)(te-ts));
	    _du_minute[MAX_BUFFER_LENGTH-1] = '\0';
	    duObj->minute = atoi(_du_minute);
	  }
	break;
	case 18:
#line 283 "mtime_iso8601.rl"
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
	break;
#line 663 "mtime_iso8601.c"
		}
	}

_again:
	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	_out: {}
	}

#line 545 "mtime_iso8601.rl"

    
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
