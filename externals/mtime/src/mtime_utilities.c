// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/*
 * strings to parse:
 *
 * time intervals <interval>
 *
 * <interval> := <start>/<end>
 * <interval> := <start>/<duration>
 * <interval> := <duration>/<end>
 * <interval> := <duration>
 *
 * repeating time intervals
 *
 * R<n>/<interval>
 *
 * ommiting n means indefinit repetitions
 *
 *  1. split string by /
 *  2. assign components to sub parts
 *  3. need to scan the repititor with strtol with repetitor+1 as start
 *  4. get the others properly scaned by the ragel parser
 */

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

#include <mtime_utilities.h>

/**
 * @brief Split ISO 8601:2004 repeated time interval strings into base components
 *
 * @param[in] recurringTimeInterval
 *         A pointer to char. The string should contain an  ISO 8601:2004 repeated time 
 *         interval string.
 *	   A string literal can be accepted.
 * @param[out] repetitor
 *         A string. Contains the repetitor part of the input string.
 * @param[out] start
 *         A string. Contains the start date part of the input string.
 * @param[out] end
 *         A string. Contains the end date part of the input string.
 * @param[out] duration
 *         A string. Contains the duration part of the input string.
 *
 * @return n/a
 *
 */

void 
splitRepetitionString(const char *recurringTimeInterval,
		      char *repetitor, char *start, char *end, char *duration)
{
  char *separator = "/";
  char *copy = NULL;
  char *token;
  char *brkb;
  
  bool dflag = false;
  
  *repetitor = '\0';
  *start = '\0';
  *end = '\0';
  *duration = '\0';
  
  // eventually remove spaces with copy

  copy = (char *) realloc(copy, (strlen(recurringTimeInterval)+1)*sizeof(char));
  strcpy(copy, recurringTimeInterval);

  token = strtok_r(copy, separator, &brkb);
  while (token != NULL)
    {
      switch (*token) 
	{
	case 'R':
	  strcpy(repetitor, token);
	  break;
	case 'P':
	  dflag = true;
	  strcpy(duration, token);
	  break;
	default:
	  if (dflag)
	    {
	      strcpy(end, token);
	    }
	  else
	    {
	      dflag = true;
	      strcpy(start, token);
	    }
	}
      token = strtok_r(NULL, separator, &brkb);
    }

  free(copy);

  return;
}

/**
 * @brief Extract number of repetitions from repetition string part.
 *
 * @param[in] repetitionString 
 *         A pointer to char. A repetition string starting with 'R'.
 *	   A string literal can be accepted.
 *
 * @return r
 *         An int representing the number of repetitions.
 */

int getRepetitions (const char *repetitionString)
{
  const char *number = repetitionString+1;
  char *ptr;

  int r = _MTIME_INDEFINETLY;
  
  r = (int) strtol(number, &ptr, 10);

  if ((r == 0) && (errno == EINVAL))
    {
      r = _MTIME_INDEFINETLY;
    }

  return (r);
}
