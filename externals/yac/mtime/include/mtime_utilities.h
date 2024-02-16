/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_utilities.h
 * 
 * @brief Support for handling ISO 8601:2004 repeated time interval strings.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 * 
 * @date December 2013
 *
 */

#ifndef _MTIME_UTILITIES_H
#define _MTIME_UTILITIES_H

#define _MTIME_INDEFINETLY -1

///provides a string length for repeated time interval components.
#define MAX_REPETITION_STR_LEN 32

void 
splitRepetitionString(const char *recurringTimeInterval,
		      char *repetitor, char *start, char *end, char *duration);

int getRepetitions (const char *repetitionString);

#endif

/**
 * @}
 */

