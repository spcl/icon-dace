// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

/* 128 bit UUID computation / Rabin fingerprinting algorithm       */

#ifndef UTIL_UUID_H
#define UTIL_UUID_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#define UUID_STRING_LENGTH 36

typedef struct
{
  unsigned char data[16];
} uuid_t;

typedef enum {
  UUID_EQUAL,
  UUID_EQUAL_LIMITED_ACCURACY,
  UUID_UNEQUAL
} cmp_UUID_t;


typedef struct 
{
  unsigned int f64[2];   // 64 bit fingerprint
  unsigned int f48[2];   // 48 bit fingerprint

  unsigned int tl64[2];  // t^l  mod P64 (needed for parallel computation)
  unsigned int tl48[2];  // t^l  mod P48 (needed for parallel computation)

  unsigned int max_zero; // largest bit position with a zero coefficient (max over all numbers)
} context_t;


void            encode_uuid(context_t *context, uuid_t* uuid);
context_t*      uuid_scan_data(const double* val, const int nval);
context_t*      concat_fingerprints(context_t* context0, context_t *context1);
void            deallocate_fingerprint(context_t* context);

void            uuid_generate(const double* val, const int nval, uuid_t* uuid);
cmp_UUID_t      compare_UUID(const uuid_t uuid_A, const uuid_t uuid_B, double* min_difference);
void            uuid_unparse(char *buffer, const uuid_t *uuid);
int             uuid_parse(uuid_t *uuid, const char *uuid_str);

#endif /* UTIL_UUID_H */
