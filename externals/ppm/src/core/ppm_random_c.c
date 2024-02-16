/**
 * @file ppm_random_c.c
 * @brief C routines to use pseudo-random number generator in Fortran
 *
 * @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "core.h"
#include "yarandom.h"
#include "ppm_random.h"
#include "ppm_extents.h"

int
PPM_irand(void)
{
  uint32_t i;
  while ((i = PPM_ya_random()) == UINT32_C(-2147483648))
    ;
  return (int)i;
}


int
PPM_irandp(void)
{
  uint32_t i = PPM_ya_random() & ((UINT32_C(1)<<31) - UINT32_C(1));
  return (int)i;
}

int
PPM_irandr(struct PPM_iinterval range)
{
  uint32_t i;
  uint32_t range_size = (uint32_t)range.last + (uint32_t)INT_MIN
    - (uint32_t)range.first + UINT32_C(1) + (uint32_t)INT_MIN;
  uint32_t range_max_multiple;
  assert(range.last >= range.first);
  range_max_multiple = (((uint32_t)-1) / range_size) * range_size;
  while ((i = PPM_ya_random()) > range_max_multiple)
    ;
  i = (i % range_size) + (uint32_t)range.first;
  return (int)i;
}



void
PPM_irand_a(int *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_irand();
}


void
PPM_irandp_a(int *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_irandp();
}


void
PPM_irandr_a(int *a, size_t n, struct PPM_iinterval range)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_irandr(range);
}

int64_t
PPM_irand8(void)
{
  uint64_t i;
  while ((i = PPM_ya_random64()) == UINT64_C(-9223372036854775808))
    ;
  return (int64_t)i;
}

int64_t
PPM_irandp8(void)
{
  uint64_t i = PPM_ya_random64() & ((UINT64_C(1)<<63) - UINT64_C(1));
  return (int64_t)i;
}


int64_t
PPM_irandr8(struct PPM_iinterval64 range)
{
  uint64_t i;
  uint64_t range_size = (uint64_t)range.last + (uint64_t)INT64_MIN
    - (uint64_t)range.first + UINT64_C(1) + (uint64_t)INT64_MIN;
  uint64_t range_max_multiple;
  assert(range.last >= range.first);
  range_max_multiple = (((uint64_t)-1) / range_size) * range_size;
  while ((i = PPM_ya_random64()) > range_max_multiple)
    ;
  i = (i % range_size) + (uint64_t)range.first;
  return (int)i;
}

void
PPM_irand8_a(int64_t *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_irand8();
}

void
PPM_irandp8_a(int64_t *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_irandp8();
}

void
PPM_irandr8_a(int64_t *a, size_t n, struct PPM_iinterval64 range)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_irandr8(range);
}

double
PPM_drand(void)
{
  return PPM_ya_fsgrandom();
}

double
PPM_drandp(void)
{
  return PPM_ya_frandom();
}

double
PPM_drandr(struct PPM_iinterval_dp range)
{
  double range_size, x;
  assert(range.last >= range.first);
  range_size = range.last - range.first;
  x = range.first + (range.last != range.first ? 1.0 : 0.0)
    * PPM_ya_frandom() * nextafter(range_size, HUGE_VAL);
  x = x <= range.last ? range.last : range.first;
  return x;
}

void
PPM_drand_a(double *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_ya_fsgrandom();
}

void
PPM_drandp_a(double *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_ya_frandom();
}

void
PPM_drandr_a(double *a, size_t n, struct PPM_iinterval_dp range)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_drandr(range);
}

float
PPM_frand(void)
{
  return PPM_ya_fsgrandomf();
}

float
PPM_frandp(void)
{
  return PPM_ya_frandomf();
}

float
PPM_frandr(struct PPM_iinterval_sp range)
{
  float range_size, x;
  assert(range.last >= range.first);
  range_size = range.last - range.first;
  x = range.first + (range.last != range.first ? 1.0f : 0.0f)
    * PPM_ya_frandomf() * nextafterf(range_size, HUGE_VALF);
  x = x <= range.last ? x : range.last;
  return x;
}

void
PPM_frand_a(float *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_ya_fsgrandomf();
}

void
PPM_frandp_a(float *a, size_t n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_ya_frandomf();
}

void
PPM_frandr_a(float *a, size_t n, struct PPM_iinterval_sp range)
{
  size_t i;
  for (i = 0; i < n; ++i)
    a[i] = PPM_frandr(range);
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
