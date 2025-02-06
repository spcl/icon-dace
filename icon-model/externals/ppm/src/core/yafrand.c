/**
 * @file yafrand.c
 * @brief compute uniformly distributed floating point numbers
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
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
#  include "config.h"
#endif

#define _ISOC99_SOURCE
#include <float.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "core/yarandom.h"
#include "core/ppm_random.h"

double
PPM_ya_frandom(void)
{
  uint64_t bits;
  int64_t mantissa;
  double ret;
  bits = (uint64_t)PPM_ya_random() + (((uint64_t)PPM_ya_random()) << 32);
  mantissa = (bits & ((UINT64_C(1) << DBL_MANT_DIG) - 1));
  ret = scalbln((double)(mantissa), -DBL_MANT_DIG);
  return ret;
}

float
PPM_ya_frandomf(void)
{
  uint32_t bits;
  int32_t mantissa;
  float ret;
  bits = PPM_ya_random();
  mantissa = (bits & ((UINT32_C(1) << FLT_MANT_DIG) - 1));
  ret = scalblnf((float)(mantissa), -FLT_MANT_DIG);
  return ret;
}

double
PPM_ya_fsgrandom(void)
{
  uint64_t bits = PPM_ya_random64();
  int64_t mantissa, sign;
  double ret;
  mantissa = (bits & ((UINT64_C(1) << DBL_MANT_DIG) - 1));
  sign = (int64_t)((bits >> (DBL_MANT_DIG - 1)) & 2) - 1;
  ret = scalbln((double)(sign * mantissa), -DBL_MANT_DIG);
  return ret;
}

float
PPM_ya_fsgrandomf(void)
{
  uint32_t bits;
  int32_t mantissa, sign;
  float ret;
  bits = PPM_ya_random();
  mantissa = (bits & ((UINT32_C(1) << FLT_MANT_DIG) - 1));
  sign = (int32_t)((bits >> (FLT_MANT_DIG - 1)) & 2) - 1;
  ret = scalblnf((float)(sign * mantissa), -FLT_MANT_DIG);
  return ret;
}

/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
