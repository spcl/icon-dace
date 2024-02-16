/**
 * @file ieee_arithmetic_emulation.c
 *
 * @brief define C functions to emulate Fortran 2003 ieee_arithmetic
 * module functionality with C99 functions
 *
 * @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
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
#include "config.h"
#endif
#include <float.h>
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#include <math.h>

#include "cfortran.h"
#include "core.h"

static int
fIsNormal(float x)
{
  return x == 0.0f || isnormal(x);
}

FCALLSCFUN1(LOGICAL,fIsNormal,PPM_IEEE_IS_NORMAL_SP,ppm_ieee_is_normal_sp,FLOAT)

static int
dIsNormal(double x)
{
  return x == 0.0 || isnormal(x);
}

FCALLSCFUN1(LOGICAL,dIsNormal,PPM_IEEE_IS_NORMAL_DP,ppm_ieee_is_normal_dp,DOUBLE)

static int
fIsNAN(float x)
{
  return isnan(x);
}

FCALLSCFUN1(LOGICAL,fIsNAN,PPM_IEEE_IS_NAN_SP,ppm_ieee_is_nan_sp,FLOAT)

static int
dIsNAN(double x)
{
  return isnan(x);
}

FCALLSCFUN1(LOGICAL,dIsNAN,PPM_IEEE_IS_NAN_DP,ppm_ieee_is_nan_dp,DOUBLE)

enum ppm_ieee_class_type
{
  ieee_signaling_nan = 1,
  ieee_quiet_nan = 2,
  ieee_negative_inf = 3,
  ieee_negative_normal = 4,
  ieee_negative_denormal = 5,
  ieee_negative_zero = 6,
  ieee_positive_zero = 7,
  ieee_positive_denormal = 8,
  ieee_positive_normal = 9,
  ieee_positive_inf = 10,
};

static float
fIeeeValue(float PPM_UNUSED(x), enum ppm_ieee_class_type v)
{
  float r;
  long e;
#ifdef HAVE_IEEE_SIGNALING_NAN
  union
  {
    float f;
    uint32_t u;
  } snan;
#endif
  switch(v)
  {
  case ieee_signaling_nan:
#ifdef HAVE_IEEE_SIGNALING_NAN
    snan.u =  UINT32_C(0x7f800001);
    return snan.f;
#else
    PPM_abort(PPM_default_comm, "unsupported ieee capability requested",
              __FILE__, __LINE__);
#endif
#ifdef NAN
  case ieee_quiet_nan:
    return NAN;
#endif
  case ieee_negative_inf:
    return -INFINITY;
  case ieee_negative_normal:
    return -1.0f;
  case ieee_negative_denormal:
    e = 0;
    do {
      r = scalblnf(-1.0f, FLT_MIN_EXP-e);
      ++e;
    } while (r < 0.0 && isnormal(r));
    return r;
  case ieee_negative_zero:
    return -0.0f;
  case ieee_positive_zero:
    return 0.0f;
  case ieee_positive_denormal:
    e = 0;
    do {
      r = scalblnf(1.0f, FLT_MIN_EXP-e);
      ++e;
    } while (r > 0.0 && isnormal(r));
    return r;
  case ieee_positive_normal:
    return 1.0f;
  case ieee_positive_inf:
    return INFINITY;
  default:
    fprintf(stderr, "v=%d\n", v);
    PPM_abort(PPM_default_comm, "invalid ieee value requested", __FILE__,
              __LINE__);
  }
}

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
FCALLSCFUN2(FLOAT,fIeeeValue,PPM_IEEE_VALUE_SP,ppm_ieee_value_sp,FLOAT,INT)
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic pop
#endif

static double
dIeeeValue(double PPM_UNUSED(x), enum ppm_ieee_class_type v)
{
  double r;
  long e;
#ifdef HAVE_IEEE_SIGNALING_NAN
  union
  {
    double f;
    uint64_t u;
  } snan;
#endif
  switch(v)
  {
  case ieee_signaling_nan:
#ifdef HAVE_IEEE_SIGNALING_NAN
    snan.u =  UINT64_C(0x7ff0000000000001);
    return snan.f;
#else
    PPM_abort(PPM_default_comm, "unsupported ieee capability requested",
              __FILE__, __LINE__);
#endif
#ifdef NAN
  case ieee_quiet_nan:
    return NAN;
#endif
  case ieee_negative_inf:
    return -INFINITY;
  case ieee_negative_normal:
    return -1.0;
  case ieee_negative_denormal:
    e = 0;
    do {
      r = scalbln(-1.0, DBL_MIN_EXP-e);
      ++e;
    } while (r < 0.0 && isnormal(r));
    return r;
  case ieee_negative_zero:
    return -0.0;
  case ieee_positive_zero:
    return 0.0;
  case ieee_positive_denormal:
    e = 0;
    do {
      r = scalbln(1.0, DBL_MIN_EXP-e);
      ++e;
    } while (r > 0.0 && isnormal(r));
    return r;
  case ieee_positive_normal:
    return 1.0;
  case ieee_positive_inf:
    return INFINITY;
  default:
    fprintf(stderr, "v=%d\n", v);
    PPM_abort(PPM_default_comm, "invalid ieee value requested", __FILE__,
              __LINE__);
  }
}

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
FCALLSCFUN2(DOUBLE,dIeeeValue,PPM_IEEE_VALUE_DP,ppm_ieee_value_dp,DOUBLE,INT)
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic pop
#endif

#ifdef HAVE_IEEE_COPYSIGN
FCALLSCFUN2(FLOAT,copysignf,PPM_COPY_SIGN_SP,ppm_copy_sign_sp,FLOAT,FLOAT)
FCALLSCFUN2(DOUBLE,copysign,PPM_COPY_SIGN_DP,ppm_copy_sign_dp,DOUBLE,DOUBLE)
#else
static float
fCopySign(float v, float s)
{
  return v * (signbit(s) ? -1.0f : 1.0f);
}

FCALLSCFUN2(FLOAT,fCopySign,PPM_COPY_SIGN_SP,ppm_copy_sign_sp,FLOAT,FLOAT)

static double
dCopySign(double v, double s)
{
  return v * (signbit(s) ? -1.0 : 1.0);
}

FCALLSCFUN2(DOUBLE,dCopySign,PPM_COPY_SIGN_DP,ppm_copy_sign_dp,DOUBLE,DOUBLE)
#endif

static int
fSupportDenormal(float x)
{
  x = fIeeeValue(x, ieee_positive_denormal);
  return x != 0.0f;
}

FCALLSCFUN1(LOGICAL,fSupportDenormal,PPM_SUPPORT_DENORMAL_SP,ppm_support_denormal_sp,FLOAT)

static int
dSupportDenormal(double x)
{
  x = dIeeeValue(x, ieee_positive_denormal);
  return x != 0.0;
}

FCALLSCFUN1(LOGICAL,dSupportDenormal,PPM_SUPPORT_DENORMAL_DP,ppm_support_denormal_dp,DOUBLE)

static int
fSupportNaN(float PPM_UNUSED(x))
{
#ifdef NAN
  return 1;
#else
  return 0;
#endif
}

FCALLSCFUN1(LOGICAL,fSupportNaN,PPM_SUPPORT_NAN_SP,ppm_support_nan_sp,FLOAT)

static int
dSupportNaN(double PPM_UNUSED(x))
{
#ifdef NAN
  return 1;
#else
  return 0;
#endif
}

FCALLSCFUN1(LOGICAL,dSupportNaN,PPM_SUPPORT_NAN_DP,ppm_support_nan_dp,DOUBLE)

static int
fSupportInf(float PPM_UNUSED(x))
{
  /* for the time being, I know no way to detect wether INFINITY is
   * an IEEE infinity or some overflowing value, so just assume it's
   * supported. */
  return 1;
}

FCALLSCFUN1(LOGICAL,fSupportInf,PPM_SUPPORT_INF_SP,ppm_support_inf_sp,FLOAT)

static int
dSupportInf(double PPM_UNUSED(x))
{
  /* for the time being, I know no way to detect wether INFINITY is
   * an IEEE infinity or some overflowing value, so just assume it's
   * supported. */
  return 1;
}

FCALLSCFUN1(LOGICAL,dSupportInf,PPM_SUPPORT_INF_DP,ppm_support_inf_dp,DOUBLE)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
