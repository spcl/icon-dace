/**
 * @file ppm_math_extensions_ddp_c.c
 * @brief C low-level functions required for ppm_math_extensions
 * DDP summation functionality
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
#include "config.h"
#endif
#include <inttypes.h>
#include <complex.h>

#include "ppm_visibility.h"
#include "ppm_math_extensions.h"
#include "xpfpa_func.h"
#include "ppm_fpu_underflow.h"

#pragma STDC CX_LIMITED_RANGE ON
#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER == 9999 && __INTEL_COMPILER_BUILD_DATE == 20110811
#pragma float_control(precise, on)
#elif __INTEL_COMPILER < 1400 || __INTEL_COMPILER >= 1600
#pragma float_control(precise, on)
#else
#pragma GCC optimize ("-fp-model=source")
#endif
#endif

double complex PPM_DSO_API_EXPORT
PPM_ddp_sum_dp(size_t n, const double *a)
{
#ifdef NEED_PRECISION_CONTROL
  uint32_t old_fpu_cw;
  xpfpa_switch_double(&old_fpu_cw);
#endif
#ifdef NEED_UNDERFLOW_CONTROL
  uint32_t old_mxcsr;
  PPM_ENABLE_DENORMALS(&old_mxcsr);
#endif
  double cr = 0.0, ci = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    double t1 = a[i] + cr,
      e = t1 - a[i],
      t2 = ((cr - e) + (a[i] - (t1 - e))) + ci;
    cr = t1 + t2;
    ci = t2 - ((t1 + t2) - t1);
  }
  double complex s = cr + ci * I;
#ifdef NEED_UNDERFLOW_CONTROL
  PPM_RESTORE_MXCSR(&old_mxcsr);
#endif
#ifdef NEED_PRECISION_CONTROL
  xpfpa_restore(old_fpu_cw);
#endif
  return s;
}

double complex PPM_DSO_API_EXPORT
PPM_ddp_add_dp_dp(double a, double b)
{
#ifdef NEED_PRECISION_CONTROL
  uint32_t old_fpu_cw;
  xpfpa_switch_double(&old_fpu_cw);
#endif
#ifdef NEED_UNDERFLOW_CONTROL
  uint32_t old_mxcsr;
  PPM_ENABLE_DENORMALS(&old_mxcsr);
#endif
  double t1 = a + b,
    e = t1 - a,
    t2 = (b - e) + (a - (t1 - e)),
    cr = t1 + t2,
    ci = t2 - ((t1 + t2) - t1);
  double complex s = cr + ci * I;
#ifdef NEED_UNDERFLOW_CONTROL
  PPM_RESTORE_MXCSR(&old_mxcsr);
#endif
#ifdef NEED_PRECISION_CONTROL
  xpfpa_restore(old_fpu_cw);
#endif
  return s;
}

double complex PPM_DSO_API_EXPORT
PPM_ddp_add_ddp_ddp(double complex a, double complex b)
{
#ifdef NEED_PRECISION_CONTROL
  uint32_t old_fpu_cw;
  xpfpa_switch_double(&old_fpu_cw);
#endif
#ifdef NEED_UNDERFLOW_CONTROL
  uint32_t old_mxcsr;
  PPM_ENABLE_DENORMALS(&old_mxcsr);
#endif
  double ar = creal(a), br = creal(b),
    t1 = ar + br,
    e = t1 - ar,
    t2 = (br - e) + (ar - (t1 - e)) + cimag(a) + cimag(b),
    cr = t1 + t2,
    ci = t2 - ((t1 + t2) - t1);
  double complex s = cr + ci * I;
#ifdef NEED_UNDERFLOW_CONTROL
  PPM_RESTORE_MXCSR(&old_mxcsr);
#endif
#ifdef NEED_PRECISION_CONTROL
  xpfpa_restore(old_fpu_cw);
#endif
  return s;
}



/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
