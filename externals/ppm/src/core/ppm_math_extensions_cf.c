/*
 * @file ppm_math_extensions_cf.c
 * @brief Fortran wrapper for ppm_math_extensions
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords:
 * @author Thomas Jahns <jahns@dkrz.de>
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
 * Commentary:
 *
 *
 *
 * Code:
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <inttypes.h>
#include <math.h>

#include "ppm_math_extensions.h"

#include "cfortran.h"


static inline void
PPM_fpu_save_cw_f(int *fpu_cw)
{
  PPM_fpu_save_cw((uint32_t *)fpu_cw);
}

FCALLSCSUB1(PPM_fpu_save_cw_f,PPM_FPU_SAVE_CW,ppm_fpu_save_cw,PINT)

static inline void
PPM_fpu_set_precision_f(int fpu_precision, int *old_fpu_cw)
{
  PPM_fpu_set_precision((enum precision)fpu_precision, (uint32_t *)old_fpu_cw);
}

FCALLSCSUB2(PPM_fpu_set_precision_f,PPM_FPU_SET_PRECISION_C,ppm_fpu_set_precision_c,INT,PINT)

static inline void
PPM_fpu_restore_cw_f(const int fpu_cw)
{
  PPM_fpu_restore_cw((uint32_t)fpu_cw);
}

FCALLSCSUB1(PPM_fpu_restore_cw_f,PPM_FPU_RESTORE_CW,ppm_fpu_restore_cw,INT)

static inline void
PPM_fpu_set_abrupt_underflow_f(int *old_mxcsr, int abrupt_underflow)
{
  PPM_fpu_set_abrupt_underflow((uint32_t *)old_mxcsr, abrupt_underflow);
}

FCALLSCSUB2(PPM_fpu_set_abrupt_underflow_f,PPM_FPU_SET_APRUPT_UNDERFLOW_C,
            ppm_fpu_set_abrupt_underflow_c,PINT,LOGICAL)

static inline void
PPM_fpu_save_mxcsr_f(int *old_mxcsr)
{
  PPM_fpu_save_mxcsr((uint32_t *)old_mxcsr);
}

FCALLSCSUB1(PPM_fpu_save_mxcsr_f,PPM_FPU_SAVE_MXCSR,
            ppm_fpu_save_mxcsr,PINT)


static inline void
PPM_fpu_restore_mxcsr_f(int old_mxcsr)
{
  PPM_fpu_restore_mxcsr((uint32_t)old_mxcsr);
}

FCALLSCSUB1(PPM_fpu_restore_mxcsr_f,PPM_FPU_RESTORE_MXCSR,
            ppm_fpu_restore_mxcsr,INT)

static inline void
PPM_assign_nan_dp(double *v)
{
  *v = NAN;
}

FCALLSCSUB1(PPM_assign_nan_dp,PPM_PPM_ASSIGN_NAN_DP,
            ppm_assign_nan_dp,PDOUBLE)

static inline void
PPM_assign_nan_sp(float *v)
{
  *v = NAN;
}

FCALLSCSUB1(PPM_assign_nan_sp,PPM_PPM_ASSIGN_NAN_SP,
            ppm_assign_nan_sp,PFLOAT)


/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
