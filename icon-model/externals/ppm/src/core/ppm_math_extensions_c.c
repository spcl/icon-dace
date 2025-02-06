/*
 * @file ppm_math_extensions_c.c
 * @brief C low-level functions required for ppm_math_extensions
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
#include <stdbool.h>

#include "ppm_math_extensions.h"
#include "xpfpa_func.h"
#include "ppm_fpu_underflow.h"

void
PPM_fpu_save_cw(uint32_t *fpu_cw)
{
  xpfpa_save(fpu_cw);
}

void
PPM_fpu_set_precision(enum precision fpu_precision, uint32_t *old_fpu_cw)
{
  switch (fpu_precision)
  {
  case PPM_FPU_PRECISION_SP:
    xpfpa_switch_single(old_fpu_cw);
    break;
  case PPM_FPU_PRECISION_DP:
    xpfpa_switch_double(old_fpu_cw);
    break;
  case PPM_FPU_PRECISION_EP:
    xpfpa_switch_double_extended(old_fpu_cw);
    break;
  }
}

void
PPM_fpu_restore_cw(const uint32_t fpu_cw)
{
  xpfpa_restore(fpu_cw);
}

void
PPM_fpu_set_abrupt_underflow(uint32_t *old_mxcsr, bool abrupt_underflow)
{
  uint32_t set_flags
    = abrupt_underflow ? 1 << PPM_FTZ_BIT | 1 << PPM_DAZ_BIT | 1 << PPM_DM_BIT
    : 1 << PPM_DM_BIT,
    clear_flags
    = abrupt_underflow ? 0U : 1 << PPM_FTZ_BIT | 1 << PPM_DAZ_BIT;
  PPM_ADJUST_MXCSR(old_mxcsr, clear_flags, set_flags);
}

void
PPM_fpu_save_mxcsr(uint32_t *old_mxcsr)
{
  PPM_SAVE_MXCSR(&old_mxcsr);
}

void
PPM_fpu_restore_mxcsr(uint32_t old_mxcsr)
{
  PPM_RESTORE_MXCSR(&old_mxcsr);
}


/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
