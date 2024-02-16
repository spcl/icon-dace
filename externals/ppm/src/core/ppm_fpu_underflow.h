/**
 * @file ppm_fpu_underflow.h
 * @brief C low-level functions required for ppm_math_extensions
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords:
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

enum {
  PPM_FTZ_BIT = 15,
  PPM_DM_BIT = 8,
  PPM_DAZ_BIT = 6
};

#ifdef HAVE_MXCSR_INLINE_ASM_X86
#define PPM_ADJUST_MXCSR(old_mxcsr, clear_flags, set_flags)              \
  do {                                                                  \
    uint32_t mxcsr;                                                     \
    __asm__ __volatile__ ("stmxcsr %0" : "=m" (*old_mxcsr));             \
    mxcsr = ((*old_mxcsr) | set_flags) & ~clear_flags;                   \
    __asm__ __volatile__ ("ldmxcsr %0" : : "m" (*&mxcsr));              \
  } while (0)

#define PPM_DISABLE_DENORMALS(old_mxcsr)                                 \
  do {                                                                  \
    uint32_t set_flags = 1 << PPM_FTZ_BIT | 1 << PPM_DAZ_BIT            \
      | 1 << PPM_DM_BIT;                                                \
    PPM_ADJUST_MXCSR((old_mxcsr), 0U, set_flags);                        \
  } while (0)

#define PPM_ENABLE_DENORMALS(old_mxcsr)                                  \
  do {                                                                  \
    uint32_t clear_flags = 1 << PPM_FTZ_BIT | 1 << PPM_DAZ_BIT,         \
    set_flags = 1 << PPM_DM_BIT;                                        \
    PPM_ADJUST_MXCSR((old_mxcsr), clear_flags, set_flags);               \
  } while (0)

#define PPM_RESTORE_MXCSR(mxcsr)                                        \
  do {                                                                  \
    __asm__ __volatile__ ("ldmxcsr %0" : : "m" (*mxcsr));               \
  } while (0)

#define PPM_SAVE_MXCSR(old_mxcsr)                                       \
  do {                                                                  \
    __asm__ __volatile__ ("stmxcsr %0" : : "m" (*old_mxcsr));           \
  } while (0)

#else
#define PPM_ADJUST_MXCSR(old_mxcsr, clear_flags, set_flags)
#define PPM_DISABLE_DENORMALS(old_mxcsr)
#define PPM_ENABLE_DENORMALS(old_mxcsr)
#define PPM_RESTORE_MXCSR(old_mxcsr)
#define PPM_SAVE_MXCSR(old_mxcsr)
#endif

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
