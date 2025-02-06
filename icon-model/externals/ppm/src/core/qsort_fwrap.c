/**
 * @file qsort_fwrap.c
 * @brief wraps around Genometools qsort variant
 *
 * @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: quicksort fortran wrapper
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

#include <stdint.h>

#include "cfortran.h"
#include "core/ppm_visibility.h"
#include "core/qsort_r_api.h"

static void
PPM_qsort_r_wrapper(void *a, int n, int es,
                    void *data, PPM_CompareWithData cmp)
{
  PPM_qsort_r(a, (size_t)n, (size_t)es, data, cmp);
}

#undef ROUTINE_5
#define ROUTINE_5 (PPM_CompareWithData)(void (*)(void))
FCALLSCSUB5(PPM_qsort_r_wrapper, PPM_QSORT_R, ppm_qsort_r,
            PVOID, INT, INT, PVOID, ROUTINE)

static void
PPM_qsort_r_mt_wrapper(void *a, int n, int es,
                       void *data, PPM_CompareWithData cmp)
{
  PPM_qsort_r_mt(a, (size_t)n, (size_t)es, data, cmp);
}


FCALLSCSUB5(PPM_qsort_r_mt_wrapper, PPM_QSORT_R_MT, ppm_qsort_r_mt,
            PVOID, INT, INT, PVOID, ROUTINE)

#undef FCALLSC_QUALIFIER
#define FCALLSC_QUALIFIER PPM_DSO_INTERNAL

static void
PPM_qsort_i4_f(int32_t *a, int n, PPM_CompareWithData cmp)
{
  PPM_qsort_r(a, (size_t)n, (size_t)4, NULL, cmp);
}

#undef ROUTINE_3
#define ROUTINE_3 (PPM_CompareWithData)(void (*)(void))
FCALLSCSUB3(PPM_qsort_i4_f, PPM_QSORT_I4, ppm_qsort_i4_f,
            INTV, INT, ROUTINE)


/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */

