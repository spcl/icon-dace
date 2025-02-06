/*
 * @file bsearch_fwrap.c
 * @brief Fortran wrappers for binary search routines.
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords: binary search Fortran
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
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stddef.h>
#include <stdlib.h>

#include "core/ppm_visibility.h"
#include "cfortran.h"

#include "core/bsearch.h"

static int
PPM_bsearch_r_f(const void *key, const void *a, int n, int es,
                void *data, PPM_CompareWithData cmp)
{
  void *match = PPM_bsearch_r(key, a, (size_t)n, (size_t)es, data, cmp);
  return (match == NULL)?0:(int)(((char *)match-(char *)a)/(ptrdiff_t)es + 1);
}

#undef ROUTINE_6
#define ROUTINE_6 (PPM_CompareWithData)(void (*)(void))
FCALLSCFUN6(INT,PPM_bsearch_r_f, PPM_BSEARCH_R, ppm_bsearch_r,
            PVOID, PVOID, INT, INT, PVOID, ROUTINE)

static int
PPM_bsearch_el_r_f(const void *key, const void *a, int n, int es,
                   void *data, PPM_CompareWithData cmp)
{
  void *match = PPM_bsearch_el_r(key, a, (size_t)n, (size_t)es, data, cmp);
  return (match == NULL)?0:(int)(((char *)match-(char *)a)/(ptrdiff_t)es + 1);
}

#undef ROUTINE_6
#define ROUTINE_6 (PPM_CompareWithData)(void (*)(void))
FCALLSCFUN6(INT,PPM_bsearch_el_r_f, PPM_BSEARCH_EL_R,
            ppm_bsearch_el_r,
            PVOID, PVOID, INT, INT, PVOID, ROUTINE)

/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
