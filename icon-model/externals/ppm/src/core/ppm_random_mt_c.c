/**
 * @file ppm_random_mt_c.c
 * @brief multi-threaded array-filling PRNG routines
 *
 * Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords: PRNG OpenMP
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
#include <limits.h>
#include <stdlib.h>
#include "core/ppm_random.h"

/* FIXME: these are a bit ugly because OpenMP 2.5 only allows for
 * signed loop counters */

void
PPM_irand_mt_a(int *restrict a, size_t n)
{
#define RHS PPM_irand()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_irandp_mt_a(int *restrict a, size_t n)
{
#define RHS PPM_irandp()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_irandr_mt_a(int *restrict a, size_t n, struct PPM_iinterval range)
{
#define RHS PPM_irandr(range)
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_drand_mt_a(double *restrict a, size_t n)
{
#define RHS PPM_ya_fsgrandom()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_drandp_mt_a(double *restrict a, size_t n)
{
#define RHS PPM_ya_frandom()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_drandr_mt_a(double *restrict a, size_t n, struct PPM_iinterval_dp range)
{
#define RHS PPM_drandr(range)
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_frand_mt_a(float *restrict a, size_t n)
{
#define RHS PPM_ya_fsgrandomf()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_frandp_mt_a(float *restrict a, size_t n)
{
#define RHS PPM_ya_frandomf()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_frandr_mt_a(float *restrict a, size_t n, struct PPM_iinterval_sp range)
{
#define RHS PPM_frandr(range)
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_irand8_mt_a(int64_t *restrict a, size_t n)
{
#define RHS PPM_irand8()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_irandp8_mt_a(int64_t *restrict a, size_t n)
{
#define RHS PPM_irandp8()
#include "ppm_omp_assign.h"
#undef RHS
}

void
PPM_irandr8_mt_a(int64_t *restrict a, size_t n, struct PPM_iinterval64 range)
{
#define RHS PPM_irandr8(range)
#include "ppm_omp_assign.h"
#undef RHS
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
