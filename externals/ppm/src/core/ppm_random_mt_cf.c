/**
 * @file ppm_random_mt_cf.c
 * @brief Fortran wrapper code for PRNG, multi-thread part
 *
 * @copyright Copyright  (C)  2021  Thomas Jahns <jahns@dkrz.de>
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
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include "core/ppm_random.h"
#include "cfortran.h"

#ifdef INT64_T_IS_LONG_LONG
#define INT64 LONGLONG
#define INT64V LONGLONGV
#define INT64VV LONGLONGVV
#define INT64VVV LONGLONGVVV
#elif defined INT64_T_IS_LONG
#define INT64 LONG
#define INT64V LONGV
#define INT64VV LONGVV
#define INT64VVV LONGVVV
#endif


static inline void
PPM_drand_mt_a_f(double *restrict a, int n)
{
  assert(n >= 0);
  PPM_drand_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_drand_mt_a_f,PPM_DRAND_MT_A,ppm_drand_mt_a,DOUBLEV,INT)
FCALLSCSUB2(PPM_drand_mt_a_f,PPM_DRAND_MT_A_2D,ppm_drand_mt_a_2d,DOUBLEVV,INT)
FCALLSCSUB2(PPM_drand_mt_a_f,PPM_DRAND_MT_A_3D,ppm_drand_mt_a_3d,DOUBLEVVV,INT)

static inline void
PPM_irand_mt_a_f(int *restrict a, int n)
{
  assert(n >= 0);
  PPM_irand_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irand_mt_a_f,PPM_IRAND_MT_A,ppm_irand_mt_a,INTV,INT)
FCALLSCSUB2(PPM_irand_mt_a_f,PPM_IRAND_MT_A_2D,ppm_irand_mt_a_2d,INTVV,INT)
FCALLSCSUB2(PPM_irand_mt_a_f,PPM_IRAND_MT_A_3D,ppm_irand_mt_a_3d,INTVVV,INT)

static inline void
PPM_irandp_mt_a_f(int *restrict a, int n)
{
  assert(n >= 0);
  PPM_irandp_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irandp_mt_a_f,PPM_IRANDP_MT_A,ppm_irandp_mt_a,INTV,INT)
FCALLSCSUB2(PPM_irandp_mt_a_f,PPM_IRANDP_MT_A_2D,ppm_irandp_mt_a_2d,INTVV,INT)
FCALLSCSUB2(PPM_irandp_mt_a_f,PPM_IRANDP_MT_A_3D,ppm_irandp_mt_a_3d,INTVVV,INT)

static inline void
PPM_irandr_mt_a_f(int *restrict a, int n, struct PPM_iinterval *restrict range)
{
  assert(n >= 0);
  PPM_irandr_mt_a(a, (size_t)n, *range);
}

FCALLSCSUB3(PPM_irandr_mt_a_f,PPM_IRANDR_MT_A,ppm_irandr_mt_a,INTV,INT,PVOID)
FCALLSCSUB3(PPM_irandr_mt_a_f,PPM_IRANDR_MT_A_2D,ppm_irandr_mt_a_2d,
            INTVV,INT,PVOID)
FCALLSCSUB3(PPM_irandr_mt_a_f,PPM_IRANDR_MT_A_3D,ppm_irandr_mt_a_3d,
            INTVVV,INT,PVOID)

static inline void
PPM_drandp_mt_a_f(double *restrict a, int n)
{
  assert(n >= 0);
  PPM_drandp_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_drandp_mt_a_f,PPM_DRANDP_MT_A,ppm_drandp_mt_a,DOUBLEV,INT)
FCALLSCSUB2(PPM_drandp_mt_a_f,PPM_DRANDP_MT_A_2D,ppm_drandp_mt_a_2d,
            DOUBLEVV,INT)
FCALLSCSUB2(PPM_drandp_mt_a_f,PPM_DRANDP_MT_A_3D,ppm_drandp_mt_a_3d,
            DOUBLEVVV,INT)

static inline void
PPM_drandr_mt_a_f(double *restrict a, int n, struct PPM_iinterval_dp *restrict range)
{
  assert(n >= 0);
  PPM_drandr_mt_a(a, (size_t)n, *range);
}

FCALLSCSUB3(PPM_drandr_mt_a_f,PPM_DRANDR_MT_A,ppm_drandr_mt_a,DOUBLEV,INT,PVOID)
FCALLSCSUB3(PPM_drandr_mt_a_f,PPM_DRANDR_MT_A_2D,ppm_drandr_mt_a_2d,
            DOUBLEVV,INT,PVOID)
FCALLSCSUB3(PPM_drandr_mt_a_f,PPM_DRANDR_MT_A_3D,ppm_drandr_mt_a_3d,
            DOUBLEVVV,INT,PVOID)

static inline void
PPM_frand_mt_a_f(float *restrict a, int n)
{
  assert(n >= 0);
  PPM_frand_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_frand_mt_a_f,PPM_FRAND_MT_A,ppm_frand_mt_a,FLOATV,INT)
FCALLSCSUB2(PPM_frand_mt_a_f,PPM_FRAND_MT_A_2D,ppm_frand_mt_a_2d,FLOATVV,INT)
FCALLSCSUB2(PPM_frand_mt_a_f,PPM_FRAND_MT_A_3D,ppm_frand_mt_a_3d,FLOATVVV,INT)

static inline void
PPM_frandp_mt_a_f(float *restrict a, int n)
{
  assert(n >= 0);
  PPM_frandp_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_frandp_mt_a_f,PPM_FRANDP_MT_A,ppm_frandp_mt_a,FLOATV,INT)
FCALLSCSUB2(PPM_frandp_mt_a_f,PPM_FRANDP_MT_A_2D,ppm_frandp_mt_a_2d,FLOATVV,INT)
FCALLSCSUB2(PPM_frandp_mt_a_f,PPM_FRANDP_MT_A_3D,ppm_frandp_mt_a_3d,
            FLOATVVV,INT)

static inline void
PPM_frandr_mt_a_f(float *restrict a, int n,
                  struct PPM_iinterval_sp *restrict range)
{
  assert(n >= 0);
  PPM_frandr_mt_a(a, (size_t)n, *range);
}

FCALLSCSUB3(PPM_frandr_mt_a_f,PPM_FRANDR_MT_A,ppm_frandr_mt_a,FLOATV,INT,PVOID)
FCALLSCSUB3(PPM_frandr_mt_a_f,PPM_FRANDR_MT_A_2D,ppm_frandr_mt_a_2d,
            FLOATVV,INT,PVOID)
FCALLSCSUB3(PPM_frandr_mt_a_f,PPM_FRANDR_MT_A_3D,ppm_frandr_mt_a_3d,
            FLOATVVV,INT,PVOID)

static inline void
PPM_irand8_mt_a_f(int64_t *restrict a, int n)
{
  assert(n >= 0);
  PPM_irand8_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irand8_mt_a_f,PPM_IRAND8_MT_A,ppm_irand8_mt_a,
            INT64V,INT)
FCALLSCSUB2(PPM_irand8_mt_a_f,PPM_IRAND8_MT_A_2D,ppm_irand8_mt_a_2d,
            INT64VV,INT)
FCALLSCSUB2(PPM_irand8_mt_a_f,PPM_IRAND8_MT_A_3D,ppm_irand8_mt_a_3d,
            INT64VVV,INT)

static inline void
PPM_irandp8_mt_a_f(int64_t *restrict a, int n)
{
  assert(n >= 0);
  PPM_irandp8_mt_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irandp8_mt_a_f,PPM_IRANDP8_MT_A,ppm_irandp8_mt_a,
            INT64V,INT)
FCALLSCSUB2(PPM_irandp8_mt_a_f,PPM_IRANDP8_MT_A_2D,ppm_irandp8_mt_a_2d,
            INT64VV,INT)
FCALLSCSUB2(PPM_irandp8_mt_a_f,PPM_IRANDP8_MT_A_3D,ppm_irandp8_mt_a_3d,
            INT64VVV,INT)

static inline void
PPM_irandr8_mt_a_f(int64_t *restrict a, int n,
                   struct PPM_iinterval64 *restrict range)
{
  assert(n >= 0);
  PPM_irandr8_mt_a(a, (size_t)n, *range);
}

FCALLSCSUB3(PPM_irandr8_mt_a_f,PPM_IRANDR8_MT_A,ppm_irandr8_mt_a,
            INT64V,INT,PVOID)
FCALLSCSUB3(PPM_irandr8_mt_a_f,PPM_IRANDR8_MT_A_2D,ppm_irandr8_mt_a_2d,
            INT64VV,INT,PVOID)
FCALLSCSUB3(PPM_irandr8_mt_a_f,PPM_IRANDR8_MT_A_3D,ppm_irandr8_mt_a_3d,
            INT64VVV,INT,PVOID)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
