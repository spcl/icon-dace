/**
 * @file ppm_random_cf.c
 * @brief Fortran wrapper code for PRNG
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

#include "ppm_random.h"


FCALLSCFUN0(INT,PPM_irand,PPM_IRAND,ppm_irand)

FCALLSCFUN0(INT,PPM_irandp,PPM_IRANDP,ppm_irandp)

static inline int
PPM_irandr_f(struct PPM_iinterval *range)
{
  return PPM_irandr(*range);
}

FCALLSCFUN1(INT,PPM_irandr_f,PPM_IRANDR,ppm_irandr,PVOID)

static inline void
PPM_irand_a_f(int *a, int n)
{
  assert(n >= 0);
  PPM_irand_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irand_a_f,PPM_IRAND_A,ppm_irand_a,INTV,INT)
FCALLSCSUB2(PPM_irand_a_f,PPM_IRAND_A_2D,ppm_irand_a_2d,INTVV,INT)
FCALLSCSUB2(PPM_irand_a_f,PPM_IRAND_A_3D,ppm_irand_a_3d,INTVVV,INT)

static inline void
PPM_irandp_a_f(int *a, int n)
{
  assert(n >= 0);
  PPM_irandp_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irandp_a_f,PPM_IRANDP_A,ppm_irandp_a,INTV,INT)
FCALLSCSUB2(PPM_irandp_a_f,PPM_IRANDP_A_2D,ppm_irandp_a_2d,INTVV,INT)
FCALLSCSUB2(PPM_irandp_a_f,PPM_IRANDP_A_3D,ppm_irandp_a_3d,INTVVV,INT)

static inline void
PPM_irandr_a_f(int *a, int n, struct PPM_iinterval *range)
{
  assert(n >= 0);
  PPM_irandr_a(a, (size_t)n, *range);
}

FCALLSCSUB3(PPM_irandr_a_f,PPM_IRANDR_A,ppm_irandr_a,INTV,INT,PVOID)
FCALLSCSUB3(PPM_irandr_a_f,PPM_IRANDR_A_2D,ppm_irandr_a_2d,INTVV,INT,PVOID)
FCALLSCSUB3(PPM_irandr_a_f,PPM_IRANDR_A_3D,ppm_irandr_a_3d,INTVVV,INT,PVOID)


FCALLSCFUN0(INT64,PPM_irand8,PPM_IRAND8,ppm_irand8)

FCALLSCFUN0(INT64,PPM_irandp8,PPM_IRANDP8,ppm_irandp8)

static inline int64_t
PPM_irandr8_f(struct PPM_iinterval64 *range)
{
  return PPM_irandr8(*range);
}

FCALLSCFUN1(INT64,PPM_irandr8_f,PPM_IRANDR8,ppm_irandr8,PVOID)

static inline void
PPM_irand8_a_f(int64_t *a, int n)
{
  assert(n >= 0);
  PPM_irand8_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irand8_a_f,PPM_IRAND8_A,ppm_irand8_a,INT64V,INT)
FCALLSCSUB2(PPM_irand8_a_f,PPM_IRAND8_A_2D,ppm_irand8_a_2d,INT64VV,INT)
FCALLSCSUB2(PPM_irand8_a_f,PPM_IRAND8_A_3D,ppm_irand8_a_3d,INT64VVV,INT)

static inline void
PPM_irandp8_a_f(int64_t *a, int n)
{
  assert(n >= 0);
  PPM_irandp8_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_irandp8_a_f,PPM_IRANDP8_A,ppm_irandp8_a,INT64V,INT)
FCALLSCSUB2(PPM_irandp8_a_f,PPM_IRANDP8_A_2D,ppm_irandp8_a_2d,INT64VV,INT)
FCALLSCSUB2(PPM_irandp8_a_f,PPM_IRANDP8_A_3D,ppm_irandp8_a_3d,INT64VVV,INT)

static inline void
PPM_irandr8_a_f(int64_t *a, int n, struct PPM_iinterval64 *range)
{
  assert(n >= 0);
  PPM_irandr8_a(a, (size_t)n, *range);
}


FCALLSCSUB3(PPM_irandr8_a_f,PPM_IRANDR8_A,ppm_irandr8_a,INT64V,INT,PVOID)
FCALLSCSUB3(PPM_irandr8_a_f,PPM_IRANDR8_A_2D,ppm_irandr8_a_2d,INT64VV,INT,PVOID)
FCALLSCSUB3(PPM_irandr8_a_f,PPM_IRANDR8_A_3D,ppm_irandr8_a_3d,INT64VVV,INT,PVOID)

FCALLSCFUN0(DOUBLE,PPM_ya_fsgrandom,PPM_DRAND,ppm_drand)

FCALLSCFUN0(DOUBLE,PPM_ya_frandom,PPM_DRANDP,ppm_drandp)

static inline double
PPM_drandr_f(struct PPM_iinterval_dp *range)
{
  return PPM_drandr(*range);
}

FCALLSCFUN1(DOUBLE,PPM_drandr_f,PPM_DRANDR,ppm_drandr,PVOID)

static inline void
PPM_drand_a_f(double *a, int n)
{
  assert(n >= 0);
  PPM_drand_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_drand_a_f,PPM_DRAND_A,ppm_drand_a,DOUBLEV,INT)
FCALLSCSUB2(PPM_drand_a_f,PPM_DRAND_A_2D,ppm_drand_a_2d,DOUBLEVV,INT)
FCALLSCSUB2(PPM_drand_a_f,PPM_DRAND_A_3D,ppm_drand_a_3d,DOUBLEVVV,INT)

static inline void
PPM_drandp_a_f(double *a, int n)
{
  assert(n >= 0);
  PPM_drandp_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_drandp_a_f,PPM_DRANDP_A,ppm_drandp_a,DOUBLEV,INT)
FCALLSCSUB2(PPM_drandp_a_f,PPM_DRANDP_A_2D,ppm_drandp_a_2d,DOUBLEVV,INT)
FCALLSCSUB2(PPM_drandp_a_f,PPM_DRANDP_A_3D,ppm_drandp_a_3d,DOUBLEVVV,INT)

static inline void
PPM_drandr_a_f(double *a, int n, struct PPM_iinterval_dp *range)
{
  assert(n >= 0);
  PPM_drandr_a(a, (size_t)n, *range);
}

FCALLSCSUB3(PPM_drandr_a_f,PPM_DRANDR_A,ppm_drandr_a,DOUBLEV,INT,PVOID)
FCALLSCSUB3(PPM_drandr_a_f,PPM_DRANDR_A_2D,ppm_drandr_a_2d,DOUBLEVV,INT,PVOID)
FCALLSCSUB3(PPM_drandr_a_f,PPM_DRANDR_A_3D,ppm_drandr_a_3d,DOUBLEVVV,INT,PVOID)

FCALLSCFUN0(FLOAT,PPM_frand,PPM_FRAND,ppm_frand)

FCALLSCFUN0(FLOAT,PPM_frandp,PPM_FRANDP,ppm_frandp)

static inline float
PPM_frandr_f(struct PPM_iinterval_sp *range)
{
  return PPM_frandr(*range);
}

FCALLSCFUN1(FLOAT,PPM_frandr_f,PPM_FRANDR,ppm_frandr,PVOID)

static inline void
PPM_frand_a_f(float *a, int n)
{
  assert(n >= 0);
  PPM_frand_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_frand_a_f,PPM_FRAND_A,ppm_frand_a,FLOATV,INT)
FCALLSCSUB2(PPM_frand_a_f,PPM_FRAND_A_2D,ppm_frand_a_2d,FLOATVV,INT)
FCALLSCSUB2(PPM_frand_a_f,PPM_FRAND_A_3D,ppm_frand_a_3d,FLOATVVV,INT)

static inline void
PPM_frandp_a_f(float *a, int n)
{
  assert(n >= 0);
  PPM_frandp_a(a, (size_t)n);
}

FCALLSCSUB2(PPM_frandp_a_f,PPM_FRANDP_A,ppm_frandp_a,FLOATV,INT)
FCALLSCSUB2(PPM_frandp_a_f,PPM_FRANDP_A_2D,ppm_frandp_a_2d,FLOATVV,INT)
FCALLSCSUB2(PPM_frandp_a_f,PPM_FRANDP_A_3D,ppm_frandp_a_3d,FLOATVVV,INT)

static inline void
PPM_frandr_a_f(float *a, int n, struct PPM_iinterval_sp *range)
{
  assert(n >= 0);
  PPM_frandr_a(a, (size_t)n, *range);
}

FCALLSCSUB3(PPM_frandr_a_f,PPM_FRANDR_A,ppm_frandr_a,FLOATV,INT,PVOID)
FCALLSCSUB3(PPM_frandr_a_f,PPM_FRANDR_A_2D,ppm_frandr_a_2d,FLOATVV,INT,PVOID)
FCALLSCSUB3(PPM_frandr_a_f,PPM_FRANDR_A_3D,ppm_frandr_a_3d,FLOATVVV,INT,PVOID)

static void
initIRand_f(MPI_Fint *comm_f, int *random_seed)
{
  MPI_Comm comm_c = MPI_COMM_NULL;
#if defined(USE_MPI)
  int flag = 0;
#  if defined (__xlC__) && defined (_AIX)
#  pragma omp critical
#  endif
  if (MPI_Initialized(&flag) == MPI_SUCCESS && flag)
    comm_c = MPI_Comm_f2c((MPI_Fint)*comm_f);
#else
  comm_c = *comm_f;
#endif
  *random_seed = (int)PPM_ya_rand_init(comm_c, *random_seed);
}

FCALLSCSUB2(initIRand_f,PPM_INITIALIZE_IRAND,ppm_initialize_irand,PVOID,PINT)

FCALLSCSUB0(PPM_ya_rand_finish,PPM_FINALIZE_IRAND,ppm_finalize_irand)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
