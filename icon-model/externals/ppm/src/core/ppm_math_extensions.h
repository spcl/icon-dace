#ifndef PPM_MATH_EXTENSIONS_H
#define PPM_MATH_EXTENSIONS_H
/**
 * @file ppm_math_extensions.h
 * @brief PPM extensions for math functionality
 *
 * @copyright Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
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
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <complex.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

/* number of decimal digits needed for string identity conversions of float */
#define PPM_FLT_DECIMAL_DIG PPM_CONF_FLT_DECIMAL_DIG
/* number of decimal exponent digits needed for string identity conversions of float */
#define PPM_FLT_EXP_DECIMAL_DIG PPM_CONF_FLT_EXP_DECIMAL_DIG
/* number of decimal mantissa digits needed for string identity conversions of double */
#define PPM_DBL_DECIMAL_DIG PPM_CONF_DBL_DECIMAL_DIG
/* number of decimal exponent digits needed for string identity conversions of double */
#define PPM_DBL_EXP_DECIMAL_DIG PPM_CONF_DBL_EXP_DECIMAL_DIG

/**
 * number of characters needed to print a float
 */
/* 2 for optional - sign and decimal separator, number of mantissa
 * digits plus 2 for e exponent symbol and exponent sign plus digits
 * for exponent
 * plus 1 for weak compilers on the Fortran side, i.e. only needed
 * for Fortran interoperability */
#define PPM_FLT_DECIMAL_WIDTH (2 + PPM_FLT_DECIMAL_DIG + 2 + PPM_FLT_EXP_DECIMAL_DIG + 1)

/**
 * number of characters needed to print a double
 */
/* 2 for optional - sign and decimal separator, number of mantissa
 * digits plus 2 for e exponent symbol and exponent sign plus digits
 * for exponent
 * plus 1 for weak compilers on the Fortran side, i.e. only needed
 * for Fortran interoperability */
#define PPM_DBL_DECIMAL_WIDTH (2 + PPM_DBL_DECIMAL_DIG + 2 + PPM_DBL_EXP_DECIMAL_DIG + 1)

void
PPM_fpu_save_cw(uint32_t *fpu_cw);

void
PPM_fpu_restore_cw(const uint32_t fpu_cw);

/*
 * On x87 FPUs we might need to control precision, one of multiple
 * alternatives exists.
 */
enum precision {
  PPM_FPU_PRECISION_SP = 1,
  PPM_FPU_PRECISION_DP = 2,
  PPM_FPU_PRECISION_EP = 3,
};
void
PPM_fpu_set_precision(enum precision fpu_precision, uint32_t *old_fpu_cw);

void
PPM_fpu_save_mxcsr(uint32_t *old_mxcsr);

void
PPM_fpu_restore_mxcsr(uint32_t old_mxcsr);

void
PPM_fpu_set_abrupt_underflow(uint32_t *old_mxcsr, bool abrupt_underflow);

double complex
PPM_ddp_sum_dp(size_t n, const double *a);

double complex
PPM_ddp_add_dp_dp(double a, double b);

double complex
PPM_ddp_add_ddp_ddp(double complex a, double complex b);

#ifdef USE_MPI
extern MPI_Op PPM_ddpdd_sum_op;

void
PPM_initialize_math_extensions_mp(void);

void
PPM_finalize_math_extensions_mp(void);
#endif

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
#endif
