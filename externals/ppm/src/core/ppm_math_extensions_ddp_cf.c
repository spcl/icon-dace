/**
 * @file ppm_math_extensions_ddp_cf.c
 * @brief Fortran wrappers for ppm_math_extensions
 * DDP summation functionality
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
#include <complex.h>

#include "ppm_math_extensions.h"
#include "cfortran.h"


static inline void
PPM_ddp_sum_dp_f2c(int n, const double *a, double *s)
{
  *(double complex *)s = PPM_ddp_sum_dp((size_t)n, a);
}

FCALLSCSUB3(PPM_ddp_sum_dp_f2c,PPM_DDP_SUM_DP,ppm_ddp_sum_dp,
            INT,DOUBLEV,DOUBLEV)

static inline void
PPM_ddp_add_dp_dp_f2c(double a, double b, double *s)
{
  *(double complex *)s = PPM_ddp_add_dp_dp(a, b);
}

FCALLSCSUB3(PPM_ddp_add_dp_dp_f2c,PPM_DDP_ADD_DP_DP,ppm_ddp_add_dp_dp,
            DOUBLE,DOUBLE,DOUBLEV)

static inline void
PPM_ddp_add_ddp_ddp_f2c(const double *a, const double *b, double *s)
{
  *(double complex *)s
    = PPM_ddp_add_ddp_ddp(*(double complex *)a, *(double complex *)b);
}

FCALLSCSUB3(PPM_ddp_add_ddp_ddp_f2c,PPM_DDP_ADD_DDP_DDP,ppm_ddp_add_ddp_ddp,
            DOUBLEV,DOUBLEV,DOUBLEV)



/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
