/**
 * @file ppm_math_extensions_ddp_mp_c.c
 * @brief C low-level functions required for ppm_math_extensions
 * DDP summation functionality for distributed systems
 *
 * Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
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
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "core/ppm_math_extensions.h"
#include "core/core.h"
#include "core/ppm_std_type_kinds_mp.h"

MPI_Op PPM_ddpdd_sum_op = MPI_OP_NULL;

static void
ddpdd(void *a_, void *b_, int *len_, MPI_Datatype *itype)
{
  size_t len = *len_ > 0 ? (size_t)*len_ : 0;
  if (*itype == PPM_DT_C_DOUBLE_PAIR_MP)
  {
    double complex *pa = a_, *pb = b_;
    for (size_t i = 0; i < len; ++i)
      pb[i] = PPM_ddp_add_ddp_ddp(pa[i], pb[i]);
  }
  else
    PPM_abort(PPM_default_comm, "invalid MPI data type in double-double"
              " summation", __FILE__, __LINE__);
}

void
PPM_initialize_math_extensions_mp(void)
{
#pragma omp master
  MPI_Op_create(ddpdd, 1, &PPM_ddpdd_sum_op);
#pragma omp barrier
}

void
PPM_finalize_math_extensions_mp(void)
{
#pragma omp master
  MPI_Op_free(&PPM_ddpdd_sum_op);
#pragma omp barrier
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
