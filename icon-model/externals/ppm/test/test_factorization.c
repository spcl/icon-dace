/*
 * @file test_factorization.c
 * @brief test C part of combinatorial routines
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords: combinatorics test
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "core/ppm_combinatorics.h"
#include "core/ppm_random.h"

int
main()
{
  unsigned random_seed = PPM_ya_rand_init(MPI_COMM_WORLD, 0);
  int tid = 0;
  uint32_t n, factors[31], product, *pfactors = factors;
  int num_factors, i, j;
  fprintf(stderr, "%d: random seed=%u\n", tid, random_seed);

  for (j = 0; j < 10; ++j)
  {
    while ((n = PPM_ya_random()) <= UINT32_C(1))
      ;
    num_factors = PPM_prime_factorization_32(n, &pfactors);
    if (num_factors < 1)
    {
      PPM_abort(MPI_COMM_WORLD, "factorization failed!", __FILE__, __LINE__);
    }
    printf("%d: n = %"PRIu32" = %"PRIu32, tid, n, factors[0]);
    product = factors[0];
    for (i = 1; i < num_factors; ++i)
    {
      product *= factors[i];
      printf(" * %"PRIu32, factors[i]);
    }
    xfputc('\n', stdout);
    assert(product == n);
  }
  return EXIT_SUCCESS;
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
