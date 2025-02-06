/**
 * @file test_combinatorics_c.c
 * @brief test C part of combinatorial routines
 *
 * @copyright (C) 2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: combinatorics test
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "core/ppm_combinatorics.h"
#include "core/ppm_random.h"

int
main()
{
  unsigned random_seed = PPM_ya_rand_init(MPI_COMM_WORLD, 0);
  int tid = 0;
  uint32_t *list;
  size_t n;
  static const struct PPM_iinterval nrange = { 5, 10 };
  static const uint32_t numRunsList[]
    = { 1000, 10000, 100000, 1000000, 10000000 };

  fprintf(stderr, "%d: random seed=%u\n", tid, random_seed);

  n = (size_t)PPM_irandr(nrange);

  uint32_t (*occCount)[n];

  list = xmalloc(sizeof (list[0]) * n);
  size_t occCountSize = sizeof (uint32_t) * n * n;
  occCount = xmalloc(occCountSize);
  size_t failures;
  for (size_t repeats = 0;
       repeats < (sizeof (numRunsList)/sizeof (numRunsList[0]));
       ++repeats)
  {
    memset(occCount, 0, occCountSize);
    uint32_t numRuns = numRunsList[repeats];
    for (uint32_t j = 0; j < numRuns; ++j)
    {
      for (size_t i = 0; i < n; ++i)
        list[i] = (uint32_t)i;
      PPM_permute_randomly(list, sizeof (list[0]), n);
      for (size_t i = 0; i < n; ++i)
        ++(occCount[i][list[i]]);
    }
    failures = 0;
    size_t occCountSum = 0;
    for (size_t j = 0; j < n; ++j)
      for (size_t i = 0; i < n; ++i)
      {
        float abs_delta
          = fabsf((float)numRuns / (float)n - (float)occCount[j][i]);
        float rel_delta = abs_delta / (float)n / (float)numRuns;
        if (rel_delta > 0.01f)
          ++failures;
        occCountSum += occCount[j][i];
      }
    if (failures > 0)
    {
      fprintf(stderr, "Number of repetitions: %"PRIu32", Failures: %zu/%zu\n",
             numRuns, failures, n * n);
      for (size_t j = 0; j < n; ++j)
      {
        printf("%11"PRIu32, occCount[j][0]);
        for (size_t i = 1; i < n; ++i)
          printf(", %11"PRIu32"", occCount[j][i]);
        fputs("\n", stdout);
      }
    }
    if (occCountSum != n * numRuns)
    {
      fprintf(stderr, "n*numRuns=%zu, occCountSum=%zu\n", n*numRuns, occCountSum);
      exit(EXIT_FAILURE);
    }

  }
  free(list);
  free(occCount);
  return failures?EXIT_FAILURE:EXIT_SUCCESS;
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */

