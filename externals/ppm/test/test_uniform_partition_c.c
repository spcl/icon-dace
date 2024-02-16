/**
 * @file test_uniform_partition_c.c
 * @brief test of uniform partitioning functions C API
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
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
#  include "config.h"
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <core/core.h>
#include <core/ppm_extents.h>
#include <core/ppm_random.h>
#include <core/ppm_xfuncs.h>
#include <ppm/ppm_set_partition_base.h>
#include <ppm/ppm_uniform_partition.h>


enum {
  MIN_DIMS = 2,
  MAX_DIMS = 9,
};

#define MAX_SIZE (INT32_C(1000000))

#define MAX(a,b)                \
  ({ __typeof (a) _a = (a);       \
    __typeof (b) _b = (b);        \
    _a > _b ? _a : _b; })

#define MIN(a,b)                \
  ({ __typeof (a) _a = (a);       \
    __typeof (b) _b = (b);        \
    _a <= _b ? _a : _b; })

static void
test1DPart(size_t numpart, struct PPM_extent parts1D[numpart],
           struct PPM_extent part_range, int symmetric);

int main()
{
  unsigned seed = PPM_ya_rand_init(MPI_COMM_WORLD, 0);
  fprintf(stderr, "random seed=%u\n", seed);
  int ndims;
  {
    static const struct PPM_iinterval dimRange = { MIN_DIMS, MAX_DIMS };
    ndims = PPM_irandr(dimRange);
  }
  struct PPM_extent domain[MAX_DIMS];
  int nparts[MAX_DIMS], symmetry[MAX_DIMS];
  struct PPM_block_decomposition pgrid[MAX_DIMS];
  for (int i = 0; i < ndims; ++i)
  {
    struct PPM_iinterval sizeRange = { 0, MAX_SIZE };
    int32_t sbase = PPM_irandr(sizeRange);
    sizeRange.last = MIN(MAX_SIZE, INT32_MAX - sbase);
    int32_t ssize = PPM_irandr(sizeRange);
    sizeRange.first = sbase;
    sizeRange.last = sbase + ssize - 1;
    nparts[i] = PPM_irandr(sizeRange);
    symmetry[i] = PPM_ya_random() & 1;
    /* can't symmetrically divide an odd number of items into an even
       number of parts */
    if (symmetry[i])
      ssize += (~(nparts[i] & 1) & (ssize & 1));
    domain[i].first = sbase;
    domain[i].size = ssize;
  }
  PPM_uniform_decomposition_nd(ndims, pgrid, domain, nparts, symmetry);
  struct PPM_extent *refPart = NULL;
  for (int i = 0; i < ndims; ++i)
  {
    refPart = xrealloc(refPart, sizeof (refPart[0]) * (size_t)nparts[i]);
    test1DPart((size_t)nparts[i], refPart, domain[i], symmetry[i]);
    assert(PPM_extents_eq((size_t)nparts[i], refPart, pgrid[i].partition));
  }
  free(refPart);
  for (int i = 0; i < ndims; ++i)
    free(pgrid[i].partition);
  return EXIT_SUCCESS;
}

static void
test1DPart(size_t nparts, struct PPM_extent parts1D[nparts],
           struct PPM_extent part_range, int symmetric)
{
  if (!nparts)
    return;
  parts1D[0] = symmetric ?
    PPM_uniform_partition_symmetric(part_range, (int)nparts, 0) :
    PPM_uniform_partition(part_range, (int)nparts, 0);
  assert(PPM_extent_start(parts1D[0]) == PPM_extent_start(part_range));
  if (symmetric)
  {
    parts1D[nparts - 1]
      = PPM_uniform_partition_symmetric(part_range, (int)nparts, (int)nparts-1);
    assert(PPM_extent_size(parts1D[0]) == PPM_extent_size(parts1D[nparts - 1]));
    size_t imax = (nparts + 1)/2;
    for (size_t i = 1; i < imax; ++i)
    {
      parts1D[i] = PPM_uniform_partition_symmetric(part_range, (int)nparts,
                                                   (int)i);
      assert(PPM_extent_end(parts1D[i - 1])
             == PPM_extent_start(parts1D[i]) - 1);
      assert(llabs((long long)PPM_extent_size(parts1D[i - 1])
                   - (long long)PPM_extent_size(parts1D[i])) < 3);
      parts1D[nparts - i - 1]
        = PPM_uniform_partition_symmetric(part_range, (int)nparts,
                                          (int)(nparts - i - 1));
      assert(PPM_extent_size(parts1D[i])
             == PPM_extent_size(parts1D[nparts - i - 1]));
    }
  }
  else
    for (size_t i = 1; i < nparts; ++i)
    {
      parts1D[i] = PPM_uniform_partition(part_range, (int)nparts, (int)i);
      assert(PPM_extent_end(parts1D[i - 1])
             == PPM_extent_start(parts1D[i]) - 1);
      assert(llabs((long long)PPM_extent_size(parts1D[i - 1])
                   - (long long)PPM_extent_size(parts1D[i])) < 2);
    }
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
