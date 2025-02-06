/**
 * @file test_distributed_array_c_mp.c
 * @brief Test distributed data structure of multiple global arrays, C version
 *
 * @copyright Copyright  (C)  2014  Thomas Jahns <jahns@dkrz.de>
 *                                  Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
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
#include <math.h>
#include <stdio.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_extents_mp.h"
#include "core/ppm_random.h"
#include "core/ppm_xfuncs.h"
#include "ppm/dist_array.h"
#include "ppm/ppm.h"

#ifdef NDEBUG
#undef NDEBUG
#endif

#define eps (1.0e-13)

static void
test_array_data(struct PPM_dist_mult_array *dm_array,
                int comm_size, int comm_rank);

int main(int argc, char **argv)
{
  enum {
    num_arrays = 2,
    max_rank = 2,
    cache_mode_auto = 0,
    cache_mode_all = 1
  };

  /* init mpi */
  xmpi(MPI_Init(&argc, &argv));
  /* init scales ppm */
  unsigned int random_seed;
  PPM_initialize(NULL, NULL, &random_seed);
  /* find rank and number of tasks */
  int comm_rank, comm_size;
  xmpi(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank));
  xmpi(MPI_Comm_size(MPI_COMM_WORLD, &comm_size));
  fprintf(stderr, "%d: random_seed=%u\n", comm_rank, random_seed);
  /* initialize distributed array */
  struct PPM_global_array_desc glob_spec[num_arrays]
    = { { .a_rank = 1, .rect = { { 0, comm_size } }, MPI_INT },
        { .a_rank = 2, .rect = { { 0, comm_size }, { 0, 2 } }, MPI_DOUBLE } };
  struct PPM_extent local_array[num_arrays][max_rank]
    = { { { comm_rank, 1 } },
        { { comm_rank, 1 }, { 0, 2 } } };
  {
    /* test local data only */
    struct PPM_dist_mult_array *dm_array = PPM_dist_mult_array_new(
      2, glob_spec, local_array[0], MPI_COMM_WORLD, 0,
      PPM_dma_sync_mode_local_only);
    test_array_data(dm_array, comm_size, comm_rank);
    PPM_dist_mult_array_delete(dm_array);
  }
  for (size_t cache_mode = cache_mode_auto;
       cache_mode <= cache_mode_all;
       ++cache_mode)
  {
    /* run with little caching at first ... */
    /* ... and no cache-evictions second */
    struct PPM_dist_mult_array *dm_array = PPM_dist_mult_array_new(
      2, glob_spec, local_array[0], MPI_COMM_WORLD,
      (cache_mode == cache_mode_auto) ? 0 : (size_t)comm_size,
      PPM_dma_sync_mode_passive_target);
    test_array_data(dm_array, comm_size, comm_rank);
    PPM_dist_mult_array_delete(dm_array);
  }
  {
    /* next try passive target mode */
    struct PPM_dist_mult_array *dm_array = PPM_dist_mult_array_new(
      2, glob_spec, local_array[0], MPI_COMM_WORLD, 0,
      PPM_dma_sync_mode_active_target);
    test_array_data(dm_array, comm_size, comm_rank);
    PPM_dist_mult_array_delete(dm_array);
  }
  /* next try tests with dm_array being switched from one mode to the next */
  {
    /* test local data only */
    struct PPM_dist_mult_array *dm_array = PPM_dist_mult_array_new(
      2, glob_spec, local_array[0], MPI_COMM_WORLD, 0,
      PPM_dma_sync_mode_local_only);
    test_array_data(dm_array, comm_size, comm_rank);
    /* switch local -> passive target */
    for (size_t cache_mode = cache_mode_auto;
         cache_mode <= cache_mode_all;
         ++cache_mode)
    {
      /* run with little caching at first ... */
      /* ... and no cache-evictions second */
      PPM_dist_mult_array_set_sync_mode(
        dm_array,
        PPM_dma_sync_mode_passive_target,
        (cache_mode == cache_mode_auto) ? 0 : (size_t)comm_size);
      test_array_data(dm_array, comm_size, comm_rank);
    }
    /* switch passive -> active */
    PPM_dist_mult_array_set_sync_mode(
      dm_array, PPM_dma_sync_mode_active_target, 0);
    test_array_data(dm_array, comm_size, comm_rank);
    /* switch active -> local */
    PPM_dist_mult_array_set_sync_mode(
      dm_array, PPM_dma_sync_mode_local_only, 0);
    test_array_data(dm_array, comm_size, comm_rank);
    /* switch local -> active */
    PPM_dist_mult_array_set_sync_mode(
      dm_array, PPM_dma_sync_mode_active_target, 0);
    test_array_data(dm_array, comm_size, comm_rank);
    /* switch active -> passive */
    for (size_t cache_mode = cache_mode_auto;
         cache_mode <= cache_mode_all;
         ++cache_mode)
    {
      /* run with little caching at first ... */
      /* ... and no cache-evictions second */
      PPM_dist_mult_array_set_sync_mode(
        dm_array,
        PPM_dma_sync_mode_passive_target,
        (cache_mode == cache_mode_auto) ? 0 : (size_t)comm_size);
      test_array_data(dm_array, comm_size, comm_rank);
    }
    /* switch passive -> local */
    PPM_dist_mult_array_set_sync_mode(
      dm_array, PPM_dma_sync_mode_local_only, 0);
    test_array_data(dm_array, comm_size, comm_rank);
    PPM_dist_mult_array_delete(dm_array);
  }
  PPM_finalize();
  xmpi(MPI_Finalize());
}

static void
test_array_data(struct PPM_dist_mult_array *dm_array,
                int comm_size, int comm_rank)
{
  int sync_mode = PPM_dist_mult_array_get_sync_mode(dm_array);
  /* fill local portion of distributed array */
  int *local_i_ref = PPM_dist_mult_array_local_ptr(dm_array, 0);
  local_i_ref[0] = comm_rank;
  /* make local part available */
  PPM_dist_mult_array_expose(dm_array);
  /* get values from other ranks */
  if (sync_mode == PPM_dma_sync_mode_passive_target)
  {
    int values_i;
    for (size_t glob_idx = 0; glob_idx < (size_t)comm_size; ++glob_idx)
    {
      int32_t coord = (int)glob_idx;
      PPM_dist_mult_array_get(dm_array, 0, &coord, &values_i);
      assert(values_i == (int)glob_idx);
    }
  }
  else if (sync_mode == PPM_dma_sync_mode_active_target)
  {
    int *values_i = xmalloc((size_t)comm_size * sizeof (*values_i));
    for (size_t glob_idx = 0; glob_idx < (size_t)comm_size; ++glob_idx)
    {
      int32_t coord = (int)glob_idx;
      PPM_dist_mult_array_get(dm_array, 0, &coord, values_i + glob_idx);
    }
    PPM_dist_mult_array_rma_sync(dm_array);
    for (size_t glob_idx = 0; glob_idx < (size_t)comm_size; ++glob_idx)
      assert(values_i[glob_idx] == (int)glob_idx);
    free(values_i);
  }
  else if (sync_mode == PPM_dma_sync_mode_local_only)
  {
    int values_i;
    int32_t coord[1] = { (int32_t)comm_rank };
    PPM_dist_mult_array_get(dm_array, 0, coord, &values_i);
    assert(values_i == comm_rank);
  }
  /* second phase: also use double array */
  PPM_dist_mult_array_unexpose(dm_array);
  double *local_d_ref
    = PPM_dist_mult_array_local_ptr(dm_array, 1);
  local_d_ref[0] = log10((double)(comm_rank + 3));
  local_d_ref[1] = exp((double)(comm_rank));
  PPM_dist_mult_array_expose(dm_array);
  if (sync_mode == PPM_dma_sync_mode_passive_target)
  {
    double values_d[2];
    for (size_t glob_idx = 0; glob_idx < (size_t)comm_size; ++glob_idx)
    {
      int32_t coord[2] = { (int32_t)glob_idx, 0 };
      PPM_dist_mult_array_get(dm_array, 1, coord, &values_d[0]);
      coord[1] = 1;
      PPM_dist_mult_array_get(dm_array, 1, coord, &values_d[1]);
      assert(fabs(values_d[0] - log10((double)(glob_idx + 3))) < eps
             && fabs(values_d[1] - exp((double)(glob_idx))) < eps);
    }
  }
  else if (sync_mode == PPM_dma_sync_mode_active_target)
  {
    double (*values_d)[2] = xcalloc((size_t)comm_size, sizeof (*values_d));
    for (size_t glob_idx = 0; glob_idx < (size_t)comm_size; ++glob_idx)
    {
      int32_t coord[2] = { (int32_t)glob_idx, 0 };
      PPM_dist_mult_array_get(dm_array, 1, coord, &values_d[glob_idx][0]);
      coord[1] = 1;
      PPM_dist_mult_array_get(dm_array, 1, coord, &values_d[glob_idx][1]);
    }
    PPM_dist_mult_array_rma_sync(dm_array);
    for (size_t glob_idx = 0; glob_idx < (size_t)comm_size; ++glob_idx)
      assert(fabs(values_d[glob_idx][0] - log10((double)(glob_idx + 3))) < eps
             && fabs(values_d[glob_idx][1] - exp((double)(glob_idx))) < eps);
    free(values_d);
  }
  else if (sync_mode == PPM_dma_sync_mode_local_only)
  {
    double values_d[2];
    int32_t coord[2] = { (int32_t)comm_rank, 0 };
    PPM_dist_mult_array_get(dm_array, 1, coord, &values_d[0]);
    coord[1] = 1;
    PPM_dist_mult_array_get(dm_array, 1, coord, &values_d[1]);
    assert(fabs(values_d[0] - log10((double)(comm_rank + 3))) < eps
           && fabs(values_d[1] - exp((double)(comm_rank))) < eps);
  }
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * license-markup: "doxygen"
 * End:
 */

