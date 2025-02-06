/**
 * @file ppm_uniform_partition_c.c
 * @brief Functions for uniform partitioning of rectilinears
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

#include <inttypes.h>

#include <core/ppm_extents.h>
#include <core/ppm_xfuncs.h>
#include <ppm/ppm_uniform_partition.h>

static inline int32_t
PPM_uniform_partition_start(struct PPM_extent set_interval, int nparts,
                            int part_idx)
{
  int32_t part_offset
    = (int32_t)(((int64_t)PPM_extent_size(set_interval)
                 * (int64_t)part_idx) / (int64_t)nparts);
  return PPM_extent_start(set_interval) + part_offset;
}

static inline int32_t
PPM_uniform_partition_start_sym(struct PPM_extent set_interval, int nparts,
                                int part_idx)
{
  int64_t sym_size = ((int64_t)(PPM_extent_size(set_interval))
                      * (int64_t)(nparts/2)) / (int64_t)nparts;
  int32_t sym_part_idx
    = (part_idx >= (nparts + 1)/2) ? (nparts - part_idx) : (part_idx);
  int32_t part_offset
    = (int32_t)((sym_size * (int64_t)sym_part_idx) / (int64_t)(nparts/2));
  if (part_idx >= (nparts + 1)/2)
  {
    part_offset = PPM_extent_size(set_interval) - part_offset;
  }
  return PPM_extent_start(set_interval) + part_offset;
}

struct PPM_extent
PPM_uniform_partition(struct PPM_extent set_interval, int nparts,
                      int part_idx)
{
  struct PPM_extent part;
  part.first = PPM_uniform_partition_start(set_interval, nparts, part_idx);
  part.size = PPM_uniform_partition_start(set_interval, nparts,
                                          part_idx + 1) - part.first;
  return part;
}

struct PPM_extent
PPM_uniform_partition_symmetric(struct PPM_extent set_interval, int nparts,
                                int part_idx)
{
  struct PPM_extent part;
  part.first = PPM_uniform_partition_start_sym(set_interval, nparts, part_idx);
  part.size = PPM_uniform_partition_start_sym(set_interval, nparts,
                                              part_idx + 1) - part.first;
  return part;
}

void
PPM_uniform_partition_nd(int ndims, const struct PPM_extent set_interval[ndims],
                         const int nparts[ndims], const int part_idx[ndims],
                         struct PPM_extent part_interval[ndims])
{
  for (int i = 0; i < ndims; ++i)
    part_interval[i] = PPM_uniform_partition(set_interval[i], nparts[i],
                                             part_idx[i]);
}

void
PPM_uniform_partition_symmetric_nd(
  int ndims, const struct PPM_extent set_interval[ndims],
  int nparts[ndims], int part_idx[ndims], int symmetry[ndims],
  struct PPM_extent part_interval[ndims])
{
  for (int i = 0; i < ndims; ++i)
    if (symmetry[i])
      part_interval[i] = PPM_uniform_partition_symmetric(
        set_interval[i], nparts[i], part_idx[i]);
    else
      part_interval[i] = PPM_uniform_partition(set_interval[i],
                                               nparts[i], part_idx[i]);
}

void
PPM_uniform_decomposition_1d(struct PPM_extent set_interval, int nparts,
                             struct PPM_extent parts[nparts],
                             int symmetric)
{
  if (symmetric)
    for (int i = 0; i < nparts; ++i)
      parts[i] = PPM_uniform_partition_symmetric(set_interval, nparts, i);
  else
    for (int i = 0; i < nparts; ++i)
      parts[i] = PPM_uniform_partition(set_interval, nparts, i);
}

void
PPM_uniform_decomposition_nd(int ndims,
                             struct PPM_block_decomposition pgrid[ndims],
                             const struct PPM_extent set_interval[ndims],
                             const int nparts[ndims],
                             const int *symmetric)
{
  if (symmetric)
  {
    for (int i = 0; i < ndims; ++i)
    {
      pgrid[i].nparts = nparts[i];
      pgrid[i].partition = xmalloc(sizeof (pgrid[i].partition[0])
                                   * (size_t)nparts[i]);
      PPM_uniform_decomposition_1d(set_interval[i], nparts[i],
                                   pgrid[i].partition, symmetric[i]);
    }
  }
  else
    for (int i = 0; i < ndims; ++i)
    {
      pgrid[i].nparts = nparts[i];
      pgrid[i].partition = xmalloc(sizeof (pgrid[i].partition[0])
                                   * (size_t)nparts[i]);
      PPM_uniform_decomposition_1d(set_interval[i], nparts[i],
                                   pgrid[i].partition, 0);
    }
}


/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
