/**
 * @file ppm_rectilinear_c.c
 * @brief Functions for rectilinear data structures.
 *
 * Compute conversions from rectilinear coordinates to logical indices
 * and neighbour coordinates or indices.
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
#include <inttypes.h>
#include <string.h>

#include <core/ppm_rectilinear.h>

int32_t
PPM_rlcoord2lidx_e(int ndims, const struct PPM_extent shape[ndims],
                   const int32_t coord[ndims])
{
  int32_t idx = 0;
  /* shape and coordinate must have identical, non-negative dimensions */
  assert(ndims >= 0);

  if (ndims > 0)
  {
    idx += coord[0] - PPM_extent_start(shape[0]);
    for (size_t i = 1; i < (size_t)ndims; ++i)
      idx += PPM_extents_size(i, shape)
        * (coord[i] - PPM_extent_start(shape[i]));
  }
  return idx;
}

int32_t
PPM_rlcoord2lidx_i(int ndims, const struct PPM_iinterval shape[ndims],
                   const int32_t coord[ndims])
{
  int32_t idx = 0;

  /* rect and coordinate must have identical, non-zero dimensions */
  assert(ndims >= 0);

  if (ndims > 0)
  {
    idx += coord[0] - PPM_iinterval_start(shape[0]);
    for (size_t i = 1; i < (size_t)ndims; ++i)
      idx += PPM_iintervals_size(i, shape)
        * (coord[i] - PPM_iinterval_start(shape[i]));
  }
  return idx;
}

void
PPM_lidx2rlcoord_e(int ndims, const struct PPM_extent shape[ndims],
                   int32_t idx, int32_t coord[ndims])
{
  int32_t rest = idx;
  for (int i = ndims - 1; i > 0; --i)
  {
    int32_t base = PPM_extents_size((size_t)i, shape);
    int32_t ofs = rest / base;
    coord[i] = ofs + PPM_extent_start(shape[i]);
    rest -= ofs * base;
  }
  coord[0] = rest + PPM_extent_start(shape[0]);
}

void
PPM_lidx2rlcoord_i(int ndims, const struct PPM_iinterval shape[ndims],
                   int32_t idx, int32_t coord[ndims])
{
  int32_t rest = idx;
  for (int i = ndims - 1; i > 0; --i)
  {
    int32_t base = PPM_iintervals_size((size_t)i, shape);
    int32_t ofs = rest / base;
    coord[i] = ofs + PPM_iinterval_start(shape[i]);
    rest -= ofs * base;
  }
  coord[0] = rest + PPM_iinterval_start(shape[0]);
}

int
PPM_num_neighbours_of_rect_elem_e(int ndims,
                                  const struct PPM_extent shape[ndims],
                                  const int32_t coord[ndims])
{
  int nnb = 0;
  for(int i = 0; i < ndims; ++i)
  {
    /* coordinate must be contained in rect */
    assert(coord[i] <= PPM_extent_end(shape[i])
           && coord[i] >= PPM_extent_start(shape[i]));
    nnb += (coord[i] > PPM_extent_start(shape[i])) &
      + (coord[i] < PPM_extent_end(shape[i]));
  }
  return nnb;
}

int
PPM_num_neighbours_of_rect_elem_i(int ndims,
                                  const struct PPM_iinterval shape[ndims],
                                  const int32_t coord[ndims])
{
  int nnb = 0;
  for(int i = 0; i < ndims; ++i)
  {
    /* coordinate must be contained in rect */
    assert(coord[i] <= PPM_iinterval_end(shape[i])
           && coord[i] >= PPM_iinterval_start(shape[i]));
    nnb += (coord[i] > PPM_iinterval_start(shape[i])) &
      + (coord[i] < PPM_iinterval_end(shape[i]));
  }
  return nnb;
}

void
PPM_lidx_nb_coords_e(int ndims, const struct PPM_extent shape[ndims],
                     int32_t idx, struct PPM_rect_coord_vec *nbcoords)
{
  int k = 0;
  int32_t coord[ndims];

  assert(ndims > 0);
  PPM_lidx2rlcoord_e(ndims, shape, idx, coord);
  nbcoords->ndims = ndims;

  int (*coords)[ndims] = (int (*)[ndims])nbcoords->coords;

  for(int i = 0; i < ndims; ++i)
  {
    if (coord[i] > PPM_extent_start(shape[i]))
    {
      memcpy(coords[k], coord, sizeof (coord));
      coords[k][i] -= 1;
      ++k;
    }
    if (coord[i] < PPM_extent_end(shape[i]))
    {
      memcpy(coords[k], coord, sizeof (coord));
      coords[k][i] = coord[i] + 1;
      ++k;
    }
  }
  nbcoords->ncoord = k;
}

void
PPM_lidx_nb_coords_i(int ndims, const struct PPM_iinterval shape[ndims],
                     int32_t idx, struct PPM_rect_coord_vec *nbcoords)
{
  int k = 0;
  int32_t coord[ndims];

  assert(ndims > 0);
  PPM_lidx2rlcoord_i(ndims, shape, idx, coord);
  nbcoords->ndims = ndims;

  int (*coords)[ndims] = (int (*)[ndims])nbcoords->coords;

  for(int i = 0; i < ndims; ++i)
  {
    if (coord[i] > PPM_iinterval_start(shape[i]))
    {
      memcpy(coords[k], coord, sizeof (coord));
      coords[k][i] -= 1;
      ++k;
    }
    if (coord[i] < PPM_iinterval_end(shape[i]))
    {
      memcpy(coords[k], coord, sizeof (coord));
      coords[k][i] = coord[i] + 1;
      ++k;
    }
  }
  nbcoords->ncoord = k;
}

int
PPM_lidx_nb_indices_e(int ndims, const struct PPM_extent shape[ndims],
                      int32_t idx, int32_t nbidx[])
{
  int k = 0;
  int32_t coord[ndims];

  PPM_lidx2rlcoord_e(ndims, shape, idx, coord);

  for (int i = 0; i < ndims; ++i)
  {
    if (coord[i] > PPM_extent_start(shape[i]))
    {
      coord[i] -= 1;
      nbidx[k] = PPM_rlcoord2lidx_e(ndims, shape, coord);
      coord[i] += 1;
      ++k;
    }
    if (coord[i] < PPM_extent_end(shape[i]))
    {
      coord[i] += 1;
      nbidx[k] = PPM_rlcoord2lidx_e(ndims, shape, coord);
      coord[i] -= 1;
      ++k;
    }
  }
  return k;
}

int
PPM_lidx_nb_indices_i(int ndims, const struct PPM_iinterval shape[ndims],
                      int32_t idx, int32_t nbidx[])
{
  int k = 0;
  int32_t coord[ndims];

  PPM_lidx2rlcoord_i(ndims, shape, idx, coord);

  for (int i = 0; i < ndims; ++i)
  {
    if (coord[i] > PPM_iinterval_start(shape[i]))
    {
      coord[i] -= 1;
      nbidx[k] = PPM_rlcoord2lidx_i(ndims, shape, coord);
      coord[i] += 1;
      ++k;
    }
    if (coord[i] < PPM_iinterval_end(shape[i]))
    {
      coord[i] += 1;
      nbidx[k] = PPM_rlcoord2lidx_i(ndims, shape, coord);
      coord[i] -= 1;
      ++k;
    }
  }
  return k;
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
