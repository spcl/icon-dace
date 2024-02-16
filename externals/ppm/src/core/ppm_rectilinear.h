/**
 * @file ppm_rectilinear.h
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

#include <core/ppm_extents.h>

int32_t
PPM_rlcoord2lidx_e(int ndims, const struct PPM_extent shape[ndims],
                   const int32_t coord[ndims]);

int32_t
PPM_rlcoord2lidx_i(int ndims, const struct PPM_iinterval shape[ndims],
                   const int32_t coord[ndims]);


void
PPM_lidx2rlcoord_e(int ndims, const struct PPM_extent shape[ndims],
                   int32_t idx, int32_t coord[ndims]);

void
PPM_lidx2rlcoord_i(int ndims, const struct PPM_iinterval shape[ndims],
                   int32_t idx, int32_t coord[ndims]);

int
PPM_num_neighbours_of_rect_elem_e(int ndims,
                                  const struct PPM_extent shape[ndims],
                                  const int32_t coord[ndims]);

int
PPM_num_neighbours_of_rect_elem_i(int ndims,
                                  const struct PPM_iinterval shape[ndims],
                                  const int32_t coord[ndims]);

struct PPM_rect_coord_vec
{
  int ndims, ncoord;
  /* contains ncoord entries of size ndims, each specifying one
     cartesian coordinate */
  int32_t coords[];
};

void
PPM_lidx_nb_coords_e(int ndims, const struct PPM_extent shape[ndims],
                     int32_t idx, struct PPM_rect_coord_vec *nbcoords);

void
PPM_lidx_nb_coords_i(int ndims, const struct PPM_iinterval shape[ndims],
                     int32_t idx, struct PPM_rect_coord_vec *nbcoords);

int
PPM_lidx_nb_indices_e(int ndims, const struct PPM_extent shape[ndims],
                      int32_t idx, int32_t nbidx[]);

int
PPM_lidx_nb_indices_i(int ndims, const struct PPM_iinterval shape[ndims],
                      int32_t idx, int32_t nbidx[]);

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
