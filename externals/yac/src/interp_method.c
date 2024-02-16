/**
 * @file interp_method.c
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include "interp_method.h"
#include "interp_weights.h"

struct interp_weights * yac_interp_method_do_search(
  struct interp_method ** method, struct interp_grid * interp_grid) {

  size_t num_src_fields = yac_interp_grid_get_num_src_fields(interp_grid);
  enum yac_location src_field_locations[num_src_fields];
  for (size_t i = 0; i < num_src_fields; ++i)
    src_field_locations[i] =
      yac_interp_grid_get_src_field_location(interp_grid, i);
  struct interp_weights * weights =
    yac_interp_weights_new(
      yac_interp_grid_get_MPI_Comm(interp_grid),
      yac_interp_grid_get_tgt_field_location(interp_grid),
      src_field_locations, num_src_fields);

  if (*method == NULL) return weights;

  size_t temp_count;
  size_t * tgt_points;
  yac_interp_grid_get_tgt_points(interp_grid, &tgt_points, &temp_count);

  size_t final_count = 0;
  while (*method != NULL) {
    final_count +=
      (*method)->vtable->do_search(
        *method, interp_grid, tgt_points + final_count, temp_count - final_count,
        weights);
    ++method;
  }

  free(tgt_points);

  return weights;
}

void yac_interp_method_delete(struct interp_method ** method) {

  while (*method != NULL) {
    struct interp_method * curr_method = *method;
    curr_method->vtable->delete(curr_method);
    ++method;
  }
}
