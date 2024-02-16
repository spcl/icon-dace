/**
 * @file interp_method_fixed.c
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

#include "interp_method_fixed.h"

static size_t do_search_fixed(struct interp_method * method,
                              struct interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct interp_weights * weights);
static void delete_fixed(struct interp_method * method);

static struct interp_method_vtable
  interp_method_fixed_vtable = {
    .do_search = do_search_fixed,
    .delete = delete_fixed};

struct interp_method_fixed {

  struct interp_method_vtable * vtable;
  double value;
};

static size_t do_search_fixed (struct interp_method * method,
                               struct interp_grid * interp_grid,
                               size_t * tgt_points, size_t count,
                               struct interp_weights * weights) {

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(interp_grid, tgt_points, count),
    .count = count};

  yac_interp_weights_add_fixed(
    weights, &tgts, ((struct interp_method_fixed *)method)->value);

  free(tgts.data);

  return count;
}

struct interp_method * yac_interp_method_fixed_new(double value) {

  struct interp_method_fixed * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_fixed_vtable;
  method->value = value;

  return (struct interp_method*)method;
}

static void delete_fixed(struct interp_method * method) {
  free(method);
}
