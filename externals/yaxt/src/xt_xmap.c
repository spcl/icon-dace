/**
 * @file xt_xmap.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_xmap.h"
#include "xt_xmap_internal.h"

MPI_Comm xt_xmap_get_communicator(Xt_xmap xmap) {

  return xmap->vtable->get_communicator(xmap);
}

int xt_xmap_get_num_destinations(Xt_xmap xmap) {

  return xmap->vtable->get_num_destinations(xmap);
}

int xt_xmap_get_num_sources(Xt_xmap xmap) {

  return xmap->vtable->get_num_sources(xmap);
}

void xt_xmap_get_destination_ranks(Xt_xmap xmap, int * ranks) {

  xmap->vtable->get_destination_ranks(xmap, ranks);
}

void xt_xmap_get_source_ranks(Xt_xmap xmap, int * ranks) {

  xmap->vtable->get_source_ranks(xmap, ranks);
}

Xt_xmap xt_xmap_copy(Xt_xmap xmap) {

  return xmap->vtable->copy(xmap);
}

void xt_xmap_delete(Xt_xmap xmap) {

  xmap->vtable->delete(xmap);
}

Xt_xmap_iter xt_xmap_get_in_iterator(Xt_xmap xmap) {

  return xmap->vtable->get_in_iterator(xmap);
}

Xt_xmap_iter xt_xmap_get_out_iterator(Xt_xmap xmap) {

  return xmap->vtable->get_out_iterator(xmap);
}

int xt_xmap_iterator_next(Xt_xmap_iter iter) {

  return iter->vtable->next(iter);
}

int xt_xmap_iterator_get_rank(Xt_xmap_iter iter) {

  return iter->vtable->get_rank(iter);
}

int const * xt_xmap_iterator_get_transfer_pos(Xt_xmap_iter iter) {

  return iter->vtable->get_transfer_pos(iter);
}

int xt_xmap_iterator_get_num_transfer_pos(Xt_xmap_iter iter) {

  return iter->vtable->get_num_transfer_pos(iter);
}

const struct Xt_pos_ext *
xt_xmap_iterator_get_transfer_pos_ext(Xt_xmap_iter iter) {
  return iter->vtable->get_transfer_pos_ext(iter);
}

int xt_xmap_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter) {
  return iter->vtable->get_num_transfer_pos_ext(iter);
}

void xt_xmap_iterator_delete(Xt_xmap_iter iter) {

  iter->vtable->delete(iter);
}

int xt_xmap_get_max_src_pos(Xt_xmap xmap) {
  return xmap->vtable->get_max_src_pos(xmap);
}

int xt_xmap_get_max_dst_pos(Xt_xmap xmap) {
  return xmap->vtable->get_max_dst_pos(xmap);
}

Xt_xmap xt_xmap_reorder(Xt_xmap xmap, enum xt_reorder_type type) {
  return xmap->vtable->reorder(xmap, type);
}

Xt_xmap xt_xmap_update_positions(
  Xt_xmap xmap, const int * src_positions, const int * dst_positions) {
  return xmap->vtable->update_pos(xmap, src_positions, dst_positions);
}

Xt_xmap xt_xmap_spread(Xt_xmap xmap, int num_repetitions,
                       const int src_displacements[num_repetitions],
                       const int dst_displacements[num_repetitions]) {
  return
    xmap->vtable->spread(
      xmap, num_repetitions, src_displacements, dst_displacements);
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
