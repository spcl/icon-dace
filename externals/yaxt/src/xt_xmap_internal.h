/**
 * @file xt_xmap_internal.h
 * @brief contains declaration for the exchange map data structure
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
#ifndef XT_XMAP_INTERNAL_H
#define XT_XMAP_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_xmap.h"

struct Xt_xmap_iter_vtable {
  int (*next)(Xt_xmap_iter iter);
  int (*get_rank)(Xt_xmap_iter iter);
  int const * (*get_transfer_pos)(Xt_xmap_iter iter);
  int (*get_num_transfer_pos)(Xt_xmap_iter iter);
  const struct Xt_pos_ext *(*get_transfer_pos_ext)(Xt_xmap_iter iter);
  int (*get_num_transfer_pos_ext)(Xt_xmap_iter iter);
  void (*delete)(Xt_xmap_iter iter);
};

struct Xt_xmap_iter_ {
  const struct Xt_xmap_iter_vtable * vtable;
};


struct Xt_xmap_vtable {

  MPI_Comm (*get_communicator)(Xt_xmap);
  int (*get_num_destinations)(Xt_xmap);
  int (*get_num_sources)(Xt_xmap);
  void (*get_destination_ranks)(Xt_xmap, int*);
  void (*get_source_ranks)(Xt_xmap, int*);
  Xt_xmap_iter (*get_out_iterator)(Xt_xmap);
  Xt_xmap_iter (*get_in_iterator)(Xt_xmap);
  Xt_xmap (*copy)(Xt_xmap);
  void (*delete)(Xt_xmap);
  int (*get_max_src_pos)(Xt_xmap);
  int (*get_max_dst_pos)(Xt_xmap);
  Xt_xmap (*reorder)(Xt_xmap xmap, enum xt_reorder_type type);
  Xt_xmap (*update_pos)(Xt_xmap xmap,
                        const int * src_positions, const int * dst_positions);
  Xt_xmap (*spread)(Xt_xmap xmap, int num_repetitions,
                    const int src_displacements[num_repetitions],
                    const int dst_displacements[num_repetitions]);
};

struct Xt_xmap_ {
  const struct Xt_xmap_vtable * vtable;
};

#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
