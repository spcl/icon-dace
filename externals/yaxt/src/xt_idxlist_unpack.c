/**
 * @file xt_idxlist_unpack.c
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

#include "core/core.h"
#include "xt/xt_mpi.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxlist_collection.h"
#include "xt/xt_idxsection.h"
#include "xt/xt_idxstripes.h"
#include "xt_idxlist_unpack.h"

typedef Xt_idxlist (*idxlist_unpack)(void*,int,int*,MPI_Comm);

static const idxlist_unpack unpack[] = {
  xt_idxempty_unpack,
  xt_idxvec_unpack,
  xt_idxlist_collection_unpack,
  xt_idxsection_unpack,
  xt_idxstripes_unpack,
};

static const unsigned xt_num_unpack_routines
= sizeof(unpack) / sizeof(unpack[0]);

Xt_idxlist xt_idxlist_unpack(void *buffer, int buffer_size,
                             int *position, MPI_Comm comm) {

  int type;
  // unpack the type of the index list to be unpacked
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &type, 1, MPI_INT, comm), comm);
  if (type >= 0 && (unsigned)type < xt_num_unpack_routines)
    return unpack[type](buffer, buffer_size, position, comm);
  Xt_abort(comm, "xt_idxlist_unpack: unknown index list type",
           __FILE__, __LINE__);

  return NULL;
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
