/**
 * @file xt_exchanger_vtable.c
 *
 * @copyright Copyright  (C)  2021 Jörg Behrens <behrens@dkrz.de>
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
#  include <config.h>
#endif

#include "core/core.h"

#include "xt_exchanger.h"
#include "xt_exchanger_mix_isend_irecv.h"
#include "xt_exchanger_simple_base.h"
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
#  include "xt_exchanger_neigh_alltoall.h"
#endif
#include "xt_exchanger_irecv_isend_packed.h"
#include "xt_exchanger_irecv_isend_ddt_packed.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_irecv_send.h"

const struct xt_exchanger_vtable *
xt_exchanger_new_get_vtable(Xt_exchanger_new exchanger_new)
{
  const struct xt_exchanger_vtable *vtab = NULL;
  if (exchanger_new == xt_exchanger_mix_isend_irecv_new)
    vtab = &xt_exchanger_mix_isend_irecv_vtable;
  else if (exchanger_new == xt_exchanger_irecv_isend_new
           || exchanger_new == xt_exchanger_irecv_isend_packed_new
           || exchanger_new == xt_exchanger_irecv_isend_ddt_packed_new
           || exchanger_new == xt_exchanger_irecv_send_new)
    vtab = &xt_exchanger_simple_base_vtable;
#if XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
  else if (exchanger_new == xt_exchanger_neigh_alltoall_new)
    vtab = &xt_exchanger_neigh_alltoall_vtable;
#endif
  else
    Xt_abort(Xt_default_comm, "unexpected exchanger constructor",
             "xt_exchanger_vtable.c", __LINE__);
  return vtab;
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
