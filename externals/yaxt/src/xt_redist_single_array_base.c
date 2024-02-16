/**
 * @file xt_redist_single_array_base.c
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>

#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_single_array_base.h"
#include "xt_redist_internal.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_request.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt_exchanger.h"
#include "xt_config_internal.h"

static Xt_redist
redist_sab_copy(Xt_redist redist);

static void
redist_sab_delete(Xt_redist redist);

static void
redist_sab_s_exchange(Xt_redist redist, int num_arrays,
                      const void **src_data, void **dst_data);

static void
redist_sab_a_exchange(Xt_redist redist, int num_arrays,
                      const void **src_data, void **dst_data,
                      Xt_request *request);

static void
redist_sab_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data);

static void
redist_sab_a_exchange1(Xt_redist redist, const void *src_data, void *dst_data,
                       Xt_request *request);

static int redist_sab_get_num_msg(Xt_redist redist,
                                  enum xt_msg_direction direction);

static MPI_Datatype
redist_sab_get_MPI_Datatype(Xt_redist redist, int rank,
                            enum xt_msg_direction direction,
                            bool do_dup);

static int
redist_sab_get_msg_ranks(Xt_redist redist,
                         enum xt_msg_direction direction,
                         int *restrict *ranks);

static MPI_Comm
redist_sab_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_sab_vtable = {
  .copy                  = redist_sab_copy,
  .delete                = redist_sab_delete,
  .s_exchange            = redist_sab_s_exchange,
  .a_exchange            = redist_sab_a_exchange,
  .s_exchange1           = redist_sab_s_exchange1,
  .a_exchange1           = redist_sab_a_exchange1,
  .get_num_msg           = redist_sab_get_num_msg,
  .get_msg_MPI_Datatype  = redist_sab_get_MPI_Datatype,
  .get_msg_ranks         = redist_sab_get_msg_ranks,
  .get_MPI_Comm          = redist_sab_get_MPI_Comm
};

typedef struct Xt_redist_sab_ *Xt_redist_sab;

struct Xt_redist_sab_ {

  const struct xt_redist_vtable *vtable;

  Xt_exchanger exchanger;

  int nmsg[2];

  MPI_Comm comm;
  int tag_offset;
};

Xt_redist xt_redist_single_array_base_new(int nsend, int nrecv,
                                          const struct Xt_redist_msg *send_msgs,
                                          const struct Xt_redist_msg *recv_msgs,
                                          MPI_Comm comm)
{
  return xt_redist_single_array_base_custom_new(nsend, nrecv, send_msgs,
                                                recv_msgs, comm,
                                                (Xt_config)&xt_default_config);
}

Xt_redist xt_redist_single_array_base_custom_new(
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs,
  const struct Xt_redist_msg *recv_msgs,
  MPI_Comm comm, Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized());

  Xt_redist_sab redist = xmalloc(sizeof (*redist));

  redist->comm = xt_mpi_comm_smart_dup(comm, &redist->tag_offset);
  Xt_exchanger_new exchanger_new
    = xt_config_get_exchange_new_by_comm(config, redist->comm);
  redist->exchanger
    = exchanger_new(nsend, nrecv, send_msgs, recv_msgs,
                    redist->comm, redist->tag_offset, config);

  redist->vtable = &redist_sab_vtable;
  redist->nmsg[SEND] = nsend;
  redist->nmsg[RECV] = nrecv;
  return (Xt_redist)redist;
}

static inline Xt_redist_sab
xrsab(void *redist)
{
  return (Xt_redist_sab)redist;
}

static Xt_redist
redist_sab_copy(Xt_redist redist)
{
  Xt_redist_sab redist_sab = xrsab(redist);
  Xt_redist_sab redist_sab_new = xmalloc(sizeof *redist_sab_new);
  redist_sab_new->vtable = redist_sab->vtable;
  for (size_t i = 0; i < 2; ++i)
    redist_sab_new->nmsg[i] = redist_sab->nmsg[i];
  redist_sab_new->comm = xt_mpi_comm_smart_dup(redist_sab->comm,
                                               &redist_sab_new->tag_offset);
  redist_sab_new->exchanger
    = xt_exchanger_copy(redist_sab->exchanger, redist_sab_new->comm,
                        redist_sab_new->tag_offset);
  return (Xt_redist)redist_sab_new;
}

static void
redist_sab_delete(Xt_redist redist) {

  Xt_redist_sab redist_sab = xrsab(redist);

  xt_exchanger_delete(redist_sab->exchanger);

  xt_mpi_comm_smart_dedup(&redist_sab->comm, redist_sab->tag_offset);

  free(redist_sab);
}

static void
redist_sab_s_exchange(Xt_redist redist, int num_arrays,
                      const void **src_data, void **dst_data)
{
  Xt_redist_sab redist_rep = xrsab(redist);
  if (num_arrays == 1)
    redist_sab_s_exchange1(redist, src_data[0], dst_data[0]);
  else
    Xt_abort(redist_rep->comm, "ERROR: multi-array s_exchange is not"
             " implemented for this xt_redist type "
             "(Xt_redist_single_array_base)", __FILE__, __LINE__);
}

static void
redist_sab_a_exchange(Xt_redist redist, int num_arrays,
                      const void **src_data, void **dst_data,
                      Xt_request *request)
{
  Xt_redist_sab redist_rep = xrsab(redist);
  if (num_arrays == 1)
    redist_sab_a_exchange1(redist, src_data[0], dst_data[0], request);
  else
    Xt_abort(redist_rep->comm, "ERROR: multi-array a_exchange is not"
             " implemented for this xt_redist type "
             "(Xt_redist_single_array_base)", __FILE__, __LINE__);
}

static void
redist_sab_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data) {

  Xt_redist_sab redist_sab = xrsab(redist);

  xt_exchanger_s_exchange(redist_sab->exchanger, src_data, dst_data);
}

static void
redist_sab_a_exchange1(Xt_redist redist, const void *src_data, void *dst_data,
                       Xt_request *request) {

  Xt_redist_sab redist_sab = xrsab(redist);

  xt_exchanger_a_exchange(redist_sab->exchanger, src_data, dst_data, request);
}

static int redist_sab_get_num_msg(Xt_redist redist,
                                  enum xt_msg_direction direction)
{
  return xrsab(redist)->nmsg[direction];
}

static MPI_Datatype
redist_sab_get_MPI_Datatype(Xt_redist redist, int rank,
                            enum xt_msg_direction direction,
                            bool do_dup)
{
  return xt_exchanger_get_MPI_Datatype(xrsab(redist)->exchanger, rank,
                                       direction, do_dup);
}

static int
redist_sab_get_msg_ranks(Xt_redist redist,
                         enum xt_msg_direction direction,
                         int *restrict *ranks)
{
  return xt_exchanger_get_msg_ranks(xrsab(redist)->exchanger, direction, ranks);
}

static MPI_Comm
redist_sab_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_sab redist_sab = xrsab(redist);

  return redist_sab->comm;
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
