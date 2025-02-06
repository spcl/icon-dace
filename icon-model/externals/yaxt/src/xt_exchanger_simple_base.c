/**
 * @file xt_exchanger_simple_base.c
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

#include <assert.h>
#include <mpi.h>
#include <stdbool.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_simple_base.h"

static Xt_exchanger
xt_exchanger_simple_base_copy(Xt_exchanger exchanger,
                              MPI_Comm newComm, int new_tag_offset);
static void xt_exchanger_simple_base_delete(Xt_exchanger exchanger);
static void xt_exchanger_simple_base_s_exchange(Xt_exchanger exchanger,
                                                const void * src_data,
                                                void * dst_data);
static void xt_exchanger_simple_base_a_exchange(Xt_exchanger exchanger,
                                                const void * src_data,
                                                void * dst_data,
                                                Xt_request *request);
static int
xt_exchanger_simple_base_get_msg_ranks(Xt_exchanger exchanger,
                                       enum xt_msg_direction direction,
                                       int *restrict *ranks);

static MPI_Datatype
xt_exchanger_simple_base_get_MPI_Datatype(Xt_exchanger exchanger,
                                          int rank,
                                          enum xt_msg_direction direction,
                                          bool do_dup);

static Xt_exchanger_omp_share
xt_exchanger_simple_base_create_omp_share(Xt_exchanger exchanger);

const struct xt_exchanger_vtable xt_exchanger_simple_base_vtable = {
  .copy = xt_exchanger_simple_base_copy,
  .delete = xt_exchanger_simple_base_delete,
  .s_exchange = xt_exchanger_simple_base_s_exchange,
  .a_exchange = xt_exchanger_simple_base_a_exchange,
  .get_msg_ranks = xt_exchanger_simple_base_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_simple_base_get_MPI_Datatype,
  .create_omp_share = xt_exchanger_simple_base_create_omp_share,
};

typedef struct Xt_exchanger_simple_base_ * Xt_exchanger_simple_base;

struct Xt_exchanger_simple_base_ {

  const struct xt_exchanger_vtable * vtable;

  int nmsg[2];
  int config_flags;
  int tag_offset;
  MPI_Comm comm;
  xt_simple_s_exchange_func s_func;
  xt_simple_a_exchange_func a_func;
  xt_simple_create_omp_share_func create_omp_share_func;
  struct Xt_redist_msg msgs[];
};

static Xt_exchanger_simple_base
xt_exchanger_simple_base_alloc(size_t nmsg)
{
  Xt_exchanger_simple_base exchanger;
  size_t header_size = sizeof(*exchanger),
    body_size = nmsg * sizeof (exchanger->msgs[0]);
  exchanger = xmalloc(header_size + body_size);
  exchanger->vtable = &xt_exchanger_simple_base_vtable;
  return exchanger;
}

static inline int
adjusted_rank(int r, int comm_rank, int comm_size)
{
  return r + (r <= comm_rank ? comm_size : 0);
}

#define XT_SORTFUNC_DECL static
#define SORT_TYPE struct Xt_redist_msg
#define SORT_TYPE_SUFFIX redist_msg
#define SORT_TYPE_CMP_LT(u,v,i,j)                       \
  (adjusted_rank((u).rank, comm_rank, comm_size)        \
   < adjusted_rank((v).rank, comm_rank, comm_size))
#define SORT_TYPE_CMP_LE(u,v,i,j)                       \
  (adjusted_rank((u).rank, comm_rank, comm_size)        \
   <= adjusted_rank((v).rank, comm_rank, comm_size))
#define SORT_TYPE_CMP_EQ(u,v,i,j) (((u).rank) == ((v).rank))
#define XT_SORT_EXTRA_ARGS_DECL , int comm_rank, int comm_size
#define XT_SORT_EXTRA_ARGS_PASS , comm_rank, comm_size
#define XT_SORT_VECSWAP_EXTRA_ARGS_DECL
#define XT_SORT_VECSWAP_EXTRA_ARGS_PASS

#include "xt_quicksort_base.h"

Xt_exchanger
xt_exchanger_simple_base_new(
  int nsend, int nrecv,
  const struct Xt_redist_msg *send_msgs,
  const struct Xt_redist_msg *recv_msgs,
  MPI_Comm comm, int tag_offset,
  xt_simple_s_exchange_func s_func,
  xt_simple_a_exchange_func a_func,
  xt_simple_create_omp_share_func create_omp_share_func,
  Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */

  if (s_func == 0)
    Xt_abort(comm, "ERROR(xt_exchanger_simple_base_new): invalid synchronous "
             "exchange function pointer", __FILE__, __LINE__);

  assert((nsend >= 0) & (nrecv >= 0));
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  Xt_exchanger_simple_base exchanger
    = xt_exchanger_simple_base_alloc(nmsg);
  exchanger->comm = comm;
  exchanger->tag_offset = tag_offset;
  exchanger->nmsg[SEND] = nsend;
  exchanger->nmsg[RECV] = nrecv;
  bool dt_dup = !(config->flags & exch_no_dt_dup);
  xt_redist_msgs_strided_copy((size_t)nsend, send_msgs, sizeof (send_msgs[0]),
                              exchanger->msgs, sizeof (exchanger->msgs[0]),
                              comm, dt_dup);
  xt_redist_msgs_strided_copy((size_t)nrecv, recv_msgs, sizeof (recv_msgs[0]),
                              exchanger->msgs + nsend,
                              sizeof (exchanger->msgs[0]),
                              comm, dt_dup);
  exchanger->s_func = s_func;
  exchanger->a_func = a_func;
  exchanger->create_omp_share_func = create_omp_share_func;

  {
    int comm_size, comm_rank, is_inter;
    xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
    xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
    int (*get_comm_size)(MPI_Comm, int *)
      = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
    xt_mpi_call(get_comm_size(comm, &comm_size), comm);
    /* In order to avoid congestion of messages, the order of send and
     * receive messages is changed. This is done by sorting the
     * messages according to the rank of the respective message
     * partner. Before the sorting to ranks that are smaller or equal
     * to the local rank the size of the communicator is added.
     *
     * example: process 5 is supposed to communicate with
     * processes: 9, 5, 2, 6, 1
     * 1. add comm_size(10): 9, 15, 12, 6, 11
     * 2. sort: 6, 9, 11, 12, 15
     * 3. subtrace comm_size again -> final order: 6, 9, 1, 2, 5
     */
    xt_quicksort_redist_msg(exchanger->msgs, (size_t)nsend,
                            comm_rank, comm_size);
    xt_quicksort_redist_msg(exchanger->msgs + (size_t)nsend, (size_t)nrecv,
                            comm_rank, comm_size);
  }

  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_simple_base_copy(Xt_exchanger exchanger,
                              MPI_Comm new_comm, int new_tag_offset)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;
  int nsend = exchanger_sb->nmsg[SEND],
    nrecv = exchanger_sb->nmsg[RECV];
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  Xt_exchanger_simple_base exchanger_copy
    = xt_exchanger_simple_base_alloc(nmsg);
  exchanger_copy->nmsg[SEND] = nsend;
  exchanger_copy->nmsg[RECV] = nrecv;
  exchanger_copy->s_func = exchanger_sb->s_func;
  exchanger_copy->a_func = exchanger_sb->a_func;
  exchanger_copy->create_omp_share_func = exchanger_sb->create_omp_share_func;
  struct Xt_redist_msg *restrict new_msgs = exchanger_copy->msgs,
    *restrict orig_msgs = exchanger_sb->msgs;
  xt_redist_msgs_strided_copy(nmsg, orig_msgs, sizeof (*orig_msgs),
                              new_msgs, sizeof (*new_msgs),
                              new_comm, true);
  exchanger_copy->comm = new_comm;
  exchanger_copy->tag_offset = new_tag_offset;
  return (Xt_exchanger)exchanger_copy;
}


static void xt_exchanger_simple_base_delete(Xt_exchanger exchanger) {

  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  size_t nmsg = (size_t)exchanger_sb->nmsg[SEND] + (size_t)exchanger_sb->nmsg[RECV];
  struct Xt_redist_msg *restrict msgs = exchanger_sb->msgs;
  xt_redist_msgs_strided_destruct(nmsg, msgs, exchanger_sb->comm,
                                  sizeof (msgs[0]));
  free(exchanger_sb);
}

static void xt_exchanger_simple_base_s_exchange(Xt_exchanger exchanger,
                                                const void * src_data,
                                                void * dst_data) {

  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  int nsend = exchanger_sb->nmsg[SEND];
  exchanger_sb->s_func(src_data, dst_data, nsend,
                       exchanger_sb->nmsg[RECV], exchanger_sb->msgs,
                       exchanger_sb->msgs + (size_t)nsend,
                       exchanger_sb->tag_offset, exchanger_sb->comm);
}

static void xt_exchanger_simple_base_a_exchange(Xt_exchanger exchanger,
                                                const void * src_data,
                                                void * dst_data,
                                                Xt_request *request) {

  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  if (exchanger_sb->a_func == NULL)
    Xt_abort(exchanger_sb->comm, "ERROR(xt_exchanger_simple_base_a_exchange): "
             "asynchronous exchange function not defined for current exchanger",
             __FILE__, __LINE__);

  int nsend = exchanger_sb->nmsg[SEND];
  exchanger_sb->a_func(src_data, dst_data, nsend,
                       exchanger_sb->nmsg[RECV], exchanger_sb->msgs,
                       exchanger_sb->msgs + (size_t)nsend,
                       exchanger_sb->tag_offset, exchanger_sb->comm, request);
}

static MPI_Datatype
xt_exchanger_simple_base_get_MPI_Datatype(Xt_exchanger exchanger,
                                          int rank,
                                          enum xt_msg_direction direction,
                                          bool do_dup)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;
  size_t nsend = (size_t)exchanger_sb->nmsg[SEND],
    nmsg = (size_t)exchanger_sb->nmsg[direction],
    ofs = direction == SEND ? 0 : nsend;
  struct Xt_redist_msg *restrict msgs = exchanger_sb->msgs + ofs;
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < nmsg; ++i)
    if (msgs[i].rank == rank) {
      if (do_dup)
        xt_mpi_call(MPI_Type_dup(msgs[i].datatype, &datatype_copy),
                    exchanger_sb->comm);
      else
        datatype_copy = msgs[i].datatype;
      break;
    }
  return datatype_copy;
}

static int
xt_exchanger_simple_base_get_msg_ranks(Xt_exchanger exchanger,
                                       enum xt_msg_direction direction,
                                       int *restrict *ranks)
{
  Xt_exchanger_simple_base exchanger_sb = (Xt_exchanger_simple_base)exchanger;
  size_t nmsg = (size_t)exchanger_sb->nmsg[direction];
  struct Xt_redist_msg *restrict msgs = exchanger_sb->msgs
    + (direction == RECV ? (size_t)exchanger_sb->nmsg[SEND] : 0);
  int *restrict ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc(nmsg * sizeof (*ranks_));
  for (size_t i = 0; i < nmsg; ++i)
    ranks_[i] = msgs[i].rank;
  return (int)nmsg;
}

static Xt_exchanger_omp_share
xt_exchanger_simple_base_create_omp_share(Xt_exchanger exchanger)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;
  int nsend = exchanger_sb->nmsg[SEND];
  return exchanger_sb->create_omp_share_func(
    nsend, exchanger_sb->nmsg[RECV],
    exchanger_sb->msgs, exchanger_sb->msgs + nsend, exchanger_sb->comm);
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
