/**
 * @file xt_exchanger_neigh_alltoall.c
 *
 * @copyright Copyright  (C)  2018 Jörg Behrens <behrens@dkrz.de>
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
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_config_internal.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_request.h"
#include "xt/xt_request_msgs.h"
#include "xt_exchanger.h"
#include "xt_exchanger_neigh_alltoall.h"

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

static Xt_exchanger
xt_exchanger_neigh_alltoall_copy(Xt_exchanger exchanger,
                                 MPI_Comm newComm, int new_tag_offset);
static void xt_exchanger_neigh_alltoall_delete(Xt_exchanger exchanger);
static void xt_exchanger_neigh_alltoall_s_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data);
static void xt_exchanger_neigh_alltoall_a_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data,
                                                   Xt_request *request);
static int
xt_exchanger_neigh_alltoall_get_msg_ranks(Xt_exchanger exchanger,
                                          enum xt_msg_direction direction,
                                          int *restrict *ranks);

static MPI_Datatype
xt_exchanger_neigh_alltoall_get_MPI_Datatype(Xt_exchanger exchanger,
                                             int rank,
                                             enum xt_msg_direction direction,
                                             bool do_dup);

struct xt_exchanger_neigh_alltoall_team_share
{
  int *ranks;
  int *one_counts;
  MPI_Aint *displs;
  MPI_Comm nb_comm;
  int nmsg[2];
};

static void xt_exchanger_neigh_alltoall_team_share_default_init(void *share)
{
  struct xt_exchanger_neigh_alltoall_team_share *team_share = share;
  team_share->ranks = NULL;
  team_share->one_counts = NULL;
  team_share->displs = NULL;
  team_share->nb_comm = MPI_COMM_NULL;
}

static void xt_exchanger_neigh_alltoall_team_share_destroy(void *share)
{
  struct xt_exchanger_neigh_alltoall_team_share *team_share = share;
  if (team_share->nb_comm != MPI_COMM_NULL)
    xt_mpi_call(MPI_Comm_free(&team_share->nb_comm), Xt_default_comm);
  free(team_share->one_counts);
  free(team_share->displs);
  free(team_share->ranks);
}

const struct xt_exchanger_vtable xt_exchanger_neigh_alltoall_vtable = {
  .copy = xt_exchanger_neigh_alltoall_copy,
  .delete = xt_exchanger_neigh_alltoall_delete,
  .s_exchange = xt_exchanger_neigh_alltoall_s_exchange,
  .a_exchange = xt_exchanger_neigh_alltoall_a_exchange,
  .get_msg_ranks = xt_exchanger_neigh_alltoall_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_neigh_alltoall_get_MPI_Datatype,
  .team_share_size = sizeof (struct xt_exchanger_neigh_alltoall_team_share),
  .team_share_default_init
  = xt_exchanger_neigh_alltoall_team_share_default_init,
  .team_share_destroy = xt_exchanger_neigh_alltoall_team_share_destroy,
};

typedef struct Xt_exchanger_neigh_alltoall_ * Xt_exchanger_neigh_alltoall;

struct Xt_exchanger_neigh_alltoall_ {

  const struct xt_exchanger_vtable * vtable;

  MPI_Datatype * datatypes;
  struct xt_exchanger_neigh_alltoall_team_share *team_share,
    team_share_[];
};

static Xt_exchanger_neigh_alltoall
xt_exchanger_neigh_alltoall_alloc(size_t nsend, size_t nrecv,
                                  void *exchanger_team_share)
{
  size_t nmsg = nsend + nrecv;
  bool need_team_share_alloc = exchanger_team_share == NULL;
  Xt_exchanger_neigh_alltoall exchanger
    = xmalloc(sizeof(*exchanger)
              + (need_team_share_alloc ? sizeof (*exchanger->team_share) : 0));
  exchanger->datatypes = xmalloc(nmsg * sizeof(*exchanger->datatypes));
  exchanger->vtable = &xt_exchanger_neigh_alltoall_vtable;
  if (need_team_share_alloc) {
    exchanger->team_share = exchanger->team_share_;
    xt_exchanger_neigh_alltoall_team_share_default_init(exchanger->team_share_);
  } else
    exchanger->team_share = exchanger_team_share;
  return exchanger;
}

static void copy_dt(size_t n,
                    const struct Xt_redist_msg *restrict msgs,
                    MPI_Datatype *restrict datatypes,
                    MPI_Comm comm, bool dt_dup) {

  for (size_t i = 0; i < n; ++i)
    if (dt_dup)
      xt_mpi_call(MPI_Type_dup(msgs[i].datatype, datatypes + i), comm);
    else
      datatypes[i] = msgs[i].datatype;
}

static void copy_ranks(size_t n,
                       const struct Xt_redist_msg *restrict msgs,
                       int *restrict ranks)
{
  for (size_t i = 0; i < n; ++i)
    ranks[i] = msgs[i].rank;
}

Xt_exchanger
xt_exchanger_neigh_alltoall_new(int nsend, int nrecv,
                                const struct Xt_redist_msg *send_msgs,
                                const struct Xt_redist_msg *recv_msgs,
                                MPI_Comm comm, int tag_offset, Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */
  (void)tag_offset;
  int flag;
  xt_mpi_call(MPI_Comm_test_inter(comm, &flag), comm);
  if (flag)
    Xt_abort(comm, "ERROR(xt_exchanger_neigh_alltoall_new): "
             "inter-communicator's are not defined for virtual topologies",
             __FILE__, __LINE__);

  assert((nsend >= 0) & (nrecv >= 0));
  Xt_exchanger_neigh_alltoall exchanger
    = xt_exchanger_neigh_alltoall_alloc((size_t)nsend, (size_t)nrecv,
                                        config->exchanger_team_share);
  bool dt_dup = !(config->flags & exch_no_dt_dup);
  copy_dt((size_t)nsend, send_msgs, exchanger->datatypes, comm, dt_dup);
  copy_dt((size_t)nrecv, recv_msgs, exchanger->datatypes + nsend, comm, dt_dup);
  struct xt_exchanger_neigh_alltoall_team_share *team_share
    = exchanger->team_share;
  if (team_share->nb_comm == MPI_COMM_NULL) {
    size_t nmsg = (size_t)nsend + (size_t)nrecv;
    size_t max_msgs = MAX((size_t)nsend, (size_t)nrecv);
    team_share->nmsg[RECV] = nrecv;
    team_share->nmsg[SEND] = nsend;
    team_share->ranks = xmalloc(nmsg * sizeof(*team_share->ranks));
    team_share->one_counts
      = xmalloc(max_msgs * sizeof(*team_share->one_counts));
    team_share->displs = xmalloc(max_msgs * sizeof(*team_share->displs));
    copy_ranks((size_t)nsend, send_msgs, team_share->ranks);
    copy_ranks((size_t)nrecv, recv_msgs, team_share->ranks+nsend);
    for (size_t i = 0; i < max_msgs; ++i) {
      team_share->one_counts[i] = 1;
      team_share->displs[i] = 0;
    }
    enum { no_reorder = 0 }; // no reordering of ranks in new comm
#if __GNUC__ >= 11 && __GNUC__ <= 13
    /* GCC 11 has no means to specify that the special value pointer
     * MPI_UNWEIGHTED does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
    xt_mpi_call(
      MPI_Dist_graph_create_adjacent(
        comm, nrecv, team_share->ranks + nsend, MPI_UNWEIGHTED, nsend,
        team_share->ranks, MPI_UNWEIGHTED, MPI_INFO_NULL, no_reorder,
        &team_share->nb_comm), comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif
  }
  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_neigh_alltoall_copy(Xt_exchanger exchanger,
                                 MPI_Comm new_comm, int new_tag_offset)
{
  (void)new_tag_offset;
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  const struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;
  size_t nsend = (size_t)team_share->nmsg[SEND],
    nrecv = (size_t)team_share->nmsg[RECV];

  Xt_exchanger_neigh_alltoall
    exchanger_copy = xt_exchanger_neigh_alltoall_alloc(nsend, nrecv, NULL);

  struct xt_exchanger_neigh_alltoall_team_share
    *team_share_copy = exchanger_copy->team_share;
  size_t nmsg = nsend + nrecv;
  size_t max_msgs = MAX(nsend, nrecv);
  team_share_copy->nmsg[SEND] = (int)nsend;
  team_share_copy->nmsg[RECV] = (int)nrecv;
  team_share_copy->ranks = xmalloc(nmsg * sizeof(*team_share_copy->ranks));
  team_share_copy->one_counts
    = xmalloc(max_msgs * sizeof(*team_share_copy->one_counts));
  team_share_copy->displs
    = xmalloc(max_msgs * sizeof(*team_share_copy->displs));
  memcpy(team_share_copy->ranks, team_share->ranks,
         nmsg * sizeof(*team_share_copy->ranks));
  for (size_t i = 0; i < max_msgs; ++i) {
    team_share_copy->one_counts[i] = 1;
    team_share_copy->displs[i] = 0;
  }
#ifndef NEED_MPICH_UNWEIGHTED_COMM_DUP_WORKAROUND
  xt_mpi_call(MPI_Comm_dup(team_share->nb_comm, &team_share_copy->nb_comm),
              new_comm);
#else
  /* MPICH up to version 3.4.2 at least cannot do MPI_Comm_dup
   * for topology communicators */
  enum { no_reorder = 0 }; // no reordering of ranks in new comm
#if __GNUC__ >= 11 && __GNUC__ <= 13
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_UNWEIGHTED does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  xt_mpi_call(
    MPI_Dist_graph_create_adjacent(
      new_comm, (int)nrecv, team_share_copy->ranks + nsend, MPI_UNWEIGHTED,
      (int)nsend, team_share_copy->ranks, MPI_UNWEIGHTED, MPI_INFO_NULL,
      no_reorder, &team_share_copy->nb_comm), new_comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif
#endif
  for (size_t i = 0; i < nmsg; ++i)
    xt_mpi_call(MPI_Type_dup(exchanger_na->datatypes[i],
                             exchanger_copy->datatypes + i), new_comm);

  return (Xt_exchanger)exchanger_copy;
}

static void xt_exchanger_neigh_alltoall_delete(Xt_exchanger exchanger) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share_na = exchanger_na->team_share;

  size_t nmsg = (size_t)team_share_na->nmsg[SEND]
    + (size_t)team_share_na->nmsg[RECV];

  for (size_t i = 0; i < nmsg; ++i)
    xt_mpi_call(MPI_Type_free(exchanger_na->datatypes + i),
                team_share_na->nb_comm);
  free(exchanger_na->datatypes);
  if (exchanger_na->team_share == exchanger_na->team_share_)
    xt_exchanger_neigh_alltoall_team_share_destroy(exchanger_na->team_share_);
  free(exchanger_na);
}

static void xt_exchanger_neigh_alltoall_s_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;

  xt_mpi_call(
    MPI_Neighbor_alltoallw(src_data, team_share->one_counts,
                           team_share->displs, exchanger_na->datatypes,
                           dst_data, team_share->one_counts,
                           team_share->displs, exchanger_na->datatypes +
                           (size_t)team_share->nmsg[SEND],
                           team_share->nb_comm),
    team_share->nb_comm);
}

static void xt_exchanger_neigh_alltoall_a_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data,
                                                   Xt_request *request) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;

#if defined OMPI_MAJOR_VERSION && OMPI_MAJOR_VERSION == 4 \
  && OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION == 2
  /* ugly work-around: Open MPI retains pointers to arguments
   * of MPI_Ineighbor_alltoallw, but the exchanger might have
   * been destroyed by the time Open MPI starts to use it.
   * Therefore, this work-around resizes the request object and
   * stores copies of some arguments in the request object tail.
   * Not pretty, I know, but Open MPI 4.0.x is too recent to ignore. */
  struct Xt_request_msgs_ {
    const struct Xt_request_vtable *vtable;
    int n;
    MPI_Comm comm;
    MPI_Request requests[];
  };
  size_t nmsg = (size_t)team_share->nmsg[SEND] + (size_t)team_share->nmsg[RECV],
    header_size = sizeof (struct Xt_request_msgs_),
    body_size = sizeof (MPI_Datatype) * nmsg + sizeof (int) * nmsg,
    dt_size = nmsg * sizeof (MPI_Datatype);
  body_size = (body_size + sizeof (MPI_Datatype) - 1)
    / sizeof (MPI_Datatype) * sizeof (MPI_Datatype);
  struct Xt_request_msgs_ *request_
    = xrealloc(xt_request_msgs_new(1, (MPI_Request[]){ MPI_REQUEST_NULL},
                                   team_share->nb_comm),
               header_size + body_size + dt_size);
  *request = (Xt_request)request_;
  MPI_Datatype *dt_copy
    = (MPI_Datatype *)(void *)((unsigned char *)request_ + header_size + body_size);
  memcpy(dt_copy, exchanger_na->datatypes, dt_size);
  xt_mpi_call(
    MPI_Ineighbor_alltoallw(src_data, team_share->one_counts,
                            team_share->displs, dt_copy,
                            dst_data, team_share->one_counts,
                            team_share->displs, dt_copy +
                            (size_t)team_share->nmsg[SEND],
                            team_share->nb_comm, request_->requests),
    team_share->nb_comm);
#else
  MPI_Request tmp_request;
  xt_mpi_call(
    MPI_Ineighbor_alltoallw(src_data, team_share->one_counts,
                            team_share->displs, exchanger_na->datatypes,
                            dst_data, team_share->one_counts,
                            team_share->displs, exchanger_na->datatypes +
                            (size_t)team_share->nmsg[SEND],
                            team_share->nb_comm, &tmp_request),
    team_share->nb_comm);

  *request = xt_request_msgs_new(1, &tmp_request, team_share->nb_comm);
#endif
}

static MPI_Datatype
xt_exchanger_neigh_alltoall_get_MPI_Datatype(Xt_exchanger exchanger,
                                             int rank,
                                             enum xt_msg_direction direction,
                                             bool do_dup)
{
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;
  size_t nsend = (size_t)team_share->nmsg[SEND],
    nmsg = (size_t)team_share->nmsg[direction],
    ofs = direction == SEND ? 0 : nsend;
  int *restrict ranks = team_share->ranks + ofs;
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < nmsg; ++i) {
    if (ranks[i] == rank) {
      if (do_dup)
        xt_mpi_call(MPI_Type_dup(exchanger_na->datatypes[i+ofs], &datatype_copy),
                    team_share->nb_comm);
      else
        datatype_copy = exchanger_na->datatypes[i+ofs];
      break;
    }
  }
  return datatype_copy;
}

static int
xt_exchanger_neigh_alltoall_get_msg_ranks(Xt_exchanger exchanger,
                                          enum xt_msg_direction direction,
                                          int *restrict *ranks)
{
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;
  size_t nsend = (size_t)team_share->nmsg[SEND],
    nmsg = (size_t)team_share->nmsg[direction],
    ofs = direction == SEND ? 0 : nsend;
  int *ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc(nmsg * sizeof(**ranks));
  memcpy(ranks_, team_share->ranks + ofs, nmsg * sizeof(**ranks));
  return (int)nmsg;
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
