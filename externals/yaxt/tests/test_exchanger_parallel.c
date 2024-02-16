/**
 * @file test_exchanger_parallel.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#include "../src/xt_exchanger.h"
#include "../src/xt_exchanger_irecv_send.h"
#include "../src/xt_exchanger_irecv_isend.h"
#include "../src/xt_exchanger_mix_isend_irecv.h"
#include "../src/xt_exchanger_irecv_isend_packed.h"
#include "../src/xt_exchanger_neigh_alltoall.h"
#include "../src/xt_redist_internal.h"
#include "../src/xt_mpi_internal.h"
#include "../src/xt_config_internal.h"
#include "core/ppm_xfuncs.h"

#include "test_redist_common.h"
#include "core/ppm_xfuncs.h"
#include "tests.h"

struct test_message {

  int rank;       // rank of communication partner
  const int *pos;   // positions to be sent/received
  int num_pos; // number of positions
};

static Xt_exchanger_new *parse_options(int *argc, char ***argv);

static void
test_bcast(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config);
static void
test_gather(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config);
static void
test_all2all(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config);
static void
test_rr(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config);
static void
test_intercomm_all2all(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config);

static int test_freq = 3;

int main(int argc, char **argv)
{

  // init mpi
  xt_mpi_call(MPI_Init(&argc, &argv), MPI_COMM_WORLD);

  xt_mpi_init();
  xt_config_defaults_init();

  Xt_exchanger_new *exchangers_new = parse_options(&argc, &argv);
  Xt_config config = xt_config_new();

  for (size_t i = 0; exchangers_new[i] != (Xt_exchanger_new)0; ++i) {
    Xt_exchanger_new exchanger_new = exchangers_new[i];

    test_bcast(MPI_COMM_WORLD, exchanger_new, config);

    test_gather(MPI_COMM_WORLD, exchanger_new, config);

    test_all2all(MPI_COMM_WORLD, exchanger_new, config);

    test_rr(MPI_COMM_WORLD, exchanger_new, config);

    test_intercomm_all2all(MPI_COMM_WORLD, exchanger_new, config);
  }
  xt_config_delete(config);
  free(exchangers_new);

  xt_mpi_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static Xt_exchanger_new *parse_options(int *argc, char ***argv)
{
  static Xt_exchanger_new exchanger_table[] = {
    [xt_exchanger_irecv_send] = xt_exchanger_irecv_send_new,
    [xt_exchanger_irecv_isend] = xt_exchanger_irecv_isend_new,
    [xt_exchanger_irecv_isend_packed] = xt_exchanger_irecv_isend_packed_new,
    [xt_exchanger_mix_isend_irecv] = xt_exchanger_mix_isend_irecv_new,
#if MPI_VERSION >= 3
    [xt_exchanger_neigh_alltoall] = xt_exchanger_neigh_alltoall_new,
#endif
  };
  enum {
    num_exchanger = sizeof (exchanger_table) / sizeof (exchanger_table[0]),
  };
  Xt_exchanger_new *exchangers_new = xmalloc(2 * sizeof (*exchangers_new));
  exchangers_new[0] = xt_exchanger_mix_isend_irecv_new;
  exchangers_new[1] = (Xt_exchanger_new)0;
  size_t cur_ex = 0;
  int opt;
  while ((opt = getopt(*argc, *argv, "m:s:")) != -1) {
    switch (opt) {
    case 'm':
      {
        int exchanger_new_id = xt_exchanger_id_by_name(optarg);
        if (exchanger_new_id == -1)
        {
          fprintf(stderr, "Unknown exchanger constructor requested: %s\n",
                  optarg);
          exit(EXIT_FAILURE);
        }
#if MPI_VERSION < 3
        else if (exchanger_new_id == xt_exchanger_neigh_alltoall)
        {
          fputs("xt_exchanger_neigh_alltoall_new requires MPI version 3.0 or "
                "higher\n", stderr);
          continue;
        }
#endif
        exchangers_new[cur_ex] = exchanger_table[exchanger_new_id];
        ++cur_ex;
        exchangers_new = xrealloc(exchangers_new,
                                  sizeof (*exchangers_new) * (cur_ex + 1));
        exchangers_new[cur_ex] = (Xt_exchanger_new)0;
      }
      break;
    case 's':
      {
        char *endptr;
        errno = 0;
        long v = strtol(optarg, &endptr, 0);
        if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
            || (errno != 0 && v == 0)) {
          perror("failed to parse argument to -s option");
          exit(EXIT_FAILURE);
        }
        if (endptr == optarg) {
          fputs("malformed argument to -s option, no digits were found\n",
                stderr);
          exit(EXIT_FAILURE);
        }
        if (v < 1 || v > INT_MAX) {
          fprintf(stderr, "value of -s option (%ld) out of range [1,%d]\n",
                  v, INT_MAX);
          exit(EXIT_FAILURE);
        }
        test_freq = (int)v;
      }
      break;
    }
  }
  return exchangers_new;
}

static void
test_bcast(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // bcast pattern
  int incr = comm_size/test_freq + (comm_size < test_freq);
  for (int i = 0; i < comm_size; i += incr) {

    // setup

    struct Xt_redist_msg *send_msgs = NULL;
    int nsend = 0;
    int nrecv = my_rank != i;

    if (my_rank == i) {
      nsend = comm_size - 1;
      send_msgs = xmalloc((size_t)nsend * sizeof (*send_msgs));
      for (size_t j = 0; j < (size_t)i; ++j)
        send_msgs[j] = (struct Xt_redist_msg){.rank=(int)j,
                                              .datatype=MPI_INT};
      for (size_t j = (size_t)i; j < (size_t)nsend; ++j)
        send_msgs[j] = (struct Xt_redist_msg){.rank=(1+(int)j)%comm_size,
                                              .datatype=MPI_INT};
    }
    struct Xt_redist_msg recv_msgs[2] =
      {{.rank=-1, .datatype=MPI_DATATYPE_NULL},
       {.rank=i, .datatype=MPI_INT}};

    Xt_exchanger exchanger = exchanger_new(nsend, nrecv,
                                           send_msgs,
                                           recv_msgs + (my_rank != i),
                                           comm, 0, config);

    // test
    int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
    for (int async = 0; async < 1 + test_async; ++async) {

      int src_data[1] = { my_rank == i ? 4711 : -1 };
      int dst_data[1] = { my_rank == i ? 4711 : -1 };

      if (async) {
        Xt_request request;
        int flag;
        xt_exchanger_a_exchange(exchanger, (void*)(src_data),
                                (void*)(dst_data), &request);
        xt_request_test(&request, &flag);
        xt_request_wait(&request);
        xt_request_test(&request, &flag);
        if (!flag) PUT_ERR("invalid flag result\n");
      } else {
        xt_exchanger_s_exchange(exchanger, (void*)(src_data),
                                (void*)(dst_data));
      }

      if (dst_data[0] != 4711) PUT_ERR("invalid data\n");
    }

    // cleanup
    free(send_msgs);
    xt_exchanger_delete(exchanger);
  }
}

static void
test_gather(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // gather pattern
  // prepare datatypes outside of loop for load-balance
  MPI_Datatype *dt_by_ofs = xmalloc((size_t)comm_size * sizeof (*dt_by_ofs));
  dt_by_ofs[0] = MPI_INT;
  for (size_t j = 1; j < (size_t)comm_size; ++j)
  {
    MPI_Type_indexed(1, (int[]){1}, (int[]){(int)j}, MPI_INT, dt_by_ofs + j);
    MPI_Type_commit(dt_by_ofs + j);
  }
  int *dst_data = xmalloc(((size_t)comm_size - 1) * sizeof (*dst_data) * 2);
  int incr = comm_size/test_freq + (comm_size < test_freq);
  for (int i = 0; i < comm_size; i += incr) {

    // setup
    int nsend = i != my_rank;
    size_t nrecv = 0;


    struct Xt_redist_msg send_msgs[1] = {{.rank=i, .datatype=MPI_INT}};
    struct Xt_redist_msg *recv_msgs = NULL;
    if (my_rank == i) {
      nrecv = (size_t)comm_size - 1;
      recv_msgs = xmalloc(nrecv * sizeof (*recv_msgs) * 2);
      for (size_t j = 0; j < nrecv; ++j) {
        recv_msgs[j].rank = (i + (int)j + 1)%comm_size;
        recv_msgs[j].datatype = dt_by_ofs[j];

        recv_msgs[j + nrecv].rank
          = (comm_size - (int)j - 1 + i)%comm_size;
        recv_msgs[j + nrecv].datatype = dt_by_ofs[j];
      }
    }

    enum { exchange_fwd, exchange_rev, num_exchanges };
    Xt_exchanger exchanger[num_exchanges];
    for (size_t j = 0; j < num_exchanges; ++j)
      exchanger[j] = exchanger_new(nsend, (int)nrecv,
                                   send_msgs, recv_msgs + j * nrecv, comm, 0,
                                   config);

    // test
    int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
    for (int async = 0; async < 1 + test_async; ++async) {
      int src_data[1] = {(my_rank+comm_size-i)%comm_size};
      for (size_t j = 0; j < 2 * nrecv; ++j)
        dst_data[j] = -1;

      for (size_t j = 0; j < num_exchanges; ++j) {
        if (async) {
          Xt_request request;
          int flag;
          xt_exchanger_a_exchange(exchanger[j], src_data,
                                  dst_data + j * nrecv, &request);
          xt_request_test(&request, &flag);
          xt_request_wait(&request);
          xt_request_test(&request, &flag);
          if (!flag) PUT_ERR("invalid flag result\n");
        } else {
          xt_exchanger_s_exchange(exchanger[j], src_data, dst_data + j * nrecv);
        }
      }

      bool mismatch = false;
      for (size_t j = 0; j < nrecv; ++j)
        mismatch |= (((size_t)dst_data[j] != j + 1)
                     || (dst_data[nrecv + j] != comm_size - (int)j - 1));
      if (mismatch)
        PUT_ERR("invalid data\n");
    }

    // cleanup
    free(recv_msgs);
    xt_exchanger_delete(exchanger[0]);
    xt_exchanger_delete(exchanger[1]);
  }
  free(dst_data);
  for (size_t j = 1; j < (size_t)comm_size; ++j)
    MPI_Type_free(dt_by_ofs + j);
  free(dt_by_ofs);
}

static void
test_all2all(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // all-to-all pattern
  // setup
  int nsend = comm_size - 1;
  int nrecv = comm_size - 1;

  struct Xt_redist_msg *send_msgs = xmalloc((size_t)nsend * sizeof (*send_msgs));
  for (size_t j = 0; j < (size_t)nsend; ++j)
    send_msgs[j].datatype = MPI_INT;
  struct Xt_redist_msg *recv_msgs
    = xmalloc((size_t)nrecv * sizeof (*recv_msgs));
  int *dst_data = xmalloc((size_t)comm_size * sizeof (*dst_data));
  MPI_Datatype *dt_by_ofs = xmalloc((size_t)comm_size * sizeof (*dt_by_ofs));
  dt_by_ofs[0] = my_rank != 0 ? MPI_INT : MPI_DATATYPE_NULL;
  for (size_t i = 1; i < (size_t)comm_size; ++i)
    if (i != (size_t)my_rank) {
      MPI_Type_indexed(1, (int[]){1}, (int[]){(int)i}, MPI_INT, dt_by_ofs + i);
      MPI_Type_commit(dt_by_ofs + i);
    } else
      dt_by_ofs[i] = MPI_DATATYPE_NULL;
  size_t comm_size_ = (size_t)comm_size,
    incr = (size_t)(comm_size/test_freq + (comm_size < test_freq));
  for (size_t i = 0; i < comm_size_ - 1; i += incr) {
    for (size_t j = 0; j < (size_t)nsend; ++j) {
      int ofs = my_rank + 1 + (int)i + (int)j;
      send_msgs[j].rank = (ofs + (ofs >= comm_size + my_rank))%comm_size;
    }
    for (size_t j = 0; j < comm_size_ - 1; j += incr) {
      for (size_t k = 0; k < (size_t)nrecv; ++k) {
        int ofs = ((int)i + (int)j + (int)k)%(comm_size - 1);
        ofs += ofs >= my_rank;
        recv_msgs[k].rank = ofs;
        recv_msgs[k].datatype = dt_by_ofs[ofs];
      }
      Xt_exchanger exchanger = exchanger_new(nsend, nrecv,
                                             send_msgs,
                                             recv_msgs,
                                             MPI_COMM_WORLD, 0, config);
      // test
      int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
      for (int async = 0; async < 1 + test_async; ++async) {
        int src_data[1] = {my_rank};
        for (size_t k = 0; k < comm_size_; ++k)
          dst_data[k] = my_rank;

        if (async) {
          Xt_request request;
          int flag;
          xt_exchanger_a_exchange(exchanger, (void*)src_data, (void*)dst_data,
                                  &request);
          xt_request_test(&request, &flag);
          xt_request_wait(&request);
          xt_request_test(&request, &flag);
          if (!flag) PUT_ERR("invalid flag result\n");
        } else {
          xt_exchanger_s_exchange(exchanger, (void*)src_data, (void*)dst_data);
        }

        bool mismatch = false;
        for (int k = 0; k < comm_size; ++k)
          mismatch |= (dst_data[k] != k);
        if (mismatch)
          PUT_ERR("invalid data\n");
      }

      // cleanup
      xt_exchanger_delete(exchanger);
    }
  }
  for (size_t i = 1; i < comm_size_; ++i)
    if (i != (size_t)my_rank)
      MPI_Type_free(dt_by_ofs + i);
  free(dt_by_ofs);
  free(dst_data);
  free(recv_msgs);
  free(send_msgs);
}

static void
test_rr(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // round robin pattern
  int incr = comm_size/test_freq + (comm_size < test_freq);
  for (int i = 1; i < comm_size; i += incr) {

    // setup
    enum { nsend = 1, nrecv = 1 };
    struct Xt_redist_msg send_msgs[nsend]
      = {{.rank=(my_rank + i)%comm_size, .datatype=MPI_INT}};
    struct Xt_redist_msg recv_msgs[nrecv]
      = {{.rank=(my_rank + comm_size - i)%comm_size, .datatype=MPI_INT}};

    Xt_exchanger exchanger = exchanger_new(nsend, nrecv, send_msgs,
                                           recv_msgs, comm, 0, config);

    // test
    int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
    for (int async = 0; async < 1 + test_async; ++async) {
      int src_data[1] = {my_rank};
      int dst_data[1] = {-1};

      if (async) {
        Xt_request request;
        int flag;
        xt_exchanger_a_exchange(exchanger, src_data, dst_data, &request);
        xt_request_test(&request, &flag);
        xt_request_wait(&request);
        xt_request_test(&request, &flag);
        if (!flag) PUT_ERR("invalid flag result\n");
      } else {
        xt_exchanger_s_exchange(exchanger, src_data, dst_data);
      }

      if (dst_data[0] != (my_rank + comm_size - i)%comm_size)
        PUT_ERR("invalid data\n");
    }

    // cleanup
    xt_exchanger_delete(exchanger);
  }
}

static void
test_intercomm_all2all(MPI_Comm comm, Xt_exchanger_new exchanger_new, Xt_config config)
{
#if MPI_VERSION >= 3
  // inter-communicator's are not defined for virtual topologies, which are
  // used by xt_exchanger_neigh_alltoall_new
  if (exchanger_new == xt_exchanger_neigh_alltoall_new) return;
#endif

  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  if (comm_size == 1) return;

  // generate intercomm
  int splitRank = (comm_size * 2) / 3;
  int group = my_rank >= splitRank;
  MPI_Comm intra_group_comm;
  MPI_Comm inter_comm; // split communicator with 2 to 1 ratio
  xt_mpi_call(MPI_Comm_split(comm, group, 0, &intra_group_comm), comm);
  xt_mpi_call(MPI_Intercomm_create(intra_group_comm, 0, comm,
                                   group ? 0 : splitRank,
                                   0, &inter_comm), comm);

  int intra_rank;
  int local_size, remote_size;
  xt_mpi_call(MPI_Comm_rank(inter_comm, &intra_rank), comm);
  xt_mpi_call(MPI_Comm_size(inter_comm, &local_size), comm);
  xt_mpi_call(MPI_Comm_remote_size(inter_comm, &remote_size), comm);

  // all-to-all pattern
  // setup
  int nsend = remote_size;
  int nrecv = remote_size;

  struct Xt_redist_msg * send_msgs =
    xmalloc((size_t)nsend * sizeof (*send_msgs));
  for (int i = 0; i < nsend; ++i) {
    send_msgs[i].rank = i;
    send_msgs[i].datatype = MPI_INT;
  }

  struct Xt_redist_msg * recv_msgs =
    xmalloc((size_t)nrecv * sizeof (*recv_msgs));
  for (int i = 0; i < nrecv; ++i) {
    recv_msgs[i].rank = i;
    recv_msgs[i].datatype = MPI_INT;
    MPI_Type_indexed(
      1, (int[]){1}, (int[]){i}, MPI_INT, &(recv_msgs[i].datatype));
    MPI_Type_commit(&(recv_msgs[i].datatype));
  }

  int *dst_data = xmalloc((size_t)nrecv * sizeof (*dst_data));

  Xt_exchanger exchanger = exchanger_new(nsend, nrecv,
                                         send_msgs,
                                         recv_msgs,
                                         inter_comm, 0, config);
  // test
  int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
  for (int async = 0; async < 1 + test_async; ++async) {

    int src_data[1] = {my_rank};
    for (int i = 0; i < nrecv; ++i) dst_data[i] = -1;

    if (async) {
      Xt_request request;
      int flag;
      xt_exchanger_a_exchange(exchanger, (void*)src_data, (void*)dst_data,
                              &request);
      xt_request_test(&request, &flag);
      xt_request_wait(&request);
      xt_request_test(&request, &flag);
      if (!flag) PUT_ERR("invalid flag result\n");
    } else {
      xt_exchanger_s_exchange(exchanger, (void*)src_data, (void*)dst_data);
    }

    int dst_data_offset = (my_rank >= splitRank)?0:splitRank;
    bool mismatch = false;
    for (int i = 0; i < nrecv; ++i)
      mismatch |= (dst_data[i] != i + dst_data_offset);
    if (mismatch)
      PUT_ERR("invalid data\n");
  }

  // cleanup
  xt_exchanger_delete(exchanger);

  for (int i = 0; i < nrecv; ++i) MPI_Type_free(&(recv_msgs[i].datatype));

  free(dst_data);
  free(recv_msgs);
  free(send_msgs);

  xt_mpi_call(MPI_Comm_free(&inter_comm), comm);
  xt_mpi_call(MPI_Comm_free(&intra_group_comm), comm);
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
