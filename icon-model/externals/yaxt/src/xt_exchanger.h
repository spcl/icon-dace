/**
 * @file xt_exchanger.h
 * @brief exchanging of data based on information provided by redist's
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
#ifndef XT_EXCHANGER_H
#define XT_EXCHANGER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <mpi.h>

#include "xt/xt_config.h"
#include "xt/xt_core.h"
#include "xt/xt_request.h"
#include "xt_redist_internal.h"
#include "core/ppm_visibility.h"


/** \example test_exchanger_parallel.c
 */

typedef struct Xt_exchanger_omp_share_ *Xt_exchanger_omp_share;

typedef struct Xt_exchanger_ *Xt_exchanger;

struct xt_exchanger_vtable {
  Xt_exchanger (*copy)(Xt_exchanger, MPI_Comm, int);
  void (*delete)(Xt_exchanger);
  void (*s_exchange)(Xt_exchanger, const void *, void *);
  void (*a_exchange)(Xt_exchanger, const void *, void *, Xt_request *request);
  int (*get_msg_ranks)(Xt_exchanger, enum xt_msg_direction, int *restrict *);
  MPI_Datatype (*get_MPI_Datatype)(Xt_exchanger, int, enum xt_msg_direction,
                                   bool);
  void (*team_share_default_init)(void *share);
  void (*team_share_destroy)(void *share);
  size_t team_share_size;
  Xt_exchanger_omp_share (*create_omp_share)(Xt_exchanger);
};

struct Xt_exchanger_ {
  struct xt_exchanger_vtable * vtable;
};

/**
 * Executes a synchronous data exchange.
 * @param[in]  exchanger exchanger object
 * @param[in]  src_data  source data
 * @param[out] dst_data  destination data
 */
PPM_DSO_INTERNAL void
xt_exchanger_s_exchange(Xt_exchanger exchanger, const void * src_data,
                        void * dst_data);

/**
 * Executes a asynchronous data exchange.
 * @param[in]  exchanger exchanger object
 * @param[in]  src_data  source data
 * @param[out] dst_data  destination data
 * @param[out] request   pointer that will reference request object
 *                       created for duration of exchange
 */
PPM_DSO_INTERNAL void
xt_exchanger_a_exchange(Xt_exchanger exchanger,
                        const void * src_data, void * dst_data,
                        Xt_request *request);

/**
 * Copies an exchange object.
 * @param[in] orig           exchanger object to be copied
 * @param[in] new_comm       communicator to be used by the new exchanger
 * @param[in] new_tag_offset new tag_offset
 */
PPM_DSO_INTERNAL Xt_exchanger
xt_exchanger_copy(Xt_exchanger orig, MPI_Comm new_comm, int new_tag_offset);

/**
 * Destructor for an exchanger.
 */
PPM_DSO_INTERNAL void
xt_exchanger_delete(Xt_exchanger);

/**
 * gets a copy of the MPI_Datatype used for a specificed message
 * @param[in] exchanger exchanger object
 * @param[in] rank      MPI rank
 * @param[in] direction specific whether the datatype of an incoming or outgoing
 *                      message is requested
 * @param[in] do_dup    mpi datatype copy will be dup if true
 * @return MPI_Datatype for the specificed message
 * @remark returns MPI_DATATYPE_NULL if there is no message matching the
 *         specificed configuration
 */
PPM_DSO_INTERNAL MPI_Datatype
xt_exchanger_get_MPI_Datatype(Xt_exchanger exchanger,
                              int rank, enum xt_msg_direction direction,
                              bool do_dup);

/**
 * Gets the ranks of all processes that receive data from/send data to the local
 * process in the specificed message.
 * @param[in]  exchanger exchanger object
 * @param[in]  direction specifices whether ranks for the outgoing or incoming
 *                       messages are requested
 * @param[out] ranks     ranks for all outgoing/incoming messages
 * @return number of outgoing/incoming message
 * @remark the user needs to ensure that array ranks is big enough to hold all
 *         ranks
 */
PPM_DSO_INTERNAL int
xt_exchanger_get_msg_ranks(Xt_exchanger exchanger,
                           enum xt_msg_direction direction,
                           int *restrict *ranks);

typedef Xt_exchanger
(*Xt_exchanger_new)(int nsend, int nrecv,
                    const struct Xt_redist_msg *send_msgs,
                    const struct Xt_redist_msg *recv_msgs,
                    MPI_Comm comm, int tag_offset,
                    Xt_config config);

/**
 * Given an exchanger constructor, query the vtable pointer
 *
 * @param[in] exchanger_new function to create new exchanger to query
 *                          vtable pointer for
 * @return vtable pointer
 */
PPM_DSO_INTERNAL const struct xt_exchanger_vtable *
xt_exchanger_new_get_vtable(Xt_exchanger_new exchanger_new);

/* team shares refer to expensive data that can be shared among
 * multiple instances of exchangers that are of the same class and
 * derived from the same xmap or communication matrix */

enum { team_share_align = sizeof (void *) };

/**
 * Given an exchanger, query the size of shared data for a team.
 *
 * @param[in] exchanger exchanger object to query team shared data size for.
 * @return size of shared data
 */
PPM_DSO_INTERNAL size_t
xt_exchanger_team_get_share_size(Xt_exchanger exchanger);

/**
 * Given an exchanger constructor, query the size of shared data for a
 * team.
 *
 * @param[in] exchanger_new function to create new exchanger to query
 *                          team shared data size for.
 * @return size of shared data
 */
PPM_DSO_INTERNAL size_t
xt_exchanger_new_team_get_share_size(Xt_exchanger_new exchanger_new);

/**
 * Given an exchanger, initialize shared data for a team to default values.
 *
 * @param[in] exchanger exchanger object to default initialize team
 * shared data for.
 * @param[out] share    object to initialize
 */
PPM_DSO_INTERNAL void
xt_exchanger_team_share_default_init(Xt_exchanger exchanger,
                                     void *share);
/**
 * Given an exchanger constructor, initialize shared data for a team
 * to default values.
 *
 * @param[in] exchanger_new exchanger constructor to default initialize team
 * shared data for.
 * @param[out] share    object to initialize
 */
PPM_DSO_INTERNAL void
xt_exchanger_new_team_share_default_init(Xt_exchanger_new exchanger_new,
                                         void *share);

/**
 * Given an exchanger, destroy shared data for a team and reset to default
 * values.
 * This call is collective for all MPI ranks in the communicator the
 * exchanger was constructed for.
 *
 * @param[in] exchanger exchanger object to destroy team shared data for.
 * @param[out] share    object to destroy
 */
PPM_DSO_INTERNAL void
xt_exchanger_team_share_destroy(Xt_exchanger exchanger, void *share);

/**
 * Given an exchanger constructor, destroy shared data for a team and
 * reset to default values.
 *
 * This call is collective for all MPI ranks
 * in the communicator the exchanger team was constructed for.
 *
 * @param[in] exchanger_new exchanger constructor to destroy team
 *                          shared data for.
 * @param[out] share        shared state object to destroy
 */
PPM_DSO_INTERNAL void
xt_exchanger_new_team_share_destroy(Xt_exchanger_new exchanger_new,
                                    void *share);


/**
 * Given an exchanger, create an object holding shared data for
 * multi-threaded invocation.
 *
 * @param[in] exchanger exchanger to prepare the shared data for.
 */
PPM_DSO_INTERNAL Xt_exchanger_omp_share
xt_exchanger_create_omp_share(Xt_exchanger exchanger);




#endif // XT_EXCHANGER_H
/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
