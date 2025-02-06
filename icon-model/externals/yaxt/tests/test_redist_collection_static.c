/**
 * @file test_redist_collection_static.c
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

#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "ctest_common.h"
#include "test_redist_common.h"

static void
simple_test(MPI_Comm comm, Xt_config config);

static void
test_repeated_redist(MPI_Comm comm, Xt_config config);

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);
  Xt_config config = redist_exchanger_option(&argc, &argv);

  simple_test(MPI_COMM_WORLD, config);
  test_repeated_redist(MPI_COMM_WORLD, config);

  xt_config_delete(config);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
simple_test(MPI_Comm comm, Xt_config config)
{ // general test with one redist
  // set up data
  enum { src_slice_len = 5, dst_slice_len = 3 };
  Xt_xmap xmap = build_odd_selection_xmap(src_slice_len, comm);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection

  static const MPI_Aint displacements[1] = {0};

  Xt_redist redist_coll
    = xt_redist_collection_static_custom_new(&redist, 1, displacements,
                                             displacements, comm, config);

  // test communicator of redist

  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);

  // test exchange
  static const double src_data[] = {1,2,3,4,5};
  static const double ref_dst_data[] = {1,3,5};
  enum { num_dst_elems = sizeof (ref_dst_data) / sizeof (ref_dst_data[0]) };

  double dst_data[num_dst_elems];
  check_redist(redist_coll, src_data,
               num_dst_elems, dst_data, fill_array_double, NULL,
               ref_dst_data, MPI_DOUBLE, MPI_DOUBLE);

  // clean up
  xt_redist_delete(redist_coll);
}

static void
test_repeated_redist(MPI_Comm comm, Xt_config config)
{ // test with one redist used three times (two exchanges)
  // set up data
  enum { src_slice_len = 5,
         dst_slice_len = (src_slice_len+1)/2 };
  Xt_xmap xmap = build_odd_selection_xmap(src_slice_len, comm);
  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection
  Xt_redist redists[3] = {redist, redist, redist};
  static const double src_data[3][src_slice_len]
    = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
  double dst_data[3][3];
  MPI_Aint src_displacements[3]
    = {0, (MPI_Aint)((size_t)(src_data[1] - src_data[0]) * sizeof (double)),
       (MPI_Aint)((size_t)(src_data[2] - src_data[0]) * sizeof(double)) };
  MPI_Aint dst_displacements[3]
    = {0, (MPI_Aint)((size_t)(dst_data[1] - dst_data[0]) * sizeof(double)),
       (MPI_Aint)((size_t)(dst_data[2] - dst_data[0]) * sizeof(double)) };

  Xt_redist redist_coll
    = xt_redist_collection_static_custom_new(redists, 3, src_displacements,
                                             dst_displacements, comm, config);

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);

  // test exchange
  static const double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};
  check_redist(redist_coll, src_data,
               sizeof (dst_data) / sizeof (dst_data[0][0]),
               dst_data, fill_array_double, NULL,
               ref_dst_data, MPI_DOUBLE, MPI_DOUBLE);

  // clean up
  xt_redist_delete(redist_coll);
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
