/**
 * @file test_idxlist_collection.c
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
#include <limits.h>
#include <stdlib.h>

#include <mpi.h>
#include <yaxt.h>

#define VERBOSE

#include "tests.h"
#include "test_idxlist_utils.h"
#include "core/ppm_xfuncs.h"

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // test packing and unpacking of collection index lists

    enum {
      num_indices = 7,
      num_ref_stripes = 2,
      num_lists = 2,
    };
    static const Xt_int index_list[num_lists][num_indices]
      = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};

    // create two index lists
    Xt_idxlist idxlists[num_lists];
    for (size_t i = 0; i < num_lists; ++i)
      idxlists[i] = xt_idxvec_new(index_list[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, num_lists);
    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, index_list[0], (int)(num_lists*num_indices));

    static const struct Xt_stripe ref_stripes[num_ref_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 7},
         {.start = 7, .stride = -1, .nstrides = 7}};
    Xt_idxlist collectionlist_copy
      = idxlist_pack_unpack_copy(collectionlist);
    check_idxlist_copy(collectionlist, collectionlist_copy,
                       2 * (size_t)num_indices, index_list[0],
                       (size_t)num_ref_stripes, ref_stripes);
    xt_idxlist_delete(collectionlist_copy);
    xt_idxlist_delete(collectionlist);
  }

  { // test copying of collection index lists

    enum {
      num_lists = 2,
      num_indices = 7,
      num_ref_stripes = 2,
    };

    static const Xt_int index_list[num_lists][num_indices]
      = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};

    Xt_idxlist idxlists[num_lists];

    // create two index lists
    for (size_t i = 0; i < num_lists; ++i)
      idxlists[i] = xt_idxvec_new(index_list[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, 2);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, &index_list[0][0], 2 * num_indices);

    // copy the index list
    Xt_idxlist collectionlist_copy = xt_idxlist_copy(collectionlist);

    // check copy of collection list
    check_idxlist(collectionlist_copy, &(index_list[0][0]),
                  2 * num_indices);

    // check the conversion to stripes

    struct Xt_stripe *stripes;
    int num_stripes;

    xt_idxlist_get_index_stripes(collectionlist, &stripes, &num_stripes);

    static const struct Xt_stripe ref_stripes[2]
      = {{.start = 1, .stride = 1, .nstrides = 7},
         {.start = 7, .stride = -1, .nstrides = 7}};

    check_stripes(stripes, num_stripes, ref_stripes, num_ref_stripes);

    // clean up
    free(stripes);

    xt_idxlist_delete(collectionlist);
    xt_idxlist_delete(collectionlist_copy);
  }

  { // test computing intersection of index list collections

    enum {
      num_indices = 7,
      num_vecs = 3
    };
    static const Xt_int index_list[num_vecs][num_indices]
      = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1},{2,6,1,4,7,3,0}},
      sorted_index_list[num_vecs * num_indices]
        = {0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,6,6,7,7,7};

    Xt_idxlist idxlists[num_vecs];

    // create two index lists
    for (size_t i = 0; i < num_vecs; ++i)
      idxlists[i] = xt_idxvec_new(index_list[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, 3);

    for (size_t i = 0; i < num_vecs; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, &(index_list[0][0]), (int)(3*num_indices));

    // test computing of intersection
    Xt_idxlist ref_idxvec = xt_idxvec_new(&(index_list[0][0]), 3*num_indices),
      intersection = xt_idxlist_get_intersection(ref_idxvec, collectionlist);

    check_idxlist(intersection, sorted_index_list, (int)(num_vecs*num_indices));

    xt_idxlist_delete(intersection);

    intersection = xt_idxlist_get_intersection(collectionlist, ref_idxvec);

    check_idxlist(intersection, sorted_index_list, (int)(num_vecs*num_indices));

    xt_idxlist_delete(intersection);
    xt_idxlist_delete(ref_idxvec);
    xt_idxlist_delete(collectionlist);
  }

  { // test using various tpyes of index lists

    enum {
      num_lists = 3,
      section_ndims = 2,
    };
    static const Xt_int index_list[] = {1,3,5,7,9,11};
    static const struct Xt_stripe stripes[]
      = {{.start = 0, .stride = 2, .nstrides = 5},
         {.start = 1, .stride = 2, .nstrides = 5}};
    static const Xt_int global_size[section_ndims] = {10,10};
    static const int local_size[section_ndims] = {5,5};
    static const Xt_int local_start[section_ndims] = {2,2};

    static const Xt_int ref_index_list[] = {1,3,5,7,9,11,
                                            0,2,4,6,8,1,3,5,7,9,
                                            22,23,24,25,26,
                                            32,33,34,35,36,
                                            42,43,44,45,46,
                                            52,53,54,55,56,
                                            62,63,64,65,66};

    Xt_idxlist idxlists[num_lists]
      = {xt_idxvec_new(index_list, sizeof(index_list)/sizeof(index_list[0])),
         xt_idxstripes_new(stripes, sizeof(stripes)/sizeof(stripes[0])),
         xt_idxsection_new(0, section_ndims, global_size, local_size, local_start)};

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, num_lists);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, ref_index_list,
                  (int)(sizeof(ref_index_list)/sizeof(ref_index_list[0])));

    xt_idxlist_delete(collectionlist);
  }

  { // check get_bounding_box
    enum { num_lists = 2 };
    Xt_idxlist idxlists[num_lists] = {xt_idxempty_new(),
                                      xt_idxempty_new()};

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, num_lists);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    enum { ndim = 3 };
    Xt_int global_size_bb[ndim];
    static const Xt_int global_start_index = 0;
    struct Xt_bounds bounds[ndim];

    for (size_t i = 0; i < ndim; ++i)
      global_size_bb[i] = 4;

    xt_idxlist_get_bounding_box(collectionlist, ndim, global_size_bb,
                                global_start_index, bounds);

    bool mismatch = false;
    for (size_t i = 0; i < ndim; ++i)
      mismatch |= (bounds[i].size != 0);

    if (mismatch)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(collectionlist);
  }

  {
    enum {
      num_lists = 2,
      num_indices = 3,
    };
    // check get_bounding_box
    static const Xt_int indices[num_lists][num_indices]
      = {{45,35,32},{32,48,33}};

    Xt_idxlist idxlists[num_lists];

    // create two index lists
    for (size_t i = 0; i < num_lists; ++i)
      idxlists[i] = xt_idxvec_new(indices[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, 2);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    enum { ndim = 3 };
    static const Xt_int global_size[ndim] = { 5, 4, 3 };
    static const Xt_int global_start_index = 1;
    struct Xt_bounds bounds[ndim];

    xt_idxlist_get_bounding_box(collectionlist, ndim, global_size,
                                global_start_index, bounds);

    static const Xt_int ref_start[3] = { 2, 2, 1 };

    bool mismatch = false;
    for (size_t i = 0; i < ndim; ++i)
      mismatch |= (bounds[i].size != 2 || bounds[i].start != ref_start[i]);

    if (mismatch)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(collectionlist);
  }

  xt_finalize();
  xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
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
