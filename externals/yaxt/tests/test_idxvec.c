/**
 * @file test_idxvec.c
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

#ifdef NDEBUG
#undef NDEBUG
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

static void
test_intersection(int num_idx_a, const Xt_int index_vector_a[],
                  int num_idx_b, const Xt_int index_vector_b[],
                  int num_idx_isect, const Xt_int index_vector_isect[]);
static void
test_idxvec_from_stripes(int num_stripes, const struct Xt_stripe stripes[],
                         int num_ref_indices, const Xt_int ref_indices[]);
static void
test_get_indices_at_positions(int num_idx, const Xt_int indices[],
                              Xt_int undef_idx,
                              int num_pos, const int pos[]);

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // test packing and unpacking of index vectors

    static const Xt_int index_vector[] = {1,2,3,4,5,6,7};
    enum { num_indices = sizeof(index_vector) / sizeof(index_vector[0]) };

    // create an index vector
    Xt_idxlist idxvector = xt_idxvec_new(index_vector, num_indices);

    // test the generated index vector
    check_idxlist(idxvector, index_vector, num_indices);

    Xt_idxlist idxvector_copy = idxlist_pack_unpack_copy(idxvector);

    // check the received index vector
    check_idxlist(idxvector_copy, index_vector, num_indices);

    // compute intersection between the two index vectors
    Xt_idxlist intersection
      = xt_idxlist_get_intersection(idxvector, idxvector_copy);

    // check the computed intersection (should be identically to the original
    // index vector)
    check_idxlist(intersection, index_vector, num_indices);

    // check the conversion to stripes
    struct Xt_stripe *stripes;
    int num_stripes;
    xt_idxlist_get_index_stripes(idxvector, &stripes, &num_stripes);

    static const struct Xt_stripe ref_stripes[]
      = {{.start = 1, .stride = 1, .nstrides = 7}};
    enum { ref_num_stripes = sizeof (ref_stripes) / sizeof (ref_stripes[0]) };

    check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

    // clean up
    free(stripes);
    xt_idxlist_delete(idxvector);
    xt_idxlist_delete(idxvector_copy);
    xt_idxlist_delete(intersection);
  }

  { // test copying of index vector

    static const Xt_int index_vector[] = {1,2,3,4,5,6,7};
    enum { num_indices = sizeof(index_vector) / sizeof(index_vector[0]) };

    // create an index vector

    Xt_idxlist idxvector = xt_idxvec_new(index_vector, num_indices);

    // test the generated index vector
    check_idxlist(idxvector, index_vector, num_indices);

    // copy the index vector
    Xt_idxlist idxvector_copy = xt_idxlist_copy(idxvector);

    // check the generated copy
    check_idxlist(idxvector_copy, index_vector, num_indices);

    // clean up
    xt_idxlist_delete(idxvector);
    xt_idxlist_delete(idxvector_copy);
  }

  { // test repeated equal indices

    static const Xt_int index_vector[] = {1,2,3,7,5,6,7,7};
    enum { num_indices = sizeof(index_vector) / sizeof(index_vector[0]) };

    // create an index vector
    Xt_idxlist idxvector = xt_idxvec_new(index_vector, num_indices);

    // test the generated index vector
    check_idxlist(idxvector, index_vector, num_indices);

    // clean up
    xt_idxlist_delete(idxvector);
  }

  { // subtest from test_ut

    enum { single_match_only = 1 };
    static const Xt_int index_vector[]
      = { 10, 15, 14, 13, 12, 15, 10, 11, 12, 13, 23, 18, 19,
          20, 21, 31, 26, 27, 28, 29 };

    static const Xt_int intersection_vector[]
      = { 12, 12, 13, 13, 14, 15, 15, 20, 21, 23, 28, 29, 31 };
    enum {
      index_num = sizeof(index_vector) / sizeof(index_vector[0]),
      intersection_num
      = sizeof(intersection_vector) / sizeof(intersection_vector[0]),
    };

    // create an index vector
    Xt_idxlist idxvector = xt_idxvec_new(index_vector, index_num);

    // get positions:
    int intersection_pos[intersection_num];
    static const int ref_intersection_pos[]
      = { 4, 8, 3, 9, 2, 1, 5, 13, 14, 10, 18, 19, 15 };
    enum {
      ref_intersection_num
      = sizeof(ref_intersection_pos) / sizeof(ref_intersection_pos[0])
    };

    assert((int)intersection_num == (int)ref_intersection_num);

    int retval = xt_idxlist_get_positions_of_indices(idxvector,
                                                     intersection_vector,
                                                     intersection_num,
                                                     intersection_pos,
                                                     single_match_only);
    assert(retval == 0);

    // check positions:
    check_offsets(intersection_num, intersection_pos, ref_intersection_pos);

    // clean up
    xt_idxlist_delete(idxvector);
  }

  {
    enum {
      num_vec = 2,
      num_idx = 3,
    };
    // check intersection
    static Xt_int index_vector[num_vec][num_idx] = {{1,2,3},{1,2,3}};

    test_intersection(num_idx, index_vector[0],
                      num_idx, index_vector[1],
                      num_idx, index_vector[0]);
  }

  {
    enum {
      num_vec = 2,
      num_idx = 3,
    };
    // check intersection
    static const Xt_int index_vector[num_vec][num_idx] = {{1,2,3},{2,3,4}};
    static const Xt_int ref_intersection_indices[] = {2,3};
    enum {
      num_ref_intersection_indices
      = sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0])
    };

    test_intersection(num_idx, index_vector[0],
                      num_idx, index_vector[1],
                      num_ref_intersection_indices, ref_intersection_indices);
  }

  {
    enum {
      num_vec = 2,
      num_idx = 3
    };
    // check intersection
    static const Xt_int index_vector[num_vec][num_idx] = {{2,3,4},{1,2,3}};
    static const Xt_int ref_intersection_indices[] = {2,3};
    enum {
      num_ref_intersection_indices
      = sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0])
    };
    test_intersection(num_idx, index_vector[0],
                      num_idx, index_vector[1],
                      num_ref_intersection_indices, ref_intersection_indices);
  }

  {
    enum {
      num_vec = 2,
      num_idx = 3
    };
    // check intersection
    static const Xt_int index_vector[num_vec][num_idx] = {{4,2,3},{3,1,2}};
    static const Xt_int ref_intersection_indices[] = {2,3};
    enum {
      num_ref_intersection_indices
      = sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0])
    };

    test_intersection(num_idx, index_vector[0],
                      num_idx, index_vector[1],
                      num_ref_intersection_indices, ref_intersection_indices);
  }

  {
    enum {
      num_vec = 2,
      num_idx = 3
    };
    // check intersection
    static const Xt_int index_vector[2][3] = {{3,1,2},{4,2,3}};
    static const Xt_int ref_intersection_indices[] = {2,3};
    enum {
      num_ref_intersection_indices
      = sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0])
    };
    test_intersection(num_idx, index_vector[0],
                      num_idx, index_vector[1],
                      num_ref_intersection_indices, ref_intersection_indices);
  }

  {
    // test constructor xt_idxvec_from_stripes_new
    static const struct Xt_stripe stripes[]
      = {{.start = 5, .stride =  1, .nstrides = 5},
         {.start = 4, .stride = -1, .nstrides = 5}};
    static const Xt_int ref_indices[] = {5,6,7,8,9,4,3,2,1,0};
    enum {
      num_stripes = sizeof (stripes) / sizeof (stripes[0]),
      num_ref_indices = sizeof (ref_indices) / sizeof (ref_indices[0]),
    };
    test_idxvec_from_stripes(num_stripes, stripes,
                             num_ref_indices, ref_indices);
  }

  {
    // test constructor xt_idxvec_from_stripes_new
    static const struct Xt_stripe stripes[]
      = {{.start = 0, .stride = 1, .nstrides = 5},
         {.start = 2, .stride = 1, .nstrides = 5},
         {.start = 4, .stride = 1, .nstrides = 5}};
    static const Xt_int ref_indices[]
      = { 0, 1, 2, 3, 4, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8 };
    enum {
      num_stripes = sizeof (stripes) / sizeof (stripes[0]),
      num_ref_indices = sizeof (ref_indices) / sizeof (ref_indices[0])
    };
    test_idxvec_from_stripes(num_stripes, stripes,
                             num_ref_indices, ref_indices);
  }

  {
    // test constructor xt_idxvec_from_stripes_new
    static const struct Xt_stripe stripes[]
      = {{.start = 2, .stride = 1, .nstrides = 5},
         {.start = 0, .stride = 1, .nstrides = 5},
         {.start = 4, .stride = 1, .nstrides = 5}};
    static const Xt_int ref_indices[]
      = { 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 4, 5, 6, 7, 8 };
    enum {
      num_stripes = sizeof (stripes) / sizeof (stripes[0]),
      num_ref_indices = sizeof (ref_indices) / sizeof (ref_indices[0])
    };
    test_idxvec_from_stripes(num_stripes, stripes,
                             num_ref_indices, ref_indices);
  }

  {
    // test constructor xt_idxvec_from_stripes_new
    static const struct Xt_stripe stripes[]
      = {{.start = 2, .stride =  1, .nstrides = 5},
         {.start = 4, .stride = -1, .nstrides = 5},
         {.start = 4, .stride =  1, .nstrides = 5}};
    static const Xt_int ref_indices[]
      = { 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 4, 5, 6, 7, 8 };
    enum {
      num_stripes = sizeof (stripes) / sizeof (stripes[0]),
      num_ref_indices = sizeof (ref_indices) / sizeof (ref_indices[0])
    };
    test_idxvec_from_stripes(num_stripes, stripes,
                             num_ref_indices, ref_indices);
  }

  {
    // test constructor xt_idxvec_from_stripes_new
    static const struct Xt_stripe stripes[]
      = {{.start =  0, .stride =  3, .nstrides = 5},
         {.start =  1, .stride =  3, .nstrides = 5},
         {.start = 14, .stride = -3, .nstrides = 5}};
    static const Xt_int ref_indices[]
      = { 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2 };
    enum {
      num_stripes = sizeof (stripes) / sizeof (stripes[0]),
      num_ref_indices = sizeof (ref_indices) / sizeof (ref_indices[0])
    };
    test_idxvec_from_stripes(num_stripes, stripes,
                             num_ref_indices, ref_indices);
  }

  {
    // test constructor xt_idxvec_from_stripes_new
    static const struct Xt_stripe stripes[]
      = {{.start =  0, .stride =  3, .nstrides = 5},
         {.start =  2, .stride =  3, .nstrides = 5},
         {.start = 14, .stride = -3, .nstrides = 5}};
    static const Xt_int ref_indices[]
      = { 0, 3, 6, 9, 12, 2, 5, 8, 11, 14, 14, 11, 8, 5, 2 };
    enum {
      num_stripes = sizeof (stripes) / sizeof (stripes[0]),
      num_ref_indices = sizeof (ref_indices) / sizeof (ref_indices[0])
    };
    test_idxvec_from_stripes(num_stripes, stripes,
                             num_ref_indices, ref_indices);
  }

  {
    // test constructor xt_idxvec_from_stripes_new
    static const struct Xt_stripe stripes[]
      = {{.start =  0, .stride = -1, .nstrides = 5},
         {.start =  1, .stride =  1, .nstrides = 5},
         {.start = -5, .stride = -1, .nstrides = 5},
         {.start =  6, .stride =  1, .nstrides = 5}};
    static const Xt_int ref_indices[]
      = {  0, -1, -2, -3, -4,  1,  2,  3,  4,  5,
          -5, -6, -7, -8, -9,  6,  7,  8,  9,  10 };
    enum {
      num_stripes = sizeof (stripes) / sizeof (stripes[0]),
      num_ref_indices = sizeof (ref_indices) / sizeof (ref_indices[0])
    };
    test_idxvec_from_stripes(num_stripes, stripes,
                             num_ref_indices, ref_indices);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: mixed valid and invalid positions
    Xt_int undef_idx = XT_INT_MIN;
    static const Xt_int indices[]
      = { 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2, 1 };
    static const int pos[]
      = { 0, 2, 7, 9, 11, 100, 11, 200, 9, 300, 7, 400, 5 };
    enum {
      num_indices = sizeof(indices) / sizeof(indices[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    test_get_indices_at_positions(num_indices, indices, undef_idx,
                                  num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: only valid positions
    Xt_int undef_idx = XT_INT_MIN;
    static const Xt_int indices[]
      = { 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2, 1 };
    static const int pos[]
      = { 0, 2, 7, 9, 11, 11, 9, 7, 5 };
    enum {
      num_indices = sizeof(indices) / sizeof(indices[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    test_get_indices_at_positions(num_indices, indices, undef_idx,
                                  num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: only invalid positions
    Xt_int undef_idx = XT_INT_MIN;
    static const Xt_int indices[]
      = { 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2, 1 };
    static const int pos[]
      = { 100, 102, 107, 109, 1011, 1011, 109, 107, 105 };
    enum {
      num_indices = sizeof(indices) / sizeof(indices[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    test_get_indices_at_positions(num_indices, indices, undef_idx,
                                  num_pos, pos);
  }

  {
    // check get_positions_of_indices
    // case: some unmatched indices
    enum { num_indices = 2 };
    static const Xt_int indices[num_indices] = { 0, 2 };

    Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);

    int position;

    if (!xt_idxlist_get_position_of_index(idxvec, 1, &position))
      PUT_ERR("xt_idxlist_get_position_of_index did not return an error\n");

    if (!xt_idxlist_get_position_of_index_off(idxvec, 1, &position, 0))
      PUT_ERR("xt_idxlist_get_position_of_index_off"
              " did not return an error\n");

    if (!xt_idxlist_get_position_of_index_off(idxvec, 0, &position, 1))
      PUT_ERR("xt_idxlist_get_position_of_index_off"
              " did not return an error\n");

    static const Xt_int selection[] = { 1, 2, 3 };
    enum { num_selection = (int)(sizeof(selection) / sizeof(*selection)) };

    int positions[num_selection];
    static const int ref_positions[] = { -1, 1, -1 };

    if (xt_idxlist_get_positions_of_indices(idxvec, selection, num_selection,
                                            positions, 0) != 2)
      PUT_ERR("xt_idxlist_get_position_of_indices did not return correct num_unmatched\n");

    bool mismatch = false;
    for (size_t i = 0; i < num_selection; i++)
      mismatch |= (positions[i] != ref_positions[i]);
    if (mismatch)
      PUT_ERR("xt_idxlist_get_position_of_indices did not return correct"
              " position\n");

    xt_idxlist_delete(idxvec);
  }

  {
    // check get_bounding_box
    static const Xt_int indices[2] = {21,42};
    int num_indices = 2;

    Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);

    enum { ndim = 3 };
    static const Xt_int global_size[ndim] = { 4, 4, 4 };
    Xt_int global_start_index = 0;
    struct Xt_bounds bounds[ndim];

    xt_idxlist_get_bounding_box(idxvec, ndim, global_size,
                                global_start_index, bounds);

    bool mismatch = false;
    for (unsigned i = 0; i < ndim; ++i)
      mismatch |= (bounds[i].size != 2 || bounds[i].start != 1);
    if (mismatch)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(idxvec);
  }

  {
    // check get_bounding_box
    Xt_int indices[5] = {45,35,32,48,33};
    int num_indices = 5;

    Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);

    enum { ndim = 3 };
    static const Xt_int global_size[ndim] = { 5, 4, 3 };
    Xt_int global_start_index = 1;
    struct Xt_bounds bounds[ndim];

    xt_idxlist_get_bounding_box(idxvec, ndim, global_size,
                                global_start_index, bounds);

    static const Xt_int ref_start[3] = {2,2,1};

    bool mismatch = false;
    for (size_t i = 0; i < ndim; ++i)
      mismatch |= (bounds[i].size != 2 || bounds[i].start != ref_start[i]);
    if (mismatch)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(idxvec);
  }

  {
    // check get_bounding_box
    static const Xt_int indices[1] = {-1};
    int num_indices = 0;

    Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);

    enum { ndim = 3 };
    static const Xt_int global_size[ndim] = { 4, 4, 4 };
    Xt_int global_start_index = 0;
    struct Xt_bounds bounds[ndim];

    xt_idxlist_get_bounding_box(idxvec, ndim, global_size,
                                global_start_index, bounds);

    bool mismatch = false;
    for (size_t i = 0; i < ndim; ++i)
      mismatch |= (bounds[i].size != 0);
    if (mismatch)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(idxvec);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static Xt_idxlist
setup_idxvec(int num_idx, const Xt_int indices[])
{
  Xt_idxlist idxvec = xt_idxvec_new(indices, num_idx);
  check_idxlist(idxvec, indices, num_idx);
  return idxvec;
}

static void
test_intersection(int num_idx_a, const Xt_int index_vector_a[],
                  int num_idx_b, const Xt_int index_vector_b[],
                  int num_idx_isect, const Xt_int index_vector_isect[])
{
  Xt_idxlist idxvec_a = setup_idxvec(num_idx_a, index_vector_a),
    idxvec_b = setup_idxvec(num_idx_b, index_vector_b),
    intersection = xt_idxlist_get_intersection(idxvec_a, idxvec_b);
  check_idxlist(intersection, index_vector_isect, num_idx_isect);
  Xt_idxlist temp[3] = { idxvec_a, idxvec_b, intersection };
  for (size_t i = 0; i < 3; ++i)
    xt_idxlist_delete(temp[i]);
}

static void
test_idxvec_from_stripes(int num_stripes, const struct Xt_stripe stripes[],
                         int num_ref_indices, const Xt_int ref_indices[])
{
  Xt_idxlist idxvec = xt_idxvec_from_stripes_new(stripes, num_stripes);
  check_idxlist(idxvec, ref_indices, num_ref_indices);
  xt_idxlist_delete(idxvec);
}

static void
test_get_indices_at_positions(int num_indices, const Xt_int *indices,
                              Xt_int undef_idx,
                              int num_pos, const int *pos)
{
  Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);
  int ref_undef_count = 0;
  num_pos = num_pos < 0 ? 0 : num_pos;
  Xt_int ref_sel_idx[num_pos];
  for (size_t i = 0; i < (size_t)num_pos; i++) {
    int p = pos[i];
    if (xt_idxlist_get_index_at_position(idxvec, p, &ref_sel_idx[i]) != 0) {
      ref_sel_idx[i] = undef_idx;
      ref_undef_count++;
    }
  }
  Xt_int sel_idx[num_pos];
  int undef_count
    = xt_idxlist_get_indices_at_positions(idxvec, pos, num_pos, sel_idx,
                                          undef_idx);
  if (undef_count != ref_undef_count)
    PUT_ERR("test_idxvec.c: undef_count != ref_undef_count\n");
  bool mismatch = false;
  for (size_t i = 0; i < (size_t)num_pos; i++)
    mismatch |= (sel_idx[i] != ref_sel_idx[i]);
  if (mismatch)
    PUT_ERR("test_idxvec.c: sel_idx != ref_sel_idx\n");
  xt_idxlist_delete(idxvec);
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
