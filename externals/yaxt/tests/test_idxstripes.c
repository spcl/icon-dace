/**
 * @file test_idxstripes.c
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
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/time.h>

#include <mpi.h>

#include <yaxt.h>

#include "core/ppm_xfuncs.h"
#include "tests.h"
#include "test_idxlist_utils.h"

static void
do_tests(Xt_idxlist idxlist, const Xt_int *ref_indices, int num_indices);

static void
stripe_test_general(int num_stripes, const struct Xt_stripe *stripes,
                    int num_ref_indices, const Xt_int *ref_indices);

static void
stripe_test_general1(void);
static void
stripe_test_general2(void);
static void
stripe_test_general3(void);
static void
stripe_test_general4(void);
static void
stripe_test_general5(void);
static void
stripe_test_general6(void);

static void
stripe_test_intersection(int num_stripes_a, const struct Xt_stripe *stripes_a,
                         int num_stripes_b, const struct Xt_stripe *stripes_b,
                         int num_ref_indices, const Xt_int *ref_indices);

static void
stripe_test_asymmetric_intersection(
  int num_stripes_a, const struct Xt_stripe *stripes_a,
  int num_stripes_b, const struct Xt_stripe *stripes_b,
  int num_ref_indices_a, const Xt_int *ref_indices_a,
  int num_ref_indices_b, const Xt_int *ref_indices_b);

static void
check_idxvec_get_indices_at_positions(int num_stripes,
                                      const struct Xt_stripe *stripes,
                                      int num_pos, const int *pos);

static void
check_idxlist_stripes_pos_ext(Xt_idxlist idxlist,
                              size_t num_stripes,
                              const struct Xt_stripe stripes[num_stripes]);

static void
check_bb(int num_stripes, const struct Xt_stripe *stripes,
         unsigned ndim, const Xt_int global_size[ndim],
         const struct Xt_bounds ref_bounds[ndim],
         Xt_int global_start_index);

static void
check_pos_ext(size_t num_stripes, const struct Xt_stripe stripes[num_stripes],
              size_t num_search_stripes,
              const struct Xt_stripe search_stripes[num_search_stripes],
              size_t ref_num_ext,
              const struct Xt_pos_ext ref_pos_ext[ref_num_ext],
              int single_match_only, int ref_unmatched,
              const char *test_desc);

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  stripe_test_general1();

  stripe_test_general2();

  stripe_test_general3();

  stripe_test_general4();

  stripe_test_general5();

  stripe_test_general6();

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 1,
      num_ref_indices = 6,
    };
    static const struct Xt_stripe
      stripes_a[num_stripes_a] = {{.start = 0, .stride = 1, .nstrides = 4},
                                  {.start = 6, .stride = 1, .nstrides = 4}},
      stripes_b[num_stripes_b] = {{.start = 1, .stride = 1, .nstrides = 8}};
    static const Xt_int ref_indices[num_ref_indices] = {1,2,3, 6,7,8};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 3,
      num_stripes_b = 2,
      num_ref_indices = 9,
    };
    static const struct Xt_stripe
      stripes_a[num_stripes_a] = {{.start =  0, .stride = 1, .nstrides = 4},
                                  {.start =  6, .stride = 1, .nstrides = 4},
                                  {.start = 11, .stride = 1, .nstrides = 4}},
      stripes_b[num_stripes_b] = {{.start =  1, .stride = 1, .nstrides = 7},
                                  {.start =  9, .stride = 1, .nstrides = 5}};
    static const Xt_int ref_indices[num_ref_indices]
      = {1,2,3, 6,7, 9, 11,12,13};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 2,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe
      stripes_a[num_stripes_a] = {{.start =  0, .stride = 1, .nstrides = 3},
                                  {.start =  8, .stride = 1, .nstrides = 3}},
      stripes_b[num_stripes_b] = {{.start =  3, .stride = 1, .nstrides = 5},
                                  {.start = 11, .stride = 1, .nstrides = 3}};
    static const Xt_int *ref_indices = NULL;

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .nstrides = 10, .stride =  1}},
      stripes_b[num_stripes_b] = {{.start = 0, .nstrides =  5, .stride =  2},
                                  {.start = 9, .nstrides =  5, .stride = -2}};
    static const Xt_int ref_indices[num_ref_indices] = {0,1,2,3,4,5,6,7,8,9};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 2,
      num_ref_indices = 6,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride =  3, .nstrides =  5},
         {.start =  1, .stride =  7, .nstrides =  5}},
      stripes_b[num_stripes_b] = {{.start =  0, .stride =  2, .nstrides =  7},
                                  {.start = 24, .stride = -1, .nstrides = 10}};
    static const Xt_int ref_indices[6] = {0,6,8,12,15,22};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .stride =  1, .nstrides = 10}},
      stripes_b[num_stripes_b] = {{.start = 5, .stride =  1, .nstrides =  5},
                                  {.start = 4, .stride = -1, .nstrides =  5}};
    static const Xt_int ref_indices[10] = {0,1,2,3,4,5,6,7,8,9};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 2,
      num_ref_indices = 7,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride = 1, .nstrides = 10},
         {.start = 20, .stride = 1, .nstrides =  5}},
      stripes_b[num_stripes_b] = {{.start =  3, .stride = 1, .nstrides =  5},
                                  {.start = 17, .stride = 1, .nstrides =  5}};
    static const Xt_int ref_indices[7] = {3,4,5,6,7,20,21};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 10,
      num_stripes_b = 5,
      num_ref_indices = 6,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride = 1, .nstrides = 2},
         {.start =  3, .stride = 1, .nstrides = 2},
         {.start =  5, .stride = 1, .nstrides = 2},
         {.start =  8, .stride = 1, .nstrides = 2},
         {.start = 10, .stride = 1, .nstrides = 2},
         {.start = 14, .stride = 1, .nstrides = 2},
         {.start = 17, .stride = 1, .nstrides = 2},
         {.start = 20, .stride = 1, .nstrides = 2},
         {.start = 23, .stride = 1, .nstrides = 2},
         {.start = 25, .stride = 1, .nstrides = 2}},
      stripes_b[num_stripes_b] =  {{.start =  5, .stride = 1, .nstrides = 3},
                                   {.start =  8, .stride = 1, .nstrides = 2},
                                   {.start = 19, .stride = 1, .nstrides = 1},
                                   {.start = 20, .stride = 1, .nstrides = 2},
                                   {.start = 30, .stride = 1, .nstrides = 2}};
    static const Xt_int ref_indices[6] = {5,6,8,9,20,21};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 3,
      num_stripes_b = 1,
      num_ref_indices_a = 7,
      num_ref_indices_b = 15,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride = 1, .nstrides =  5},
         {.start =  1, .stride = 1, .nstrides =  5},
         {.start =  2, .stride = 1, .nstrides =  5}},
      stripes_b[num_stripes_b] = {{.start = -2, .stride = 1, .nstrides = 10}};
    static const Xt_int ref_indices_a[num_ref_indices_a] = {0,1,2,3,4,5,6},
      ref_indices_b[num_ref_indices_b] = {0,1,1,2,2,2,3,3,3,4,4,4,5,5,6};

    stripe_test_asymmetric_intersection(num_stripes_a, stripes_a,
                                        num_stripes_b, stripes_b,
                                        num_ref_indices_a, ref_indices_a,
                                        num_ref_indices_b, ref_indices_b);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .stride = 2, .nstrides = 5}},
      stripes_b[num_stripes_b] = {{.start = 1, .stride = 2, .nstrides = 5}};

    static const Xt_int ref_indices[1] = {0};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 3,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .stride = 5, .nstrides = 20}},
      stripes_b[num_stripes_b] = {{.start = 1, .stride = 7, .nstrides = 15}};
    static const Xt_int ref_indices[3] = {15,50,85};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test, both ranges overlap in range but have no
    // indices in common because of stride
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 34, .stride = 29, .nstrides = 12}},
      stripes_b[num_stripes_b] = {{.start = 36, .stride = 7, .nstrides = 2}};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, NULL);
  }

  { // intersection test, same as before but with negative stride
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 353, .stride = -29, .nstrides = 12}},
      stripes_b[num_stripes_b] = {{.start = 36, .stride = 7, .nstrides = 2}};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, NULL);
  }

  { // intersection test, both ranges overlap in range but have no
    // indices in common because of stride
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = -923, .stride = 1678, .nstrides = 2}},
      stripes_b[num_stripes_b] = {{.start = 844, .stride = -1779, .nstrides = 2}};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, NULL);
  }


  {
    /* test case where overlap equals start of one and end of other stripe */
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 1,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 95, .stride = -29, .nstrides = 2}},
      stripes_b[num_stripes_b] = {{.start = 81, .stride = 14, .nstrides = 2}};
    static const Xt_int ref_indices[num_ref_indices] = {95};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  {
    /* test case where overlap equals end of both stripes */
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 1,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 546, .stride = 14, .nstrides = 2}},
      stripes_b[num_stripes_b] = {{.start = 354, .stride = 206, .nstrides = 2}};
    static const Xt_int ref_indices[num_ref_indices] = {560};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

#if 328669 <= XT_INT_MAX
  {
    /* test case where overlap equals start of one and end of other
     * stripe and strides are huge */
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 1,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 50201107, .stride = 99947, .nstrides = 2}},
      stripes_b[num_stripes_b]
      = {{.start = 110701013, .stride = -60399959, .nstrides = 2}};
    static const Xt_int ref_indices[num_ref_indices] = {50301054};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }
#endif

  {
    // generate idxvec from stripes
    static const struct Xt_stripe stripes[]
      = {{.start = 4, .stride = 1, .nstrides = 1},
         {.start = 5, .stride = 1, .nstrides = 1},
         {.start = 10, .stride = -10, .nstrides = 2}};
    int num_stripes = sizeof(stripes) / sizeof(stripes[0]);
    Xt_idxlist idxvec_a = xt_idxvec_from_stripes_new(stripes, num_stripes);

    // generate idxvec containing single elemente that is also with idxvec_a
    Xt_int index_vector[] = {5};
    int num_indices = sizeof(index_vector) / sizeof(index_vector[0]);
    Xt_idxlist idxvec_b = xt_idxvec_new(index_vector, num_indices);

    // intersect both index vectors
    Xt_idxlist intersection =
      xt_idxlist_get_intersection(idxvec_a, idxvec_b);

    // test result
    if (xt_idxlist_get_num_indices(intersection) != 1) {
      PUT_ERR("wrong number of indices in intersection\n");
    } else {
      Xt_int idx;
      xt_idxlist_get_index_at_position(intersection, 0, &idx);
      if (idx != index_vector[0])
        PUT_ERR("wrong index in intersection\n");
    }
    xt_idxlist_delete(intersection);
    xt_idxlist_delete(idxvec_b);
    xt_idxlist_delete(idxvec_a);
  }

  {
    static const Xt_int index_vector[] = {
      3375, 3376,
      3379, 3380, 3381,
      3387, 3388, 3389, 3390, 3391, 3392, 3393,
      3421, 3422, 3423, 3424, 3425, 3426, 3427,
      3444,
      3458, 3459,
      3461, 3462, 3463, 3464, 3465, 3466, 3467, 3468, 3469, 3470,
        3471, 3472, 3473, 3474, 3475, 3476, 3477, 3478, 3479, 3480,
      3529,
      3606, 3607, 3608,
      3611, 3612, 3613, 3614,
      3617,
      3620, 3621, 3622, 3623, 3624, 3625, 3626, 3627, 3628, 3629,
        3630, 3631,
      3684, 3685, 3686, 3687, 3688, 3689, 3690, 3691, 3692, 3693,
        3694, 3695, 3696, 3697, 3698, 3699, 3700, 3701, 3702, 3703,
        3704, 3705, 3706, 3707, 3708, 3709,
      3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722,
        3723, 3724, 3725, 3726, 3727, 3728, 3729, 3730, 3731,
      3741, 3742,
      3931, 3932,
      3374,
      3382,
      3385,
      3394,
      3404,
      3408,
      3412,
      3440,
      3443,
      3457,
      3481,
      3483,
      3527,
      3619,
      3735,
      3743,
      3925,
      3930,
      3377, 3378,
      3383, 3384,
      3386,
      3395,
      3397, 3398,
      3400,
      3402, 3403,
      3407,
      3409, 3410,
      3413,
      3420,
      3441, 3442,
      3445,
      3448, 3449,
      3451,
      3460,
      3482,
      3519, 3520,
      3526,
      3528,
      3530,
      3592, 3593,
      3595, 3596, 3597,
      3609, 3610,
      3615, 3616,
      3618,
      3644,
      3710, 3711, 3712,
      3732, 3733,
      3736, 3737,
      3748, 3749,
      3753, 3754,
      3759, 3760,
      3766, 3767,
      3919, 3920,
      3924,
      3926,
      3933, 3934,
      2589,
      2602,
      2680,
      3326,
      3340, 3341,
      3396,
      3401,
      3411,
      3414,
      3418,
      3446, 3447,
      3450,
      3515,
      3521,
      3525,
      3582,
      3590, 3591,
      3594,
      3642,
      3734,
      3738,
      3747,
      3750,
      3761,
      3765,
      3865,
      3918,
      3923,
      3935
    };
    enum { num_indices = sizeof(index_vector) / sizeof(index_vector[0]) };
    Xt_idxlist idxlist = xt_idxvec_new(index_vector, num_indices);

    static const struct Xt_stripe stripes[] = {
      {.start = 3326, .stride = 14, .nstrides = 2 },
      {.start = 3341, .stride = 33, .nstrides = 1 },
      {.start = 3374, .stride = 1 , .nstrides = 25},
      {.start = 3400, .stride = 1 , .nstrides = 5 },
      {.start = 3407, .stride = 1 , .nstrides = 8 },
      {.start = 3418, .stride = 2 , .nstrides = 1 },
      {.start = 3420, .stride = 1 , .nstrides = 8 },
      {.start = 3440, .stride = 1 , .nstrides = 12},
      {.start = 3457, .stride = 1 , .nstrides = 27},
      {.start = 3515, .stride = 4 , .nstrides = 1 },
      {.start = 3519, .stride = 1 , .nstrides = 3 },
      {.start = 3525, .stride = 1 , .nstrides = 6 },
      {.start = 3582, .stride = 8 , .nstrides = 1 },
      {.start = 3590, .stride = 1 , .nstrides = 8 },
      {.start = 3606, .stride = 1 , .nstrides = 26},
      {.start = 3642, .stride = 2 , .nstrides = 2 },
      {.start = 3684, .stride = 1 , .nstrides = 55},
      {.start = 3741, .stride = 1 , .nstrides = 3 },
      {.start = 3747, .stride = 1 , .nstrides = 4 },
      {.start = 3753, .stride = 1 , .nstrides = 2 },
      {.start = 3759, .stride = 1 , .nstrides = 3 },
      {.start = 3765, .stride = 1 , .nstrides = 3 },
      {.start = 3865, .stride = 53, .nstrides = 1 },
      {.start = 3918, .stride = 1 , .nstrides = 3 },
      {.start = 3923, .stride = 1 , .nstrides = 4 },
      {.start = 3930, .stride = 1 , .nstrides = 6 }
    };
    enum { num_stripes = sizeof(stripes) / sizeof(stripes[0]) };

    check_idxlist_stripes_pos_ext(idxlist, (size_t)num_stripes, stripes);

    xt_idxlist_delete(idxlist);
  }

#if 328669 <= XT_INT_MAX
  {
    static const Xt_int index_vector[] = {
      328669, 30608, 38403
    };
    enum { num_indices = sizeof (index_vector) / sizeof (index_vector[0]) };
    Xt_idxlist idxlist = xt_idxvec_new(index_vector, num_indices);

    static const struct Xt_stripe stripes[] = {
      { .start = 30608, .stride = 7795, .nstrides = 2 },
    };
    enum { num_stripes = sizeof(stripes) / sizeof(stripes[0]) };

    check_idxlist_stripes_pos_ext(idxlist, (size_t)num_stripes, stripes);

    xt_idxlist_delete(idxlist);
  }

  {
    static const Xt_int index_vector[] = {
      679605, 726349, 726346 };
    enum { num_indices = sizeof (index_vector) / sizeof (index_vector[0]) };
    Xt_idxlist idxlist = xt_idxvec_new(index_vector, num_indices);

    static const struct Xt_stripe stripes[] = {
      {.start = 679605, .stride = 46741, .nstrides = 2  },
    };
    enum { num_stripes = sizeof (stripes) / sizeof (stripes[0]) };

    check_idxlist_stripes_pos_ext(idxlist, (size_t)num_stripes, stripes);

    xt_idxlist_delete(idxlist);
  }
#endif

  {
    enum {
      NUM_ITERATIONS=128,
      MAX_NUM_INDICES=1024,
      MAX_INDEX=1024,
    };

    /* demonstrates problem of invalid intersection
     * get_stripe_intersection({916,-287,2}, {603,300,2})
     * gives an intersection with nstrides == 1 even though we'd need
     * nstrides == 0
     */
    enum {
       num_seeds_max = 3
    };
    unsigned int seeds[num_seeds_max];
    seeds[0] = 6544668;
    seeds[1] = 1561743077;
    const char *renv = getenv("YAXT_FULLY_RANDOM_TESTS");
    bool fully_random_test
      = (renv && renv[0] && ((!renv[1] && (renv[0] == '1' || renv[0] == 'y'
                                           || renv[0] == 'Y'))
                             || !strcasecmp(renv, "yes")));
    size_t num_seeds = num_seeds_max - (size_t)!fully_random_test;
    if (fully_random_test)
    {
      struct timeval tv;
      int rc = gettimeofday(&tv, NULL);
      if (rc) {
        perror("cannot query time of day");
        xt_mpi_call(MPI_Abort(MPI_COMM_WORLD, 1), MPI_COMM_WORLD);
      }
      seeds[num_seeds_max-1]
        = (unsigned)((unsigned)tv.tv_sec ^ (unsigned)tv.tv_usec);
    }
    Xt_int *restrict indices = NULL;
    for (size_t seed_idx = 0; seed_idx < num_seeds; ++seed_idx) {
      srand(seeds[seed_idx]);
      if (seed_idx > 0)
        fprintf(stderr, "seed=%u\n", seeds[seed_idx]);
      for (int iteration = 0; iteration < NUM_ITERATIONS; ++iteration) {

        int num_indices = rand()%MAX_NUM_INDICES;

        indices
          = xrealloc(indices, (size_t)num_indices * sizeof(*indices));

        for (int i = 0; i < num_indices; ++i)
          indices[i] = (Xt_int)(rand()%(2*MAX_INDEX)-MAX_INDEX);

        Xt_idxlist idxlist = xt_idxvec_new(indices, num_indices);

        struct Xt_stripe * stripes;
        int num_stripes;

        xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);

        check_idxlist_stripes_pos_ext(idxlist, (size_t)num_stripes, stripes);

        xt_idxlist_delete(idxlist);
        free(stripes);
      }
    }
    free(indices);
  }

  {
    static const Xt_int index_vector[] = {
      178, 179, 180, 181, 182, 183, 184,
      186, 187, 188, 189, 190,
      194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212,
      217,
      223,
      426,
      428, 429, 430,
      434, 435, 436, 437, 438, 439, 440,
      442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455,
        456, 457, 458,
      670, 671, 672, 673, 674, 675, 676, 677,
      682,
      684, 685, 686, 687, 688, 689, 690,
      692,
      695,
      703, 704, 705, 706, 707,
      894, 895, 896, 897, 898, 899, 900, 901,
      906, 907, 908,
      913,
      915,
      921, 922, 923, 924, 925, 926, 927,
      1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103,
      1107, 1108, 1109, 1110, 1111,
      1113, 1114,
      1119, 1120, 1121,
      2095, 2096, 2097, 2098,
      2100,
      2102, 2103, 2104, 2105,
      2107, 2108, 2109, 2110,
      2112,
      2118,
      2120, 2121, 2122, 2123, 2124, 2125,
      2127, 2128, 2129, 2130,
      2134,
      2140, 2141, 2142, 2143,
      2145,
      2148, 2149,
      2151, 2152, 2153, 2154, 2155, 2156,
      683,
      691,
      903,
      914,
      1105,
      1115,
      2099,
      2106,
      2111,
      2115,
      2126,
      2132,
      2139,
      2144,
      2147,
      2150,
      2305,
      427,
      465, 466,
      678,
      693,
      902,
      909,
      1104,
      1112,
      2101,
      2113, 2114,
      2116, 2117,
      2119,
      2131,
      2136,
      2138,
      2146,
      2297,
      2302,
      2304,
      2307
    };
    enum { num_indices = sizeof(index_vector) / sizeof(index_vector[0]) };
    Xt_idxlist idxlist = xt_idxvec_new(index_vector, num_indices);

    static const struct  Xt_stripe stripes[] = {
      {.start = 670, .stride = 1, .nstrides = 9},
      {.start = 682, .stride = 1, .nstrides = 12},
      {.start = 695, .stride = 8, .nstrides = 1},
      {.start = 703, .stride = 1, .nstrides = 5},
      {.start = 894, .stride = 1, .nstrides = 10},
      {.start = 906, .stride = 1, .nstrides = 4},
      {.start = 913, .stride = 1, .nstrides = 3},
      {.start = 921, .stride = 1, .nstrides = 7}
    };
    enum { num_stripes = sizeof(stripes) / sizeof(stripes[0]) };

    check_idxlist_stripes_pos_ext(idxlist, num_stripes, stripes);

    xt_idxlist_delete(idxlist);
  }

  {
    static const Xt_int index_vector[] = {
2055, 2056, 2060, 2193, 2199, 2203, 2211, 2212, 2278, 2281, 2311, 2312, 2316,
2317, 2322, 2332, 2447, 2448, 2452, 2585, 2591, 2595, 2603, 2604, 2670, 2673,
2703, 2704, 2708, 2709, 2714, 2724, 2839, 2840, 2844, 2977, 2983, 2987, 2995,
2996, 3062, 3065, 3095, 3096, 3100, 3101, 3106, 3116, 3231, 3232, 3236, 3369,
3375, 3379, 3387, 3388, 3454, 3457, 3487, 3488, 3492, 3493, 3498, 3508, 3623,
3624, 3628, 3761, 3767, 3771, 3779, 3780, 3846, 3849, 3879, 3880, 3884, 3885,
3890, 3900, 3997, 4001, 4002, 4053, 4057, 4084, 4085, 4092, 4102, 4188, 4192,
4201, 4373, 4377, 4378, 4429, 4433, 4460, 4461, 4468, 4478, 4564, 4568, 4577,
4749, 4753, 4754, 4805, 4809, 4836, 4837, 4844, 4854, 4945, 4953, 5125, 5129,
5130, 5181, 5185, 5212, 5213, 5220, 5230, 5321, 5329, 5501, 5505, 5506, 5557,
5561, 5588, 5589, 5596, 5606, 5697, 5705, 162, 163, 166, 168, 171, 172, 173,
177, 181, 362, 363, 367, 369, 375, 378, 382, 383, 386, 570, 571, 574, 576, 579,
580, 581, 585, 589, 758, 759, 763, 765, 769, 774, 775, 778, 962, 963, 966, 968,
971, 972, 973, 977, 981, 1150, 1151, 1155, 1157, 1161, 1166, 1167, 1170, 1354,
1355, 1358, 1360, 1363, 1364, 1365, 1369, 1373, 1542, 1543, 1547, 1549, 1553,
1558, 1559, 1562, 1746, 1747, 1750, 1752, 1755, 1756, 1757, 1761, 1918, 1919,
1923, 1925, 1929, 1934, 1935, 1938, 1988, 1989, 2024, 2025, 2032, 2033, 2036,
2038, 2039, 2048, 2049, 2053, 2054, 2057, 2058, 2061, 2076, 2077, 2091, 2092,
2093, 2095, 2097, 2126, 2127, 2144, 2145, 2149, 2150, 2156, 2198, 2204, 2205,
2207, 2245, 2253, 2254, 2256, 2268, 2269, 2277, 2279, 2280, 2283, 2287, 2298,
2299, 2307, 2308, 2309, 2310, 2333, 2334, 2380, 2381, 2416, 2417, 2424, 2425,
2428, 2430, 2431, 2440, 2441, 2445, 2446, 2449, 2450, 2453, 2468, 2469, 2483,
2484, 2485, 2487, 2489, 2518, 2519, 2536, 2537, 2541, 2542, 2548, 2590, 2596,
2597, 2599, 2637, 2645, 2646, 2648, 2660, 2661, 2669, 2671, 2672, 2675, 2679,
2690, 2691, 2699, 2700, 2701, 2702, 2725, 2726, 2772, 2773, 2808, 2809, 2816,
2817, 2820, 2822, 2823, 2832, 2833, 2837, 2838, 2841, 2842, 2845, 2860, 2861,
2875, 2876, 2877, 2879, 2881, 2910, 2911, 2928, 2929, 2933, 2934, 2940, 2982,
2988, 2989, 2991, 3029, 3037, 3038, 3040, 3052, 3053, 3061, 3063, 3064, 3067,
3071, 3082, 3083, 3091, 3092, 3093, 3094, 3117, 3118, 3164, 3165, 3200, 3201,
3208, 3209, 3212, 3214, 3215, 3224, 3225, 3229, 3230, 3233, 3234, 3237, 3252,
3253, 3267, 3268, 3269, 3271, 3273, 3302, 3303, 3320, 3321, 3325, 3326, 3332,
3374, 3380, 3381, 3383, 3421, 3429, 3430, 3432, 3444, 3445, 3453, 3455, 3456,
3459, 3463, 3474, 3475, 3483, 3484, 3485, 3486, 3509, 3510, 3556, 3557, 3592,
3593, 3600, 3601, 3604, 3606, 3607, 3616, 3617, 3621, 3622, 3625, 3626, 3629,
3644, 3645, 3659, 3660, 3661, 3663, 3665, 3694, 3695, 3712, 3713, 3717, 3718,
3724, 3766, 3772, 3773, 3775, 3813, 3821, 3822, 3824, 3836, 3837, 3845, 3847,
3848, 3851, 3855, 3866, 3867, 3875, 3876, 3877, 3878, 3901, 3902, 3948, 3949,
3984, 3985, 3992, 3993, 3996, 3998, 3999, 4008, 4009, 4013, 4014, 4017, 4018,
4021, 4036, 4037, 4051, 4052, 4054, 4055, 4058, 4090, 4091, 4093, 4108, 4109,
4112, 4113, 4114, 4158, 4164, 4165, 4193, 4199, 4200, 4212, 4213, 4222, 4223,
4225, 4227, 4231, 4242, 4243, 4250, 4251, 4271, 4272, 4274, 4324, 4325, 4360,
4361, 4368, 4369, 4372, 4374, 4375, 4384, 4385, 4389, 4390, 4393, 4394, 4397,
4412, 4413, 4427, 4428, 4430, 4431, 4434, 4466, 4467, 4469, 4484, 4485, 4488,
4489, 4490, 4534, 4540, 4541, 4569, 4575, 4576, 4588, 4589, 4598, 4599, 4601,
4603, 4607, 4618, 4619, 4626, 4627, 4647, 4648, 4650, 4700, 4701, 4736, 4737,
4744, 4745, 4748, 4750, 4751, 4760, 4761, 4765, 4766, 4769, 4770, 4773, 4788,
4789, 4803, 4804, 4806, 4807, 4810, 4842, 4843, 4845, 4860, 4861, 4864, 4865,
4866, 4910, 4916, 4917, 4951, 4952, 4964, 4965, 4974, 4975, 4977, 4979, 4983,
4994, 4995, 5002, 5003, 5023, 5024, 5026, 5076, 5077, 5112, 5113, 5120, 5121,
5124, 5126, 5127, 5136, 5137, 5141, 5142, 5145, 5146, 5149, 5164, 5165, 5179,
5180, 5182, 5183, 5186, 5218, 5219, 5221, 5236, 5237, 5240, 5241, 5242, 5286,
5292, 5293, 5327, 5328, 5340, 5341, 5350, 5351, 5353, 5355, 5359, 5370, 5371,
5378, 5379, 5399, 5400, 5402, 5452, 5453, 5488, 5489, 5496, 5497, 5500, 5502,
5503, 5512, 5513, 5517, 5518, 5521, 5522, 5525, 5540, 5541, 5555, 5556, 5558,
5559, 5562, 5594, 5595, 5597, 5612, 5613, 5616, 5617, 5618, 5662, 5668, 5669,
5703, 5704, 5716, 5717, 5726, 5727, 5729, 5731, 5735, 5746, 5747, 5754, 5755,
5775, 5776, 5778, 5958, 5959, 5962, 5964, 5967, 5968, 5971, 5973, 6154, 6155,
6159, 6161, 6167, 6170, 6172, 6173, 6350, 6351, 6354, 6356, 6359, 6360, 6363,
6530, 6531, 6535, 6537, 6543, 6546, 6548, 6549, 6726, 6727, 6730, 6732, 6735,
6736, 6739, 6906, 6907, 6911, 6913, 6919, 6922, 6924, 6925, 7102, 7103, 7106,
7108, 7111, 7112, 7115, 7282, 7283, 7287, 7289, 7295, 7298, 7300, 7301, 7478,
7479, 7482, 7484, 7487, 7488, 7491, 7646, 7647, 7651, 7653, 7657, 7660, 7661,
130, 161, 169, 170, 336, 361, 366, 384, 538, 569, 577, 578, 736, 757, 762, 776,
930, 961, 969, 970, 1128, 1149, 1154, 1168, 1322, 1353, 1361, 1362, 1520, 1541,
1546, 1560, 1714, 1745, 1753, 1754, 1896, 1917, 1922, 1936, 1985, 2019, 2031,
2035, 2040, 2044, 2052, 2059, 2062, 2071, 2087, 2090, 2094, 2140, 2148, 2153,
2157, 2206, 2257, 2263, 2267, 2284, 2288, 2293, 2295, 2305, 2306, 2377, 2411,
2423, 2427, 2432, 2436, 2444, 2451, 2454, 2463, 2479, 2482, 2486, 2532, 2540,
2545, 2549, 2598, 2649, 2655, 2659, 2676, 2680, 2685, 2687, 2697, 2698, 2769,
2803, 2815, 2819, 2824, 2828, 2836, 2843, 2846, 2855, 2871, 2874, 2878, 2924,
2932, 2937, 2941, 2990, 3041, 3047, 3051, 3068, 3072, 3077, 3079, 3089, 3090,
3161, 3195, 3207, 3211, 3216, 3220, 3228, 3235, 3238, 3247, 3263, 3266, 3270,
3316, 3324, 3329, 3333, 3382, 3433, 3439, 3443, 3460, 3464, 3469, 3471, 3481,
3482, 3553, 3587, 3599, 3603, 3608, 3612, 3620, 3627, 3630, 3639, 3655, 3658,
3662, 3708, 3716, 3721, 3725, 3774, 3825, 3831, 3835, 3852, 3856, 3861, 3863,
3873, 3874, 3945, 3979, 3991, 3995, 4000, 4004, 4012, 4019, 4022, 4031, 4033,
4047, 4050, 4104, 4106, 4115, 4207, 4221, 4228, 4232, 4237, 4249, 4252, 4321,
4355, 4367, 4371, 4376, 4380, 4388, 4395, 4398, 4407, 4409, 4423, 4426, 4480,
4482, 4491, 4583, 4597, 4604, 4608, 4613, 4625, 4628, 4697, 4731, 4743, 4747,
4752, 4756, 4764, 4771, 4774, 4783, 4785, 4799, 4802, 4856, 4858, 4867, 4959,
4973, 4980, 4984, 4989, 5001, 5004, 5073, 5107, 5119, 5123, 5128, 5132, 5140,
5147, 5150, 5159, 5161, 5175, 5178, 5232, 5234, 5243, 5335, 5349, 5356, 5360,
5365, 5377, 5380, 5449, 5483, 5495, 5499, 5504, 5508, 5516, 5523, 5526, 5535,
5537, 5551, 5554, 5608, 5610, 5619, 5711, 5725, 5732, 5736, 5741, 5753, 5756,
5930, 5957, 5965, 5966, 6128, 6153, 6158, 6174, 6322, 6349, 6357, 6358, 6504,
6529, 6534, 6550, 6698, 6725, 6733, 6734, 6880, 6905, 6910, 6926, 7074, 7101,
7109, 7110, 7256, 7281, 7286, 7302, 7450, 7477, 7485, 7486, 7624, 7645, 7650,
7662
    };
    enum { num_indices = sizeof(index_vector) / sizeof(index_vector[0]) };
    Xt_idxlist idxlist = xt_idxvec_new(index_vector, num_indices);

    static const struct  Xt_stripe stripes[] = {
{.start = 173, .stride = 408, .nstrides = 2},
{.start = 973, .stride = 392, .nstrides = 3},
{.start = 1985, .stride = 4, .nstrides = 2},
{.start = 2044, .stride = 4, .nstrides = 2},
{.start = 2049, .stride = 3, .nstrides = 1},
{.start = 2052, .stride = 1, .nstrides = 9},
{.start = 2062, .stride = 131, .nstrides = 2},
{.start = 2198, .stride = 1, .nstrides = 2},
{.start = 2203, .stride = 1, .nstrides = 5},
{.start = 2211, .stride = 1, .nstrides = 2},
{.start = 2263, .stride = 4, .nstrides = 1},
{.start = 2267, .stride = 1, .nstrides = 3},
{.start = 2277, .stride = 1, .nstrides = 5},
{.start = 2283, .stride = 1, .nstrides = 2},
{.start = 2287, .stride = 1, .nstrides = 2},
{.start = 2293, .stride = 2, .nstrides = 2},
{.start = 2298, .stride = 1, .nstrides = 2},
{.start = 2305, .stride = 1, .nstrides = 8},
{.start = 2316, .stride = 1, .nstrides = 2},
{.start = 2322, .stride = 10, .nstrides = 1},
{.start = 2332, .stride = 1, .nstrides = 3},
{.start = 2377, .stride = 4, .nstrides = 2},
{.start = 2436, .stride = 4, .nstrides = 2},
{.start = 2441, .stride = 3, .nstrides = 1},
{.start = 2444, .stride = 1, .nstrides = 9},
{.start = 2454, .stride = 131, .nstrides = 2},
{.start = 2590, .stride = 1, .nstrides = 2},
{.start = 2595, .stride = 1, .nstrides = 5},
{.start = 2603, .stride = 1, .nstrides = 2},
{.start = 2655, .stride = 4, .nstrides = 1},
{.start = 2659, .stride = 1, .nstrides = 3},
{.start = 2669, .stride = 1, .nstrides = 5},
{.start = 2675, .stride = 1, .nstrides = 2},
{.start = 2679, .stride = 1, .nstrides = 2},
{.start = 2685, .stride = 2, .nstrides = 2},
{.start = 2690, .stride = 1, .nstrides = 2},
{.start = 2697, .stride = 1, .nstrides = 8},
{.start = 2708, .stride = 1, .nstrides = 2},
{.start = 2714, .stride = 10, .nstrides = 1},
{.start = 2724, .stride = 1, .nstrides = 3},
{.start = 2769, .stride = 4, .nstrides = 2},
{.start = 2828, .stride = 4, .nstrides = 2},
{.start = 2833, .stride = 3, .nstrides = 1},
{.start = 2836, .stride = 1, .nstrides = 9},
{.start = 2846, .stride = 131, .nstrides = 2},
{.start = 2982, .stride = 1, .nstrides = 2},
{.start = 2987, .stride = 1, .nstrides = 5},
{.start = 2995, .stride = 1, .nstrides = 2},
{.start = 3047, .stride = 4, .nstrides = 1},
{.start = 3051, .stride = 1, .nstrides = 3},
{.start = 3061, .stride = 1, .nstrides = 5},
{.start = 3067, .stride = 1, .nstrides = 2},
{.start = 3071, .stride = 1, .nstrides = 2},
{.start = 3077, .stride = 2, .nstrides = 2},
{.start = 3082, .stride = 1, .nstrides = 2},
{.start = 3089, .stride = 1, .nstrides = 8},
{.start = 3100, .stride = 1, .nstrides = 2},
{.start = 3106, .stride = 10, .nstrides = 1},
{.start = 3116, .stride = 1, .nstrides = 3},
{.start = 3161, .stride = 4, .nstrides = 2},
{.start = 3220, .stride = 4, .nstrides = 2},
{.start = 3225, .stride = 3, .nstrides = 1},
{.start = 3228, .stride = 1, .nstrides = 9},
{.start = 3238, .stride = 131, .nstrides = 2},
{.start = 3374, .stride = 1, .nstrides = 2},
{.start = 3379, .stride = 1, .nstrides = 5},
{.start = 3387, .stride = 1, .nstrides = 2},
{.start = 3439, .stride = 4, .nstrides = 1},
{.start = 3443, .stride = 1, .nstrides = 3},
{.start = 3453, .stride = 1, .nstrides = 5},
{.start = 3459, .stride = 1, .nstrides = 2},
{.start = 3463, .stride = 1, .nstrides = 2},
{.start = 3469, .stride = 2, .nstrides = 2},
{.start = 3474, .stride = 1, .nstrides = 2},
{.start = 3481, .stride = 1, .nstrides = 8},
{.start = 3492, .stride = 1, .nstrides = 2},
{.start = 3498, .stride = 10, .nstrides = 1},
{.start = 3508, .stride = 1, .nstrides = 3},
{.start = 3553, .stride = 4, .nstrides = 2},
{.start = 3612, .stride = 4, .nstrides = 2},
{.start = 3617, .stride = 3, .nstrides = 1},
{.start = 3620, .stride = 1, .nstrides = 9},
{.start = 3630, .stride = 131, .nstrides = 2},
{.start = 3766, .stride = 1, .nstrides = 2},
{.start = 3771, .stride = 1, .nstrides = 5},
{.start = 3779, .stride = 1, .nstrides = 2},
{.start = 3831, .stride = 4, .nstrides = 1},
{.start = 3835, .stride = 1, .nstrides = 3},
{.start = 3845, .stride = 1, .nstrides = 5},
{.start = 3851, .stride = 1, .nstrides = 2},
{.start = 3855, .stride = 1, .nstrides = 2},
{.start = 3861, .stride = 2, .nstrides = 2},
{.start = 3866, .stride = 1, .nstrides = 2},
{.start = 3873, .stride = 1, .nstrides = 8},
{.start = 3884, .stride = 1, .nstrides = 2},
{.start = 3890, .stride = 10, .nstrides = 1},
{.start = 3900, .stride = 1, .nstrides = 3},
{.start = 3945, .stride = 3, .nstrides = 2},
{.start = 3979, .stride = 5, .nstrides = 2},
{.start = 3985, .stride = 6, .nstrides = 1},
{.start = 3991, .stride = 1, .nstrides = 3},
{.start = 3995, .stride = 2, .nstrides = 1},
{.start = 3997, .stride = 1, .nstrides = 6},
{.start = 4031, .stride = 2, .nstrides = 2},
{.start = 4036, .stride = 1, .nstrides = 2},
{.start = 4047, .stride = 3, .nstrides = 1},
{.start = 4050, .stride = 1, .nstrides = 6},
{.start = 4057, .stride = 1, .nstrides = 2},
{.start = 4084, .stride = 1, .nstrides = 2},
{.start = 4090, .stride = 1, .nstrides = 4},
{.start = 4102, .stride = 2, .nstrides = 4},
{.start = 4109, .stride = 3, .nstrides = 1},
{.start = 4112, .stride = 1, .nstrides = 4},
{.start = 4188, .stride = 4, .nstrides = 2},
{.start = 4193, .stride = 6, .nstrides = 1},
{.start = 4199, .stride = 1, .nstrides = 3},
{.start = 4321, .stride = 3, .nstrides = 2},
{.start = 4355, .stride = 5, .nstrides = 2},
{.start = 4361, .stride = 6, .nstrides = 1},
{.start = 4367, .stride = 1, .nstrides = 3},
{.start = 4371, .stride = 2, .nstrides = 1},
{.start = 4373, .stride = 1, .nstrides = 6},
{.start = 4407, .stride = 2, .nstrides = 2},
{.start = 4412, .stride = 1, .nstrides = 2},
{.start = 4423, .stride = 3, .nstrides = 1},
{.start = 4426, .stride = 1, .nstrides = 6},
{.start = 4433, .stride = 1, .nstrides = 2},
{.start = 4460, .stride = 1, .nstrides = 2},
{.start = 4466, .stride = 1, .nstrides = 4},
{.start = 4478, .stride = 2, .nstrides = 4},
{.start = 4485, .stride = 3, .nstrides = 1},
{.start = 4488, .stride = 1, .nstrides = 4},
{.start = 4564, .stride = 4, .nstrides = 2},
{.start = 4569, .stride = 6, .nstrides = 1},
{.start = 4575, .stride = 1, .nstrides = 3},
{.start = 4697, .stride = 3, .nstrides = 2},
{.start = 4731, .stride = 5, .nstrides = 2},
{.start = 4737, .stride = 6, .nstrides = 1},
{.start = 4743, .stride = 1, .nstrides = 3},
{.start = 4747, .stride = 2, .nstrides = 1},
{.start = 4749, .stride = 1, .nstrides = 6},
{.start = 4783, .stride = 2, .nstrides = 2},
{.start = 4788, .stride = 1, .nstrides = 2},
{.start = 4799, .stride = 3, .nstrides = 1},
{.start = 4802, .stride = 1, .nstrides = 6},
{.start = 4809, .stride = 1, .nstrides = 2},
{.start = 4836, .stride = 1, .nstrides = 2},
{.start = 4842, .stride = 1, .nstrides = 4},
{.start = 4854, .stride = 2, .nstrides = 4},
{.start = 4861, .stride = 3, .nstrides = 1},
{.start = 4864, .stride = 1, .nstrides = 4},
{.start = 4945, .stride = 6, .nstrides = 1},
{.start = 4951, .stride = 1, .nstrides = 3},
{.start = 5107, .stride = 5, .nstrides = 2},
{.start = 5113, .stride = 6, .nstrides = 1},
{.start = 5119, .stride = 1, .nstrides = 3},
{.start = 5123, .stride = 2, .nstrides = 1},
{.start = 5125, .stride = 1, .nstrides = 6},
{.start = 5159, .stride = 2, .nstrides = 2},
{.start = 5164, .stride = 1, .nstrides = 2},
{.start = 5175, .stride = 3, .nstrides = 1},
{.start = 5178, .stride = 1, .nstrides = 6},
{.start = 5185, .stride = 1, .nstrides = 2},
{.start = 5212, .stride = 1, .nstrides = 2},
{.start = 5218, .stride = 1, .nstrides = 4},
{.start = 5230, .stride = 2, .nstrides = 4},
{.start = 5237, .stride = 3, .nstrides = 1},
{.start = 5240, .stride = 1, .nstrides = 4},
{.start = 5321, .stride = 6, .nstrides = 1},
{.start = 5327, .stride = 1, .nstrides = 3},
{.start = 5483, .stride = 5, .nstrides = 2},
{.start = 5489, .stride = 6, .nstrides = 1},
{.start = 5495, .stride = 1, .nstrides = 3},
{.start = 5499, .stride = 2, .nstrides = 1},
{.start = 5501, .stride = 1, .nstrides = 6},
{.start = 5535, .stride = 2, .nstrides = 2},
{.start = 5540, .stride = 1, .nstrides = 2},
{.start = 5551, .stride = 3, .nstrides = 1},
{.start = 5554, .stride = 1, .nstrides = 6},
{.start = 5561, .stride = 1, .nstrides = 2},
{.start = 5588, .stride = 1, .nstrides = 2},
{.start = 5594, .stride = 1, .nstrides = 4},
{.start = 5606, .stride = 2, .nstrides = 4},
{.start = 5613, .stride = 3, .nstrides = 1},
{.start = 5616, .stride = 1, .nstrides = 4},
{.start = 5697, .stride = 6, .nstrides = 1},
{.start = 5703, .stride = 1, .nstrides = 3}
    };
    enum { num_stripes = sizeof(stripes) / sizeof(stripes[0]) };

    check_idxlist_stripes_pos_ext(idxlist, num_stripes, stripes);

    xt_idxlist_delete(idxlist);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: mixed valid and invalid positions
    static const struct Xt_stripe stripes[]
      = {{.start=0, .nstrides=5, .stride=1},
         {.start=10, .nstrides=5, .stride=1},
         {.start=20, .nstrides=5, .stride=-1}};
    static const int pos[] = {0,2,7,9,11,100,11,200,9,300,18,400,5};
    enum {
      num_pos = sizeof(pos) / sizeof(pos[0]),
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: only valid positions
    static const struct Xt_stripe stripes[]
      = {{.start= 0, .nstrides=3, .stride= 1},
         {.start=10, .nstrides=2, .stride= 1},
         {.start=20, .nstrides=6, .stride=-1},
         {.start=30, .nstrides=7, .stride=-1}};
    static const int pos[] = {-1,0,1,2,3,4,23,5,6,7,8,9,10,11,12, 0,2,100,2000};
    enum {
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: complete permutation
    static const struct Xt_stripe stripes[]
      = {{.start= 0, .nstrides=3, .stride= 1},
         {.start=10, .nstrides=2, .stride= 1},
         {.start=20, .nstrides=6, .stride=-1},
         {.start=30, .nstrides=7, .stride=-1}};
    static const int pos[] = {4,7,2,5,9,0,10,6,11,8,12,1,3};
    enum {
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: only invalid positions
    static const struct Xt_stripe stripes[]
      = {{.start=0, .nstrides=5, .stride=1},
         {.start=10, .nstrides=5, .stride=1},
         {.start=20, .nstrides=5, .stride=-1}};
    static const int pos[] = {-10,200,700,90,90,18,141};
    enum {
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  { // test with overlapping stripes
    enum {
      num_stripes = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 0, .stride = 1, .nstrides = 5},
         {.start = 1, .stride = 1, .nstrides = 5}};
    static const Xt_int ref_indices[num_ref_indices] = {0,1,2,3,4, 1,2,3,4,5};

    stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
  }

  { // check get_bounding_box
    enum {
      num_stripes = 0,
      num_ref_indices = 10,
      ndim = 3,
    };
    static const struct Xt_stripe stripes[1]
      = {{.start = -1, .stride = -1, .nstrides = -1}};
    static const Xt_int global_size[ndim] = { 4, 4, 4 };
    static const struct Xt_bounds ref_bounds[ndim]
      = { { .start = -1, .size = 0 },
          { .start = -1, .size = 0 },
          { .start = -1, .size = 0 } };
    Xt_int global_start_index = 0;

    check_bb(num_stripes, stripes,
             ndim, global_size, ref_bounds, global_start_index);
  }

  { // check get_bounding_box
    enum {
      num_stripes = 3,
      ndim = 3,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 47, .stride = -12, .nstrides = 2},
         {.start = 32, .stride =  12, .nstrides = 2},
         {.start = 36, .stride =  12, .nstrides = 2}};
    static const Xt_int global_size[ndim] = { 5, 4, 3 };
    static const struct Xt_bounds ref_bounds[ndim]
      = { { .start = 2, .size = 2 },
          { .start = 2, .size = 2 },
          { .start = 1, .size = 2 } };
    Xt_int global_start_index = 1;

    check_bb(num_stripes, stripes,
             ndim, global_size, ref_bounds, global_start_index);
  }

  {
    enum {
      num_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 0,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10}},
      search_stripe = {.start = 10, .stride = -1, .nstrides = 5 };

    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
      = {{.start = 9, .size = -5}};

    check_pos_ext(num_stripes, stripes,
                  1, &search_stripe, num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "simple inverted stripe");
  }

  {
    enum {
      num_stripes = 1,
      num_search_stripes = 2,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 5,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 10, .stride = -1, .nstrides = 5 },
           {.start = 10, .stride = -1, .nstrides = 5 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 9, .size = -5}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "simple inverted stripe");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 10, .stride = 1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 9, .size = 2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search inc stripe over inc gap");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 25, .stride = -1, .nstrides = 11},
         {.start = 10, .stride = -1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 10, .stride = 1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 11, .size = -2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search inc stripe over dec gap");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 25, .stride = -1, .nstrides = 11},
         {.start = 10, .stride = -1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 15, .stride = -1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 10, .size = 2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over dec gap");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 15, .stride = -1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 10, .size = -2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over inc gap");
  }

  {
    enum {
      num_stripes = 3,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 8,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10},
         {.start = 29, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 32, .stride = -1, .nstrides = 30 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 23, .size = -22}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over 2 inc gap");
  }

  {
    enum {
      num_stripes = 5,
      num_search_stripes = 1,
      num_ref_pos_ext = 5,
      num_ref_unmatched = 0,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10},
         {.start = 29, .stride = 1, .nstrides = 10},
         {.start = 14, .stride = -1, .nstrides = 4},
         {.start = 28, .stride = -1, .nstrides = 4}},
      search_stripes[num_search_stripes]
        = {{.start = 32, .stride = -1, .nstrides = 30 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
      = {{.start = 23, .size = -4},{.start = 34, .size = 4},
         {.start = 19, .size = -10},{.start = 30, .size = 4},
         {.start = 9, .size = -8}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over jumbled stripes");
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
do_tests(Xt_idxlist idxlist, const Xt_int *ref_indices, int num_indices) {

  check_idxlist(idxlist, ref_indices, num_indices);

  struct Xt_stripe * stripes;
  int num_stripes;
  Xt_idxlist temp_idxlist;

  xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);
  temp_idxlist = xt_idxvec_from_stripes_new(stripes, num_stripes);

  check_idxlist(temp_idxlist, ref_indices, num_indices);

  xt_idxlist_delete(temp_idxlist);

  free(stripes);

  {
    // test packing and unpacking
    Xt_idxlist idxlist_copy
      = idxlist_pack_unpack_copy(idxlist);

    // check copy
    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up
    xt_idxlist_delete(idxlist_copy);
  }

  {
    // test copying
    Xt_idxlist idxlist_copy = xt_idxlist_copy(idxlist);

    // check copy
    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up
    xt_idxlist_delete(idxlist_copy);
  }
}

static void
stripe_test_general(int num_stripes, const struct Xt_stripe *stripes,
                    int num_ref_indices, const Xt_int *ref_indices)
{
  Xt_idxlist idxstripes = xt_idxstripes_new(stripes, num_stripes);
  do_tests(idxstripes, ref_indices, num_ref_indices);
  int num_ext;
  struct Xt_pos_ext *pos_ext;
  xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, num_stripes, stripes,
                                           &num_ext, &pos_ext, 1);
  bool mismatch = false;
  size_t num_pos = 0;
  for (int i = 0; i < num_ext; ++i)
  {
    size_t ext_size = (size_t)pos_ext[i].size;
    mismatch |= (num_pos != (size_t)pos_ext[i].start);
    num_pos += ext_size;
  }
  if (mismatch)
    PUT_ERR("position/start mismatch\n");
  if (num_pos != (size_t)xt_idxlist_get_num_indices(idxstripes))
    PUT_ERR("index list length/positions overlap mismatch\n");

  free(pos_ext);
  xt_idxlist_delete(idxstripes);
  Xt_idxlist idxvec = xt_idxvec_new(ref_indices, num_ref_indices);
  idxstripes = xt_idxstripes_from_idxlist_new(idxvec);
  check_idxlist(idxstripes, ref_indices, num_ref_indices);
  xt_idxlist_delete(idxvec);
  xt_idxlist_delete(idxstripes);
}

static void
stripe_test_general1(void)
{
  enum {
    num_stripes = 3,
    num_ref_indices = 15,
  };
  static const struct Xt_stripe stripes[num_stripes]
    = {{.start =  0, .stride = 1, .nstrides = 5},
       {.start = 10, .stride = 1, .nstrides = 5},
       {.start = 20, .stride = 1, .nstrides = 5}};
  static const Xt_int ref_indices[num_ref_indices]
    = {0,1,2,3,4, 10,11,12,13,14, 20,21,22,23,24};
  stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
}

static void
stripe_test_general2(void)
{
  enum {
    num_stripes = 3,
    num_ref_indices = 15,
  };
  static const struct Xt_stripe stripes[num_stripes]
    = {{.start =  0, .stride = 1, .nstrides = 5},
       {.start = 10, .stride = 2, .nstrides = 5},
       {.start = 20, .stride = 3, .nstrides = 5}};
  static const Xt_int ref_indices[num_ref_indices]
    = {0,1,2,3,4, 10,12,14,16,18, 20,23,26,29,32};

  stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
}

static void
stripe_test_general3(void)
{
  enum {
    num_stripes = 2,
    num_ref_indices = 10,
  };
  static const struct Xt_stripe stripes[num_stripes]
    = {{.start = 0, .stride = 6, .nstrides = 5},
       {.start = 1, .stride = 3, .nstrides = 5}};
  static const Xt_int ref_indices[num_ref_indices]
    = {0,6,12,18,24, 1,4,7,10,13};
  stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
}

static void
stripe_test_general4(void)
{
  enum {
    num_stripes = 2,
    num_ref_indices = 10,
  };
  static const struct Xt_stripe stripes[num_stripes]
    = {{.start = 0, .stride = -1, .nstrides = 5},
       {.start = 1, .stride =  1, .nstrides = 5}};
  static const Xt_int ref_indices[num_ref_indices]
    = {0,-1,-2,-3,-4, 1,2,3,4,5};
  stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
}

static void
stripe_test_general5(void)
{
  enum {
    num_stripes = 2,
    num_ref_indices = 10,
  };
  static const struct Xt_stripe stripes[num_stripes]
    = {{.start = 9, .stride = -2, .nstrides = 5},
       {.start = 0, .stride =  2, .nstrides = 5}};
  static const Xt_int ref_indices[num_ref_indices]
    = {9,7,5,3,1, 0,2,4,6,8};
  stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
}

static void
stripe_test_general6(void)
{
  enum {
    num_stripes = 1,
    num_ref_indices = 0,
  };
  static const struct Xt_stripe stripes[num_stripes]
    = {{.start = 179, .stride = -2, .nstrides = 0}};
  stripe_test_general(num_stripes, stripes, num_ref_indices, NULL);
}

static void
stripe_test_intersection(int num_stripes_a, const struct Xt_stripe *stripes_a,
                         int num_stripes_b, const struct Xt_stripe *stripes_b,
                         int num_ref_indices, const Xt_int *ref_indices)
{
  Xt_idxlist idxstripes_a = xt_idxstripes_new(stripes_a, num_stripes_a),
    idxstripes_b = xt_idxstripes_new(stripes_b, num_stripes_b);

  // compute intersections
  Xt_idxlist intersection[2]
    = { xt_idxlist_get_intersection(idxstripes_a, idxstripes_b),
        xt_idxlist_get_intersection(idxstripes_b, idxstripes_a) };

  // check intersection
  do_tests(intersection[0], ref_indices, num_ref_indices);
  do_tests(intersection[1], ref_indices, num_ref_indices);

  // clean up
  xt_idxlist_delete(idxstripes_a);
  xt_idxlist_delete(idxstripes_b);
  xt_idxlist_delete(intersection[0]);
  xt_idxlist_delete(intersection[1]);
}

static void
stripe_test_asymmetric_intersection(
  int num_stripes_a, const struct Xt_stripe *stripes_a,
  int num_stripes_b, const struct Xt_stripe *stripes_b,
  int num_ref_indices_a, const Xt_int *ref_indices_a,
  int num_ref_indices_b, const Xt_int *ref_indices_b)
{
  Xt_idxlist idxstripes_a = xt_idxstripes_new(stripes_a, num_stripes_a),
    idxstripes_b = xt_idxstripes_new(stripes_b, num_stripes_b);

  // compute intersection
  Xt_idxlist intersection[2]
    = { xt_idxlist_get_intersection(idxstripes_a, idxstripes_b),
        xt_idxlist_get_intersection(idxstripes_b, idxstripes_a) };

  // check intersection
  do_tests(intersection[0], ref_indices_a, num_ref_indices_a);
  do_tests(intersection[1], ref_indices_b, num_ref_indices_b);

  // clean up
  xt_idxlist_delete(idxstripes_a);
  xt_idxlist_delete(idxstripes_b);
  xt_idxlist_delete(intersection[0]);
  xt_idxlist_delete(intersection[1]);
}

/* test whether
 *  xt_idxlist_get_index_at_position and
 *  xt_idxlist_get_indices_at_positions
 * give the same results
 */
static void
check_idxvec_get_indices_at_positions(int num_stripes,
                                      const struct Xt_stripe *stripes,
                                      int num_pos, const int *pos)
{
  static const Xt_int undef_idx = XT_INT_MIN;
  Xt_idxlist idxlist = xt_idxstripes_new(stripes, num_stripes);
  Xt_int ref_sel_idx[num_pos];
  int ref_undef_count = 0;
  for (int i=0; i<num_pos; i++) {
    int p = pos[i];
    if (xt_idxlist_get_index_at_position(idxlist, p, &ref_sel_idx[i]) != 0) {
      ref_sel_idx[i] = undef_idx;
      ref_undef_count++;
    }
  }
  Xt_int sel_idx[num_pos];
  int undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos, num_pos,
                                                        sel_idx, undef_idx);
  if (undef_count != ref_undef_count)
    PUT_ERR("test_idxstripes.c: (undef_count != ref_undef_count)\n");
  bool mismatch = false;
  for (int i=0; i<num_pos; i++)
    mismatch |= (sel_idx[i] != ref_sel_idx[i]);
  if (mismatch)
    PUT_ERR("test_idxstripes.c: sel_idx mismatch\n");

  xt_idxlist_delete(idxlist);
}

static struct Xt_pos_ext *
get_idxlist_pos_exts_of_index_stripes(
  Xt_idxlist idxlist,
  size_t num_stripes, const struct Xt_stripe stripes[num_stripes],
  size_t *num_ext, int single_match_only)
{
  int num_ext_;
  struct Xt_pos_ext *pos_ext_;
  int retval = xt_idxlist_get_pos_exts_of_index_stripes(
    idxlist, (int)num_stripes, stripes, &num_ext_, &pos_ext_,
    single_match_only);
  if (retval != 0 || num_ext_ < 0)
    PUT_ERR("error in xt_idxlist_get_pos_exts_of_index_stripes\n");
  *num_ext = (size_t)num_ext_;
  return pos_ext_;
}


static void
check_idxlist_stripes_pos_ext(Xt_idxlist idxlist,
                              size_t num_stripes,
                              const struct Xt_stripe stripes[num_stripes])
{
  size_t num_ext;
  int single_match_only = 1;
  struct Xt_pos_ext *restrict pos_ext
    = get_idxlist_pos_exts_of_index_stripes(idxlist,
                                            num_stripes, stripes,
                                            &num_ext, single_match_only);

  { // testing of results
    Xt_idxlist intersection
      = xt_idxvec_from_stripes_new(stripes, (int)num_stripes);
    for (size_t i = 0, k = 0; i < num_ext; ++i) {
      int abs_pos_ext_size = abs(pos_ext[i].size),
        jsign = pos_ext[i].size >= 0 ? 1 : -1;
      for (int j = 0; j < abs_pos_ext_size; ++j, ++k) {
        Xt_int intersection_index;
        xt_idxlist_get_index_at_position(intersection, (int)k,
                                         &intersection_index);
        int send_pos = pos_ext[i].start + jsign * j;
        Xt_int orig_index;
        xt_idxlist_get_index_at_position(idxlist, send_pos, &orig_index);
        if (intersection_index != orig_index) {
          PUT_ERR("error in xt_idxlist_get_pos_exts_of_index_stripes\n");
          fprintf(stderr, "intersection pos %zu index %"XT_INT_FMT
                  " orig pos %d index %"XT_INT_FMT"\n",
                  k, intersection_index, send_pos, orig_index);
        }
      }
    }
    xt_idxlist_delete(intersection);
  }
  free(pos_ext);
}


static void
check_bb(int num_stripes, const struct Xt_stripe *stripes,
         unsigned ndim, const Xt_int global_size[ndim],
         const struct Xt_bounds ref_bounds[ndim],
         Xt_int global_start_index)
{
  Xt_idxlist idxstripes = xt_idxstripes_new(stripes, num_stripes);
  struct Xt_bounds bounds[ndim];

  xt_idxlist_get_bounding_box(idxstripes, ndim, global_size,
                              global_start_index, bounds);

  bool mismatch = false;
  for (size_t i = 0; i < ndim; ++i)
    mismatch |= ((bounds[i].size != ref_bounds[i].size) |
                 ((ref_bounds[i].size != 0)
                  & (bounds[i].start != ref_bounds[i].start)));
  if (mismatch)
    PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

  xt_idxlist_delete(idxstripes);
}

static void
check_pos_ext(size_t num_stripes, const struct Xt_stripe stripes[num_stripes],
              size_t num_search_stripes,
              const struct Xt_stripe search_stripes[num_search_stripes],
              size_t num_ref_pos_ext,
              const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext],
              int single_match_only, int ref_unmatched,
              const char *test_desc)
{
  Xt_idxlist idxstripes = xt_idxstripes_new(stripes, (int)num_stripes);
  int num_ext;
  struct Xt_pos_ext *pos_ext;
  int unmatched
    = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes,
                                               (int)num_search_stripes,
                                               search_stripes,
                                               &num_ext, &pos_ext,
                                               single_match_only);
  (void)test_desc;
  if (unmatched != ref_unmatched)
    PUT_ERR("error in number of unmatched indices for %s", test_desc);
  if (num_ext < 0 || (size_t)num_ext != num_ref_pos_ext)
    PUT_ERR("error finding %s\n", test_desc);
  bool mismatched_start = false,
    mismatched_size = false;
  for (size_t i = 0; i < num_ref_pos_ext; ++i) {
    mismatched_start |= (pos_ext[i].start != ref_pos_ext[i].start);
    mismatched_size |= (pos_ext[i].size != ref_pos_ext[i].size);
  }
  if (mismatched_start)
    PUT_ERR("incorrect starting position found in %s\n", test_desc);
  if (mismatched_size)
    PUT_ERR("incorrect position extent length found in %s\n", test_desc);
  free(pos_ext);
  xt_idxlist_delete(idxstripes);
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
