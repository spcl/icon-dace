! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE test_mo_octree
  USE FORTUTF
  USE ISO_C_BINDING, ONLY: c_double

CONTAINS

  SUBROUTINE TEST_octree_1

    USE mo_octree, ONLY: t_range_octree, octree_init, octree_finalize, &
      &                       octree_count_point, octree_query_point
    USE helpers, ONLY: rrand
    IMPLICIT NONE

    INTEGER, PARAMETER :: SUCCESS = 0
    INTEGER, PARAMETER :: wp = c_double
    INTEGER, PARAMETER :: NOBJECTS = 10000
    INTEGER, PARAMETER :: NQUERIES = 50

    INTEGER, PARAMETER :: TEST_RESULT1(422) = &
                          [6454, 6450, 6452, 6310, 6324, 6325, 6388, 6351, 6579, 6625, 6493, &
                           6547, 6097, 6092, 6125, 6135, 6109, 6115, 6009, 6014, 6042, 6029, &
                           6293, 6269, 6277, 6141, 6153, 6220, 7097, 7134, 7013, 7023, 7000, &
                           7073, 7067, 7050, 7269, 7251, 7246, 7319, 7330, 7286, 7155, 7222, &
                           6757, 6769, 6732, 6649, 6710, 6917, 6966, 6944, 6961, 6829, 6877, &
                           6850, 6870, 5113, 5141, 5147, 5012, 5055, 5066, 5033, 5046, 5281, &
                           5186, 5167, 5171, 5207, 5214, 4774, 4746, 4753, 4792, 4666, 4716, &
                           4732, 4699, 4712, 4988, 4922, 4851, 4848, 4861, 5745, 5798, 5815, &
                           5788, 5791, 5680, 5726, 5700, 5952, 5847, 5857, 5819, 5436, 5423, &
                           5422, 5447, 5340, 5596, 5614, 5505, 5528, 5539, 9098, 9095, 9132, &
                           9032, 9039, 9008, 9080, 9043, 9049, 9336, 9171, 9189, 9250, 9211, &
                           9223, 8770, 8756, 8829, 8798, 8691, 8708, 8741, 8711, 8954, 8859, &
                           8872, 8907, 8908, 8667, 9781, 9784, 9800, 9811, 9810, 9676, 9731, &
                           9922, 9978, 9871, 9867, 9888, 9436, 9490, 9496, 9493, 9605, 9595, &
                           9626, 9525, 9512, 7764, 7759, 7809, 7811, 7688, 7675, 7674, 7725, &
                           7662, 7937, 7946, 7952, 7931, 7978, 7960, 7850, 7860, 7829, 7881, &
                           7435, 7432, 7441, 7424, 7337, 7350, 7400, 7380, 7376, 7614, 7612, &
                           7594, 7521, 7501, 7564, 7546, 8421, 8416, 8413, 8495, 8463, 8374, &
                           8385, 8585, 8569, 8579, 8547, 8145, 8019, 8016, 8007, 8064, 8063, &
                           8041, 8277, 8199, 8175, 8227, 8241, 406, 2008, 1985, 2018, 2203, 2177, &
                           2179, 2228, 2214, 2133, 2113, 2151, 1794, 1785, 1780, 1817, 1806, &
                           1808, 1739, 1933, 1915, 1950, 1858, 1869, 435, 1893, 2581, 2586, 2631, &
                           2512, 2544, 2738, 293, 2717, 2711, 2713, 2720, 2772, 2756, 2666, 2662, &
                           2646, 363, 2317, 2362, 2340, 2348, 2272, 2298, 2283, 335, 337, 340, &
                           2499, 2385, 2383, 591, 570, 560, 1174, 1150, 562, 1067, 1079, 1063, &
                           1116, 735, 758, 766, 674, 654, 680, 636, 684, 854, 907, 891, 882, 615, &
                           816, 1530, 1573, 1453, 1668, 489, 1607, 1612, 1276, 1271, 1302, 1296, &
                           1195, 1232, 1233, 1227, 1335, 3931, 3972, 3965, 3875, 3883, 3858, &
                           3890, 4061, 84, 3995, 4032, 3673, 3695, 3625, 109, 111, 3789, 3792, &
                           3741, 3742, 3719, 3725, 3753, 3760, 4479, 4481, 4461, 4499, 4401, &
                           4451, 1, 4632, 4538, 4571, 4569, 4581, 4576, 4577, 4578, 4192, 4203, &
                           4232, 64, 4157, 4133, 4320, 4382, 4358, 41, 4271, 4270, 4296, 3119, &
                           3115, 3121, 3142, 3126, 224, 3240, 3236, 3254, 3250, 3257, 3212, 210, &
                           2855, 2880, 2795, 2793, 2807, 2991, 3004, 258, 2913, 3585, 3486, 3472, &
                           3344, 3493, 3348, 3338, 3412, 165, 3294, 3543, 3573, 3287, 3331, 3532, &
                           3314, 157, 3908, 9379, 6021, 229, 6941, 4613, 6041, 2455, 6512, 4742, &
                           281, 9191, 4586, 1684, 8617, 2985, 1722, 7057, 8548, 4991]

    TYPE(t_range_octree) :: octree ! octree data structure
    REAL(wp) :: brange(2, 3) ! box range (min/max, dim=1,2,3)
    REAL(wp), DIMENSION(NOBJECTS, 3) :: pmin, pmax ! dim=(1,...,nobjects, x/y/z)  : range (corners)
    INTEGER :: iindex(NOBJECTS) ! object index (optional)
    REAL(wp) :: p(3) ! point to insert
    INTEGER :: num_obj ! number of found objects
    INTEGER :: obj_list(NOBJECTS) ! result: list of objects in traversed boxes.
    LOGICAL :: lassert, tmp(NOBJECTS)
    INTEGER :: i, j, seed
    REAL(c_double) :: r1, r2
    CHARACTER(len=100) :: label

    ! generate some test data
    seed = 5
    DO i = 1, nobjects
      DO j = 1, 3
        CALL rrand(seed, r1)
        CALL rrand(seed, r2)
        pmin(i, j) = MIN(r1, r2)
        pmax(i, j) = MAX(r1, r2)
      END DO
    END DO
    brange(1, 1:3) = [-1.0_wp, -1.0_wp, -1.0_wp]
    brange(2, 1:3) = [1.0_wp, 1.0_wp, 1.0_wp]
    ! constructor
    iindex = [(i, i=1, nobjects)]
    CALL octree_init(octree, brange, pmin, pmax, iindex)

    DO i = 1, NQUERIES
      ! test point
      CALL rrand(seed, p(1))
      CALL rrand(seed, p(2))
      CALL rrand(seed, p(3))

      ! query number of box-shaped objects that are hit by `p`
      num_obj = octree_count_point(octree, p)
      ! query list of objects that are hit by `p`
      obj_list(:) = -1
      CALL octree_query_point(octree, p, obj_list)

      ! check: length of result array matches `num_obj`
      WRITE (label, '(A,I3,A)') 'TEST_octree', i, '_1'
      CALL TAG_TEST(label)
      CALL ASSERT_EQUAL(num_obj == COUNT(obj_list(:) /= -1), .TRUE.)

      ! check: result array does not contain duplicates
      tmp(:) = .FALSE.; tmp(obj_list(1:num_obj)) = .TRUE.
      WRITE (label, '(A,I3,A)') 'TEST_octree', i, '_2'
      CALL TAG_TEST(label)
      CALL ASSERT_EQUAL(COUNT(tmp) == num_obj, .TRUE.)
      ! check: result list for 1st query point
      WRITE (label, '(A,I3,A)') 'TEST_octree', i, '_3'
      CALL TAG_TEST(label)
      CALL ASSERT_EQUAL((i /= 1) .OR. (ALL(TEST_RESULT1 == obj_list(1:num_obj))), .TRUE.)
    END DO

    ! destructor
    CALL octree_finalize(octree)
  END SUBROUTINE

END MODULE test_mo_octree
