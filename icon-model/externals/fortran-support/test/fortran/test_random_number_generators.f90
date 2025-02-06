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

MODULE TEST_mo_random_number_generators
  USE FORTUTF
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  USE helpers, ONLY: assert_statistics, calculate_mean_wp, calculate_variance_wp

CONTAINS
  SUBROUTINE TEST_CALCULATE_STATISTICS_WP
    REAL(wp) :: array(8) = (/2.0_wp, 4.0_wp, 4.0_wp, 4.0_wp, 5.0_wp, 5.0_wp, 7.0_wp, 9.0_wp/)
    CALL TAG_TEST("TEST_calculate_mean_wp")
    CALL ASSERT_ALMOST_EQUAL(calculate_mean_wp(array), 5.0_wp)
    CALL TAG_TEST("TEST_calculate_variance_wp")
    CALL ASSERT_ALMOST_EQUAL(calculate_variance_wp(array), 4.0_wp)
  END SUBROUTINE

  SUBROUTINE TEST_CLIP_WP
    USE mo_random_number_generators
    REAL(wp) :: array(6) = (/4.1_wp, -4.1_wp, 2.0_wp, -2.0_wp, 1.0_wp, -1.0_wp/)
    CALL TAG_TEST("TEST_clip_wp")
    CALL clip(2.0_wp, array)
    CALL ASSERT_EQUAL(array, (/2.0_wp, -2.0_wp, 2.0_wp, -2.0_wp, 1.0_wp, -1.0_wp/))
  END SUBROUTINE

  SUBROUTINE TEST_UNIFORM_RANDOM_NUMBER_WP_BASIC
    USE mo_random_number_generators
    INTEGER, PARAMETER :: n = 1e5
    INTEGER :: i
    INTEGER(int32) :: rng_state
    REAL(wp) :: array(n)

    rng_state = initialize_random_number_generator(11, 17)
    DO i = 1, n
      array(i) = generate_uniform_random_number(rng_state)
    END DO

    CALL assert_statistics("TEST_UNIFORM_RANDOM_NUMBER_WP_BASIC", &
                           array, mean=0.5_wp, variance=(1.0_wp/12.0_wp), tol=1.0e-2_wp, max=1.0_wp, min=0.0_wp)

  END SUBROUTINE

  SUBROUTINE TEST_UNIFORM_RANDOM_NUMBER_WP_ISEED
    USE mo_random_number_generators
    INTEGER, PARAMETER :: n = 1e5
    INTEGER :: i
    INTEGER(int32) :: rng_state
    REAL(wp) :: array(n)

    DO i = 1, n
      rng_state = initialize_random_number_generator(i, 11)
      array(i) = generate_uniform_random_number(rng_state)
    END DO

    CALL assert_statistics("TEST_UNIFORM_RANDOM_NUMBER_WP_ISEED", &
                           array, mean=0.5_wp, variance=(1.0_wp/12.0_wp), tol=1.0e-2_wp, max=1.0_wp, min=0.0_wp)

  END SUBROUTINE

  SUBROUTINE TEST_UNIFORM_RANDOM_NUMBER_WP_JSEED
    USE mo_random_number_generators
    INTEGER, PARAMETER :: n = 1e5
    INTEGER :: i
    INTEGER(int32) :: rng_state
    REAL(wp) :: array(n)

    DO i = 1, n
      rng_state = initialize_random_number_generator(11, i)
      array(i) = generate_uniform_random_number(rng_state)
    END DO

    CALL assert_statistics("TEST_UNIFORM_RANDOM_NUMBER_WP_JSEED", &
                           array, mean=0.5_wp, variance=(1.0_wp/12.0_wp), tol=1.0e-2_wp, max=1.0_wp, min=0.0_wp)

  END SUBROUTINE

  SUBROUTINE TEST_RANDOM_NORMAL_VALUES_WP
    USE mo_random_number_generators
    INTEGER, PARAMETER :: n = 1e5 + 1
    REAL(wp) :: array(n)
    REAL(wp) :: mean = -5.0_wp
    REAL(wp) :: variance = 7.1_wp
    REAL(wp) :: value_range = 1.0e32_wp

    CALL random_normal_values(seed=11, values_range=value_range, values=array, mean=mean, stdv=SQRT(variance))

    CALL assert_statistics("TEST_RANDOM_NORMAL_VALUES_WP", &
                           array, mean=mean, variance=variance, tol=1.0e-2_wp, max=value_range, min=-value_range)

  END SUBROUTINE

END MODULE TEST_mo_random_number_generators
