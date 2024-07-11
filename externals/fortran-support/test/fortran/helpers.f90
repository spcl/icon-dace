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

MODULE helpers
  USE ISO_C_BINDING, ONLY: c_double
  USE mo_io_units, ONLY: find_next_free_unit
  USE mo_exception, ONLY: message
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64

CONTAINS

  ! deterministic sequence of "random" numbers.
  ! (rval in [-1,1], linear congruential generator)
  SUBROUTINE rrand(seed, rval)
    INTEGER, INTENT(INOUT) :: seed
    REAL(c_double), INTENT(OUT) :: rval
    !
    INTEGER, PARAMETER :: rand_m = ISHFT(2, 15) + 1 ! 2**16+1
    INTEGER, PARAMETER :: rand_a = 75
    INTEGER, PARAMETER :: rand_c = 74
    seed = MOD((rand_a*seed + rand_c), rand_m)
    rval = 2.0_c_double*REAL(seed, c_double)/REAL(rand_m - 1, c_double) - 1.0_c_double
  END SUBROUTINE rrand

  SUBROUTINE custom_exit()
    STOP
  END SUBROUTINE

  SUBROUTINE custom_exit_dummy()
  END SUBROUTINE

  SUBROUTINE open_new_logfile(nerr, file)
    CHARACTER(len=*), INTENT(IN) :: file
    INTEGER, INTENT(OUT) :: nerr

    nerr = find_next_free_unit(10, 20)
    OPEN (unit=nerr, file=file, status='replace', action='write')

  END SUBROUTINE open_new_logfile

  SUBROUTINE open_logfile(nerr, file)
    CHARACTER(len=*), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: nerr

    OPEN (unit=nerr, file=file, status='old', action='read')

  END SUBROUTINE open_logfile

  FUNCTION calculate_mean_wp(array) RESULT(mean)
    REAL(wp), INTENT(IN) :: array(:)
    REAL(wp)             :: mean
    INTEGER :: i, n

    mean = 0.0_wp
    variance = 0.0_wp

    n = SIZE(array)

    DO i = 1, SIZE(array)
      mean = mean + array(i)
    END DO
    mean = mean/n

  END FUNCTION

  FUNCTION calculate_variance_wp(array) RESULT(variance)
    REAL(wp), INTENT(IN) :: array(:)
    REAL(wp)             :: mean, variance
    INTEGER :: i, n

    mean = 0.0_wp
    variance = 0.0_wp

    n = SIZE(array)
    mean = calculate_mean_wp(array)

    DO i = 1, SIZE(array)
      variance = variance + (array(i) - mean)**2
    END DO
    variance = variance/n

  END FUNCTION

  SUBROUTINE assert_statistics(test_name, array, mean, variance, tol, max, min)
    USE FORTUTF
    CHARACTER(LEN=*), INTENT(IN) :: test_name
    REAL(wp), INTENT(IN) :: array(:)
    REAL(wp), INTENT(IN) :: mean, variance, tol, max, min

    CALL TAG_TEST(test_name//"__mean")
    CALL ASSERT_ALMOST_EQUAL(calculate_mean_wp(array), mean, tol)
    CALL TAG_TEST(test_name//"__variance")
    CALL ASSERT_ALMOST_EQUAL(calculate_variance_wp(array), variance, tol)
    CALL TAG_TEST(test_name//"__max")
    CALL ASSERT_LESS_THAN_EQUAL(MAXVAL(array), max)
    CALL TAG_TEST(test_name//"__min")
    CALL ASSERT_GREATER_THAN_EQUAL(MINVAL(array), min)

  END SUBROUTINE

END MODULE helpers
