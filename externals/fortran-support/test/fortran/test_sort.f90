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

MODULE TEST_mo_util_sort
  USE FORTUTF
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64

CONTAINS
  SUBROUTINE TEST_quicksort_real
    USE mo_util_sort, ONLY: quicksort
    REAL(wp) :: to_sort(6) = (/144.4, 58.6, 4.3, 7.8, 10.0, 11.0/)
    CALL TAG_TEST("TEST_quicksort_real_before")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_after")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_real2
    USE mo_util_sort, ONLY: quicksort
    REAL(wp) :: to_sort(6) = (/144.4, 11.0, 4.3, 58.6, 10.0, 7.8/)
    CALL TAG_TEST("TEST_quicksort_real2_before")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real2_after")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_real_random
    USE mo_util_sort, ONLY: quicksort
    REAL(wp) :: to_sort(6)
    ! Generate random numbers geq 0.0 and < 256.0
    CALL RANDOM_NUMBER(to_sort)
    to_sort = to_sort*256

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_real_random")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_permutation_real
    USE mo_util_sort, ONLY: quicksort
    REAL(wp) :: to_sort(6) = (/144.4, 58.6, 4.3, 7.8, 10.0, 11.0/)
    INTEGER :: idx_permutation(6) = (/1, 2, 3, 4, 5, 6/)
    CALL TAG_TEST("TEST_quicksort_permutation_real_before")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .FALSE.)

    CALL quicksort(to_sort, idx_permutation)

    CALL TAG_TEST("TEST_quicksort_permutation_real_after")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .TRUE.)
    CALL TAG_TEST("TEST_quicksort_permutation_real_permutation")
    CALL ASSERT_EQUAL(has_same_values_int(idx_permutation, (/3, 4, 5, 6, 2, 1/)), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_int
    USE mo_util_sort, ONLY: quicksort
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    CALL TAG_TEST("TEST_quicksort_int_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_int_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_int2
    USE mo_util_sort, ONLY: quicksort
    INTEGER :: to_sort(6) = (/58, 4, 144, 10, 7, 11/)
    CALL TAG_TEST("TEST_quicksort_int2_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_int2_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_int_random
    USE mo_util_sort, ONLY: quicksort
    INTEGER :: to_sort(6)
    REAL(wp) :: random_wp(6)
    ! Generate random numbers between 0 and 255
    CALL RANDOM_NUMBER(random_wp)
    to_sort = FLOOR(random_wp*256)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_int_random")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_permutation_int
    USE mo_util_sort, ONLY: quicksort
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    INTEGER :: idx_permutation(6) = (/1, 2, 3, 4, 5, 6/)
    CALL TAG_TEST("TEST_quicksort_permutation_int_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL quicksort(to_sort, idx_permutation)

    CALL TAG_TEST("TEST_quicksort_permutation_int_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
    CALL TAG_TEST("TEST_quicksort_permutation_int_permutation")
    CALL ASSERT_EQUAL(has_same_values_int(idx_permutation, (/3, 4, 5, 6, 2, 1/)), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string
    USE mo_util_sort, ONLY: quicksort
    CHARACTER :: to_sort(6) = (/'A', 'C', 'Y', 'E', 'S', 'H'/)
    CALL TAG_TEST("TEST_quicksort_string_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string2
    USE mo_util_sort, ONLY: quicksort
    CHARACTER :: to_sort(6) = (/'Y', 'H', 'A', 'S', 'E', 'C'/)
    CALL TAG_TEST("TEST_quicksort_string2_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string2_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string3
    USE mo_util_sort, ONLY: quicksort
    CHARACTER :: to_sort(6) = (/'P', 'M', 'W', 'G', 'K', 'D'/)
    CALL TAG_TEST("TEST_quicksort_string3_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string3_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string4
    USE mo_util_sort, ONLY: quicksort
    CHARACTER :: to_sort(6) = (/'B', 'L', 'Q', 'S', 'Z', 'T'/)
    CALL TAG_TEST("TEST_quicksort_string4_before")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)

    CALL quicksort(to_sort)

    CALL TAG_TEST("TEST_quicksort_string4_after")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_insertion_sort_int
    USE mo_util_sort, ONLY: insertion_sort
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    CALL TAG_TEST("TEST_insertion_sort_int_before")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)

    CALL insertion_sort(to_sort)

    CALL TAG_TEST("TEST_insertion_sort_int_after")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_insertion_sort_int_random
    USE mo_util_sort, ONLY: insertion_sort
    INTEGER :: to_sort(6)
    REAL(wp) :: random_wp(6)
    ! Generate random numbers between 0 and 255
    CALL RANDOM_NUMBER(random_wp)
    to_sort = FLOOR(random_wp*256)

    CALL insertion_sort(to_sort)

    CALL TAG_TEST("TEST_insertion_sort_int_random")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  LOGICAL FUNCTION has_same_values_int(array, ref)
    INTEGER, INTENT(IN) :: array(:), ref(:)
    INTEGER :: i

    has_same_values_int = .TRUE.
    DO i = 1, SIZE(array)
      IF (array(i) /= ref(i)) THEN
        has_same_values_int = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION has_same_values_int

  LOGICAL FUNCTION is_sorted_real(array)
    REAL(wp), INTENT(IN) :: array(:)
    INTEGER :: i

    is_sorted_real = .TRUE.
    DO i = 1, SIZE(array) - 1
      IF (array(i) > array(i + 1)) THEN
        is_sorted_real = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION is_sorted_real

  LOGICAL FUNCTION is_sorted_int(array)
    INTEGER, INTENT(IN) :: array(:)
    INTEGER :: i

    is_sorted_int = .TRUE.
    DO i = 1, SIZE(array) - 1
      IF (array(i) > array(i + 1)) THEN
        is_sorted_int = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION is_sorted_int

  LOGICAL FUNCTION is_sorted_string(array)
    CHARACTER, INTENT(IN) :: array(:)
    INTEGER :: i

    is_sorted_string = .TRUE.
    DO i = 1, SIZE(array) - 1
      IF (array(i) > array(i + 1)) THEN
        is_sorted_string = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION is_sorted_string

END MODULE TEST_mo_util_sort
