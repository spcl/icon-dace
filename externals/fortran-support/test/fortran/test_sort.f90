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
    USE mo_util_sort
    REAL(wp) :: to_sort(6) = (/144.4, 58.6, 4.3, 7.8, 10.0, 11.0/)
    CALL TAG_TEST("TEST_quicksort_real_1")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .FALSE.)
    CALL quicksort(to_sort)
    CALL TAG_TEST("TEST_quicksort_real_2")
    CALL ASSERT_EQUAL(is_sorted_real(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_int
    USE mo_util_sort
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    CALL TAG_TEST("TEST_quicksort_int_1")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)
    CALL quicksort(to_sort)
    CALL TAG_TEST("TEST_quicksort_int_2")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_permutation_int
    USE mo_util_sort
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    INTEGER :: idx_permutation(6) = (/1, 2, 3, 4, 5, 6/)
    CALL TAG_TEST("TEST_quicksort_permutation_int_1")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)
    CALL quicksort(to_sort, idx_permutation)
    CALL TAG_TEST("TEST_quicksort_permutation_int_2")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .TRUE.)
    CALL TAG_TEST("TEST_quicksort_permutation_int_3")
    CALL ASSERT_EQUAL(has_same_values_int(idx_permutation, (/3, 4, 5, 6, 2, 1/)), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_quicksort_string
    USE mo_util_sort
    CHARACTER :: to_sort(6) = (/'A', 'C', 'Y', 'E', 'S', 'H'/)
    CALL TAG_TEST("TEST_quicksort_string_1")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .FALSE.)
    CALL quicksort(to_sort)
    CALL TAG_TEST("TEST_quicksort_string_2")
    CALL ASSERT_EQUAL(is_sorted_string(to_sort), .TRUE.)
  END SUBROUTINE

  SUBROUTINE TEST_insertion_sort_int
    USE mo_util_sort
    INTEGER :: to_sort(6) = (/144, 58, 4, 7, 10, 11/)
    CALL TAG_TEST("TEST_insertion_sort_int_1")
    CALL ASSERT_EQUAL(is_sorted_int(to_sort), .FALSE.)
    CALL quicksort(to_sort)
    CALL TAG_TEST("TEST_insertion_sort_int_2")
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
