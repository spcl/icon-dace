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

!! Module containing basic sorting algorithms (required until
!! a Fortran STL has been invented).

MODULE mo_util_sort

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: quicksort



  PUBLIC :: insertion_sort

  ! generic interface for in-situ QuickSort sorting routine
  INTERFACE quicksort
    MODULE PROCEDURE quicksort_real
    MODULE PROCEDURE quicksort_int
    MODULE PROCEDURE quicksort_permutation_int
    MODULE PROCEDURE quicksort_string
  END INTERFACE quicksort






  INTERFACE swap
    MODULE PROCEDURE swap_permutation_int
    MODULE PROCEDURE swap_int
  END INTERFACE swap

  INTERFACE insertion_sort
    MODULE PROCEDURE insertion_sort_int
  END INTERFACE insertion_sort

  ! A simple tool to reorder DATA.
  CHARACTER(*), PARAMETER :: modname = "mo_util_sort"

CONTAINS

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of REAL values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_real(a, permutation, l_in, r_in)
    REAL(wp), INTENT(INOUT)           :: a(:) !< array for in-situ sorting
    INTEGER, INTENT(INOUT), OPTIONAL :: permutation(:) !< (optional) permutation of indices
    INTEGER, INTENT(IN), OPTIONAL :: l_in, r_in !< left, right partition indices
    ! local variables
    INTEGER  :: i, j, l, r, t_p
    REAL(wp) :: t, v

    IF (PRESENT(l_in)) THEN
      l = l_in
    ELSE
      l = 1
    END IF
    IF (PRESENT(r_in)) THEN
      r = r_in
    ELSE
      r = SIZE(a, 1)
    END IF
    IF (r > l) THEN
      v = a(r)
      i = l - 1
      j = r
      LOOP: DO
        CNTLOOP1: DO
          i = i + 1
          IF (a(i) >= v) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2: DO
          j = j - 1
          IF ((a(j) <= v) .OR. (j == 1)) EXIT CNTLOOP2
        END DO CNTLOOP2
        t = a(i)
        a(i) = a(j)
        a(j) = t
        IF (PRESENT(permutation)) THEN
          t_p = permutation(i)
          permutation(i) = permutation(j)
          permutation(j) = t_p
        END IF
        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t
      IF (PRESENT(permutation)) THEN
        permutation(j) = permutation(i)
        permutation(i) = permutation(r)
        permutation(r) = t_p
      END IF
      CALL quicksort(a, permutation, l, i - 1)
      CALL quicksort(a, permutation, i + 1, r)
    END IF
  END SUBROUTINE quicksort_real

  SUBROUTINE swap_int(a, i, j)
    !> array for in-situ sorting
    INTEGER, INTENT(INOUT)           :: a(:)
    !> indices to be exchanged
    INTEGER, INTENT(IN)              :: i, j
    ! local variables
    INTEGER :: t

    t = a(i)
    a(i) = a(j)
    a(j) = t
  END SUBROUTINE swap_int

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_int(a, l_in, r_in)
    INTEGER, INTENT(INOUT)           :: a(:) !< array for in-situ sorting
    INTEGER, INTENT(IN), OPTIONAL :: l_in, r_in !< left, right partition indices
    ! local variables
    INTEGER :: i, j, l, r, t, v, m

    IF (PRESENT(l_in)) THEN
      l = l_in
    ELSE
      l = 1
    END IF
    IF (PRESENT(r_in)) THEN
      r = r_in
    ELSE
      r = SIZE(a, 1)
    END IF
    IF (r > l) THEN
      i = l - 1
      j = r

      ! median-of-three selection of partitioning element
      IF ((r - l) > 3) THEN
        m = (l + r)/2
        IF (a(l) > a(m)) CALL swap(a, l, m)
        IF (a(l) > a(r)) THEN
          CALL swap(a, l, r)
        ELSE IF (a(r) > a(m)) THEN
          CALL swap(a, r, m)
        END IF
      END IF

      v = a(r)
      LOOP: DO
        CNTLOOP1: DO
          i = i + 1
          IF (a(i) >= v) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2: DO
          j = j - 1
          IF ((a(j) <= v) .OR. (j == 1)) EXIT CNTLOOP2
        END DO CNTLOOP2
        t = a(i)
        a(i) = a(j)
        a(j) = t
        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t
      CALL quicksort(a, l, i - 1)
      CALL quicksort(a, i + 1, r)
    END IF
  END SUBROUTINE quicksort_int
  SUBROUTINE swap_permutation_int(a, i, j, permutation)
    !> array for in-situ sorting
    INTEGER, INTENT(INOUT)           :: a(:)
    !> indices to be exchanged
    INTEGER, INTENT(IN)              :: i, j
    !> (optional) permutation of indices
    INTEGER, INTENT(INOUT)           :: permutation(:)
    ! local variables
    INTEGER :: t, t_p

    t = a(i)
    a(i) = a(j)
    a(j) = t
    t_p = permutation(i)
    permutation(i) = permutation(j)
    permutation(j) = t_p
  END SUBROUTINE swap_permutation_int

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_permutation_int(a, permutation, l_in, r_in)
    INTEGER, INTENT(INOUT)           :: a(:) !< array for in-situ sorting
    INTEGER, INTENT(INOUT)           :: permutation(:) !< (optional) permutation of indices
    INTEGER, INTENT(IN), OPTIONAL :: l_in, r_in !< left, right partition indices
    ! local variables
    INTEGER :: i, j, l, r, t_p, t, v, m

    IF (PRESENT(l_in)) THEN
      l = l_in
    ELSE
      l = 1
    END IF
    IF (PRESENT(r_in)) THEN
      r = r_in
    ELSE
      r = SIZE(a, 1)
    END IF
    IF (r > l) THEN
      i = l - 1
      j = r

      ! median-of-three selection of partitioning element
      IF ((r - l) > 3) THEN
        m = (l + r)/2
        IF (a(l) > a(m)) CALL swap(a, l, m, permutation)
        IF (a(l) > a(r)) THEN
          CALL swap(a, l, r, permutation)
        ELSE IF (a(r) > a(m)) THEN
          CALL swap(a, r, m, permutation)
        END IF
      END IF

      v = a(r)
      LOOP: DO
        CNTLOOP1: DO
          i = i + 1
          IF (a(i) >= v) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2: DO
          j = j - 1
          IF ((a(j) <= v) .OR. (j == 1)) EXIT CNTLOOP2
        END DO CNTLOOP2
        t = a(i)
        a(i) = a(j)
        a(j) = t

        t_p = permutation(i)
        permutation(i) = permutation(j)
        permutation(j) = t_p

        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t

      permutation(j) = permutation(i)
      permutation(i) = permutation(r)
      permutation(r) = t_p
      CALL quicksort(a, permutation, l, i - 1)
      CALL quicksort(a, permutation, i + 1, r)
    END IF
  END SUBROUTINE quicksort_permutation_int

  SUBROUTINE swap_string(a, b)
    CHARACTER(LEN=*), INTENT(INOUT)  :: a, b !< strings for in-situ swap
    CHARACTER(LEN=LEN(a)) :: dummy
    dummy = a
    a = b
    b = dummy
  END SUBROUTINE swap_string

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_string(a, l_in, r_in)
    CHARACTER(LEN=*), INTENT(INOUT)           :: a(:) !< array for in-situ sorting
    INTEGER, INTENT(IN), OPTIONAL :: l_in, r_in !< left, right partition indices
    ! local variables
    INTEGER :: i, j, l, r, m
    CHARACTER(len=:), ALLOCATABLE :: v, t

    IF (PRESENT(l_in)) THEN
      l = l_in
    ELSE
      l = 1
    END IF
    IF (PRESENT(r_in)) THEN
      r = r_in
    ELSE
      r = SIZE(a, 1)
    END IF
    IF (r > l) THEN
      ALLOCATE (CHARACTER(LEN(a(1))) :: v, t)

      i = l - 1
      j = r

      ! median-of-three selection of partitioning element
      IF ((r - l) > 3) THEN
        m = (l + r)/2
        IF (a(l) > a(m)) CALL swap_string(a(l), a(m))
        IF (a(l) > a(r)) THEN
          CALL swap_string(a(l), a(r))
        ELSE IF (a(r) > a(m)) THEN
          CALL swap_string(a(r), a(m))
        END IF
      END IF

      v = a(r)
      LOOP: DO
        CNTLOOP1: DO
          i = i + 1
          IF (a(i) >= v) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2: DO
          j = j - 1
          IF ((a(j) <= v) .OR. (j == 1)) EXIT CNTLOOP2
        END DO CNTLOOP2
        t = a(i)
        a(i) = a(j)
        a(j) = t
        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t
      CALL quicksort(a, l, i - 1)
      CALL quicksort(a, i + 1, r)
      DEALLOCATE (v, t)
    END IF

  END SUBROUTINE quicksort_string

  SUBROUTINE insertion_sort_int(a)
    INTEGER, INTENT(inout) :: a(:)

    INTEGER :: t, h
    INTEGER :: i, n

    n = SIZE(a)

    DO i = 2, n
      t = a(i)
      DO h = i - 1, 1, -1
        IF (t >= a(h)) EXIT
        a(h + 1) = a(h)
      END DO
      a(h + 1) = t
    END DO
  END SUBROUTINE insertion_sort_int

END MODULE mo_util_sort
