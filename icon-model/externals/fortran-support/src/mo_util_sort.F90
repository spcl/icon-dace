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
#ifdef __SX__
  PUBLIC :: radixsort, radixsort_int
#endif
  PUBLIC :: insertion_sort

  ! generic interface for in-situ QuickSort sorting routine
  INTERFACE quicksort
    MODULE PROCEDURE quicksort_real
    MODULE PROCEDURE quicksort_int
    MODULE PROCEDURE quicksort_permutation_int
    MODULE PROCEDURE quicksort_string
  END INTERFACE quicksort
#ifdef __SX__
  INTERFACE radixsort
    MODULE PROCEDURE radixsort_real
    MODULE PROCEDURE radixsort_int
  END INTERFACE
#endif
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

#ifdef __SX__
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Description:
  !!    Sorts an array of type dp with a radix sort using two bits at once.
  !!    Doubles the required memory, but saves up to 28% percent performance
  !!    when compared to the 1bt-version.
  !!    Scales linearly with the number of elements for large number of elements.
  !!    The difference between initial and final state is recorded by a permutation state perm
  !! Variables:
  !!    array: The dp-array to be sorted
  !!    perm: permutation state to record the sorting procedure
  !!          Must have the same size as array!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE radixsort_real(array, perm)

    IMPLICIT NONE

    REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: perm
    INTEGER(KIND=wp), DIMENSION(SIZE(array)) :: iarray

    INTEGER(KIND=wp), DIMENSION(SIZE(array)) :: bucket00, bucket01, bucket10, bucket11
    INTEGER, DIMENSION(SIZE(array)) :: pbucket00, pbucket01, pbucket10, pbucket11
    INTEGER :: idx00, idx01, idx10, idx11

    INTEGER(KIND=wp), PARAMETER :: bitmaskr2i = -1
    ! REAL(KIND=wp), PARAMETER :: bitmaski2r = Z'FFFFFFFFFFFFFFFF'

    INTEGER :: i, j
    INTEGER :: n, offs
    INTEGER(KIND=wp) :: tmp

    n = SIZE(array)

    ! Transfer the bit pattern of the real array to an i8 array
    iarray = TRANSFER(array, bitmaskr2i, n)

    ! Loop over the bits
    DO i = 0, STORAGE_SIZE(tmp) - 2, 2
      idx00 = 0
      idx01 = 0
      idx10 = 0
      idx11 = 0
      ! sort numbers into buckets based on their bit i
      DO j = 1, n
        ! extract the i-th bit
        tmp = IBITS(iarray(j), i, 2)
        IF (tmp < 2) THEN
          IF (tmp == 0) THEN
            !increase the counter of the selected bucket
            idx00 = idx00 + 1
            ! add number to bucket
            bucket00(idx00) = iarray(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket00(idx00) = perm(j)
          ELSE
            !increase the counter of the selected bucket
            idx01 = idx01 + 1
            ! add number to bucket
            bucket01(idx01) = iarray(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket01(idx01) = perm(j)
          END IF
        ELSE
          IF (tmp == 2) THEN
            !increase the counter of the selected bucket
            idx10 = idx10 + 1
            ! add number to bucket
            bucket10(idx10) = iarray(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket10(idx10) = perm(j)
          ELSE
            !increase the counter of the selected bucket
            idx11 = idx11 + 1
            ! add number to bucket
            bucket11(idx11) = iarray(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket11(idx11) = perm(j)
          END IF
        END IF
      END DO

      ! copy the presorted numbers back onto the array
      DO j = 1, idx00
        iarray(j) = bucket00(j)
        IF (PRESENT(perm)) perm(j) = pbucket00(j)
      END DO
      offs = idx00
      DO j = 1, idx01
        iarray(offs + j) = bucket01(j)
        IF (PRESENT(perm)) perm(offs + j) = pbucket01(j)
      END DO
      offs = offs + idx01
      DO j = 1, idx10
        iarray(offs + j) = bucket10(j)
        IF (PRESENT(perm)) perm(offs + j) = pbucket10(j)
      END DO
      offs = offs + idx10
      DO j = 1, idx11
        iarray(offs + j) = bucket11(j)
        IF (PRESENT(perm)) perm(offs + j) = pbucket11(j)
      END DO
    END DO

    ! sort by sign
    idx00 = 0
    idx01 = 0
    idx10 = 0
    idx11 = 0
    i = STORAGE_SIZE(tmp) - 2
    ! sort numbers into buckets based on their first bit and
    DO j = 1, n
      tmp = IBITS(iarray(j), i, 2)
      IF (tmp < 2) THEN
        IF (tmp == 0) THEN
          !increase the counter of the selected bucket
          idx00 = idx00 + 1
          ! add number to bucket
          bucket00(idx00) = iarray(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket00(idx00) = perm(j)
        ELSE
          !increase the counter of the selected bucket
          idx01 = idx01 + 1
          ! add number to bucket
          bucket01(idx01) = iarray(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket01(idx01) = perm(j)
        END IF
      ELSE
        IF (tmp == 2) THEN
          !increase the counter of the selected bucket
          idx10 = idx10 + 1
          ! add number to bucket
          bucket10(idx10) = iarray(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket10(idx10) = perm(j)
        ELSE
          !increase the counter of the selected bucket
          idx11 = idx11 + 1
          ! add number to bucket
          bucket11(idx11) = iarray(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket11(idx11) = perm(j)
        END IF
      END IF
    END DO

    ! copy the presorted numbers back onto the array
    ! the half inverse order is a result of the ieee 745 float sign bit standard
    DO j = 1, idx11
      iarray(idx11 - j + 1) = bucket11(j)
      IF (PRESENT(perm)) perm(idx11 - j + 1) = pbucket11(j)
    END DO
    offs = idx11
    DO j = 1, idx10
      iarray(offs + idx10 - j + 1) = bucket10(j)
      IF (PRESENT(perm)) perm(offs + idx10 - j + 1) = pbucket10(j)
    END DO
    offs = offs + idx10
    DO j = 1, idx00
      iarray(offs + j) = bucket00(j)
      IF (PRESENT(perm)) perm(offs + j) = pbucket00(j)
    END DO
    offs = offs + idx00
    DO j = 1, idx01
      iarray(offs + j) = bucket01(j)
      IF (PRESENT(perm)) perm(offs + j) = pbucket01(j)
    END DO

    ! Transfer the bit pattern of the i8 array to the r8 array
    array = TRANSFER(iarray, Z'FFFFFFFFFFFFFFFF', n)

    RETURN
  END SUBROUTINE radixsort_real
#endif
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
#ifdef __SX__
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Description:
  !!    Sorts an array of type int with a radix sort using two bits at once.
  !!    Doubles the required memory, but saves up to 28% percent performance
  !!    when compared to the 1bt-version.
  !!    Scales linearly with the number of elements for large number of elements.
  !!    The difference between initial and final state is recorded by a permutation state perm
  !! Variables:
  !!    array: The int-array to be sorted
  !!    perm: permutation state to record the sorting procedure
  !!          Must have the same size as array!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE radixsort_int(array, perm)

    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: perm

    INTEGER, DIMENSION(SIZE(array)) :: bucket00, bucket01, bucket10, bucket11
    INTEGER, DIMENSION(SIZE(array)) :: pbucket00, pbucket01, pbucket10, pbucket11
    INTEGER :: idx00, idx01, idx10, idx11

    INTEGER :: i, j
    INTEGER :: n, offs
    INTEGER :: tmp

    n = SIZE(array)

    ! Loop over the bits
    DO i = 0, STORAGE_SIZE(tmp) - 1, 2
      idx00 = 0
      idx01 = 0
      idx10 = 0
      idx11 = 0
      ! sort numbers into buckets based on their bit i
      DO j = 1, n
        tmp = IBITS(array(j), i, 2)
        IF (tmp < 2) THEN
          IF (tmp == 0) THEN
            !increase the counter of the selected bucket
            idx00 = idx00 + 1
            ! add number to bucket
            bucket00(idx00) = array(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket00(idx00) = perm(j)
          ELSE
            !increase the counter of the selected bucket
            idx01 = idx01 + 1
            ! add number to bucket
            bucket01(idx01) = array(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket01(idx01) = perm(j)
          END IF
        ELSE
          IF (tmp == 2) THEN
            !increase the counter of the selected bucket
            idx10 = idx10 + 1
            ! add number to bucket
            bucket10(idx10) = array(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket10(idx10) = perm(j)
          ELSE
            !increase the counter of the selected bucket
            idx11 = idx11 + 1
            ! add number to bucket
            bucket11(idx11) = array(j)
            ! same operation on the permutation
            IF (PRESENT(perm)) pbucket11(idx11) = perm(j)
          END IF
        END IF
      END DO

      ! copy the presorted numbers back onto the array
      DO j = 1, idx00
        array(j) = bucket00(j)
        IF (PRESENT(perm)) perm(j) = pbucket00(j)
      END DO
      offs = idx00
      DO j = 1, idx01
        array(offs + j) = bucket01(j)
        IF (PRESENT(perm)) perm(offs + j) = pbucket01(j)
      END DO
      offs = offs + idx01
      DO j = 1, idx10
        array(offs + j) = bucket10(j)
        IF (PRESENT(perm)) perm(offs + j) = pbucket10(j)
      END DO
      offs = offs + idx10
      DO j = 1, idx11
        array(offs + j) = bucket11(j)
        IF (PRESENT(perm)) perm(offs + j) = pbucket11(j)
      END DO
    END DO

    ! sort by sign
    idx00 = 0
    idx01 = 0
    idx10 = 0
    idx11 = 0
    i = STORAGE_SIZE(tmp) - 2
    ! sort numbers into buckets based on their first bit
    DO j = 1, n
      tmp = IBITS(array(j), i, 2)
      IF (tmp < 2) THEN
        IF (tmp == 0) THEN
          !increase the counter of the selected bucket
          idx00 = idx00 + 1
          ! add number to bucket
          bucket00(idx00) = array(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket00(idx00) = perm(j)
        ELSE
          !increase the counter of the selected bucket
          idx01 = idx01 + 1
          ! add number to bucket
          bucket01(idx01) = array(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket01(idx01) = perm(j)
        END IF
      ELSE
        IF (tmp == 2) THEN
          !increase the counter of the selected bucket
          idx10 = idx10 + 1
          ! add number to bucket
          bucket10(idx10) = array(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket10(idx10) = perm(j)
        ELSE
          !increase the counter of the selected bucket
          idx11 = idx11 + 1
          ! add number to bucket
          bucket11(idx11) = array(j)
          ! same operation on the permutation
          IF (PRESENT(perm)) pbucket11(idx11) = perm(j)
        END IF
      END IF
    END DO

    ! copy the presorted numbers back onto the array
    ! the half inverse order is a result of the ieee 745 float sign bit standard
    DO j = 1, idx10
      array(j) = bucket10(j)
      IF (PRESENT(perm)) perm(j) = pbucket10(j)
    END DO
    offs = idx10
    DO j = 1, idx11
      array(offs + j) = bucket11(j)
      IF (PRESENT(perm)) perm(offs + j) = pbucket11(j)
    END DO
    offs = offs + idx11
    DO j = 1, idx00
      array(offs + j) = bucket00(j)
      IF (PRESENT(perm)) perm(offs + j) = pbucket00(j)
    END DO
    offs = offs + idx00
    DO j = 1, idx01
      array(offs + j) = bucket01(j)
      IF (PRESENT(perm)) perm(offs + j) = pbucket01(j)
    END DO

    RETURN
  END SUBROUTINE radixsort_int
#endif
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
    INTEGER, INTENT(INOUT) :: a(:)

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
