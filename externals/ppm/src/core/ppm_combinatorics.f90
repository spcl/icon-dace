!> @file ppm_combinatorics.f90
!! @brief simple functions and transformations from combinatorics
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords:
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!> gathers some base routines for combinatorial problems
#include "fc_feature_defs.inc"
MODULE ppm_combinatorics
  USE ppm_base, ONLY: assertion
  USE ppm_std_type_kinds, ONLY: i4, i8
  USE ppm_compare, ONLY: cmp_i4, rcmp_i4
  USE ppm_random, ONLY: irand, irandr
  USE ppm_extents, ONLY: iinterval, extent_size, extent_start
  USE ppm_sort, ONLY: sorted_insertion, qsort_i4
  IMPLICIT NONE
  PRIVATE
  !> randomly permute a given array
  INTERFACE permute
    MODULE PROCEDURE permute_randomly_i4
    MODULE PROCEDURE permute_randomly_i8
  END INTERFACE permute

  !> produce random selection
  INTERFACE selection
    MODULE PROCEDURE random_selection_i4
  END INTERFACE selection

  !> produce random selection and its complement
  INTERFACE combination
    MODULE PROCEDURE random_combination_i4
  END INTERFACE combination

  !> establish whether an array is a permutation of another array or a range
  INTERFACE is_permutation
    MODULE PROCEDURE is_permutation_of_range_i4
    MODULE PROCEDURE is_permutation_of_array_i4
  END INTERFACE is_permutation

  !> compute list of prime factors (sorted, but duplicates included)
  INTERFACE prime_factorization
    MODULE PROCEDURE prime_factorization_i4
  END INTERFACE prime_factorization

  PUBLIC :: permute, is_permutation, selection, combination
  PUBLIC :: prime_factorization

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_combinatorics.f90'
CONTAINS
  !> randomly permute a given array
  !! @param a array to permute
  SUBROUTINE permute_randomly_i4(a)
    INTEGER(i4), INTENT(inout) :: a(:)

    INTEGER :: n, i, s, left
    INTEGER(i4) :: temp

    n = SIZE(a)

    left = n
    DO i = 1, n-1
      temp = a(i)
      s = MODULO(irand(), left) + i
      a(i) = a(s)
      a(s) = temp
      left = left - 1
    END DO
  END SUBROUTINE permute_randomly_i4

  !> randomly permute a given array
  !! @param a array to permute
  SUBROUTINE permute_randomly_i8(a)
    INTEGER(i8), INTENT(inout) :: a(:)

    INTEGER :: n, i, s, left
    INTEGER(i8) :: temp

    n = SIZE(a)

    left = n
    DO i = 1, n
      temp = a(i)
      s = MODULO(irand(), left) + i
      a(i) = a(s)
      a(s) = temp
      left = left - 1
    END DO
  END SUBROUTINE permute_randomly_i8

  !> establish whether one array of integers is a permutation of a range
  !! @param a array
  !! @param r range
  !! @return \c .TRUE. if \c a is permutation of \c r, \c .FALSE. if not
  PURE FUNCTION is_permutation_of_range_i4(a, r) RESULT(p)
    INTEGER(i4), INTENT(in) :: a(:)
    TYPE(iinterval), INTENT(in) :: r
    LOGICAL :: p

    INTEGER :: i, n
    INTEGER(i4) :: lb, ub, ofs
    LOGICAL, ALLOCATABLE :: seen_mask(:)

    p = .TRUE.
    n = SIZE(a)
    ALLOCATE(seen_mask(n))
    seen_mask = .FALSE.
    lb = r%first
    ub = r%last
    DO i = 1, n
      ofs = a(i) - lb + 1
      IF (a(i) < lb .OR. a(i) > ub .OR. seen_mask(ofs)) THEN
        p = .FALSE.
        EXIT
      END IF
      seen_mask(ofs) = .TRUE.
    END DO
    DEALLOCATE(seen_mask)
  END FUNCTION is_permutation_of_range_i4

  !> establish whether one array of integers is a permutation of another
  !! @param a first array
  !! @param b second array
  !! @return \c .TRUE. if \c a is permutation of \c b
#ifndef __G95__
  PURE &
#endif
       FUNCTION is_permutation_of_array_i4(a, b) RESULT(p)
    LOGICAL :: p
    INTEGER(i4), INTENT(in) :: a(:), b(:)

    INTEGER(i4), ALLOCATABLE :: sorted(:,:)
    INTEGER :: n

    p = .FALSE.
    n = SIZE(a)
    IF (n == SIZE(b)) THEN
      ALLOCATE(sorted(n, 2))
      sorted(:, 1) = a
      sorted(:, 2) = b
      CALL qsort_i4(sorted(:,1), cmp_i4)
      CALL qsort_i4(sorted(:,2), cmp_i4)
      p = ALL(sorted(:, 1) == sorted(:, 2))
      DEALLOCATE(sorted)
    END IF
  END FUNCTION is_permutation_of_array_i4

  !> produce random selection from range
  !!
  !! normally this procedure is only efficient when SIZE(selected)
  !! is significantly smaller than extent_size(range), when SIZE(selected)
  !! approaches extent_size(range) it's cheaper to produce a permutation
  !! and use the first SIZE(selected) elements
  !! @param[out] selected array to fill with random selection from range,
  !! where SIZE(selected) <= extent_size(range)
  !! @param range set of integers to selected from
  SUBROUTINE random_selection_i4(selected, range)
    INTEGER(i4), INTENT(out) :: selected(:)
    TYPE(iinterval), INTENT(in) :: range

    INTEGER :: i, j, n, t
    TYPE(iinterval) :: range_left

    n = SIZE(selected)
    CALL assertion(n < extent_size(range), filename, __LINE__, &
         "SIZE(selected) > extent_size(range)")

    range_left = range
    DO i = 1, n
      t = irandr(range_left)
      DO j = i - 1, 1, -1
        IF (selected(j) <= t) THEN
          t = t + 1
        ELSE
          EXIT
        END IF
      END DO
      selected(i) = t
      range_left%last = range_left%last - 1
      CALL sorted_insertion(selected, i, 4, 0, rcmp_i4)
    END DO

  END SUBROUTINE random_selection_i4

  !> produce random selection and its complement from range
  !! @param range range to select from
  !! @param[out] selected array where a random selection of integers
  !! from range are written to
  !! @param[out] not_selected integers from range not written to
  !! selected are stored here
  SUBROUTINE random_combination_i4(selected, not_selected, range)
    INTEGER(i4), INTENT(out) :: selected(:), not_selected(:)
    TYPE(iinterval), INTENT(in) :: range

    INTEGER :: i, j, n, ssize, nssize, left, ofs, s, t

    n = extent_size(range)
    ssize = SIZE(selected)
    nssize = SIZE(not_selected)
    CALL assertion(ssize + nssize == n, filename, __LINE__, &
         "selected and not_selected must have complementary size")

    ofs = extent_start(range)
    DO i = 1, ssize
      selected(i) = i + ofs - 1
    END DO
    DO i = 1, nssize
      not_selected(i) = i + ssize + ofs - 1
    END DO
    left = n
    DO i = 1, ssize
      s = irandr(iinterval(1, left)) + i - 1
      IF (s > ssize) THEN
        j = s - ssize
        t = not_selected(j)
        not_selected(j) = selected(i)
        selected(i) = t
      ELSE
        t = selected(s)
        selected(s) = selected(i)
        selected(i) = t
      END IF
      left = left - 1
    END DO
  END SUBROUTINE random_combination_i4

  !> compute list of prime factors (sorted, but duplicates included)
  !! @param n integer to factor
  !! @param factors integers, where factor(i)
  !! is prime and PRODUCT(factors)==n, factors is allocated to the
  !! necessary size
  PURE SUBROUTINE prime_factorization_i4(n, factors)
    INTEGER(i4), INTENT(in) :: n
    INTEGER(i4), ALLOCATABLE, INTENT(out) :: factors(:)

    ! a 32 bit integer cannot have more than 32 factors
    INTEGER(i4) :: temp(31), num_factors
    INTERFACE
      PURE SUBROUTINE ppm_prime_factorization_i4(n, num_factors, factors)
        USE ppm_std_type_kinds, ONLY: i4
        INTEGER(i4), INTENT(in) :: n
        INTEGER(i4), INTENT(out) :: num_factors, factors(31)
      END SUBROUTINE ppm_prime_factorization_i4
    END INTERFACE

    CALL ppm_prime_factorization_i4(n, num_factors, temp)
    ALLOCATE(factors(num_factors))
    factors = temp(1:num_factors)
  END SUBROUTINE prime_factorization_i4

END MODULE ppm_combinatorics
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
