! mo_array_utils.f90 - Routines to work with structured data
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

!> @file mo_array_utils.f90
!> @author Thomas Jahns <jahns@dkrz.de>
!> @brief general routines to work with array structured data
!>
!> contains helper routines useful in parallel environments
!>
!> Mai 2011 - Veronika Gayler: added function inc_monotonic_interval_closest_midpoint
!>            using the example of dec_monotonic_interval_closest_midpoint (T. Jahns)
!> Mai 2011 - Veronika Gayler: renamed functions to meet the fortran standard of 31 
!>            characters maximum
!>
!> This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
!> Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
!> Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
!> version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
!> doi: 10.1029/2018MS001400.

MODULE mo_array_utils

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE
                                              !   former names (more than 31 char)
  PUBLIC :: inc_regular_bisection, &          ! inc_regular_interval_bisection
            inc_monotonic_bisection, &        ! inc_monotonic_interval_bisection
            inc_monotonic_closest_midpoint, & ! inc_monotonic_interval_closest_midpoint
            dec_monotonic_closest_midpoint    ! dec_monotonic_interval_closest_midpoint

CONTAINS

  !> search a for index i such that x element of [ a(i), a(i + 1) )
  !> a satisfies a(i) <= a(i + 1)
  PURE FUNCTION inc_monotonic_bisection(a, x, first_k) RESULT(i)
    REAL(dp), INTENT(in)           :: a(:), x
    INTEGER,  INTENT(in), OPTIONAL :: first_k

    INTEGER :: i, j, k, n

    n = SIZE(a)
    i = 1
    j = SIZE(a)
    IF (PRESENT(first_k)) THEN
      k = first_k
      IF (k > n .OR. k < 1) THEN
        ! error: invalid inital k, needs to be in [sidx,eidx]
        i = - 1
        RETURN
      END IF
    ELSE
      k = (i + j)/2
    END IF
    DO
      IF (x < a(k)) THEN
        j = k
      ELSE
        i = k
      END IF
      IF (i + 1 >= j) EXIT
      k = (i + j)/2
    END DO
  END FUNCTION inc_monotonic_bisection

  PURE FUNCTION inc_monotonic_bisection_slice(sidx, eidx, a, x, first_k) RESULT(i)
    INTEGER,  INTENT(in)           :: sidx, eidx
    REAL(dp), INTENT(in)           :: a(sidx:eidx), x
    INTEGER,  INTENT(in), OPTIONAL :: first_k

    INTEGER :: i

    IF (PRESENT(first_k)) THEN
      i = sidx - 1 + &
           inc_monotonic_bisection(a(sidx:eidx), x, &
           first_k = first_k - sidx + 1)
    ELSE
      i = sidx - 1 + &
           inc_monotonic_bisection(a(sidx:eidx), x)
    END IF
  END FUNCTION inc_monotonic_bisection_slice

  PURE FUNCTION inc_regular_bisection_slice(sidx, eidx, a, x) RESULT(i)
    INTEGER,  INTENT(in) :: sidx, eidx
    REAL(dp), INTENT(in) :: a(sidx:eidx), x

    INTEGER :: i

    i = sidx - 1 + inc_regular_bisection(a(sidx:eidx), x)
  END FUNCTION inc_regular_bisection_slice

  PURE FUNCTION inc_regular_bisection(a, x) RESULT(i)
    REAL(dp), INTENT(in) :: a(:), x

    INTEGER :: first_k, n, i

    n = SIZE(a)
    first_k = MIN(INT((x - a(1))/(a(n) - a(1)) * n), &
         n)
    i = inc_monotonic_bisection(a, x, first_k)
  END FUNCTION inc_regular_bisection

  !> search for i such that (a(i - 1) + a(i))/2 < x < (a(i) + a(i + 1))/2
  !> a satisfies a(i) <= a(i + 1)
  !> also considers alb <= x < (a(1) + a(2))/2
  !> and (a(eidx - 1) + a(eidx))/2 < x < (a(eidx) + aub)/2
  !> if alb or aub are given respectively
  PURE FUNCTION inc_monotonic_closest_midpoint(a, x, alb, aub, first_i) RESULT(i)
    REAL(dp), INTENT(in)           :: a(:), x
    REAL(dp), INTENT(in), OPTIONAL :: alb, aub
    INTEGER,  INTENT(in), OPTIONAL :: first_i

    INTEGER :: i, m, n

    n = SIZE(a)
    IF (n < 2) THEN
      i = -1
      RETURN
    END IF
    IF (PRESENT(alb)) THEN
      IF (alb <= x .AND. x < (a(1) + a(2)) * 0.5_dp) THEN
        i = 1
        RETURN
      ELSE IF(x < (a(1) + a(2)) * 0.5_dp) THEN
        i = -1
        RETURN
      END IF
    ELSE IF(x < (a(1) + a(2)) * 0.5_dp) THEN
      i = -1
      RETURN
    END IF
    IF (PRESENT(aub)) THEN
      IF ((a(n - 1) + a(n)) * 0.5_dp <= x .AND. x < aub) THEN
        i = n
        RETURN
      ELSE IF((a(n - 1) + a(n)) * 0.5_dp <= x) THEN
        i = -n
        RETURN
      END IF
    ELSE IF((a(n - 1) + a(n)) * 0.5_dp <= x) THEN
      i = -n
      RETURN
    END IF
    IF (n < 3) THEN
      i =  -1
      RETURN
    END IF
    ! at this point the following holds:
    ! 1. (a(1) + a(2))*0.5_dp < x and (a(n - 1) + a(n))*0.5_dp >= x
    ! 2. n > 2
    ! therefore it is guaranteed, we can find an index satisfying
    ! above condition
    m = 1
    i = (n + 1)/ 2
    IF (PRESENT(first_i)) THEN
      IF (first_i < n .AND. first_i > 1) i = first_i
    END IF
    DO
      IF ((a(i - 1) + a(i))*0.5_dp <= x .AND. x < (a(i) + a(i + 1))/2) THEN
        RETURN
      ELSE IF ((a(i - 1) + a(i))*0.5_dp < x) THEN
        m = i
      ELSE IF (x <= (a(i) + a(i + 1))*0.5_dp) THEN
        n = i
      END IF
      i = (m + n + 1)/2
    END DO

  END FUNCTION inc_monotonic_closest_midpoint

  !> search for i such that (a(i - 1) + a(i))/2 < x < (a(i) + a(i + 1))/2
  !> a satisfies a(i) <= a(i + 1)
  !> also considers (alb + a(sidx))/2 < x < (a(sidx) + a(sidx + 1))/2
  !> and (a(eidx - 1) + a(eidx))/2 < x < (a(eidx) + aub)/2
  !> if alb or aub are given respectively
  PURE FUNCTION inc_monotonic_closest_mid_slice(sidx, eidx, a, x, alb, aub, first_i) RESULT(i)
    INTEGER,  INTENT(in) :: sidx, eidx
    REAL(dp), INTENT(in) :: a(sidx:eidx), x
    REAL(dp), INTENT(in), OPTIONAL :: alb, aub
    INTEGER,  INTENT(in), OPTIONAL :: first_i

    INTEGER :: i

    i = -1
  END FUNCTION inc_monotonic_closest_mid_slice

  !> search for i such that (a(i - 1) + a(i))/2 >= x > (a(i) + a(i + 1))/2
  !> a satisfies a(i) >= a(i + 1)
  !> also considers aub >= x > (a(1) + a(2))/2
  !> and (a(size(a) - 1) + a(size(a)))/2 >= x > alb
  !> if alb or aub are given respectively
  !> returns -1 if no index satisfies the condition
  PURE FUNCTION dec_monotonic_closest_midpoint(a, x, alb, aub, first_i) RESULT(i)
    REAL(dp), INTENT(in)           :: a(:), x
    REAL(dp), INTENT(in), OPTIONAL :: alb, aub
    INTEGER,  INTENT(in), OPTIONAL :: first_i

    INTEGER :: i, m, n

    n = SIZE(a)
    IF (n < 2) THEN
      i = -1
      RETURN
    END IF
    IF (PRESENT(aub)) THEN
      IF (aub >= x .AND. x > (a(1) + a(2)) * 0.5_dp) THEN
        i = 1
        RETURN
      ELSE IF(x > (a(1) + a(2)) * 0.5_dp) THEN
        i = -1
        RETURN
      END IF
    ELSE IF(x > (a(1) + a(2)) * 0.5_dp) THEN
      i = -1
      RETURN
    END IF
    IF (PRESENT(alb)) THEN
      IF ((a(n - 1) + a(n)) * 0.5_dp >= x .AND. x > alb) THEN
        i = n
        RETURN
      ELSE IF((a(n - 1) + a(n)) * 0.5_dp >= x) THEN
        i = -n
        RETURN
      END IF
    ELSE IF((a(n - 1) + a(n)) * 0.5_dp >= x) THEN
      i = -n
      RETURN
    END IF
    IF (n < 3) THEN
      i =  -1
      RETURN
    END IF
    ! at this point the following holds:
    ! 1. (a(1) + a(2))*0.5_dp < x and (a(n - 1) + a(n))*0.5_dp >= x
    ! 2. n > 2
    ! therefore it is guaranteed, we can find an index satisfying
    ! above condition
    m = 1
    i = (n + 1)/ 2
    IF (PRESENT(first_i)) THEN
      IF (first_i < n .AND. first_i > 1) i = first_i
    END IF
    DO
      IF ((a(i - 1) + a(i))*0.5_dp < x) THEN
        n = i
      ELSE IF (x <= (a(i) + a(i + 1))*0.5_dp) THEN
        m = i
      ELSE ! IF ((a(i - 1) + a(i))*0.5_dp >= x .AND. x > (a(i) + a(i + 1))/2)
        RETURN
      END IF
      i = (m + n + 1)/2
    END DO
  END FUNCTION dec_monotonic_closest_midpoint


  !> search for i such that (a(i - 1) + a(i))/2 >= x > (a(i) + a(i + 1))/2
  !> a satisfies a(i) >= a(i + 1)
  !> also considers aub >= x > (a(sidx) + a(sidx + 1))/2
  !> and (a(eidx - 1) + a(eidx))/2 >= x > alb
  !> if alb or aub are given respectively
  PURE FUNCTION dec_monotonic_closest_mid_slice(sidx, eidx, a, x, alb, aub, first_i) RESULT(i)
    INTEGER,  INTENT(in)           :: sidx, eidx
    REAL(dp), INTENT(in)           :: a(sidx:eidx), x
    REAL(dp), INTENT(in), OPTIONAL :: alb, aub
    INTEGER,  INTENT(in), OPTIONAL :: first_i

    INTEGER :: i

    IF (PRESENT(first_i)) THEN
      i = sidx - 1 + &
           dec_monotonic_closest_midpoint(a(sidx:eidx), x, &
           aub = aub, alb = alb, first_i = first_i - sidx + 1)
    ELSE
      i = sidx - 1 + &
           dec_monotonic_closest_midpoint(a(sidx:eidx), x, &
           aub = aub, alb = alb)
    END IF
  END FUNCTION dec_monotonic_closest_mid_slice

END MODULE mo_array_utils
