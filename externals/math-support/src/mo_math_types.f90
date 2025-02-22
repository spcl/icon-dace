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

!>
!! Contains basic math types
!!

MODULE mo_math_types
  USE ISO_C_BINDING, ONLY: c_int64_t
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64, &
                                                                              sp => real32, &
                                                                              dp => real64

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_cartesian_coordinates
  PUBLIC :: t_geographical_coordinates
  PUBLIC :: t_line
  PUBLIC :: t_tangent_vectors
  PUBLIC :: t_Statistics

  ! cartesian coordinate class
  TYPE t_cartesian_coordinates
    REAL(wp) :: x(3)
  END TYPE t_cartesian_coordinates

  ! geographical coordinate class
  TYPE t_geographical_coordinates
    REAL(wp) :: lon
    REAL(wp) :: lat
  END TYPE t_geographical_coordinates

  ! the two coordinates on the tangent plane
  TYPE t_tangent_vectors
    REAL(wp) :: v1
    REAL(wp) :: v2
  END TYPE t_tangent_vectors

  ! line class
  TYPE t_line
    TYPE(t_geographical_coordinates) :: p1
    TYPE(t_geographical_coordinates) :: p2
  END TYPE t_line

  TYPE :: t_Statistics
    INTEGER(c_int64_t) :: sampleCount
    REAL(wp) :: MIN, mean, MAX
  CONTAINS
    PROCEDURE :: reset => statistics_reset
    GENERIC :: add => addData_s1d, addData_d1d, addStatistics
    ! scan a given array AND update the statistics accordingly
    PROCEDURE :: addData_s1d => statistics_addData_s1d
    PROCEDURE :: addData_d1d => statistics_addData_d1d
    PROCEDURE :: addStatistics => statistics_addStatistics ! update the statistics with the contents of another t_Statistics object
  END TYPE t_Statistics

CONTAINS

  SUBROUTINE statistics_reset(me)
    CLASS(t_Statistics), INTENT(INOUT) :: me

    me%sampleCount = 0_c_int64_t
    me%MIN = HUGE(me%MIN)
    me%mean = 0.0_wp
    me%MAX = -HUGE(me%MAX)
  END SUBROUTINE statistics_reset

  SUBROUTINE statistics_addData_s1d(me, DATA)
    CLASS(t_Statistics), INTENT(INOUT) :: me
    REAL(sp), INTENT(IN) :: DATA(:)
    INTEGER :: i, icount
    REAL(wp) :: data_sum, data_max, data_min, data_wp

    TYPE(t_Statistics) :: newStatistics

    CALL newStatistics%reset()

    icount = 0
    data_sum = 0._wp
    data_max = -HUGE(DATA)
    data_min = HUGE(DATA)
!$OMP PARALLEL DO PRIVATE(data_wp), REDUCTION(+:data_sum,icount), REDUCTION(MAX:data_max), REDUCTION(MIN:data_min)
    DO i = 1, SIZE(DATA)
      icount = icount + 1
      data_wp = REAL(DATA(i), wp)
      data_sum = data_sum + data_wp
      data_max = MAX(data_max, data_wp)
      data_min = MIN(data_min, data_wp)
    END DO
!$OMP END PARALLEL DO
    newStatistics%sampleCount = icount
    newStatistics%MIN = data_min
    newStatistics%MAX = data_max
    IF (icount > 0) THEN
      newStatistics%mean = data_sum/REAL(icount, wp)
    END IF
    CALL me%add(newStatistics)
  END SUBROUTINE statistics_addData_s1d

  SUBROUTINE statistics_addData_d1d(me, DATA)
    CLASS(t_Statistics), INTENT(INOUT) :: me
    REAL(dp), INTENT(IN) :: DATA(:)
    INTEGER :: i, icount
    REAL(wp) :: data_sum, data_max, data_min

    TYPE(t_Statistics) :: newStatistics

    CALL newStatistics%reset()

    icount = 0
    data_sum = 0._wp
    data_max = -HUGE(DATA)
    data_min = HUGE(DATA)
!$OMP PARALLEL DO REDUCTION(+:data_sum,icount), REDUCTION(MAX:data_max), REDUCTION(MIN:data_min)
    DO i = 1, SIZE(DATA)
      icount = icount + 1
      data_sum = data_sum + DATA(i)
      data_max = MAX(data_max, DATA(i))
      data_min = MIN(data_min, DATA(i))
    END DO
!$OMP END PARALLEL DO
    newStatistics%sampleCount = icount
    newStatistics%MIN = data_min
    newStatistics%MAX = data_max
    IF (icount > 0) THEN
      newStatistics%mean = data_sum/REAL(icount, wp)
    END IF

    CALL me%add(newStatistics)
  END SUBROUTINE statistics_addData_d1d

  SUBROUTINE statistics_addStatistics(me, other)
    CLASS(t_Statistics), INTENT(INOUT) :: me
    CLASS(t_Statistics), INTENT(IN) :: other

    INTEGER(c_int64_t) :: newSampleCount

    newSampleCount = me%sampleCount + other%sampleCount
    me%MIN = MIN(me%MIN, other%MIN)
    me%mean = (me%mean*REAL(me%sampleCount, wp) + other%mean*REAL(other%sampleCount, wp))/REAL(newSampleCount, wp)
    me%MAX = MAX(me%MAX, other%MAX)
    me%sampleCount = newSampleCount
  END SUBROUTINE statistics_addStatistics

END MODULE mo_math_types
