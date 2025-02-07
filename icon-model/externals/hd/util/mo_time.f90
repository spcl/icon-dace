! mo_time.f90 - Utilities for handling time information 
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_time

  !
  ! Authors:
  !
  ! S. Hagemann, HZG-IfK, june 2020, original source

!!  use mo_grid,         ONLY: 
!
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_timestep, leapyear, next_timestep    ! calc_date

  ! days per month
  INTEGER, DIMENSION(12) :: ndaymon  = [ 31,28,31, 30,31,30, 31,31,30, 31,30,31 ]


  !-----------------------------------------------------------------------------
CONTAINS
! *******************************************************************************
  SUBROUTINE next_timestep(ires, iday, imon, iyear, lendofyear)
! *******************************************************************************

    INTEGER, INTENT(in) :: ires     ! Temporal resolution: 1: daily, 2: monthly
    INTEGER, INTENT(inout) :: iday 
    INTEGER, INTENT(inout) :: imon
    INTEGER, INTENT(inout) :: iyear 
    LOGICAL, INTENT(out) :: lendofyear
!
    INTEGER  :: ndays                   ! days per year
!
    ndays = 365 ; ndaymon(2) = 28
    lendofyear = .FALSE.
    IF (leapyear(iyear)) THEN 
      ndays=366
      ndaymon(2) = 29
    ENDIF
    IF (ires.EQ.1) THEN
      IF (imon.eq.0) imon = 1
      iday = iday + 1
      IF (imon.EQ.12 .AND.iday.EQ.32) THEN
        iday = 1 ; imon = 1 
        iyear = iyear + 1
        WRITE(*,*) ' The year is now A.D. ', iyear
      ELSE IF (iday.GT.ndaymon(imon)) THEN
        iday = 1 ; imon = imon + 1
      ENDIF
      IF (imon.EQ.12 .AND. iday.EQ.31) lendofyear = .TRUE.
    ELSE IF (ires.EQ.2) THEN
      imon = imon + 1
      IF (imon.EQ.13) THEN
        imon = 1 
        iyear = iyear + 1
      ENDIF
      IF (imon.EQ.12) lendofyear = .TRUE.
    ENDIF

  END SUBROUTINE next_timestep

! *******************************************************************************
  SUBROUTINE calc_timestep(ires, iday, imon, iyear, jahr1, istep)
! *******************************************************************************

    INTEGER, INTENT(in) :: ires     ! Temporal resolution: 1: daily, 2: monthly
    INTEGER, INTENT(in) :: iday 
    INTEGER, INTENT(in) :: imon
    INTEGER, INTENT(in) :: iyear 
    INTEGER, INTENT(in) :: jahr1    ! Start Year
    INTEGER, INTENT(out) :: istep   ! Time step since 1.1. of jahr1
!
    INTEGER :: i
    INTEGER  :: ndays               ! days per year
!
    ndays = 365 ; ndaymon(2) = 28
    IF (leapyear(iyear)) THEN 
      ndays=366
      ndaymon(2) = 29
    ENDIF
    istep = 1
    IF (ires.EQ.1) THEN
      IF (iyear.GT.jahr1) THEN
        DO i=jahr1, iyear-1
          ndays=365
          IF (leapyear(i)) ndays=366
          istep = istep + ndays
        ENDDO
      ENDIF
      IF (imon.GT.1) THEN
        DO i=1, imon-1
          IF (leapyear(iyear)) ndaymon(2)=29
          istep = istep + ndaymon(i)
        ENDDO
      ENDIF
      IF (iday.GT.1) THEN
        istep = istep + iday-1
      ENDIF
    ELSE IF (ires.EQ.2) THEN
      IF (iyear.GT.jahr1) THEN
        istep = istep + (iyear-jahr1) * 12
      ENDIF
      IF (imon.GT.1) THEN
        istep = istep + imon -1
      ENDIF
    ENDIF

  END SUBROUTINE calc_timestep

! *******************************************************************************
  FUNCTION leapyear(iyear) RESULT (leap)
! *******************************************************************************
! *** checks year iyear whether it is a leap year
!
    INTEGER, INTENT(in)          :: iyear            ! Year
    LOGICAL :: leap

        leap=.FALSE.
        IF (MOD(iyear, 100) /= 0) THEN
          IF (MOD(iyear, 4) == 0) leap = .true.
        ELSE IF (MOD(iyear, 400) == 0) THEN
          leap = .true.
        ENDIF

  END FUNCTION leapyear

END MODULE mo_time

