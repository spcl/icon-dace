!> Contains mathematical constants
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_math_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_math_constants,        ONLY: pi, deg2rad

  IMPLICIT NONE
  PUBLIC


  ! QUINCY simplified time control
  REAL(wp), PARAMETER   :: one_year  = 365._wp                !< length of one year in days
  REAL(wp), PARAMETER   :: one_day   = 86400._wp              !< length of one day in seconds
  ! conversion of degree to radians
  REAL(wp), PARAMETER   :: degcircle = 360._wp                !< degrees per circle

  REAL(wp), PARAMETER   :: hourday         = 24.0_wp   !< hours per day
  REAL(wp), PARAMETER   :: seconds_per_day = 86400._wp !< seconds per day

  ! QUINCY math helpers
  REAL(wp), PARAMETER   :: eps12 = 0.000000000001_wp
  REAL(wp), PARAMETER   :: eps8  = 0.00000001_wp
  REAL(wp), PARAMETER   :: eps4  = 0.0001_wp
  REAL(wp), PARAMETER   :: eps1  = 0.1_wp
  REAL(wp), PARAMETER   :: zero  = 0.0_wp           !< just zero


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_math_constants'

!CONTAINS

#endif
END MODULE mo_jsb_math_constants
