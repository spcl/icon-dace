!> Contains constants for the natural disturbances processes
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
MODULE mo_disturbance_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC
  ! Submodel identifiers for subroutine disturbed_fract and relocate_disturbed_carbon
  ! to tell the subroutines which kind of parametrizations they should use.
  ! However, only two of them are used in JS4.
  !INTEGER , PARAMETER, PUBLIC  :: DIST_FIRE_WOOD        =  1 !! Calculation wrt. burning woody pfts
  !INTEGER , PARAMETER, PUBLIC  :: DIST_FIRE_GRASS       =  2 !! Calculation wrt. burning grass pfts
  INTEGER , PARAMETER, PUBLIC  :: DIST_FIRE             =  3 !! Calculation wrt. any burning pfts
  !INTEGER , PARAMETER, PUBLIC  :: DIST_WINDBREAK_WOOD   =  4 !! Calculation wrt. windbreak of woody types
  !INTEGER , PARAMETER, PUBLIC  :: DIST_WINDBREAK_GRASS  =  8 !! Calculation wrt. windbreak of grass types
  !                                                           !! (ok: 8 is unlikely to ever be used - maybe for extreme
  !                                                           !! precipitation events?!)
  INTEGER , PARAMETER, PUBLIC  :: DIST_WINDBREAK        = 12 !! Calculation wrt. windbreak of any pft types
  !INTEGER , PARAMETER, PUBLIC  :: FIRE_NONE             =  0
  INTEGER , PARAMETER, PUBLIC  :: FIRE_JSBACH           =  1
  !INTEGER , PARAMETER, PUBLIC  :: FIRE_ARORA            =  2
  !INTEGER , PARAMETER, PUBLIC  :: FIRE_THONICKE         =  3

  REAL(wp), PARAMETER          :: persist_rel_hum  = 0.95_wp   ! factor concerning the smoothing of relative air humidity (~14 days)
  REAL(wp), PARAMETER          :: persist_wind_10m = 0.9995_wp ! factor concerning the smoothing of daily maximum wind speed

#endif
END MODULE mo_disturbance_constants
