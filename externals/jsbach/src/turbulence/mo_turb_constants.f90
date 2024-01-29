!> Contains constants for the turbulence processes
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
MODULE mo_turb_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! REAL(wp), PARAMETER ::          &
  !   !
  !   ! Parameters for roughness length
  !   & blending_height          = 100._wp,  & !< Blending height [m]
  !   & roughness_snow           = 0.001_wp, & !< Roughness length of snow-covered surfaces [m]
  !   & roughness_bare           = 0.005_wp, & !< Roughness length of bare land surfaces [m]
  !   & roughness_lai_saturation = 0.4_wp      !< factor in roughness length dependence on LAI (should be between 0.1 and 1)

#endif
END MODULE mo_turb_constants
