!> Contains parameter for implementation of some specific model code
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
MODULE mo_jsb_impl_constants
 
  USE mo_kind,                  ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! general parameters
  INTEGER, PARAMETER :: SHORT_NAME_LEN = 10

  ! using the functionality of LOGICAL var with REAL var
  REAL(wp), PARAMETER :: false            = 0.0_wp    !< zero
  REAL(wp), PARAMETER :: true             = 1.0_wp    !< one
  REAL(wp), PARAMETER :: test_false_true  = 0.5_wp    !< 0.5, inbetween false and true; 'IF(var < test_false_true) -> FALSE'

  INTEGER,  PARAMETER :: ifs_nsoil         = 4                                  !< Number of soil levels in IFS initial file
  REAL(wp), PARAMETER :: ifs_soil_depth(4) = (/ 0.07_wp,0.28_wp,1._wp,2.89_wp/) !< Soil layer depths in IFS initial file [m]

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_impl_constants'

! CONTAINS



END MODULE mo_jsb_impl_constants
