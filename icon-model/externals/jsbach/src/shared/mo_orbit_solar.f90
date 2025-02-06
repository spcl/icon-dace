!> Contains some functions for orbital and solar parameters
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
MODULE mo_orbit_solar
#ifndef __NO_JSBACH__

  USE mo_jsb_orbit_solar_iface, ONLY: inquire_declination, compute_cos_zenith_angle

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: inquire_declination, compute_cos_zenith_angle

  CHARACTER(len=*), PARAMETER :: modname = 'mo_orbit_solar'

! CONTAINS

#endif
END MODULE mo_orbit_solar
