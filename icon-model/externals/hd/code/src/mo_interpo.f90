! mo_interpo.f90 - Weighting factors and indices for time interpolation
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_interpo

  ! Weighting factors and indices for time interpolation
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  USE mo_kind,  ONLY: dp

  IMPLICIT NONE

  PUBLIC

  REAL(dp):: wgt1, wgt2
  INTEGER :: nmw1, nmw2, nmw1cl, nmw2cl

  REAL(dp):: wgtd1, wgtd2
  INTEGER :: ndw1, ndw2

  ! weightings and indicies for radiation calculation time step

  REAL(dp):: wgt1_m, wgt2_m
  INTEGER :: nmw1_m, nmw2_m
 
END MODULE mo_interpo
