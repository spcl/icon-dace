!> Contains constants for the carbon processes
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
MODULE mo_carbon_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_jsb_varlist,            ONLY: VARNAME_LEN
  USE mo_jsb_physical_constants, ONLY: amco2  !< Molar weight of CO2 [g/mol]

  IMPLICIT NONE
  PUBLIC

  ! Yasso requires tile dependant lctlib constants for carbon relocation, accessed via indices
  INTEGER, PARAMETER  :: i_lctlib_acid           = 1 !< Index of constants for acid yasso pools in the lctlib
  INTEGER, PARAMETER  :: i_lctlib_water          = 2 !< Index of constants for water yasso pools in the lctlib
  INTEGER, PARAMETER  :: i_lctlib_ethanol        = 3 !< Index of constants for ethanol yasso pools in the lctlib
  INTEGER, PARAMETER  :: i_lctlib_nonsoluble     = 4 !< Index of constants for nonsoluble yasso pools in the lctlib
  INTEGER, PARAMETER  :: i_lctlib_humus          = 5 !< Index of constants for humus yasso pools in the lctlib

  ! Further constants
  REAL(wp), PARAMETER :: molarMassC_kg           = 12.01_wp * 1.e-3_wp  !< 12.01e-3_wp   -  Mass of 1 mol C   in kg
  REAL(wp), PARAMETER :: molarMassCO2_kg         = amco2 * 1.e-3_wp     !< 44.0095e-3_wp -  Mass of 1 mol CO2 in kg

  REAL(wp), PARAMETER :: FireFracWood2Atmos      =  0.2_wp              !< Fraction of wood carbon emitted to the atm by fire.
                                                                        !< In JS3 called def_fire_fract_wood_2_atmos.
  REAL(wp), PARAMETER :: fract_wood_aboveGround  = 0.7_wp               !< Fraction of C above ground in wood pool (for separation
                                                                        !  of woody litter into above and below ground litter pools)
  REAL(wp), PARAMETER :: fract_green_aboveGround = 0.5_wp               !< Fraction of C above ground in green pool (for separation
                                                                        !  of green litter into above and below ground litter pools)
  REAL(wp), PARAMETER :: sec_per_day             = 86400._wp            !< seconds per day
  REAL(wp), PARAMETER :: days_per_year           = 365.25_wp
  REAL(wp), PARAMETER :: sec_per_year            = days_per_year * sec_per_day

  ! Definitions for lcc

  !> list of potential active variables for the carbon process: For these the carbon interface provides active relocation
  CHARACTER(len=VARNAME_LEN) :: carbon_potential_active_vars(11)  = [character(len=VARNAME_LEN) :: &
     & 'c_green', 'c_woods', 'c_reserve',                                                          &
     & 'c_acid_ag1', 'c_water_ag1', 'c_ethanol_ag1', 'c_nonsoluble_ag1',                           &
     & 'c_acid_ag2', 'c_water_ag2', 'c_ethanol_ag2', 'c_nonsoluble_ag2' ]

  !> list of required passive vars if vars are active (because these are used to relocate the active variables)
  INTEGER, PARAMETER :: nr_of_yasso_pools = 18

  CHARACTER(len=VARNAME_LEN) :: carbon_required_passive_vars(nr_of_yasso_pools) = [character(len=VARNAME_LEN) :: &
    & 'c_acid_ag1', 'c_water_ag1', 'c_ethanol_ag1', 'c_nonsoluble_ag1',               &
    & 'c_acid_bg1', 'c_water_bg1', 'c_ethanol_bg1', 'c_nonsoluble_bg1', 'c_humus_1',  &
    & 'c_acid_ag2', 'c_water_ag2', 'c_ethanol_ag2', 'c_nonsoluble_ag2',               &
    & 'c_acid_bg2', 'c_water_bg2', 'c_ethanol_bg2', 'c_nonsoluble_bg2', 'c_humus_2' ]

  !> helper for automatization (order needs to match var name list in carbon_required_passive_vars)
  INTEGER :: coefficient_ind_of_passive_vars(nr_of_yasso_pools) = (/ &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble,                 &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble, i_lctlib_humus, &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble,                 &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble, i_lctlib_humus  /)

#endif
END MODULE mo_carbon_constants
