!> Contains constants for the soil and snow energy process
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
MODULE mo_sse_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_physical_constants,  ONLY: ks, ki, clw, cs, ci, rhoh2o, dens_snow, dens_ice

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                    &
    ! Parameters for solid soil particulate
    & hcond_dry_soil     = 0.05_wp,         & !< Dry soil heat conductivity (O'Donnel et al 2009)             [W/m/K]
    ! Parameters for liquid water
    & vol_hcap_water     = clw * rhoh2o,    & !< Volumetric heat capacity of liquid water                     [J/m3/K]
    & hcond_water        = 0.57_wp,         & !< Thermal conductivity of liquid water (Hillel 1982)           [W/K/m]
    ! Parameters for frozen water in soil (ice)
    ! @todo: vol_hcap_ice is 2090000 in jsbach3
    & vol_hcap_ice       = ci * dens_ice,   & !< =1931202 Volumetric heat capacity of ice                     [J/m3/K]
    & hcond_ice          = ki,              & !< =2.1656  thermal conductivity of ice                         [W/K/m]
    ! Parameters for organic soil component
    ! @todo these values need to be double-checked and potentially adapted for top and deep soil layers
    & vol_hcap_org_top       = 2.5e+06_wp,  & !< Organic layer volumetric heat capacity (Beringer et al 2001) [J/m3/K]
    & hcond_org_top          = 0.25_wp,     & !< Organic layer heat conductivity (Beringer et al 2001)        [W/m/K]
    & hcond_dry_org_top      = 0.05_wp,     & !< Dry soil heat conductivity (Farouki 1981)                    [W/m/K]
    & vol_hcap_org_below     = 2.5e+06_wp,  & !< Organic layer volumetric heat capacity (Beringer et al 2001) [J/m3/K]
    & hcond_org_below        = 0.25_wp,     & !< Organic layer heat conductivity (Beringer et al 2001)        [W/m/K]
    & hcond_dry_org_below    = 0.05_wp,     & !< Dry soil heat conductivity (Farouki 1981)                    [W/m/K]
    ! Parameters for bedrock
    & vol_hcap_bedrock   = 2.e+06_wp,       & !< Volumetric heat capacity of bedrock (Bonan 2002)             [J/m3/K]
    & hcond_bedrock      = 2._wp,           & !< Thermal conductivity of bedrock (Tarnocai 2004)              [W/m/K]
    ! Parameters for free soil pores (air)
    & vol_hcap_air       = 1300._wp,        & !< Volumetric heat capacity of free soil pores (from VIC model) [J/m3/K]
    ! Parameters for snow
    & snow_depth_min     = 0.02_wp,         & !< Minimum snow depth below which no snow layer is considered   [m]
    ! @todo: vol_hcap_snow is 634500 in jsbach3
    & vol_hcap_snow      = cs * dens_snow,  & !< =627000 Volumetric heat capacity of snow                     [J/m3/K]
    & vol_hcap_snow_min  = 526500._wp,      & !< Minimum snow heat capacity (fresh snow)                      [J/m3/K]
    & hcond_snow         = ks,              & !< =0.31 Thermal conductivity of snow                           [W/m/K]
    ! @todo: hcond_snow_min is 0.1 in jsbach3
    & hcond_snow_min     = 0.03_wp            !< Minimum thermal conductivity fo snow (fresh snow)            [W/m/K]

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sse_constants'

#endif
END MODULE mo_sse_constants
