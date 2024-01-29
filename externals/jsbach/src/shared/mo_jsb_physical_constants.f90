!> Contains physical constants
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
MODULE mo_jsb_physical_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_physical_constants_iface

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER :: &
    & dens_soil          = 2700._wp,   & !< Density of soil particulate (Hillel 1980, Johansen 1977)           [kg/m3]
    & dens_org           = 1300._wp,   & !< Density of organic part of soil                                    [kg/m3]
    & dens_ice           = rhoi,       & !< Density of frozen water in soil (ice)                              [kg/m3]
    & dens_snow          = rhos,       & !< =300.0  Density of snow                                            [kg/m3]
    & dens_snow_min      = 50._wp,     & !< Minimum density of snow (for fresh snow)                           [kg/m3]
    & dens_snow_max      = 300._wp       !< Maximum density of snow                                            [kg/m3]

  REAL(wp), PARAMETER :: &
    molarMassDryAir = amd   * 1.e-3_wp, & ! 28.970e-3_wp. Mass [kg] per amount of substance [mol] of dry air
    molarMassCO2    = amco2 * 1.e-3_wp    ! 44.0095e-3_wp. Mass [kg] per amount of substance [mol] of CO2

  !! quincy
  ! physical constants
  REAL(wp), PARAMETER :: &
    r_gas        = 8.3144_wp        ,  &      !< universal gas constant (J / mol / K) 
    r_gas_dryair = 287.058_wp       ,  &      !< specific gas constant for dry air (J / kg / K) 
    r_gas_vapour = 461.51_wp        ,  &      !< specific gas constant for water vapour (J / kg / K) 
    Dwv          = 0.0000242_wp     ,  &      !< diffusion coefficient for water vapour in air (m2 s-1) 
    Dco2         = 0.0000151_wp               !< diffusion coefficient for CO2 in air (m2 s-1) 
    ! Tzero == tmelt -> used from ECHAM echam/src/echam/jsbach4_util_mods.f90 via mo_physical_constants_iface
  ! molar masses
  REAL(wp), PARAMETER :: &
    molar_mass_C    = 12.0112_wp    ,  &      !< g C / mol C
    molar_mass_C12  = 12.0_wp       ,  &      !< g C / mol C
    molar_mass_C13  = 13.0_wp       ,  &      !< g C / mol C
    molar_mass_C14  = 14.0_wp       ,  &      !< g C / mol C
    molar_mass_N    = 14.0072_wp    ,  &      !< g N / mol N
    molar_mass_N14  = 14.0_wp       ,  &      !< g N / mol N
    molar_mass_N15  = 15.0_wp       ,  &      !< g N / mol N
    molar_mass_P    = 30.997_wp               !< g P / mol P
  ! isotope standards
  REAL(wp), PARAMETER :: &
    PDB_STANDARD_C13           = 0.0112372_wp                    , &  !< 13C/12C PDB standard
    AAB_STANDARD_C14           = 0.2261_wp                       , &  !< 0.95·0.238 Bq/gC is the absolute specific activity of 
                                                                      !! the radiocarbon standard (i.e. 95% of the activity of the OxA-I standard)
    ATM_STANDARD_N15           = 0.0036765_wp                    , &  !< 15N/14N atmospheric standard, Robinson 2001, TREE
    AvogadroC                  = 6.022_wp * 1.e23_wp             , &  !< Avogadro Number
    lambda_C14                 = 3.8332_wp * 1.e-12_wp           , &  !< decay constant of C14 (1/s)
    AvogadroC_times_lambda_C14 = 6.022_wp * 3.8332_wp * 1.e11_wp , &  !< Avogadro Number * decay constant of C14 to 
                                                                      !! reduce the use of the very large/small AvogadroC/lambda_C14 
    scale_correction_c14       = 1.e12_wp                             !< conversion from mol to pico-mol for C14 (for numerical precision) 
#endif
END MODULE mo_jsb_physical_constants
