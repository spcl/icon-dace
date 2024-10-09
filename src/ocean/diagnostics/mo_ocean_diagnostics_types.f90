! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Contains the variables to set up the ocean model.

#include "iconfor_dsl_definitions.inc"

MODULE mo_ocean_diagnostics_types

  USE mo_kind,                ONLY: wp, sp
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    & success, max_char_length, min_dolic,               &
    & full_coriolis, beta_plane_coriolis,                &
    & f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_math_types,          ONLY: t_cartesian_coordinates,      &
    & t_geographical_coordinates

  PUBLIC :: t_ocean_monitor

  PUBLIC :: t_ocean_regions
  PUBLIC :: t_ocean_region_volumes
  PUBLIC :: t_ocean_region_areas
  PUBLIC :: t_ocean_basins


  ! diagnostic variables
  TYPE t_ocean_monitor
    REAL(wp), POINTER :: amoc26n(:)
    REAL(wp), POINTER :: volume(:)
    REAL(wp), POINTER :: kin_energy(:)
    REAL(wp), POINTER :: pot_energy(:)
    REAL(wp), POINTER :: total_energy(:)
    REAL(wp), POINTER :: global_heat_content(:)
    REAL(wp), POINTER :: global_heat_content_solid(:)
    REAL(wp), POINTER :: total_salt(:)
    REAL(wp), POINTER :: total_saltinseaice(:)
    REAL(wp), POINTER :: total_saltinliquidwater(:)
    REAL(wp), POINTER :: vorticity(:)
    REAL(wp), POINTER :: enstrophy(:)
    REAL(wp), POINTER :: ssh_global(:)
    REAL(wp), POINTER :: sst_global(:)
    REAL(wp), POINTER :: sss_global(:)
    REAL(wp), POINTER :: potential_enstrophy(:)
    REAL(wp), POINTER :: HeatFlux_ShortWave(:)
    REAL(wp), POINTER :: HeatFlux_LongWave(:)
    REAL(wp), POINTER :: HeatFlux_Sensible(:)
    REAL(wp), POINTER :: HeatFlux_Latent(:)
    REAL(wp), POINTER :: HeatFlux_Total(:)
    REAL(wp), POINTER :: FrshFlux_Precipitation(:)
    REAL(wp), POINTER :: FrshFlux_SnowFall(:)
    REAL(wp), POINTER :: FrshFlux_Evaporation(:)
    REAL(wp), POINTER :: FrshFlux_Runoff(:)
    REAL(wp), POINTER :: FrshFlux_TotalSalt(:)
    REAL(wp), POINTER :: FrshFlux_TotalOcean(:)
    REAL(wp), POINTER :: FrshFlux_TotalIce(:)
    REAL(wp), POINTER :: FrshFlux_VolumeIce(:)
    REAL(wp), POINTER :: FrshFlux_VolumeTotal(:)
    REAL(wp), POINTER :: HeatFlux_Relax(:)
    REAL(wp), POINTER :: FrshFlux_Relax(:)
    REAL(wp), POINTER :: TempFlux_Relax(:)
    REAL(wp), POINTER :: SaltFlux_Relax(:)
    REAL(wp), POINTER :: totalsnowfall(:)

    REAL(wp), POINTER :: ice_volume_nh(:)!                                                           [km3]
    REAL(wp), POINTER :: ice_volume_sh(:)!                                                           [km3]
    REAL(wp), POINTER :: ice_extent_nh(:)!                                                           [km2]
    REAL(wp), POINTER :: ice_extent_sh(:)!                                                           [km2]
    ! ice transport through {{{
    REAL(wp), POINTER :: ice_framStrait(:) !                                                          [Sv]
    ! }}}
    ! throug, POINTER  flows {{{
    REAL(wp), POINTER :: gibraltar(:)     ! though flow                                               [Sv]
    REAL(wp), POINTER :: denmark_strait(:)! though flow                                               [Sv]
    REAL(wp), POINTER :: drake_passage(:) ! though flow                                               [Sv]
    REAL(wp), POINTER :: indonesian_throughflow(:) !                                                  [Sv]
    REAL(wp), POINTER :: scotland_iceland(:) !                                                        [Sv]
    REAL(wp), POINTER :: mozambique(:)
    REAL(wp), POINTER :: framStrait(:)
    REAL(wp), POINTER :: beringStrait(:)
    REAL(wp), POINTER :: barentsOpening(:)
    REAL(wp), POINTER :: agulhas(:)
    REAL(wp), POINTER :: agulhas_long(:)
    REAL(wp), POINTER :: agulhas_longer(:)
    REAL(wp), POINTER :: florida_strait(:)
    ! }}}
    REAL(wp), ALLOCATABLE :: tracer_content(:)
  END TYPE t_ocean_monitor

    
  
  !----------------------------------------------------------------------------
  !
  ! Ocean areas/regions:
  !  0 = land point
  !  1 = Greenland-Iceland-Norwegian Sea
  !  2 = Arctic Ocean
  !  3 = Labrador Sea
  !  4 = North Atlantic Ocean
  !  5 = Tropical Atlantic Ocean
  !  6 = Southern Ocean
  !  7 = Indian Ocean
  !  8 = Tropical Pacific Ocean
  !  9 = North Pacific Ocean
  !
  !-----------------------------
  TYPE t_ocean_regions
    INTEGER :: &
      & land                            = 0,&
      & greenland_iceland_norwegian_sea = 1,&
      & arctic_ocean                    = 2,&
      & labrador_sea                    = 3,&
      & north_atlantic                  = 4,&
      & tropical_atlantic               = 5,&
      & southern_ocean                  = 6,&
      & indian_ocean                    = 7,&
      & tropical_pacific                = 8,&
      & north_pacific                   = 9,&
      & caribbean                       = -33
  END TYPE t_ocean_regions
  TYPE t_ocean_region_volumes
    REAL(wp)            :: &
      & land                            = 0.0_wp,&
      & greenland_iceland_norwegian_sea = 0.0_wp,&
      & arctic_ocean                    = 0.0_wp,&
      & labrador_sea                    = 0.0_wp,&
      & north_atlantic                  = 0.0_wp,&
      & tropical_atlantic               = 0.0_wp,&
      & southern_ocean                  = 0.0_wp,&
      & indian_ocean                    = 0.0_wp,&
      & tropical_pacific                = 0.0_wp,&
      & north_pacific                   = 0.0_wp,&
      & caribbean                       = 0.0_wp,&
      & total                           = 0.0_wp
  END TYPE t_ocean_region_volumes
  TYPE t_ocean_region_areas
    REAL(wp)            :: &
      & land                            = 0.0_wp,&
      & greenland_iceland_norwegian_sea = 0.0_wp,&
      & arctic_ocean                    = 0.0_wp,&
      & labrador_sea                    = 0.0_wp,&
      & north_atlantic                  = 0.0_wp,&
      & tropical_atlantic               = 0.0_wp,&
      & southern_ocean                  = 0.0_wp,&
      & indian_ocean                    = 0.0_wp,&
      & tropical_pacific                = 0.0_wp,&
      & north_pacific                   = 0.0_wp,&
      & caribbean                       = 0.0_wp,&
      & total                           = 0.0_wp
  END TYPE t_ocean_region_areas
  !-----------------------------
  !
  ! Ocean basins:
  !  1: Atlantic; 3: Pacific, for Indean and pacific the area values ara used
  !
  !-----------------------------
  TYPE t_ocean_basins
    INTEGER :: &
      & atlantic = 1, pacific = 3
  END TYPE t_ocean_basins
  !-----------------------------
  
END MODULE mo_ocean_diagnostics_types

