!> QUINCY soil-physics constants
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and define soil-physics-quincy constants
!>
MODULE mo_spq_constants
#ifndef __NO_QUINCY__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC


  REAL(wp), SAVE :: &
    soil_heat_cap,  &           !< soil heat capacity (temporarily for SLM offline) J/kg/K
    !soil_frozen_heat_cap, &     !< soil thermal conductivity frozen soil (W m-1 K-1)
    !rho_bulk_soil,  &          !< soil bulk denisty (kg / m3)
    soil_therm_cond, &          !< soil thermal conductivity (W m-1 K-1)
    soil_frozen_therm_cond, &   !< soil thermal conductivity frozen soil(W m-1 K-1)
! constants in simple multi-layer soil moisture scheme to describe the
    InterceptionEfficiency, &   !< Efficiency of interception of precipitation as rain
    w_skin_max, &               !< maximum water storage per unit LAI (m)
    drain_frac, &               !< preferential flow fraction of infiltrating water per meter soil
    fdrain_srf_runoff, &        !< fraction of surface runoff that goes to drainage (preferential flow)
    percolation_up_max, &       !< maximum upwards   directed material transport across soil layers, i.e., 'percolation_sl'
    percolation_down_max, &     !< maximum downwards directed material transport across soil layers, i.e., 'percolation_sl'
    frac_w_lat_loss_max, &      !< maximum of 'frac_w_lat_loss_sl', i.e., fraction of lateral (horizontal) water loss
                                !!  of 'w_soil_sl_old' (prev. timestep) to limit nutrient loss at extreme rainfall events
! constants for water freezing and ice melting
    latent_heat_fusion, &       !< latent heat release/uptake for water fusion freezing/thawing [J kg-1]
    w_density, &                !< Density of water   [kg m-3]
    ice_density, &              !< Density of ice     [kg m-3]
    water2ice_density_ratio, &  !< water to ice density ratio                 [-]
    fact_water_supercooled      !< factor to modify the amount of liquid water at soil temperatures below Tzero (supercooled water)  [-]

 REAL(wp), SAVE :: &
    height_wind                 !<

 ! constants in pedo-transfer functions according to Saxton & Rawls 2006
 REAL(wp) :: &
    wpot_pwp, &                 !< water potential at permanent wilitng point (Pa)
    wpot_fc, &                  !< water potential at field capacity (Pa)
    wpot_min, &                 !< minimum water potential (MPa) for numerical reasons
    k_pwp_a, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_c, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_s, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_sc, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_at, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_bt, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    k_fc_a, &                   !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_c, &                   !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_s, &                   !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_sc, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_at, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_bt, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_ct, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_sat_a, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_c, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_s, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_sc, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_at, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_bt, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_ct, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_dt, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    kdiff_sat_max               !< maximum diffusivity of soil at saturation (m/s)

  ! snow
  REAL(wp) :: &
    snow_dens,      &           !< Snow density                                                 [kg/m3]
    snow_height_min, &          !< Minimum snow height for thermal conductivity calculations
                                ! below this amount, no snow in invididual layer is assumed     [m]
    snow_heat_capa,  &          !< Volumetric heat capacity of snow                             [kJ/m3/K]
    snow_therm_cond,     &      !< Heat conductivity snow                                       [W/m/K]
    albedo_snow,         &      !< Snow albedo                                                  [-]
    w_snow_min,          &      !> minimum snow water content                                   [m]
    w_snow_max                  !> maximum water content of snow layer, should be function of snow grid depth attribute [m]
  ! constants used as factor for multiplication with namelist parameters (spq_config) when "parameter sensitivity" is enabled
  ! see spq_config & mo_qs_set_parameters:Modify_quincy_namelist_parameter
  REAL(wp) :: &
    f_psensi_soil_sand, &
    f_psensi_soil_silt, &
    f_psensi_soil_awc_prescribe, &
    f_psensi_bulk_density

  ! surface energy balance
  REAL(wp) :: &
    t_soil_max_increase, &      !< max increase in soil temperature per layer and timestep [K]
    t_soil_max_decrease         !< max decrease in soil temperature per layer and timestep; a negative number [K]

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_spq_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from spq
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_spq_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    !
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_spq_constants'


    soil_heat_cap           = 1480.0_wp     ! textbook for wet sandy soil                      [J / K]
    !soil_frozen_heat_cap    = 2000.0_wp     ! Heat capacity frozen soil, textbook knowledge    [J / K]
    !rho_bulk_soil           = 1695.86_wp    ! Bernhard Ahrens for Hainich -> using site specific value now
    soil_therm_cond         = 2.0_wp        ! Text book for wet sandy soil                     [J / K / s]
    soil_frozen_therm_cond  = 1.4_wp        ! Approximation for frozen Siberian soil, should be generally applicable to frozen soil ,Xu et al. (2020)   [J / K / s]

    InterceptionEfficiency  = 0.25_wp       ! JSBACH4
    w_skin_max              = 0.0002_wp     ! JSBACH4
    drain_frac              = 0.01_wp       ! tuned
    fdrain_srf_runoff       = 0.95_wp       ! set from ORCHIDEE
    percolation_up_max      = -0.75_wp      !< pure assumption
    percolation_down_max    = 0.75_wp       !< pure assumption
    frac_w_lat_loss_max     = 0.75_wp       !< simply following percolation_down_max
    height_wind             = 10.0_wp       ! JSBACH ?
    latent_heat_fusion      = 333550.0_wp   ! chemistry textbook                 [J / kg]
    w_density               = 998.0_wp      ! water density at 20 degrees C (e.g. Tanaka et al., 2001)                 [kg / m3]
    ice_density             = 916.7_wp      ! density of ice, textbook value                                           [kg / m3]
    water2ice_density_ratio = w_density / ice_density !                                                  [-]
    fact_water_supercooled  = 1.0_wp        ! fraction of soil water relative to the PWP that remains liquid below 0degC         [-]
                                            ! set to 1.2 for now to be safe
    ! Saxton & Rawls, 2006, includes conversion from fractional to mass% texture
    wpot_pwp        = 1500._wp
    wpot_fc         = 33.0_wp
    wpot_min        = -15.0_wp  ! no plant can survive 8 MPa...
    k_pwp_a         = 0.031_wp
    k_pwp_c         = 0.487_wp
    k_pwp_s         = -0.024_wp
    k_pwp_sc        = 0.068_wp
    k_pwp_at        = -0.02_wp
    k_pwp_bt        = 0.14_wp
    k_fc_a          = 0.299_wp
    k_fc_c          = 0.195_wp
    k_fc_s          = -0.251_wp
    k_fc_sc         = 0.452_wp
    k_fc_at         = -0.015_wp
    k_fc_bt         = -0.373_wp
    k_fc_ct         = 1.293_wp
    k_sat_a         = 0.078_wp
    k_sat_c         = 0.034_wp
    k_sat_s         = 0.278_wp
    k_sat_sc        = -0.584_wp
    k_sat_at        = -0.107_wp
    k_sat_bt        = 0.6360_wp
    k_sat_ct        = 0.097_wp
    k_sat_dt        = 0.043_wp
    kdiff_sat_max   = 1930.0_wp / 1000.0_wp / 3600._wp ! m/s

    ! snow
    snow_dens           = 150.0_wp                                !> Literature research, 30-800 kg/m-3, adapted based on calibration with FI-So2 & Siberian sites
                                                                  !! TODO: could vary with age (see JSBACH4) or temperature
    snow_height_min     = 0.02_wp                                 !< same as jsbach4
    w_snow_min          = snow_height_min * snow_dens / w_density !< minimum snow amount for thermal conductivity calculations
    w_snow_max          = 0.05_wp * snow_dens / w_density         !< determined by snow grid, can the 0.05_wp be replaced by the grid information? Jan?
    snow_heat_capa      = 3100.0_wp                               !< constant from bried literature review, otherwise check jsbach4 range, adapted based on calibration with FI-So2 & Siberian sites
    snow_therm_cond     = 0.06_wp                                  !< constant from brief literature search, range 0.03-0.8 d, adapted based on calibration with FI-So2 & Siberian sites
                                                                  !!  depending on snow properties, within range used in jsbach4


    ! constants used as factor for multiplication with namelist parameters (spq_config) when "parameter sensitivity" is enabled
    f_psensi_soil_sand            = 1.0_wp
    f_psensi_soil_silt            = 1.0_wp
    f_psensi_soil_awc_prescribe   = 1.0_wp
    f_psensi_bulk_density         = 1.0_wp

    ! surface energy balance
    t_soil_max_increase           = 5.0_wp        !< best guess (https://gitlab.dkrz.de/jsbach/jsbach/-/issues/190#note_204075)
    t_soil_max_decrease           = -5.0_wp       !< best guess (https://gitlab.dkrz.de/jsbach/jsbach/-/issues/190#note_204075)

  END SUBROUTINE init_spq_constants

#endif
END MODULE mo_spq_constants
