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

! This module provides definition of types of the ICON ocean surface module,
! i.e. the sea ice and coupling between the atmopshere and the ocean model.
!
! Types are grouped as follows:
! - Sea ice:                                   t_sea_ice, t_sea_ice_budgets
! - Surface fluxes (categories):               t_atmos_fluxes

MODULE mo_sea_ice_types
  USE mo_kind,                ONLY: wp
  USE mo_math_types,          ONLY: t_cartesian_coordinates

  IMPLICIT NONE
  PRIVATE

  PUBLIC  :: t_sea_ice
  PUBLIC  :: t_sea_ice_budgets

  PUBLIC  :: t_atmos_fluxes

!! ---------------- OBSOLETE, replaced by p_oce_sfc ---------------------------------
!  PUBLIC  :: t_sfc_flx
!  PUBLIC  :: t_ptr2d

! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------

  ! -----------------
  ! Sea ice types
  ! -----------------
  !
  TYPE t_sea_ice_budgets
  ! accumulated fields of the sea-ice state
    REAL(wp), POINTER :: &
      & salt_00    (:,:)       !,   &
! --- not currently in use ---
!      & salt_01    (:,:)       ,   &
!      & salt_02    (:,:)       ,   &
!      & salt_03    (:,:)
  END TYPE t_sea_ice_budgets

  TYPE t_sea_ice

  ! The description of the sea-ice state, defined on cell-centers
  ! dimension: (nproma, k, nblks_c)

    INTEGER ::  kice                   ! Number of ice-thickness classes

    REAL(wp), POINTER :: &
      & hi         (:,:,:)       ,   & ! Ice thickness                                          [m]
      & hs         (:,:,:)       ,   & ! Snow thickness                                         [m]
      & conc       (:,:,:)       ,   & ! Ice concentration in each ice class
      & concSum    (:,:)         ,   & ! Total ice concentration within a grid cell
      & vol        (:,:,:)       ,   & ! Ice volume                                             [m^3]
      & vols       (:,:,:)       ,   & ! Snow volume                                            [m^3]
      & delhi      (:,:,:)       ,   & ! Change in mean ice thickness due to growth/melt        [m]
      & delhs      (:,:,:)       ,   & ! Change in mean snow thickness due to growth/melt       [m]
      ! functionality of hiold and hsold was replaced by delhi and delhs
      ! to be deleted together with the old thermodynamics routines
      & hiold      (:,:,:)       ,   & ! Ice thickness at the begining of ice time step         [m]
      & hsold      (:,:,:)       ,   & ! Snow thickness at the begining of ice time step        [m]
      ! thermodynamics, fast
      & Tsurf      (:,:,:)       ,   & ! Surface temperature                                    [C]
      & Tfw        (:,:)         ,   & ! Ocean surface freezing temperature                     [C]
      & Qtop       (:,:,:)       ,   & ! Energy flux available for surface melting              [W/m2]
      & Qbot       (:,:,:)       ,   & ! Conductive heat flux at ice-ocean interface (fast)     [W/m2]
      & alb        (:,:,:)       ,   & ! Albedo of snow-ice system
      ! thermodynamics, slow
      ! Note: variables like zHeatOceI and heatOceI SHOULD NOT have kice dimension
      !       to be changed once old thermodynamics routines are deleted
      !       delhi, delhs, snow_to_ice SHOULD NOT have kice dim and should represent cell-mean (mult by conc)
      !       just like totalsnowfall doesn't
      & zHeatOceI  (:,:,:)       ,   & ! Heat flux that goes into ice from below                [W/m^2]
      & heatOceI   (:,:,:)       ,   & ! Heat flux to ocean due to the ice growth               [W/m^2]
      & heatOceW   (:,:)         ,   & ! Heat flux to ocean from the atmosphere                 [W/m^2]
      & Qbot_slow  (:,:,:)       ,   & ! Engery flux at ice-ocean interface (slow)              [W/m^2]
      & snow_to_ice(:,:,:)       ,   & ! Amount of snow that is transformed to ice              [m]
      & newice     (:,:)         ,   & ! New ice growth in open water                           [m]
      & totalsnowfall(:,:)       ,   & ! Snow fall on ice-covered part of cell (water equiv.)   [m]
      & draft      (:,:,:)       ,   & ! Water equivalent of ice and snow over ice covered area [m]
      & draftave   (:,:)         ,   & ! Averaged water equivalent of ice and snow over grid area [m]
      & zUnderIce  (:,:)         ,   & ! water in upper ocean grid cell below ice               [m]
      & draftave_old(:,:)         ,   & ! water in upper ocean grid cell below ice               [m]
      ! thermodynamics, fast (winton scheme only -- not currently functional)
      & T1         (:,:,:)       ,   & ! Temperature upper layer                                [C]
      & T2         (:,:,:)       ,   & ! Temperature lower layer                                [C]
      & E1         (:,:,:)       ,   & ! Energy content upper layer                             [Jm/kg]
      & E2         (:,:,:)       ,   & ! Energy content lower layer                             [Jm/kg]
      ! thermodynamics, slow (winton scheme only -- not currently functional)
      & surfmelt   (:,:,:)       ,   & ! Surface melt water running into ocean                  [m]
      & surfmeltT  (:,:,:)       ,   & ! Mean temperature of surface melt water                 [C]
      & evapwi     (:,:,:)             ! Amount of evaporated water if no ice left              [kg/m2]

      ! dynamics
    REAL(wp), POINTER :: &
      & u_prog     (:,:)         ,   & ! Zonal velocity (prognostic, rotated grid)              [m/s]
      & v_prog     (:,:)         ,   & ! Meridional velocity (prognostic, rotated)              [m/s]
      & u          (:,:)         ,   & ! Zonal velocity on cell centre (diagnostic)             [m/s]
      & v          (:,:)         ,   & ! Meridional velocity on cell centre (diagn.)            [m/s]
      & vn_e       (:,:)         ,   &  ! Edge normal velocity(prognostic)                      [m/s]
      & vt_e       (:,:)         ,   & ! Edge tangential velocity prognostic                     [m/s]
      & Delta      (:,:)           ! Delta                                               [1/s]   
    ! not currently used categorywise limiter
    REAL(wp), ALLOCATABLE :: hi_lim(:) ! Thickness limit                                        [m]

    TYPE(t_sea_ice_budgets) :: budgets

  END TYPE t_sea_ice

  ! global type variables
  TYPE(t_sea_ice),PUBLIC, SAVE, TARGET :: v_sea_ice

! ---------------------------------------------------------------------------------------

  ! --------------------------------
  ! Surface flux types for sea ice
  ! --------------------------------
  !
  TYPE t_atmos_fluxes

    ! #slo# 2015-02: all fluxes used in sea ice model now defined as pointer, sized via add_var
    ! TODO: move to class t_sea_ice and/or reorganize/cleanup
    ! dimension: (nproma, {k,} nblks_c)

    ! at ice surface
    REAL(wp), POINTER ::  &
      & sens        (:,:,:),           & ! Sensible heat flux at ice surface           [W/m2]
      & lat         (:,:,:),           & ! Latent heat flux at ice surface             [W/m2]
      & LWnet       (:,:,:),           & ! net LW radiation flux at ice surface        [W/m2]
      & SWnet       (:,:,:),           & ! net SW radiation flux over ice              [W/m2]
      & dsensdT     (:,:,:),           & ! d sensible Flux / d T_surf                  [W/m2/K]
      & dlatdT      (:,:,:),           & ! d latent Flux / d T_surf                    [W/m2/K]
      & dLWdT       (:,:,:),           & ! d radiation Flux / d T_surf                 [W/m2/K]
      & stress_x    (:,:),             & ! Wind stress at the ice surface              [Pa]
      & stress_y    (:,:)                ! Wind stress at the ice surface              [Pa]

    ! at open ocean surface
    REAL(wp), POINTER ::   &
      & sensw       (:,:),             & ! Sensible heat flux over water               [W/m2]
      & latw        (:,:),             & ! Latent heat flux over water                 [W/m2]
      & LWnetw      (:,:),             & ! net LW radiation flux over water            [W/m2]
      & SWnetw      (:,:),             & ! net SW radiation flux over water            [W/m2]
      & stress_xw   (:,:),             & ! Wind stress at the ocean surface            [Pa]
      & stress_yw   (:,:)                ! Wind stress at the ocean surface            [Pa]

    ! precipitation over entire gridcell
    REAL(wp), POINTER ::   &
      & rprecw      (:,:),             & ! liquid precipitation rate                   [m/s]
      & rpreci      (:,:)                ! solid  precipitation rate                   [m/s]

!   Lwin and LWout is a useless nomenclature, only LWnet is calculated in the bulk formulas
!      & LWin        (:,:),             & ! incoming LW radiation flux                  [W/m2]
!      & LWout       (:,:,:),           & ! outgoing LW radiation flux at ice surface   [W/m2]
!      & LWoutw      (:,:),             & ! outgoing LW radiation flux over water       [W/m2]
!   bot was never used. Heat exchange at ice bottom is given by zHeatOceI
!      & bot         (:,:,:),           & ! Ocean heat flux at ice bottom               [W/m2]

    ! Albedos
    REAL(wp), POINTER ::     &
      & albvisdir   (:,:,:),           & ! VIS direct/paralell (ice)
      & albvisdif   (:,:,:),           & ! VIS diffuse (ice)
      & albnirdir   (:,:,:),           & ! NIR direct/paralell (ice)
      & albnirdif   (:,:,:),           & ! NIR diffuse (ice)
      & albvisdirw  (:,:),             & ! VIS direct/paralell (ocean)
      & albvisdifw  (:,:),             & ! VIS diffuse (ocean)
      & albnirdirw  (:,:),             & ! NIR direct/paralell (ocean)
      & albnirdifw  (:,:)                ! NIR diffuse (ocean)

    INTEGER ::     counter

! ------------------------------------------------------------------------------------------------------------------
! Variables below are replaced by the respective variables in t_ocean_surface.
! Only used in the old formulation of sea ice thermodynamics and update_ocean_surface
! Everything below is to be deleted once old thermodynamics routines are deleted
! ------------------------------------------------------------------------------------------------------------------
    REAL(wp), POINTER ::   &
      ! relaxaton
!      &  topBoundCond_Temp_vdiff  (:,:),     & ! forcing of temperature in vertical diffusion equation     [K*m/s]
!      &  topBoundCond_Salt_vdiff  (:,:),     & ! forcing of salinity in vertical diffusion equation        [psu*m/s]
!      &  data_surfRelax_Temp      (:,:),     & ! contains data to which temperature is relaxed             [K]
!      &  data_surfRelax_Salt      (:,:),     & ! contains data to which salinity is relaxed                [psu]
!      &  HeatFlux_Relax           (:,:),     & ! surface heat flux due to relaxation                       [W/m2]
!      &  FrshFlux_Relax           (:,:),     & ! surface freshwater flux due to relaxation                 [m/s]
!      &  TempFlux_Relax           (:,:),     & ! temperature tracer flux due to relaxation                 [K/s]
!      &  SaltFlux_Relax           (:,:),     & ! salinity tracer flux due to relaxation                    [psu/s]
      ! heat fluxes
      &  HeatFlux_ShortWave       (:,:),     & ! surface short wave heat flux                              [W/m2]
      &  HeatFlux_LongWave        (:,:),     & ! surface long wave heat flux                               [W/m2]
      &  HeatFlux_Sensible        (:,:),     & ! surface sensible heat flux                                [W/m2]
      &  HeatFlux_Latent          (:,:),     & ! surface latent heat flux                                  [W/m2]
      &  HeatFlux_Total           (:,:),     & ! sum of forcing surface heat flux                          [W/m2]
      ! fresh water + salt
      &  FrshFlux_Precipitation   (:,:),     & ! total precipitation flux                                  [m/s]
      &  FrshFlux_SnowFall        (:,:),     & ! total snow flux                                           [m/s]
      &  FrshFlux_Evaporation     (:,:),     & ! evaporation flux                                          [m/s]
      &  FrshFlux_Runoff          (:,:),     & ! river runoff flux                                         [m/s]
      &  FrshFlux_TotalSalt       (:,:),     & ! sum of forcing surface freshwater flux from BC            [m/s]
      &  FrshFlux_TotalOcean      (:,:),     & ! forcing surface freshwater flux at open ocean             [m/s]
      &  FrshFlux_TotalIce        (:,:),     & ! forcing surface freshwater flux under sea ice             [m/s]
      &  FrshFlux_VolumeIce       (:,:),     & ! forcing volume flux for height equation under sea ice     [m/s]
      &  FrshFlux_VolumeTotal     (:,:),     & ! sum of forcing volume flux including relaxation           [m/s]
!      &  cellThicknessUnderIce    (:,:),     &
      !  surface wind stress
      &  topBoundCond_windStress_u(:,:),     & ! forcing of zonal component of velocity equation,
      &  topBoundCond_windStress_v(:,:)        ! forcing of meridional component of velocity equation,

    TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
      & ALLOCATABLE :: topBoundCond_windStress_cc(:,:)

  END TYPE t_atmos_fluxes

END MODULE mo_sea_ice_types
