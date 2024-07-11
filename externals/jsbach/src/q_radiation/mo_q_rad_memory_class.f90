!> QUINCY radiation process memory
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
!>#### definition and init of (memory) variables for the radiation process
!>
MODULE mo_q_rad_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: BARE_TYPE, VEG_TYPE, LAND_TYPE, GLACIER_TYPE, LAKE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  dsl4jsb_Use_processes Q_RAD_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_rad_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 90

  ! ======================================================================================================= !
  !>Type definition for radiation memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_rad_memory

    TYPE(t_jsb_var_real2d)  :: &
      & sw_srf_net,            & !< net shortwave radiation at surface [W m-2]
      & swvis_srf_net,         & !< net shortwave radiation at surface in visible range [W m-2]
      & swnir_srf_net,         & !< net shortwave radiation at surface in NIR range [W m-2]
      & lw_srf_net,            & !< net longwave radiation at surface [W m-2]
      & rad_srf_net              !< net radiation at surface [W m-2]

    TYPE(t_jsb_var_real2d)  :: &
      & alb_vis,               & !< Albedo of the whole tile in the VIS (visible) range [unitless]
      & alb_nir,               & !< Albedo of the whole tile in the NIR (near-infrared) range [unitless]
      & alb_vis_snow,          & !< Albedo of snow in the VIS range [unitless]
      & alb_nir_snow,          & !< Albedo of snow in the NIR range [unitless]
      & alb_vis_pond,          & !< Albedo of ponds in the VIS range [unitless]
      & alb_nir_pond             !< Albedo of ponds in the NIR range [unitless]
    ! soil
    TYPE(t_jsb_var_real2d) :: &
      & alb_vis_soil,         & !< Albedo in VIS range of soil [unitless]
                                !<  Initialized from bc_land_phys.nc file.
                                !<  If use_albedo_soil_scheme =
                                !<    .TRUE.:  then it is calculated and may incl. litter layer and carbon within soil dependent
                                !<             on the use_alb_soil_organic_C and use_alb_soil_litter switches.
                                !<    .FALSE.: it is not calculated but remains as taken from bc_land_phys.nc file during
                                !<             initialization. alb_vis_soil of the bc_land_phys.nc file represents the soil albedo
                                !<             including carbon in the soil as well as litter on the soil, but without the canopy
                                !<             above. It was calculated from MODIS data (aboard the TERRA satellite) using linear
                                !<             regression to project the albedo when fapar=0 (LAI/Canopy =0) (for details see
                                !<             JSBACH3 documentation 4.1.1 or ask Thomas Raddatz).
      & alb_nir_soil            !< Same as alb_vis_soil but for NIR.

    ! vegetation
    TYPE(t_jsb_var_real2d) :: &
      & alb_vis_can,          & !< Albedo of canopy (without snow cover) in VIS range. [unitless]
                                !<   is calculated at runtime
      & alb_nir_can             !< Albedo of canopy (without snow cover) in NIR range. [unitless]
                                !<   is calculated at runtime

    TYPE(t_jsb_var_real2d)  :: &
      & net_radiation                   !< surface net radiation [W m-2]

    TYPE(t_jsb_var_real2d)  :: &
      & albedo,                &        !< total surface albedo [unitless]
      & albedo_noon                     !< albedo at local solar noon, value is modified only every 24h [unitless]

    TYPE(t_jsb_var_real2d)  :: &
      & arad_vis_soil,                & !< radiation in VIS range absorbed by the soil [W m-2]
      & arad_nir_soil                   !< radiation in NIR range absorbed by the soil [W m-2]

    TYPE(t_jsb_var_real2d)  :: &
      & arad_vis_can,                 & !< radiation in VIS range absorbed by the canopy [W m-2]
      & arad_nir_can,                 & !< radiation in NIR range absorbed by the canopy [W m-2]
      & appfd,                        & !< absorbed Photosynthetically Active Photon Flux Density (a.k.a. faPAR) [micro-mol m-2 s-1]
      & rfr_ratio_boc,                & !< red to far-red ratio at the bottom of the canopy [unitless]
      & rfr_ratio_boc_tvegdyn_mavg      !< red to far-red ratio at the bottom of the canopy [unitless]

    TYPE(t_jsb_var_real3d)  :: &   ! DIMENSION(ncanopy)
      & ppfd_sunlit_cl                , & !< Photosynthetically Active Photon Flux Density on sunlit leaves [micro-mol m-2 s-1]
      & ppfd_shaded_cl                , & !< Photosynthetically Active Photon Flux Density on shaded leaves [micro-mol m-2 s-1]
      & ppfd_sunlit_daytime_cl        , & !< average PPFD on sunlit leaves of the previous day [micro-mol m-2 s-1]
      & ppfd_shaded_daytime_cl        , & !< average PPFD on shaded leaves of the previous day [micro-mol m-2 s-1]
      & ppfd_sunlit_daytime_dacc_cl   , & !< daily accumulated PPFD on sunlit leaves [micro-mol m-2 d-1]
      & ppfd_shaded_daytime_dacc_cl   , & !< daily accumulated PPFD on shaded leaves of the previous day [micro-mol m-2 d-1]
      & ppfd_sunlit_tfrac_mavg_cl     , & !< average PPFD on sunlit leaves at tfrac-timescale [micro-mol m-2 s-1]
      & ppfd_shaded_tfrac_mavg_cl     , & !< average PPFD on shaded leaves at tfrac-timescale [micro-mol m-2 s-1]
      & ppfd_sunlit_tcnl_mavg_cl      , & !< average PPFD on sunlit leaves at tcnl-timescale [micro-mol m-2 s-1]
      & ppfd_shaded_tcnl_mavg_cl          !< average PPFD on shaded leaves at tcnl-timescale [micro-mol m-2 s-1]

    TYPE(t_jsb_var_real3d)  :: &   ! DIMENSION(ncanopy)
      & arad_sunlit_vis_cl          , & !< absorbed radiation in visible band on sunlit leaves [W m-2]
      & arad_shaded_vis_cl          , & !< absorbed radiation in visible band on shaded leaves [W m-2]
      & fleaf_sunlit_vis_cl         , & !< fraction of leaves that are sunlit in VIS band [unitless]
      & arad_sunlit_nir_cl          , & !< absorbed radiation in NIR band on sunlit leaves [W m-2]
      & arad_shaded_nir_cl          , & !< absorbed radiation in NIR band on shaded leaves [W m-2]
      & fleaf_sunlit_nir_cl             !< fraction of leaves that are sunlit in NIR band [unitless]

    ! 3D Radiative transfer variables (appended by _sp for spectral_bands)
    TYPE(t_jsb_var_real3d)  :: &   ! DIMENSION(nspec)
      & canopy_reflectance_sp           !< top of the canopy directional reflectance factor [unitless]

  CONTAINS
    PROCEDURE :: Init => Init_q_rad_memory
  END TYPE t_q_rad_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_rad_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: radiation
  !>
  SUBROUTINE Init_q_rad_memory(mem, prefix, suffix, lct_ids, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_rad_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),   POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface                      ! Vertical grid
    TYPE(t_jsb_vgrid),  POINTER :: vgrid_canopy_q_assimi        ! Vertical grid
    TYPE(t_jsb_vgrid),  POINTER :: vgrid_spectral               ! Vertical grid used to display the sprectal bands of light
    INTEGER                     :: table                        ! ...
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_rad_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    table                 = tables(1)
    hgrid                 => Get_grid(mem%grid_id)
    surface               => Get_vgrid('surface')
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
#ifdef __QUINCY_STANDALONE__
    vgrid_spectral => Get_vgrid('spectral_bands')
#endif

    CALL mem%Add_var('sw_srf_net', mem%sw_srf_net, &
      & hgrid, surface, &
      & t_cf('sw_srf_net', 'W m-2', 'net shortwave radiation at surface'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('swvis_srf_net', mem%swvis_srf_net, &
      & hgrid, surface, &
      & t_cf('swvis_srf_net', 'W m-2', 'net shortwave radiation at surface in visible range'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('swnir_srf_net', mem%swnir_srf_net, &
      & hgrid, surface, &
      & t_cf('swnir_srf_net', 'W m-2', 'net shortwave radiation at surface in NIR range'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('lw_srf_net', mem%lw_srf_net, &
      & hgrid, surface, &
      & t_cf('lw_srf_net', 'W m-2', 'net longwave radiation at surface'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('rad_srf_net', mem%rad_srf_net, &
      & hgrid, surface, &
      & t_cf('rad_srf_net', 'W m-2', 'net radiation at surface'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_vis', mem%alb_vis, &
      & hgrid, surface, &
      & t_cf('alb_vis', 'unitless', 'Albedo of the whole tile in the VIS (visible) range'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_nir', mem%alb_nir, &
      & hgrid, surface, &
      & t_cf('alb_nir', 'unitless', 'Albedo of the whole tile in the NIR (near-infrared) range'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_vis_snow', mem%alb_vis_snow, &
      & hgrid, surface, &
      & t_cf('alb_vis_snow', 'unitless', 'Albedo of snow in the VIS range'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_nir_snow', mem%alb_nir_snow, &
      & hgrid, surface, &
      & t_cf('alb_nir_snow', 'unitless', 'Albedo of snow in the NIR range'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    ! Additional variables for soil
    IF ( (One_of(BARE_TYPE, lct_ids(:)) > 0 .OR. &
      &   One_of(VEG_TYPE, lct_ids(:)) > 0 .OR. &
      &   One_of(LAND_TYPE, lct_ids(:)) > 0)                                 &
      & ) THEN

      CALL mem%Add_var('alb_vis_soil', mem%alb_vis_soil, &
        & hgrid, surface, &
        & t_cf('alb_vis_soil', 'unitless', 'Albedo in VIS range of soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('alb_nir_soil', mem%alb_nir_soil, &
        & hgrid, surface, &
        & t_cf('alb_nir_soil', 'unitless', 'Albedo in VIS range of soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
    END IF ! soil

    ! Additional variables for vegetation
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      CALL mem%Add_var('alb_vis_can', mem%alb_vis_can, &
        & hgrid, surface, &
        & t_cf('alb_vis_can', 'unitless', 'Albedo of canopy (without snow cover) in VIS range'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('alb_nir_can', mem%alb_nir_can, &
        & hgrid, surface, &
        & t_cf('alb_nir_can', 'unitless', 'Albedo of canopy (without snow cover) in NIR range'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_radiation', mem%net_radiation, &
        & hgrid, surface, &
        & t_cf('net_radiation', 'W m-2', 'surface net radiation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp )

      CALL mem%Add_var('albedo', mem%albedo, &
        & hgrid, surface, &
        & t_cf('albedo', 'unitless', 'total surface albedo'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

        CALL mem%Add_var('albedo_noon', mem%albedo_noon, &
        & hgrid, surface, &
        & t_cf('albedo_noon', 'unitless', 'albedo at local solar noon'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_vis_soil', mem%arad_vis_soil, &
        & hgrid, surface, &
        & t_cf('arad_vis_soil', 'W m-2', 'radiation in VIS range absorbed by the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_nir_soil', mem%arad_nir_soil, &
        & hgrid, surface, &
        & t_cf('arad_nir_soil', 'W m-2', 'radiation in NIR range absorbed by the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_vis_can', mem%arad_vis_can, &
        & hgrid, surface, &
        & t_cf('arad_vis_can', 'W m-2', 'radiation in VIS range absorbed by the canopy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_nir_can', mem%arad_nir_can, &
        & hgrid, surface, &
        & t_cf('arad_nir_can', 'W m-2', 'radiation in NIR range absorbed by the canopy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('appfd', mem%appfd, &
        & hgrid, surface, &
        & t_cf('appfd', 'micro-mol m-2 s-1', 'absorbed Photosynthetically Active Photon Flux Density (faPAR)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rfr_ratio_boc', mem%rfr_ratio_boc, &
        & hgrid, surface, &
        & t_cf('rfr_ratio_boc', 'unitless', 'red to far-red ratio at the bottom of the canopy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rfr_ratio_boc_tvegdyn_mavg', mem%rfr_ratio_boc_tvegdyn_mavg, &
        & hgrid, surface, &
        & t_cf('rfr_ratio_boc_tvegdyn_mavg', 'unitless', 'red to far-red ratio at the bottom of the canopy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_sunlit_cl', mem%ppfd_sunlit_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_sunlit_cl', 'micro-mol m-2 s-1', 'Photosynthetically Active Photon Flux Density on sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_shaded_cl', mem%ppfd_shaded_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_shaded_cl', 'micro-mol m-2 s-1', 'Photosynthetically Active Photon Flux Density on shaded leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_sunlit_daytime_cl', mem%ppfd_sunlit_daytime_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_sunlit_daytime_cl', 'micro-mol m-2 s-1', 'average PPFD on sunlit leaves of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_shaded_daytime_cl', mem%ppfd_shaded_daytime_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_shaded_daytime_cl', 'micro-mol m-2 s-1', 'average PPFD on shaded leaves of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_sunlit_daytime_dacc_cl', mem%ppfd_sunlit_daytime_dacc_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_sunlit_daytime_dacc_cl', 'micro-mol m-2 d-1', 'daily accumulated PPFD on sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_shaded_daytime_dacc_cl', mem%ppfd_shaded_daytime_dacc_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_shaded_daytime_dacc_cl', 'micro-mol m-2 d-1', &
        &      'daily accumulated PPFD on shaded leaves of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_sunlit_tfrac_mavg_cl', mem%ppfd_sunlit_tfrac_mavg_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_sunlit_tfrac_mavg_cl', 'micro-mol m-2 s-1', 'average PPFD on sunlit leaves at tfrac-timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_shaded_tfrac_mavg_cl', mem%ppfd_shaded_tfrac_mavg_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_shaded_tfrac_mavg_cl', 'micro-mol m-2 s-1', 'average PPFD on shaded leaves at tfrac-timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_sunlit_tcnl_mavg_cl', mem%ppfd_sunlit_tcnl_mavg_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_sunlit_tcnl_mavg_cl', 'micro-mol m-2 s-1', 'average PPFD on sunlit leaves at tcnl-timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ppfd_shaded_tcnl_mavg_cl', mem%ppfd_shaded_tcnl_mavg_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('ppfd_shaded_tcnl_mavg_cl', 'micro-mol m-2 s-1', 'average PPFD on shaded leaves at tcnl-timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_sunlit_vis_cl', mem%arad_sunlit_vis_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('arad_sunlit_vis_cl', 'W m-2', 'absorbed radiation in visible band on sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_shaded_vis_cl', mem%arad_shaded_vis_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('arad_shaded_vis_cl', 'W m-2', 'absorbed radiation in visible band on shaded leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_vis_cl', mem%fleaf_sunlit_vis_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fleaf_sunlit_vis_cl', 'unitless', 'fraction of leaves that are sunlit in VIS band'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_sunlit_nir_cl', mem%arad_sunlit_nir_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('arad_sunlit_nir_cl', 'W m-2', 'absorbed radiation in NIR band on sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('arad_shaded_nir_cl', mem%arad_shaded_nir_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('arad_shaded_nir_cl', 'W m-2', 'absorbed radiation in NIR band on shaded leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_nir_cl', mem%fleaf_sunlit_nir_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fleaf_sunlit_nir_cl', 'unitless', 'fraction of leaves that are sunlit in NIR band'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

#ifdef __QUINCY_STANDALONE__
      CALL mem%Add_var('canopy_reflectance_sp', mem%canopy_reflectance_sp, &
        & hgrid, vgrid_spectral, &
        & t_cf('canopy_reflectance_sp', 'unitless', 'top of the canopy directional reflectance factor'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
#endif

    END IF
  END SUBROUTINE Init_q_rad_memory

#endif
END MODULE mo_q_rad_memory_class
