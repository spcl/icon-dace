!> QUINCY soil-physics process memory
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
!>#### definition and init of (memory) variables for the soil-physics-quincy process
!>
MODULE mo_spq_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_spq_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars      = 80


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for spq memory
  !!
  !! @par includes: \n
  !!    spq variables spq\%VAR
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_spq_memory


    TYPE(t_jsb_var_real2d)          :: &
                                    soil_depth      , & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
                                    root_depth      , & !< actual rooting depth [m]
                                    t_srf_old       , & !< previous surface temperature [K]
                                    zril_old        , & !< previous timestep's Reynold's number
                                    elevation           !< currently 0.0 [m]

    TYPE(t_jsb_var_real2d)          :: &
      &   num_sl_above_bedrock                          !< number of soil layers above bedrock, i.e., with a width larger eps8

    ! JSBACH SOIL
    TYPE(t_jsb_var_real3d)          :: &
                                    soil_depth_sl, &            !< Width ! of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_width_sl, &        !< width of each soil layer that can be saturated with water, identical with soil_depth_sl (until bedrock) [m]
                                    soil_lay_depth_center_sl, & !< depth at the center of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_ubound_sl, & !< depth at the upper bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_lbound_sl    !< depth at the lower bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]

    ! JSBACH HYDRO MEMORY (temporary here to get rid of vegsoil construct)
    TYPE(t_jsb_var_real2d)          :: &
                                    interception,     & !< surface interception evaporation                                         [kg m-2 s-1]
                                    pet,              & !< potential evaporation                                                    [kg m-2 s-1]
                                    evaporation,      & !< surface evaporation                                                      [kg m-2 s-1]
                                    evaporation_snow, & !< snow surface evaporation                                                 [kg m-2 s-1] unit correct ?
                                    srf_runoff,       & !< surface runoff                                                           [kg m-2 s-1]
                                    drainage,         & !< drainage                                                                 [kg m-2 s-1]
                                    drainage_fraction   !< water mass fraction lost from the soil column by drainage                [1 s-1]

    TYPE(t_jsb_var_real2d)          :: &
      & fact_q_air                      , &  !< factor representing air at actual water content [unitless]
      & fact_qsat_srf                        !< factor representing air at saturating water content [unitless]

    ! Additional variables for spq
    TYPE(t_jsb_var_real2d)          :: &
                                    w_skin,           & !< Water content in skin reservoir                                          [m]
                                    w_soil_root,      & !< Water content in root zone of the soil                                   [m]
                                    w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
                                    w_soil_root_pwp,  & !< Water content at permanent wilting point in root zone of the soil        [m]
                                    t_soil_root         !< Temperature of fine roots in the soil                                    [K]

    TYPE(t_jsb_var_real3d)          :: &
      & frac_w_lat_loss_sl, &   !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep) [unitless]
                                    gw_runoff_sl,     & !< lateral (horizontal) soilwater to groundwater runoff per layer           [kg m-2 yr-1]
                                    percolation_sl,   & !< fraction of material transferred to below layer                          [fraction s-1]
                                    drainage_sl,      & !< vertical water flow (drainage) per soil layer across soil layers         [kg m-2 s-1]
                                    w_soil_sl,        & !< Water content in soil layers                                             [m]
                                    w_soil_sat_sl,    & !< Water content of soil layers at saturation
                                                        !! (derived from volumetric capacity at saturation)                         [m]
                                    w_soil_fc_sl,     & !< Water content of soil layers at field capacity
                                                        !! (derived from volumetric field capacity)                                 [m]
                                    w_soil_pwp_sl,    & !< Water content of soil layers at permanent wilting point
                                                        !! (derived from vol. perm wilt point)                                      [m]
                                    saxtonA,          & !< coefficient in moisture-tension relationship
                                    saxtonB,          & !< coefficient in moisture-tension relationship
                                    saxtonC,          & !< exponent of moisture-conductivity relationship
                                    kdiff_sat_sl,     & !< saturated hydraulic conductivity                                         [m s-1]
                                    w_ice_sl            !< ice content of soil layer                                                [m]

    ! heat fluxes
    TYPE(t_jsb_var_real2d)          :: &
                                    sensible_heat_flx     , &  !< sensible heat flux [W m-2]
                                    latent_heat_flx       , &  !< latent heat flux [W m-2]
                                    ground_heat_flx       , &  !< ground heat flux [W m-2]
                                    ground_heat_flx_old        !< ground heat flux of previous timestep [W m-2]

    ! transpiration
    TYPE(t_jsb_var_real2d)          :: &
                                    transpiration       !< transpiration in [kg m-2 s-1]

    !Soil Water addons for SLM in HYDRO
    TYPE(t_jsb_var_real2d)          :: &
                                    w_soil_root_theta, & !< root zone integrated plant available water, fraction of maximum
                                    w_soil_root_pot      !< soil water potential [MPa]

    TYPE(t_jsb_var_real3d)          :: &
                                    w_soil_pot_sl        !< soil water potential per layer [MPa]

    ! soil energy
    TYPE(t_jsb_var_real3d)          :: &
                                    t_soil_sl            !< soil temperature [K]

    TYPE(t_jsb_var_real3d)          :: &
      & w_soil_freeze_flux, &                            !< soil water flux from liquid water to ice       [m]
      & w_soil_melt_flux                                 !< soil water flux from ice to liquid water       [m]

    ! soil physical properties (not sure where they should belong to, but they are needed in the long-term)
    TYPE(t_jsb_var_real3d)          :: &
                                    heat_capa_sl,         &  !< soil heat capacity [J kg-1]
                                    bulk_dens_sl,         &  !< soil bulk density [kg m-3]
                                    therm_cond_sl            !< soil thermal conductivity [W m-2 K-1]
    TYPE(t_jsb_var_real3d)          :: &
                                    sand_sl             , &  !< soil sand proportion of mineral soil
                                    silt_sl             , &  !< soil silt proportion of mineral soil
                                    clay_sl             , &  !< soil clay proportion of mineral soil(re-calculated as: clay = 1.0 -sand -silt)
                                    volume_min_sl       , &  !< soil mineral volumetric fraction
                                    qmax_org_min_sl     , &  !< maximum OM sorption capacity of mineral soil, [mol C kg-1 mineral soil]
                                    qmax_po4_min_sl     , &  !< maximum po4 sorption capacity of mineral soil, [mol P kg-1 mineral soil]
                                    qmax_nh4_min_sl     , &  !< maximum nh4 sorption capacity of mineral soil, [mol N kg-1 mineral soil]
                                    qmax_po4_om_sl           !< maximum po4 sorption capacity of SOM, [mol P kg-1 SOM]

    ! variables for snow
    TYPE(t_jsb_var_real2d)          :: snow_height, &           !< height (thickness) of all snow layers                                 [m]
                                       snow_soil_heat_flux, &   !< heat flux between snow and soil                          [J m-2]
                                       snow_srf_heat_flux,  &   !< ..
                                       snow_melt_to_soil        !< Snow melt water to soil flux                              [kg m-2 s-1]

    TYPE(t_jsb_var_real3d)          :: snow_present_snl       , &  !< variable to check if there is any snow for each layer     [-]
                                       t_snow_snl             , &  !< temperature of snow layer                                 [K]
                                       w_snow_snl             , &  !< water content of snow layer                               [m]
                                       snow_lay_thickness_snl , &  !< snow layer thickness                                      [m]
                                       w_snow_max_snl         , &  !< maximum water content of snow layer                       [m]
                                       w_liquid_snl                !< Liquid water in snow layers                               [m]


  CONTAINS
    PROCEDURE :: Init => Init_spq_memory
  END TYPE t_spq_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SPQ_ process
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_spq_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC , MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_spq_memory), INTENT(inout), TARGET  :: mem             !> spq memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_snow_spq               ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_spq_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model          => Get_model(model_id)
    table          =  tables(1)
    hgrid          => Get_grid(model%grid_id)
    surface        => Get_vgrid('surface')
    vgrid_soil_sb  => Get_vgrid('soil_layer_sb')
    vgrid_snow_spq => Get_vgrid('snow_layer_spq')


    ! ----------------------------------------------------------------------------------------------------- !
    ! add memory only for LAND & PFT tiles
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      CALL mem%Add_var('soil_depth', mem%soil_depth, &
        & hgrid, surface, &
        & t_cf('soil_depth', 'm', 'Soil depth derived from textures (bedrock), for hydrology'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('root_depth', mem%root_depth, &
        & hgrid, surface, &
        & t_cf('root_depth', 'm', 'actual rooting depth'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_srf_old', mem%t_srf_old, &
        & hgrid, surface, &
        & t_cf('t_srf_old', 'K', 'previous surface temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('zril_old', mem%zril_old, &
        & hgrid, surface, &
        & t_cf('zril_old', '', 'previous timesteps Reynolds number'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('elevation', mem%elevation, &
        & hgrid, surface, &
        & t_cf('elevation', 'm', 'Mean orography'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart=.FALSE., &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var('num_sl_above_bedrock', mem%num_sl_above_bedrock, &
        & hgrid, surface, &
        & t_cf('num_sl_above_bedrock', '[unitless]', 'number of soil layers above bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=0.0_wp)

      CALL mem%Add_var('soil_depth_sl', mem%soil_depth_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_depth_sl', 'm', 'Width of each soil layer that can be saturated with water, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_width_sl', mem%soil_lay_width_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_width_sl', 'm', 'Width of each soil layer that can be saturated with water, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_depth_center_sl', mem%soil_lay_depth_center_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_depth_center_sl', 'm', 'Depth at the center of each soil layer, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_depth_ubound_sl', mem%soil_lay_depth_ubound_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_depth_ubound_sl', 'm', 'Depth at the upper boundary of each soil layer, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_depth_lbound_sl', mem%soil_lay_depth_lbound_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_depth_lbound_sl', 'm', 'Depth at the lower boundary of each soil layer, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('interception', mem%interception, &
        & hgrid, surface, &
        & t_cf('interception', 'kg m-2 s-1', 'surface interception evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('pet', mem%pet, &
        & hgrid, surface, &
        & t_cf('pet', 'kg m-2 s-1', 'potential evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('evaporation', mem%evaporation, &
        & hgrid, surface, &
        & t_cf('evaporation', 'kg m-2 s-1', 'surface evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('evaporation_snow', mem%evaporation_snow, &
        & hgrid, surface, &
        & t_cf('evaporation_snow', 'kg m-2 s-1', 'snow surface evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('srf_runoff', mem%srf_runoff, &
        & hgrid, surface, &
        & t_cf('srf_runoff', 'kg m-2 s-1', 'surface runoff'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('drainage', mem%drainage, &
        & hgrid, surface, &
        & t_cf('drainage', 'kg m-2 s-1', 'drainage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('drainage_fraction', mem%drainage_fraction, &
        & hgrid, surface, &
        & t_cf('drainage_fraction', '1 s-1', 'water mass fraction lost from the soil column by drainage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fact_q_air', mem%fact_q_air, &
        & hgrid, surface, &
        & t_cf('fact_q_air', 'unitless', 'factor representing air at actual water content'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fact_qsat_srf', mem%fact_qsat_srf, &
        & hgrid, surface, &
        & t_cf('fact_qsat_srf', 'unitless', 'factor representing air at saturating water content'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_skin', mem%w_skin, &
        & hgrid, surface, &
        & t_cf('w_skin', 'm', 'Water content in skin reservoir'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root', mem%w_soil_root, &
        & hgrid, surface, &
        & t_cf('w_soil_root', 'm', 'Water content in root zone of the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_fc', mem%w_soil_root_fc, &
        & hgrid, surface, &
        & t_cf('w_soil_root_fc', 'm', 'Water content at field capacity in root zone of the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_pwp', mem%w_soil_root_pwp, &
        & hgrid, surface, &
        & t_cf('w_soil_root_pwp', 'm', 'Water content at permanent wilting point in root zone of the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_soil_root', mem%t_soil_root, &
        & hgrid, surface, &
        & t_cf('t_soil_root', 'K', 'Temperature of fine roots in the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('frac_w_lat_loss_sl', mem%frac_w_lat_loss_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('frac_w_lat_loss_sl', '', 'fraction of lateral (horizontal) water loss of w_soil_sl_old'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('gw_runoff_sl', mem%gw_runoff_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('gw_runoff_sl', 'kg m-2 s-1', 'lateral groundwater runoff per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('percolation_sl', mem%percolation_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('percolation_sl', 'fraction s-1', 'fraction of mass transferred to below layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('drainage_sl', mem%drainage_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('drainage_sl', 'kg m-2 s-1', 'drainage of soil layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_sl', mem%w_soil_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_sl', 'm', 'Water content in soil layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_sat_sl', mem%w_soil_sat_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_sat_sl', 'm', 'Water content of soil layers at saturation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_fc_sl', mem%w_soil_fc_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_fc_sl', 'm', 'Water content of soil layers at field capacity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_pwp_sl', mem%w_soil_pwp_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_pwp_sl', 'm', 'Water content of soil layers at permanent wilting point'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonA', mem%saxtonA, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('saxtonA', '', 'coefficient in moisture-tension relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonB', mem%saxtonB, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('saxtonB', '', 'coefficient in moisture-tension relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonC', mem%saxtonC, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('saxtonC', '', 'exponent of moisture-conductivity relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('kdiff_sat_sl', mem%kdiff_sat_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('kdiff_sat_sl', 'm s-1', 'saturated hydraulic conductivity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_ice_sl', mem%w_ice_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_ice_sl', 'm', 'Soil layer ice content'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sensible_heat_flx', mem%sensible_heat_flx, &
        & hgrid, surface, &
        & t_cf('sensible_heat_flx', 'W m-2', 'sensible heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('latent_heat_flx', mem%latent_heat_flx, &
        & hgrid, surface, &
        & t_cf('latent_heat_flx', 'W m-2', 'latent heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ground_heat_flx', mem%ground_heat_flx, &
        & hgrid, surface, &
        & t_cf('ground_heat_flx', 'W m-2', 'ground heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ground_heat_flx_old', mem%ground_heat_flx_old, &
        & hgrid, surface, &
        & t_cf('ground_heat_flx_old', 'W m-2', 'ground heat flux of previous timestep'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transpiration', mem%transpiration, &
        & hgrid, surface, &
        & t_cf('transpiration', 'kg m-2 s-1', 'transpiration in'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_theta', mem%w_soil_root_theta, &
        & hgrid, surface, &
        & t_cf('w_soil_root_theta', '', 'root zone integrated plant available water, fraction of maximum '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_pot', mem%w_soil_root_pot, &
        & hgrid, surface, &
        & t_cf('w_soil_root_pot', 'MPa', 'soil water potential'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_pot_sl', mem%w_soil_pot_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_pot_sl', 'MPa', 'soil water potential per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_soil_sl', mem%t_soil_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('t_soil_sl', 'K', 'soil temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_freeze_flux', mem%w_soil_freeze_flux, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_freeze_flux', '-', 'Water flux from liquid to ice from freezing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_melt_flux', mem%w_soil_melt_flux, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_melt_flux', '-', 'Water transfer from ice to liquid from melting'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('heat_capa_sl', mem%heat_capa_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('heat_capa_sl', 'J kg-1', 'soil heat capacity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('bulk_dens_sl', mem%bulk_dens_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('bulk_dens_sl', 'kg m-3', 'soil bulk density'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('therm_cond_sl', mem%therm_cond_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('therm_cond_sl', 'W m-2 K-1', 'soil thermal conductivity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sand_sl', mem%sand_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sand_sl', '', 'soil sand proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('silt_sl', mem%silt_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('silt_sl', '', 'soil silt proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('clay_sl', mem%clay_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('clay_sl', '', 'soil clay proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('volume_min_sl', mem%volume_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('volume_min_sl', '', 'soil mineral volumetric proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_org_min_sl', mem%qmax_org_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_org_min_sl', 'mol C kg-1 mineral soil', 'maximum OM sorption capacity of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_po4_min_sl', mem%qmax_po4_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_po4_min_sl', 'mol P kg-1 mineral soil', 'maximum po4 sorption capacity of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_nh4_min_sl', mem%qmax_nh4_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_nh4_min_sl', 'mol N kg-1 mineral soil', 'maximum po4 sorption capacity of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_po4_om_sl', mem%qmax_po4_om_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_po4_om_sl', 'mol P kg-1 SOM', 'maximum po4 sorption capacity of SOM'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('snow_height', mem%snow_height, &
          &  hgrid, surface, &
          &  t_cf('snow_height','m', 'snow depth'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = BASIC, &
          &  loutput = .TRUE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var( 'snow_soil_heat_flux', mem%snow_soil_heat_flux, &
          &  hgrid, surface, &
          &  t_cf('snow_soil_heat_flux','J m-2', 'snow soil heat flux'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_melt_to_soil', mem%snow_melt_to_soil, &
          &  hgrid, surface, &
          &  t_cf('snow_melt_to_soil','m', 'snow melt water to soil'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_srf_heat_flux', mem%snow_srf_heat_flux, &
          &  hgrid, surface, &
          &  t_cf('snow_srf_heat_flux','m', 'snow srf heat flux'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_present_snl', mem%snow_present_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('snow_present_snl','-', 'check if there is snow'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('t_snow_snl', mem%t_snow_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('t_snow_snl','m', 'temperature snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('w_snow_snl', mem%w_snow_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_snow_snl','m', 'water content snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_lay_thickness_snl', mem%snow_lay_thickness_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('snow_lay_thickness_snl','m', 'snow layer thickness'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('w_snow_max_snl', mem%w_snow_max_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_snow_max_snl','m', 'maximum water content snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('w_liquid_snl', mem%w_liquid_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_liquid_snl','m', 'w_liquid_snl'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

    END IF ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)
  END SUBROUTINE Init_spq_memory

#endif
END MODULE mo_spq_memory_class
