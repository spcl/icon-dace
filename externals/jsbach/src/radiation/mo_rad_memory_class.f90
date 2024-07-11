!> Contains the memory class for radiation.
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
MODULE mo_rad_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_util,             ONLY: One_of
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,    ONLY: BARE_TYPE, VEG_TYPE, LAND_TYPE, GLACIER_TYPE, LAKE_TYPE

  dsl4jsb_Use_processes SEB_
  dsl4jsb_Use_config(SEB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_rad_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 90

  TYPE, EXTENDS(t_jsb_memory) :: t_rad_memory

    TYPE(t_jsb_var_real2d)  :: &
      & sw_srf_net,            & !< net shortwave radiation at surface [W m-2]
      & swvis_srf_net,         & !< net shortwave radiation at surface in visible range [W m-2]
      & swnir_srf_net,         & !< net shortwave radiation at surface in NIR range [W m-2]
      & lw_srf_net,            & !< net longwave radiation at surface [W m-2]
      & rad_srf_net              !< net radiation at surface [W m-2]

    TYPE(t_jsb_var_real2d) :: &
      & alb_background,        & !< Background albedo (without snow or canopy) of tile
                                 !< Note, alb_background is ONLY used for FUNCTION get_surface_albedo_simple
      & alb_vis,               & !< Albedo of the whole tile in the VIS (visible) range [unitless]
      & alb_nir,               & !< Albedo of the whole tile in the NIR (near-infrared) range [unitless]
      & alb_vis_snow,          & !< Albedo of snow in the VIS range [unitless]
      & alb_nir_snow,          & !< Albedo of snow in the NIR range [unitless]
      & alb_vis_pond,          & !< Albedo of ponds in the VIS range [unitless]
      & alb_nir_pond             !< Albedo of ponds in the NIR range [unitless]


    !
    ! Additional variables for soil
    TYPE(t_jsb_var_real2d) :: &
      & alb_vis_soil,         & ! Albedo in VIS range of soil
                                ! Initialized from bc_land_phys.nc file.
                                ! If use_albedo_soil_scheme =
                                !   .TRUE.:  then it is calculated and may incl. litter layer and carbon within soil dependent
                                !            on the use_alb_soil_organic_C and use_alb_soil_litter switches.
                                !   .FALSE.: it is not calculated but remains as taken from bc_land_phys.nc file during
                                !            initialization. alb_vis_soil of the bc_land_phys.nc file represents the soil albedo
                                !            including carbon in the soil as well as litter on the soil, but without the canopy
                                !            above. It was calculated from MODIS data (aboard the TERRA satellite) using linear
                                !            regression to project the albedo when fapar=0 (LAI/Canopy =0) (for details see
                                !            JSBACH3 documentation 4.1.1 or ask Thomas Raddatz).
      & alb_nir_soil,         & ! Same as alb_vis_soil but for NIR.
      ! R: Albedo of mineral soil is needed when Freya scheme is used. Not in the bc_lnd_phy.nc file yet!
      & alb_vis_mineralsoil,  & !< Albedo of soil without carbon or litter.
                                !! from bc_land_phys.nc file, necessary when use_albedo_soil_scheme = .TRUE.
                                !! for SOC-dependent albedo calculation
      & alb_nir_mineralsoil,  & ! from bc_land_phys.nc file, necessary when use_albedo_soil_scheme = .TRUE. for
                                ! SOC dependend albedo calculation
      & alb_vis_lnd,          & ! Albedo on the non-lake part of tile in the VIS (visible) range
      & alb_nir_lnd             ! Albedo of the non-lake part of tile in the NIR (near-infrared) range

    !
    ! Additional variables for vegetation
    TYPE(t_jsb_var_real2d) :: &
      & alb_vis_can,          & ! Albedo of canopy (without snow cover) in VIS range.
                                ! Never calculated. Initialized from bc_land_phys.nc file.
                                ! If use_alb_canopy =
                                !   .TRUE.:  remains as taken from bc_land_phys.nc file during initialization. These data were
                                !            calculated from MODIS data (aboard the TERRA satellite) using linear regression to project
                                !            the albedo when fapar=1 (LAI/Canopy >=1) (for details see JSBACH3 documentation 4.1.1
                                !            or ask Thomas Raddatz).
                                !   .FALSE.: taken from lctlib file (AlbedoCanopyVIS)
      & alb_nir_can,          & !< Albedo of canopy (without snow cover) in NIR range. Is never calculated:
                                !< If use_alb_canopy=.TRUE. taken from bc_land_phys.nc file, else from lctlib-File (AlbedoCanopyNIR)
      & sky_view_fract,       & ! Fraction of bare ground below vegetation (without accounting for stems) that is "visible" for the sun
      & sky_view_fract_stem,  & ! Fraction of soil below vegetation which is (for the sun) covered by stems but not by LAI
      & snow_age,             &
      !
      & par,                  & ! sw_par_down from atmosphere in [W/m^2]
      & par_down_mol,         & ! Downward PAR flux in mol(photons)/(m^2 s) but only for current time step/"delta_time"
      & faPAR,                & ! fraction of PAR absorbed by the canopy and per leaf area
      ! & par_down_mol_nacc,    & ! incoming PAR accumulated (over output timestep) in mol(photons)/(m^2)
      ! & par_down_mol_tavg,    & ! incoming PAR (mean value over output timestep) [MOL PHOTONS/M^2 S]
      & soil_reflectivity_par,& ! Soil reflectivity in the PAR region but only for current time step/"delta_time"
      & fract_par_direct        ! Fraction of direct PAR
      ! & apar_nacc

    TYPE(t_jsb_var_real3d) :: & ! alt: (nland, ncanopy, ntiles)
      & lai_cl,               & ! LAI per canopy layer; preliminarily, should later be put to another process...
      & faPAR_cl,             & ! Absorbed PAR per canopy layer and per leaf area [fraction of irradiance]
                                ! As fraction of incoming radiation at the top of the canopy layers (=1).
      & apar_per_lai_cl         ! Absorbed PAR per canopy layer LAI (lai_cl). Only for current time step/dtime.
                                ! [(absorbed photons) / (m^2(leaf area) * s)]
      ! & apar_per_lai_cl_nacc, & ! PAR absorbed by the canopy accumulated (over output timestep) in mol(photons)/m^2
      ! & apar_per_lai_cl_tavg    ! PAR absorbed by the canopy (mean value over output timestep)  in mol(photons)/(m^2 * s)

    !
    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & rad_net_lwtr,         & !< Net radiation over lake water
      & rad_net_lice,         & !< Net radiation over lake ice
      & albedo_lwtr,          & !< Albedo of open water on lake
      & albedo_lice             !< Albedo of ice on lake

  CONTAINS
    PROCEDURE :: Init => Init_rad_memory
  END TYPE t_rad_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: radiation
  !>
  SUBROUTINE Init_rad_memory(mem, prefix, suffix, lct_ids, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    USE mo_rad_constants,     ONLY: AlbedoLakeWater, AlbedoLakeIceMin

    CLASS(t_rad_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)            :: prefix
    CHARACTER(len=*),    INTENT(in)            :: suffix
    INTEGER,             INTENT(in)            :: lct_ids(:)
    INTEGER,             INTENT(in)            :: model_id

    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER  :: model                        !< model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface, canopy_layer        !< Vertical grids
    INTEGER                     :: table                        !< ...

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_rad_memory'

    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model        => Get_model(model_id)
    table        = tables(1)
    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
    canopy_layer => Get_vgrid('canopy_layer')

    dsl4jsb_Get_config(SEB_)

    CALL mem%Add_var( 'sw_srf_net', mem%sw_srf_net,                        &
      & hgrid, surface,                                                    &
      & t_cf('sw_srf_net', '', 'Net surface shortwave radiation'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE.,                                                  &
      & output_level=BASIC,                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'swvis_srf_net', mem%swvis_srf_net,                    &
      & hgrid, surface,                                                      &
      & t_cf('swvis_srf_net', '', 'Net surface shortwave radiation in VIS'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=MEDIUM,                                                 &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'swnir_srf_net', mem%swnir_srf_net,                    &
      & hgrid, surface,                                                      &
      & t_cf('swnir_srf_net', '', 'Net surface shortwave radiation in NIR'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=MEDIUM,                                                 &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'lw_srf_net', mem%lw_srf_net,                        &
      & hgrid, surface,                                                    &
      & t_cf('lw_srf_net', '', 'Net surface longwave radiation'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE.,                                                  &
      & output_level=BASIC,                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'rad_srf_net', mem%rad_srf_net,                      &
      & hgrid, surface,                                                    &
      & t_cf('rad_srf_net', '', 'Net surface radiation'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE.,                                                  &
      & output_level=BASIC,                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alb_background', mem%alb_background,                &
      & hgrid, surface,                                                    &
      & t_cf('surface_albedo_background', '', ''),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE.,                                                  &
      & output_level=MEDIUM,                                               &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'alb_vis', mem%alb_vis,                              &
      & hgrid, surface,                                                    &
      & t_cf('gridcell_surface_albedo_visible', '', ''),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.TRUE.,                                                   &
      & output_level=BASIC,                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alb_nir', mem%alb_nir,                              &
      & hgrid, surface,                                                    &
      & t_cf('gridcell_surface_albedo_nir', '', ''),                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.TRUE.,                                                   &
      & output_level=BASIC,                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alb_vis_snow', mem%alb_vis_snow,                    &
      & hgrid, surface,                                                    &
      & t_cf('albedo_vis_snow', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & output_level=MEDIUM,                                               &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'alb_nir_snow', mem%alb_nir_snow,                    &
      & hgrid, surface,                                                    &
      & t_cf('albedo_nir_snow', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & output_level=MEDIUM,                                               &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)

    ! Additional variables for land tiles
    IF ( (One_of(BARE_TYPE, lct_ids(:)) > 0 .OR. One_of(VEG_TYPE, lct_ids(:)) > 0           &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. One_of(GLACIER_TYPE, lct_ids(:)) > 0) &
      & ) THEN

      CALL mem%Add_var( 'alb_vis_lnd', mem%alb_vis_lnd,                             &
        & hgrid, surface,                                                           &
        & t_cf('surface_albvis_lnd', '', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.TRUE.,                                                          &
        & output_level=MEDIUM, &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'alb_nir_lnd', mem%alb_nir_lnd,                             &
        & hgrid, surface,                                                           &
        & t_cf('surface_albnir_lnd', '', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.TRUE.,                                                          &
        & output_level=MEDIUM, &
        & initval_r=0.0_wp )

    END IF ! land

    ! Additional variables for soil
    IF ( (One_of(BARE_TYPE, lct_ids(:)) > 0 .OR. One_of(VEG_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0)                                 &
      & ) THEN

      CALL mem%Add_var( 'alb_vis_soil', mem%alb_vis_soil,                    &
        & hgrid, surface,                                                    &
        & t_cf('soil_albedo_with_carbon_and_litter', '', ''),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=MEDIUM,                                               &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

      CALL mem%Add_var( 'alb_nir_soil', mem%alb_nir_soil,                    &
        & hgrid, surface,                                                    &
        & t_cf('soil_albedo_with_carbon_and_litter', '', ''),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=MEDIUM,                                               &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

      CALL mem%Add_var( 'alb_vis_mineralsoil', mem%alb_vis_mineralsoil,      &
        & hgrid, surface,                                                    &
        & t_cf('alb_vis_mineralsoil', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=MEDIUM,                                               &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var( 'alb_nir_mineralsoil', mem%alb_nir_mineralsoil,      &
        & hgrid, surface,                                                    &
        & t_cf('alb_nir_mineralsoil', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=MEDIUM,                                               &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var( 'alb_vis_pond', mem%alb_vis_pond,                    &
        & hgrid, surface,                                                    &
        & t_cf('alb_vis_pond', '', ''),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.TRUE.,                                                   &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'alb_nir_pond', mem%alb_nir_pond,                    &
        & hgrid, surface,                                                    &
        & t_cf('alb_nir_pond', '', ''),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.TRUE.,                                                   &
        & initval_r=0.0_wp )
    END IF

    ! Additional variables for vegetation
    IF ( (One_of(VEG_TYPE, lct_ids(:)) > 0 .OR. One_of(LAND_TYPE, lct_ids(:)) > 0) &
      & ) THEN

      CALL mem%Add_var( 'alb_vis_can', mem%alb_vis_can,                      &
        & hgrid, surface,                                                    &
        & t_cf('albedo_vis_canopy', '', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE.,                                                  &
        & output_level=MEDIUM,                                               &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'alb_nir_can', mem%alb_nir_can,                      &
        & hgrid, surface,                                                    &
        & t_cf('albedo_nir_canpy', '', ''),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE.,                                                  &
        & output_level=MEDIUM,                                               &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'snow_age', mem%snow_age,                            &
        & hgrid, surface,                                                    &
        & t_cf('snow_age', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

      CALL mem%Add_var( 'sky_view_fract', mem%sky_view_fract,                &
        & hgrid, surface,                                                    &
        & t_cf('sky_view_fract', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var( 'sky_view_fract_stem', mem%sky_view_fract_stem,      &
        & hgrid, surface,                                                    &
        & t_cf('sky_view_fract_stem', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var( 'par', mem%par,                                      &
        & hgrid, surface,                                                    &
        & t_cf('Downward PAR flux in [W/m^2]', '', ''),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL mem%Add_var( 'par_down_mol', mem%par_down_mol,                    &
        & hgrid, surface,                                                    &
        & t_cf('Downward PAR flux in mol (photons)', '', ''),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! R: initval =0

      CALL mem%Add_var( 'fapar', mem%faPAR,                                  &
        & hgrid, surface,                                                    &
        & t_cf('Fraction of PAR absorbed by the canopy', '', ''),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=BASIC,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! CALL mem%Add_var( 'par_down_mol_nacc', mem%par_down_mol_nacc,          &
      !   & hgrid, surface,                                                    &
      !   & t_cf('Incoming PAR accumulated (over output timestep)', '', ''),   &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      !   & prefix, suffix,                                                    &
      !   & lrestart=.TRUE., initval_r=0.0_wp )

      ! CALL mem%Add_var( 'par_down_mol_tavg', mem%par_down_mol_tavg,          &
      !   & hgrid, surface,                                                    &
      !   & t_cf('Incoming PAR (mean value over output timestep)', '', ''),    &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      !   & prefix, suffix,                                                    &
      !   & lrestart=.FALSE., initval_r=0.0_wp ) ! R: initval =0

      CALL mem%Add_var( 'soil_reflectivity_par', mem%soil_reflectivity_par,      &
        & hgrid, surface,                                                        &
        & t_cf('Soil reflectivity in the PAR region', '', ''),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp )

      CALL mem%Add_var( 'fract_par_direct', mem%fract_par_direct,            &
        & hgrid, surface,                                                    &
        & t_cf('Fraction of direct PAR', '', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.FALSE., initval_r=0.5_wp )

      CALL mem%Add_var('lai_cl', mem%lai_cl,                                 &
        & hgrid, canopy_layer,                                               &
        & t_cf('LAI for each canopy layer', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL mem%Add_var('fapar_cl', mem%faPAR_cl,                                  &
        & hgrid, canopy_layer,                                                    &
        & t_cf('Absorbed PAR per canopy layer [fraction of irradiance]', '', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
        & prefix, suffix,                                                         &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL mem%Add_var('apar_per_lai_cl', mem%apar_per_lai_cl,               &
        & hgrid, canopy_layer,                                               &
        & t_cf('Absorbed PAR per canopy layer (LAI area).', '', ''),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp )

      !  CALL mem%Add_var('apar_per_lai_cl_nacc', mem%apar_per_lai_cl_nacc,                 &
      !    & hgrid, canopy_layer,                                                           &
      !    & t_cf('PAR absorbed by the canopy (accumulated over output timestep)', '', ''), &
      !    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      !    & prefix, suffix,                                                                &
      !    & lrestart=.TRUE., initval_r=0.0_wp )

      !  CALL mem%Add_var('apar_per_lai_cl_tavg', mem%apar_per_lai_cl_tavg,                &
      !    & hgrid, canopy_layer,                                                          &
      !    & t_cf('PAR absorbed by the canopy (mean value over output timestep)', '', ''), &
      !    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
      !    & prefix, suffix,                                                               &
      !    & lrestart=.FALSE., initval_r=0.0_wp )

      !  CALL mem%Add_var('apar_nacc', mem%apar_nacc,                                       &
      !    & hgrid, surface,                                                                &
      !    & t_cf('PAR absorbed by the canopy (accumulated over output timestep)', '', ''), &
      !    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      !    & prefix, suffix,                                                                &
      !    & lrestart=.TRUE., initval_r=0.0_wp )
    END IF

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 'rad_net_lwtr', mem%rad_net_lwtr,                     &
        & hgrid, surface,                                                     &
        & t_cf('rad_net_lwtr', '', 'Net surface radiation over lake water'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & lrestart=.FALSE.,                                                   &
        & output_level=MEDIUM,                                                &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'albedo_lwtr', mem%albedo_lwtr,                       &
        & hgrid, surface,                                                     &
        & t_cf('albedo_lwtr', '', 'Albedo of open water on lakes'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & output_level=MEDIUM,                                                &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=AlbedoLakeWater )

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN

        CALL mem%Add_var( 'rad_net_lice', mem%rad_net_lice,                     &
          & hgrid, surface,                                                     &
          & t_cf('rad_net_lice', '', 'Net surface radiation over lake ice'),    &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
          & prefix, suffix,                                                     &
          & lrestart=.FALSE.,                                                   &
          & output_level=MEDIUM,                                                &
          & initval_r=0.0_wp )

        CALL mem%Add_var( 'albedo_lice', mem%albedo_lice,                       &
          & hgrid, surface,                                                     &
          & t_cf('albedo_lice', '', 'Albedo of ice on lake water'),             &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
          & prefix, suffix,                                                     &
          & output_level=MEDIUM,                                                &
          & loutput=.TRUE., lrestart=.TRUE., initval_r=AlbedoLakeIceMin )

      END IF ! l_ice_on_lakes
    END IF ! lakes

  END SUBROUTINE Init_rad_memory

#endif
END MODULE mo_rad_memory_class
