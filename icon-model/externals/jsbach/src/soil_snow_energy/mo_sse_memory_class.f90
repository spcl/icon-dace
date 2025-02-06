!> Contains the memory class for the soil and snow energy process.
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
MODULE mo_sse_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes SSE_

  ! Use process configurations
  dsl4jsb_Use_config(SSE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sse_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  TYPE, EXTENDS(t_jsb_memory) :: t_sse_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      & vol_heat_cap,         & !< Volumetric heat capacity of soil solids (from input file) [J m-3 K-1]
      & heat_cond,            & !< Thermal conductivity of soil solids (from input file) [J m-1 s-1 K-1]
      & vol_heat_cap_snow,    & !< (dynamic) Volumetric heat capacity of snow (total depth) [J m-3 K-1]
      & heat_cond_snow,       & !< Thermal conductivity of snow [J m-1 s-1 K-1]
      & thermal_diffusivity     !< Thermal diffusivity of soil [m2 s-1]
    !
    TYPE(t_jsb_var_real3d) :: &
      & vol_heat_cap_sl,      & !< (dynamic) Volumetric heat capacity of dry soil [J m-3 K-1]
      & heat_cond_sl,         & !< (dynamic) Thermal conductivity of soil [J m-1 s-1 K-1]
      & t_soil_sl,            & !< Temperature of soil layers [K]
      & t_soil_acoef,         & !< Soil A coefficient of Richtmeyer and Morton scheme
      & t_soil_bcoef,         & !< Soil B coefficient of Richtmeyer and Morton scheme
      & snow_depth_sl,        & !< Depth of snow layers [m]
      & t_snow,               & !< Temperature of snow layers [K]
      & t_snow_acoef,         & !< Snow A coefficient of Richtmeyer and Morton scheme
      & t_snow_bcoef            !< Snow B coefficient of Richtmeyer and Morton scheme
    TYPE(t_jsb_var_real2d) :: &
      & hcap_grnd,            & !< Heat capacity of the ground (uppermost layer)        [J m-2 K-1]
      & hcap_grnd_old,        & !< Heat capacity of the ground (uppermost_layer), old   [J m-2 K-1]
      & grnd_hflx,            & !< ground heat flux
      & grnd_hflx_old,        & !< ground heat flux (old)
      & thaw_depth              !< Thawing depth [m]

  CONTAINS
    PROCEDURE :: Init => Init_sse_memory
  END TYPE t_sse_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_memory_class'

CONTAINS

  SUBROUTINE Init_sse_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM !, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_sse_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)            :: prefix
    CHARACTER(len=*),    INTENT(in)            :: suffix
    INTEGER,             INTENT(in)            :: lct_ids(:)
    INTEGER,             INTENT(in)            :: model_id

    dsl4jsb_Def_config(SSE_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_e, snow_e      ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sse_memory'

    IF (SIZE(lct_ids) < 1) CALL finish(routine, 'lcd_ids not specified')

    model => Get_model(model_id)

    dsl4jsb_Get_Config(SSE_)

    table = tables(1)

    hgrid => Get_grid(mem%grid_id)

    surface => Get_vgrid('surface')
    soil_e  => Get_vgrid('soil_depth_energy')
    IF( dsl4jsb_Config(SSE_)%l_snow) THEN
      snow_e  => Get_vgrid('snow_depth_energy')
    END IF

    CALL mem%Add_var( 'vol_heat_cap', mem%vol_heat_cap,                 &
      & hgrid, surface,                                                      &
      & t_cf('vol_heat_capacity_soil', 'J m-3 K-1', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'vol_heat_cap_sl', mem%vol_heat_cap_sl,           &
      & hgrid, soil_e,                                                       &
      & t_cf('vol_heat_capacity_soil', 'J m-3 K-1', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'vol_heat_cap_snow', mem%vol_heat_cap_snow,       &
      & hgrid, surface,                                                      &
      & t_cf('vol_heat_capacity_snow', 'J m-3 K-1', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'heat_cond', mem%heat_cond,                       &
      & hgrid, surface,                                                      &
      & t_cf('heat_cond_soil', 'J m-1 s-1 K-1', ''),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'heat_cond_sl', mem%heat_cond_sl,                 &
      & hgrid, soil_e,                                                       &
      & t_cf('heat_cond_soil', 'J m-1 s-1 K-1', ''),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'heat_cond_snow', mem%heat_cond_snow,             &
      & hgrid, surface,                                                      &
      & t_cf('heat_cond_snow', 'J m-1 s-1 K-1', ''),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    IF (      .NOT. dsl4jsb_Config(SSE_)%l_heat_cond_map                     &
      & .AND. .NOT. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      CALL mem%Add_var( 'thermal_diffusivity', mem%thermal_diffusivity,        &
        & hgrid, surface,                                                      &
        & t_cf('thermal_diffusivity_soil', 'm2 s-1', ''),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )
    END IF

    CALL mem%Add_var( 't_soil_sl', mem%t_soil_sl,                            &
      & hgrid, soil_e,                                                       &
      & t_cf('soil_temperature', 'K', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=280.0_wp )

    CALL mem%Add_var( 't_soil_acoef', mem%t_soil_acoef,                      &
      & hgrid, soil_e,                                                       &
      & t_cf('soil_temperature_acoef', '', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 't_soil_bcoef', mem%t_soil_bcoef,                      &
      & hgrid, soil_e,                                                       &
      & t_cf('soil_temperature_bcoef', '', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      CALL mem%Add_var( 'snow_depth_sl', mem%snow_depth_sl,                    &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_depth_sl', 'm', 'Depth of snow layers'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.TRUE.,                                                     &
        & output_level=MEDIUM,                                                 &
        & initval_r=0._wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 't_snow', mem%t_snow,                                  &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_temperature', 'K', ''),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level=BASIC,                                                  &
        & initval_r=273.15_wp )

      CALL mem%Add_var( 't_snow_acoef', mem%t_snow_acoef,                      &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_temperature_acoef', '', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 't_snow_bcoef', mem%t_snow_bcoef,                      &
        & hgrid, snow_e,                                                       &
        & t_cf('snow_temperature_bcoef', '', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )
    END IF

    CALL mem%Add_var( 'hcap_grnd', mem%hcap_grnd,                            &
      & hgrid, surface,                                                      &
      & t_cf('heat_capacity_ground', 'J m-2 K-1', 'Ground heat capacity'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=1.3e5_wp )

    CALL mem%Add_var( 'hcap_grnd_old', mem%hcap_grnd_old,                            &
      & hgrid, surface,                                                              &
      & t_cf('heat_capacity_ground_old', 'J m-2 K-1', 'Ground heat capacity (old)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC,                                                          &
      & initval_r=1.3e5_wp )

    CALL mem%Add_var( 'grnd_hflx', mem%grnd_hflx,                            &
      & hgrid, surface,                                                      &
      & t_cf('grnd_hflx', 'J m-2 s-1', 'Ground heat flux'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=10.0_wp )

    CALL mem%Add_var( 'grnd_hflx_old', mem%grnd_hflx_old,                    &
      & hgrid, surface,                                                      &
      & t_cf('grnd_hflx_old', 'J m-2 s-1', 'Ground heat flux (old)'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=10.0_wp )

    CALL mem%Add_var( 'thaw_depth', mem%thaw_depth,                          &
      & hgrid, surface,                                                      &
      & t_cf('thaw_depth', 'm', 'Thawing depth'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=MEDIUM, &
      & initval_r=10.0_wp, l_aggregate_all=.TRUE. )

  END SUBROUTINE Init_sse_memory

#endif
END MODULE mo_sse_memory_class
