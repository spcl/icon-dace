!> Contains types for the atmosphere-land interface memory.
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
MODULE mo_l2a_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  !USE mo_exception,         ONLY: message
  USE mo_util,              ONLY: One_of

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,     ONLY: LAKE_TYPE

  dsl4jsb_Use_processes CARBON_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_l2a_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 30

  TYPE, EXTENDS(t_jsb_memory) :: t_l2a_memory
    TYPE(t_jsb_var_real2d) :: &
      t_srf       ,   &
      t_eff_srf   ,   &
      qsat_srf    ,   &
      fact_q_air  ,   &
      fact_qsat_srf , &
      zh_srf      ,   &
      zm_srf      ,   &
      alb_vis_dif ,   &
      alb_nir_dif ,   &
      alb_vis_dir ,   &
      alb_nir_dir ,   &
      latent_hflx_wtr    , &
      sensible_hflx_wtr  , &
      latent_hflx_ice    , &
      sensible_hflx_ice  , &
      carbon_conservation_test, &
      yday_c_state_sum
  CONTAINS
    PROCEDURE :: Init => Init_l2a_memory
  END TYPE t_l2a_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_memory_class'

CONTAINS

  SUBROUTINE Init_l2a_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_l2a_memory),      INTENT(inout), TARGET :: mem
    CHARACTER(len=*),         INTENT(in)    :: prefix
    CHARACTER(len=*),         INTENT(in)    :: suffix
    INTEGER,                  INTENT(in)    :: lct_ids(:)
    INTEGER,                  INTENT(in)    :: model_id ! necessary to get modelconfig

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_jsb_grid),        POINTER :: hgrid
    TYPE(t_jsb_vgrid),       POINTER :: surface

    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_l2a_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL mem%Add_var( 't_srf', mem%t_srf,                                    &
      & hgrid, surface,                                                      &
      & t_cf('surface_temperature', 'K', ''),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 't_eff_srf', mem%t_eff_srf,                            &
      & hgrid, surface,                                                      &
      & t_cf('effective_surface_temperature', 'K', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'qsat_srf', mem%qsat_srf,                              &
      & hgrid, surface,                                                      &
      & t_cf('saturated_surface_specific_humidity', '', ''),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'fact_q_air', mem%fact_q_air,                          &
      & hgrid, surface,                                                      &
      & t_cf('factor_air_moisture', '', ''),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'fact_qsat_srf', mem%fact_qsat_srf,                    &
      & hgrid, surface,                                                      &
      & t_cf('factor_qsat_srf', '', ''),                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'zh_srf', mem%zh_srf,                                  &
      & hgrid, surface,                                                      &
      & t_cf('surface_roughness_humidity', '', ''),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'zm_srf', mem%zm_srf,                                  &
      & hgrid, surface,                                                      &
      & t_cf('surface_roughness_momentum', '', ''),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alb_vis_dif', mem%alb_vis_dif,                        &
      & hgrid, surface,                                                      &
      & t_cf('albedo_visible_diffuse', '', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alb_nir_dif', mem%alb_nir_dif,                        &
      & hgrid, surface,                                                      &
      & t_cf('albedo_near_infrared_diffuse', '', ''),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alb_vis_dir', mem%alb_vis_dir,                        &
      & hgrid, surface,                                                      &
      & t_cf('albedo_visible_direct', '', ''),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alb_nir_dir', mem%alb_nir_dir,                        &
      & hgrid, surface,                                                      &
      & t_cf('albedo_near_infrared_dir', '', ''),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%latent_hflx_wtr%ptr)) THEN

      CALL mem%Add_var( 'latent_hflx_wtr', mem%latent_hflx_wtr,                &
        & hgrid, surface,                                                      &
        & t_cf('latent_hflx_wtr', '', 'Latent heat flux over lake water'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'sensible_hflx_wtr', mem%sensible_hflx_wtr,            &
        & hgrid, surface,                                                      &
        & t_cf('sensible_hflx_wtr', '', 'Sensible heat flux over lake water'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'latent_hflx_ice', mem%latent_hflx_ice,                &
        & hgrid, surface,                                                      &
        & t_cf('latent_hflx_ice', '', 'Latent heat flux over lake ice'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'sensible_hflx_ice', mem%sensible_hflx_ice,            &
        & hgrid, surface,                                                      &
        & t_cf('sensible_hflx_ice', '', 'Sensible heat flux over lake ice'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

    END IF

    IF (model%processes(CARBON_)%p%config%active) THEN
      CALL mem%Add_var( 'carbon_conservation_test', mem%carbon_conservation_test,      &
        & hgrid, surface,                                                              &
        & t_cf('carbon_conservation_test', '', 'Carbon conservation test'),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & loutput = .TRUE., output_level=BASIC, lrestart=.TRUE.,                       &
        & initval_r=-999.0_wp )

      CALL mem%Add_var( 'yday_c_state_sum', mem%yday_c_state_sum,                                  &
        & hgrid, surface,                                                                          &
        & t_cf('yday_c_state_sum', '', 'yday c state sum required for carbon conservation test'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & prefix, suffix,                                                                          &
        & output_level=BASIC, loutput = .TRUE., lrestart=.TRUE.,                                   &
        & initval_r=0.0_wp )
    END IF

  END SUBROUTINE Init_l2a_memory

#endif
END MODULE mo_l2a_memory_class
