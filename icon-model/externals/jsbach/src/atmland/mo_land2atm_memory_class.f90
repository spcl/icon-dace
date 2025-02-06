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

  USE mo_kind,                   ONLY: wp
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: LAKE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d

#ifdef __QUINCY_STANDALONE__
#else
  dsl4jsb_Use_processes CARBON_
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_l2a_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 30

  ! ======================================================================================================= !
  !>Type definition for land2atm memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_l2a_memory

    TYPE(t_jsb_var_real2d) :: &
      & t_srf                           , &  !< surface temperature [K]
      & t_eff_srf                       , &  !< effective surface temperature [K]
      & qsat_srf                        , &  !< surface saturated specific humidity [g g-1]
      & zh_srf                          , &  !< roughnesslength for heat transfer [m]
      & zm_srf                          , &  !< roughnesslength for momentum transfer [m]
      & alb_vis_dif                     , &  !< surface albedo for diffuse radiation in the visible band [unitless]
      & alb_nir_dif                     , &  !< surface albedo for diffuse radiation in the near infrared band [unitless]
      & alb_vis_dir                     , &  !< surface albedo for direct radiation in the visible band [unitless]
      & alb_nir_dir                          !< surface albedo for direct radiation in the near infrared band [unitless]

    ! C in jsbach
    TYPE(t_jsb_var_real2d) :: &
      & carbon_conservation_test        , &  !< difference of change in total C sum and Cflux
      & yday_c_state_sum                     !< total C sum

    ! C in Quincy
    TYPE(t_jsb_var_real2d) :: &
      & q_total_c,            & !< total C, sum of all quincy carbon pools
      & q_c_conservation_test   !< difference in quincy total C of last timestep and nbp and total C of this timestep

    ! For lakes
    TYPE(t_jsb_var_real2d) :: &
      & latent_hflx_wtr                 , &  !< latent heat flux over water [W m-2]
      & sensible_hflx_wtr               , &  !< sensible heat flux over water [W m-2]
      & latent_hflx_ice                 , &  !< latent heat flux over ice [W m-2]
      & sensible_hflx_ice                    !< sensible heat flux over ice [W m-2]

  CONTAINS
    PROCEDURE :: Init => Init_l2a_memory
  END TYPE t_l2a_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_l2a_memory_class'

CONTAINS


  ! ======================================================================================================= !
  !>initialize memory (variables) for the process: land2atm
  !>
  SUBROUTINE Init_l2a_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_l2a_memory),       INTENT(inout), TARGET :: mem
    CHARACTER(len=*),          INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),          INTENT(in)            :: suffix          !> tile name
    INTEGER,                   INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,                   INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*),PARAMETER :: routine = TRIM(modname)//':Init_l2a_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(mem%grid_id)
    surface       => Get_vgrid('surface')

    ! ----------------------------------------------------------------------------------------------------- !
    !>add L2A_ variables to memory
    !>
    CALL mem%Add_var('t_srf', mem%t_srf, &
      & hgrid, surface, &
      & t_cf('t_srf', 'K', 'surface temperature'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('t_eff_srf', mem%t_eff_srf, &
      & hgrid, surface, &
      & t_cf('t_eff_srf', 'K', 'effective surface temperature'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('qsat_srf', mem%qsat_srf, &
      & hgrid, surface, &
      & t_cf('qsat_srf', 'g g-1', 'surface saturated specific humidity'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('zh_srf', mem%zh_srf, &
      & hgrid, surface, &
      & t_cf('zh_srf', 'm', 'roughnesslength for heat transfer'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('zm_srf', mem%zm_srf, &
      & hgrid, surface, &
      & t_cf('zm_srf', 'm', 'roughnesslength for momentum transfer'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_vis_dif', mem%alb_vis_dif, &
      & hgrid, surface, &
      & t_cf('alb_vis_dif', 'unitless', 'surface albedo for diffuse radiation in the visible band'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_nir_dif', mem%alb_nir_dif, &
      & hgrid, surface, &
      & t_cf('alb_nir_dif', 'unitless', 'surface albedo for diffuse radiation in the near infrared band'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_vis_dir', mem%alb_vis_dir, &
      & hgrid, surface, &
      & t_cf('alb_vis_dir', 'unitless', 'surface albedo for direct radiation in the visible band'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('alb_nir_dir', mem%alb_nir_dir, &
      & hgrid, surface, &
      & t_cf('alb_nir_dir', 'unitless', 'surface albedo for direct radiation in the near infrared band'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)

#ifndef __NO_QUINCY__
    CALL mem%Add_var('q_total_c', mem%q_total_c, &
      & hgrid, surface, &
      & t_cf('q_total_c', '[mol m-2]', 'total C, sum of all quincy carbon pools'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)

    CALL mem%Add_var('q_c_conservation_test', mem%q_c_conservation_test, &
      & hgrid, surface, &
      & t_cf('q_c_conservation_test', '[mol m-2]', 'difference in change in quincy total C and nbp'), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & output_level = BASIC, &
      & loutput = .TRUE., &
      & lrestart = .FALSE., &
      & initval_r = 0.0_wp)
#endif
#ifdef __QUINCY_STANDALONE__
#else
    ! carbon
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
#endif

    ! lake
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0) THEN

      CALL mem%Add_var('latent_hflx_wtr', mem%latent_hflx_wtr, &
        & hgrid, surface, &
        & t_cf('latent_hflx_wtr', 'W m-2', 'latent heat flux over water'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sensible_hflx_wtr', mem%sensible_hflx_wtr, &
        & hgrid, surface, &
        & t_cf('sensible_hflx_wtr', 'W m-2', 'sensible heat flux over water'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('latent_hflx_ice', mem%latent_hflx_ice, &
        & hgrid, surface, &
        & t_cf('latent_hflx_ice', 'W m-2', 'latent heat flux over ice'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sensible_hflx_ice', mem%sensible_hflx_ice, &
        & hgrid, surface, &
        & t_cf('sensible_hflx_ice', 'W m-2', 'sensible heat flux over ice'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

    END IF ! lake

  END SUBROUTINE Init_l2a_memory

#endif
END MODULE mo_l2a_memory_class
