!> Contains the memory class for the surface.
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
MODULE mo_hd_memory_class
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind, ONLY: wp

  !USE mo_jsb_model_class,  ONLY: t_jsb_model
  !USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,      ONLY: BASIC !, MEDIUM, FULL

  ! Use of prcesses in this module
  !dsl4jsb_Use_processes HD_

  ! Use process configurations
  !dsl4jsb_Use_config(HD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hd_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type for integer array; there's no t_jsb_var_ for integers
  TYPE t_int3d
    INTEGER, POINTER :: ptr(:,:,:)
  END TYPE t_int3d

  TYPE, EXTENDS(t_jsb_memory) :: t_hd_memory

    TYPE(t_jsb_var_real2d) :: &
      hd_mask         , &  !< HD mask
      coast_ocean  ,    &  !< mask for coastal cells of the ocean
      ret_overlflow   , &  !< retention constant of overland flow
      ret_baseflow    , &  !< retention constant of baseflow
      ret_riverflow   , &  !< retention constant of riverflow
      nres_overlflow  , &  !< number of overland flow reservoirs
      nres_baseflow   , &  !< number of baseflow reservoirs
      nres_riverflow  , &  !< number of riverflow reservoirs
      nsplit          , &  !< number of channels river splits into
      outflow_runoff  , &  !< grid cell outflow due to runoff
      outflow_drainage, &  !< grid cell outflow due to drainage
      outflow_rivers  , &  !< grid cell outflow due to rivers
      local_budget    , &  !< total HD reservoir content
      local_fluxes    , &  !< grid cell in- and outflow
      local_wbal_error, &  !< grid cell water balance error
      water_budget    , &  !< HD water budget
      water_budget_change, & !< water budget change within the time step
      outflow_resid   , &  !< residual outflow that's not routed to discharge - only used for debugging (see mo_hd_interface.f90)
      outflow_count        !< count of outflows from a given cell - only used for debugging (see mo_hd_interface.f90)

    TYPE(t_jsb_var_real3d) :: &
      overlflow_res,          &  !< content of overland flow reservoir cascade
      baseflow_res,           &  !< content of baseflow reservoir cascade
      riverflow_res              !< content of groundflow reservoir cascade

    TYPE(t_int3d) ::  &
      nidx_upstream,         &
      bidx_upstream

    ! Diagnostic global land means/sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) :: &
      water_error_gsum           !< Total HD water imbalance during time step       [m3]

    INTEGER ::               &
      nneigh,                &  !< maximum number of neighboring (inflow) cells
      nres_o_max,            &  !< maximum number of overland flow reservoirs
      nres_b_max,            &  !< maximum number of baseflow reservoirs
      nres_r_max                !< maximum number of riverflow reservoirs
  CONTAINS
    PROCEDURE :: Init => Init_hd_memory
  END TYPE t_hd_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_memory_class'

CONTAINS

  SUBROUTINE Init_hd_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_hd_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),   INTENT(in)    :: prefix
    CHARACTER(len=*),   INTENT(in)    :: suffix
    INTEGER,            INTENT(in)    :: lct_ids(:)
    INTEGER,            INTENT(in)    :: model_id

    !dsl4jsb_Def_config(HD_)

    TYPE(t_jsb_grid),  POINTER :: hgrid                      ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, hd_o, hd_b, hd_r  ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hd_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => get_vgrid('surface')

    !vg   wie kann ich die level dynamisch definieren??
    hd_o => Get_vgrid('hd_nres_overlflow')
    hd_b => get_vgrid('hd_nres_baseflow')
    hd_r => get_vgrid('hd_nres_riverflow')

    !dsl4jsb_Get_config(HD_)

    CALL mem%Add_var( 'hd_mask', mem%hd_mask,                                                   &
      & hgrid, surface,                                                                         &
      & t_cf('hd_mask', '', 'HD_ mask: -1 ocean, 0 ocean inflow, 1 land, 2 internal drainage'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
      & prefix, suffix,                                                                         &
      & lrestart=.FALSE.,                                                                       &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'coast_ocean', mem%coast_ocean,                                    &
      & hgrid, surface,                                                                  &
      & t_cf('coast_ocean', '', 'coastal grid cells of the ocean '),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( prefix//'ret_overlflow'//suffix, mem%ret_overlflow,                &
      & hgrid, surface,                                                                  &
      & t_cf('ret_overlflow', 'day', 'overland flow retention constant'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'ret_baseflow', mem%ret_baseflow,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('ret_baseflow', 'day', 'baseflow retention constant'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'ret_riverflow', mem%ret_riverflow,                                &
      & hgrid, surface,                                                                  &
      & t_cf('ret_riverflow', 'day', 'riverflow retention constant'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_overlflow', mem%nres_overlflow,                              &
      & hgrid, surface,                                                                  &
      & t_cf('nres_overlflow', '', 'number of overland flow reservoirs'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_baseflow', mem%nres_baseflow,                                &
      & hgrid, surface,                                                                  &
      & t_cf('nres_baseflow', '', 'number of baseflow reservoirs'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_riverflow', mem%nres_riverflow,                              &
      & hgrid, surface,                                                                  &
      & t_cf('nres_riverflow', '', 'number of riverflow reservoirs'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nsplit', mem%nsplit,                                              &
      & hgrid, surface,                                                                  &
      & t_cf('nsplit', '', 'number of channels river splits into'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'outflow_runoff', mem%outflow_runoff,                              &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_runoff', 'm3/s', 'grid cell outflow due to runoff'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_drainage', mem%outflow_drainage,                          &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_drainage', 'm3/s', 'grid cell outflow due to drainage'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_rivers', mem%outflow_rivers,                              &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_rivers', 'm3/s', 'grid cell outflow due to rivers'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_resid', mem%outflow_resid,                                &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_resid', 'm3/s', 'grid cell outflow residual: should be zero.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_count', mem%outflow_count,                                &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_count', 'm3/s', 'count of outflows from grid cell'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'overlflow_res', mem%overlflow_res,                                &
      & hgrid, hd_o,                                                                     &
      & t_cf('overlflow_res', '', 'content of the overland flow reservoir'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'baseflow_res', mem%baseflow_res,                                  &
      & hgrid, hd_b,                                                                     &
      & t_cf('baseflow_res', '', 'content of the baseflow reservoir'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'riverflow_res', mem%riverflow_res,                                &
      & hgrid, hd_r,                                                                     &
      & t_cf('riverflow_res', '', 'content of the riverflow reservoir'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'local_budget', mem%local_budget,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('local_budget', 'm3/s', 'Total HD reservoir content'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'local_fluxes', mem%local_fluxes,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('local_fluxes', 'm3/s', 'Fluxes for local water balance check'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'local_wbal_error', mem%local_wbal_error,                          &
      & hgrid, surface,                                                                  &
      & t_cf('local_wbal_error', 'm3', 'Local water balance error'),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & loutput=.TRUE., output_level=BASIC,                                              &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget', mem%water_budget,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_budget', 'm3', 'HD water budget'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget_change', mem%water_budget_change,                    &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_budget_change', 'm3/timestep', 'HD water budget change within time step'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

#ifdef __ICON__
    ! Diagnostic 1d global land variables for experiment monitoring
    !       (1d stream variables are not supported with echam)
    !
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var( 'water_error_gsum', mem%water_error_gsum,                          &
        & hgrid, surface,                                                                  &
        & t_cf('hd_water_error_gsum', 'm3', 'Total HD water imbalance during time step'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix,                                                                  &
        & lrestart=.FALSE.,                                                                &
        & initval_r=0.0_wp )
    ENDIF
#endif

  END SUBROUTINE Init_hd_memory

#endif
END MODULE mo_hd_memory_class
