! Provide an implementation of the initialization of the sea-ice.
!
! Provide an implementation of the subroutines used for initialization of
! the surface module.
!
!
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

MODULE mo_ice_init_thermo
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_dynamics_config,     ONLY: nold
  USE mo_coupling_config,     ONLY: is_coupled_to_atmo
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D !, t_patch_vert
  USE mo_exception,           ONLY: finish, message
  USE mo_impl_constants,      ONLY: success, max_char_length, sea_boundary

  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, ki, ks, Tf, mu
  USE mo_ocean_nml,           ONLY: no_tracer
  USE mo_sea_ice_nml,         ONLY: i_ice_dyn, i_ice_advec, &
    &                                use_IceInitialization_fromTemperature, use_constant_tfreez, &
    &                               init_analytic_conc_param, init_analytic_hi_param, &
    &                               init_analytic_hs_param, init_analytic_temp_under_ice, &
    &                               albi, albedoW_sim, initialize_seaice_fromfile, &
    &                               i_ice_therm
  USE mo_ocean_nml,           ONLY: limit_seaice, seaice_limit
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_state,         ONLY: v_base, ocean_restart_list, ocean_default_list
  USE mo_var_list,            ONLY: add_var, t_var_list_ptr
  USE mo_var_groups,          ONLY: groups
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, GRID_UNSTRUCTURED
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL,      &
    &                               GRID_UNSTRUCTURED_VERT, GRID_VERTEX,    &
    &                               GRID_UNSTRUCTURED_EDGE, GRID_EDGE
  USE mo_zaxis_type,          ONLY: ZA_GENERIC_ICE, ZA_SURFACE
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_io_config,           ONLY: lnetcdf_flt64_output

  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes, t_sea_ice_budgets

  USE mo_ice_fem_init,            ONLY: ice_init_fem
  USE mo_ice_fem_icon_init,        ONLY: init_fem_wgts, init_fem_wgts_extra, destruct_fem_wgts,  &
    &                               ice_fem_grid_init, ice_fem_grid_post
  USE mo_ocean_initial_conditions, ONLY: init_cell_2D_variable_fromFile

#include "add_var_acc_macro.inc"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ice_init
  PUBLIC :: construct_sea_ice, destruct_sea_ice
  PUBLIC :: construct_atmos_fluxes

  CHARACTER(len=12)           :: str_module    = 'IceInitTherm'  ! Output of module for 1 line debug

CONTAINS


  !-------------------------------------------------------------------------
  !
  !> ice_init
  !!
  SUBROUTINE ice_init( patch_3D, p_os, ice, cellThicknessUnderIce)
    TYPE(t_patch_3D), TARGET              :: patch_3D
    TYPE(t_hydro_ocean_state)             :: p_os
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    REAL(wp), DIMENSION(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks), &
      & INTENT(OUT)                       :: cellThicknessUnderIce

    !local variables
    REAL(wp), DIMENSION(nproma,ice%kice, patch_3D%p_patch_2D(1)%alloc_cell_blocks) :: &
      & Tinterface  ! temperature at snow-ice interface

    TYPE(t_patch), POINTER                :: p_patch
    LOGICAL  :: has_missValue
    REAL(wp) :: missValue, hi_max
    REAL(wp) :: read_in(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    !INTEGER i,j,k      ! counter for loops
    INTEGER k !, jb, jc, i_startidx_c, i_endidx_c! counter for loops
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_init_thermo:ice_init'

  INTEGER                     :: idt_src       = 1         ! Level of detail for 1 line debug
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch => patch_3D%p_patch_2D(1)

    ! Calculate the sea surface freezing temperature                        [C]
    IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
      ice%Tfw(:,:) = Tf
    ELSE
      ice%Tfw(:,:) = -mu * p_os%p_prog(nold(1))%tracer(:,1,:,2)
    ENDIF

    ice% Tsurf(:,:,:)  = Tf
    ice% T1   (:,:,:)  = Tf
    ice% T2   (:,:,:)  = Tf
    ice% conc (:,:,:)  = 0.0_wp

    IF (initialize_seaice_fromfile) THEN
      hi_max = seaice_limit * v_base%del_zlev_m(1)

      CALL init_cell_2D_variable_fromFile(patch_3d, read_in, "seaice_hi_2D", has_missValue, missValue)
      IF (limit_seaice) THEN
        ice%hi(:,1,:) = MIN(read_in(:,:),hi_max)
      ELSE
         ice%hi(:,1,:) = read_in(:,:)
      ENDIF
      CALL init_cell_2D_variable_fromFile(patch_3d, read_in, "seaice_hs_2D", has_missValue, missValue)
      ice%hs(:,1,:) = read_in(:,:)
      CALL init_cell_2D_variable_fromFile(patch_3d, read_in, "seaice_conc_2D", has_missValue, missValue)
      ice%conc(:,1,:) = read_in(:,:)

    ELSE
      ! Stupid initialisation trick for Levitus initialisation
      IF (use_IceInitialization_fromTemperature) THEN
        WHERE (p_os%p_prog(nold(1))%tracer(:,1,:,1) <= init_analytic_temp_under_ice .and. v_base%lsm_c(:,1,:) <= sea_boundary )
          ice%hi(:,1,:)   = init_analytic_hi_param
          ice%hs(:,1,:)   = init_analytic_hs_param
          ice%conc(:,1,:) = init_analytic_conc_param
        ENDWHERE
      ! or constant initialization for ice, snow and concentration
      ELSE
        WHERE (v_base%lsm_c(:,1,:) <= sea_boundary )
          ice%hi(:,1,:)    = init_analytic_hi_param
          ice%hs(:,1,:)    = init_analytic_hs_param
          ice%conc(:,1,:)  = init_analytic_conc_param
        ENDWHERE
      ENDIF
    ENDIF


    DO k=1,ice%kice
      WHERE(ice% hi(:,k,:) > 0.0_wp)
        ice% Tsurf (:,k,:) = ice%Tfw(:,:)
        ice% T1    (:,k,:) = ice%Tfw(:,:)
        ice% T2    (:,k,:) = ice%Tfw(:,:)
        Tinterface (:,k,:) = (ice%Tfw(:,:) * (ki/ks * ice%hs(:,k,:)/ice%hi(:,k,:)) &
            &                + ice%Tsurf(:,k,:)) / (1.0_wp+ki/ks * ice%hs(:,k,:)/ice%hi(:,k,:))
        ice% conc  (:,k,:) = ice%conc(:,k,:)/REAL(ice%kice,wp)
        ice% T1    (:,k,:) = ice%Tfw(:,:) + 2._wp/3._wp*(Tinterface(:,k,:)-ice%Tfw(:,:))
        ice% T2    (:,k,:) = ice%Tfw(:,:) + 1._wp/3._wp*(Tinterface(:,k,:)-ice%Tfw(:,:))
        ice% draft (:,k,:) = (rhos * ice%hs(:,k,:) + rhoi * ice%hi(:,k,:)) / rho_ref
      END WHERE
    ENDDO

    ice%concSum(:,:)   = SUM(ice%conc(:,:,:), 2)
    ice%draftave (:,:) = sum(ice%draft(:,:,:) * ice%conc(:,:,:),2)
    ice%zUnderIce(:,:) = v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) - ice%draftave(:,:)

    cellThicknessUnderIce (:,:) = ice%zUnderIce(:,:)

    IF ( i_ice_dyn == 1 ) THEN ! AWI dynamics
      CALL ice_fem_grid_init(patch_3D)
      CALL ice_init_fem
      CALL ice_fem_grid_post(p_patch)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('IceInit: hi       ' ,ice%hi       ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: hs       ' ,ice%hs       ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: conc     ' ,ice%conc     ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: draft    ' ,ice%draft    ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: draftave ' ,ice%draftave ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: zUnderIce' ,ice%zUnderIce,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: Tfw      ' ,ice%Tfw      ,str_module, idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE ice_init

  !-------------------------------------------------------------------------
  !
  !> Constructor of sea-ice model, allocates all components and assigns zero.
  !!
  SUBROUTINE construct_sea_ice(patch_3D, p_ice, i_no_ice_thick_class)
    TYPE(t_patch_3D),TARGET,INTENT(IN)    :: patch_3D
    TYPE (t_sea_ice),       INTENT(INOUT) :: p_ice
    INTEGER,                INTENT(IN)    :: i_no_ice_thick_class

    !Local variables
    INTEGER :: alloc_cell_blocks, nblks_v, nblks_e, ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_init_thermo:construct_sea_ice'

    TYPE(t_patch),POINTER    :: p_patch
    INTEGER                  :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch           => patch_3D%p_patch_2D(1)
    alloc_cell_blocks =  p_patch%alloc_cell_blocks
    nblks_v           =  p_patch%nblks_v
    nblks_e           =  p_patch%nblks_e

    p_ice%kice = i_no_ice_thick_class

    !$ACC ENTER DATA COPYIN(p_ice)

    CALL add_var(ocean_restart_list, 'hi', p_ice%hi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hi', 'm', 'ice thickness', datatype_flt),&
      &          grib2_var(10, 2, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"), lopenacc=.TRUE.)
    __acc_attach(p_ice%hi)

    CALL add_var(ocean_restart_list, 'hs', p_ice%hs ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hs', 'm', 'snow thickness', datatype_flt),&
      &          grib2_var(10, 2, 192, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"), lopenacc=.TRUE.)
    __acc_attach(p_ice%hs)

    CALL add_var(ocean_restart_list, 'conc', p_ice%conc ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('conc', '', 'ice concentration in each ice class', datatype_flt),&
      &          grib2_var(10, 2, 0, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"), lopenacc=.TRUE.)
    __acc_attach(p_ice%conc)
!fixme why is this needed
    CALL add_var(ocean_restart_list, 'concSum', p_ice%concSum ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('concSum', '', 'total ice concentration', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), lopenacc=.TRUE.)
    __acc_attach(p_ice%concSum)

    CALL add_var(ocean_default_list, 'vol', p_ice%vol ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('vol', 'm^3', 'ice volume', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%vol)

    CALL add_var(ocean_default_list, 'vols', p_ice%vols ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('vols', 'm^3', 'snow volume', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%vols)

    CALL add_var(ocean_default_list, 'delhi', p_ice%delhi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('delhi', 'm', 'Change in ice mean thickness due to thermodynamic effects', datatype_flt),&
      &          grib2_var(10, 2, 6, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%delhi)

    CALL add_var(ocean_default_list, 'delhs', p_ice%delhs ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('delhs', 'm', 'Change in mean snow thickness due to thermodynamic melting', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%delhs)

      ! functionality of hiold and hsold was replaced by delhi and delhs
      ! to be deleted together with the old thermodynamics routines
    CALL add_var(ocean_default_list, 'hiold', p_ice%hiold ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hiold', 'm', 'ice thickness (last timstep)', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%hiold)

    CALL add_var(ocean_default_list, 'hsold', p_ice%hsold ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hsold', 'm', 'snow thickness (last timstep)', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%hsold)
!fixme
    ! thermodynamics, fast
    CALL add_var(ocean_restart_list, 'Tsurf', p_ice%Tsurf ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Tsurf', '', 'surface temperature', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%Tsurf)
!fixme not needed when tfw = tfreeze
    CALL add_var(ocean_restart_list, 'Tfw', p_ice%Tfw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Tfw', '', 'ocean surface freezing temperature', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%Tfw)

    CALL add_var(ocean_restart_list, 'Qtop', p_ice%Qtop ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Qtop', 'W/m^2', 'Energy flux available for surface melting', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%Qtop)

    CALL add_var(ocean_restart_list, 'Qbot', p_ice%Qbot ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Qbot', 'W/m^2', 'Conductive heat flux at ice-ocean interface', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%Qbot)

    CALL add_var(ocean_default_list, 'alb', p_ice%alb ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('alb', '', 'albedo of snow-ice system', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%alb)

    ! thermodynamics, slow
    CALL add_var(ocean_default_list, 'Qbot_slow', p_ice%Qbot_slow ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Qbot_slow', 'W/m^2', 'Energy flux at ice-ocean interface', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%Qbot_slow)

    CALL add_var(ocean_default_list, 'zHeatOceI', p_ice%zHeatOceI,GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('zHeatOceI', 'W/m^2', 'Oceanic Heat flux', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%zHeatOceI)

    CALL add_var(ocean_default_list, 'heatOceI', p_ice%heatOceI ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('heatOceI', 'W/m^2', 'Heat flux to ocean from the ice growth', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%heatOceI)

    CALL add_var(ocean_default_list, 'heatOceW', p_ice%heatOceW ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('heatOceW', 'W/m^2', 'Heat flux to ocean from the atmosphere', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%heatOceW)

    CALL add_var(ocean_default_list, 'snow_to_ice', p_ice%snow_to_ice ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('snow_to_ice', 'm', 'amount of snow that is transformed to ice', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%snow_to_ice)

    CALL add_var(ocean_default_list, 'newice', p_ice%newice ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('newice', 'm', 'new ice growth in open water', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%newice)

    CALL add_var(ocean_default_list, 'totalsnowfall', p_ice%totalsnowfall ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('totalsnowfall', 'm', 'total snow fall on sea ice', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%totalsnowfall)
!fixme
    CALL add_var(ocean_restart_list, 'draft', p_ice%draft ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('draft', 'm', 'water equiv. of ice and snow on ice covered area', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%draft)
!fixme
    CALL add_var(ocean_restart_list, 'draftave', p_ice%draftave ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('draftave', 'm', 'avrg. water equiv. of ice and snow on grid area', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%draftave)
!fixme
    CALL add_var(ocean_restart_list, 'zUnderIce', p_ice%zUnderIce ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('zUnderIce', 'm', 'water in upper ocean grid cell below ice', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%zUnderIce)
!fixme is this usefull ?
    CALL add_var(ocean_default_list, 'draft_old', p_ice%draftave_old ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('draftave_old', 'm', 'avrg. water equiv. of ice and snow on grid area', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%draftave_old)

    !fixme winton scheme -- currently not functional
    IF ( i_ice_therm == 2 ) THEN
    ! thermodynamics, fast (winton scheme only -- not currently functional)
    CALL add_var(ocean_restart_list, 'T1', p_ice%T1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T1', 'C', 'Temperature upper layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%T1)

    CALL add_var(ocean_restart_list, 'T2', p_ice%T2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T2', 'C', 'Temperature lower layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%T2)

    CALL add_var(ocean_restart_list, 'E1', p_ice%E1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E1', 'Jm/kg', 'Energy content upper layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%E1)

    CALL add_var(ocean_restart_list, 'E2', p_ice%E2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E2', 'Jm/kg', 'Energy content lower layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%E2)

    CALL add_var(ocean_default_list, 'surfmelt', p_ice%surfmelt ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('surfmelt', 'm', 'surface melt water running into ocean', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%surfmelt)

    CALL add_var(ocean_default_list, 'surfmeltT', p_ice%surfmeltT ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('surfmeltT', 'C', 'Mean temperature of surface melt water', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%surfmeltT)

    CALL add_var(ocean_default_list, 'evapwi', p_ice%evapwi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('evapwi', 'kg m-2', 'amount of evaporated water if no ice left', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%evapwi)
  ELSE
    CALL add_var(ocean_default_list, 'T1', p_ice%T1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T1', 'C', 'Temperature upper layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%T1)

    CALL add_var(ocean_default_list, 'T2', p_ice%T2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T2', 'C', 'Temperature lower layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%T2)

    CALL add_var(ocean_default_list, 'E1', p_ice%E1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E1', 'Jm/kg', 'Energy content upper layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%E1)

    CALL add_var(ocean_default_list, 'E2', p_ice%E2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E2', 'Jm/kg', 'Energy content lower layer', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(p_ice%E2)

  ENDIF

    ! dynamics

    CALL add_var(ocean_default_list, 'ice_u', p_ice%u ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_u', 'm s-1', 'zonal velocity', datatype_flt),&
      &          grib2_var(10, 2, 4, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_default"), lopenacc=.TRUE., initval=0.0_wp)
    __acc_attach(p_ice%u)

    CALL add_var(ocean_default_list, 'ice_v', p_ice%v ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_v', 'm s-1', 'meridional velocity', datatype_flt),&
      &          grib2_var(10, 2, 5, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_default"), lopenacc=.TRUE., initval=0.0_wp)
    __acc_attach(p_ice%v)

    IF ( i_ice_dyn == 1 ) THEN

      CALL add_var(ocean_restart_list, 'ice_u_prog', p_ice%u_prog ,&
           &          GRID_UNSTRUCTURED_VERT, ZA_SURFACE, &
           &          t_cf_var('ice_u_prog', 'm s-1', 'zonal velocity', datatype_flt),&
           &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_VERTEX),&
           &          ldims=(/nproma,nblks_v/),lrestart_cont=.TRUE., lopenacc=.TRUE.)
      __acc_attach(p_ice%u_prog)

      CALL add_var(ocean_restart_list, 'ice_v_prog', p_ice%v_prog ,&
           &          GRID_UNSTRUCTURED_VERT, ZA_SURFACE, &
           &          t_cf_var('ice_v', 'm s-1', 'meridional velocity', datatype_flt),&
           &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_VERTEX),&
           &          ldims=(/nproma,nblks_v/),lrestart_cont=.TRUE., lopenacc=.TRUE.)
      __acc_attach(p_ice%v_prog)

      CALL add_var(ocean_default_list, 'ice_vn', p_ice%vn_e ,&
           &          GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
           &          t_cf_var('ice_vn', 'm s-1', 'zonal velocity', datatype_flt),&
           &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_EDGE),&
           &          ldims=(/nproma,nblks_e/), lopenacc=.TRUE.)
      __acc_attach(p_ice%vn_e)

    ELSE IF ( i_ice_dyn == 2 ) THEN

      CALL add_var(ocean_restart_list, 'ice_vn', p_ice%vn_e ,&
        &          GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
        &          t_cf_var('ice_vn', 'm s-1', 'zonal velocity', datatype_flt),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_EDGE),&
        &          ldims=(/nproma,nblks_e/), lopenacc=.TRUE.)
      __acc_attach(p_ice%vn_e)

      CALL add_var(ocean_restart_list, 'ice_vt', p_ice%vt_e ,&
        &          GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
        &          t_cf_var('ice_vt', 'm/s', 'zonal velocity', datatype_flt),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_EDGE),&
        &          ldims=(/nproma,nblks_e/),lrestart_cont=.TRUE.)

      CALL add_var(ocean_default_list, 'Delta', p_ice%Delta ,&
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('Delta', '1/s', 'Defomation', datatype_flt),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_default"),lrestart_cont=.TRUE.)

    ENDIF

    ! not currently used categorywise limiter
    ALLOCATE(p_ice%hi_lim(i_no_ice_thick_class), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for hi_lim failed')
    END IF

    IF(p_ice%kice==1)THEN
      p_ice%hi_lim = 0.0_wp
    ELSEIF(p_ice%kice==8)THEN
      p_ice%hi_lim(:)=(/ 0.0_wp, 0.1_wp, 0.3_wp, 0.7_wp, 1.1_wp, 1.5_wp, 2.0_wp, 2.5_wp /)
    ENDIF

    IF ( i_ice_dyn > 0 ) THEN ! AWI dynamics
      CALL init_fem_wgts(patch_3D)
      IF (i_ice_advec == 1) THEN ! AWI advection
        CALL init_fem_wgts_extra(p_patch)
      ENDIF
    ENDIF

    ! add accumulated fields
    CALL message(TRIM(routine), 'end' )

    CALL construct_sea_ice_budgets(patch_3D,p_ice%budgets, ocean_default_list)
  END SUBROUTINE construct_sea_ice


  SUBROUTINE construct_sea_ice_budgets(patch_3d,budgets, varlist)
    TYPE(t_patch_3d) :: patch_3d
    TYPE(t_sea_ice_budgets) :: budgets
    TYPE(t_var_list_ptr)        :: varlist

    INTEGER :: alloc_cell_blocks
    TYPE(t_patch), POINTER :: patch
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    patch             => patch_3D%p_patch_2D(1)
    alloc_cell_blocks =  patch%alloc_cell_blocks

    CALL add_var(varlist, 'salt_00', budgets%salt_00 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('salt_00', 'kg', '', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_budgets"))
  END SUBROUTINE construct_sea_ice_budgets

  !-------------------------------------------------------------------------
  !
  !> Destructor of sea-ice model, deallocates all components.
  !!
  SUBROUTINE destruct_sea_ice(p_ice)
    TYPE (t_sea_ice),  INTENT (INOUT) :: p_ice
    !Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_init_thermo:destruct_sea_ice'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    DEALLOCATE(p_ice%hi_lim, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for hi_lim failed')
    END IF

    IF ( i_ice_dyn == 1 ) THEN ! AWI dynamics
      CALL destruct_fem_wgts
    ENDIF

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_sea_ice

  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmos fluxes for hydrostatic ocean
  !!
  SUBROUTINE construct_atmos_fluxes(p_patch, atmos_fluxes, i_no_ice_thick_class)
    !
    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_atmos_fluxes ), INTENT(INOUT) :: atmos_fluxes
    INTEGER,               INTENT(IN)    :: i_no_ice_thick_class
    ! Local variables
    INTEGER :: alloc_cell_blocks

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_init_thermo:construct_atmos_fluxes'
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    ! TODO: cleanup unused fluxes, sort in ice_diag
    ! Variables with 3 dimensions: fluxes over ice-covered part, second dimension is ice class

    !$ACC ENTER DATA COPYIN(atmos_fluxes)

    CALL add_var(ocean_default_list, 'atmos_fluxes_lat', atmos_fluxes%lat,                            &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_lat', 'W m-2', 'atmos_fluxes_lat', datatype_flt),            &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%lat)

    CALL add_var(ocean_default_list, 'atmos_fluxes_sens', atmos_fluxes%sens,                          &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_sens', 'W m-12', 'atmos_fluxes_sens', datatype_flt),          &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%sens)

    CALL add_var(ocean_default_list, 'atmos_fluxes_LWnet', atmos_fluxes%LWnet,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_LWnet', 'W m-2', 'atmos_fluxes_LWnet', datatype_flt),        &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%LWnet)

    CALL add_var(ocean_default_list, 'atmos_fluxes_SWnet', atmos_fluxes%SWnet,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_SWnet', 'W m-12', 'atmos_fluxes_SWnet', datatype_flt),        &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%SWnet)

    CALL add_var(ocean_default_list, 'atmos_fluxes_dsensdT', atmos_fluxes%dsensdT,                    &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_dsensdT', 'W m-2 K-1', 'atmos_fluxes_dsensdT', datatype_flt),  &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%dsensdT)

    CALL add_var(ocean_default_list, 'atmos_fluxes_dlatdT', atmos_fluxes%dlatdT,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_dlatdT', 'W m-2 K-1', 'atmos_fluxes_dlatdT', datatype_flt),    &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%dlatdT)

    CALL add_var(ocean_default_list, 'atmos_fluxes_dLWdT', atmos_fluxes%dLWdT,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_dLWdT', 'W m-2 K-1', 'atmos_fluxes_dLWdT', datatype_flt),      &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%dLWdT)

    ! Variables with 2 dimensions: fluxes over icefree part in sea ice model

    CALL add_var(ocean_default_list, 'atmos_fluxes_latw', atmos_fluxes%latw,                          &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_latw', 'W m-2', 'atmos_fluxes_latw', datatype_flt),          &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%latw)

    CALL add_var(ocean_default_list, 'atmos_fluxes_sensw', atmos_fluxes%sensw,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_sensw', 'W m-12', 'atmos_fluxes_sensw', datatype_flt),        &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%sensw)

    CALL add_var(ocean_default_list, 'atmos_fluxes_LWnetw', atmos_fluxes%LWnetw,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_LWnetw', 'W m-2', 'atmos_fluxes_LWnetw', datatype_flt),      &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%LWnetw)

    CALL add_var(ocean_default_list, 'atmos_fluxes_SWnetw', atmos_fluxes%SWnetw,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_SWnetw', 'W m-2', 'atmos_fluxes_SWnetw', datatype_flt),      &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%SWnetw)

    CALL add_var(ocean_default_list, 'atmos_fluxes_rprecw', atmos_fluxes%rprecw,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_rprecw', 'm s-1', 'atmos_fluxes_rprecw', datatype_flt),       &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%rprecw)

    CALL add_var(ocean_default_list, 'atmos_fluxes_rpreci', atmos_fluxes%rpreci,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_rpreci', 'm s-1', 'atmos_fluxes_rpreci', datatype_flt),       &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%rpreci)

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_x', atmos_fluxes%stress_x,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('surface_downward_eastward_stress', 'Pa',   'atmos_fluxes_stress_x', datatype_flt),  &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
      &          lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_y', atmos_fluxes%stress_y,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('surface_downward_northward_stress', 'Pa',   'atmos_fluxes_stress_y', datatype_flt),  &
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
      &          lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_xw', atmos_fluxes%stress_xw,                &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('surface_downward_eastward_stress', 'Pa',   'atmos_fluxes_stress_xw', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
      &          lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_yw', atmos_fluxes%stress_yw,                &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('surface_downward_northward_stress', 'Pa',   'atmos_fluxes_stress_yw', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
      &          lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    IF (.NOT.is_coupled_to_atmo()) THEN
      __acc_attach(atmos_fluxes%stress_x)
      __acc_attach(atmos_fluxes%stress_y)
      __acc_attach(atmos_fluxes%stress_xw)
      __acc_attach(atmos_fluxes%stress_yw)
    END IF  !  coupled

!fixme check if this can be avoided
    CALL add_var(ocean_default_list, 'albvisdirw', atmos_fluxes%albvisdirw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albvisdirw', '1', 'albvisdirw', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albvisdirw)

    CALL add_var(ocean_default_list, 'albvisdifw', atmos_fluxes%albvisdifw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albvisdifw', '1', 'albvisdifw', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albvisdifw)

    CALL add_var(ocean_default_list, 'albnirdirw', atmos_fluxes%albnirdirw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albnirdirw', '1', 'albnirdirw', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albnirdirw)

    CALL add_var(ocean_default_list, 'albnirdifw', atmos_fluxes%albnirdifw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albnirdifw', '1', 'albnirdifw', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albnirdifw)

    CALL add_var(ocean_restart_list, 'albvisdir', atmos_fluxes%albvisdir ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albvisdir', '1', 'albvisdir', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albvisdir)

    CALL add_var(ocean_restart_list, 'albvisdif', atmos_fluxes%albvisdif ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albvisdif', '1', 'albvisdif', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albvisdif)

    CALL add_var(ocean_restart_list, 'albnirdir', atmos_fluxes%albnirdir ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albnirdir', '1', 'albnirdir', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albnirdir)

    CALL add_var(ocean_restart_list, 'albnirdif', atmos_fluxes%albnirdif ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albnirdif', '1', 'albnirdif', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"), lopenacc=.TRUE.)
    __acc_attach(atmos_fluxes%albnirdif)

    ! Initialize with zero
    atmos_fluxes%counter = 0

    ! Initialise the albedos sensibly
    atmos_fluxes%albvisdir (:,:,:) = albi
    atmos_fluxes%albvisdif (:,:,:) = albi
    atmos_fluxes%albnirdir (:,:,:) = albi
    atmos_fluxes%albnirdif (:,:,:) = albi
    atmos_fluxes%albvisdirw(:,:) = albedoW_sim
    atmos_fluxes%albvisdifw(:,:) = albedoW_sim
    atmos_fluxes%albnirdirw(:,:) = albedoW_sim
    atmos_fluxes%albnirdifw(:,:) = albedoW_sim

    ! forcing of zonal component of velocity equation,
    CALL add_var(ocean_default_list,'atmos_fluxes_topBC_windStress_u', atmos_fluxes%topBoundCond_windStress_u,     &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_topBC_windStress_u', 'Pa', 'atmos_fluxes_topBoundCond_windStress_u', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! forcing of meridional component of velocity equation,
    CALL add_var(ocean_default_list,'atmos_fluxes_topBC_windStress_v', atmos_fluxes%topBoundCond_windStress_v,     &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_topBC_windStress_v', 'Pa', 'atmos_fluxes_topBoundCond_windStress_v', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! surface short wave heat flux                              [W/m2]
    CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_ShortWave', atmos_fluxes%HeatFlux_ShortWave,         &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_HeatFlux_ShortWave', 'W m-2', 'atmos_fluxes_HeatFlux_ShortWave', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! surface long wave heat flux                               [W/m2]
    CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_LongWave', atmos_fluxes%HeatFlux_LongWave,           &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_HeatFlux_LongWave', 'W m-2', 'atmos_fluxes_HeatFlux_LongWave', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! surface sensible heat flux                                [W/m2]
    CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_Sensible', atmos_fluxes%HeatFlux_Sensible,           &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_HeatFlux_Sensible', 'W m-2', 'atmos_fluxes_HeatFlux_Sensible', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! surface latent heat flux                                  [W/m2]
    CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_Latent', atmos_fluxes%HeatFlux_Latent,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_HeatFlux_Latent', 'W m-2', 'atmos_fluxes_HeatFlux_Latent', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! total heat flux                                  [W/m2]
    CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_Total', atmos_fluxes%HeatFlux_Total,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_HeatFlux_Total', 'W m-2', 'atmos_fluxes_HeatFlux_Total', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! total precipitation flux                                  [m/s]
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_Precipitation', atmos_fluxes%FrshFlux_Precipitation, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_Precipitation', 'm s-1', 'atmos_fluxes_FrshFlux_Precipitation', &
      &          datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! total snow flux                                           [m/s]
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_SnowFall', atmos_fluxes%FrshFlux_SnowFall,           &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_SnowFall', 'm s-1', 'atmos_fluxes_FrshFlux_SnowFall', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! evaporation flux                                          [m/s]
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_Evaporation', atmos_fluxes%FrshFlux_Evaporation,     &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_Evaporation', 'm s-1', 'atmos_fluxes_FrshFlux_Evaporation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),     &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
     
    ! river runoff flux                                         [m/s]
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_Runoff', atmos_fluxes%FrshFlux_Runoff, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_Runoff', 'm s-1', 'atmos_fluxes_FrshFlux_Runoff', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_TotalSalt', atmos_fluxes%FrshFlux_TotalSalt, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_TotalSalt', 'm s-1]', 'atmos_fluxes_FrshFlux_TotalSalt', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_TotalOcean', atmos_fluxes%FrshFlux_TotalOcean, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_TotalOcean', 'm s-1', 'atmos_fluxes_FrshFlux_TotalOcean', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_TotalIce', atmos_fluxes%FrshFlux_TotalIce, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_TotalIce', 'm s-1', 'atmos_fluxes_FrshFlux_TotalIce', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_VolumeIce', atmos_fluxes%FrshFlux_VolumeIce, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_VolumeIce', 'm s-1', 'atmos_fluxes_FrshFlux_VolumeIce', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    !  atmos_fluxes%FrshFlux_VolumeTotal is zero due to disassociated pointer in ocean_bulk
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_VolumeTotal', atmos_fluxes%FrshFlux_VolumeTotal, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_VolumeTotal', 'm s-1', 'atmos_fluxes_FrshFlux_VolumeTotal', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
!   CALL add_var(ocean_default_list,'atmos_flux_cellThicknessUnderIce', atmos_fluxes%cellThicknessUnderIce, &
!     &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
!     &          t_cf_var('atmos_flux_cellThicknessUnderIce', 'm', 'atmos_flux_cellThicknessUnderIce', datatype_flt),&
!     &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
!     &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL message(TRIM(routine), 'end' )
  END SUBROUTINE construct_atmos_fluxes


END MODULE mo_ice_init_thermo
