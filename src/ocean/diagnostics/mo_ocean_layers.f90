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

! Contains a diagnostic that estimates the transport between isopycnal layers
! and that also estimates diapycnal velocities.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_layers
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: &
    & vert_cor_type,                                          &
    & n_zlev, bottom_drag_coeff,                              &
    & n_dlev, rho_lev_in, mode_layers,                        & ! parameter from namelist
    & OceanReferenceDensity,                                  &
    & ReferencePressureIndbars,                               &
    & EOS_TYPE, &
    & LinearThermoExpansionCoefficient, &
    & LinearHalineContractionCoefficient, &
    & use_lbound_dirichlet!,        &

  USE mo_ocean_physics_types, ONLY: t_ho_params, v_params
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_impl_constants,      ONLY: success, max_char_length, min_dolic, sea
  USE mo_var_list,            ONLY: add_var
  USE mo_grib2,               ONLY: grib2_var, t_grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32 => CDI_DATATYPE_FLT32, &
    &                               DATATYPE_FLT64 => CDI_DATATYPE_FLT64, &
    &                               DATATYPE_INT8 => CDI_DATATYPE_INT8, &
    &                               DATATYPE_PACK16 => CDI_DATATYPE_PACK16, &
    &                               tstep_constant, GRID_LONLAT, GRID_UNSTRUCTURED
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell, grid_unstructured_edge, &
    &                               grid_unstructured_vert, grid_vertex, GRID_ZONAL
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_var_groups,          ONLY: groups, max_groups
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_onEdges_Pointer_3d_wp, t_onCells_HalfLevels_Pointer_wp, t_operator_coeff, t_hydro_ocean_diag
  USE mo_ocean_state,         ONLY: oce_config, ocean_default_list
  USE mo_physical_constants,  ONLY: grav
  USE mo_cf_convention
  USE mo_zaxis_type,          ONLY: &
    & za_depth_below_sea, za_depth_below_sea_half, za_surface, &
    & za_oce_layer_interface, za_oce_layer_centre
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max, sync_patch_array_mult
  USE mo_ocean_thermodyn,     ONLY: calculate_density_onColumn
  !USE mo_ocean_math_operators,ONLY: div_oce_3d
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, &
    & timer_extra10, timer_extra11
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_physical_constants,  ONLY: grav, clw
  USE mo_ocean_surface_types, ONLY: t_ocean_surface
  USE mo_ocean_thermodyn,     ONLY: calc_neutralslope_coeff_func_onColumn

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_layers
  PUBLIC :: calc_layers

  ! --- namelist parameters 
  !(if paras should not be read hear and not in mo_ocean_namelist; does not work)
  !INTEGER :: n_dlev = 5          ! number of density layers
  !INTEGER :: mode_layers = 1                     ! mode to derive layer transport
  !REAL(wp) :: rho_lev_in(1024)   ! density levels rho_lev_in(n_dlev+1)

  ! --- other global variables
  REAL(wp), DIMENSION(:), ALLOCATABLE :: rho_lev
  REAL(wp), DIMENSION(:), ALLOCATABLE :: rho_lev_cent
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: layer_thickness_c_old
  INTEGER :: tstep_count


CONTAINS

!To do (optional):
!   * check whether namelist reading can be done here instead of in main namelist
!   * rename rho_lev to layer_int and rho_lev_cent to layer_cent 
!   * rename sigma_2 to layer_scalar
!   * make switch to allow other level-types like temp, salt or sigma0

  SUBROUTINE init_layers(patch_3d, ocean_state_diag)
    TYPE(t_patch_3D), TARGET, INTENT(in)      :: patch_3d
    TYPE(t_hydro_ocean_diag), INTENT(inout)   :: ocean_state_diag

    INTEGER :: jd
    CHARACTER(LEN = *), PARAMETER :: routine = 'mo_ocean_layers:init_layers'
    INTEGER :: status,i, alloc_cell_blocks, nblks_e, nblks_v
    INTEGER :: datatype_flt

    !INTEGER :: iunit, i_status
    !NAMELIST/ocean_layers_nml/ n_dlev, rho_lev_in

    !CALL position_nml ('ocean_layers_nml', status=i_status)
    !IF (my_process_is_stdio()) THEN
    !  iunit = temp_defaults()
    !  WRITE(iunit, ocean_layers_nml)    ! write defaults to temporary text file
    !END IF
    !SELECT CASE (i_status)
    !CASE (positioned)
    !  READ (nnml, ocean_layers_nml)                            ! overwrite default settings
    !  IF (my_process_is_stdio()) THEN
    !    iunit = temp_settings()
    !    WRITE(iunit, ocean_layers_nml)    ! write settings to temporary text file
    !  END IF
    !END SELECT

    ALLOCATE(rho_lev(1:n_dlev+1))
    ALLOCATE(rho_lev_cent(1:n_dlev))

    ! --- level interfaces taken from namelist
    rho_lev(:) = -10.0_wp
    rho_lev(1:n_dlev+1) = rho_lev_in(1:n_dlev+1)
    !rho_lev(:) = rho_lev(:) !+ 1000.0_wp
    
    ! --- if n_dlev was set to a too large value in namelist there will appear negative values which is bad for the algorithm
    DO jd = 1,n_dlev+1
      IF (rho_lev(jd)<0.) THEN
        rho_lev(jd)=50.
      ENDIF
    ENDDO

    rho_lev_cent(:) = 0.5*(rho_lev(1:n_dlev)+rho_lev(2:n_dlev+1))

    ! Example layer definition
    !rho_lev(1:38) = (/ &
    !   & 20.  , 21.  , 28.1 , 28.9 , 29.7 , 30.5 , 30.95, 31.5 , 32.05, &
    !   & 32.6 , 33.15, 33.7 , 34.25, 34.75, 35.15, 35.5 , 35.8 , 36.04, &
    !   & 36.2 , 36.38, 36.52, 36.62, 36.7 , 36.77, 36.83, 36.89, 36.97, &
    !   & 37.02, 37.06, 37.09, 37.11, 37.13, 37.15, 37.2 , 37.3 , 37.4 , &
    !   & 39.  , 40. &
    !             & /)

    ! --- some information to print out as log
    CALL message (TRIM(routine), "Layers package setup:")
    WRITE(message_text,'(a,i4)') "mode_layers = ", mode_layers
    CALL message (TRIM(routine), message_text)
    WRITE(message_text,'(a,i4)') "n_dlev = ", n_dlev
    CALL message (TRIM(routine), message_text)
    DO jd = 1, n_dlev+1
      WRITE(message_text,'(a,i4,a,f10.3)') "rho_lev(", jd, ") = ", rho_lev(jd)
      CALL message (TRIM(routine), message_text)
    END DO
    DO jd = 1, n_dlev
      WRITE(message_text,'(a,i4,a,f10.3)') "rho_lev_cent(", jd, ") = ", rho_lev_cent(jd)
      CALL message (TRIM(routine), message_text)
    END DO

    ! --- here starts diagnostics
    alloc_cell_blocks = patch_3d%p_patch_2d(1)%alloc_cell_blocks
    nblks_e = patch_3d%p_patch_2d(1)%nblks_e
    nblks_v = patch_3d%p_patch_2d(1)%nblks_v
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

    CALL add_var(ocean_default_list,'mass_flux_lay', ocean_state_diag%mass_flux_lay, &
      & grid_unstructured_edge,&
      & za_oce_layer_centre, &
      & t_cf_var('mass_flux_lay', 'm2 s-1', 'mass flux in isopycnal layer', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_dlev,nblks_e/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
    !__acc_attach(ocean_state_diag%mass_flux_lay)

    CALL add_var(ocean_default_list,'layer_thickness_e', ocean_state_diag%layer_thickness_e, &
      & grid_unstructured_edge,&
      & za_oce_layer_centre, &
      & t_cf_var('layer_thickness_e', 'm', 'thickness of isopycnal layer on edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_dlev,nblks_e/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
    !__acc_attach(ocean_state_diag%layer_thickness_e)

    CALL add_var(ocean_default_list,'layer_thickness_c', ocean_state_diag%layer_thickness_c, &
      & grid_unstructured_cell,&
      & za_oce_layer_centre, &
      & t_cf_var('layer_thickness_c', 'm', 'thickness of isopycnal layer on cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_dlev,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
    !__acc_attach(ocean_state_diag%layer_thickness_c)

    IF (mode_layers == 1) THEN
      CALL add_var(ocean_default_list,'dhdt_tot', ocean_state_diag%dhdt_tot, &
        & grid_unstructured_cell,&
        & za_oce_layer_centre, &
        & t_cf_var('dhdt_tot', 'm/s', 'total layer thickness change', datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_dlev,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
      !__acc_attach(ocean_state_diag%dhdt_tot)

      CALL add_var(ocean_default_list,'dhdt_srf', ocean_state_diag%dhdt_srf, &
        & grid_unstructured_cell,&
        & za_oce_layer_centre, &
        & t_cf_var('dhdt_srf', 'm/s', 'layer thickness change by surface density flux', datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_dlev,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
      !__acc_attach(ocean_state_diag%dhdt_srf)

      CALL add_var(ocean_default_list,'dhdt_hfl', ocean_state_diag%dhdt_hfl, &
        & grid_unstructured_cell,&
        & za_oce_layer_centre, &
        & t_cf_var('dhdt_hfl', 'm/s', 'layer thickness change by horizontal flux', datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_dlev,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
      !__acc_attach(ocean_state_diag%dhdt_hfl)
    ENDIF

    CALL add_var(ocean_default_list,'div_mass_flux_lay', ocean_state_diag%div_mass_flux_lay, &
      & grid_unstructured_cell,&
      & za_oce_layer_centre, &
      & t_cf_var('div_mass_flux_lay', 'm/s', 'divergence of mass flux within layer', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_dlev,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
    !__acc_attach(ocean_state_diag%div_mass_flux_lay)

    CALL add_var(ocean_default_list, 'diapycnal_velocity', ocean_state_diag%diapycnal_velocity, &
      & grid_unstructured_cell,&
      & za_oce_layer_interface, &
      & t_cf_var('diapycnal_velocity','m s-1','diapycnal velocity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_dlev+1,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
    !__acc_attach(ocean_state_diag%diapycnal_velocity)

    CALL add_var(ocean_default_list, 'sigma2', ocean_state_diag%sigma2, &
      & grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('sigma2','kg/m^3','potential density ref. to 2000m', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
    !__acc_attach(ocean_state_diag%sigma2)

    IF (mode_layers == 1) THEN
      CALL add_var(ocean_default_list, 'sflx_dens', ocean_state_diag%sflx_dens, &
        & grid_unstructured_cell,&
        & za_surface, &
        & t_cf_var('sflx_dens','kg/s/m^2','potential density ref. to 2000m', datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
      !__acc_attach(ocean_state_diag%sflx_dens)

      ! drho = -rho0 alphaT dT + rho0 betaS dS

      CALL add_var(ocean_default_list, 'alphaT', ocean_state_diag%alphaT, &
        & grid_unstructured_cell,&
        & za_surface,&
        & t_cf_var('alphaT','1/K','thermal expansion coefficient', datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
      !__acc_attach(ocean_state_diag%alphaT)

      CALL add_var(ocean_default_list, 'betaS', ocean_state_diag%betaS, &
        & grid_unstructured_cell,&
        & za_surface,&
        & t_cf_var('betaS','m^3/kg','haline expansion coefficient', datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
      !__acc_attach(ocean_state_diag%betaS)

      CALL add_var(ocean_default_list, 'weight_e_sum', ocean_state_diag%weight_e_sum, &
        & grid_unstructured_edge,&
        & za_depth_below_sea,&
        & t_cf_var('weight_e_sum','','weight_e_sum', datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_layers"), lopenacc=.TRUE.)
      !__acc_attach(ocean_state_diag%weight_e_sum)
    ENDIF

    ocean_state_diag%layer_thickness_c = 0.0_wp
    ALLOCATE(layer_thickness_c_old(nproma,1:n_dlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks))
    layer_thickness_c_old = 0.0_wp
    tstep_count = 0

  END SUBROUTINE init_layers

  SUBROUTINE calc_layers(patch_3d, ocean_state, p_oce_sfc, op_coeffs, params_oce, stretch_c, stretch_e)
    TYPE(t_patch_3d ), TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in) :: ocean_state
    TYPE(t_operator_coeff), INTENT(in)        :: op_coeffs
    TYPE(t_ho_params), TARGET, INTENT(inout)  :: params_oce
    REAL(wp), INTENT(IN), OPTIONAL :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN), OPTIONAL :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)

    TYPE(t_patch), POINTER                    :: patch_2d
    TYPE(t_subset_range), POINTER             :: edges_in_domain, cells_in_domain, all_edges, all_cells
    TYPE (t_ocean_surface), INTENT(IN)        :: p_oce_sfc
    INTEGER,  DIMENSION(:,:,:),POINTER        :: idx_c, blk_c

    REAL(wp), DIMENSION(:,:,:), POINTER       :: dhdt_tot
    REAL(wp), DIMENSION(:,:,:), POINTER       :: dhdt_srf
    REAL(wp), DIMENSION(:,:,:), POINTER       :: dhdt_hfl
    REAL(wp), DIMENSION(:,:,:), POINTER       :: div_mass_flux_lay
    REAL(wp), DIMENSION(:,:,:), POINTER       :: sigma2
    REAL(wp), DIMENSION(:,:,:), POINTER       :: mass_flux_lay
    REAL(wp), DIMENSION(:,:,:), POINTER       :: layer_thickness_e
    REAL(wp), DIMENSION(:,:,:), POINTER       :: layer_thickness_c
    REAL(wp), DIMENSION(:,:,:), POINTER       :: diapycnal_velocity
    REAL(wp), DIMENSION(:,:,:), POINTER       :: weight_e_sum
    REAL(wp), DIMENSION(:,:), POINTER         :: sflx_dens
    REAL(wp), DIMENSION(:,:), POINTER         :: alphaT 
    REAL(wp), DIMENSION(:,:), POINTER         :: betaS
    INTEGER, DIMENSION(:,:), POINTER          :: dolic_e
    INTEGER, DIMENSION(:,:), POINTER           :: dolic_c
    REAL(wp), DIMENSION(n_zlev)                :: rho_c
    REAL(wp), DIMENSION(n_zlev+1)              :: rho_ci
    REAL(wp), DIMENSION(n_zlev)                :: rho_e
    REAL(wp), DIMENSION(n_zlev+1)              :: rho_ei
    REAL(wp), DIMENSION(n_zlev)                :: dzc
    REAL(wp), DIMENSION(n_zlev)                :: dze
    REAL(wp), DIMENSION(n_zlev)               :: ones

    INTEGER,  DIMENSION(:,:,:), POINTER       :: idx, blk
    INTEGER :: level
    INTEGER :: levels

    INTEGER                                   :: jc, je, jk, jd, jds, jde, blockNo, tracer_index, kk
    INTEGER                                   :: start_index, end_index
    REAL(wp) :: neutral_coeff(1:n_zlev, 2)
    REAL(wp) :: sflx_temp, sflx_salt
    REAL(wp) :: dmin, dmax, ddif, weight
    REAL(wp) :: mfe, div
    
    LOGICAL, SAVE                             :: firstcall=.true.

    ! --- for debugging
    !INTEGER                                   :: bli, jei
    !REAL(wp)                                  :: tmp1, tmp2
    !write(*,*) "Layers diagnostic."

    ones = 1.0_wp

    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    all_cells       => patch_2d%cells%ALL
    all_edges       => patch_2d%edges%ALL

    idx_c => patch_3d%p_patch_2d(1)%edges%cell_idx
    blk_c => patch_3d%p_patch_2d(1)%edges%cell_blk
    dolic_e            => patch_3d%p_patch_1d(1)%dolic_e
    dolic_c            => patch_3d%p_patch_1d(1)%dolic_c

    sigma2             => ocean_state%p_diag%sigma2
    mass_flux_lay      => ocean_state%p_diag%mass_flux_lay
    layer_thickness_e  => ocean_state%p_diag%layer_thickness_e
    layer_thickness_c  => ocean_state%p_diag%layer_thickness_c
    div_mass_flux_lay  => ocean_state%p_diag%div_mass_flux_lay
    diapycnal_velocity => ocean_state%p_diag%diapycnal_velocity

    sigma2 = 0.0_wp
    mass_flux_lay = 0.0_wp
    layer_thickness_e = 0.0_wp
    layer_thickness_c = 0.0_wp
    div_mass_flux_lay = 0.0_wp
    diapycnal_velocity = 0.0_wp

    ! --- calculate sigma_2000 (potential density referenced to 2000m)
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        ! FIXME: nnew or nold?
        ocean_state%p_diag%sigma2(jc,:,blockNo)  = (calculate_density_onColumn(&
            & ocean_state%p_prog(nnew(1))%tracer(jc,:,blockNo,1), & ! temp
            & ocean_state%p_prog(nnew(1))%tracer(jc,:,blockNo,2), & ! salt
            & 2000.0_wp*ones, n_zlev) - 1000._wp) * patch_3d%wet_c(jc,:,blockNo)
      ENDDO
    ENDDO

    IF (mode_layers == 1) THEN
      sflx_dens          => ocean_state%p_diag%sflx_dens
      alphaT             => ocean_state%p_diag%alphaT
      betaS              => ocean_state%p_diag%betaS
      dhdt_tot           => ocean_state%p_diag%dhdt_tot
      dhdt_srf           => ocean_state%p_diag%dhdt_srf
      dhdt_hfl           => ocean_state%p_diag%dhdt_hfl
      weight_e_sum       => ocean_state%p_diag%weight_e_sum

      sflx_dens = 0.0_wp
      alphaT = 0.0_wp
      betaS = 0.0_wp
      dhdt_tot = 0.0_wp
      dhdt_srf = 0.0_wp
      dhdt_hfl = 0.0_wp
      weight_e_sum = 0.0_wp

      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        DO je = start_index, end_index
          ! FIXME: Check whether ==sea also includes wet point at boundary
          IF (patch_3d%lsm_e(je,1,blockNo)==sea) THEN
            IF ( vert_cor_type == 0 ) THEN
              dze = patch_3d%p_patch_1d(1)%prism_thick_e(je,:,blockNo)
            ELSEIF ( vert_cor_type == 1 ) THEN
              dze = patch_3d%p_patch_1d(1)%prism_thick_e(je,:,blockNo)*stretch_e(je,blockNo)
            ENDIF

            ! FIXME: How best treating case when one of the neighbours is land?
            ! FIXME: *\ interpolate vertically to get density on interfaces
            !        *\ extrapolate boundary values
            !        *\ pay attention when extrapolating land values; use kbot

            ! interpolate horizontally from center to edges
            rho_e = 0.0_wp
            DO jk = 1, dolic_e(je,blockNo)
              rho_e(jk) = patch_3d%wet_e(je,jk,blockNo) * 0.5_wp*( &
                &   ocean_state%p_diag%sigma2(idx_c(je,blockNo,1),jk,blk_c(je,blockNo,1)) &
                & + ocean_state%p_diag%sigma2(idx_c(je,blockNo,2),jk,blk_c(je,blockNo,2)) )
            ENDDO

            ! interpolate vertically from center to interfaces
            rho_ei = 0.0_wp
            DO jk = 2, dolic_e(je,blockNo)
              rho_ei(jk) = (rho_e(jk-1)*dze(jk) + rho_e(jk)*dze(jk-1)) &
                & / (dze(jk)+dze(jk-1))
            ENDDO
            jk = 1
            rho_ei(jk) = (rho_ei(jk+1)-rho_ei(jk+2))*dze(jk)/dze(jk+1) + rho_ei(jk+1)
            !dens(nzmin)    =dens(nzmin+1)      +(dens(nzmin+1)-dens(nzmin+2))            *helem(nzmin, elem)/helem(nzmin+1,elem)
            !!if (blockNo==100 .and. je==2) then
            !!  write(*,*) 'jk = ', jk
            !!  write(*,*) 'rho_ei: jk, jk+1, jk+2', rho_ei(jk), rho_ei(jk+1), rho_ei(jk+2)
            !!  write(*,*) 'dze: jk, jk+1', dze(jk), dze(jk+1)
            !!endif
            jk = dolic_e(je,blockNo)+1
            rho_ei(jk) = (rho_ei(jk-1)-rho_ei(jk-2))*dze(jk-1)/dze(jk-2) + rho_ei(jk-1)

            ! main loop
            DO jk = 1, dolic_e(je,blockNo)
            
              dmin = minval(rho_ei(jk:jk+1))
              dmax = maxval(rho_ei(jk:jk+1))
              ddif = abs(rho_ei(jk)-rho_ei(jk+1))

              jds = n_dlev ! if all levels are above dmin, we need jds>jde
              DO jd = 2, n_dlev
                IF (rho_lev(jd) > dmin) THEN
                  jds = jd
                  EXIT
                ENDIF 
              ENDDO 

              jde = 1 ! if all levels are below dmax, we need jds>jde
              DO jd = n_dlev-1,1,-1
                IF (rho_lev(jd) < dmax) THEN
                  jde = jd
                  EXIT
                ENDIF
              ENDDO 

              !if (jds>jde) jds=jde    ! whole cell is within two density levels
              !if (jds==jde) jde=jds+1 ! one density level is within cell
              !if (jds<jde) do nothing -> normal case
              !!if (rho_lev(jds)>=dmax) jds=jde
              !!if (rho_lev(jde)<=dmin) jde=jds ! necessary?
               
              mfe = ocean_state%p_diag%mass_flx_e(je,jk,blockNo)
              IF (jds<=jde) THEN
                weight = abs(rho_lev(jds)-dmin)/ddif 
                mass_flux_lay(je,jds-1,blockNo) = mass_flux_lay(je,jds-1,blockNo) + mfe*weight
                layer_thickness_e(je,jds-1,blockNo) = layer_thickness_e(je,jds-1,blockNo) + dze(jk)*weight
                !if (blockNo==100 .and. je==2 .and. jk==24) then
                !  write(*,*) 'weight = ', weight
                !endif
                weight_e_sum(je,jk,blockNo) = weight_e_sum(je,jk,blockNo)+weight
                DO jd = jds,jde-1 ! more than one density level is within cell
                  weight = abs(rho_lev(jd+1)-rho_lev(jd))/ddif
                  mass_flux_lay(je,jd,blockNo) = mass_flux_lay(je,jd,blockNo) + mfe*weight
                  layer_thickness_e(je,jd,blockNo) = layer_thickness_e(je,jd,blockNo) + dze(jk)*weight
                !if (blockNo==100 .and. je==2 .and. jk==24) then
                !  write(*,*) 'weight = ', weight
                !endif
                  weight_e_sum(je,jk,blockNo) = weight_e_sum(je,jk,blockNo)+weight
                ENDDO
                weight = abs(dmax-rho_lev(jde))/ddif 
                mass_flux_lay(je,jde,blockNo) = mass_flux_lay(je,jde,blockNo) + mfe*weight
                layer_thickness_e(je,jde,blockNo) = layer_thickness_e(je,jde,blockNo) + dze(jk)*weight
                !if (blockNo==100 .and. je==2 .and. jk==24) then
                !  write(*,*) 'weight = ', weight
                !endif
                weight_e_sum(je,jk,blockNo) = weight_e_sum(je,jk,blockNo)+weight
              ELSE ! jds>jde whole cell is within two density levels
                ! use jde since it is the upper index in this case
                mass_flux_lay(je,jde,blockNo) = mass_flux_lay(je,jde,blockNo) + mfe
                layer_thickness_e(je,jde,blockNo) = layer_thickness_e(je,jde,blockNo) + dze(jk)
                weight_e_sum(je,jk,blockNo) = 1.0_wp
                !!if (blockNo==100 .and. je==2) then
                !!  write(*,*) 'Dumped all at once!'
                !!endif
              ENDIF

              !! - debugging
              !!if (blockNo==100 .and. je==2) then
              !!  write(*,*) 'jk = ', jk
              !!  write(*,*) 'mass_flx_e = ', mfe
              !!!  write(*,*) 'blockNo = ', blockNo, 'je = ', je, 'jk = ', jk
              !!!  !write(*,*) 'rho_lev = ', rho_lev
              !!!  write(*,*) 'dmin = ', dmin, 'dmax = ', dmax, 'ddif = ', ddif
              !!  write(*,*) 'jds = ', jds, 'jde = ', jde
              !!  write(*,*) 'rho_ei(jk) = ', rho_ei(jk), ' rho_ei(jk+1) = ', rho_ei(jk+1)
              !!  write(*,*) 'mass_flux_lay(jds-1) = ', mass_flux_lay(je,jds-1,blockNo)
              !!  write(*,*) 'mass_flux_lay(jds) = ', mass_flux_lay(je,jds,blockNo)
              !!  write(*,*) 'mass_flux_lay(jde) = ', mass_flux_lay(je,jde,blockNo)
              !!  write(*,*) 'weight_e_sum = ', weight_e_sum(je,jk,blockNo)
              !!!  write(*,*) 'weight = ', weight
              !!!  !write(*,*) 'rho_e(jk) = ', rho_e(jk), ' rho_lev(kk) = ', rho_lev(kk), ' kk = ', kk, ' jk = ', jk
              !!!  !write(*,*) 'mass_flux_e(jk) = ', ocean_state%p_diag%mass_flx_e(je,jk,blockNo)
              !!  write(*,*) '--------------------------------------'
              !!end if
            ENDDO ! jk
          ENDIF ! patch_3d%lsm_e(je,1,blockNo)==sea)
        ENDDO ! je
      ENDDO ! blockNo

      ! --- 
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_index, end_index)
        idx => patch_2d%cells%edge_idx
        blk => patch_2d%cells%edge_blk
        DO jc = start_index, end_index
          IF (patch_3d%lsm_c(jc,1,blockNo)==sea) THEN
            dzc = patch_3d%p_patch_1d(1)%prism_thick_c(jc,:,blockNo)

            ! get density on cell
            DO jk = 1, dolic_c(jc,blockNo)
              rho_c(jk) = patch_3d%wet_c(jc,jk,blockNo) * ocean_state%p_diag%sigma2(jc,jk,blockNo)
            ENDDO

            ! interpolate vertically from center to interfaces
            DO jk = 2, dolic_c(jc,blockNo)
              rho_ci(jk) = (rho_c(jk-1)*dzc(jk) + rho_c(jk)*dzc(jk-1)) &
                & / (dzc(jk)+dzc(jk-1))
            ENDDO
            jk = 1
            rho_ci(jk) = (rho_ci(jk+1)-rho_ci(jk+2))*dzc(jk)/dzc(jk+1) + rho_ci(jk+1)
            jk = dolic_c(jc,blockNo)+1
            rho_ci(jk) = (rho_ci(jk-1)-rho_ci(jk-2))*dzc(jk-1)/dzc(jk-2) + rho_ci(jk-1)

            ! main loop
            DO jk = 1, dolic_c(jc,blockNo)
            
              dmin = minval(rho_ci(jk:jk+1))
              dmax = maxval(rho_ci(jk:jk+1))
              ddif = abs(rho_ci(jk)-rho_ci(jk+1))

              jds = n_dlev ! if all levels are above dmin, we need jds>jde
              DO jd = 2, n_dlev
                IF (rho_lev(jd) > dmin) THEN
                  jds = jd
                  EXIT
                ENDIF 
              ENDDO 

              jde = 1 ! if all levels are below dmax, we need jds>jde
              DO jd = n_dlev-1,1,-1
                IF (rho_lev(jd) < dmax) THEN
                  jde = jd
                  EXIT
                ENDIF
              ENDDO

              !if (jds>jde) jds=jde    ! whole cell is within two density levels
              !if (jds==jde) jds=jde-1 ! one density level is within cell
              !if (jds<jde) do nothing -> normal case
              !!if (rho_lev(jds)>=dmax) jds=jde
              !!if (rho_lev(jde)<=dmin) jde=jds ! necessary?
               
              div =  &
                & ocean_state%p_diag%mass_flx_e(idx(jc,blockNo,1),jk,blk(jc,blockNo,1)) * op_coeffs%div_coeff(jc,jk,blockNo,1) + &
                & ocean_state%p_diag%mass_flx_e(idx(jc,blockNo,2),jk,blk(jc,blockNo,2)) * op_coeffs%div_coeff(jc,jk,blockNo,2) + &
                & ocean_state%p_diag%mass_flx_e(idx(jc,blockNo,3),jk,blk(jc,blockNo,3)) * op_coeffs%div_coeff(jc,jk,blockNo,3)
              IF (jds<=jde) THEN
                weight = (rho_lev(jds)-dmin)/ddif 
                dhdt_hfl(jc,jds-1,blockNo) = dhdt_hfl(jc,jds-1,blockNo) + div*weight
                layer_thickness_c(jc,jds-1,blockNo) = layer_thickness_c(jc,jds-1,blockNo) + dzc(jk)*weight
                DO jd = jds,jde-1 ! more than one density level is within cell
                  weight = (rho_lev(jd+1)-rho_lev(jd))/ddif
                  dhdt_hfl(jc,jd,blockNo) = dhdt_hfl(jc,jd,blockNo) + div*weight
                  layer_thickness_c(jc,jd,blockNo) = layer_thickness_c(jc,jd,blockNo) + dzc(jk)*weight
                ENDDO
                weight = (dmax-rho_lev(jde))/ddif 
                dhdt_hfl(jc,jde,blockNo) = dhdt_hfl(jc,jde,blockNo) + div*weight
                layer_thickness_c(jc,jde,blockNo) = layer_thickness_c(jc,jde,blockNo) + dzc(jk)*weight
              ELSE ! jds>jde whole cell is within two density levels
                dhdt_hfl(jc,jds,blockNo) = dhdt_hfl(jc,jds,blockNo) + div
                layer_thickness_c(jc,jds,blockNo) = layer_thickness_c(jc,jds,blockNo) + dzc(jk)
              ENDIF
              ! ! - debugging
              ! if (blockNo==100 .and. jc==2 .and. jk==9) then
              !   write(*,*) 'blockNo = ', blockNo, 'jc = ', jc, 'jk = ', jk
              !   write(*,*) 'rho_lev = ', rho_lev
              !   write(*,*) 'dmin = ', dmin, 'dmax = ', dmax, 'ddif = ', ddif
              !   write(*,*) 'jds = ', jds, 'jde = ', jde
              !   !write(*,*) 'rho_e(jk) = ', rho_e(jk), ' rho_lev(kk) = ', rho_lev(kk), ' kk = ', kk, ' jk = ', jk
              !   !write(*,*) 'mass_flux_e(jk) = ', ocean_state%p_diag%mass_flx_e(je,jk,blockNo)
              !   write(*,*) '--------------------------------------'
              ! end if
            ENDDO ! jk
          ENDIF ! patch_3d%lsm_c(jc,1,blockNo)==sea
        ENDDO ! jc
      ENDDO ! blockNo

      ! --- derive thickness change of layer
      ! write(*,*) dtime
      ! write(*,*) 'layer_thickness_c_old = ', layer_thickness_c_old(2,:,100)
      ! write(*,*) 'layer_thickness_c = ', layer_thickness_c(2,:,100)
      IF (firstcall) THEN
        dhdt_tot(:,:,:) = 0.0_wp
      ELSE
        dhdt_tot(:,:,:) = (layer_thickness_c(:,:,:) - layer_thickness_c_old(:,:,:)) / dtime
      ENDIF
      firstcall = .FALSE.
      layer_thickness_c_old = layer_thickness_c
      !if (blockNo==100 .and. jc==2 .and. jk==9) then
      ! write(*,*) 'layer_thickness_c_old = ', layer_thickness_c_old(2,:,100)
      ! write(*,*) 'dhdt_tot = ', dhdt_tot(2,:,100)
      ! write(*,*) '------------------------'

      ! --- surface transformation
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          levels       = dolic_e(jc, blockNo)   ! 

          ! only for ocean cells
          IF (dolic_c(jc,blockNo)>0) THEN

          ! FIXME: Try to derive neutral_coeff only for the surface
          ! calculate alpha and beta
          IF(EOS_TYPE/=1)THEN
            !Nonlinear EOS, slope coefficients are calculated via the McDougall-method
            ! FIXME: use salinity instead of salinityColumn?
            neutral_coeff = calc_neutralslope_coeff_func_onColumn( &
              & ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,1), &
              & ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2), &
              !& patch_3D%p_patch_1d(1)%depth_cellinterface(jc,2:levels+1,blockNo), &
              & 2000.0_wp*ones(1:levels), &
              & levels) 
          ELSEIF(EOS_TYPE==1)THEN
            !Linear EOS: slope coefficients are equal to EOS-coefficients
            neutral_coeff(:,1) = LinearThermoExpansionCoefficient
            neutral_coeff(:,2) = LinearHalineContractionCoefficient
          ENDIF
          sflx_temp = p_oce_sfc%HeatFlux_Total(jc,blockNo) / clw 
          sflx_salt = p_oce_sfc%FrshFlux_TotalOcean(jc,blockNo) * OceanReferenceDensity

          sflx_dens(jc,blockNo) = ( &
            &   neutral_coeff(1,1) * sflx_temp &
            & - neutral_coeff(1,2) * sflx_salt )
          alphaT(jc,blockNo) = neutral_coeff(1,1)
          betaS(jc,blockNo)  = neutral_coeff(1,2)

          jd = minloc(abs(rho_lev_cent - ocean_state%p_diag%sigma2(jc,1,blockNo)),1)
          dhdt_srf(jc,jd,blockNo) = dhdt_srf(jc,jd,blockNo) + sflx_dens(jc,blockNo)/(rho_lev(jd+1)-rho_lev(jd))
          !DO jd = 1, n_dlev
          !  IF (ocean_state%p_diag%sigma2(jc,1,blockNo)>rho_lev(jd) & 
          !    &  .and. ocean_state%p_diag%sigma2(jc,1,blockNo)<=rho_lev(jd+1)) then
          !    dhdt_srf(jc,jd,blockNo) = dhdt_srf(jc,jd,blockNo) & 
          !      & + sflx_dens(jc,blockNo) &
          !      !& * cell_area(jc,blockNo) &
          !      & / (rho_lev(jd+1)-rho_lev(jd))
          !  ENDIF
          !ENDDO
          !if (blockNo==100 .and. jc==2) then
          !  write(*,*) 'jd = ', jd, 'rho_lev_cent = ', rho_lev_cent(jd)
          !  write(*,*) 'sigma2 = ', ocean_state%p_diag%sigma2(jc,1,blockNo)
          !  write(*,*) 'drho_lev = ', (rho_lev(jd+1)-rho_lev(jd))
          !  write(*,*) 'sflx_dens(jc,blockNo) = ', sflx_dens(jc,blockNo)
          !  write(*,*) '------------------------'
          !ENDIF

          ENDIF !(kbot>0)
        ENDDO
      ENDDO
    ELSE  !(mode_layers != 1)
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        DO je = start_index, end_index
          IF (patch_3d%lsm_e(je,1,blockNo)==sea) THEN
            IF ( vert_cor_type == 0 ) THEN
              dze = patch_3d%p_patch_1d(1)%prism_thick_e(je,:,blockNo)
            ELSEIF ( vert_cor_type == 1 ) THEN
              dze = patch_3d%p_patch_1d(1)%prism_thick_e(je,:,blockNo)*stretch_e(je,blockNo)
            ENDIF
            ! interpolate horizontally from center to edges
            rho_e = 0.0_wp
            DO jk = 1, dolic_e(je,blockNo)
              rho_e(jk) = patch_3d%wet_e(je,jk,blockNo) * 0.5_wp*( &
                &   ocean_state%p_diag%sigma2(idx_c(je,blockNo,1),jk,blk_c(je,blockNo,1)) &
                & + ocean_state%p_diag%sigma2(idx_c(je,blockNo,2),jk,blk_c(je,blockNo,2)) )
            ENDDO

            ! main loop
            DO jk = 1, dolic_e(je,blockNo)
              jds = 1 ! This implies: Add all the transport to the first layer, that cannot be attributed to any other layer."
              DO jd = 1, n_dlev
                IF (rho_e(jk)>rho_lev(jd) .AND. rho_e(jk)<=rho_lev(jd+1)) THEN
                  jds = jd
                ENDIF            
              ENDDO
!              jds = n_dlev ! if all levels are above dmin, we need jds>jde
!              DO jd = 2, n_dlev 
!                IF (rho_lev(jd) > dmin) THEN
!                  jds = jd
!                  EXIT
!                ENDIF 
!              ENDDO 
              ! --- layer depth on edges
              layer_thickness_e(je,jds,blockNo) = layer_thickness_e(je,jds,blockNo) &
                & + patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo)
              ! --- accumulate mass flux of this layer (derive mass flux as for wvel)
              mass_flux_lay(je,jds,blockNo) = mass_flux_lay(je,jds,blockNo) & 
                & + ocean_state%p_diag%mass_flx_e(je,jk,blockNo)
              ! debugging
              !if (blockNo==100 .and. je==2) then
              !  write(*,*) 'jk = ', jk
              !  write(*,*) 'mass_flx_e = ', ocean_state%p_diag%mass_flx_e(je,jk,blockNo)
              !  write(*,*) 'jds = ', jds
              !  write(*,*) 'mass_flux_lay = ', mass_flux_lay(je,jds,blockNo)
              !  write(*,*) '------------------------'
              !ENDIF
            ENDDO
          ENDIF ! (patch_3d%lsm_e(je,1,blockNo)==sea)
        ENDDO
      ENDDO

      ! ---
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_index, end_index)
        idx => patch_2d%cells%edge_idx
        blk => patch_2d%cells%edge_blk
        DO jc = start_index, end_index
          IF (patch_3d%lsm_c(jc,1,blockNo)==sea) THEN
            IF ( vert_cor_type == 0 ) THEN
              dzc = patch_3d%p_patch_1d(1)%prism_thick_c(jc,:,blockNo)
            ELSEIF ( vert_cor_type == 1 ) THEN
              dzc = patch_3d%p_patch_1d(1)%prism_thick_c(jc,:,blockNo)*stretch_c(jc,blockNo)
            ENDIF
            ! get density on cell
            DO jk = 1, dolic_c(jc,blockNo)
              rho_c(jk) = patch_3d%wet_c(jc,jk,blockNo) * ocean_state%p_diag%sigma2(jc,jk,blockNo)
            ENDDO

            ! main loop
            DO jk = 1, dolic_c(jc,blockNo)
              jds = 1 ! dump all transport which cannot be attributed to first layer
              DO jd = 1, n_dlev
                IF (rho_c(jk)>rho_lev(jd) .AND. rho_c(jk)<=rho_lev(jd+1)) THEN
                  jds = jd
                ENDIF            
              ENDDO
              ! --- layer depth on cells
              layer_thickness_c(jc,jds,blockNo) = layer_thickness_c(jc,jds,blockNo) &
                & + patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,blockNo)
            ENDDO
          ENDIF ! (patch_3d%lsm_e(je,1,blockNo)==sea)
        ENDDO
      ENDDO
    ENDIF !(mode_layers == x)

    ! --- divergence of isopycnal flow
    level = 1
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      idx => patch_2d%cells%edge_idx
      blk => patch_2d%cells%edge_blk
      DO jd = 1, n_dlev
        DO jc = start_index, end_index
          !write(*,*) 'blockNo = ', blockNo, ' jd = ', jd, ' jc = ', jc ! debugging
          div_mass_flux_lay(jc,jd,blockNo) =  &
            & mass_flux_lay(idx(jc,blockNo,1),jd,blk(jc,blockNo,1)) * op_coeffs%div_coeff(jc,level,blockNo,1) + &
            & mass_flux_lay(idx(jc,blockNo,2),jd,blk(jc,blockNo,2)) * op_coeffs%div_coeff(jc,level,blockNo,2) + &
            & mass_flux_lay(idx(jc,blockNo,3),jd,blk(jc,blockNo,3)) * op_coeffs%div_coeff(jc,level,blockNo,3)
        ENDDO
      ENDDO
    ENDDO

    ! --- diapycnal velocity
    ! (integrate from densest layer upward and assume no diap transport through
    !  lowest layer)
    diapycnal_velocity = 0.0_wp
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        DO jd = n_dlev, 1, -1
          diapycnal_velocity(jc,jd,blockNo) = &
            & diapycnal_velocity(jc,jd+1,blockNo) - div_mass_flux_lay(jc,jd,blockNo)
        ENDDO
      ENDDO
    ENDDO

    tstep_count = tstep_count + 1
  END SUBROUTINE calc_layers

END MODULE mo_ocean_layers
