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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_idemix
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: &
    & n_zlev, bottom_drag_coeff,                              &
    & HarmonicViscosity_reference, velocity_VerticalDiffusion_background,                 &
    & Temperature_VerticalDiffusion_background, Salinity_VerticalDiffusion_background, no_tracer,                       &
    & tracer_convection_MixingCoefficient,                                     &
    & BiharmonicViscosity_scaling, HarmonicViscosity_scaling, &
    & VelocityDiffusion_order,                                &
    & BiharmonicViscosity_reference,                          &
    & tracer_RichardsonCoeff, velocity_RichardsonCoeff,                    &
    & PPscheme_type,                                &
    & PPscheme_Constant_type,                       &
    !& PPscheme_ICON_type,                        &
    & PPscheme_ICON_Edge_type,                   &
    & PPscheme_ICON_Edge_vnPredict_type,         &
    & use_wind_mixing,                                        &
    & HorizontalViscosity_SmoothIterations,                   &
    & convection_InstabilityThreshold,                        &
    & RichardsonDiffusion_threshold,                          &
    & use_reduced_mixing_under_ice,                           &
    & k_tracer_dianeutral_parameter,                          &
    & k_tracer_isoneutral_parameter, k_tracer_GM_kappa_parameter,    &
    & GMRedi_configuration,GMRedi_combined,                   &
    & GM_only,Redi_only,                                      &
    & laplacian_form,                                         &
    & HorizontalViscosity_SpatialSmoothFactor,                &
    & VerticalViscosity_TimeWeight, OceanReferenceDensity,    &
    & tracer_TopWindMixing, WindMixingDecayDepth,             &
    & velocity_TopWindMixing, TracerHorizontalDiffusion_scaling, &
    &  Temperature_HorizontalDiffusion_Background,            &
    &  Temperature_HorizontalDiffusion_Reference,             &
    &  Salinity_HorizontalDiffusion_Background,               &
    &  Salinity_HorizontalDiffusion_Reference,                &
    &  HarmonicViscosity_background,                          &
    &  BiharmonicViscosity_background,                        &
    &  LeithHarmonicViscosity_background, LeithHarmonicViscosity_reference,    &
    &  LeithHarmonicViscosity_scaling,                                         &
    &  LeithBiharmonicViscosity_background, LeithBiharmonicViscosity_reference,&
    &  LeithBiharmonicViscosity_scaling,                       &
    &  LeithClosure_order,   LeithClosure_form, &
    &  TracerDiffusion_LeithWeight, &!Salinity_ConvectionRestrict, &
    !&  max_turbulenece_TracerDiffusion_amplification, &
    &  ReferencePressureIndbars,                              &
    &  tau_v,                       &
    &  tau_h,                       &
    &  gamma,                       &
    &  jstar,                       &
    &  mu0,                         &
    &  l_use_idemix_forcing,        &
    &  l_idemix_osborn_cox_kv,      &
    &  n_hor_iwe_prop_iter,         &
    &  fpath_iwe_surforc,           &
    &  name_iwe_surforc,            &
    &  fpath_iwe_botforc,           &
    &  name_iwe_botforc!,            &

  USE mo_ocean_physics_types, ONLY: t_ho_params, v_params, WindMixingDecay, WindMixingLevel
   !, l_convection, l_pp_scheme
  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_impl_constants,      ONLY: success, max_char_length, min_dolic, sea
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge,                           &
    &                               grid_unstructured_edge, grid_unstructured_cell
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_onEdges_Pointer_3d_wp, t_onCells_HalfLevels_Pointer_wp, t_operator_coeff
  USE mo_ocean_state,         ONLY: oce_config
  USE mo_physical_constants,  ONLY: grav, sitodbar,sal_ref
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_dynamics_config,     ONLY: nold!, nnew
  USE mo_run_config,          ONLY: dtime
  USE mo_var_list,            ONLY: add_var,                  &
    & add_ref, t_var_list_ptr
  USE mo_cf_convention
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: datatype_pack16, DATATYPE_FLT32, DATATYPE_FLT64, filetype_nc2, &
    &                               GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: &
    & za_depth_below_sea, za_depth_below_sea_half, za_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max, sync_patch_array_mult
  USE  mo_ocean_thermodyn,    ONLY: calculate_density_onColumn
  USE mo_ocean_math_operators,ONLY: div_oce_3d
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, &
    & timer_extra10, timer_extra11
  USE mo_statistics,          ONLY: global_minmaxmean
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_ocean_idemix_base,   ONLY: init_idemix, coeffs_idemix!, integrate_idemix
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_read_interface,      ONLY: read_2D_1Time, on_cells, t_stream_id, &
    & read_netcdf_broadcast_method, openInputFile, closeFile

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_idemix
  PUBLIC :: setup_idemix

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: iwe_surf_forc
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: iwe_bott_forc

CONTAINS


  SUBROUTINE setup_idemix(patch_3d)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_stream_id) :: stream_id 
    
    !! namelist variables
    !REAL(wp) :: &
    !  tau_v                     = 86400.0             ,& 
    !  tau_h                     = 1296000.0           ,&
    !  gamma                     = 1.570               ,&
    !  jstar                     = 10.0                ,&
    !  mu0                       = 1.33333333
    !! FIXME: these need to go to calc_idemix
    !LOGICAL :: &
    !  l_idemix_avo_dvo_direct   = .false.
    !INTEGER :: &
    !  n_hor_iwe_prop_iter       = 1
    !CHARACTER(LEN=120) :: fpath_iwe_surforc=''
    !CHARACTER(LEN=120) :: fpath_iwe_botforc=''
    !CHARACTER(LEN=40) :: name_iwe_surforc=''
    !CHARACTER(LEN=40) :: name_iwe_botforc=''

    ! other variables
    !REAL(wp) :: iwe_surf_forc(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: iwe_bott_forc(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL  :: has_missValue = .false.
    REAL(wp) :: missValue = -123456789999999999.0

    patch_2d   => patch_3d%p_patch_2d(1)

    ALLOCATE(iwe_surf_forc(nproma, patch_2d%alloc_cell_blocks))
    ALLOCATE(iwe_bott_forc(nproma, patch_2d%alloc_cell_blocks))
    iwe_surf_forc(:,:) = 0.0_wp
    iwe_bott_forc(:,:) = 0.0_wp

    fpath_iwe_surforc = 'idemix_surface_forcing.nc'
    !fpath_iwe_surforc = '/mnt/lustre01/work/mh0033/m300602/proj_vmix/icon/idemix_forcing/fourier_smooth_2005_cfsr_inert_OceanOnly_Icos_0158km_etopo40.nc'
    name_iwe_surforc = 'niw_forc' 

    fpath_iwe_botforc = 'idemix_bottom_forcing.nc'
    name_iwe_botforc = 'wave_dissipation' 

    ! initialise IDEMIX parameters
    CALL init_idemix(tau_v, tau_h, gamma, jstar, mu0)

    if ( l_use_idemix_forcing ) then 
      ! read internal wave surface forcing
      CALL openInputFile(stream_id, fpath_iwe_surforc, patch_2d, &
        &                read_netcdf_broadcast_method)
      CALL read_2D_1Time( stream_id=stream_id, location=on_cells, &
        & variable_name=name_iwe_surforc, fill_array=iwe_surf_forc,         &
        & has_missValue=has_missValue, missValue=missValue)
      CALL closeFile(stream_id)

      ! read internal wave bottom forcing
      CALL openInputFile(stream_id, fpath_iwe_botforc, patch_2d, &
        &                read_netcdf_broadcast_method)
      CALL read_2D_1Time( stream_id=stream_id, location=on_cells, &
        & variable_name=name_iwe_botforc, fill_array=iwe_bott_forc,         &
        & has_missValue=has_missValue, missValue=missValue)
      CALL closeFile(stream_id)
    else
      iwe_surf_forc(:,:) = 0.01_wp
      iwe_bott_forc(:,:) = 0.01_wp
    endif

    ! convert from W/m^2 to m^3/s^3
    ! (only 20% of the niw-input are available to penetrate into the deeper ocean)
    iwe_surf_forc(:,:) = iwe_surf_forc(:,:) / OceanReferenceDensity * 0.2
    iwe_bott_forc(:,:) = iwe_bott_forc(:,:) / OceanReferenceDensity

  END SUBROUTINE setup_idemix

  SUBROUTINE calc_idemix(patch_3d, ocean_state, params_oce, op_coeffs, atmos_fluxes)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_atmos_fluxes)                  :: atmos_fluxes
    !REAL(wp),          INTENT(in)         :: fu10   (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) ! t_atmos_for_ocean%fu10
    TYPE(t_ho_params), INTENT(inout)      :: params_oce
    TYPE(t_operator_coeff),INTENT(in)     :: op_coeffs
    !REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10

    ! pointer for convenience 
    REAL(wp), POINTER :: dz(:,:,:)
    REAL(wp), POINTER :: dzi(:,:,:)
    REAL(wp), POINTER :: dzw(:,:,:)
    REAL(wp), POINTER :: vned(:,:,:,:)
    !REAL(wp), POINTER :: vvel(:,:,:)
    REAL(wp), POINTER :: temp(:,:,:)
    REAL(wp), POINTER :: salt(:,:,:)
    REAL(wp), POINTER :: dens(:,:,:)
    REAL(wp), POINTER :: Av_old(:,:,:)
    REAL(wp), POINTER :: kv_old(:,:,:)
    REAL(wp), POINTER :: iwe(:,:,:)

    INTEGER, POINTER :: kbot(:,:)

    ! loop variables
    INTEGER :: jc, je, jk, blockNo, tracer_index, jtrc
    INTEGER :: start_index, end_index
    INTEGER :: levels

    !INTEGER :: cell_1_idx, cell_1_block, cell_2_idx, cell_2_block

    REAL(wp) :: rho_up(n_zlev), rho_down(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    REAL(wp) :: Nsqr(n_zlev+1), Ssqr(n_zlev+1)

    INTEGER :: tstep_count

    ! FIXME: what to do with that?
    REAL(wp) :: iwe_Av(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: iwe_kv(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    ! FIXME: this needs to be read in
    !REAL(wp) :: iwe_surf_forc(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: iwe_bott_forc(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    ! variables for horizontal iwe propagation
    REAL(wp) :: v0_e(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: v0_iwe_c(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flx_e(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_v0_iwe_e
    REAL(wp) :: div_flx_c(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    INTEGER :: n

    INTEGER,  DIMENSION(:,:,:),   POINTER :: idx_e, blk_e
    INTEGER,  DIMENSION(:,:,:),   POINTER :: idx_c, blk_c

    REAL(wp) :: dummy_zeros(n_zlev+1)
    LOGICAL  :: debug

    dummy_zeros = 0.0
    debug = .false.

    iwe => params_oce%vmix_params%iwe(:,:,:)
    iwe_kv = 0.0
    iwe_Av = 0.0

    dz  => patch_3d%p_patch_1d(1)%prism_center_dist_c
    dzi => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
    dzw => patch_3d%p_patch_1d(1)%prism_thick_c
    !uvel => ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x
    temp => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    ! FIXME: use sal_ref in case of temp is only tracer
    salt => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)

    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    ! renaming stuff
    patch_2D   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    all_cells  => patch_2D%cells%ALL
    levels = n_zlev

    !salinity(1:levels) = sal_ref
    rho_up(:)=0.0_wp    
    rho_down(:)=0.0_wp

    !write(*,*) "TKE before:"
    !write(*,*) tke(8,:,10)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        ! calculate N2
        pressure(1:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 1:levels, blockNo) * ReferencePressureIndbars
        !rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
        !  & salinity(1:levels-1), pressure(2:levels), levels-1)
        !rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
        !  & salinity(2:levels), pressure(2:levels), levels-1)
        rho_up(1:levels-1)  = calculate_density_onColumn(&
            & temp(jc,1:levels-1,blockNo), &
            & salt(jc,1:levels-1,blockNo), &
            & pressure(2:levels), levels-1)
        rho_down(2:levels)  = calculate_density_onColumn(&
            & temp(jc,2:levels,blockNo), &
            & salt(jc,2:levels,blockNo), &
            & pressure(2:levels), levels-1)
        Nsqr = 0.
        DO jk = 2, n_zlev 
          Nsqr(jk) = grav/OceanReferenceDensity * (rho_down(jk) - rho_up(jk-1)) *  dzi(jc,jk,blockNo)
        ENDDO

      !if (jc==8 .and. blockNo==10) then
      !  write(*,*) 'jc = ', jc, 'blockNo = ', blockNo, 'tstep_count = ', tstep_count
      !  !write(*,*) 'kbot(jc,blockNo) = ', patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
      !  write(*,*) 'kbot(jc,blockNo) = ', kbot(jc,blockNo)
      !  write(*,*) 'temp(jc,:,blockNo) = ', temp(jc,:,blockNo)
      !endif

      ! main call to calculate idemix
      if (kbot(jc,blockNo)>0) then
      CALL coeffs_idemix(                                                            &
                   ! parameter
                   dzw             = dzw(jc,:,blockNo),                              &
                   dzt             = dz(jc,:,blockNo),                               &
                   nlev            = kbot(jc,blockNo),                               &
                   max_nlev        = n_zlev,                                         & 
                   dtime           = dtime,                                          &
                   ! FIXME: get name of Coriolis param
                   coriolis        = patch_2D%cells%f_c(jc,blockNo),                 &
                   ! essentials 
                   iwe_old         = iwe(jc,:,blockNo),                              & ! in
                   iwe_new         = iwe(jc,:,blockNo),                              & ! out
                   forc_iw_surface = iwe_surf_forc(jc,blockNo),                      & ! in
                   forc_iw_bottom  = iwe_bott_forc(jc,blockNo),                      & ! in
                   ! FIXME: nils: better output IDEMIX Ri directly
                   alpha_c         = params_oce%vmix_params%iwe_alpha_c(jc,:,blockNo), & ! out (for Ri IMIX)
                   ! only for Osborn shortcut 
                   ! FIXME: nils: put this to tke
                   KappaM_out      = iwe_Av(jc,:,blockNo),                           & ! out
                   KappaH_out      = iwe_kv(jc,:,blockNo),                           & ! out
                   Nsqr            = Nsqr(:),                                        & ! in
                   ! diagnostics
                   iwe_Ttot        = params_oce%vmix_params%iwe_Ttot(jc,:,blockNo), &
                   iwe_Tdif        = params_oce%vmix_params%iwe_Tdif(jc,:,blockNo), &
                   iwe_Thdi        = params_oce%vmix_params%iwe_Thdi(jc,:,blockNo), &
                   iwe_Tdis        = params_oce%vmix_params%iwe_Tdis(jc,:,blockNo), &
                   iwe_Tsur        = params_oce%vmix_params%iwe_Tsur(jc,:,blockNo), &
                   iwe_Tbot        = params_oce%vmix_params%iwe_Tbot(jc,:,blockNo), &
                   c0              = params_oce%vmix_params%iwe_c0(jc,:,blockNo),   &
                   v0              = params_oce%vmix_params%iwe_v0(jc,:,blockNo),   &
                   ! debugging
                   debug           = debug,                                          &
                   !i = jc,                                                           &
                   !j = blockNo,                                                      &
                   !tstep_count = tstep_count,                                        &
                   int_1     = params_oce%vmix_params%vmix_dummy_1(jc,:,blockNo), &
                   int_2     = params_oce%vmix_params%vmix_dummy_2(jc,:,blockNo), &
                   int_3     = params_oce%vmix_params%vmix_dummy_3(jc,:,blockNo)  &
                   )

      end if
    !if (jc==8 .and. blockNo==10) then
    !!  write(*,*) params_oce%vmix_params%vmix_dummy_1(jc,:,blockNo)
    !  write(*,*) 'dzw = ', dzw(jc,:,blockNo)
    !  write(*,*) 'dzt = ', dz(jc,:,blockNo)
    !  stop
    !end if
      ENDDO
    ENDDO

! start: iwe hor prop 
! continue here
    ! add contribution from horizontal wave propagation
    if (n_hor_iwe_prop_iter>0) then

      CALL sync_patch_array(sync_c, patch_2D, iwe)
      CALL sync_patch_array(sync_c, patch_2D, params_oce%vmix_params%iwe_v0)

      ! temporarily store old iwe values for diag
      params_oce%vmix_params%iwe_Thdi(:,:,:) = iwe(:,:,:) 

      ! restrict iwe_v0 to fullfill stability criterium
      !DO j=1,je
      !  DO i=1,ie-1
      !    DO k=1,ke+1
      !      fxa = sqrt(   0.2_wp*min(dlxu(i,j),dlyu(i,j))**2 &
      !                  / (tau_h*dt/n_hor_iwe_prop_iter) )
      !      iwe_v0(i,j,k) = min(iwe_v0(i,j,k), fxa)
      !    ENDDO
      !  ENDDO
      !ENDDO

      ! interpolate hor. group velocity cell center to edges
      v0_e(:,:,:) = 0.0_wp
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      !DO blockNo = patch_2D%edges%all%start_block, patch_2D%edges%all%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        ! FIXME: in or outside block loop?
        ! FIXME: Shouldn't this be done once outside of all loops?
        ! FIXME: probably idx_e and idx_c necessary
        idx_c => patch_3D%p_patch_2D(1)%edges%cell_idx
        blk_c => patch_3D%p_patch_2D(1)%edges%cell_blk
        DO je = start_index, end_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)+1
            v0_e(je,jk,blockNo) = &
              & 0.5_wp * (    params_oce%vmix_params%iwe_v0(  &
              &   idx_c(je,blockNo,1),jk,blk_c(je,blockNo,1) ) &
              &             + params_oce%vmix_params%iwe_v0(  &
              &   idx_c(je,blockNo,2),jk,blk_c(je,blockNo,2) ) )
          ENDDO
        ENDDO
      ENDDO

      ! sub-timestepping of horizontal iwe prop.
      DO n=1,n_hor_iwe_prop_iter

        ! derive product v0*iwe
        v0_iwe_c(:,:,:) = params_oce%vmix_params%iwe_v0(:,:,:)*iwe(:,:,:)
        !v0_iwe_c(:,:,:) = 0.0_wp
        !DO blockNo = all_cells%start_block, all_cells%end_block
        !  CALL get_index_range(all_cells, blockNo, start_index, end_index)
        !  DO jc = start_index, end_index
        !    DO jk = 1, n_zlev+1
        !        v0_iwe_c(jc,jk,blockNo) = v0(jc,jk,blockNo)*iwe(jc,jk,blockNo)
        !    ENDDO
        !  ENDDO
        !ENDDO

        ! derive flux
        !grad_v0_iwe_e(:,:,:) = 0.0_wp
        flx_e(:,:,:) = 0.0_wp
        ! FIXME: in or outside block loop?
        ! FIXME: Shouldn't this be done once outside of all loops?
        ! FIXME: probably idx_e and idx_c necessary
        idx_c => patch_3D%p_patch_2D(1)%edges%cell_idx
        blk_c => patch_3D%p_patch_2D(1)%edges%cell_blk
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        !DO blockNo = patch_2D%edges%all%start_block, patch_2D%edges%all%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
          DO je = start_index, end_index
            ! leave out last vert. iwe layer since grad_coeff is not defined there
            DO jk = 1, patch_3D%p_patch_1d(1)%dolic_e(je,blockNo)
              ! FIXME: how can we include correct dz here?
              grad_v0_iwe_e =  op_coeffs%grad_coeff(je,jk,blockNo) *                 &
                & (  v0_iwe_c(idx_c(je,blockNo,2),jk,blk_c(je,blockNo,2))        &
                &  - v0_iwe_c(idx_c(je,blockNo,1),jk,blk_c(je,blockNo,1)) )
              ! FIXME: How to asure that flux through solid boundaries is zero?
              flx_e(je,jk,blockNo) = tau_h*v0_e(je,jk,blockNo)*grad_v0_iwe_e
            ENDDO
          ENDDO
        ENDDO

        ! sync (better save than sorry)
        CALL sync_patch_array(sync_e, patch_2D, flx_e)

        !! derive flux
        !flx_e(:,:,:) = 0.0_wp
        !DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        !  CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        !  DO je = start_index, end_index
        !    DO jk = 1, patch_3D%p_patch_1d(1)%dolic_e(je,blockNo)+1
        !      ! FIXME: How to asure that flux through solid boundaries is zero?
        !      flx_e(je,jk,blockNo) = tau_h*v0_e(je,jk,blockNo)*grad_v0_times_iwe(je,jk,blockNo)
        !    ENDDO
        !  ENDDO
        !ENDDO

        ! divergence of flux
        div_flx_c(:,:,:) = 0.0_wp
        ! FIXME: in or outside block loop?
        ! FIXME: Shouldn't this be done once outside of all loops?
        ! FIXME: probably idx_e and idx_c necessary
        idx_e => patch_3D%p_patch_2D(1)%cells%edge_idx
        blk_e => patch_3D%p_patch_2D(1)%cells%edge_blk
        DO blockNo = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, blockNo, start_index, end_index)
          DO jc = start_index, end_index
            ! leave out last vert. iwe layer since grad_coeff is not defined there
            DO jk = 1, patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
              ! FIXME: how can we include correct dz here?
              div_flx_c(jc,jk,blockNo) =  &
                & flx_e(idx_e(jc,blockNo,1),jk,blk_e(jc,blockNo,1)) * op_coeffs%div_coeff(jc,jk,blockNo,1) + &
                & flx_e(idx_e(jc,blockNo,2),jk,blk_e(jc,blockNo,2)) * op_coeffs%div_coeff(jc,jk,blockNo,2) + &
                & flx_e(idx_e(jc,blockNo,3),jk,blk_e(jc,blockNo,3)) * op_coeffs%div_coeff(jc,jk,blockNo,3)
            ENDDO
          ENDDO
        ENDDO

        ! apply tendency
        ! FIXME: Do we need to divide by cell volume?
        iwe(:,:,:) = iwe(:,:,:) + dtime/n_hor_iwe_prop_iter * div_flx_c(:,:,:)
        CALL sync_patch_array(sync_c, patch_2D, iwe)

      ENDDO ! n=1,n_hor_iwe_prop_iter

      ! derive diagnostic of hor. wave propagation and total tendency
      params_oce%vmix_params%iwe_Thdi(:,:,:) = (iwe(:,:,:) - params_oce%vmix_params%iwe_Thdi(:,:,:))/dtime
      params_oce%vmix_params%iwe_Ttot(:,:,:) = params_oce%vmix_params%iwe_Ttot(:,:,:) + params_oce%vmix_params%iwe_Thdi(:,:,:)
    end if
! end: iwe hor prop 

    if ( l_idemix_osborn_cox_kv ) then
      ! write tke vert. diffusivity to vert tracer diffusivities
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          ! FIXME: nils: make loop over all tracer
          DO jtrc = 1, no_tracer
            params_oce%a_tracer_v(jc,:,blockNo,jtrc) = iwe_kv(jc,:,blockNo)
          ENDDO
        ENDDO
      ENDDO
    end if

    !write(*,*) "Stopping..."
    !stop
    !write(*,*) "TKE after:"
    !write(*,*) tke(8,:,10)
    
    !params_oce%a_tracer_v = 1e-5
    !write(*,*) 'a_tracer_v = ', params_oce%a_tracer_v(8,:,10,1) 
    !write(*,*) 'kbot = ', kbot(8,10) 
    !write(*,*) 'tke_kv = ', tke_kv(8,:,10) 
    !write(*,*) 'tke = ', tke(8,:,10) 
    !write(*,*) 'fu10 = ', fu10(8,10) 
  END SUBROUTINE calc_idemix
  
END MODULE mo_ocean_idemix
