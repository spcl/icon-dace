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

!NEC$ options "-O0"

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_tke
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
    & vert_mix_type, vmix_pp, vmix_tke, vmix_idemix_tke,      & ! by_nils
    & HorizontalViscosity_SmoothIterations,                   &
    & convection_InstabilityThreshold,                        &
    & RichardsonDiffusion_threshold,                          &
    & lambda_wind,                                            &
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
    &  ReferencePressureIndbars,    &
    ! tke parameter
    &  c_k,                         &
    &  c_eps,                       &
    &  alpha_tke,                   &
    &  mxl_min,                     &
    &  use_Kappa_min,               &
    &  KappaM_min,                  &
    &  KappaH_min,                  &
    &  KappaM_max,                  &
    &  cd,                          &
    &  tke_min,                     &
    &  tke_mxl_choice,              &
    &  tke_surf_min,                &
    &  only_tke,                    &
    &  l_lc,                        & !by_Oliver
    &  clc,                         & !by_Oliver 
    &  use_ubound_dirichlet,        &
    &  use_lbound_dirichlet,        &
    &  vert_cor_type

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
  USE mo_math_constants,      ONLY: dbl_eps,pi
  USE mo_dynamics_config,     ONLY: nold!, nnew
  USE mo_run_config,          ONLY: dtime
  USE mo_cf_convention
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: datatype_pack16, DATATYPE_FLT32, DATATYPE_FLT64, filetype_nc2, &
    &                               GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: &
    & za_depth_below_sea, za_depth_below_sea_half, za_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max, sync_patch_array_mult
  USE  mo_ocean_thermodyn,    ONLY: calculate_density_onColumn, calculate_density_onColumn_elem
  USE mo_ocean_math_operators,ONLY: div_oce_3d
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, &
    & timer_extra10, timer_extra11
  USE mo_statistics,          ONLY: global_minmaxmean
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_ocean_tke_base,      ONLY: init_tke, coeffs_tke!, integrate_tke
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  !USE turb_data,              ONLY: rhon ! air density
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_tke
  PUBLIC :: setup_tke

  INTERFACE calc_tke
#if defined(__LVECTOR__) && !defined(__LVEC_BITID__)
    MODULE PROCEDURE calc_tke_vector
#else
    MODULE PROCEDURE calc_tke_scalar
#endif
  END INTERFACE

  CHARACTER(LEN=*), PARAMETER :: module_name = 'mo_ocean_tke'


CONTAINS


  SUBROUTINE setup_tke()
    ! FIXME: So far this routine is not called!!!
    ! If it is called delet definitions from mo_ocean_tke_base
    !REAL(wp) :: c_k
    !REAL(wp) :: c_eps
    !REAL(wp) :: alpha_tke
    !REAL(wp) :: mxl_min
    !LOGICAL  :: use_Kappa_min
    !REAL(wp) :: KappaM_min
    !REAL(wp) :: KappaH_min
    !REAL(wp) :: KappaM_max
    !REAL(wp) :: cd
    !REAL(wp) :: tke_min
    !INTEGER  :: tke_mxl_choice
    !REAL(wp) :: tke_surf_min
    !LOGICAL  :: only_tke
    !LOGICAL  :: use_ubound_dirichlet
    !LOGICAL  :: use_lbound_dirichlet

    !! FIXME: In the end this needs to go to the namelist
    !c_k        = 0.1_wp
    !c_eps      = 0.7_wp
    !alpha_tke  = 30.0_wp
    !mxl_min    = 1.d-8_wp
    !use_Kappa_min = .false.
    !KappaM_min = 1.d-4_wp
    !KappaH_min = 1.d-5_wp
    !KappaM_max = 100.0_wp
    !cd         = 3.75_wp
    !tke_min    = 1.d-6_wp
    !tke_mxl_choice = 2
    !tke_surf_min = 1.d-4_wp
    !only_tke = .true.
    !use_ubound_dirichlet = .false.
    !use_lbound_dirichlet = .false.

    !write(*,*) c_k
    !write(*,*) c_eps
    !write(*,*) alpha_tke
    !write(*,*) mxl_min
    !write(*,*) use_Kappa_min
    !write(*,*) KappaM_min
    !write(*,*) KappaH_min
    !write(*,*) KappaM_max
    !write(*,*) cd
    !write(*,*) tke_min
    !write(*,*) tke_mxl_choice
    !write(*,*) tke_surf_min
    !write(*,*) only_tke
    !write(*,*) use_ubound_dirichlet
    !write(*,*) use_lbound_dirichlet

    CALL init_tke(c_k            = c_k,            &
                  c_eps          = c_eps,          &
                  cd             = cd,             &
                  alpha_tke      = alpha_tke,      &
                  mxl_min        = mxl_min,        &
                  use_Kappa_min  = use_Kappa_min,  &
                  KappaM_min     = KappaM_min,     &
                  KappaH_min     = KappaH_min,     &
                  KappaM_max     = KappaM_max,     &
                  tke_mxl_choice = tke_mxl_choice, &
                  use_ubound_dirichlet = use_ubound_dirichlet, &
                  use_lbound_dirichlet = use_lbound_dirichlet, &
                  only_tke       = only_tke,       &
                  l_lc           = l_lc,           &
                  clc            = clc,            &
                  tke_min        = tke_min,        &
                  tke_surf_min   = tke_surf_min    )
                  !tke_userdef_constants=tke_userdef_constants)
  END SUBROUTINE setup_tke

#ifdef _OPENACC
  SUBROUTINE calc_tke_scalar(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_atmos_fluxes)                  :: atmos_fluxes
    TYPE(t_ho_params), INTENT(inout)      :: params_oce
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum   (:,:) ! sea ice concentration
    LOGICAL, INTENT(IN), OPTIONAL        :: lacc

    ! pointer for convenience 
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells

    REAL(wp), POINTER :: prism_thick_c(:,:,:)
    INTEGER,  POINTER :: dolic_c(:,:)

    REAL(wp), POINTER :: dz(:,:,:)
    REAL(wp), POINTER :: dzi(:,:,:)
    REAL(wp), POINTER :: vned(:,:,:,:)
    REAL(wp), POINTER :: temp(:,:,:)
    REAL(wp), POINTER :: salt(:,:,:)
    REAL(wp), POINTER :: dens(:,:,:)
    REAL(wp), POINTER :: Av_old(:,:,:)
    REAL(wp), POINTER :: kv_old(:,:,:)
    REAL(wp), POINTER :: tke(:,:,:)
 
    ! Langmuir turbulence   
    REAL(wp), POINTER :: tke_plc(:,:,:)
    REAL(wp), POINTER :: tke_plc_ptr(:,:)
    REAL(wp), POINTER :: wlc(:,:,:)
    REAL(wp), POINTER :: hlc(:,:)
    REAL(wp), POINTER :: u_stokes(:,:)
    REAL(wp), POINTER :: depth_CellInterface(:,:,:)


    ! loop variables
    INTEGER :: jc, blockNo, je,jk, tracer_index, jtrc
    INTEGER :: start_index, end_index
    INTEGER :: levels

    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block

    REAL(wp) :: rho_up(nproma,n_zlev), rho_down(nproma,n_zlev)
    REAL(wp) :: pressure(nproma,n_zlev)
    REAL(wp) :: Nsqr(nproma,n_zlev+1), Ssqr(nproma,n_zlev+1), tmp, tke_old(nproma,n_zlev+1)

    ! Parameters for zstar
    REAL(wp) :: s_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    ! stretching factor


    INTEGER :: tstep_count, level

    ! put this later to a global place
    REAL(wp) :: tke_Av(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_kv(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iw_alpha_c(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iwe(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iwe_forcing(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: forc_tke_surf_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: forc_rho_surf_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: bottom_fric_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tmp_dzw(nproma, n_zlev)
    REAL(wp) :: tmp_dzt(nproma, n_zlev+1)
    REAL(wp) :: tau_abs

    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    tke => params_oce%vmix_params%tke(:,:,:)

    !$ACC DATA CREATE(tke_Av, tke_kv, tke_iw_alpha_c, tke_iwe, tke_iwe_forcing, forc_tke_surf_2D) &
    !$ACC   CREATE(forc_rho_surf_2D, bottom_fric_2D, tmp_dzw, tmp_dzt, s_c) &
    !$ACC   CREATE(rho_up, rho_down, pressure, Nsqr, Ssqr, tke_old) &
    !$ACC   IF(lzacc)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    tke_kv(:,:,:) = 0.0
    tke_Av(:,:,:) = 0.0
    tke_iw_alpha_c(:,:,:) = 0.0
    tke_iwe(:,:,:) = 0.0
    tke_iwe_forcing(:,:,:) = 0.0
    !$ACC END KERNELS

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    forc_rho_surf_2D(:,:) = 0.0
    bottom_fric_2D(:,:) = 0.0
    !$ACC END KERNELS
    !$ACC WAIT(1)

    IF(l_lc) THEN
      ! Langmuir turbulence variables
      tke_plc  => params_oce%vmix_params%tke_plc(:,:,:)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      tke_plc(:,:,:) = 0.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)
      hlc      => params_oce%vmix_params%hlc(:,:)
      wlc      => params_oce%vmix_params%wlc(:,:,:)
      u_stokes => params_oce%vmix_params%u_stokes(:,:)
      depth_CellInterface => patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,:)
    ENDIF

    dz  => patch_3d%p_patch_1d(1)%prism_center_dist_c
    dzi => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
    prism_thick_c => patch_3d%p_patch_1d(1)%prism_thick_c
    !uvel => ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x
    temp => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    ! FIXME: use sal_ref in case of temp is only tracer
    salt => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)

    Av_old => params_oce%a_veloc_v(:,:,:)
    ! Use a_tracer_v of temperature here
    kv_old => params_oce%a_tracer_v(:,:,:,1)
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    ! renaming stuff
    patch_2D   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    all_cells  => patch_2D%cells%ALL

    ! special settings if IDEMIX is used together with TKE
    IF ( vert_mix_type==vmix_idemix_tke ) THEN
      ! use iwe dissipation as forcing for tke
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      tke_iwe_forcing(:,:,:) = -1.0_wp * params_oce%vmix_params%iwe_Tdis(:,:,:)
      !$ACC END KERNELS
      !$ACC WAIT(1)
    ENDIF

    ! set zstar related parameters
    IF (vert_cor_type == 1) THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      s_c(:,:) = ocean_state%p_prog(nold(1))%stretch_c(:,:)
      !$ACC END KERNELS
      !$ACC WAIT(1)
    ELSE
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      s_c(:,:) = 1.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)
    ENDIF

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO jc = 1,nproma
      pressure(jc,1) = 1.0_wp
      pressure(jc,2:n_zlev) = patch_3d%p_patch_1d(1)%zlev_i(2:n_zlev) * ReferencePressureIndbars
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        levels = dolic_c(jc,blockNo)
        IF (dolic_c(jc,blockNo) > 0) THEN

          tke_old(jc,:) = tke(jc,:,blockNo)

          rho_up(jc,:)=0.0_wp
          rho_down(jc,:)=0.0_wp
          Nsqr(jc,:) = 0.0_wp
          Ssqr(jc,:) = 0.0_wp

          ! wind stress for tke surface forcing (reduced under sea ice)
          IF (use_reduced_mixing_under_ice) THEN
            tau_abs = (1.0_wp - concsum(jc,blockNo))**2 &
                 &    * SQRT((atmos_fluxes%stress_xw(jc,blockNo))**2 &
                 &    + (atmos_fluxes%stress_yw(jc,blockNo))**2 )
          ELSE
            tau_abs = (1.0_wp - concsum(jc,blockNo)) &
                 &    * SQRT((atmos_fluxes%stress_xw(jc,blockNo))**2 &
                 &    + (atmos_fluxes%stress_yw(jc,blockNo))**2 )
          ENDIF

          forc_tke_surf_2D(jc,blockNo) = tau_abs / OceanReferenceDensity

          ! calculate N2
          DO level = 1,levels-1
            rho_up(jc,level) = calculate_density_onColumn_elem(temp(jc,level,blockNo), salt(jc,level,blockNo), &
                                                pressure(jc,level+1))
          END DO

          DO level = 2,levels
            rho_down(jc,level) = calculate_density_onColumn_elem(temp(jc,level,blockNo), salt(jc,level,blockNo), &
                                                pressure(jc,level))
          END DO

          DO jk = 2, levels
            Nsqr(jc,jk) = grav/OceanReferenceDensity * (rho_down(jc,jk) - rho_up(jc,jk-1)) *  dzi(jc,jk,blockNo) / s_c(jc,blockNo)
          ENDDO

          ! calculate shear
          DO jk = 2, levels
            Ssqr(jc,jk) = SUM(((  ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x   &
                 &              - ocean_state%p_diag%p_vn(jc,jk,  blockNo)%x ) &
                 &            * dzi(jc,jk,blockNo)/s_c(jc,blockNo) )**2)
          ENDDO

          IF (l_lc) THEN
            ! calculate Langmuir cell additional term after Axell (2002)

            ! calculate Stoke's drift
            ! Approximation if there is no information about the wave field
            ! As done in Nemo
            ! FIXME: do we need to divide tau by rho?

            ! Option used in NEMO model (https://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf, p.197) see also Breivik et al. (2015)
            ! They assume rhoair=1.2 kg/m3 and cd=1.5e-03:
            ! u_stokes = 0.016/(1.2 * 1.5e-03)^0.5 * |tau|^0.5; although they seem to use rhoair=1.2 kg/m3
            !u_stokes(jc,blockNo) = 0.377_wp * SQRT(tau_abs)                       ! [tau]=N2/m2
            u_stokes(jc,blockNo) = 0.016_wp/SQRT(1.2_wp * 1.5e-03_wp)*SQRT(tau_abs)           ! [tau]=N2/m2, rhoair=1.2, cd=1.5*10e-03

            ! This is done in Coulevard et al (2020, doi:10.5194/gmd-13-3067-2020), see Fig.2
            ! u_stokes(jc,blockNo) = 0.377_wp * SQRT(forc_tke_surf_2D(jc,blockNo))

            ! other option from Li and Garrett (1993)
            !u_stokes(jc,blockNo) = 0.016_wp * fu10(jc,blockNo)

            ! or original version from Axell (2002)
            !LLC = 0.12_wp*(u10**2/g)
            !u_stokes(jc,:,blockNo) = 0.0016*u10*EXP(depth/LLC)


            ! find depth of langmuir cell (hlc). hlc is the depth to which a water
            ! parcel with kinetic energy 0.5*u_stokes**2 can reach on its own by
            ! converting its kinetic energy to potential energy.
            hlc(jc,blockNo) = 0.0_wp
            DO jk = 2, levels
              tmp = SUM( Nsqr(jc,2:jk)*depth_CellInterface(jc,2:jk,blockNo) &  ! see axell (2002) eq.47
                           *prism_thick_c(jc,2:jk,blockNo)*s_c(jc,blockNo) )

              IF(tmp > 0.5_wp*u_stokes(jc,blockNo)**2.0_wp) THEN
                hlc(jc,blockNo) = depth_CellInterface(jc,jk,blockNo)*s_c(jc,blockNo)
                EXIT
              ENDIF
            ENDDO

            ! calculate langmuir cell velocity scale (wlc)
            ! Note: Couvelard et al (2020) set clc=0.3 instead of default 0.15 from
            ! Axell (2002); results in deeper MLDs and better spatial MLD pattern.
            DO jk = 2, levels
              IF ( ( depth_CellInterface(jc,jk,blockNo)*s_c(jc,blockNo) <= hlc(jc,blockNo) )  &
                   .AND. ( patch_3D%wet_c(jc,jk,blockNo) .EQ. 1 ) ) THEN
                wlc(jc,jk,blockNo) = clc * u_stokes(jc,blockNo)*SIN(pi*depth_CellInterface(jc,jk,blockNo)/hlc(jc,blockNo))
              ELSE
                wlc(jc,jk,blockNo) = 0.0_wp
              ENDIF
            ENDDO

            ! calculate langmuir turbulence term (tke_plc)
            IF (hlc(jc,blockNo) > 0.0_wp) THEN
              tke_plc(jc,1:levels,blockNo) = wlc(jc,1:levels,blockNo)**3.0_wp / hlc(jc,blockNo)
            ELSE
              tke_plc(jc,1:levels,blockNo) = 0.0_wp
            ENDIF

          ENDIF
        ENDIF
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        tmp_dzw(jc,:) = prism_thick_c(jc,:,blockNo)*s_c(jc,blockNo)
        tmp_dzt(jc,:) = dz(jc,:,blockNo)*s_c(jc,blockNo)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
      IF (l_lc) THEN
        tke_plc_ptr => tke_plc(:,:,blockNo)
      ELSE
        tke_plc_ptr => NULL()
      END IF
      CALL coeffs_tke(                                                                        &
                            start_index = start_index,                                        &
                            end_index = end_index,                                            &
                            nproma = nproma,                                                  &
                            levels = dolic_c(:,blockNo),                                      &
                            tke_old = tke_old(:,:),                                           &
                            tke_new = tke(:,:,blockNo),                                       &
                            KappaM_out = tke_Av(:,:,blockNo),                                 &
                            KappaH_out = tke_kv(:,:,blockNo),                                 &
                            vmix_int_1 = params_oce%vmix_params%vmix_dummy_1(:,:,blockNo),    &
                            vmix_int_2 = params_oce%vmix_params%vmix_dummy_2(:,:,blockNo),    &
                            vmix_int_3 = params_oce%vmix_params%vmix_dummy_3(:,:,blockNo),    &
                            dzw = tmp_dzw(:,:),                                               &
                            dzt = tmp_dzt(:,:),                                               &
                            max_nlev = n_zlev,                                                &
                            Ssqr = Ssqr(:,:),                                                 &
                            Nsqr = Nsqr(:,:),                                                 &
                            tke_Tbpr = params_oce%vmix_params%tke_Tbpr(:,:,blockNo),          &
                            tke_Tspr = params_oce%vmix_params%tke_Tspr(:,:,blockNo),          &
                            tke_Tdif = params_oce%vmix_params%tke_Tdif(:,:,blockNo),          &
                            tke_Tdis = params_oce%vmix_params%tke_Tdis(:,:,blockNo),          &
                            tke_Twin = params_oce%vmix_params%tke_Twin(:,:,blockNo),          &
                            tke_Tiwf = params_oce%vmix_params%tke_Tiwf(:,:,blockNo),          &
                            tke_Tbck = params_oce%vmix_params%tke_Tbck(:,:,blockNo),          &
                            tke_Ttot = params_oce%vmix_params%tke_Ttot(:,:,blockNo),          &
                            tke_Lmix = params_oce%vmix_params%tke_Lmix(:,:,blockNo),          &
                            tke_Pr   = params_oce%vmix_params%tke_Pr(:,:,blockNo),            &
                            tke_plc  = tke_plc_ptr,                                           & !by_Oliver
                            forc_tke_surf = forc_tke_surf_2D(:,blockNo),                      &
                            E_iw = tke_iwe(:,:,blockNo),                                      &
                            dtime = dtime,                                                    &
                            iw_diss = tke_iwe_forcing(:,:,blockNo),                           &
                            forc_rho_surf = forc_rho_surf_2D(:,blockNo),                      &
                            rho_ref = OceanReferenceDensity,                                  &
                            grav = grav,                                                      &
                            alpha_c = tke_iw_alpha_c(:,:,blockNo),                            &
                            lacc = lzacc                                                      &
                            )
    END DO

    ! interpolate vert. visosity from cell center to edges
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      params_oce%a_veloc_v(:,:,blockNo) = 0.0
      !$ACC END KERNELS

      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_index, end_index
        levels       = patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
        cell_1_idx   = patch_2D%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
        cell_2_idx   = patch_2D%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
        DO jk = 2, levels
          params_oce%a_veloc_v(je,jk,blockNo) = &
               & 0.5_wp * (    tke_Av(cell_1_idx,jk,cell_1_block) &
               &             + tke_Av(cell_2_idx,jk,cell_2_block) )
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)

    ! write tke vert. diffusivity to vert tracer diffusivities
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        ! FIXME: nils: make loop over all tracer
        DO jtrc = 1, no_tracer
          params_oce%a_tracer_v(jc,:,blockNo,jtrc) = tke_kv(jc,:,blockNo)
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)

    !$ACC END DATA

  END SUBROUTINE calc_tke_scalar

#else

  SUBROUTINE calc_tke_scalar(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_atmos_fluxes)                  :: atmos_fluxes
    !REAL(wp),          INTENT(in)         :: fu10   (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) ! t_atmos_for_ocean%fu10
    TYPE(t_ho_params), INTENT(inout)      :: params_oce
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum   (:,:) ! sea ice concentration

    ! pointer for convenience 
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells

    REAL(wp), POINTER :: prism_thick_c(:,:,:)
    INTEGER,  POINTER :: dolic_c(:,:)

    REAL(wp), POINTER :: dz(:,:,:)
    REAL(wp), POINTER :: dzi(:,:,:)
    REAL(wp), POINTER :: vned(:,:,:,:)
    !REAL(wp), POINTER :: vvel(:,:,:)
    REAL(wp), POINTER :: temp(:,:,:)
    REAL(wp), POINTER :: salt(:,:,:)
    REAL(wp), POINTER :: dens(:,:,:)
    REAL(wp), POINTER :: Av_old(:,:,:)
    REAL(wp), POINTER :: kv_old(:,:,:)
    REAL(wp), POINTER :: tke(:,:,:)
 
    ! Langmuir turbulence   
    REAL(wp), POINTER :: tke_plc(:,:,:)
    REAL(wp), POINTER :: tke_plc_ptr(:)
    REAL(wp), POINTER :: wlc(:,:,:)
    REAL(wp), POINTER :: hlc(:,:)
    REAL(wp), POINTER :: u_stokes(:,:)
    REAL(wp), POINTER :: depth_CellInterface(:,:,:)


    ! loop variables
    INTEGER :: jc, blockNo, je,jk, tracer_index, jtrc
    INTEGER :: start_index, end_index
    INTEGER :: levels

    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block

    REAL(wp) :: rho_up(n_zlev), rho_down(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    REAL(wp) :: Nsqr(n_zlev+1), Ssqr(n_zlev+1), tmp, tke_old(n_zlev+1)
    !REAL(wp), POINTER :: vert_density_grad(:,:,:)

    ! Parameters for zstar
    REAL(wp) :: s_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    ! stretching factor


    INTEGER :: tstep_count

    ! put this later to a global place
    REAL(wp) :: tke_Av(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_kv(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tbpr(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tspr(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tdif(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tdis(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Twin(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tiwf(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tbck(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Ttot(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Lmix(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Pr(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: dummy_1(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: dummy_2(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: dummy_3(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iw_alpha_c(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iwe(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iwe_forcing(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: forc_tke_surf_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: forc_rho_surf_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: bottom_fric_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tau_abs

!    REAL(wp) :: minmaxmean(3)

    REAL(wp) :: dummy_zeros(n_zlev+1)
    dummy_zeros = 0.0

    !tke = 0.0
    tke => params_oce%vmix_params%tke(:,:,:)
    !tke_Tbpr = 0.0
    !tke_Tspr = 0.0
    !tke_Tdif = 0.0
    !tke_Tdis = 0.0
    !tke_Twin = 0.0
    !tke_Tiwf = 0.0
    !tke_Tbck = 0.0
    !tke_Ttot = 0.0
    !tke_Lmix = 0.0
    !tke_Pr = 0.0
    !dummy_1 = 0.0
    !dummy_2 = 0.0
    !dummy_3 = 0.0
    tke_kv = 0.0
    tke_Av = 0.0
    tke_iw_alpha_c = 0.0
    tke_iwe = 0.0
    tke_iwe_forcing = 0.0
    forc_rho_surf_2D = 0.0
    bottom_fric_2D = 0.0

    IF(l_lc) THEN
      ! Langmuir turbulence variables
      tke_plc  => params_oce%vmix_params%tke_plc(:,:,:)
      tke_plc = 0.0_wp
      hlc      => params_oce%vmix_params%hlc(:,:)
      wlc      => params_oce%vmix_params%wlc(:,:,:)
      u_stokes => params_oce%vmix_params%u_stokes(:,:)
      depth_CellInterface => patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,:)
    ENDIF

    dz  => patch_3d%p_patch_1d(1)%prism_center_dist_c
    dzi => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
    prism_thick_c => patch_3d%p_patch_1d(1)%prism_thick_c
    !uvel => ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x
    temp => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    ! FIXME: use sal_ref in case of temp is only tracer
    salt => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)

    Av_old => params_oce%a_veloc_v(:,:,:)
    ! Use a_tracer_v of temperature here
    kv_old => params_oce%a_tracer_v(:,:,:,1)
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    ! renaming stuff
    patch_2D   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    all_cells  => patch_2D%cells%ALL

    ! special settings if IDEMIX is used together with TKE
    IF ( vert_mix_type==vmix_idemix_tke ) THEN
      ! use iwe dissipation as forcing for tke
      tke_iwe_forcing(:,:,:) = -1.0_wp * params_oce%vmix_params%iwe_Tdis(:,:,:)
    ENDIF

    ! set zstar related parameters
    IF (vert_cor_type == 1) THEN
      s_c(:,:) = ocean_state%p_prog(nold(1))%stretch_c(:,:)
    ELSE
      s_c(:,:) = 1.0_wp
    ENDIF

    !write(*,*) "TKE before:"
    !write(*,*) tke(8,:,10)

!    minmaxmean(:) = global_minmaxmean(values=concsum(:,:), in_subset=all_cells)
!    CALL debug_print_MaxMinMean('concsum',   minmaxmean, 'calc_tke', 1)

!ICON_OMP_PARALLEL PRIVATE(pressure)

    pressure(1) = 1.0_wp
    pressure(2:n_zlev) = patch_3d%p_patch_1d(1)%zlev_i(2:n_zlev) * ReferencePressureIndbars

!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, levels, &
!ICON_OMP  tau_abs, Nsqr, Ssqr, rho_up, rho_down, tstep_count, tmp, tke_old, &
!ICON_OMP blockNo, tke_plc_ptr)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        levels = dolic_c(jc,blockNo)
        IF (dolic_c(jc,blockNo) > 0) THEN

          tke_old(:) = tke(jc,:,blockNo)

          rho_up(:)=0.0_wp
          rho_down(:)=0.0_wp
          Nsqr(:) = 0.0_wp
          Ssqr(:) = 0.0_wp

          ! wind stress for tke surface forcing (reduced under sea ice)

          IF (use_reduced_mixing_under_ice) THEN
            tau_abs = (1.0_wp - concsum(jc,blockNo))**2 &
                 &    * SQRT((atmos_fluxes%stress_xw(jc,blockNo))**2 &
                 &    + (atmos_fluxes%stress_yw(jc,blockNo))**2 )
          ELSE
            tau_abs = (1.0_wp - concsum(jc,blockNo)) &
                 &    * SQRT((atmos_fluxes%stress_xw(jc,blockNo))**2 &
                 &    + (atmos_fluxes%stress_yw(jc,blockNo))**2 )
          ENDIF

          forc_tke_surf_2D(jc,blockNo) = tau_abs / OceanReferenceDensity

          ! calculate N2

          rho_up(1:levels-1)  = calculate_density_onColumn(&
               & temp(jc,1:levels-1,blockNo), &
               & salt(jc,1:levels-1,blockNo), &
               & pressure(2:levels), levels-1)
          rho_down(2:levels)  = calculate_density_onColumn(&
               & temp(jc,2:levels,blockNo), &
               & salt(jc,2:levels,blockNo), &
               & pressure(2:levels), levels-1)

          DO jk = 2, levels
            Nsqr(jk) = grav/OceanReferenceDensity * (rho_down(jk) - rho_up(jk-1)) *  dzi(jc,jk,blockNo) / s_c(jc,blockNo)
          ENDDO

          ! calculate shear
          DO jk = 2, levels
            Ssqr(jk) = SUM(((  ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x   &
                 &              - ocean_state%p_diag%p_vn(jc,jk,  blockNo)%x ) &
                 &            * dzi(jc,jk,blockNo)/s_c(jc,blockNo) )**2)
          ENDDO


          
      !if (jc==8 .and. blockNo==10) then
      !  write(*,*) 'jc = ', jc, 'blockNo = ', blockNo, 'tstep_count = ', tstep_count
      !  !write(*,*) 'dolic_c(jc,blockNo) = ', patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
      !  write(*,*) 'dolic_c(jc,blockNo) = ', dolic_c(jc,blockNo)
      !  write(*,*) 'temp(jc,:,blockNo) = ', temp(jc,:,blockNo)
      !endif

          IF (l_lc) THEN
            ! calculate Langmuir cell additional term after Axell (2002)

            ! calculate Stoke's drift
            ! Approximation if there is no information about the wave field
            ! As done in Nemo
            ! FIXME: do we need to divide tau by rho?

            ! Option used in NEMO model (https://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf, p.197) see also Breivik et al. (2015)
            ! They assume rhoair=1.2 kg/m3 and cd=1.5e-03:
            ! u_stokes = 0.016/(1.2 * 1.5e-03)^0.5 * |tau|^0.5; although they seem to use rhoair=1.2 kg/m3
            !u_stokes(jc,blockNo) = 0.377_wp * SQRT(tau_abs)                       ! [tau]=N2/m2
            u_stokes(jc,blockNo) = 0.016_wp/SQRT(1.2_wp * 1.5e-03_wp)*SQRT(tau_abs)           ! [tau]=N2/m2, rhoair=1.2, cd=1.5*10e-03

            ! This is done in Coulevard et al (2020, doi:10.5194/gmd-13-3067-2020), see Fig.2
            ! u_stokes(jc,blockNo) = 0.377_wp * SQRT(forc_tke_surf_2D(jc,blockNo))

            ! other option from Li and Garrett (1993)
            !u_stokes(jc,blockNo) = 0.016_wp * fu10(jc,blockNo)

            ! or original version from Axell (2002)
            !LLC = 0.12_wp*(u10**2/g)
            !u_stokes(jc,:,blockNo) = 0.0016*u10*EXP(depth/LLC)


            ! find depth of langmuir cell (hlc). hlc is the depth to which a water
            ! parcel with kinetic energy 0.5*u_stokes**2 can reach on its own by
            ! converting its kinetic energy to potential energy.
            hlc(jc,blockNo) = 0.0_wp
            DO jk = 2, levels
              tmp = SUM( Nsqr(2:jk)*depth_CellInterface(jc,2:jk,blockNo) &  ! see axell (2002) eq.47
                           *prism_thick_c(jc,2:jk,blockNo)*s_c(jc,blockNo) )

              IF(tmp > 0.5_wp*u_stokes(jc,blockNo)**2.0_wp) THEN
                hlc(jc,blockNo) = depth_CellInterface(jc,jk,blockNo)*s_c(jc,blockNo)
                EXIT
              ENDIF
            ENDDO

            ! calculate langmuir cell velocity scale (wlc)
            ! Note: Couvelard et al (2020) set clc=0.3 instead of default 0.15 from
            ! Axell (2002); results in deeper MLDs and better spatial MLD pattern.
            DO jk = 2, levels
!              IF ( depth_CellInterface(jc,jk,blockNo)*s_c(jc,blockNo) <= hlc(jc,blockNo) ) THEN
              IF ( ( depth_CellInterface(jc,jk,blockNo)*s_c(jc,blockNo) <= hlc(jc,blockNo) )  &
                   .AND. ( patch_3D%wet_c(jc,jk,blockNo) .EQ. 1 ) ) THEN
                wlc(jc,jk,blockNo) = clc * u_stokes(jc,blockNo)*SIN(pi*depth_CellInterface(jc,jk,blockNo)/hlc(jc,blockNo))
              ELSE
                wlc(jc,jk,blockNo) = 0.0_wp
              ENDIF
            ENDDO

            ! calculate langmuir turbulence term (tke_plc)
            IF (hlc(jc,blockNo) > 0.0_wp) THEN
              tke_plc(jc,1:levels,blockNo) = wlc(jc,1:levels,blockNo)**3.0_wp / hlc(jc,blockNo)
            ELSE
              tke_plc(jc,1:levels,blockNo) = 0.0_wp
            ENDIF

          ENDIF

      ! main vmix call to calculate tke

          IF (l_lc) THEN
            tke_plc_ptr => tke_plc(jc,:,blockNo)
          ELSE
            tke_plc_ptr => NULL()
          ENDIF

          CALL coeffs_tke(                                                          &
              ! parameter                                                           
              i = jc,                                                               &
              j = blockNo,                                                          &
              tstep_count = tstep_count,                                            &
              tke_old      = tke_old(:),                                            & ! in
              tke_new      = tke(jc,:,blockNo),                                     & ! out
              KappaM_out   = tke_Av(jc,:,blockNo),                                  & ! out
              KappaH_out   = tke_kv(jc,:,blockNo),                                  & ! out
              vmix_int_1  = params_oce%vmix_params%vmix_dummy_1(jc,:,blockNo),      & !
              vmix_int_2  = params_oce%vmix_params%vmix_dummy_2(jc,:,blockNo),      & !
              vmix_int_3  = params_oce%vmix_params%vmix_dummy_3(jc,:,blockNo),      & !
              dzw          = prism_thick_c(jc,:,blockNo)*s_c(jc,blockNo),           &
              dzt          = dz(jc,:,blockNo)*s_c(jc,blockNo),                      &
              nlev         = dolic_c(jc,blockNo),                                   &
              max_nlev     = n_zlev,                                                &
              Ssqr         = Ssqr(:),                                               & ! in
              Nsqr         = Nsqr(:),                                               & ! in
              tke_Tbpr     = params_oce%vmix_params%tke_Tbpr(jc,:,blockNo),         &
              tke_Tspr     = params_oce%vmix_params%tke_Tspr(jc,:,blockNo),         &
              tke_Tdif     = params_oce%vmix_params%tke_Tdif(jc,:,blockNo),         &
              tke_Tdis     = params_oce%vmix_params%tke_Tdis(jc,:,blockNo),         &
              tke_Twin     = params_oce%vmix_params%tke_Twin(jc,:,blockNo),         &
              tke_Tiwf     = params_oce%vmix_params%tke_Tiwf(jc,:,blockNo),         &
              tke_Tbck     = params_oce%vmix_params%tke_Tbck(jc,:,blockNo),         &
              tke_Ttot     = params_oce%vmix_params%tke_Ttot(jc,:,blockNo),         &
              tke_Lmix     = params_oce%vmix_params%tke_Lmix(jc,:,blockNo),         &
              tke_Pr       = params_oce%vmix_params%tke_Pr(jc,:,blockNo),           &
              tke_plc      = tke_plc_ptr,                                           & !by_Oliver
              forc_tke_surf= forc_tke_surf_2D(jc,blockNo),                          &
              E_iw         = tke_iwe(jc,:,blockNo),                                 & ! for IDEMIX Ri
              dtime        = dtime,                                                 &
              bottom_fric  = bottom_fric_2D(jc,blockNo),                            &
              iw_diss      = tke_iwe_forcing(jc,:,blockNo),                         &
              forc_rho_surf= forc_rho_surf_2D(jc,blockNo),                          &
              rho_ref      = OceanReferenceDensity,                                 &
              grav         = grav,                                                  &
              ! essentials
              ! FIXME: nils: better calc IDEMIX Ri directly in ! VMIX/IDEMIX
              alpha_c      = tke_iw_alpha_c(jc,:,blockNo)                           & ! for IDEMIX Ri
              ! forcing
              ! diagnostics
              ! debugging
              !tke          = tke(i,j,:),             &
              !tke_diss_out = tke_diss(i,j,:),        & !
              !KappaM_out   = avo(i,j,:),             & !
              !KappaH_out   = dvo(i,j,:),             & !
              !old_tke_diss = tke_diss(i,j,:),        &
              !Kappa_GM     = Kappa_GM,             &
              !tke_userdef_constants = tke_userdef_constants)
              )
        ENDIF ! dolic > 0


      ENDDO
    ENDDO
!ICON_OMP_END_DO

    IF (.TRUE.) THEN
      ! interpolate vert. visosity from cell center to edges
!ICON_OMP_DO PRIVATE(start_index, end_index, je, levels, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP jk)
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        params_oce%a_veloc_v(:,:,blockNo) = 0.0

        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        DO je = start_index, end_index
          levels       = patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          cell_1_idx   = patch_2D%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
          cell_2_idx   = patch_2D%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
          DO jk = 2, levels
            params_oce%a_veloc_v(je,jk,blockNo) = &
                 & 0.5_wp * (    tke_Av(cell_1_idx,jk,cell_1_block) &
                 &             + tke_Av(cell_2_idx,jk,cell_2_block) )
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (.TRUE.) THEN

      ! write tke vert. diffusivity to vert tracer diffusivities
!ICON_OMP_DO PRIVATE(start_index, end_index, jc)
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          ! FIXME: nils: make loop over all tracer
          DO jtrc = 1, no_tracer
            params_oce%a_tracer_v(jc,:,blockNo,jtrc) = tke_kv(jc,:,blockNo)
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_DO
    ENDIF
!ICON_OMP_END_PARALLEL

    !write(*,*) "Stopping..."
    !stop
    !write(*,*) "TKE after:"
    !write(*,*) tke(8,:,10)

    !params_oce%a_tracer_v = 1e-5
    !write(*,*) 'a_tracer_v = ', params_oce%a_tracer_v(8,:,10,1)
    !write(*,*) 'dolic_c = ', dolic_c(8,10)
    !write(*,*) 'tke_kv = ', tke_kv(8,:,10)
    !write(*,*) 'tke = ', tke(8,:,10)
    !write(*,*) 'fu10 = ', fu10(8,10)
  END SUBROUTINE calc_tke_scalar

#endif

  SUBROUTINE calc_tke_vector(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(IN) :: ocean_state
    TYPE(t_atmos_fluxes), INTENT(IN) :: atmos_fluxes
    TYPE(t_ho_params), INTENT(INOUT) :: params_oce
    REAL(wp), TARGET, INTENT(IN) :: fu10(:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET, INTENT(IN) :: concsum(:,:) !< sea ice concentration
    LOGICAL, INTENT(IN), OPTIONAL        :: lacc

    ! pointer for convenience
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells
    REAL(wp), POINTER :: dz(:,:,:)
    REAL(wp), POINTER :: dzi(:,:,:)
    REAL(wp), POINTER :: dzw(:,:,:)
    !REAL(wp), POINTER :: vvel(:,:,:)
    REAL(wp), POINTER :: temp(:,:,:)
    REAL(wp), POINTER :: salt(:,:,:)
    REAL(wp), POINTER :: Av_old(:,:,:)
    REAL(wp), POINTER :: kv_old(:,:,:)
    REAL(wp), POINTER :: tke(:,:,:)

    ! Langmuir turbulence
    REAL(wp), POINTER, CONTIGUOUS :: tke_plc(:,:,:)
    REAL(wp), POINTER, CONTIGUOUS :: tke_plc_ptr(:,:)

    REAL(wp), POINTER :: wlc(:,:,:)
    REAL(wp), POINTER :: hlc(:,:)
    REAL(wp), POINTER :: u_stokes(:,:)
    REAL(wp), POINTER :: depth_CellInterface(:,:,:)


    INTEGER, POINTER :: dolic_c(:,:)

    ! loop variables
    INTEGER :: jc, blockNo, je, jk, jtrc
    INTEGER :: start_index, end_index
    INTEGER :: max_levels

    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block

    REAL(wp) :: rho_up, rho_down
    REAL(wp) :: pressure(2:n_zlev)
    REAL(wp) :: Nsqr(nproma, n_zlev+1), Ssqr(nproma, n_zlev+1), tmp(nproma)
    !REAL(wp), POINTER :: vert_density_grad(:,:,:)

    ! Parameters for zstar
    REAL(wp) :: s_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    ! stretching factor
    REAL(wp) :: e_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    ! surface height

    REAL(wp) :: dzw_stretched(nproma, n_zlev)
    REAL(wp) :: dzt_stretched(nproma, n_zlev+1)

    LOGICAL :: active(nproma)
    LOGICAL :: lzacc

    INTEGER :: tstep_count

    ! put this later to a global place
    REAL(wp) :: tke_Av(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_kv(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_old(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tbpr(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tspr(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tdif(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tdis(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Twin(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tiwf(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Tbck(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Ttot(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Lmix(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: tke_Pr(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: dummy_1(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: dummy_2(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: dummy_3(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iw_alpha_c(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iwe(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tke_iwe_forcing(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: forc_tke_surf_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: forc_rho_surf_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: bottom_fric_2D(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: tau_abs

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish('calc_tke_vector', 'OpenACC version currently not implemented')
#endif

    ! Suppress Unused-argument warning
    IF (SIZE(fu10) < 0) THEN; END IF

    !tke = 0.0
    tke => params_oce%vmix_params%tke(:,:,:)
    tke_old = tke
    !tke_Tbpr = 0.0
    !tke_Tspr = 0.0
    !tke_Tdif = 0.0
    !tke_Tdis = 0.0
    !tke_Twin = 0.0
    !tke_Tiwf = 0.0
    !tke_Tbck = 0.0
    !tke_Ttot = 0.0
    !tke_Lmix = 0.0
    !tke_Pr = 0.0
    !dummy_1 = 0.0
    !dummy_2 = 0.0
    !dummy_3 = 0.0
    tke_kv = 0.0
    tke_Av = 0.0
    tke_iw_alpha_c = 0.0
    tke_iwe = 0.0
    tke_iwe_forcing = 0.0
    forc_rho_surf_2D = 0.0
    bottom_fric_2D = 0.0

    dzw_stretched = 0.0
    dzt_stretched = 0.0

    IF(l_lc) THEN
      ! Langmuir turbulence variables
      tke_plc  => params_oce%vmix_params%tke_plc(:,:,:)
      hlc      => params_oce%vmix_params%hlc(:,:)
      wlc      => params_oce%vmix_params%wlc(:,:,:)
      u_stokes => params_oce%vmix_params%u_stokes(:,:)
      depth_CellInterface => patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,:)  !interface depths(jc,1:levels,blockNo)
    ELSE
      NULLIFY(tke_plc, hlc, wlc, u_stokes, depth_CellInterface)
    ENDIF

    dz  => patch_3d%p_patch_1d(1)%prism_center_dist_c
    dzi => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
    dzw => patch_3d%p_patch_1d(1)%prism_thick_c
    !uvel => ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x
    temp => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    ! FIXME: use sal_ref in case of temp is only tracer
    salt => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)

    Av_old => params_oce%a_veloc_v(:,:,:)
    ! Use a_tracer_v of temperature here
    kv_old => params_oce%a_tracer_v(:,:,:,1)
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    ! renaming stuff
    patch_2D   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    all_cells  => patch_2D%cells%ALL

    max_levels = n_zlev

    ! special settings if IDEMIX is used together with TKE
    IF ( vert_mix_type==vmix_idemix_tke ) THEN
      ! use iwe dissipation as forcing for tke
      tke_iwe_forcing(:,:,:) = -1.0_wp * params_oce%vmix_params%iwe_Tdis(:,:,:)
    ENDIF

    ! set zstar related parameters
    IF (vert_cor_type == 1) THEN
      s_c(:,:) = ocean_state%p_prog(nold(1))%stretch_c(:,:)
      e_c(:,:) = ocean_state%p_prog(nold(1))%eta_c(:,:)
    ELSE
      s_c(:,:) = 1.0_wp
      e_c(:,:) = 0.0_wp
    ENDIF

    pressure(2:n_zlev) = patch_3d%p_patch_1d(1)%zlev_i(2:n_zlev) * ReferencePressureIndbars

    !write(*,*) "TKE before:"
    !write(*,*) tke(8,:,10)
!ICON_OMP_PARALLEL PRIVATE(rho_up, rho_down)

!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, max_levels, &
!ICON_OMP  tau_abs, Nsqr, Ssqr, tstep_count, tmp, blockNo, tke_plc_ptr, active)

 ! error?! !ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk,                     &
 ! error?! !ICON_OMP  tau_abs, Nsqr, Ssqr, tstep_count, blockNo, forc_tke_surf_2D,  &
 ! error?! !ICON_OMP  forc_rho_surf_2D, u_stokes, bottom_fric_2D, hlc, tmp, wlc,    &
 ! error?! !ICON_OMP  Av_old, kv_old, tke_plc, tke, tke_Av, tke_kv, tke_iwe,        &
 ! error?! !ICON_OMP  tke_iwe_forcing, tke_iw_alpha_c)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      max_levels = MAXVAL(dolic_c(start_index:end_index, blockNo))

      Nsqr = 0.
      Ssqr = 0.

      DO jc = start_index, end_index
        ! wind stress for tke surface forcing
        tau_abs = (1._wp - concsum(jc,blockNo)) &
            &     * SQRT((atmos_fluxes%stress_xw(jc,blockNo))**2 &
            &     + (atmos_fluxes%stress_yw(jc,blockNo))**2 )
        forc_tke_surf_2D(jc,blockNo) = tau_abs / OceanReferenceDensity

        IF (l_lc) THEN
          ! calculate Langmuir cell additional term after Axell (2002)

          ! calculate Stoke's drift
          ! Approximation if there is no information about the wave field
          ! As done in Nemo
          ! FIXME: do we need to divide tau by rho?

          ! Option used in NEMO model (https://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf, p.197) see also Breivik et al. (2015)
          ! They assume rhoair=1.2 kg/m3 and cd=1.5e-03:
          ! u_stokes = 0.016/(1.2 * 1.5e-03)^0.5 * |tau|^0.5; although they seem to use rhoair=1.2 kg/m3
          !u_stokes(jc,blockNo) = 0.377_wp * SQRT(tau_abs)                       ! [tau]=N2/m2
          u_stokes(jc,blockNo) = 0.016_wp/SQRT(1.2_wp * 1.5e-03_wp)*SQRT(tau_abs)           ! [tau]=N2/m2, rhoair=1.2, cd=1.5*10e-03

          ! This is done in Coulevard et al (2020, doi:10.5194/gmd-13-3067-2020), see Fig.2
          ! u_stokes(jc,blockNo) = 0.377_wp * SQRT(forc_tke_surf_2D(jc,blockNo))

          ! other option from Li and Garrett (1993)
          !u_stokes(jc,blockNo) = 0.016_wp * fu10(jc,blockNo)

          ! or original version from Axell (2002)
          !LLC = 0.12_wp*(u10**2/g)
          !u_stokes(jc,:,blockNo) = 0.0016*u10*EXP(depth/LLC)
        END IF
      END DO

      ! Loop over internal interfaces, surface (jk=1) and bottom (jk=kbot+1) excluded
      DO jk = 2, max_levels
        !NEC$ ivdep
        DO jc = start_index, end_index
          IF (jk <= dolic_c(jc,blockNo)) THEN
            ! calculate N2
            !NEC$ inline_complete
            rho_up = calculate_density_onColumn_elem(&
                & temp(jc,jk-1,blockNo), &
                & salt(jc,jk-1,blockNo), &
                & pressure(jk))
            !NEC$ inline_complete
            rho_down = calculate_density_onColumn_elem(&
                & temp(jc,jk,blockNo), &
                & salt(jc,jk,blockNo), &
                & pressure(jk))

            Nsqr(jc,jk) = grav/OceanReferenceDensity * (rho_down - rho_up) *  dzi(jc,jk,blockNo) / s_c(jc,blockNo)

            ! calculate shear
            Ssqr(jc,jk) = SUM(((  ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x   &
              &              - ocean_state%p_diag%p_vn(jc,jk,  blockNo)%x ) &
              &            * dzi(jc,jk,blockNo)/s_c(jc,blockNo) )**2)
          END IF
        END DO
      END DO

      DO jk = 1, max_levels
        DO jc = start_index, end_index
          dzw_stretched(jc,jk) = dzw(jc,jk,blockNo)* s_c(jc,blockNo)
        END DO
      END DO
      DO jk = 1, max_levels + 1
        DO jc = start_index, end_index
          dzt_stretched(jc,jk) = dz(jc,jk,blockNo) * s_c(jc,blockNo)
        END DO
      END DO

      !if (jc==8 .and. blockNo==10) then
      !  write(*,*) 'jc = ', jc, 'blockNo = ', blockNo, 'tstep_count = ', tstep_count
      !  !write(*,*) 'kbot(jc,blockNo) = ', patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
      !  write(*,*) 'kbot(jc,blockNo) = ', kbot(jc,blockNo)
      !  write(*,*) 'temp(jc,:,blockNo) = ', temp(jc,:,blockNo)
      !endif
      IF (l_lc) THEN

        ! find depth of langmuir cell (hlc). hlc is the depth to which a water
        ! parcel with kinetic energy 0.5*u_stokes**2 can reach on its own by
        ! converting its kinetic energy to potential energy.
        DO jc = start_index, end_index
          hlc(jc,blockNo) = 0.0_wp
          tmp(jc) = 0.0_wp
          active(jc) = .TRUE.
        END DO

        DO jk = 2, max_levels
          DO jc = start_index, end_index
            IF(active(jc) .AND. jk <= dolic_c(jc,blockNo)) THEN
              tmp(jc) = tmp(jc) + Nsqr(jc,jk+1)*(depth_CellInterface(jc,jk+1,blockNo)*s_c(jc,blockNo) - e_c(jc,blockNo))

              IF(tmp(jc) > 0.5_wp*u_stokes(jc,blockNo)**2.0_wp) THEN
                hlc(jc,blockNo) = depth_CellInterface(jc,jk,blockNo)*s_c(jc,blockNo) - e_c(jc,blockNo)
                active(jc) = .FALSE.
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        ! calculate langmuir cell velocity scale (wlc)
        ! Note: Couvelard et al (2020) set clc=0.3 instead of default 0.15 from
        ! Axell (2002); results in deeper MLDs and better spatial MLD pattern.
        DO jk = 2, max_levels
          DO jc = start_index, end_index
            IF ( ( depth_CellInterface(jc,jk,blockNo)*s_c(jc,blockNo) - e_c(jc,blockNo) <= hlc(jc,blockNo) )  &
                .AND. ( patch_3D%wet_c(jc,jk,blockNo) .EQ. 1 ) ) then
              wlc(jc,jk,blockNo) = clc * u_stokes(jc,blockNo)*SIN(pi*depth_CellInterface(jc,jk,blockNo)/hlc(jc,blockNo))
            ELSE
              wlc(jc,jk,blockNo) = 0.0_wp
            ENDIF
          ENDDO
        ENDDO

        ! calculate langmuir turbulence term (tke_plc)
        DO jk = 1, max_levels + 1
          DO jc = start_index, end_index
            IF (hlc(jc,blockNo) > 0.0_wp) THEN
              tke_plc(jc,jk,blockNo) = wlc(jc,jk,blockNo)**3.0_wp / hlc(jc,blockNo)
            ELSE
              tke_plc(jc,jk,blockNo) = 0.0_wp
            ENDIF
          ENDDO
        ENDDO
      END IF

      ! main call to calculate tke
      IF (max_levels > 0) THEN

        IF (l_lc) THEN
          tke_plc_ptr => tke_plc(:,:,blockNo)
        ELSE
          tke_plc_ptr => NULL()
        END IF

        CALL coeffs_tke(                                                            &
            ! parameter
            nproma = nproma,                                                        &
            si = start_index,                                                       &
            ei = end_index,                                                         &
            j = blockNo,                                                            &
            tstep_count = tstep_count,                                              &
            tke_old      = tke_old(:,:,blockNo),                                    & ! in
            tke_new      = tke(:,:,blockNo),                                        & ! out
            KappaM_out   = tke_Av(:,:,blockNo),                                     & ! out
            KappaH_out   = tke_kv(:,:,blockNo),                                     & ! out
            vmix_int_1  = params_oce%vmix_params%vmix_dummy_1(:,:,blockNo),         & !
            vmix_int_2  = params_oce%vmix_params%vmix_dummy_2(:,:,blockNo),         & !
            vmix_int_3  = params_oce%vmix_params%vmix_dummy_3(:,:,blockNo),         & !
            dzw          = dzw_stretched(:,:),                                      &
            dzt          = dzt_stretched(:,:),                                      &
            nlev         = dolic_c(:,blockNo),                                      &
            max_nlev     = n_zlev,                                                  &
            Ssqr         = Ssqr(:,:),                                               & ! in
            Nsqr         = Nsqr(:,:),                                               & ! in
            tke_Tbpr     = params_oce%vmix_params%tke_Tbpr(:,:,blockNo),            &
            tke_Tspr     = params_oce%vmix_params%tke_Tspr(:,:,blockNo),            &
            tke_Tdif     = params_oce%vmix_params%tke_Tdif(:,:,blockNo),            &
            tke_Tdis     = params_oce%vmix_params%tke_Tdis(:,:,blockNo),            &
            tke_Twin     = params_oce%vmix_params%tke_Twin(:,:,blockNo),            &
            tke_Tiwf     = params_oce%vmix_params%tke_Tiwf(:,:,blockNo),            &
            tke_Tbck     = params_oce%vmix_params%tke_Tbck(:,:,blockNo),            &
            tke_Ttot     = params_oce%vmix_params%tke_Ttot(:,:,blockNo),            &
            tke_Lmix     = params_oce%vmix_params%tke_Lmix(:,:,blockNo),            &
            tke_Pr       = params_oce%vmix_params%tke_Pr(:,:,blockNo),              &
            tke_plc      = tke_plc_ptr,                                             & !by_Oliver, NULL = not present
            forc_tke_surf= forc_tke_surf_2D(:,blockNo),                             &
            E_iw         = tke_iwe(:,:,blockNo),                                    & ! for IDEMIX Ri
            dtime        = dtime,                                                   &
            bottom_fric  = bottom_fric_2D(:,blockNo),                               &
            old_kappaM   = Av_old(:,:,blockNo),                                     & ! in
            old_KappaH   = kv_old(:,:,blockNo),                                     & ! in
            iw_diss      = tke_iwe_forcing(:,:,blockNo),                            &
            forc_rho_surf= forc_rho_surf_2D(:,blockNo),                             &
            rho_ref      = OceanReferenceDensity,                                   &
            grav         = grav,                                                    &
            ! essentials
            ! FIXME: nils: better calc IDEMIX Ri directly in ! VMIX/IDEMIX
            alpha_c      = tke_iw_alpha_c(:,:,blockNo)                              & ! for IDEMIX Ri
            ! forcing
            ! diagnostics
            ! debugging
            !tke          = tke(i,j,:),             &
            !tke_diss_out = tke_diss(i,j,:),        & !
            !KappaM_out   = avo(i,j,:),             & !
            !KappaH_out   = dvo(i,j,:),             & !
            !old_tke_diss = tke_diss(i,j,:),        &
            !Kappa_GM     = Kappa_GM,             &
            !tke_userdef_constants = tke_userdef_constants)
          )
      END IF
    ENDDO
!ICON_OMP_END_DO

    if (.true.) then
    ! interpolate vert. visosity from cell center to edges
!ICON_OMP_DO PRIVATE(start_index, end_index, je, max_levels, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP jk)
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      params_oce%a_veloc_v(:,:,blockNo) = 0.0

      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      max_levels = MAXVAL(patch_3d%p_patch_1d(1)%dolic_e(start_index:end_index, blockNo))

      DO jk = 2, max_levels
        DO je = start_index, end_index
          IF (jk <= patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)) THEN
            cell_1_idx   = patch_2D%edges%cell_idx(je,blockNo,1)
            cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
            cell_2_idx   = patch_2D%edges%cell_idx(je,blockNo,2)
            cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
            params_oce%a_veloc_v(je,jk,blockNo) = &
              & 0.5_wp * (    tke_Av(cell_1_idx,jk,cell_1_block) &
              &             + tke_Av(cell_2_idx,jk,cell_2_block) )
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    end if

    if (.true.) then

    ! write tke vert. diffusivity to vert tracer diffusivities
!ICON_OMP_DO PRIVATE(start_index, end_index, jc)
   DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jtrc = 1, no_tracer
        params_oce%a_tracer_v(start_index:end_index,:,blockNo,jtrc) = tke_kv(start_index:end_index,:,blockNo)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
    end if
!ICON_OMP_END_PARALLEL

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
  END SUBROUTINE calc_tke_vector

END MODULE mo_ocean_tke
