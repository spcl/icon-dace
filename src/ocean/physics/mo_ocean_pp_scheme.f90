! Provide an implementation of the ocean physics.
!
! Provide an implementation of the physical parameters and characteristics
! for the hydrostatic ocean model.
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

!----------------------------
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
MODULE mo_ocean_pp_scheme
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: &
    & n_zlev,                             &
!     & HarmonicViscosity_reference, velocity_VerticalDiffusion_background,                 &
    & Temperature_VerticalDiffusion_background, Salinity_VerticalDiffusion_background, no_tracer,  &
    & tracer_convection_MixingCoefficient,                                     &
!     & BiharmonicViscosity_scaling, HarmonicViscosity_scaling, &
!     & VelocityDiffusion_order,                                &
!     & BiharmonicViscosity_reference,                          &
    & tracer_RichardsonCoeff, velocity_RichardsonCoeff,        &
    & PPscheme_type,                             &
    & PPscheme_Constant_type,                    &
    & PPscheme_ICON_Edge_type,                   &
    & PPscheme_ICON_Edge_vnPredict_type,         &
    & use_wind_mixing,                                        &
!     & HorizontalViscosity_SmoothIterations,                   &
    & convection_InstabilityThreshold,                        &
    & RichardsonDiffusion_threshold,                          &
    & lambda_wind,                                            &
    & use_reduced_mixing_under_ice,                           &
    & k_tracer_dianeutral_parameter,                          &
    & k_tracer_isoneutral_parameter, k_tracer_GM_kappa_parameter,    &
    & GMRedi_configuration,GMRedi_combined,                   &
    & GM_only,Redi_only,  Cartesian_Mixing,                   &
!     & laplacian_form,                                         &
!     & HorizontalViscosity_SpatialSmoothFactor,                &
    & VerticalViscosity_TimeWeight, OceanReferenceDensity,    &
    & tracer_TopWindMixing, WindMixingDecayDepth,             &
    & velocity_TopWindMixing, &
!     &  Temperature_HorizontalDiffusion_Background,            &
!     &  Temperature_HorizontalDiffusion_Reference,             &
!     &  Salinity_HorizontalDiffusion_Background,               &
!     &  Salinity_HorizontalDiffusion_Reference,                &
    &  ReferencePressureIndbars

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
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = 'mo_pp_scheme'
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug

  PUBLIC :: update_PP_scheme
  PUBLIC :: update_PP_scheme_zstar
  PUBLIC :: ICON_PP_Edge_vnPredict_scheme

  REAL(wp), POINTER :: WindAmplitude_at10m(:,:)  
  REAL(wp), POINTER :: SeaIceConcentration(:,:)
  PUBLIC :: calculate_rho4GMRedi
  

CONTAINS

 !-------------------------------------------------------------------------
  !>
  !! 
  !!
  !<Optimize:inUse:done>
  SUBROUTINE calculate_rho4GMRedi(patch_3d, temperature, salinity, rho_GM) 

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    !TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), INTENT(IN)                 :: temperature(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN)                 :: salinity(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(INOUT)              :: rho_GM(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !TYPE(t_ho_params), INTENT(inout)     :: params_oce

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), sal(n_zlev)
    !REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    !edges_in_domain => patch_2D%edges%in_domain
    !cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    !z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    levels = n_zlev


!     IF (ltimer) CALL timer_start(timer_extra10)

!ICON_OMP_PARALLEL PRIVATE(salinity,z_rho_up, z_rho_down)
    !sal(1:levels) = sal_ref

    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp
    rho_GM(:,:,:)=0.0_wp

!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, &
!ICON_OMP  tracer_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (levels < 2) CYCLE

        !--------------------------------------------------------
        sal(1:levels)      = salinity(jc,1:levels,blockNo)

        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * OceanReferenceDensity * sitodbar

        z_rho_up(1:levels-1)  &
        &= calculate_density_onColumn(temperature(jc,1:levels-1,blockNo),&
                                    & salinity   (jc,1:levels-1,blockNo), pressure(2:levels), levels-1)
        z_rho_down(2:levels)  &
        &= calculate_density_onColumn(temperature(jc, 2:levels, blockNo),&
                                    & salinity   (jc, 2:levels, blockNo), pressure(2:levels), levels-1)


        DO jk = 2, levels
          rho_GM(jc,jk,blockNo)=0.5_wp*(z_rho_up(jk)+z_rho_down(jk))
        
          !z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
          !  & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
        END DO ! levels
        !ocean_state%p_diag%grad_rho_PP_vert(jc,2:levels,blockNo)=z_vert_density_grad_c(jc,2:levels,blockNo)
        rho_GM(jc,1,blockNo)=rho_GM(jc,2,blockNo)
 
      END DO ! index


    END DO ! blocks
!ICON_OMP_END_DO

!ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra10)
!     IF (ltimer) CALL timer_start(timer_extra11)
! !ICON_OMP_PARALLEL

!     IF (ltimer) CALL timer_stop(timer_extra11)
    DO jk=1,10
      CALL dbg_print('calc_rho4GMRedi: rho_GM',rho_GM(:,jk,:),&
        & module_name, idt_src, in_subset=all_cells) 
    END DO          

  END SUBROUTINE calculate_rho4GMRedi
  !-------------------------------------------------------------------------


!<Optimize:inUse:done>
  SUBROUTINE update_PP_scheme(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum(:,:) ! t_sea_ice%concsum
    TYPE(t_ho_params), INTENT(inout)     :: params_oce
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs

    INTEGER :: tracer_index
    !-------------------------------------------------------------------------
    WindAmplitude_at10m => fu10
    SeaIceConcentration => concsum

    SELECT CASE (PPscheme_type)
    CASE (PPscheme_Constant_type)
      !nothing to do!In sbr init_ho_params (see above)
      !tracer mixing coefficient params_oce%A_tracer_v(:,:,:, tracer_index) is already
      !initialzed with params_oce%A_tracer_v_back(tracer_index)
      !and velocity diffusion coefficient
      RETURN


    CASE (PPscheme_ICON_Edge_type)
      CALL ICON_PP_Edge_scheme2(patch_3d, ocean_state, params_oce)

    CASE (PPscheme_ICON_Edge_vnPredict_type)
!       CALL update_PhysicsParameters_ICON_PP_Tracer(patch_3d, ocean_state)
      CALL ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce)
      ! the velovity friction will be updated during dynamics

    CASE default
      CALL finish("update_ho_params", "unknown PPscheme_type")
    END SELECT

  END SUBROUTINE update_PP_scheme
  !-------------------------------------------------------------------------

!<Optimize:inUse:done>
  SUBROUTINE update_PP_scheme_zstar(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs, &
      & eta_c, stretch_c, stretch_e) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum(:,:) ! t_sea_ice%concsum
    TYPE(t_ho_params), INTENT(inout)     :: params_oce
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    REAL(wp), INTENT(IN) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 

    INTEGER :: tracer_index
 
    !-------------------------------------------------------------------------
    WindAmplitude_at10m => fu10
    SeaIceConcentration => concsum

    SELECT CASE (PPscheme_type)
    CASE (PPscheme_Constant_type)
      !nothing to do!In sbr init_ho_params (see above)
      !tracer mixing coefficient params_oce%A_tracer_v(:,:,:, tracer_index) is already
      !initialzed with params_oce%A_tracer_v_back(tracer_index)
      !and velocity diffusion coefficient
      RETURN

    !! FIXME: scheme2 has no zstar port
    CASE (PPscheme_ICON_Edge_type)
      CALL ICON_PP_Edge_scheme_zstar(patch_3d, ocean_state, params_oce, &
        & eta_c, stretch_c, stretch_e)


    CASE (PPscheme_ICON_Edge_vnPredict_type)
      CALL ICON_PP_Edge_scheme_zstar(patch_3d, ocean_state, params_oce, &
        & eta_c, stretch_c, stretch_e)

    CASE default
      CALL finish("update_ho_params", "unknown PPscheme_type")
    END SELECT
    !-------------------------------------------------------------------------


  END SUBROUTINE update_PP_scheme_zstar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! As in the ICON_PP_scheme, but
  !! velocity gradients for the vertical viscocity are clculated on edges
  !!
  !<Optimize:inUse:done>
  SUBROUTINE ICON_PP_Edge_scheme2(patch_3d, ocean_state, params_oce) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_ho_params), INTENT(inout)            :: params_oce

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp), POINTER :: z_ri_cell(:,:,:)
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)

    !Below is a set of variables and parameters for tracer and velocity
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
    REAL(wp) :: diffusion_weight
    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge
    REAL(wp), POINTER :: tracer_windMixing(:,:), velocity_windMixing(:,:)
    REAL(wp) :: z_vert_density_grad_e(n_zlev+1)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    z_ri_cell =>  ocean_state%p_diag%Richardson_Number
    levels = n_zlev

    !-------------------------------------------------------------------------
    z_grav_rho                   = grav/OceanReferenceDensity
    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
    !-------------------------------------------------------------------------
!     IF (ltimer) CALL timer_start(timer_extra10)
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, &
!ICON_OMP tracer_index, diffusion_weight, tracer_windMixing) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      tracer_windMixing => params_oce%tracer_windMixing(:,:,blockNo)

      IF (use_wind_mixing) THEN
        DO jc = start_index, end_index

          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          ! wind-mixing: Marsland et al., 2003
          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
          tracer_windMixing(jc,1) = tracer_TopWindMixing  * WindAmplitude_at10m(jc,blockNo)**3 * &
            & (1.0_wp - SeaIceConcentration(jc,blockNo))

          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
          DO jk = 2, levels
           tracer_windMixing(jc,jk) =  tracer_windMixing(jc,jk-1) * WindMixingDecay(jk) * &
             & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_c(jc,jk,blockNo),0.0_wp))

          END DO! levels

        END DO ! index

      END IF  ! use_wind_mixing
      !-----------------------------------------------------------

      DO tracer_index = 1, no_tracer
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          DO jk = 2, levels
            ! by default use the richardson formula, no convection
            params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = MIN( &
               params_oce%a_tracer_v_back(tracer_index) +   &
               tracer_RichardsonCoeff / ((1.0_wp + z_c1_t * z_ri_cell(jc, jk,blockNo))**3) &
               + tracer_windMixing(jc,jk), &
               tracer_convection_MixingCoefficient)

            IF (z_vert_density_grad_c(jc,jk,blockNo) <= convection_InstabilityThreshold) THEN
              ! convection
              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = tracer_convection_MixingCoefficient
            ELSE
              IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold) THEN
                ! interpolate between convection and richardson diffusion
                diffusion_weight =  &
                  & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
                  & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
                params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
                  & tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
                  & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
              ENDIF
            ENDIF
            
          ENDDO ! levels
        ENDDO !  block index
      ENDDO ! tracer_index]

    END DO ! blocks
!ICON_OMP_END_DO

! !ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra10)
!     IF (ltimer) CALL timer_start(timer_extra11)
! !ICON_OMP_PARALLEL
    !--------------------------------------------
    ! Calculate params_oce%A_veloc_v:
!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP jk, dz, density_differ_edge, z_shear_edge, richardson_edge,z_vert_density_grad_e,velocity_windMixing) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      velocity_windMixing => params_oce%velocity_windMixing(:,:,blockNo)

      DO je = start_index, end_index

        cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
        cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          z_vert_density_grad_e(jk) =  0.5_wp * &
            & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
            &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
        ENDDO

        IF (use_wind_mixing) THEN
          velocity_windMixing(je, 1) =           &
            & velocity_TopWindMixing &
            & *  (0.5_wp * &
            &    (WindAmplitude_at10m(cell_1_idx,cell_1_block) + WindAmplitude_at10m(cell_2_idx,cell_2_block)))**3 &
            & * (1.0_wp - 0.5_wp *       &
            &    (SeaIceConcentration(cell_1_idx,cell_1_block) + SeaIceConcentration(cell_2_idx,cell_2_block)))

            ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
            DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
              velocity_windMixing(je, jk) =  velocity_windMixing(je, jk-1) * WindMixingDecay(jk) * &
                & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_e(jk),0.0_wp))
            END DO! levels
        END IF  ! use_wind_mixing

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          ! TODO: the following expect equally sized cells
          ! compute density gradient at edges
          dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
            &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo))
          density_differ_edge = z_vert_density_grad_e(jk) * dz
          z_shear_edge = dbl_eps + &
            & (ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
            &  ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)   )**2

          richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)

          params_oce%a_veloc_v(je,jk,blockNo) =                 &
            & params_oce%a_veloc_v_back * dz +             &
            & velocity_RichardsonCoeff /                           &
            & ((1.0_wp + z_c1_v * richardson_edge)**2) +   &
            & velocity_windMixing(je, jk)

        END DO ! jk = 2, levels
      ENDDO ! je = start_index, end_index
    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra11)

  END SUBROUTINE ICON_PP_Edge_scheme2
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! As in the ICON_PP_scheme, but
  !! velocity gradients for the vertical viscocity are clculated on edges
  !!
  !<Optimize:inUse:done>
  SUBROUTINE ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_ho_params), INTENT(inout)            :: params_oce

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    REAL(wp) :: z_shear_cell
    REAL(wp) :: z_ri_cell(nproma, n_zlev)
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)

    !Below is a set of variables and parameters for tracer and velocity
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
    REAL(wp) :: diffusion_weight
    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge
    REAL(wp), POINTER :: tracer_windMixing(:,:), velocity_windMixing(:,:)
    REAL(wp) :: z_vert_density_grad_e(n_zlev+1)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    levels = n_zlev

    !-------------------------------------------------------------------------
    z_grav_rho                   = grav/OceanReferenceDensity
    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
    !-------------------------------------------------------------------------
!     IF (ltimer) CALL timer_start(timer_extra10)
!ICON_OMP_PARALLEL PRIVATE(salinity, z_rho_up, z_rho_down, pressure, z_ri_cell, &
!ICON_OMP tracer_windMixing, z_vert_density_grad_e,velocity_windMixing)
    salinity(1:levels) = sal_ref
    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp
    pressure(:) = 0._wp
    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, &
!ICON_OMP z_shear_cell, tracer_index, diffusion_weight) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      z_ri_cell(:,:) = 0.0_wp
      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (levels < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        ! pressure in dbars
        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * ReferencePressureIndbars
        z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
          & salinity(1:levels-1), pressure(2:levels), levels-1)
        z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
          & salinity(2:levels), pressure(2:levels), levels-1)

        IF(GMRedi_configuration/=Cartesian_Mixing) THEN
          ocean_state%p_diag%rho_GM(jc,2:levels-1,blockNo)=0.5_wp*(z_rho_up(2:levels-1)+z_rho_down(2:levels-1))
          ocean_state%p_diag%rho_GM(jc,1,blockNo)=0.5_wp*(z_rho_up(2)+z_rho_down(2))
          ocean_state%p_diag%rho_GM(jc,levels:n_zlev,blockNo)=ocean_state%p_diag%rho_GM(jc,levels-1,blockNo)
        ENDIF
        
        DO jk = 2, levels
                    
          z_shear_cell = dbl_eps + &
            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
          z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
          z_ri_cell(jc, jk) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
            & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
                                                                     ! this is canceled out in this formula
        END DO ! levels
        ocean_state%p_diag%grad_rho_PP_vert(jc,2:levels,blockNo)=z_vert_density_grad_c(jc,2:levels,blockNo)
      END DO ! index
      !-----------------------------------------------------------
      tracer_windMixing => params_oce%tracer_windMixing(:,:,blockNo)
!       tracer_windMixing(:,:) = 0.0_wp
      IF (use_wind_mixing) THEN
        DO jc = start_index, end_index

          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          ! wind-mixing: Marsland et al., 2003
          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
          tracer_windMixing(jc,1) = tracer_TopWindMixing  * WindAmplitude_at10m(jc,blockNo)**3 * &
            & (1.0_wp - SeaIceConcentration(jc,blockNo))

          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
          DO jk = 2, levels
           tracer_windMixing(jc,jk) =  tracer_windMixing(jc,jk-1) * WindMixingDecay(jk) * &
             & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_c(jc,jk,blockNo),0.0_wp))

          END DO! levels

        END DO ! index

      END IF  ! use_wind_mixing
      !-----------------------------------------------------------

      DO tracer_index = 1, no_tracer
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          DO jk = 2, levels
            ! by default use the richardson formula, no convection
            params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = MIN( &
              & params_oce%a_tracer_v_back(tracer_index) +   &
              & tracer_RichardsonCoeff / ((1.0_wp + z_c1_t * z_ri_cell(jc, jk))**3) + tracer_windMixing(jc,jk), &
              & tracer_convection_MixingCoefficient)

            IF (z_vert_density_grad_c(jc,jk,blockNo) <= convection_InstabilityThreshold) THEN
              ! convection
              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = tracer_convection_MixingCoefficient
            ELSE
              IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold) THEN
                ! interpolate between convection and richardson diffusion
                diffusion_weight =  &
                  & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
                  & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
                params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
                  & tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
                  & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
              ENDIF
            ENDIF
            
          ENDDO ! levels
        ENDDO !  block index
      ENDDO ! tracer_index]

    END DO ! blocks
!ICON_OMP_END_DO

! !ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra10)
!     IF (ltimer) CALL timer_start(timer_extra11)
! !ICON_OMP_PARALLEL
    !--------------------------------------------
    ! Calculate params_oce%A_veloc_v:
!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP jk, dz, density_differ_edge, z_shear_edge, richardson_edge) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      velocity_windMixing => params_oce%velocity_windMixing(:,:,blockNo)

      DO je = start_index, end_index

        cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
        cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          z_vert_density_grad_e(jk) =  0.5_wp * &
            & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
            &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
        ENDDO

        IF (use_wind_mixing) THEN
          velocity_windMixing(je, 1) =           &
            & velocity_TopWindMixing &
            & *  (0.5_wp * &
            &    (WindAmplitude_at10m(cell_1_idx,cell_1_block) + WindAmplitude_at10m(cell_2_idx,cell_2_block)))**3 &
            & * (1.0_wp - 0.5_wp *       &
            &    (SeaIceConcentration(cell_1_idx,cell_1_block) + SeaIceConcentration(cell_2_idx,cell_2_block)))

            ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
            DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
              velocity_windMixing(je, jk) =  velocity_windMixing(je, jk-1) * WindMixingDecay(jk) * &
                & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_e(jk),0.0_wp))
            END DO! levels
        END IF  ! use_wind_mixing

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          ! TODO: the following expect equally sized cells
          ! compute density gradient at edges
          dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
            &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo))
          density_differ_edge = z_vert_density_grad_e(jk) * dz
          z_shear_edge = dbl_eps + &
            & (ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
            &  ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)   )**2

          richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)

          params_oce%a_veloc_v(je,jk,blockNo) =                 &
            & params_oce%a_veloc_v_back * dz +             &
            & velocity_RichardsonCoeff /                           &
            & ((1.0_wp + z_c1_v * richardson_edge)**2) +   &
            & velocity_windMixing(je, jk)

        END DO ! jk = 2, levels
      ENDDO ! je = start_index, end_index
    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra11)

  END SUBROUTINE ICON_PP_Edge_scheme
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
!! As in the ICON_PP_scheme, but
  !! velocity gradients for the vertical viscocity are clculated on edges
  !!
  !<Optimize:inUse:done>
  SUBROUTINE ICON_PP_Edge_scheme_zstar(patch_3d, ocean_state, params_oce, &
      & eta_c, stretch_c, stretch_e) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_ho_params), INTENT(inout)            :: params_oce
    REAL(wp), INTENT(IN) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    REAL(wp) :: z_shear_cell
    REAL(wp) :: z_ri_cell(nproma, n_zlev)
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)

    !Below is a set of variables and parameters for tracer and velocity
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
    REAL(wp) :: diffusion_weight
    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge
    REAL(wp), POINTER :: tracer_windMixing(:,:), velocity_windMixing(:,:)
    REAL(wp) :: z_vert_density_grad_e(n_zlev+1)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    levels = n_zlev

    !-------------------------------------------------------------------------
    z_grav_rho                   = grav/OceanReferenceDensity
    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
    !-------------------------------------------------------------------------
!     IF (ltimer) CALL timer_start(timer_extra10)
!ICON_OMP_PARALLEL PRIVATE(salinity, z_rho_up, z_rho_down, pressure, z_ri_cell, &
!ICON_OMP tracer_windMixing, z_vert_density_grad_e,velocity_windMixing)
    salinity(1:levels) = sal_ref
    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp
    pressure(:) = 0._wp
    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, &
!ICON_OMP z_shear_cell, tracer_index, diffusion_weight) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      z_ri_cell(:,:) = 0.0_wp
      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (levels < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        ! pressure in dbars
        pressure(2:levels) = (patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * stretch_c(jc,blockNo) &
          & - eta_c(jc,blockNo)) * ReferencePressureIndbars !+ ocean_state%p_aux%bc_total_top_potential(jc,blockNo)
        z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
          & salinity(1:levels-1), pressure(2:levels), levels-1)
        z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
          & salinity(2:levels), pressure(2:levels), levels-1)

        IF(GMRedi_configuration/=Cartesian_Mixing) THEN
          ocean_state%p_diag%rho_GM(jc,2:levels-1,blockNo)=0.5_wp*(z_rho_up(2:levels-1)+z_rho_down(2:levels-1))
          ocean_state%p_diag%rho_GM(jc,1,blockNo)=0.5_wp*(z_rho_up(2)+z_rho_down(2))
          ocean_state%p_diag%rho_GM(jc,levels:n_zlev,blockNo)=ocean_state%p_diag%rho_GM(jc,levels-1,blockNo)
        ENDIF
        
        DO jk = 2, levels
                    
          z_shear_cell = dbl_eps + &
            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
          z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo) / stretch_c(jc,blockNo)
          z_ri_cell(jc, jk) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * stretch_c(jc,blockNo) * z_grav_rho * &
            & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
                                                                     ! this is canceled out in this formula
        END DO ! levels
        ocean_state%p_diag%grad_rho_PP_vert(jc,2:levels,blockNo)=z_vert_density_grad_c(jc,2:levels,blockNo)
      END DO ! index
      !-----------------------------------------------------------
      tracer_windMixing => params_oce%tracer_windMixing(:,:,blockNo)
!       tracer_windMixing(:,:) = 0.0_wp
      IF (use_wind_mixing) THEN
        DO jc = start_index, end_index

          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          ! wind-mixing: Marsland et al., 2003
          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
          tracer_windMixing(jc,1) = tracer_TopWindMixing  * WindAmplitude_at10m(jc,blockNo)**3 * &
            & (1.0_wp - SeaIceConcentration(jc,blockNo))

          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
          DO jk = 2, levels
           tracer_windMixing(jc,jk) =  tracer_windMixing(jc,jk-1) * WindMixingDecay(jk) * &
             & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_c(jc,jk,blockNo),0.0_wp))

          END DO! levels

        END DO ! index

      END IF  ! use_wind_mixing
      !-----------------------------------------------------------

      DO tracer_index = 1, no_tracer
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          DO jk = 2, levels
            ! by default use the richardson formula, no convection
            params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = MIN( &
              & params_oce%a_tracer_v_back(tracer_index) +   &
              & tracer_RichardsonCoeff / ((1.0_wp + z_c1_t * z_ri_cell(jc, jk))**3) + tracer_windMixing(jc,jk), &
              & tracer_convection_MixingCoefficient)

            IF (z_vert_density_grad_c(jc,jk,blockNo) <= convection_InstabilityThreshold) THEN
              ! convection
              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = tracer_convection_MixingCoefficient
            ELSE
              IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold) THEN
                ! interpolate between convection and richardson diffusion
                diffusion_weight =  &
                  & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
                  & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
                params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
                  & tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
                  & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
              ENDIF
            ENDIF
            
          ENDDO ! levels
        ENDDO !  block index
      ENDDO ! tracer_index]

    END DO ! blocks
!ICON_OMP_END_DO

! !ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra10)
!     IF (ltimer) CALL timer_start(timer_extra11)
! !ICON_OMP_PARALLEL
    !--------------------------------------------
    ! Calculate params_oce%A_veloc_v:
!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP jk, dz, density_differ_edge, z_shear_edge, richardson_edge) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      velocity_windMixing => params_oce%velocity_windMixing(:,:,blockNo)

      DO je = start_index, end_index

        cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
        cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          z_vert_density_grad_e(jk) =  0.5_wp * &
            & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
            &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
        ENDDO

        IF (use_wind_mixing) THEN
          velocity_windMixing(je, 1) =           &
            & velocity_TopWindMixing &
            & *  (0.5_wp * &
            &    (WindAmplitude_at10m(cell_1_idx,cell_1_block) + WindAmplitude_at10m(cell_2_idx,cell_2_block)))**3 &
            & * (1.0_wp - 0.5_wp *       &
            &    (SeaIceConcentration(cell_1_idx,cell_1_block) + SeaIceConcentration(cell_2_idx,cell_2_block)))

            ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
            DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
              velocity_windMixing(je, jk) =  velocity_windMixing(je, jk-1) * WindMixingDecay(jk) * &
                & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_e(jk),0.0_wp))
            END DO! levels
        END IF  ! use_wind_mixing

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          ! TODO: the following expect equally sized cells
          ! compute density gradient at edges
          dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
            &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo)) * stretch_e(je,blockNo)
          density_differ_edge = z_vert_density_grad_e(jk) * dz
          z_shear_edge = dbl_eps + &
            & (ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
            &  ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)   )**2

          richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)

          params_oce%a_veloc_v(je,jk,blockNo) =                 &
            & params_oce%a_veloc_v_back * dz +             &
            & velocity_RichardsonCoeff /                           &
            & ((1.0_wp + z_c1_v * richardson_edge)**2) +   &
            & velocity_windMixing(je, jk)

        END DO ! jk = 2, levels
      ENDDO ! je = start_index, end_index
    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra11)

  END SUBROUTINE ICON_PP_Edge_scheme_zstar
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! This calculate vertical viscocity for the predicted vn
  !! and time-weights it with the previous viscocity
  !!
  !<Optimize:inUse:done>
  SUBROUTINE ICON_PP_Edge_vnPredict_scheme(patch_3d, &
    & blockNo, start_index, end_index, ocean_state, vn_predict, lacc) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    INTEGER, INTENT(in) :: blockNo, start_index, end_index
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    REAL(wp) :: vn_predict(:,:)
    LOGICAL, INTENT(in), OPTIONAL     :: lacc

    ! Local variables
    INTEGER :: je,jk
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: levels
    LOGICAL :: lzacc

    !Below is a set of variables and parameters for tracer and velocity
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp) :: diffusion_weight
    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge, vn_diff, new_velocity_friction
    !-------------------------------------------------------------------------
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
    REAL(wp), POINTER :: wind_mixing(:,:)
    REAL(wp) :: z_vert_density_grad_e(1:n_zlev+1)
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_ho_params), POINTER :: params_oce
    CHARACTER(LEN=*), PARAMETER :: routine='ICON_PP_Edge_vnPredict_scheme'

    CALL set_acc_host_or_device(lzacc, lacc)


    !-------------------------------------------------------------------------
    params_oce      => v_params
    patch_2D        => patch_3d%p_patch_2d(1)
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho ! already calculated
    wind_mixing     => params_oce%velocity_windMixing(:,:,blockNo) ! already calculated
    levels = n_zlev

    !-------------------------------------------------------------------------
    z_grav_rho                   = grav/OceanReferenceDensity
    z_inv_OceanReferenceDensity  = 1.0_wp/OceanReferenceDensity
    !-------------------------------------------------------------------------

    !$ACC DATA CREATE(z_vert_density_grad_e) IF(lzacc)

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO je = start_index, end_index

      cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
      cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
      cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
      cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

      DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
        z_vert_density_grad_e(jk) =  0.5_wp * &
          & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
          &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
      ENDDO

      DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
        ! TODO: the following expect equally sized cells
        ! compute density gradient at edges
        dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
          &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo))
        density_differ_edge = z_vert_density_grad_e(jk) * dz
!         vn_diff = MAX( &
!           & ABS(ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
!           &      ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)), &
!           & ABS(vn_predict(je,jk) - &
!           &     vn_predict(je,jk-1)))
        vn_diff =  &
          & ABS(vn_predict(je,jk) - &
          &     vn_predict(je,jk-1))

        z_shear_edge = dbl_eps + vn_diff**2

        richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)

        new_velocity_friction = &
          & params_oce%a_veloc_v_back * dz +                              &
          & velocity_RichardsonCoeff / ((1.0_wp + z_c1_v * richardson_edge)**2)+  &
          & wind_mixing(je,jk)

        ! the average of the calculated velocity friction based on the old velocity and the predicted one
        params_oce%a_veloc_v(je,jk,blockNo) = &
!           & MAX(params_oce%a_veloc_v(je,jk,blockNo), new_velocity_friction )
           & VerticalViscosity_TimeWeight * params_oce%a_veloc_v(je,jk,blockNo) + &
           & (1.0_wp - VerticalViscosity_TimeWeight) * new_velocity_friction

      END DO ! jk = 2, levels
    END DO ! je = start_index, end_index
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    !$ACC END DATA

  END SUBROUTINE ICON_PP_Edge_vnPredict_scheme
  !-------------------------------------------------------------------------


END MODULE mo_ocean_pp_scheme
