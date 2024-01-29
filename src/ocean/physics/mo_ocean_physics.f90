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
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_physics
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: &
#ifdef __LVECTOR__
    & eos_type,                                               &
#endif
    & n_zlev, bottom_drag_coeff,                              &
    & HarmonicViscosity_reference, velocity_VerticalDiffusion_background,                 &
    & Temperature_VerticalDiffusion_background, Salinity_VerticalDiffusion_background, no_tracer,                       &
    & tracer_convection_MixingCoefficient,                                     &
    & BiharmonicViscosity_scaling, HarmonicViscosity_scaling, &
    & VelocityDiffusion_order,                                &
    & n_points_in_munk_layer,                                 &
    & BiharmonicViscosity_reference,                          &
    & tracer_RichardsonCoeff, velocity_RichardsonCoeff,                    &
    & use_wind_mixing,                                        &
    & vert_mix_type, vmix_pp, vmix_tke, vmix_idemix_tke,         & ! by_nils / by_ogut
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
    & OceanReferenceDensity,                                  &
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
    &  TracerDiffusion_LeithWeight,             &
    &  max_turbulenece_TracerDiffusion,                       &
    &  LeithViscosity_SmoothIterations, LeithViscosity_SpatialSmoothFactor
  USE mo_ocean_surface_types,            ONLY:  &
      t_ocean_surface                             ! contains p_oce_sfc

  USE mo_ocean_physics_types, ONLY: t_ho_params, v_params, WindMixingDecay, WindMixingLevel
   !, l_convection, l_pp_scheme
  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_impl_constants,      ONLY: success, max_char_length, min_dolic, sea
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_edge, grid_unstructured_cell
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_onEdges_Pointer_3d_wp, t_onCells_HalfLevels_Pointer_wp, t_operator_coeff
  USE mo_ocean_state,         ONLY: oce_config
  USE mo_physical_constants,  ONLY: grav, sitodbar,sal_ref
  USE mo_math_constants,      ONLY: dbl_eps, pi, rad2deg
  USE mo_dynamics_config,     ONLY: nold!, nnew
  USE mo_run_config,          ONLY: dtime
  USE mo_var_list,            ONLY: t_var_list_ptr
!  USE mo_cf_convention
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: datatype_pack16, DATATYPE_FLT32, DATATYPE_FLT64, filetype_nc2, &
    &                               GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: &
    & za_depth_below_sea, za_depth_below_sea_half, za_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max, sync_patch_array_mult
  USE mo_ocean_thermodyn,    ONLY: calculate_density_onColumn, calculate_density_onColumn_elem, calculate_density_mpiom_onColumn
  USE mo_ocean_math_operators,ONLY: div_oce_3d
  USE mo_timer,               ONLY: timers_level, timer_start, timer_stop, timer_upd_phys
  USE mo_statistics,          ONLY: global_minmaxmean
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_ocean_pp_scheme,     ONLY: update_PP_scheme, update_PP_scheme_zstar
  USE mo_ocean_tke,           ONLY: calc_tke, setup_tke
  USE mo_ocean_idemix,        ONLY: calc_idemix, setup_idemix
  USE mo_ocean_physics_types, ONLY: t_ho_params, v_params, &
   & WindMixingDecay, WindMixingLevel
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_velocity_diffusion, ONLY: veloc_diff_harmonic_div_grad
  USE mo_ocean_GM_Redi,       ONLY: init_GMRedi
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_physics'
  CHARACTER(LEN=12)           :: str_module    = 'ocePhysics  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug

  ! Public interface

  !PUBLIC :: init_ho_physics
  PUBLIC :: init_ho_params
  PUBLIC :: update_ho_params
  PUBLIC :: update_ho_params_zstar
!   PUBLIC :: calc_characteristic_physical_numbers
  PUBLIC :: scale_horizontal_diffusion, copy2Dto3D

  ! variables
  TYPE (t_var_list_ptr), PUBLIC :: ocean_params_list

  CHARACTER(LEN=*), PARAMETER :: module_name = 'mo_ocean_physics'

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Initialisation of ocean physics
  !!
!<Optimize:inUse:initOnly>
  SUBROUTINE init_ho_params(  patch_3d, physics_param, fu10 )
    TYPE(t_patch_3d ),POINTER, INTENT(in) :: patch_3d
    TYPE (t_ho_params)                          :: physics_param
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10

    ! Local variables
    INTEGER :: i, i_no_trac
    INTEGER :: je, blockNo, jk,jc, il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_edges, owned_edges
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: tracer_basisCoeff(nproma, patch_3D%p_patch_2d(1)%nblks_e)

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ho_params'

    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_edges => patch_2D%edges%ALL
    owned_edges => patch_2D%edges%owned
    !-------------------------------------------------------------------------
    !Init from namelist
    physics_param%a_veloc_v_back = velocity_VerticalDiffusion_background
    physics_param%a_veloc_v      = velocity_VerticalDiffusion_background

    IF(GMRedi_configuration==GMRedi_combined&
      &.OR.GMRedi_configuration==GM_only.OR.GMRedi_configuration==Redi_only)THEN

      CALL init_GMRedi(patch_3d)

      physics_param%k_tracer_isoneutral = k_tracer_isoneutral_parameter
      physics_param%k_tracer_dianeutral = k_tracer_dianeutral_parameter

      physics_param%k_tracer_GM_kappa = k_tracer_GM_kappa_parameter
    ENDIF

    CALL dbg_print('edge area:', patch_2D%edges%area_edge, str_module,0, &
      & in_subset=patch_2D%edges%owned)
    CALL dbg_print('edge lentgh:', patch_2D%edges%primal_edge_length, str_module,0, &
      & in_subset=patch_2D%edges%owned)
    CALL dbg_print('dual edge lentgh:', patch_2D%edges%dual_edge_length, str_module,0, &
      & in_subset=patch_2D%edges%owned)


    IF (VelocityDiffusion_order == 1 .OR. VelocityDiffusion_order == 21 ) THEN
      ! harmonic or harmonic+biharmonic
      CALL scale_horizontal_diffusion(patch_3D=patch_3D, DiffusionScaling=HarmonicViscosity_scaling, &
        & DiffusionReferenceValue=HarmonicViscosity_reference, &
        & DiffusionBackgroundValue=HarmonicViscosity_background,  &
        & out_DiffusionCoefficients=physics_param%HarmonicViscosity_BasisCoeff)
      CALL dbg_print('HarmonicVisc:'     ,physics_param%HarmonicViscosity_BasisCoeff,str_module,0, &
        & in_subset=patch_2D%edges%owned)
      IF (laplacian_form == 2) THEN
        ! divide by dual_edge_length to be used as a coeeficient directly to the grad operator
        CALL divideBy(physics_param%HarmonicViscosity_BasisCoeff, patch_2D%edges%dual_edge_length, &
          all_edges)
      ENDIF
      CALL copy2Dto3D(physics_param%HarmonicViscosity_BasisCoeff, physics_param%HarmonicViscosity_coeff, &
        all_edges)
    ENDIF

    IF (VelocityDiffusion_order == 2 .OR. VelocityDiffusion_order == 21 ) THEN
      ! biharmonic or harmonic+biharmonic
      CALL scale_horizontal_diffusion(patch_3D=patch_3D, DiffusionScaling=BiharmonicViscosity_scaling, &
        & DiffusionReferenceValue=BiharmonicViscosity_reference, &
        & DiffusionBackgroundValue=BiharmonicViscosity_background,  &
        & out_DiffusionCoefficients=physics_param%BiharmonicViscosity_BasisCoeff)
      CALL dbg_print('BiharmonicVisc:'     ,physics_param%BiharmonicViscosity_BasisCoeff,str_module,0, &
        & in_subset=patch_2D%edges%owned)
      CALL copy2Dto3D(physics_param%BiharmonicViscosity_BasisCoeff, physics_param%BiharmonicViscosity_coeff, &
        all_edges)
    ENDIF

    IF (LeithClosure_order == 1 .or.  LeithClosure_order == 21) THEN
      CALL scale_horizontal_diffusion(patch_3D=patch_3D, DiffusionScaling = LeithHarmonicViscosity_scaling, &
        & DiffusionReferenceValue  = LeithHarmonicViscosity_reference, &
        & DiffusionBackgroundValue = LeithHarmonicViscosity_background,  &
        & out_DiffusionCoefficients= physics_param%LeithHarmonicViscosity_BasisCoeff)
      CALL dbg_print('LeithHarmVisc:'     ,physics_param%LeithHarmonicViscosity_BasisCoeff,str_module,0, &
        & in_subset=patch_2D%edges%owned)
    ENDIF
    IF (LeithClosure_order == 2 .or.  LeithClosure_order == 21) THEN
      CALL scale_horizontal_diffusion(patch_3D=patch_3D, DiffusionScaling = LeithBiharmonicViscosity_scaling, &
        & DiffusionReferenceValue  = LeithBiharmonicViscosity_reference, &
        & DiffusionBackgroundValue = LeithBiharmonicViscosity_background,  &
        & out_DiffusionCoefficients= physics_param%LeithBiharmonicViscosity_BasisCoeff)
      CALL dbg_print('LeithBiharmVisc:'     ,physics_param%LeithBiharmonicViscosity_BasisCoeff,str_module,0, &
        & in_subset=patch_2D%edges%owned)
    ENDIF

    DO i=1,no_tracer

      IF(i==1)THEN!temperature
        physics_param%Tracer_HorizontalDiffusion_Background(i) = Temperature_HorizontalDiffusion_Background
        physics_param%Tracer_HorizontalDiffusion_Reference(i) = Temperature_HorizontalDiffusion_Reference
        physics_param%a_tracer_v_back(i) = Temperature_VerticalDiffusion_background
      ELSEIF(i==2)THEN!salinity
        physics_param%Tracer_HorizontalDiffusion_Background(i) = Salinity_HorizontalDiffusion_Background
        physics_param%Tracer_HorizontalDiffusion_Reference(i) = Salinity_HorizontalDiffusion_Reference
        physics_param%a_tracer_v_back(i) = Salinity_VerticalDiffusion_background
      ELSE

        CALL finish (method_name,'number of tracers exceeds number of background values')
      ENDIF

      CALL scale_horizontal_diffusion(patch_3D=patch_3D, &
        & DiffusionScaling=TracerHorizontalDiffusion_scaling, &
        & DiffusionReferenceValue=physics_param%Tracer_HorizontalDiffusion_Reference(i), &
        & DiffusionBackgroundValue=physics_param%Tracer_HorizontalDiffusion_Background(i), &
        & out_DiffusionCoefficients=physics_param%TracerDiffusion_BasisCoeff(:,:,i))
      CALL copy2Dto3D(physics_param%TracerDiffusion_BasisCoeff(:,:,i), physics_param%TracerDiffusion_coeff(:,:,:,i), &
        all_edges)
      CALL dbg_print('Tracer Diff Basis:', physics_param%TracerDiffusion_BasisCoeff(:,:,i),str_module,0, &
        & in_subset=patch_2D%edges%owned)
      CALL dbg_print('Tracer Diff:', physics_param%TracerDiffusion_coeff(:,:,:,i),str_module,0, &
        & in_subset=patch_2D%edges%owned)


      physics_param%a_tracer_v(:,:,:,i) = physics_param%a_tracer_v_back(i)

    END DO

    physics_param%bottom_drag_coeff = bottom_drag_coeff
    !$ACC UPDATE DEVICE(physics_param%A_veloc_v, physics_param%a_tracer_v, physics_param%bottom_drag_coeff)

    ! precalculate exponential wind mixing decay with depth
    DO jk=2,n_zlev
      WindMixingDecay(jk) = EXP(-patch_3d%p_patch_1d(1)%del_zlev_m(jk-1)/WindMixingDecayDepth)
      WindMixingLevel(jk) = lambda_wind * patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk-1)
    ENDDO

    ! setup tke scheme
    SELECT CASE(vert_mix_type)
    CASE(vmix_pp)
!      write(*,*) ''
      CALL message(method_name,'Setup vmix_pp scheme.')

    CASE(vmix_tke) ! by_nils
!      write(*,*) 'Setup tke scheme.'
      CALL message(method_name,'Setup tke scheme.')
      CALL setup_tke()
    CASE(vmix_idemix_tke) ! by_nils
!      write(*,*) 'Setup idemix_tke scheme.'
      CALL message(method_name,'Setup idemix_tke scheme.')
      CALL setup_idemix(patch_3d)
      CALL setup_tke()
    CASE default
!      write(*,*) "Unknown vert_mix_type!"
!      stop
      CALL finish(method_name, 'Unknown vert_mix_type!')

    END SELECT

  END SUBROUTINE init_ho_params
  !-------------------------------------------------------------------------

   !-------------------------------------------------------------------------
  REAL(wp) FUNCTION edge_area(patch_2D, je,jb)
    TYPE(t_patch), POINTER :: patch_2D
    INTEGER, INTENT(in) :: je, jb

    INTEGER :: c1_idx, c1_blk,  c2_idx, c2_blk

    c1_idx = patch_2D%edges%cell_idx(je,jb,1)
    c1_blk = patch_2D%edges%cell_blk(je,jb,1)
    c2_idx = patch_2D%edges%cell_idx(je,jb,2)
    c2_blk = patch_2D%edges%cell_blk(je,jb,2)

    edge_area = (patch_2D%cells%area(c1_idx, c1_blk) + patch_2D%cells%area(c2_idx, c2_blk)) * 0.5_wp

!     edge_area = patch_2D%edges%primal_edge_length(je,jb)*patch_2D%edges%dual_edge_length(je,jb))
  END FUNCTION edge_area
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE scale_horizontal_diffusion(patch_3D, &
    & DiffusionScaling, DiffusionReferenceValue, DiffusionBackgroundValue, out_DiffusionCoefficients)
    TYPE(t_patch_3d), POINTER :: patch_3D
    INTEGER, INTENT(in) :: DiffusionScaling
    REAL(wp), INTENT(in) :: DiffusionReferenceValue, DiffusionBackgroundValue
    REAL(wp) ::  out_DiffusionCoefficients(:,:)

    TYPE(t_patch), POINTER :: patch_2D
    INTEGER :: i, je, jb, jk
    INTEGER :: start_index, end_index
    REAL(wp) :: z_lower_bound_diff, C_MPIOM
    REAL(wp) :: z_diff_multfac, z_diff_efdt_ratio
    REAL(wp) :: points_in_munk_layer
    REAL(wp) :: minmaxmean_length(3), reference_scale
    REAL(wp) :: minDualEdgeLength, minEdgeLength, meanDualEdgeLength, meanEdgeLength
    REAL(wp) :: maxDualEdgeLength, maxEdgeLength
    REAL(wp) :: minCellArea, meanCellArea, maxCellArea
    TYPE(t_subset_range), POINTER :: all_edges, owned_edges
    REAL(wp):: length_scale, dual_length_scale, prime_length_scale

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':scale_horizontal_diffusion'

    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_edges => patch_2D%edges%ALL
    owned_edges => patch_2D%edges%owned
    !-------------------------------------------------------------------------
    points_in_munk_layer = REAL(n_points_in_munk_layer,wp)

    minmaxmean_length = global_minmaxmean(patch_2D%edges%dual_edge_length, owned_edges)
    minDualEdgeLength = minmaxmean_length(1)
    meanDualEdgeLength = minmaxmean_length(3)
    maxDualEdgeLength = minmaxmean_length(2)
    minmaxmean_length = global_minmaxmean(patch_2D%edges%primal_edge_length, owned_edges)
    minEdgeLength = minmaxmean_length(1)
    meanEdgeLength = minmaxmean_length(3)
    maxEdgeLength = minmaxmean_length(2)
    minmaxmean_length = global_minmaxmean(patch_2D%cells%area, patch_2D%cells%owned)
    minCellArea = minmaxmean_length(1)
    meanCellArea = minmaxmean_length(3)
    maxCellArea = minmaxmean_length(2)

    SELECT CASE(DiffusionScaling)

    CASE(0)
      out_DiffusionCoefficients(:,:) = 0.0_wp

    CASE(1)
      out_DiffusionCoefficients(:,:) = DiffusionReferenceValue + DiffusionBackgroundValue

!     CASE(2)!calculate coefficients based on Leith closure. Start value is the same as case(2).
!       out_DiffusionCoefficients(:,:) = 3.82E-12_wp * &
!         & (points_in_munk_layer * maxEdgeLength)**3


   CASE(2)
      ! linear scale
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + DiffusionReferenceValue * &
            & SQRT(edge_area(patch_2D,je,jb))

        END DO
      END DO

 !   CASE(3)! calculate coefficients for each location based on MUNK layer
 !     DO jb = all_edges%start_block, all_edges%end_block
 !       CALL get_index_range(all_edges, jb, start_index, end_index)
 !       out_DiffusionCoefficients(:,jb) = 0.0_wp
 !       DO je = start_index, end_index
 !         !calculate lower bound for diffusivity
 !         !The factor cos(lat) is omitted here, because of equatorial reference (cf. Griffies, eq. (18.29))
 !         out_DiffusionCoefficients(je,jb) = 3.82E-12_wp  &
 !           & *(points_in_munk_layer * patch_2D%edges%primal_edge_length(je,jb))**3
 !       END DO
 !     END DO

   CASE(3)
      ! **2
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + DiffusionReferenceValue * &
            & edge_area(patch_2D,je,jb)

        END DO
      END DO

    ! **3
    CASE(4)
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + DiffusionReferenceValue * &
            & SQRT(edge_area(patch_2D,je,jb))**3

        END DO
      END DO

    CASE(7)
      ! multiply DiffusionReferenceValue by dual_edge_length**3
      ! recommended values:
      !  Harmonic viscosity: 1.5E-11
      !  Biharmonic viscosicity: 3.15e-3
      !  Tracer diffusion: 5.0E-13
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + DiffusionReferenceValue * patch_2D%edges%dual_edge_length(je,jb)**3

        END DO
      END DO

    CASE(8)
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + DiffusionReferenceValue * patch_2D%edges%primal_edge_length(je,jb)**3

        END DO
      END DO

    CASE(9)
      ! multiply DiffusionReferenceValue by sqrt(dual_edge_length**3)
      ! recommended values:
      !  Biharmonic viscosicity:
      !  Tracer diffusion: 1.25E-5
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + DiffusionReferenceValue * SQRT(patch_2D%edges%dual_edge_length(je,jb)**3)

        END DO
      END DO

    CASE(10)
      ! **4
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + DiffusionReferenceValue * &
            & edge_area(patch_2D,je,jb)**2

        END DO
      END DO

!     CASE(11)
!       !The number that controls all that the "z_diff_efdt_ratio"
!       !is different. Higher z_diff_efdt_ratio decreases the final
!       !diffusion coefficient
!       z_diff_efdt_ratio = 10000.0_wp * DiffusionReferenceValue
!       z_diff_multfac = (1._wp/ (z_diff_efdt_ratio*64._wp))/3._wp
!       DO jb = all_edges%start_block, all_edges%end_block
!         CALL get_index_range(all_edges, jb, start_index, end_index)
!         DO je = start_index, end_index
!           out_DiffusionCoefficients(je,jb) = &
!           & patch_2D%edges%area_edge(je,jb)*patch_2D%edges%area_edge(je,jb)*z_diff_multfac
!         END DO
!       END DO
!       !          z_diff_multfac = 0.0045_wp*dtime/3600.0_wp
!       !         DO jb = all_edges%start_block, all_edges%end_block
!       !            CALL get_index_range(all_edges, jb, start_index, end_index)
!       !            DO je = start_index, end_index
!       !              physics_param%HarmonicViscosity_reference(je,:,jb) = z_diff_multfac*&
!       !              &maxval(patch_2D%edges%primal_edge_length)**4
!       !            END DO
!       !          END DO

   CASE(12)
      ! Simple scaling of the backgound diffusion by the dual edge length^4
      ! This is meant to be used with the non-uniform grids
      !This follows the MPI-OM convention
      C_MPIOM = DiffusionReferenceValue*dtime/3600.0_wp
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        DO je = start_index, end_index

          length_scale = &
          & sqrt(edge_area(patch_2D,je,jb))

          out_DiffusionCoefficients(je,jb)=C_MPIOM*length_scale**2

        END DO
      END DO

    CASE(13)
      ! Simple scaling of the constant diffusion by prime edge lenght + cell area (estimated)
      reference_scale = 4.0_wp / (3.0_wp * maxEdgeLength * maxCellArea)
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index

          length_scale = &
            & patch_2D%edges%primal_edge_length(je,jb)**2 * patch_2D%edges%dual_edge_length(je,jb) &
            & * reference_scale

            out_DiffusionCoefficients(je,jb) = &
              & DiffusionBackgroundValue +        &
              & DiffusionReferenceValue * length_scale

        END DO
      END DO

    CASE(14)
      ! Simple scaling of the constant diffusion using the dual edge
      ! this is the default scaling
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index
          dual_length_scale = patch_2D%edges%dual_edge_length(je,jb) / maxDualEdgeLength
          length_scale = dual_length_scale**3

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + &
            & DiffusionReferenceValue * length_scale

        END DO
      END DO

    CASE(15)
      ! Simple scaling of the constant diffusion using the prime edge
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        out_DiffusionCoefficients(:,jb) = 0.0_wp
        DO je = start_index, end_index
          prime_length_scale = patch_2D%edges%primal_edge_length(je,jb) / maxEdgeLength
          length_scale = prime_length_scale**3

          out_DiffusionCoefficients(je,jb) = &
            & DiffusionBackgroundValue + &
            & DiffusionReferenceValue * length_scale

        END DO
      END DO

    CASE DEFAULT
        CALL finish (method_name, 'uknown DiffusionScaling')

    END SELECT

    ! smooth if requested
    DO i=1, HorizontalViscosity_SmoothIterations
      CALL smooth_lapl_diff( patch_3d, &
        & out_DiffusionCoefficients, HorizontalViscosity_SpatialSmoothFactor )
    ENDDO

  END SUBROUTINE scale_horizontal_diffusion
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Calculation of a lower bound for horizontal velocity diffusion of laplacian type ! that is
  !! required to have N (default =1) points in Munk layer. The lower bound is calculated ! with
  !! respect to the equator.
  !! The code is based on  Griffies, Fundamentals of ocean climate modeling, sect 18, p. 413.  !
  !! The lower bound is given in units [m^2/s].
  !!
  !!
  !
  SUBROUTINE calc_lower_bound_veloc_diff(  patch_2D, lower_bound_diff )
    TYPE(t_patch), TARGET, INTENT(in)  :: patch_2D
    REAL(wp), INTENT(inout)              :: lower_bound_diff

    ! Local variables
    REAL(wp) :: points_in_munk_layer
    REAL(wp)            :: z_largest_edge_length
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-------------------------------------------------------------------------
    edges_in_domain => patch_2D%edges%in_domain

    ! Get the largest edge length globally
    z_largest_edge_length = global_max(MAXVAL(patch_2D%edges%primal_edge_length))

    !calculate lower bound for diffusivity: The factor cos(lat) is omitted here, because of
    !equatorial reference (cf. Griffies, eq.  (18.29))
    points_in_munk_layer = REAL(n_points_in_munk_layer,wp)
    lower_bound_diff = 3.82E-12_wp*(points_in_munk_layer*z_largest_edge_length)**3

  END SUBROUTINE calc_lower_bound_veloc_diff
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse:initOnly>
  SUBROUTINE smooth_lapl_diff( patch_3d, k_h, smoothFactor )
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)    :: k_h(:,:)
    REAL(wp), INTENT(in)       ::smoothFactor
    ! Local variables
    TYPE(t_patch), POINTER  :: patch_2D
    INTEGER :: je,jv,jk, jev, ile, ibe
    INTEGER :: il_v1,ib_v1, il_v2,ib_v2
    INTEGER :: start_index, end_index, blockNo
    INTEGER :: i_startidx_v, i_endidx_v
    REAL(wp), POINTER :: z_k_ave_v(:,:)
    INTEGER  :: sea_edges_onLevel
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER ::all_edges, verts_in_domain
    !-------------------------------------------------------------------------
    patch_2D => patch_3d%p_patch_2D(1)
    all_edges => patch_2D%edges%all
    verts_in_domain => patch_2D%verts%in_domain

    ALLOCATE(z_k_ave_v(nproma,patch_2D%nblks_v))
    z_k_ave_v(:,:) = 0.0_wp

    DO blockNo = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockNo, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        DO jev = 1, patch_2D%verts%num_edges(jv,blockNo)
          ile = patch_2D%verts%edge_idx(jv,blockNo,jev)
          ibe = patch_2D%verts%edge_blk(jv,blockNo,jev)
            z_k_ave_v(jv,blockNo)= z_k_ave_v(jv,blockNo) + k_h(ile,ibe)
        END DO ! jev = 1, patch_2D%verts%num_edges(jv,blockNo)
        z_k_ave_v(jv,blockNo) = z_k_ave_v(jv,blockNo) / REAL(patch_2D%verts%num_edges(jv,blockNo),wp)
      ENDDO
    END DO

    ! we do need to sync here
    CALL sync_patch_array(sync_v, patch_2D, z_k_ave_v)

    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index, end_index)
      DO je = start_index, end_index

        il_v1 = patch_2D%edges%vertex_idx(je,blockNo,1)
        ib_v1 = patch_2D%edges%vertex_blk(je,blockNo,1)
        il_v2 = patch_2D%edges%vertex_idx(je,blockNo,2)
        ib_v2 = patch_2D%edges%vertex_blk(je,blockNo,2)

        k_h(je,blockNo)= 0.5_wp * smoothFactor * (z_k_ave_v(il_v1,ib_v1) + z_k_ave_v(il_v2,ib_v2)) + &
          & (1.0_wp - smoothFactor) * k_h(je,blockNo)
      ENDDO
    END DO

    ! we do not need to sync edge coefficients
    ! CALL sync_patch_array(sync_e, patch_2D, k_h)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=1  ! output print levels - 0: print in any case
    CALL dbg_print('smoothed Laplac Diff.'     ,k_h                     ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
    DEALLOCATE(z_k_ave_v)

  END SUBROUTINE smooth_lapl_diff
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse:initOnly>
  SUBROUTINE smooth_lapl_diff_3D( patch_3d, k_h, smoothFactor )
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)    :: k_h(:,:,:)
    REAL(wp), INTENT(in)       ::smoothFactor
    ! Local variables
    TYPE(t_patch), POINTER  :: patch_2D
    INTEGER :: je,jv,jk, jev, ile, ibe, level
    INTEGER :: il_c1,ib_c1, il_c2,ib_c2
    INTEGER :: start_index, end_index, blockNo
    REAL(wp), POINTER :: z_k_ave_c(:,:,:)
    INTEGER  :: sea_edges_onLevel(n_zlev)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER ::all_edges, cells_in_domain
    !-------------------------------------------------------------------------
    patch_2D => patch_3d%p_patch_2D(1)
    all_edges => patch_2D%edges%all
    cells_in_domain => patch_2D%cells%in_domain

    ALLOCATE(z_k_ave_c(nproma,n_zlev, patch_2D%alloc_cell_blocks))

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, start_index, end_index, jv, jev, ile, ibe, sea_edges_onLevel, level)
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      z_k_ave_c(:,:,blockNo) = 0.0_wp

      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      DO jv = start_index, end_index
        sea_edges_onLevel(:) = 0
        DO jev = 1, patch_2D%cells%num_edges(jv,blockNo)
          ile = patch_2D%cells%edge_idx(jv,blockNo,jev)
          ibe = patch_2D%cells%edge_blk(jv,blockNo,jev)
          DO level=1,patch_3D%p_patch_1D(1)%dolic_e(ile,ibe)
            z_k_ave_c(jv,level,blockNo)= z_k_ave_c(jv,level,blockNo) + k_h(ile,level,ibe)
            sea_edges_onLevel(level) = sea_edges_onLevel(level) + 1
          END DO
        END DO
        DO level=1,n_zlev
          IF (sea_edges_onLevel(level) > 1) &
            z_k_ave_c(jv,level,blockNo) = z_k_ave_c(jv,level,blockNo) / REAL(sea_edges_onLevel(level),wp)
        ENDDO
      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! we do need to sync here
    CALL sync_patch_array(sync_c, patch_2D, z_k_ave_c)

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, start_index, end_index, je, il_c1, ib_c1, il_c2, ib_c2, level)
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index, end_index)
      DO je = start_index, end_index

        il_c1 = patch_2D%edges%cell_idx(je,blockNo,1)
        ib_c1 = patch_2D%edges%cell_blk(je,blockNo,1)
        il_c2 = patch_2D%edges%cell_idx(je,blockNo,2)
        ib_c2 = patch_2D%edges%cell_blk(je,blockNo,2)

        DO level=1,patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
          k_h(je,level,blockNo)= 0.5_wp * smoothFactor * (z_k_ave_c(il_c1,level,ib_c1) + z_k_ave_c(il_c2,level,ib_c2)) + &
            & (1.0_wp - smoothFactor) * k_h(je,level,blockNo)
        ENDDO
      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! we do not need to sync edge coefficients
    ! CALL sync_patch_array(sync_e, patch_2D, k_h)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=1  ! output print levels - 0: print in any case
    CALL dbg_print('smoothed Laplac Diff.'     ,k_h                     ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
    DEALLOCATE(z_k_ave_c)

  END SUBROUTINE smooth_lapl_diff_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse:initOnly>
  SUBROUTINE smooth_lapl_diff_verts_3D( patch_3d, k_h, smoothFactor )
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)    :: k_h(:,:,:)
    REAL(wp), INTENT(in)       ::smoothFactor
    ! Local variables
    TYPE(t_patch), POINTER  :: patch_2D
    INTEGER :: je,jv,jk, jev, ile, ibe, level
    INTEGER :: il_v1,ib_v1, il_v2,ib_v2
    INTEGER :: start_index, end_index, blockNo
    INTEGER :: i_startidx_v, i_endidx_v
    REAL(wp), POINTER :: z_k_ave_v(:,:,:)
    INTEGER  :: sea_edges_onLevel(n_zlev)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER ::all_edges, verts_in_domain
    !-------------------------------------------------------------------------
    patch_2D => patch_3d%p_patch_2D(1)
    all_edges => patch_2D%edges%all
    verts_in_domain => patch_2D%verts%in_domain

    ALLOCATE(z_k_ave_v(nproma,n_zlev, patch_2D%nblks_v))
    z_k_ave_v(:,:,:) = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, i_startidx_v, i_endidx_v, jv, jev, ile, ibe, sea_edges_onLevel, level)
    DO blockNo = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockNo, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        sea_edges_onLevel(:) = 0
        DO jev = 1, patch_2D%verts%num_edges(jv,blockNo)
          ile = patch_2D%verts%edge_idx(jv,blockNo,jev)
          ibe = patch_2D%verts%edge_blk(jv,blockNo,jev)
          DO level=1,patch_3D%p_patch_1D(1)%dolic_e(ile,ibe)
            z_k_ave_v(jv,level,blockNo)= z_k_ave_v(jv,level,blockNo) + k_h(ile,level,ibe)
            sea_edges_onLevel(level) = sea_edges_onLevel(level) + 1
          END DO
        END DO
        DO level=1,n_zlev
          IF (sea_edges_onLevel(level) > 1) &
            z_k_ave_v(jv,level,blockNo) = z_k_ave_v(jv,level,blockNo) / REAL(sea_edges_onLevel(level),wp)
        ENDDO
      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! we do need to sync here
    CALL sync_patch_array(sync_v, patch_2D, z_k_ave_v)

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, start_index, end_index, je, il_v1, ib_v1, il_v2, ib_v2, level)
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index, end_index)
      DO je = start_index, end_index

        il_v1 = patch_2D%edges%vertex_idx(je,blockNo,1)
        ib_v1 = patch_2D%edges%vertex_blk(je,blockNo,1)
        il_v2 = patch_2D%edges%vertex_idx(je,blockNo,2)
        ib_v2 = patch_2D%edges%vertex_blk(je,blockNo,2)

        DO level=1,patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
          k_h(je,level,blockNo)= 0.5_wp * smoothFactor * (z_k_ave_v(il_v1,level,ib_v1) + z_k_ave_v(il_v2,level,ib_v2)) + &
            & (1.0_wp - smoothFactor) * k_h(je,level,blockNo)
        ENDDO
      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! we do not need to sync edge coefficients
    ! CALL sync_patch_array(sync_e, patch_2D, k_h)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=1  ! output print levels - 0: print in any case
    CALL dbg_print('smoothed Laplac Diff.'     ,k_h                     ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
    DEALLOCATE(z_k_ave_v)

  END SUBROUTINE smooth_lapl_diff_verts_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Update of ocean physics parameters
  !!
  !! Update of ocean physics: This routine is used used only if time-dependent
  !! changes of physical parametrizations.
  !! Currently vertical mixing coefficients for tracers and vertical diffusivity are updated.
  !! Dependent on the local Richardson number the diffusivity are calculated
  !! (Large & Gent JPO 29, (1999), 449-464).
  !! The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean
  !! Modelling 5, 2003).
  !! The notational convention is also taken from this paper( cf. eqs (14) and (19)).
  !! What is missing is the fractional ice cover (see eqs. (15-16)).
  !! Eq. (18) is the Redi part that is not implemented, yet
  !!
!<Optimize:inUse:done>
  SUBROUTINE update_ho_params(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs, atmos_fluxes, p_oce_sfc, lacc)
    !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum(:,:) ! t_sea_ice%concsum
    TYPE(t_ho_params), INTENT(inout)     :: params_oce
    TYPE (t_ocean_surface), INTENT(IN)   :: p_oce_sfc
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    TYPE(t_atmos_fluxes)                 :: atmos_fluxes
    LOGICAL, INTENT(IN), OPTIONAL        :: lacc

    INTEGER :: tracer_index
    LOGICAL :: lzacc
    !INTEGER :: vert_mix_type=2 ! by_nils ! FIXME: make this a namelist parameter
    !-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    start_timer(timer_upd_phys,1)

!    CALL calc_characteristic_physical_numbers(patch_3d, ocean_state)

!   Calculate the vertical density gradient on the interfaces (zgrad_rho)
!   and the Richardson Number ; shall be used in PP and possibly TKE
    CALL calc_vertical_stability(patch_3d, ocean_state, lacc=lzacc)

    SELECT CASE(vert_mix_type)
    CASE(vmix_pp)
      CALL update_PP_scheme(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs)
    CASE(vmix_tke)
      !write(*,*) 'Do calc_tke...'
#ifdef _OPENACC
      CALL calc_tke(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum, lacc=lzacc)
#else
      CALL calc_tke(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum)
#endif
    CASE(vmix_idemix_tke)
      !write(*,*) 'Do calc_idemix...'
      CALL calc_idemix(patch_3d, ocean_state, params_oce, op_coeffs, atmos_fluxes)
      !CALL calc_tke(patch_3d, ocean_state, params_oce, atmos_fluxes)
      CALL calc_tke(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum)
    CASE default
      write(*,*) "Unknown vert_mix_type!"
    END SELECT

    IF (LeithClosure_order == 1 .or.  LeithClosure_order == 21) THEN
      IF (LeithClosure_form == 1) THEN
        CALL calculate_LeithClosure_harmonic_vort(patch_3d, ocean_state, params_oce, op_coeffs)
      ELSEIF (LeithClosure_form == 2) THEN
        CALL calculate_LeithClosure_harmonic_vort_w(patch_3d, ocean_state, params_oce, op_coeffs)
      ELSEIF (LeithClosure_form == 3) THEN
        CALL calculate_LeithClosure_harmonic_vort_div(patch_3d, ocean_state, params_oce, op_coeffs)
      ELSEIF (LeithClosure_form == 4) THEN
        CALL calculate_LeithClosure_harmonicDivGrad_VortDiv(patch_3d, ocean_state, params_oce, op_coeffs)
      ELSEIF (LeithClosure_form == 5) THEN
        CALL calculate_LeithClosure_harmonic_w(patch_3d, ocean_state, params_oce, op_coeffs)
      ELSEIF (LeithClosure_form == 6) THEN
        CALL calculate_LeithClosure_harmonic_w2(patch_3d, ocean_state, params_oce, op_coeffs)
      ENDIF
    ENDIF
    IF (LeithClosure_order == 2 .or.  LeithClosure_order == 21) THEN
      CALL calculate_LeithClosure_biharmonic(patch_3d, ocean_state, params_oce, op_coeffs)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    !idt_src=4  ! output print levels (1-5, fix)
    !CALL dbg_print('UpdPar: p_vn%x(1)'         ,ocean_state%p_diag%p_vn%x(1)    ,str_module,idt_src, &
    !  & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    !CALL dbg_print('UpdPar: p_vn%x(2)'         ,ocean_state%p_diag%p_vn%x(2)    ,str_module,idt_src, &
    !  & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     idt_src=2  ! output print levels (1-5, fix)
!     DO tracer_index = 1, no_tracer
!       CALL dbg_print('UpdPar FinalTracerMixing'  ,params_oce%a_tracer_v(:,:,:,tracer_index), str_module,idt_src, &
!         & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     ENDDO
!     CALL dbg_print('UpdPar FinalVelocMixing'   ,params_oce%a_veloc_v     ,str_module,idt_src, &
!       & in_subset=patch_3d%p_patch_2d(1)%edges%owned)
    !---------------------------------------------------------------------
    stop_timer(timer_upd_phys,1)

  END SUBROUTINE update_ho_params
  !-------------------------------------------------------------------------


  SUBROUTINE update_ho_params_zstar(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs, atmos_fluxes, p_oce_sfc, &
      & eta_c, stretch_c, stretch_e)
    !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum(:,:) ! t_sea_ice%concsum
    TYPE(t_ho_params), INTENT(inout)     :: params_oce
    TYPE (t_ocean_surface), INTENT(IN)   :: p_oce_sfc
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    TYPE(t_atmos_fluxes)                 :: atmos_fluxes
    REAL(wp), INTENT(IN) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor

    INTEGER :: tracer_index
    !INTEGER :: vert_mix_type=2 ! by_nils ! FIXME: make this a namelist parameter
    !-------------------------------------------------------------------------
    start_timer(timer_upd_phys,1)

!    CALL calc_characteristic_physical_numbers(patch_3d, ocean_state)

!   Calculate the vertical density gradient on the interfaces (zgrad_rho)
!   and the Richardson Number ; shall be used in PP and possibly TKE
     CALL calc_vertical_stability_zstar(patch_3d, ocean_state, eta_c, stretch_c)


    SELECT CASE(vert_mix_type)
    CASE(vmix_pp)
      CALL update_PP_scheme_zstar(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs, &
                                  & eta_c, stretch_c, stretch_e)
    CASE(vmix_tke)
      !write(*,*) 'Do calc_tke...'
      ! tke does not need a special routine for zstar
      CALL calc_tke(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum)
    CASE(vmix_idemix_tke)
      !write(*,*) 'Do calc_idemix...'
      CALL calc_idemix(patch_3d, ocean_state, params_oce, op_coeffs, atmos_fluxes)
      CALL calc_tke(patch_3d, ocean_state, params_oce, atmos_fluxes, fu10, concsum)
    CASE default
      write(*,*) "Unknown vert_mix_type!"
    END SELECT



    stop_timer(timer_upd_phys,1)

  END SUBROUTINE update_ho_params_zstar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the viscosity coefficients following the Leith closure.
  !! Implemented is the pure leith closure and a modified version of it. Both options
  !! are available for harmonic as well as for biharmonic diffusion.
  !!
  !! Liteature: Fox-Kemper, Menemenlis, Can Large Eddy Simulations improve Mesoscale Rich Ocean Models ?
  !!
  SUBROUTINE calculate_LeithClosure_harmonic_vort(patch_3d, ocean_state, param, operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: operators_coeff

    !Local variables
    INTEGER :: jk, blockNo, je, jc,jb, i
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
!     INTEGER :: LEITH_EXPONENT
    TYPE(t_subset_range), POINTER ::edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: div_e, grad_vort_abs
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain

    start_level = 1
   !-------------------------------------------------------------------------------
   !leith closure for Laplacian(harmonic) viscosity

!      LEITH_EXPONENT=3
    !1) calculation leith closure or modified LeithClosure_type
!     IF(LeithClosure_type==1)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level,end_level, &
!ICON_OMP vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk, grad_vort_abs, i) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index
        end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)
        DO level=start_level,end_level

          grad_vort_abs = &
            &  ABS((ocean_state%p_diag%vort(vertex1_idx,level,vertex1_blk) - &
            &       ocean_state%p_diag%vort(vertex2_idx,level,vertex2_blk))  &
            &  * patch_2D%edges%inv_primal_edge_length(je,blockNo))

          param%HarmonicViscosity_coeff(je,level,blockNo) = &
            & param%HarmonicViscosity_BasisCoeff(je,blockNo) + &
            & grad_vort_abs * param%LeithHarmonicViscosity_BasisCoeff(je,blockNo)

          DO i=1,no_tracer
            param%TracerDiffusion_coeff(je,level,blockNo,i) = &
              & param%TracerDiffusion_BasisCoeff(je,blockNo,i) +  &
              & MIN(grad_vort_abs * TracerDiffusion_LeithWeight * &
              &     param%LeithHarmonicViscosity_BasisCoeff(je,blockNo), &
              &     param%TracerDiffusion_BasisCoeff(je,blockNo,i) * max_turbulenece_TracerDiffusion)
          END DO

        END DO
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO
   ! this sync most probably is not needed
!    CALL sync_patch_array(sync_e, patch_2D, param%HarmonicViscosity_coeff)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('LeithClosure: viscosity',param%HarmonicViscosity_coeff,&
      & str_module,idt_src, in_subset=edges_in_domain)
    IF (TracerDiffusion_LeithWeight > 0.0_wp) THEN
      DO i=1,no_tracer
        CALL dbg_print('LeithClosure: tracer diff',param%TracerDiffusion_coeff(:,:,:,i),&
          &str_module,idt_src, in_subset=edges_in_domain)
      END DO
    ENDIF
!   !---------------------------------------------------------------------
    ! smooth if requested
    DO i=1, LeithViscosity_SmoothIterations
      CALL smooth_lapl_diff_3D( patch_3d, &
        & param%HarmonicViscosity_coeff, LeithViscosity_SpatialSmoothFactor )
    ENDDO

  END SUBROUTINE calculate_LeithClosure_harmonic_vort
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the viscosity coefficients following the Leith closure.
  !! Implemented is the pure leith closure and a modified version of it. Both options
  !! are available for harmonic as well as for biharmonic diffusion.
  !!
  !! Liteature: Fox-Kemper, Menemenlis, Can Large Eddy Simulations improve Mesoscale Rich Ocean Models ?
  !!
  SUBROUTINE calculate_LeithClosure_harmonic_vort_w(patch_3d, ocean_state, param, operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: operators_coeff

    !Local variables
    INTEGER :: jk, blockNo, je, jc,jb, i
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
!     INTEGER :: LEITH_EXPONENT
    TYPE(t_subset_range), POINTER ::edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: grad_w_e, grad_vort, LeithCoeff
!     REAL(wp):: div_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain

    start_level = 1
   !-------------------------------------------------------------------------------
   !leith closure for Laplacian(harmonic) viscosity

!      LEITH_EXPONENT=3
    !1) calculation leith closure or modified LeithClosure_type
!     ELSEIF(LeithClosure_type==2)THEN
!     CALL dbg_print('LeithClosure: ptp_vn',ocean_state%p_diag%ptp_vn,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!
!     CALL div_oce_3d( ocean_state%p_diag%ptp_vn, patch_3D, operators_coeff%div_coeff, div_c, &
!       & subset_range=patch_2d%cells%all)

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level,end_level, &
!ICON_OMP vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk, grad_vort, grad_w_e, &
!ICON_OMP cell1_idx, cell1_blk, cell2_idx, cell2_blk, LeithCoeff, i) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index
        end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)
        cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
        cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
        cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
        cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

        DO level=start_level,end_level

          grad_vort = &
            &  (ocean_state%p_diag%vort(vertex1_idx,level,vertex1_blk) - &
            &   ocean_state%p_diag%vort(vertex2_idx,level,vertex2_blk))  &
            &    * patch_2D%edges%inv_primal_edge_length(je,blockNo)

          grad_w_e = &
            & (  ocean_state%p_diag%w(cell2_idx,level,  cell2_blk) &
            &  + ocean_state%p_diag%w(cell2_idx,level+1,cell2_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level,  cell1_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level+1,cell1_blk)) &
            &    * 0.5_wp * patch_2D%edges%inv_dual_edge_length(je,blockNo)

          LeithCoeff = param%LeithHarmonicViscosity_BasisCoeff(je,blockNo) *  &
            &   sqrt(grad_vort**2  + grad_w_e**2)

          param%HarmonicViscosity_coeff(je,level,blockNo)  = &
            & param%HarmonicViscosity_BasisCoeff(je,blockNo) + &
            & LeithCoeff

!           param%TracerDiffusion_coeff(je,level,blockNo,1) = grad_vort_abs
!           param%TracerDiffusion_coeff(je,level,blockNo,2) = div_e

          DO i=1,no_tracer
            param%TracerDiffusion_coeff(je,level,blockNo,i) = &
              & param%TracerDiffusion_BasisCoeff(je,blockNo,i) + &
              & MIN(LeithCoeff * TracerDiffusion_LeithWeight,    &
              &     max_turbulenece_TracerDiffusion)
          END DO

        END DO
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO
   ! this sync most probably is not needed
!    CALL sync_patch_array(sync_e, patch_2D, param%HarmonicViscosity_coeff)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('Leith 2: Harm.viscos.',param%HarmonicViscosity_coeff,&
      & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: coeff.',param%LeithHarmonicViscosity_BasisCoeff,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: vort', ocean_state%p_diag%vort,&
!       & str_module,idt_src, in_subset=patch_2D%verts%in_domain)
!     CALL dbg_print('Leith 2: w',ocean_state%p_diag%w,&
!       & str_module,idt_src, in_subset=patch_2D%cells%in_domain)


!     CALL dbg_print('LeithClosure: grad_vort_abs',param%TracerDiffusion_coeff(:,:,:,1),&
!       &str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('LeithClosure: div_e',param%TracerDiffusion_coeff(:,:,:,2),&
!       &str_module,idt_src, in_subset=edges_in_domain)
    IF (TracerDiffusion_LeithWeight > 0.0_wp) THEN
!       DO i=1,no_tracer
        CALL dbg_print('LeithClosure: tracer diff',param%TracerDiffusion_coeff(:,:,:,1),&
          &str_module,idt_src, in_subset=edges_in_domain)
!       END DO
    ENDIF
!   !---------------------------------------------------------------------
    ! smooth if requested
    DO i=1, LeithViscosity_SmoothIterations
      CALL smooth_lapl_diff_3D( patch_3d, &
        & param%HarmonicViscosity_coeff, LeithViscosity_SpatialSmoothFactor )
    ENDDO

  END SUBROUTINE calculate_LeithClosure_harmonic_vort_w
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_LeithClosure_harmonic_w(patch_3d, ocean_state, param, operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: operators_coeff

    !Local variables
    INTEGER :: jk, blockNo, je, jc,jb, i
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
!     INTEGER :: LEITH_EXPONENT
    TYPE(t_subset_range), POINTER ::edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: grad_w_e,  LeithCoeff
!     REAL(wp):: div_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain

    start_level = 1
   !-------------------------------------------------------------------------------
   !leith closure for Laplacian(harmonic) viscosity

!      LEITH_EXPONENT=3
    !1) calculation leith closure or modified LeithClosure_type
!     ELSEIF(LeithClosure_type==2)THEN
!     CALL dbg_print('LeithClosure: ptp_vn',ocean_state%p_diag%ptp_vn,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!
!     CALL div_oce_3d( ocean_state%p_diag%ptp_vn, patch_3D, operators_coeff%div_coeff, div_c, &
!       & subset_range=patch_2d%cells%all)

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level,end_level, &
!ICON_OMP grad_w_e, &
!ICON_OMP cell1_idx, cell1_blk, cell2_idx, cell2_blk, LeithCoeff, i) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index
        end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

        cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
        cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
        cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
        cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

        DO level=start_level,end_level

          grad_w_e = &
            & (  ocean_state%p_diag%w(cell2_idx,level,  cell2_blk) &
            &  + ocean_state%p_diag%w(cell2_idx,level+1,cell2_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level,  cell1_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level+1,cell1_blk)) &
            &    * 0.5_wp * patch_2D%edges%inv_dual_edge_length(je,blockNo)

          LeithCoeff = param%LeithHarmonicViscosity_BasisCoeff(je,blockNo) *  &
            &   ABS(grad_w_e)

          param%HarmonicViscosity_coeff(je,level,blockNo)  = &
            & param%HarmonicViscosity_BasisCoeff(je,blockNo) + &
            & LeithCoeff

!           param%TracerDiffusion_coeff(je,level,blockNo,1) = grad_vort_abs
!           param%TracerDiffusion_coeff(je,level,blockNo,2) = div_e

        END DO
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO
   ! this sync most probably is not needed
!    CALL sync_patch_array(sync_e, patch_2D, param%HarmonicViscosity_coeff)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('Leith 5: Harm.viscos.',param%HarmonicViscosity_coeff,&
      & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: coeff.',param%LeithHarmonicViscosity_BasisCoeff,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: vort', ocean_state%p_diag%vort,&
!       & str_module,idt_src, in_subset=patch_2D%verts%in_domain)
!     CALL dbg_print('Leith 2: w',ocean_state%p_diag%w,&
!       & str_module,idt_src, in_subset=patch_2D%cells%in_domain)


!     CALL dbg_print('LeithClosure: grad_vort_abs',param%TracerDiffusion_coeff(:,:,:,1),&
!       &str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('LeithClosure: div_e',param%TracerDiffusion_coeff(:,:,:,2),&
!       &str_module,idt_src, in_subset=edges_in_domain)
    !---------------------------------------------------------------------
    ! smooth if requested
    DO i=1, LeithViscosity_SmoothIterations
      CALL smooth_lapl_diff_3D( patch_3d, &
        & param%HarmonicViscosity_coeff, LeithViscosity_SpatialSmoothFactor )
    ENDDO

    !---------------------------------------------------------------------
    ! tracer diffusion
!     IF (TracerDiffusion_LeithWeight > 0.0_wp) THEN
!       DO i=1,no_tracer
!       param%TracerDiffusion_coeff(je,level,blockNo,i) = &
!         & param%TracerDiffusion_BasisCoeff(je,blockNo,i) + &
!         & MIN(param%HarmonicViscosity_coeff(je,1,blockNo) * TracerDiffusion_LeithWeight,    &
!         &     max_turbulenece_TracerDiffusion)
!       END DO
! !       DO i=1,no_tracer
!       CALL dbg_print('LeithClosure: tracer diff',param%TracerDiffusion_coeff(:,:,:,1),&
!       &str_module,idt_src, in_subset=edges_in_domain)
! !       END DO
!     ENDIF

  END SUBROUTINE calculate_LeithClosure_harmonic_w
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  SUBROUTINE calculate_LeithClosure_harmonic_w2(patch_3d, ocean_state, param, operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: operators_coeff

    !Local variables
    INTEGER :: jk, blockNo, je, jc,jb, i
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
!     INTEGER :: LEITH_EXPONENT
    TYPE(t_subset_range), POINTER ::edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: grad_w_e,  LeithCoeff
!     REAL(wp):: div_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain

    start_level = 1
   !-------------------------------------------------------------------------------
   !leith closure for Laplacian(harmonic) viscosity

!      LEITH_EXPONENT=3
    !1) calculation leith closure or modified LeithClosure_type
!     ELSEIF(LeithClosure_type==2)THEN
!     CALL dbg_print('LeithClosure: ptp_vn',ocean_state%p_diag%ptp_vn,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!
!     CALL div_oce_3d( ocean_state%p_diag%ptp_vn, patch_3D, operators_coeff%div_coeff, div_c, &
!       & subset_range=patch_2d%cells%all)

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level,end_level, &
!ICON_OMP grad_w_e, &
!ICON_OMP cell1_idx, cell1_blk, cell2_idx, cell2_blk, LeithCoeff, i) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index
        end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

        cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
        cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
        cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
        cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

        DO level=start_level,end_level

          grad_w_e = &
            & (  ocean_state%p_diag%w(cell2_idx,level,  cell2_blk) &
            &  + ocean_state%p_diag%w(cell2_idx,level+1,cell2_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level,  cell1_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level+1,cell1_blk)) &
            &    * 0.5_wp * patch_2D%edges%inv_dual_edge_length(je,blockNo)

          LeithCoeff = param%LeithHarmonicViscosity_BasisCoeff(je,blockNo) *  &
            &   (grad_w_e**2)

          param%HarmonicViscosity_coeff(je,level,blockNo)  = &
            & param%HarmonicViscosity_BasisCoeff(je,blockNo) + &
            & LeithCoeff

!           param%TracerDiffusion_coeff(je,level,blockNo,1) = grad_vort_abs
!           param%TracerDiffusion_coeff(je,level,blockNo,2) = div_e

        END DO
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO
   ! this sync most probably is not needed
!    CALL sync_patch_array(sync_e, patch_2D, param%HarmonicViscosity_coeff)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('Leith 6: Harm.viscos.',param%HarmonicViscosity_coeff,&
      & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: coeff.',param%LeithHarmonicViscosity_BasisCoeff,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: vort', ocean_state%p_diag%vort,&
!       & str_module,idt_src, in_subset=patch_2D%verts%in_domain)
!     CALL dbg_print('Leith 2: w',ocean_state%p_diag%w,&
!       & str_module,idt_src, in_subset=patch_2D%cells%in_domain)


!     CALL dbg_print('LeithClosure: grad_vort_abs',param%TracerDiffusion_coeff(:,:,:,1),&
!       &str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('LeithClosure: div_e',param%TracerDiffusion_coeff(:,:,:,2),&
!       &str_module,idt_src, in_subset=edges_in_domain)
    !---------------------------------------------------------------------
    ! smooth if requested
    DO i=1, LeithViscosity_SmoothIterations
      CALL smooth_lapl_diff_3D( patch_3d, &
        & param%HarmonicViscosity_coeff, LeithViscosity_SpatialSmoothFactor )
    ENDDO

    !---------------------------------------------------------------------
    ! tracer diffusion
!     IF (TracerDiffusion_LeithWeight > 0.0_wp) THEN
!       DO i=1,no_tracer
!   param%TracerDiffusion_coeff(je,level,blockNo,i) = &
!     & param%TracerDiffusion_BasisCoeff(je,blockNo,i) + &
!     & MIN(param%HarmonicViscosity_coeff(je,1,blockNo) * TracerDiffusion_LeithWeight,    &
!     &     max_turbulenece_TracerDiffusion)
!       END DO
! !       DO i=1,no_tracer
!       CALL dbg_print('LeithClosure: tracer diff',param%TracerDiffusion_coeff(:,:,:,1),&
!   &str_module,idt_src, in_subset=edges_in_domain)
! !       END DO
!     ENDIF

  END SUBROUTINE calculate_LeithClosure_harmonic_w2
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>  !! !  SUBROUTINE calculates the viscosity coefficients following the Leith closure.
  !! Implemented is the pure leith closure and a modified version of it. Both options
  !! are available for harmonic as well as for biharmonic diffusion.
  !! this is the original equation
  !!
  !! Liteature: Fox-Kemper, Menemenlis, Can Large Eddy Simulations improve Mesoscale Rich Ocean Models ?
  !!
  SUBROUTINE calculate_LeithClosure_harmonic_vort_div(patch_3d, ocean_state, param, operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: operators_coeff

    !Local variables
    INTEGER :: jk, blockNo, je, jc,jb, i
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
!     INTEGER :: LEITH_EXPONENT
    TYPE(t_subset_range), POINTER ::edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: grad_w_e, grad_vort, LeithCoeff
!     REAL(wp):: div_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain

    start_level = 1
   !-------------------------------------------------------------------------------
   !leith closure for Laplacian(harmonic) viscosity

!      LEITH_EXPONENT=3
    !1) calculation leith closure or modified LeithClosure_type
!     ELSEIF(LeithClosure_type==2)THEN
!     CALL dbg_print('LeithClosure: ptp_vn',ocean_state%p_diag%ptp_vn,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!
!     CALL div_oce_3d( ocean_state%p_diag%ptp_vn, patch_3D, operators_coeff%div_coeff, div_c, &
!       & subset_range=patch_2d%cells%all)

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level,end_level, &
!ICON_OMP vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk, grad_vort, grad_w_e, &
!ICON_OMP cell1_idx, cell1_blk, cell2_idx, cell2_blk, LeithCoeff, i) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index
        end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)
        cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
        cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
        cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
        cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

        DO level=start_level,end_level

          grad_vort = &
            &  (ocean_state%p_diag%vort(vertex1_idx,level,vertex1_blk) - &
            &   ocean_state%p_diag%vort(vertex2_idx,level,vertex2_blk))  &
            &    * patch_2D%edges%inv_primal_edge_length(je,blockNo)

          grad_w_e = &
            & (  ocean_state%p_diag%w(cell2_idx,level+1,  cell2_blk) &
            &  - ocean_state%p_diag%w(cell2_idx,level,cell2_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level+1,  cell1_blk) &
            &  + ocean_state%p_diag%w(cell1_idx,level+1,cell1_blk)) &
            &    * 0.5_wp * patch_2D%edges%inv_dual_edge_length(je,blockNo)

          LeithCoeff = param%LeithHarmonicViscosity_BasisCoeff(je,blockNo) *  &
            &   sqrt(grad_vort**2  + grad_w_e**2)

          param%HarmonicViscosity_coeff(je,level,blockNo)  = &
            & param%HarmonicViscosity_BasisCoeff(je,blockNo) + &
            & LeithCoeff

!           param%TracerDiffusion_coeff(je,level,blockNo,1) = grad_vort_abs
!           param%TracerDiffusion_coeff(je,level,blockNo,2) = div_e

          DO i=1,no_tracer
            param%TracerDiffusion_coeff(je,level,blockNo,i) = &
              & param%TracerDiffusion_BasisCoeff(je,blockNo,i) + &
              & MIN(LeithCoeff * TracerDiffusion_LeithWeight,    &
              &     max_turbulenece_TracerDiffusion)
          END DO

        END DO
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO
   ! this sync most probably is not needed
!    CALL sync_patch_array(sync_e, patch_2D, param%HarmonicViscosity_coeff)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('Leith 3: Harm.viscos.',param%HarmonicViscosity_coeff,&
      & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: coeff.',param%LeithHarmonicViscosity_BasisCoeff,&
!       & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('Leith 2: vort', ocean_state%p_diag%vort,&
!       & str_module,idt_src, in_subset=patch_2D%verts%in_domain)
!     CALL dbg_print('Leith 2: w',ocean_state%p_diag%w,&
!       & str_module,idt_src, in_subset=patch_2D%cells%in_domain)


!     CALL dbg_print('LeithClosure: grad_vort_abs',param%TracerDiffusion_coeff(:,:,:,1),&
!       &str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('LeithClosure: div_e',param%TracerDiffusion_coeff(:,:,:,2),&
!       &str_module,idt_src, in_subset=edges_in_domain)
    IF (TracerDiffusion_LeithWeight > 0.0_wp) THEN
!       DO i=1,no_tracer
        CALL dbg_print('LeithClosure: tracer diff',param%TracerDiffusion_coeff(:,:,:,1),&
          &str_module,idt_src, in_subset=edges_in_domain)
!       END DO
    ENDIF
!   !---------------------------------------------------------------------
    ! smooth if requested
    DO i=1, LeithViscosity_SmoothIterations
      CALL smooth_lapl_diff_3D( patch_3d, &
        & param%HarmonicViscosity_coeff, LeithViscosity_SpatialSmoothFactor )
    ENDDO

  END SUBROUTINE calculate_LeithClosure_harmonic_vort_div
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the viscosity coefficients following the Leith closure.
  !! Implemented is the pure leith closure and a modified version of it. Both options
  !! are available for harmonic as well as for biharmonic diffusion.
  !!
  !! Liteature: Fox-Kemper, Menemenlis, Can Large Eddy Simulations improve Mesoscale Rich Ocean Models ?
  !!
  SUBROUTINE calculate_LeithClosure_harmonicDivGrad_VortDiv(patch_3d, ocean_state, param, operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: operators_coeff

    !Local variables
    INTEGER :: jk, blockNo, je, jc,jb, i
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
!     INTEGER :: LEITH_EXPONENT
    TYPE(t_subset_range), POINTER ::edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: grad_w_e, LeithCoeff
    REAL(wp):: laplacian_vn_out(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain

    start_level = 1
   !-------------------------------------------------------------------------------
   !leith closure for Laplacian(harmonic) viscosity
    CALL veloc_diff_harmonic_div_grad( patch_3D, operators_coeff%grad_coeff, ocean_state%p_diag, &
      & operators_coeff, laplacian_vn_out)

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level,end_level, grad_w_e, &
!ICON_OMP cell1_idx, cell1_blk, cell2_idx, cell2_blk, LeithCoeff, i) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index
        end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

        cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
        cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
        cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
        cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

        DO level=start_level,end_level


          grad_w_e = &
            & (  ocean_state%p_diag%w(cell2_idx,level,  cell2_blk) &
            &  + ocean_state%p_diag%w(cell2_idx,level+1,cell2_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level,  cell1_blk) &
            &  - ocean_state%p_diag%w(cell1_idx,level+1,cell1_blk)) &
            &    * 0.5_wp * patch_2D%edges%inv_dual_edge_length(je,blockNo)

          LeithCoeff = param%LeithHarmonicViscosity_BasisCoeff(je,blockNo) *  &
            &   sqrt(laplacian_vn_out(je,level,blockNo)**2  + grad_w_e**2)

          param%HarmonicViscosity_coeff(je,level,blockNo)  = &
            & param%HarmonicViscosity_BasisCoeff(je,blockNo) + &
            & LeithCoeff

!           param%TracerDiffusion_coeff(je,level,blockNo,1) = grad_vort_abs
!           param%TracerDiffusion_coeff(je,level,blockNo,2) = div_e

          DO i=1,no_tracer
            param%TracerDiffusion_coeff(je,level,blockNo,i) = &
              & param%TracerDiffusion_BasisCoeff(je,blockNo,i) + &
              & MIN(LeithCoeff * TracerDiffusion_LeithWeight,    &
              &     param%TracerDiffusion_BasisCoeff(je,blockNo,i) * max_turbulenece_TracerDiffusion)
          END DO

        END DO
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO
   ! this sync most probably is not needed
!    CALL sync_patch_array(sync_e, patch_2D, param%HarmonicViscosity_coeff)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('LeithClosure: viscosity',param%HarmonicViscosity_coeff,&
      & str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('LeithClosure: grad_vort_abs',param%TracerDiffusion_coeff(:,:,:,1),&
!       &str_module,idt_src, in_subset=edges_in_domain)
!     CALL dbg_print('LeithClosure: div_e',param%TracerDiffusion_coeff(:,:,:,2),&
!       &str_module,idt_src, in_subset=edges_in_domain)
    IF (TracerDiffusion_LeithWeight > 0.0_wp) THEN
!       DO i=1,no_tracer
        CALL dbg_print('Tracer diff BasisCoeff',param%TracerDiffusion_BasisCoeff(:,:,1),&
          &str_module,idt_src, in_subset=edges_in_domain)
        CALL dbg_print('LeithClosure: tracer diff',param%TracerDiffusion_coeff(:,:,:,1),&
          &str_module,idt_src, in_subset=edges_in_domain)
!       END DO
    ENDIF
!   !---------------------------------------------------------------------
    ! smooth if requested
    DO i=1, LeithViscosity_SmoothIterations
      CALL smooth_lapl_diff_3D( patch_3d, &
        & param%HarmonicViscosity_coeff, LeithViscosity_SpatialSmoothFactor )
    ENDDO

  END SUBROUTINE calculate_LeithClosure_harmonicDivGrad_VortDiv
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the viscosity coefficients following the Leith closure.
  !! Implemented is the pure leith closure and a modified version of it. Both options
  !! are available for harmonic as well as for biharmonic diffusion.
  !!
  !! Liteature: Fox-Kemper, Menemenlis, Can Large Eddy Simulations improve Mesoscale Rich Ocean Models ?
  !!
  !!
  SUBROUTINE calculate_LeithClosure_biharmonic(patch_3d, ocean_state, param, operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: operators_coeff

    !Local variables
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    INTEGER :: jk, blockNo, je, jc,jb, i
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk, vertex3_idx, vertex3_blk
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
!     INTEGER :: LEITH_EXPONENT
    REAL(wp):: div_e(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp):: div_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: vort_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: vort_e(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp):: laplacian_vort(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: laplacian_div(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: laplacian_vort_e, laplacian_div_e
    !REAL(wp):: size_grad_S_horz_vec(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain

    idt_src=1  ! output print level (1-5, fix)
    start_level = 1

!     div_e        (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
!     div_c        (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     vort_e        (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
!     vort_c       (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     laplacian_vort(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     laplacian_div(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
   !-------------------------------------------------------------------------------

!      LEITH_EXPONENT=6

    !1) calculation of vertical vorticity at cell centers

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, jc, level, vertex1_idx, &
!ICON_OMP vertex1_blk, vertex2_idx, vertex2_blk, vertex3_idx, vertex3_blk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      DO jc=start_cell_index, end_cell_index


        vertex1_idx = patch_2d%cells%vertex_idx(jc,blockNo,1)
        vertex1_blk = patch_2d%cells%vertex_blk(jc,blockNo,1)
        vertex2_idx = patch_2d%cells%vertex_idx(jc,blockNo,2)
        vertex2_blk = patch_2d%cells%vertex_blk(jc,blockNo,2)
        vertex3_idx = patch_2d%cells%vertex_idx(jc,blockNo,3)
        vertex3_blk = patch_2d%cells%vertex_blk(jc,blockNo,3)

        DO level=start_level, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)

          vort_c(jc,level,blockNo) = &
            &  ocean_state%p_diag%vort(vertex1_idx,level,vertex1_blk) &
            &+ ocean_state%p_diag%vort(vertex2_idx,level,vertex2_blk) &
            &+ ocean_state%p_diag%vort(vertex3_idx,level,vertex3_blk)

        END DO
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

!     CALL sync_patch_array(sync_c, patch_2D, vort_c)
!     CALL dbg_print('LeithBiharm:vort_c',vort_c,&
!       & str_module,idt_src, in_subset=cells_in_domain)

    !2) calculate laplacian of vertical velocity
    CALL tracer_diffusion_horz_local(patch_3D, vort_c, ocean_state, vort_e)!,subset_range = edges_in_domain)
    CALL sync_patch_array(sync_e, patch_2D, vort_e)

!     CALL dbg_print('LeithBiharm:vort_e', vort_e,&
!       & str_module,idt_src, in_subset=edges_in_domain)

    CALL div_oce_3d( vort_e, patch_3D, operators_coeff%div_coeff, laplacian_vort)
    CALL sync_patch_array(sync_c, patch_2D, laplacian_vort)
!     CALL dbg_print('LeithBiharm:laplacian_vort',laplacian_vort,&
!       & str_module,idt_src, in_subset=cells_in_domain)


   !3a) In case of pure Leith, we have all to calculate the coefficient
    IF(LeithClosure_form==1)THEN
     !Now aggregate the final parameter
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, level, cell1_idx, &
!ICON_OMP cell1_blk, cell2_idx, cell2_blk, laplacian_vort_e) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je=start_edge_index, end_edge_index
!           end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

          cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
          cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
          cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
          cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)
          DO level=start_level, patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

            laplacian_vort_e &
              & = 0.5_wp*(laplacian_vort(cell1_idx,level,cell1_blk) + laplacian_vort(cell2_idx,level,cell2_blk))
!             laplacian_vort_e = laplacian_vort_e*laplacian_vort_e

            !Add the leith viscosity coefficient
            param%BiharmonicViscosity_coeff(je,level,blockNo) = &
              & param%BiharmonicViscosity_BasisCoeff(je,blockNo)  &
              & + param%LeithBiharmonicViscosity_BasisCoeff(je,blockNo) * ABS(laplacian_vort_e)

          END DO
        END DO
      END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

!    DO level= start_level, 5!end_level-1
!     write(*,*)'SIZES Leith-Viscosity',level,maxval(param%HarmonicViscosity_BasisCoeff(:,level,:)),&
!     & minval(param%HarmonicViscosity_BasisCoeff(:,level,:)), maxval(laplacian_vort(:,level,:)),minval(laplacian_vort(:,level,:)),&
!      &maxval(length_scale)**LEITH_EXPONENT,minval(length_scale)**LEITH_EXPONENT
!     END DO

    !3b) In case of modified Leith, we need additionally the Laplacian of divergence
    ELSEIF(LeithClosure_form==2)THEN

      CALL div_oce_3d(ocean_state%p_diag%vn_time_weighted, patch_3D, operators_coeff%div_coeff, div_c)
      CALL sync_patch_array(sync_c, patch_2D, div_c)
!       CALL dbg_print('LeithBiharm:div_c',div_c,&
!         & str_module,idt_src, in_subset=cells_in_domain)
      !CALL div_oce_3d( ocean_state%p_diag%ptp_vn, patch_3D, operators_coeff%div_coeff, div_c)
      !  The next two calls calculate the laplacian of the divergence without any diffusion parameter
      CALL tracer_diffusion_horz_local(patch_3D, div_c, ocean_state, div_e)!,subset_range = edges_in_domain)
      CALL sync_patch_array(sync_e, patch_2D, div_e)
!       CALL dbg_print('LeithBiharm:div_e',div_e,&
!         & str_module,idt_src, in_subset=edges_in_domain)
      CALL div_oce_3d( div_e, patch_3D, operators_coeff%div_coeff, laplacian_div)
      CALL sync_patch_array(sync_c, patch_2D, laplacian_div)
!       CALL dbg_print('LeithBiharm:laplacian_div',laplacian_div,&
!         & str_module,idt_src, in_subset=cells_in_domain)

     !Now aggregate the final parameter
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, level, cell1_idx, &
!ICON_OMP cell1_blk, cell2_idx, cell2_blk, laplacian_div_e, laplacian_vort_e) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je=start_edge_index, end_edge_index
!           end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

          cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
          cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
          cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
          cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

          DO level=start_level,  patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

            ! maybe we need to weight them
            laplacian_div_e = &
              & 0.5_wp*(laplacian_div(cell1_idx,level,cell1_blk) + laplacian_div(cell2_idx,level,cell2_blk))
!             laplacian_div_e = laplacian_div_e * laplacian_div_e

            laplacian_vort_e = &
              & 0.5_wp*(laplacian_vort(cell1_idx,level,cell1_blk) + laplacian_vort(cell2_idx,level,cell2_blk))
!             laplacian_vort_e = laplacian_vort_e * laplacian_vort_e

            !The modified leith viscosity coefficient
            param%BiharmonicViscosity_coeff(je,level,blockNo) = &
              & param%BiharmonicViscosity_BasisCoeff(je,blockNo)  &
              & + param%LeithBiharmonicViscosity_BasisCoeff(je,blockNo) *  &
              & sqrt(laplacian_div_e**2 + laplacian_vort_e**2)
          END DO
        END DO
      END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

!    DO level= start_level, 5!end_level-1
!     write(*,*)'SIZES Leith-Viscosity',level,maxval(param%HarmonicViscosity_BasisCoeff(:,level,:)),&
!     & minval(param%HarmonicViscosity_BasisCoeff(:,level,:)), &
!     &maxval(laplacian_vort(:,level,:)),minval(laplacian_vort(:,level,:)),&
!     &maxval(laplacian_div(:,level,:)),minval(laplacian_div(:,level,:)),&
!     &maxval(length_scale)**LEITH_EXPONENT,minval(length_scale)**LEITH_EXPONENT
!     END DO

   ENDIF

   !CALL sync_patch_array(sync_e, patch_2D, grad_vort_abs)
   !Note: this sync most probably is nor needed

!idt_src=1  ! output print level (1-5, fix)
!      CALL debug_print_MaxMinMean('after ocean_gmres: h-new', minmaxmean, str_module, idt_src)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
     idt_src=1  ! output print level (1-5, fix)
     CALL dbg_print('LeithClosure_type: biharm visc',param%BiharmonicViscosity_coeff,&
       &str_module,idt_src, in_subset=edges_in_domain)
!   !---------------------------------------------------------------------

!    DO level= start_level, 5!end_level-1
!    write(*,*)'Leith-Viscosity',level,maxval(param%HarmonicViscosity_BasisCoeff(:,level,:)),&
!    & minval(param%HarmonicViscosity_BasisCoeff(:,level,:))!, maxval(grad_vort_abs(:,level,:)),minval(grad_vort_abs(:,level,:))!,&
!    !&maxval(grad_vort_abs(:,level,:)),minval(grad_vort_abs(:,level,:))
!    END DO
    ! smooth if requested
    DO i=1, LeithViscosity_SmoothIterations
      CALL smooth_lapl_diff_3D( patch_3d, &
        & param%BiharmonicViscosity_coeff, LeithViscosity_SpatialSmoothFactor )
    ENDDO

  END SUBROUTINE calculate_LeithClosure_biharmonic
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates importan physical numbers: Richardson number, Buoancy frequency, baroclinic wave speed, Rossby radius.
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_characteristic_physical_numbers(patch_3d, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    !TYPE(t_ho_params),                 INTENT(inout) :: param
    !TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff

    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level, blockNo
    !INTEGER :: start_edge_index, end_edge_index, je
    TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: lambda
    REAL(wp) :: depth_scale, depth
    REAL(wp) :: Coriolis_abs,z_grav_rho
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
    REAL(wp), PARAMETER :: beta_param=2.29E-11
    REAL(wp), PARAMETER :: latitude_threshold=5.0_wp
    REAL(wp) :: Rossby_radius_equator(nproma)
    REAL(wp) :: Rossby_radius_offequator(nproma)


    REAL(wp) :: z_shear_cell
    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)

    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    start_level=1
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    !-------------------------------------------------------------------------------
    z_grav_rho = grav/OceanReferenceDensity

!ICON_OMP_PARALLEL PRIVATE(salinity)
    salinity(1:n_zlev) = sal_ref
!ICON_OMP_DO PRIVATE(start_cell_index, end_cell_index, cell_index, end_level, level, pressure, z_rho_up, z_rho_down, &
!ICON_OMP z_shear_cell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      ocean_state%p_diag%Richardson_Number(:,:, blockNo) = 0.0_wp
      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
      DO cell_index = start_cell_index, end_cell_index

        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        IF (end_level < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:end_level) = ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        pressure(2:end_level) = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, 2:end_level, blockNo)&
        & * OceanReferenceDensity * sitodbar

        z_rho_up(1:end_level-1)  &
        &= calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level-1,blockNo,1), &
                                    & salinity(1:end_level-1), pressure(2:end_level),end_level-1)

        z_rho_down(2:end_level) &
        &= calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,2:end_level,blockNo,1), &
                                    & salinity(2:end_level), pressure(2:end_level), end_level-1)


        DO level = 2, end_level

          z_shear_cell = dbl_eps + &
            & SUM((ocean_state%p_diag%p_vn(cell_index,level-1,blockNo)%x - ocean_state%p_diag%p_vn(cell_index,level,blockNo)%x)**2)

          z_vert_density_grad_c(cell_index,level,blockNo) = (z_rho_down(level) - z_rho_up(level-1)) *  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)

          !adjusted vertical derivative (follows MOM, see Griffies-book, (p. 332, eq. (15.15)) or MOM-5 manual (sect. 23.7.1.1)
          !z_vert_density_grad_c(cell_index,level,blockNo)=min(z_vert_density_grad_c(cell_index,level,blockNo),-dbl_eps)
          ocean_state%p_diag%Richardson_Number(cell_index, level, blockNo) &
          &= MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(cell_index,level,blockNo) * z_grav_rho * &
          & (z_rho_down(level) - z_rho_up(level-1)) / z_shear_cell, 0.0_wp)
        END DO ! levels
        ocean_state%p_diag%grad_rho_PP_vert(cell_index,2:end_level,blockNo)=z_vert_density_grad_c(cell_index,2:end_level,blockNo)

      END DO ! index
    END DO
!ICON_OMP_END_DO


!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, &
!ICON_OMP Rossby_radius_equator, Rossby_radius_offequator) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      Rossby_radius_equator   (:)               =0.0_wp
      Rossby_radius_offequator(:)               =0.0_wp
      ocean_state%p_diag%Wavespeed_baroclinic(:,blockNo)=0.0_wp

      DO cell_index = start_cell_index, end_cell_index

          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          DO level = start_level, end_level

            !calculate buoyancy frequency
            ocean_state%p_diag%Buoyancy_Freq(cell_index,level,blockNo)&
            & = -z_grav_rho*z_vert_density_grad_c(cell_index,level,blockNo)

            !calculate baroclinic wave speed
            ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)&
            &=ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)&
            &+(1.0_wp/pi)*ocean_state%p_diag%Buoyancy_Freq(cell_index,level,blockNo)&
            &*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(cell_index,level,blockNo)

          END DO

          IF(abs(patch_2d%cells%center(cell_index, blockNo)%lat)*rad2deg>=latitude_threshold)THEN
            Rossby_radius_offequator(cell_index)&
            &=ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)&
            &/abs(patch_2d%cells%f_c(cell_index,blockNo))
          ENDIF
          Rossby_radius_equator(cell_index)&
            &=ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)/(2.0_wp*beta_param)
          ocean_state%p_diag%Rossby_Radius(cell_index,blockNo)&
            &=min(Rossby_radius_offequator(cell_index),Rossby_radius_equator(cell_index))
      END DO
    END DO
!ICON_OMP_END_DO


!ICON_OMP_END_PARALLEL

!   Do level=1,n_zlev
!   write(0,*)'max-min',level,&
!    &maxval( ocean_state%p_diag%Richardson_Number(:,level,:)),&
!    &minval( ocean_state%p_diag%Richardson_Number(:,level,:)),&
!    &maxval( ocean_state%p_diag%Buoyancy_Freq(:,level,:)),&
!    &minval( ocean_state%p_diag%Buoyancy_Freq(:,level,:))
!   End do
!   write(0,*)'max-min',&
!    &maxval( ocean_state%p_diag%Wavespeed_baroclinic),&
!    &minval( ocean_state%p_diag%Wavespeed_baroclinic),&
!    &maxval( ocean_state%p_diag%Rossby_Radius),&
!    &maxval(Rossby_radius_offequator), maxval(Rossby_radius_equator)

!  stop

  END SUBROUTINE calc_characteristic_physical_numbers
  !-------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE calc_vertical_stability(patch_3d, ocean_state, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    LOGICAL, INTENT(in), OPTIONAL                    :: lacc

    !Local variables
    INTEGER :: start_index, end_index, cell_index,level,end_level, blockNo
    LOGICAL :: lzacc

    TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D

    REAL(wp) :: z_grav_rho
    REAL(wp) :: z_shear_cell
    CHARACTER(len=*), PARAMETER :: routine = 'calc_vertical_stability'

#if defined(__LVECTOR__) && !defined(__LVEC_BITID__)
    REAL(wp) :: z_rho_up(nproma,n_zlev), z_rho_down(nproma,n_zlev) !, density(n_zlev)
    REAL(wp) :: pressure(nproma,n_zlev), salinity(nproma,n_zlev)!

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish('calc_vertical_stability', 'OpenACC version for LVECTOR currently not implemented')
#endif

    IF (eos_type /= 2) THEN
     write(0,*) "Vector version for eos_type =",eos_type," not yet implemented."
     stop
    ENDIF
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    !-------------------------------------------------------------------------------

    z_grav_rho = grav/OceanReferenceDensity

    !$ACC DATA CREATE(z_rho_up, pressure, z_rho_down, salinity) IF(lzacc)

    !ICON_OMP_PARALLEL PRIVATE(salinity, z_rho_up, z_rho_down, pressure)
    salinity = sal_ref
    z_rho_up=0.0_wp
    z_rho_down=0.0_wp
    pressure = 0._wp

    !ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, end_level, level, &
    !ICON_OMP z_shear_cell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      ocean_state%p_diag%Richardson_Number(:, :, blockNo) = 0.0_wp
      ocean_state%p_diag%zgrad_rho(:,:, blockNo) = 0.0_wp

      DO level = 2, MAXVAL(patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,blockNo))
        DO cell_index = start_index, end_index

          end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
          IF (end_level < 2 .OR. level > end_level) CYCLE

          pressure(cell_index,level) = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, level, blockNo) &
             * OceanReferenceDensity * sitodbar
        ENDDO
      ENDDO

      IF(no_tracer >= 2) THEN
        DO level = 1, MAXVAL(patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,blockNo))
          DO cell_index = start_index, end_index

            end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
            IF (end_level < 2 .OR. level > end_level) CYCLE

            salinity(cell_index,level) = ocean_state%p_prog(nold(1))%tracer(cell_index,level,blockNo,2)
          ENDDO
        ENDDO
      ENDIF
      DO level = 2, MAXVAL(patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,blockNo))
        DO cell_index = start_index, end_index

          end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
          IF (end_level < 2 .OR. level > end_level-1) CYCLE
          z_rho_up(cell_index,level) = &
             calculate_density_mpiom_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,level,blockNo,1), &
             salinity(cell_index,level), pressure(cell_index,level+1) )

          z_rho_down(cell_index,level) = &
             calculate_density_mpiom_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,level,blockNo,1), &
             salinity(cell_index,level), pressure(cell_index,level) )
        ENDDO
      ENDDO
      DO cell_index = start_index, end_index
        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        IF (end_level < 2 ) CYCLE
        z_rho_up(cell_index,1) = &
             calculate_density_mpiom_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,1,blockNo,1), &
             salinity(cell_index,1), pressure(cell_index,2) )

        z_rho_down(cell_index,end_level) = &
             calculate_density_mpiom_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,end_level,blockNo,1), &
             salinity(cell_index,end_level), pressure(cell_index,end_level) )
      ENDDO

      DO level = 2, MAXVAL(patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,blockNo))
        DO cell_index = start_index, end_index
          end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
          IF (end_level < 2 .OR. level > end_level) CYCLE

          z_shear_cell = dbl_eps + &
               SUM((ocean_state%p_diag%p_vn(cell_index,level-1,blockNo)%x &
               - ocean_state%p_diag%p_vn(cell_index,level,blockNo)%x)**2)

          ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo) = (z_rho_down(cell_index,level) - z_rho_up(cell_index,level-1)) *  &
               patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)

          !adjusted vertical derivative (follows MOM, see Griffies-book,
          ! (p. 332, eq. (15.15)) or MOM-5 manual (sect. 23.7.1.1)
          !ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo)= &
          !     MIN(ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo),-dbl_eps)

          ocean_state%p_diag%Richardson_Number(cell_index, level, blockNo) &
               = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(cell_index,level,blockNo) * z_grav_rho * &
               (z_rho_down(cell_index,level) - z_rho_up(cell_index,level-1)) / z_shear_cell, 0.0_wp)
        END DO ! index
      END DO ! levels
#else
    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev) !, density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)

    CALL set_acc_host_or_device(lzacc, lacc)

    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    !-------------------------------------------------------------------------------

    z_grav_rho = grav/OceanReferenceDensity

    !$ACC DATA CREATE(salinity, z_rho_up, z_rho_down, pressure) IF(lzacc)

    !ICON_OMP_PARALLEL PRIVATE(salinity, z_rho_up, z_rho_down, pressure)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    salinity(1:n_zlev) = sal_ref
    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp
    pressure(:) = 0._wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

    !ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, end_level, level, &
    !ICON_OMP z_shear_cell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      ocean_state%p_diag%Richardson_Number(:, :, blockNo) = 0.0_wp
      ocean_state%p_diag%zgrad_rho(:,:, blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC   PRIVATE(salinity, z_rho_up, z_rho_down, pressure) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO cell_index = start_index, end_index

        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        IF (end_level < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:end_level) = ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        pressure(2:end_level) = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, 2:end_level, blockNo) &
             * OceanReferenceDensity * sitodbar

#ifdef _OPENACC
        DO level = 1, end_level-1
          z_rho_up(level) = calculate_density_onColumn_elem(ocean_state%p_prog(nold(1))%tracer(cell_index,level,blockNo,1), &
                                              salinity(level), pressure(level+1))
        END DO
#else
        z_rho_up(1:end_level-1) = &
             calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level-1,blockNo,1), &
             salinity(1:end_level-1), pressure(2:end_level),end_level-1)
#endif

#ifdef _OPENACC
        DO level = 2, end_level
          z_rho_down(level) = calculate_density_onColumn_elem(ocean_state%p_prog(nold(1))%tracer(cell_index,level,blockNo,1), &
                                              salinity(level), pressure(level))
        END DO
#else
        z_rho_down(2:end_level) = &
             calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,2:end_level,blockNo,1), &
             salinity(2:end_level), pressure(2:end_level), end_level-1)
#endif

        DO level = 2, end_level

          z_shear_cell = dbl_eps + &
               SUM((ocean_state%p_diag%p_vn(cell_index,level-1,blockNo)%x &
               - ocean_state%p_diag%p_vn(cell_index,level,blockNo)%x)**2)

          ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo) = (z_rho_down(level) - z_rho_up(level-1)) *  &
               patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)

          !adjusted vertical derivative (follows MOM, see Griffies-book,
          ! (p. 332, eq. (15.15)) or MOM-5 manual (sect. 23.7.1.1)
          !ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo)= &
          !     MIN(ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo),-dbl_eps)

          ocean_state%p_diag%Richardson_Number(cell_index, level, blockNo) &
               = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(cell_index,level,blockNo) * z_grav_rho * &
               (z_rho_down(level) - z_rho_up(level-1)) / z_shear_cell, 0.0_wp)
        END DO ! levels
      END DO ! index
      !$ACC END PARALLEL LOOP
#endif
    END DO
    !$ACC WAIT(1)
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    !$ACC END DATA

  END SUBROUTINE calc_vertical_stability
  !-------------------------------------------------------------------------


!<Optimize:inUse>
#ifdef __LVECTOR__
  SUBROUTINE calc_vertical_stability_zstar(patch_3d, ocean_state, eta_c, stretch_c, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    REAL(wp), INTENT(IN) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL                    :: lacc

    !Local variables
    INTEGER :: start_index, end_index, cell_index,level,end_level, blockNo

    TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D

    REAL(wp) :: z_grav_rho
    REAL(wp) :: z_shear_cell
    REAL(wp) :: z_rho_up, z_rho_down !, density(n_zlev)
    REAL(wp) :: pressure, salinity_up, salinity_down
    LOGICAL :: lzacc

    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    !-------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish('calc_vertical_stability_zstar', 'OpenACC version for LVECTOR currently not implemented')
#endif

    z_grav_rho = grav/OceanReferenceDensity

    !ICON_OMP_PARALLEL PRIVATE(salinity_up, salinity_down, z_rho_up, z_rho_down, pressure)
    !ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, end_level, level, &
    !ICON_OMP z_shear_cell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      ocean_state%p_diag%Richardson_Number(:, :, blockNo) = 0.0_wp
      ocean_state%p_diag%zgrad_rho(:,:, blockNo) = 0.0_wp

      DO level = 2, MAXVAL(patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,blockNo))
        DO cell_index = start_index, end_index

          IF (level <= patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)) THEN

            IF(no_tracer >= 2) THEN
              salinity_up = ocean_state%p_prog(nold(1))%tracer(cell_index,level-1,blockNo,2)
              salinity_down = ocean_state%p_prog(nold(1))%tracer(cell_index,level,blockNo,2)
            ELSE
              salinity_up = sal_ref
              salinity_down = sal_ref
            ENDIF

            !--------------------------------------------------------
            pressure = (patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, level, blockNo) &
              & * stretch_c(cell_index, blockNo) - eta_c(cell_index,blockNo))  &
              & * OceanReferenceDensity * sitodbar

            z_rho_up = &
                calculate_density_onColumn_elem(ocean_state%p_prog(nold(1))%tracer(cell_index,level-1,blockNo,1), &
                salinity_up, pressure)

            z_rho_down = &
                calculate_density_onColumn_elem(ocean_state%p_prog(nold(1))%tracer(cell_index,level,blockNo,1), &
                salinity_down, pressure)

            z_shear_cell = dbl_eps + &
                SUM((ocean_state%p_diag%p_vn(cell_index,level-1,blockNo)%x &
                - ocean_state%p_diag%p_vn(cell_index,level,blockNo)%x)**2)

            ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo) = (z_rho_down - z_rho_up) *  &
                patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo) / stretch_c(cell_index, blockNo)

            !adjusted vertical derivative (follows MOM, see Griffies-book,
            ! (p. 332, eq. (15.15)) or MOM-5 manual (sect. 23.7.1.1)
            !ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo)= &
            !     MIN(ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo),-dbl_eps)

            ocean_state%p_diag%Richardson_Number(cell_index, level, blockNo) &
                = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(cell_index,level,blockNo) * stretch_c(cell_index,blockNo) * z_grav_rho * &
                (z_rho_down - z_rho_up) / z_shear_cell, 0.0_wp)
          END IF
        END DO ! index
      END DO ! levels
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

  END SUBROUTINE calc_vertical_stability_zstar
#else
  SUBROUTINE calc_vertical_stability_zstar(patch_3d, ocean_state, eta_c, stretch_c, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    REAL(wp), INTENT(IN) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL                    :: lacc

    !Local variables
    INTEGER :: start_index, end_index, cell_index,level,end_level, blockNo

    TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D

    REAL(wp) :: z_grav_rho
    REAL(wp) :: z_shear_cell
    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev) !, density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    LOGICAL :: lzacc

    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    !-------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish('calc_vertical_stability_zstar', 'OpenACC version currently not implemented')
#endif

    z_grav_rho = grav/OceanReferenceDensity

    !ICON_OMP_PARALLEL PRIVATE(salinity, z_rho_up, z_rho_down, pressure)
    salinity(1:n_zlev) = sal_ref
    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp
    pressure(:) = 0._wp

    !ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, end_level, level, &
    !ICON_OMP z_shear_cell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      ocean_state%p_diag%Richardson_Number(:, :, blockNo) = 0.0_wp
      ocean_state%p_diag%zgrad_rho(:,:, blockNo) = 0.0_wp

      DO cell_index = start_index, end_index

        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        IF (end_level < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:end_level) = ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        pressure(2:end_level) = (patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, 2:end_level, blockNo) &
           & * stretch_c(cell_index, blockNo) - eta_c(cell_index,blockNo))  &
           & * OceanReferenceDensity * sitodbar

        z_rho_up(1:end_level-1) = &
             calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level-1,blockNo,1), &
             salinity(1:end_level-1), pressure(2:end_level),end_level-1)

        z_rho_down(2:end_level) = &
             calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,2:end_level,blockNo,1), &
             salinity(2:end_level), pressure(2:end_level), end_level-1)

        DO level = 2, end_level

          z_shear_cell = dbl_eps + &
               SUM((ocean_state%p_diag%p_vn(cell_index,level-1,blockNo)%x &
               - ocean_state%p_diag%p_vn(cell_index,level,blockNo)%x)**2)

          ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo) = (z_rho_down(level) - z_rho_up(level-1)) *  &
               patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo) / stretch_c(cell_index, blockNo)

          !adjusted vertical derivative (follows MOM, see Griffies-book,
          ! (p. 332, eq. (15.15)) or MOM-5 manual (sect. 23.7.1.1)
          !ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo)= &
          !     MIN(ocean_state%p_diag%zgrad_rho(cell_index,level,blockNo),-dbl_eps)

          ocean_state%p_diag%Richardson_Number(cell_index, level, blockNo) &
               = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(cell_index,level,blockNo) * stretch_c(cell_index,blockNo) * z_grav_rho * &
               (z_rho_down(level) - z_rho_up(level-1)) / z_shear_cell, 0.0_wp)
        END DO ! levels
      END DO ! index
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

  END SUBROUTINE calc_vertical_stability_zstar
#endif
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
  SUBROUTINE tracer_diffusion_horz_local(patch_3D, trac_in, p_os, diff_flx, k_t, subset_range)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp), OPTIONAL                :: k_t(:,:,:) !mixing coefficient for tracer
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    !
    !Local variables
    INTEGER :: level, blockNo, edge_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2D%cells%in_domain
    !-------------------------------------------------------------------------------

    IF(PRESENT(k_t))THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        diff_flx(:,:,blockNo) = 0.0_wp
        DO edge_index = start_edge_index, end_edge_index
          !Get indices of two adjacent triangles
          il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

          DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

            diff_flx(edge_index,level,blockNo) = &
              &   k_t(edge_index,level,blockNo) &
              & * patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
              & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
              & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)

          ENDDO

        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
    ENDIF

    ELSEIF(.NOT.PRESENT(k_t))THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        diff_flx(:,:,blockNo) = 0.0_wp
        DO edge_index = start_edge_index, end_edge_index
          !Get indices of two adjacent triangles
          il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

          DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

            diff_flx(edge_index,level,blockNo) = &
              & patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
              & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
              & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)

          ENDDO

        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      IF (PRESENT(subset_range)) THEN
        IF (.NOT. subset_range%is_in_domain) &
          & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
      ENDIF

    ENDIF

!     CALL dbg_print('LeithDiff:trac_in',trac_in,&
!       & str_module,1, in_subset=cells_in_domain)
!     CALL dbg_print('LeithDiff:diff_flx',diff_flx,&
!       & str_module,1, in_subset=edges_in_domain)
!     CALL dbg_print('LeithDiff:thick_e',patch_3D%p_patch_1d(1)%prism_thick_e,&
!       & str_module,1, in_subset=edges_in_domain)
!     CALL dbg_print('LeithDiff:dual_edge_length',patch_2D%edges%inv_dual_edge_length,&
!       & str_module,1, in_subset=edges_in_domain)

  END SUBROUTINE tracer_diffusion_horz_local
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE divideBy(nominator, denominator, in_subset)
    REAL(wp) :: nominator(:,:) ! INTENT(inout)
    REAL(wp) :: denominator(:,:) ! INTENT(in)
    TYPE(t_subset_range), TARGET :: in_subset

    INTEGER :: block, start_index, end_index, idx

!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index, idx)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          nominator(idx, block) = nominator(idx, block) / denominator(idx, block)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE divideBy
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE copy2Dto3D (from, to, in_subset)
    REAL(wp) :: from(:,:) ! INTENT(in)
    REAL(wp) :: to(:,:,:) ! INTENT(in)
    TYPE(t_subset_range), TARGET :: in_subset

    INTEGER :: block, level, start_index, end_index, idx

!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index, idx, level)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = 1, in_subset%vertical_levels(idx,block)
            to(idx, level, block) = from(idx, block)
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE copy2Dto3D
  !-----------------------------------------------------------------------


END MODULE mo_ocean_physics
