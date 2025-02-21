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

! DSL definitions 












!---------------------
! block definitions




!---------------------
! mappings



















!---------------------
! connectivity









!---------------------
! generic types




!---------------------
! shortcuts







































!---------------------
! Upper-lower case
!----------------------------
MODULE mo_ocean_physics_types
  !-------------------------------------------------------------------------
  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: &
    & n_zlev, bottom_drag_coeff,                    &
    & no_tracer,                                    &
    & BiharmonicViscosity_scaling, HarmonicViscosity_scaling, &
    & VelocityDiffusion_order,                                &
    & BiharmonicViscosity_reference,                          &
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
    & OceanReferenceDensity,    &
    & vert_mix_type, vmix_pp, vmix_tke, vmix_idemix_tke,      &
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
    &  TracerDiffusion_LeithWeight,l_lc

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
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_dynamics_config,     ONLY: nold!, nnew
  USE mo_run_config,          ONLY: dtime
  USE mo_var_list,            ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_groups,          ONLY: groups
  USE mo_cf_convention
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: datatype_pack16, DATATYPE_FLT32, DATATYPE_FLT64, filetype_nc2, &
    &                               GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: &
    & za_depth_below_sea, za_depth_below_sea_half, za_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max, sync_patch_array_mult
  USE mo_io_config,           ONLY: lnetcdf_flt64_output

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




  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_physics_types'
!   CHARACTER(LEN=12)           :: str_module    = 'ocePhysics  '  ! Output of module for 1 line debug

  ! Public interface
  PUBLIC :: t_ho_params, ocean_params_list

  !PUBLIC :: init_ho_physics
  PUBLIC :: construct_ho_params
  PUBLIC :: destruct_ho_params

  ! variables
  TYPE (t_var_list_ptr) :: ocean_params_list

  TYPE t_vmix_params
    REAL(wp),POINTER ::            &
      ! start by_nils
      & vmix_dummy_1(:,:,:)      ,& !
      & vmix_dummy_2(:,:,:)      ,& !
      & vmix_dummy_3(:,:,:)      ,& !
      & tke(:,:,:)                ,& ! turbulent kinetic energy
      & tke_Tbpr(:,:,:)           ,& ! tke tend by bpr
      & tke_Tspr(:,:,:)           ,& ! tke tend by spr
      & tke_Tdif(:,:,:)           ,& ! tke tend by dif
      & tke_Tdis(:,:,:)           ,& ! tke tend by dis
      & tke_Twin(:,:,:)           ,& ! tke tend by win
      & tke_Tiwf(:,:,:)           ,& ! tke tend by iwf
      & tke_Tbck(:,:,:)           ,& ! tke tend by bck
      & tke_Ttot(:,:,:)           ,& ! tke tend by tot
      & tke_Lmix(:,:,:)           ,& ! tke mixing length
      & tke_Pr(:,:,:)             ,& ! tke Prandtl number
      & iwe(:,:,:)                ,& ! internal wave energy
      & iwe_Ttot(:,:,:)           ,& ! iwe tend by tot
      & iwe_Tdif(:,:,:)           ,& ! iwe tend by dif
      & iwe_Thdi(:,:,:)           ,& ! iwe tend by hdi
      & iwe_Tdis(:,:,:)           ,& ! iwe tend by dis
      & iwe_Tsur(:,:,:)           ,& ! iwe tend by sur
      & iwe_Tbot(:,:,:)           ,& ! iwe tend by bot
      & iwe_alpha_c(:,:,:)        ,& ! iwe dissip coef
      & iwe_c0(:,:,:)             ,& ! vertical group velocity
      & iwe_v0(:,:,:)             ,& ! horizontal group velocity
      & tke_plc(:,:,:)            ,& ! langmuir turbulence
      & u_stokes(:,:)             ,& ! Stokes drift (m/s)
      & wlc(:,:,:)                ,& ! Langmuir turbulence velocity scale
      & hlc(:,:)                     ! depth of langmuir cell (m)
      ! end by_nils
  END TYPE t_vmix_params


  ! Parameters below appear directly in the ocean model/equation. They are eventually
  ! dynamically updated by using the "ocean-physics" structure. #slo# - not yet
  TYPE t_ho_params

    ! diffusion coefficients for horizontal velocity, temp. and salinity, dim=(nproma,n_zlev,nblks_e)
    REAL(wp),POINTER ::     &
      & HarmonicViscosity_BasisCoeff(:,:),   & ! coefficient of horizontal velocity harmonic diffusion
      & BiharmonicViscosity_BasisCoeff(:,:), & ! coefficient of horizontal velocity biharmonic diffusion
      & LeithHarmonicViscosity_BasisCoeff(:,:), & ! coefficient of Leith scaled basis coefficients, use in calculating the HarmonicViscosity_coeff
      & LeithBiharmonicViscosity_BasisCoeff(:,:), & ! coefficient of Leith scaled basis coefficients, use in calculating the HarmonicViscosity_coeff
      & HarmonicViscosity_coeff(:,:,:),       & ! coefficient of total diffusion
      & BiharmonicViscosity_coeff(:,:,:),     & ! coefficient of total diffusion
      & TracerDiffusion_BasisCoeff(:,:,:),    &  ! coefficient of horizontal tracer diffusion
      & TracerDiffusion_coeff(:,:,:,:)  ! coefficient of horizontal tracer diffusion
!       & TracerDiffusion_coeff(:,:,:,:)  ! coefficient of horizontal tracer diffusion
    TYPE(t_onEdges_Pointer_3d_wp),ALLOCATABLE :: tracer_h_ptr(:)

    ! diffusion coefficients for vertical velocity, temp. and salinity, dim=(nproma,n_zlev+1,nblks_e)
    REAL(wp),POINTER ::     &
      & a_veloc_v(:,:,:),  & ! coefficient of vertical velocity diffusion
      & a_tracer_v(:,:,:,:)  ! coefficient of vertical tracer diffusion
    TYPE(t_onCells_HalfLevels_Pointer_wp), ALLOCATABLE :: tracer_v_ptr(:)

    !constant background values of coefficients above
    REAL(wp) :: &
      & a_veloc_v_back   ! coefficient of vertical velocity diffusion

    REAL(wp), POINTER, DIMENSION(:,:,:) :: tracer_windMixing
    REAL(wp), POINTER, DIMENSION(:,:,:) :: velocity_windMixing 

    REAL(wp),ALLOCATABLE ::         &
      & Tracer_HorizontalDiffusion_Reference(:),    &
      & Tracer_HorizontalDiffusion_Background(:),   & ! coefficient of horizontal tracer diffusion dim=no_tracer
      & a_tracer_v_back(:)       ! coefficient of vertical tracer diffusion dim=no_tracer

    REAL(wp),POINTER ::     &
      & k_tracer_isoneutral(:,:,:),  & ! coefficient of isoneutral tracer diffusion diffusion at cells
      & k_tracer_dianeutral(:,:,:)  ! coefficient of dianeutral tracer diffusion
   REAL(wp),POINTER ::     &
      & k_tracer_GM_kappa(:,:,:) ! coefficient of Gent-McWilliams mesoscale eddyparametrizations
    REAL(wp) :: bottom_drag_coeff

    ! by_nils
    TYPE(t_vmix_params) :: vmix_params
    ! by_nils

  END TYPE t_ho_params

  TYPE(t_ho_params),PUBLIC,TARGET :: v_params

  ! used in only for the winf mixing 
  REAL(wp), POINTER, PUBLIC :: WindMixingDecay(:), WindMixingLevel(:)

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Construction of arrays for ocean physics
  !!
  !! Construction of arrays for ocean physics ...
  !!
  !
  !
!<Optimize:inUse:initOnly>
  SUBROUTINE construct_ho_params(patch_2D, params_oce, ocean_restart_list)

    TYPE(t_patch),      INTENT(IN)    :: patch_2D
    TYPE (t_ho_params), INTENT(INOUT) :: params_oce
    TYPE (t_var_list_ptr),  INTENT(INOUT) :: ocean_restart_list

    ! Local variables
    INTEGER :: ist, i,jtrc
    INTEGER :: alloc_cell_blocks, nblks_e
    INTEGER :: datatype_flt

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = this_mod_name//':construct_ho_physics'

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'construct hydro ocean physics')

    CALL vlr_add(ocean_params_list, 'ocean_params_list', patch_id=patch_2D%id, &
      & lrestart=.FALSE., model_type=TRIM(get_my_process_name()))

    ! determine size of arrays
    alloc_cell_blocks = patch_2D%alloc_cell_blocks
    nblks_e = patch_2D%nblks_e

    !$ACC ENTER DATA COPYIN(params_oce, params_oce%vmix_params)

    IF (VelocityDiffusion_order == 1 .OR. VelocityDiffusion_order == 21) THEN
! cfname: difmxylo2d, m2 s-1, ocean_momentum_xy_laplacian_diffusivity
      CALL add_var(ocean_params_list, 'HarmonicViscosity_BasisCoeff', &
        & params_oce%HarmonicViscosity_BasisCoeff , grid_unstructured_edge,&
        & za_surface, &
        & t_cf_var('HarmonicViscosity_BasisCoeff', 'm2 s-1', 'ocean_momentum_xy_laplacian_diffusivity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,nblks_e/),in_group=groups("oce_physics"))

! cfname: difmxylo, m2 s-1, ocean_momentum_xy_laplacian_diffusivity
      CALL add_var(ocean_params_list, 'HarmonicViscosity_coeff', &
        & params_oce%HarmonicViscosity_coeff , grid_unstructured_edge,&
        & za_depth_below_sea, &
        & t_cf_var('HarmonicViscosity_coeff', 'm2 s-1', 'ocean_momentum_xy_laplacian_diffusivity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"), lopenacc=.TRUE.)
      
    ENDIF

    IF (VelocityDiffusion_order == 2 .OR. VelocityDiffusion_order == 21) THEN
!cfname: difmxybo2d, m4 s-1, ocean_momentum_xy_biharmonic_diffusivity
      CALL add_var(ocean_params_list, 'BiharmonicViscosity_BasisCoeff', &
        & params_oce%BiharmonicViscosity_BasisCoeff , grid_unstructured_edge,&
        & za_surface, &
        & t_cf_var('BiharmonicViscosity_BasisCoeff', 'm4 s-1', 'ocean_momentum_xy_biharmonic_diffusivity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,nblks_e/),in_group=groups("oce_physics"))

! cfname: difmxybo, m4 s-1, ocean_momentum_xy_biharmonic_diffusivity
      CALL add_var(ocean_params_list, 'BiharmonicViscosity_coeff', &
        & params_oce%BiharmonicViscosity_coeff , grid_unstructured_edge,&
        & za_depth_below_sea, &
        & t_cf_var('BiharmonicViscosity_coeff', 'm4 s-1', 'ocean_momentum_xy_biharmonic_diffusivity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))
   ENDIF
   
    IF (LeithClosure_order == 1 .or.  LeithClosure_order == 21) THEN
      CALL add_var(ocean_params_list, 'LeithHarmonicViscosity_BasisCoeff', &
        & params_oce%LeithHarmonicViscosity_BasisCoeff , grid_unstructured_edge,&
        & za_surface, &
        & t_cf_var('LeithHarmonicViscosity_BasisCoeff', 'm2 s-1', 'horizontal velocity diffusion', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,nblks_e/),in_group=groups("oce_physics"))
    ENDIF
    IF (LeithClosure_order == 2 .or.  LeithClosure_order == 21) THEN
      CALL add_var(ocean_params_list, 'LeithBiharmonicViscosity_BasisCoeff', &
        & params_oce%LeithBiharmonicViscosity_BasisCoeff , grid_unstructured_edge,&
        & za_surface, &
        & t_cf_var('LeithBiharmonicViscosity_BasisCoeff', 'm4 s-1', 'horizontal velocity diffusion', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,nblks_e/),in_group=groups("oce_physics"))
    ENDIF

    CALL add_var(ocean_params_list, 'A_veloc_v', params_oce%a_veloc_v , grid_unstructured_edge,&
      & za_depth_below_sea_half, &
      & t_cf_var('ocean_vertical_momentum_diffusivity', 'm2 s-1', 'vertical velocity diffusion', datatype_flt),&
      & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev+1,nblks_e/),in_group=groups("oce_physics","oce_diag"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'velocity_windMixing', params_oce%velocity_windMixing , grid_unstructured_edge,&
      & za_depth_below_sea_half, &
      & t_cf_var('velocity_windMixing', 'm2 s-1', 'vertical velocity windMixing', datatype_flt),&
      & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev+1,nblks_e/),in_group=groups("oce_physics","oce_diag"))

    ! start by_nils
    ! --- vmix dummy variables
    !IF (vert_mix_type .EQ. vmix_tke) THEN
    IF (vert_mix_type .EQ. vmix_tke .or. vert_mix_type .EQ.  vmix_idemix_tke) THEN
    CALL add_var(ocean_params_list, 'vmix_dummy_1', params_oce%vmix_params%vmix_dummy_1, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('vmix_dummy_1', 'fixme', 'vmix_dummy_1', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'vmix_dummy_2', params_oce%vmix_params%vmix_dummy_2, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('vmix_dummy_2', 'fixme', 'vmix_dummy_2', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'vmix_dummy_3', params_oce%vmix_params%vmix_dummy_3, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('vmix_dummy_3', 'fixme', 'vmix_dummy_3', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    ! --- TKE variables
    CALL add_var(ocean_restart_list, 'tke', params_oce%vmix_params%tke, &
        & grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('specific_turbulent_kinetic_energy_of_sea_water', 'm2 s-2', 'turbulent kinetic energy', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
        & lrestart_cont=.TRUE., in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    
 
    CALL add_var(ocean_params_list, 'tke_Tbpr', params_oce%vmix_params%tke_Tbpr, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Tbpr', 'm2 s-3', 'TKE tend bpr', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Tspr', params_oce%vmix_params%tke_Tspr, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Tspr', 'm2 s-3', 'TKE tend spr', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Tdif', params_oce%vmix_params%tke_Tdif, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Tdif', 'm2 s-3', 'TKE tend dif', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Tdis', params_oce%vmix_params%tke_Tdis, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Tdis', 'm2 s-3', 'TKE tend dis', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Twin', params_oce%vmix_params%tke_Twin, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Twin', 'm2 s-3', 'TKE tend win', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Tiwf', params_oce%vmix_params%tke_Tiwf, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Tiwf', 'm2 s-3', 'TKE tend iwf', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Tbck', params_oce%vmix_params%tke_Tbck, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Tbck', 'm2 s-3', 'TKE tend bck', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Ttot', params_oce%vmix_params%tke_Ttot, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Ttot', 'm2 s-3', 'TKE tend tot', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Lmix', params_oce%vmix_params%tke_Lmix, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Lmix', 'm', 'TKE mixing length', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    CALL add_var(ocean_params_list, 'tke_Pr', params_oce%vmix_params%tke_Pr, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('tke_Pr', '', 'TKE Prandtl number', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
    

    IF (l_lc) THEN
      CALL add_var(ocean_params_list, 'tke_plc', params_oce%vmix_params%tke_plc, &
         & grid_unstructured_cell, za_depth_below_sea_half, &
         & t_cf_var('tke_plc', 'm2 s-3', 'TKE langmuir turbulence', datatype_flt),&
         & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
         & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
         & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
      

      CALL add_var(ocean_params_list, 'wlc',params_oce%vmix_params%wlc, &
         & grid_unstructured_cell, za_depth_below_sea_half, &
         & t_cf_var('wlc', 'm s-1', 'Langmuir turbulence velocity scale', datatype_flt),&
         & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
         & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
         & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
      

      CALL add_var(ocean_params_list, 'hlc',params_oce%vmix_params%hlc, &
         & grid_unstructured_cell, za_surface, &
         & t_cf_var('hlc', 'm', 'Depth of langmuir cell', datatype_flt),&
         & grib2_var(255, 255, 255, datatype_pack16,GRID_UNSTRUCTURED,grid_cell),&
         & ldims=(/nproma,alloc_cell_blocks/), &
         & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
      

      CALL add_var(ocean_params_list, 'u_stokes',params_oce%vmix_params%u_stokes, &
         & grid_unstructured_cell, za_surface, &
         & t_cf_var('u_stokes', 'm s-1', 'Stokes drift', datatype_flt),&
         & grib2_var(255, 255, 255, datatype_pack16,GRID_UNSTRUCTURED,grid_cell),&
         & ldims=(/nproma,alloc_cell_blocks/), &
         & in_group=groups("oce_vmix_tke"), lopenacc=.TRUE.)
      

    ENDIF

    ENDIF

    ! --- IWE variables
    IF (vert_mix_type .EQ. vmix_idemix_tke) THEN
      CALL add_var(ocean_restart_list, 'iwe', params_oce%vmix_params%iwe, &
         & grid_unstructured_cell, za_depth_below_sea_half, &
         & t_cf_var('iwe', 'm2 s-2', 'internal wave energy', datatype_flt),&
         & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
         & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
         & lrestart_cont=.TRUE., in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_Ttot', params_oce%vmix_params%iwe_Ttot, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_Ttot', 'm2 s-3', 'IWE tend tot', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_Tdif', params_oce%vmix_params%iwe_Tdif, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_Tdif', 'm2 s-3', 'IWE tend dif', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_Thdi', params_oce%vmix_params%iwe_Thdi, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_Thdi', 'm2 s-3', 'IWE tend hdi', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_Tdis', params_oce%vmix_params%iwe_Tdis, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_Tdis', 'm2 s-3', 'IWE tend dis', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_Tsur', params_oce%vmix_params%iwe_Tsur, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_Tsur', 'm2 s-3', 'IWE tend sur', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_Tbot', params_oce%vmix_params%iwe_Tbot, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_Tbot', 'm2 s-3', 'IWE tend bot', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_alpha_c', params_oce%vmix_params%iwe_alpha_c, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_alpha_c', 's m-2', 'IWE dissip coef', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_c0', params_oce%vmix_params%iwe_c0, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_c0', 'm2 s-3', 'IWE ver group vel', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    CALL add_var(ocean_params_list, 'iwe_v0', params_oce%vmix_params%iwe_v0, &
       & grid_unstructured_cell, za_depth_below_sea_half, &
       & t_cf_var('iwe_v0', 'm2 s-3', 'IWE hor group vel', datatype_flt),&
       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
       & in_group=groups("oce_vmix_iwe"))

    ! end by_nils
    ENDIF

    !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL add_var(ocean_params_list, 'TracerDiffusion_coeff', params_oce%TracerDiffusion_coeff , &
        & grid_unstructured_edge, za_surface, &
        & t_cf_var('TracerDiffusion_coeff', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e,no_tracer/), &
        & lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE.)
      CALL add_var(ocean_params_list, 'TracerDiffusion_BasisCoeff', params_oce%TracerDiffusion_BasisCoeff , &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('TracerDiffusion_BasisCoeff', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,nblks_e,no_tracer/), &
        & lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE.)
      CALL add_var(ocean_params_list, 'A_tracer_v', params_oce%a_tracer_v , &
        & grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('A_tracer_v', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev+1,alloc_cell_blocks,no_tracer/), &
        & lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE., lopenacc=.TRUE.)
      

      CALL add_var(ocean_params_list, 'tracer_windMixing', params_oce%tracer_windMixing , &
        & grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('tracer_windMixing', '', 'tracer_windMixing', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
        & loutput=.TRUE., lrestart=.FALSE.,in_group=groups("oce_physics","oce_diag"))

      ! Reference to individual tracer, for I/O
      ALLOCATE(params_oce%tracer_h_ptr(no_tracer))
      ALLOCATE(params_oce%tracer_v_ptr(no_tracer))
      DO jtrc = 1,no_tracer
        CALL add_ref( ocean_params_list, 'TracerDiffusion_coeff',&
          & 'K_tracer_h_'//TRIM(oce_config%tracer_shortnames(jtrc)),     &
          & params_oce%tracer_h_ptr(jtrc)%p,                             &
          & grid_unstructured_edge, za_depth_below_sea,               &
          & t_cf_var('K_tracer_h_'//TRIM(oce_config%tracer_shortnames(jtrc)), &
          & 'kg/kg', &
          & TRIM(oce_config%tracer_longnames(jtrc))//'(K_tracer_h_)', &
          & datatype_flt), &
          & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
          & ref_idx=jtrc, &
          & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))
        CALL add_ref( ocean_params_list, 'A_tracer_v',&
          & 'A_tracer_v_'//TRIM(oce_config%tracer_shortnames(jtrc)),     &
          & params_oce%tracer_h_ptr(jtrc)%p,                             &
          & grid_unstructured_cell, za_depth_below_sea_half,            &
          & t_cf_var('A_tracer_v_'//TRIM(oce_config%tracer_shortnames(jtrc)), &
          & 'kg/kg', &
          & TRIM(oce_config%tracer_longnames(jtrc))//'(A_tracer_v)', &
          & datatype_flt), &
          & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
          & ref_idx=jtrc, &
          & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_physics"))

      END DO
      !TODO     use the following code, if add_var support 1d arrays:
      !TODO     CALL add_var(ocean_params_list, 'K_tracer_h_back', params_oce%K_tracer_h_back , &
      !TODO     &            GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
      !TODO     &            t_cf_var('K_tracer_h_back', '', '1:temperature 2:salinity', datatype_flt),&
      !TODO     &            grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_EDGE),&
      !TODO     &            ldims=(/ no_tracer /))
      !TODO     CALL add_var(ocean_params_list, 'A_tracer_v_back', params_oce%A_tracer_v_back , &
      !TODO     &            GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      !TODO     &            t_cf_var('A_tracer_v_back', '', '1:temperature 2:salinity', datatype_flt),&
      !TODO     &            grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      !TODO     &            ldims=(/no_tracer/))
    ENDIF ! no_tracer > 0


    ALLOCATE(params_oce%Tracer_HorizontalDiffusion_Background(no_tracer),  &
      & params_oce%Tracer_HorizontalDiffusion_Reference(no_tracer),                  &
      & params_oce%a_tracer_v_back(no_tracer),            &
      & stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for horizontal background tracer diffusion failed')
    END IF

   IF(GMRedi_configuration==GMRedi_combined &
   &.OR.GMRedi_configuration==GM_only.OR.GMRedi_configuration==Redi_only)THEN

     CALL add_var(ocean_params_list, 'k_tracer_isoneutral', params_oce%k_tracer_isoneutral, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('k_tracer_isoneutral at edges', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(ocean_params_list, 'k_tracer_dianeutral', params_oce%k_tracer_dianeutral, &
        & grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('k_tracer_dianeutral', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


     CALL add_var(ocean_params_list, 'k_tracer_GM_kappa', params_oce%k_tracer_GM_kappa, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('k_tracer_GM_kappa at cells', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ENDIF

    ALLOCATE(WindMixingDecay(1:n_zlev+1), WindMixingLevel(1:n_zlev+1))
    
  END SUBROUTINE construct_ho_params
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Destruction of arrays for ocean physics
  !!
  !
!<Optimize:inUse:initOnly>
  SUBROUTINE destruct_ho_params(params_oce)

    TYPE (t_ho_params), INTENT(inout) :: params_oce

    ! Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = this_mod_name//':destruct_ho_physics'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'destruct hydro ocean physics')

    CALL vlr_del(ocean_params_list)

    DEALLOCATE(params_oce%a_tracer_v_back,               &
      & params_oce%Tracer_HorizontalDiffusion_Reference, &
      & params_oce%Tracer_HorizontalDiffusion_Background,&
      & stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for tracer Diffusion Background failed')
    END IF
    
    DEALLOCATE(WindMixingDecay, WindMixingLevel)
    
  END SUBROUTINE destruct_ho_params
  !-------------------------------------------------------------------------


END MODULE mo_ocean_physics_types
