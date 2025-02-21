! Contains the ocean variables that hamocc uses
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
!=============================================================================================
MODULE mo_ocean_hamocc_couple_state

  USE mo_kind,                ONLY: wp, sp
  USE mo_model_domain,        ONLY: t_patch_3d, t_patch
  USE mo_ocean_tracer_transport_types, ONLY: t_ocean_transport_state
  USE mo_ocean_nml,           ONLY: n_zlev, TracerHorizontalDiffusion_scaling, &
    &  Salinity_HorizontalDiffusion_Background,  &
    &  Salinity_HorizontalDiffusion_Reference,   &
    &  TracerHorizontalDiffusion_scaling
  USE mo_impl_constants,      ONLY: success, max_char_length, TLEV_NNEW
  USE mo_parallel_config,     ONLY: nproma
  USE mo_var_groups,          ONLY: groups 
  USE mo_var_list,            ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_metadata,        ONLY: get_timelevel_string
  USE mo_grib2,               ONLY: grib2_var, t_grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32 => CDI_DATATYPE_FLT32, &
    &                               datatype_FLT64 => CDI_datatype_FLT64, &
    &                               DATATYPE_INT8 => CDI_DATATYPE_INT8, &
    &                               DATATYPE_PACK16 => CDI_DATATYPE_PACK16, &
    &                               tstep_constant, GRID_UNSTRUCTURED
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell, grid_unstructured_edge, &
      &                             grid_unstructured_vert, grid_vertex
  USE mo_zaxis_type,          ONLY: za_depth_below_sea, za_depth_below_sea_half, za_surface
  USE mo_ocean_physics,       ONLY: scale_horizontal_diffusion, copy2Dto3D
  USE mo_cf_convention
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, t_ocean_transport_state, t_hamocc_ocean_state
  PUBLIC :: hamocc_ocean_state, hamocc_ocean_state_list
  PUBLIC :: construct_hamocc_ocean_state, destruct_hamocc_ocean_state
  
  !----------------------------------------------
  TYPE t_ocean_to_hamocc_state
    REAL(wp), POINTER, DIMENSION(:,:) :: top_dilution_coeff
    REAL(wp), POINTER, DIMENSION(:,:) :: h_old
    REAL(wp), POINTER, DIMENSION(:,:) :: h_new
    REAL(wp), POINTER, DIMENSION(:,:) :: h_old_withIce
    REAL(wp), POINTER, DIMENSION(:,:) :: ice_concentration_sum
    
    REAL(wp),  POINTER, DIMENSION(:,:,:)    :: temperature
    REAL(wp),  POINTER, DIMENSION(:,:,:)    :: salinity
    REAL(wp),  POINTER, DIMENSION(:,:,:)    :: press_hyd   ! (agg)

    ! get the from the ocean the salinity diffusion coefficients
    REAL(wp), POINTER, DIMENSION(:,:,:)            :: hor_diffusion_coeff ! this is actually constant, needs to be initialized, not communicated
    REAL(wp), POINTER, DIMENSION(:,:,:) :: ver_diffusion_coeff
    
    ! thiese are actually from the atmosphere, but for the moment we will keep them here
    REAL(wp), POINTER, DIMENSION(:,:) :: short_wave_flux
    REAL(wp), POINTER, DIMENSION(:,:) :: wind10m
    REAL(wp), POINTER, DIMENSION(:,:) :: co2_mixing_ratio   

    ! thses are for the zstar
!     REAL(wp), POINTER, DIMENSION(:,:) :: eta_c  not used
    REAL(wp), POINTER, DIMENSION(:,:) :: stretch_c
    REAL(wp), POINTER, DIMENSION(:,:) :: stretch_c_new
    REAL(wp), POINTER, DIMENSION(:,:) :: draftave

  END TYPE t_ocean_to_hamocc_state
  !-------------------------_state---------------------
  
  !----------------------------------------------
  TYPE t_hamocc_to_ocean_state
  
   ! this is actually to the atmosphere, but for the moment we keep it here
   REAL(wp), POINTER, DIMENSION(:,:) :: co2_flux

   REAL(wp),  POINTER, DIMENSION(:,:,:) ::  swr_fraction 
    
  END TYPE t_hamocc_to_ocean_state
  !----------------------------------------------
  
    !----------------------------------------------
  TYPE t_hamocc_ocean_state
  
    TYPE(t_patch_3D), POINTER     :: patch_3D
    TYPE(t_ocean_to_hamocc_state) :: ocean_to_hamocc_state
    TYPE(t_hamocc_to_ocean_state) :: hamocc_to_ocean_state
    TYPE(t_ocean_transport_state), POINTER :: ocean_transport_state
  
  END TYPE t_hamocc_ocean_state
  !----------------------------------------------
  
  TYPE(t_hamocc_ocean_state), TARGET :: hamocc_ocean_state
  TYPE(t_var_list_ptr)           :: hamocc_ocean_state_list
  TYPE(t_ocean_transport_state), TARGET :: ocean_transport_state
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE construct_hamocc_ocean_state(patch_3d)
    TYPE(t_patch_3D), POINTER, INTENT(in) :: patch_3D

    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: alloc_cell_blocks, nblks_e !, nblks_v
    REAL(wp), ALLOCATABLE :: hor_diffusion_coeff_2D(:,:)
    INTEGER :: datatype_flt
    
    
    patch_2d => patch_3d%p_patch_2d(1)
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    hamocc_ocean_state%patch_3D => patch_3D
    ocean_transport_state%patch_3D => patch_3D
    hamocc_ocean_state%ocean_transport_state => ocean_transport_state
    
    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF
 
    CALL vlr_add(hamocc_ocean_state_list, 'hamocc_ocean_state_list', &
      & patch_id=patch_2d%id, lrestart=.FALSE., loutput=.TRUE.,           &
      & model_type='hamocc')

    ! transport state, ocean to hamocc
    CALL add_var(hamocc_ocean_state_list,'vn', ocean_transport_state%vn, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn', '', 'vn', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("hamocc_ocean_state"))
    ocean_transport_state%vn = 0.0_wp

    CALL add_var(hamocc_ocean_state_list,'mass_flux_e', ocean_transport_state%mass_flux_e, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('mass flux e', '', 'mass flux at edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("hamocc_ocean_state"))
    ocean_transport_state%vn = 0.0_wp
 
    CALL add_var(hamocc_ocean_state_list, 'w', ocean_transport_state%w, &
      & grid_unstructured_cell, za_depth_below_sea_half, &
      & t_cf_var('w','m/s','vertical velocity at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("hamocc_ocean_state"))
    ocean_transport_state%w = 0.0_wp
 
    CALL add_var(hamocc_ocean_state_list, 'h_old', ocean_transport_state%h_old , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('h_old', 'm', 'h_old', datatype_flt,'h_old'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    ocean_transport_state%h_old = 0.0_wp

    CALL add_var(hamocc_ocean_state_list, 'h_new', ocean_transport_state%h_new , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('h_new', 'm', 'h_new', datatype_flt,'h_new'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    ocean_transport_state%h_new = 0.0_wp


    ! zstar variables
    
    ! not used
!     CALL add_var(hamocc_ocean_state_list, 'eta_c', hamocc_ocean_state%ocean_to_hamocc_state%eta_c , &
!       & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
!       & t_cf_var('eta_c', 'm', 'eta_c', datatype_flt,'eta_c'),&
!       & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
!       & ldims=(/nproma,alloc_cell_blocks/),&
!       & in_group=groups("hamocc_ocean_state"))
!     hamocc_ocean_state%ocean_to_hamocc_state%eta_c = 0.0_wp

    CALL add_var(hamocc_ocean_state_list, 'stretch_c', hamocc_ocean_state%ocean_to_hamocc_state%stretch_c , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('stretch_c', 'unitless', 'stretch_c', datatype_flt,'stretch_c'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%stretch_c = 0.0_wp

    CALL add_var(hamocc_ocean_state_list, 'stretch_c_new', hamocc_ocean_state%ocean_to_hamocc_state%stretch_c_new , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('stretch_c_new', 'unitless', 'stretch_c_new', datatype_flt,'stretch_c_new'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%stretch_c_new = 0.0_wp

    CALL add_var(hamocc_ocean_state_list, 'draftave', hamocc_ocean_state%ocean_to_hamocc_state%draftave , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('draftave', 'm', 'draftave', datatype_flt,'draftave'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%draftave = 0.0_wp

    CALL add_var(hamocc_ocean_state_list, 'h_old_withIce', hamocc_ocean_state%ocean_to_hamocc_state%h_old_withIce , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('h_old_withIce', 'm', 'h_old_withIce', datatype_flt,'h_old_withIce'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%h_old_withIce = 0.0_wp
    
    ! ocean to hamocc
    ! just add pointers for the h_old, h_new to the transport
    hamocc_ocean_state%ocean_to_hamocc_state%h_old => hamocc_ocean_state%ocean_transport_state%h_old
    hamocc_ocean_state%ocean_to_hamocc_state%h_new => hamocc_ocean_state%ocean_transport_state%h_new
    
    CALL add_var(hamocc_ocean_state_list, 'top_dilution_coeff', hamocc_ocean_state%ocean_to_hamocc_state%top_dilution_coeff , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('top_dilution_coeff', 'unitless', 'top_dilution_coeff', datatype_flt,'top_dilution_coeff'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%top_dilution_coeff = 1.0_wp
     
    CALL add_var(hamocc_ocean_state_list, 'ice_concentration_sum', &
      & hamocc_ocean_state%ocean_to_hamocc_state%ice_concentration_sum , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('ice_concentration_sum', 'm', 'ice_concentration_sum', datatype_flt,'ice_concentration_sum'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%ice_concentration_sum = 0.0_wp
     
     CALL add_var(hamocc_ocean_state_list, 'temperature', hamocc_ocean_state%ocean_to_hamocc_state%temperature, &
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('temperature', 'C', 'temperature', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%temperature = 10.0_wp
      
     CALL add_var(hamocc_ocean_state_list, 'salinity', hamocc_ocean_state%ocean_to_hamocc_state%salinity, &
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('salinity', '', 'salinity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%salinity = 10.0_wp
      
    !! tracer horizontal diffusion
    CALL add_var(hamocc_ocean_state_list,'hor_diffusion_coeff',hamocc_ocean_state%ocean_to_hamocc_state%hor_diffusion_coeff, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn', '', 'tracer hor_diffusion_coeff', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("hamocc_ocean_state"))

   !  (agg)
     CALL add_var(hamocc_ocean_state_list, 'press_hyd', hamocc_ocean_state%ocean_to_hamocc_state%press_hyd, &
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('press_hyd', '', 'press_hyd', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%press_hyd = 1.0_wp    ! hydrostatic pressure, unit [m] (see ocean/dynamics/mo_ocean_types.f90)

 
    ! this is initialized as in the ocean, 
    ! in general, it does not need to be communicated, 
    ! unless we employ some tracer turbulance diffusion scheme 
    ALLOCATE(hor_diffusion_coeff_2D(nproma,nblks_e))    
    CALL scale_horizontal_diffusion(patch_3D=patch_3D, &
      & DiffusionScaling=TracerHorizontalDiffusion_scaling, &
      & DiffusionReferenceValue=Salinity_HorizontalDiffusion_Reference, &
      & DiffusionBackgroundValue=Salinity_HorizontalDiffusion_Background, &
      & out_DiffusionCoefficients=hor_diffusion_coeff_2D)
    CALL copy2Dto3D(hor_diffusion_coeff_2D, hamocc_ocean_state%ocean_to_hamocc_state%hor_diffusion_coeff, patch_2d%edges%all)
    DEALLOCATE(hor_diffusion_coeff_2D)    
     
     
    CALL add_var(hamocc_ocean_state_list, 'ver_diffusion_coeff', hamocc_ocean_state%ocean_to_hamocc_state%ver_diffusion_coeff, &
      & grid_unstructured_cell, za_depth_below_sea_half, &
      & t_cf_var('ver_diffusion_coeff','','ver_diffusion_coeff', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%ver_diffusion_coeff = 0.0_wp
   
    CALL add_var(hamocc_ocean_state_list, 'short_wave_flux', hamocc_ocean_state%ocean_to_hamocc_state%short_wave_flux , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('short_wave_flux', '', 'short_wave_flux', datatype_flt,'short_wave_flux'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%short_wave_flux = 0.0_wp

     CALL add_var(hamocc_ocean_state_list, 'wind10m', hamocc_ocean_state%ocean_to_hamocc_state%wind10m , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('wind10m', 'm/s', 'wind10m', datatype_flt,'wind10m'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%wind10m = 0.0_wp

     CALL add_var(hamocc_ocean_state_list, 'co2_mixing_ratio', hamocc_ocean_state%ocean_to_hamocc_state%co2_mixing_ratio , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('co2_mixing_ratio', '', 'co2_mixing_ratio', datatype_flt,'co2_mixing_ratio'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%ocean_to_hamocc_state%co2_mixing_ratio = 0.0_wp

    ! hamocc to ocean
    CALL add_var(hamocc_ocean_state_list, 'co2_flux', hamocc_ocean_state%hamocc_to_ocean_state%co2_flux , &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
      & t_cf_var('co2_flux', '', 'co2_flux', datatype_flt,'co2_flux'),&
      & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups("hamocc_ocean_state"))
    hamocc_ocean_state%hamocc_to_ocean_state%co2_flux = 0.0_wp
    
   ! relative swr absorption factor from hamocc (LFB_BGC_OCE)
    CALL add_var(hamocc_ocean_state_list, 'swr_fraction', hamocc_ocean_state%hamocc_to_ocean_state%swr_fraction,&
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('swr_fraction','1','swr_fraction', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("hamocc_ocean_state"))
   hamocc_ocean_state%hamocc_to_ocean_state%swr_fraction = 1.0_wp 
    
   END SUBROUTINE construct_hamocc_ocean_state
   !-------------------------------------------------------------------------

   !-------------------------------------------------------------------------
   SUBROUTINE destruct_hamocc_ocean_state()
   
      CALL vlr_del(hamocc_ocean_state_list)
   
   END SUBROUTINE destruct_hamocc_ocean_state
   !-------------------------------------------------------------------------
   
    
    
END MODULE mo_ocean_hamocc_couple_state

