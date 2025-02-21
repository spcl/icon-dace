!  Contains the data structures for the hydrostatic ocean model.
!
!  Contains the data structures to store the hydrostatic & boussinesq ocean model state.
!  Implementation is based on ICON-Shallow-Water model
!  to store the shallow water model state and other auxiliary variables.
!  Constructors and destructors for these data structures are also defined here.
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
MODULE mo_icon_output_variables
  !-------------------------------------------------------------------------
  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_kind,                ONLY: wp, sp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: success, max_char_length, TLEV_NNEW
  USE mo_run_config,          ONLY: test_mode
  USE mo_mpi,                 ONLY: get_my_global_mpi_id, global_mpi_barrier,my_process_is_mpi_test
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,     &
    &                               success, max_char_length, MIN_DOLIC,                  &
    &                               full_coriolis, beta_plane_coriolis,                   &
    &                               f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell,         &
    &                               grid_unstructured_edge, grid_unstructured_vert,       &
    &                               grid_vertex
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch,t_patch_3d, t_grid_cells, t_grid_edges, p_patch_local_parent
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, grid_sphere_radius, grid_angular_velocity, &
    & use_dummy_cell_closure
  USE mo_dynamics_config,     ONLY: nnew, nold, nnow
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_var_list,            ONLY: add_var,                  &
    &                               add_ref, t_var_list_ptr
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_metadata,        ONLY: get_timelevel_string
  USE mo_var_groups,          ONLY: groups 
  USE mo_cf_convention
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_grib2,               ONLY: grib2_var, t_grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32 => CDI_DATATYPE_FLT32, &
    &                               DATATYPE_FLT64 => CDI_DATATYPE_FLT64, &
    &                               DATATYPE_INT8 => CDI_DATATYPE_INT8, &
    &                               DATATYPE_PACK16 => CDI_DATATYPE_PACK16, &
    &                               tstep_constant, GRID_LONLAT, GRID_UNSTRUCTURED
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell, grid_unstructured_edge, &
    &                               grid_unstructured_vert, grid_vertex, GRID_ZONAL
  USE mo_zaxis_type,          ONLY: za_depth_below_sea, za_depth_below_sea_half, za_surface
  !  USE mo_ocean_config,        ONLY: ignore_land_points
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_alloc_patches,       ONLY: destruct_patches

  IMPLICIT NONE
  PRIVATE

  !public interface
  !
  ! subroutines
  PUBLIC :: construct_icon_output_variables
  PUBLIC :: destruct_icon_output_variables
  PUBLIC :: zlevels, dz_full_level
  PUBLIC :: patch_3d
 
  !----------------------------------------------------------------------------
  INTEGER :: zlevels = 0
  INTEGER, PARAMETER :: max_allocated_levels = 1024
  REAL(wp) :: dz_full_level(max_allocated_levels) = 0  ! namelist input of layer thickness
  
  TYPE t_output_collection
    REAL(sp),  POINTER, DIMENSION(:,:,:) :: output_variable
  END TYPE t_output_collection

  ! variables
  TYPE(t_patch_3d), POINTER :: patch_3d => NULL()
  TYPE(t_var_list_ptr)  :: output_default_list
  TYPE(t_output_collection) :: myOutputCollection

CONTAINS

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE construct_icon_output_variables
  
    CHARACTER(LEN=max_char_length) :: listname
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(len=64) :: model_name
    INTEGER :: alloc_cell_blocks
    INTEGER :: datatype_flt
    
    ! creat a var_list
    model_name=get_my_process_name()
    patch_2d => patch_3d%p_patch_2d(1)
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    datatype_flt = DATATYPE_FLT32

    
    ! IMO the number of variable lists should be as small as possible
    ! default list: elements can be written to disk, but not to the restart file
    WRITE(listname,'(a)')  'output_default_list'
    CALL vlr_add(output_default_list, listname, patch_id=patch_2d%id, &
      & lrestart=.FALSE.,model_type=TRIM(model_name), loutput=.TRUE.)

    ! add an output variable
    CALL add_var(output_default_list, 'out_var', myOutputCollection%output_variable, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('out_var','-','output variable', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,zlevels,alloc_cell_blocks/),in_group=groups("output_model"))
    myOutputCollection%output_variable = 0.0_wp
    
      
  END SUBROUTINE construct_icon_output_variables
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE destruct_icon_output_variables
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'destruct_icon_output_variables'
 
    !-------------------------------------------------------------------------
    CALL message(TRIM(method_name), 'starting...')

    CALL vlr_del(output_default_list)

    !The 3D-icon_output version of previous calls
    CALL destruct_patches( patch_3d%p_patch_2d )
    CALL destruct_patches( p_patch_local_parent )
    NULLIFY( patch_3d%p_patch_2d )
!     CALL destruct_patch_3d( patch_3d )

  END SUBROUTINE destruct_icon_output_variables
  !-------------------------------------------------------------------------

END MODULE mo_icon_output_variables
