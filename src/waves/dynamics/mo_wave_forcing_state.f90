! Allocation/deallocation of wave forcing state
!
! This module contains routines for setting up the external data state.
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

MODULE mo_wave_forcing_state
  USE mo_exception,            ONLY: message, finish
  USE mo_var_list,             ONLY: add_var, t_var_list_ptr
  USE mo_wave_forcing_types,   ONLY: t_wave_forcing
  USE mo_grid_config,          ONLY: n_dom
  USE mo_impl_constants,       ONLY: max_char_length, SUCCESS
  USE mo_model_domain,         ONLY: t_patch
  USE mo_zaxis_type,           ONLY: ZA_SURFACE, ZA_HEIGHT_10M
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_parallel_config,      ONLY: nproma
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, &
       &                             GRID_UNSTRUCTURED
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_var_groups,           ONLY: groups
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL, &
       &                             GRID_UNSTRUCTURED_EDGE, GRID_EDGE


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_forcing_state'

  PUBLIC :: construct_wave_forcing_state
  PUBLIC :: destruct_wave_forcing_state
  PUBLIC :: wave_forcing_state
  PUBLIC :: wave_forcing_state_list

  TYPE(t_wave_forcing),       TARGET, ALLOCATABLE :: wave_forcing_state(:)
  TYPE(t_var_list_ptr),       TARGET, ALLOCATABLE :: wave_forcing_state_list(:)

CONTAINS

  SUBROUTINE construct_wave_forcing_state(p_patch)

    TYPE(t_patch),                      INTENT(IN)    :: p_patch(:)

    CHARACTER(len=max_char_length) :: listname
    CHARACTER(len=*), PARAMETER :: routine = modname//'::construct_wave_forcing_state'

    INTEGER :: jg, ist

    !-------------------------------------------------------------------------
    ALLOCATE(wave_forcing_state(n_dom), wave_forcing_state_list(n_dom), stat=ist)
    IF (ist /= SUCCESS) CALL finish(routine, &
         &  'allocation of wave forcing state array and list failed')

    DO jg = 1, n_dom
      ! Build forcing state list
      ! includes memory allocation
      WRITE(listname,'(a,i2.2)') 'wave_forcing_state_of_domain_',jg
      CALL new_wave_forcing_state_list(p_patch(jg), wave_forcing_state(jg), wave_forcing_state_list(jg), TRIM(listname))
    END DO

  END SUBROUTINE construct_wave_forcing_state


  SUBROUTINE new_wave_forcing_state_list(p_patch, p_forcing, p_forcing_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_forcing),  INTENT(INOUT) :: p_forcing
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_forcing_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: nblks_c, nblks_e
    INTEGER :: shape2d_c(2), shape2d_e(2)

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    shape2d_c      = (/nproma, nblks_c/)
    shape2d_e      = (/nproma, nblks_e/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(p_forcing_list, listname, patch_id=p_patch%id, lrestart=.TRUE.)

    !forcing group
    !wind 10
    cf_desc    = t_cf_var('u10m', 'm s-1 ','U-Component of wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'u10m', p_forcing%u10m,                         &
         &        GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('v10m', 'm s-1 ','V-Component of wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'v10m', p_forcing%v10m,                         &
         &        GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('sp10m', 'm s-1 ','Wind speed (SP_10M)', datatype_flt)
    grib2_desc = grib2_var(0, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'sp10m', p_forcing%sp10m,                      &
         &        GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('dir10m', 'deg ','Wind direction (DD_10M)', datatype_flt)
    grib2_desc = grib2_var(0, 2, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'dir10m', p_forcing%dir10m,                       &
         &        GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    !sea ice
    cf_desc    = t_cf_var('sea_ice_c', 'frac','sea ice fraction at cells', datatype_flt)
    grib2_desc = grib2_var(10, 2, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'sea_ice_c', p_forcing%sea_ice_c,              &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('sea_ice_e', 'frac','sea ice fraction at edges', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_forcing_list, 'sea_ice_e', p_forcing%sea_ice_e,              &
         &        GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_e, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('ice_free_mask_c', '-', 'ice-free mask at cells', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'ice_free_mask_c', p_forcing%ice_free_mask_c,  &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_forcing") )

    !sea level
    cf_desc    = t_cf_var('sea_level_c', 'm','sea level height at cells', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'sea_level_c', p_forcing%sea_level_c,          &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('sea_level_e', 'm','sea level height at edges', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_forcing_list, 'sea_level_e', p_forcing%sea_level_e,          &
         &        GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_e, in_group=groups("wave_forcing") )

    !ocean surface currents
    cf_desc    = t_cf_var('usoce_c', 'm s-1 ','zonal ocean surface current at cells', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'usoce_c', p_forcing%usoce_c,                  &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('usoce_e', 'm s-1 ','zonal ocean surface current at edges', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_forcing_list, 'usoce_e', p_forcing%usoce_e,                  &
         &        GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_e, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('vsoce_c', 'm s-1 ','meridional ocean surface current at cells', datatype_flt)
    grib2_desc = grib2_var(10, 1, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'vsoce_c', p_forcing%vsoce_c,                  &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('vsoce_e', 'm s-1 ','meridional ocean surface current at edges', datatype_flt)
    grib2_desc = grib2_var(10, 1, 3, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_forcing_list, 'vsoce_e', p_forcing%vsoce_e,                  &
         &        GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_e, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('sp_soce_c', 'm s-1 ','ocean surface current velocity at cells', datatype_flt)
    grib2_desc = grib2_var(10, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'sp_soce_c', p_forcing%sp_soce_c,              &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

    cf_desc    = t_cf_var('dir_soce_c', 'deg','ocean surface current direction at cells', datatype_flt)
    grib2_desc = grib2_var(10, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_forcing_list, 'dir_soce_c', p_forcing%dir_soce_c,            &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         &        ldims=shape2d_c, in_group=groups("wave_forcing") )

  END SUBROUTINE new_wave_forcing_state_list

  !>
  !! Destruction of wave forcing variable lists and memory deallocation
  !!
  SUBROUTINE destruct_wave_forcing_state()

    INTEGER :: jg, ist

    CHARACTER(len=*), PARAMETER :: routine = modname//'::destruct_wave_forcing_state'

    DO jg = 1, n_dom
      ! delete wave forcing  state list elements
      CALL vlr_del(wave_forcing_state_list(jg))
    END DO

    !-------------------------------------------------------------------------
    DEALLOCATE(wave_forcing_state, wave_forcing_state_list, stat=ist)
    IF (ist /= SUCCESS) CALL finish(routine, &
         &  'deallocation of wave forcing state array and list failed')

    CALL message (TRIM(routine), 'Destruction of wave forcing state finished')

  END SUBROUTINE destruct_wave_forcing_state

END MODULE mo_wave_forcing_state
