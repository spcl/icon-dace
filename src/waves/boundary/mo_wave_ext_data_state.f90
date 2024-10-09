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

! Allocation/deallocation of external parameter state
!
! This module contains routines for setting up the external data state.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_ext_data_state

  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_ext_data_types, ONLY: t_external_wave
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, VNAME_LEN, SUCCESS
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                               GRID_CELL, GRID_EDGE
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_list,            ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_grid_config,         ONLY: n_dom
  USE mo_parallel_config,     ONLY: nproma
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, &
    &                               GRID_UNSTRUCTURED
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_zaxis_type,          ONLY: ZA_SURFACE


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_ext_data_state'

  ! variables
  PUBLIC :: wave_ext_data
  PUBLIC :: wave_ext_data_list

  ! subroutines
  PUBLIC :: construct_wave_ext_data_state
  PUBLIC :: destruct_wave_ext_data_state

  TYPE(t_external_wave), ALLOCATABLE :: wave_ext_data(:)
  TYPE(t_var_list_ptr),  ALLOCATABLE :: wave_ext_data_list(:)

CONTAINS

  !>
  !! Constructor for wave external data state and list.
  !!
  SUBROUTINE construct_wave_ext_data_state (p_patch)
    TYPE(t_patch),                      INTENT(IN)    :: p_patch(:)

    INTEGER :: jg
    INTEGER :: ist
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
         routine = modname//':construct_wave_ext_data'

    !-------------------------------------------------------------------------
    ALLOCATE(wave_ext_data(n_dom), wave_ext_data_list(n_dom), stat=ist)
    IF (ist /= SUCCESS) CALL finish(routine, &
      &  'allocation of wave ext_data state array and list failed')

    ! Build external data list for constant-in-time fields for the wave model
    DO jg = 1, n_dom
       WRITE(listname,'(a,i2.2)') 'ext_data_wave_D',jg
       CALL new_ext_data_wave_list(p_patch(jg), wave_ext_data(jg),       &
         &                          wave_ext_data_list(jg), TRIM(listname))
    END DO

    CALL message (routine, 'Construction of wave ext_data state finished')

  END SUBROUTINE construct_wave_ext_data_state


  !>
  !! Destructor for wave external data state and list.
  !!
  SUBROUTINE destruct_wave_ext_data_state

    INTEGER :: jg
    INTEGER :: ist

    CHARACTER(len=*), PARAMETER :: &
         routine = modname//':destruct_wave_ext_data'

    DO jg = 1,n_dom
      ! Delete list of wave elements
      CALL vlr_del(wave_ext_data_list(jg))
    END DO

    DEALLOCATE(wave_ext_data, wave_ext_data_list, stat=ist)
    IF (ist/=SUCCESS) CALL finish (routine,&
      & 'deallocation of wave ext_data state array and list failed')

    CALL message (TRIM(routine), 'Destruction of wave ext_data state finished')

  END SUBROUTINE destruct_wave_ext_data_state



  !>
  !! Allocation of components of wave ext_data state
  !!
  SUBROUTINE new_ext_data_wave_list( p_patch, ext_data_wave, ext_data_wave_list, listname)

    TYPE(t_patch), INTENT(IN)            :: & !< current patch
      &  p_patch

    TYPE(t_external_wave), INTENT(INOUT) :: & !< current external data structure
      &  ext_data_wave

    TYPE(t_var_list_ptr), INTENT(INOUT)  :: & !< current external data list
      &  ext_data_wave_list

    CHARACTER(len=*), INTENT(IN)         :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, & !< number of cell blocks to allocate
      &        nblks_e    !< number of edge blocks to allocate

    INTEGER :: shape2d_c(2), shape2d_e(2), shape3d_c_2(3)
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt
    INTEGER :: ist           ! status
    INTEGER :: comp

    CHARACTER(len=VNAME_LEN) :: out_name

    CHARACTER(len=*), PARAMETER :: &
         routine = modname//':new_ext_data_wave_list'

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !--------------------------------------------------------------

    ! determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ibits = DATATYPE_PACK16 ! packing accuracy of horizontal slice

    ! predefined array shapes
    shape2d_c   = (/ nproma, nblks_c /)
    shape2d_e   = (/ nproma, nblks_e /)
    shape3d_c_2 = (/2, nproma, nblks_c /)

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(                             &
      &  ext_data_wave%bathymetry_c,     &
      &  ext_data_wave%bathymetry_e,     &
      &  ext_data_wave%geo_depth_grad_c, &
      &  ext_data_wave%depth_c,          &
      &  ext_data_wave%depth_e)


    CALL vlr_add(ext_data_wave_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &          model_type=get_my_process_name())


    ! bathymetric height at cell center
    cf_desc    = t_cf_var('Model bathymetry at cell center', 'm', &
      &                   'Model bathymetry', datatype_flt)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( ext_data_wave_list, 'bathymetry_c', ext_data_wave%bathymetry_c,  &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      &           lrestart=.FALSE., loutput=.TRUE., ldims=shape2d_c )


    ! bathymetric height at cell edge
    !
    ! bathymetry_e  ext_data_wave%bathymetry_e(nproma,nblks_e)
    cf_desc    = t_cf_var('Model bathymetry at cell edges', 'm', &
                          'Model bathymetry', datatype_flt)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( ext_data_wave_list, 'bathymetry_e', ext_data_wave%bathymetry_e,  &
      &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,         &
      &           lrestart=.FALSE., loutput=.TRUE., ldims=shape2d_e )


    ! bathymetry geographical gradient
    !
    ! geo_depth_grad_c  ext_data_wave%geo_depth_grad_c(2,nproma,nblks_c)
    cf_desc    = t_cf_var('geo_depth_grad_c', 'm/m', 'depth geo-gradient at cell center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(ext_data_wave_list, 'geo_depth_grad_c', ext_data_wave%geo_depth_grad_c, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      &          lrestart=.FALSE., loutput=.FALSE.,                                     &
      &          ldims=shape3d_c_2,lcontainer=.TRUE.)

    ALLOCATE(ext_data_wave%grad_ptr(2),  STAT=ist)
    IF (ist/=SUCCESS) CALL finish(routine,               &
          'allocation of grad_ptr failed')

    comp=1 ! U-component
    out_name = 'u_depth_grad_c'
    CALL add_ref(ext_data_wave_list, 'geo_depth_grad_c',                                   &
           & TRIM(out_name), ext_data_wave%grad_ptr(comp)%p_2d,                            &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                           &
           & t_cf_var(TRIM(out_name), 'm/m', 'component of depth gradient', datatype_flt), &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),                &
           & ref_idx=comp, ldims=shape2d_c, lrestart=.FALSE., loutput=.TRUE., opt_var_ref_pos=1)

    comp=2 ! V-component
    out_name = 'v_depth_grad_c'
    CALL add_ref(ext_data_wave_list, 'geo_depth_grad_c', &
           & TRIM(out_name), ext_data_wave%grad_ptr(comp)%p_2d, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                            &
           & t_cf_var(TRIM(out_name), 'm/m','component of depth gradient', datatype_flt),        &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),                 &
           & ref_idx=comp, ldims=shape2d_c, lrestart=.FALSE., loutput=.TRUE., opt_var_ref_pos=1)

    ! water depth
    !
    ! depth_c  ext_data_wave%depth_c(nproma,nblks_c)
    cf_desc    = t_cf_var('depth_c', 'm', 'Water depth at cell center', datatype_flt)
    grib2_desc = grib2_var(10, 4, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(ext_data_wave_list, 'depth_c', ext_data_wave%depth_c,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE., ldims=shape2d_c)

    ! depth_e  ext_data_wave%depth_e(nproma,nblks_e)
    cf_desc    = t_cf_var('depth_e', 'm', 'Water depth at cell edges', datatype_flt)
    grib2_desc = grib2_var(10, 4, 14, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(ext_data_wave_list, 'depth_e', ext_data_wave%depth_e,  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE., ldims=shape2d_e)

  END SUBROUTINE new_ext_data_wave_list

END MODULE mo_wave_ext_data_state
