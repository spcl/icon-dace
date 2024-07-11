! Allocation/deallocation of external parameter state
!
! This module contains routines for setting up the external data state.
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
!----------------------------

MODULE mo_wave_ext_data_state

  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_ext_data_types, ONLY: t_external_wave
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                               GRID_CELL, GRID_EDGE
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_list,            ONLY: add_var, t_var_list_ptr
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
  SUBROUTINE construct_wave_ext_data_state (p_patch, wave_ext_data, wave_ext_data_list)
    TYPE(t_patch),                      INTENT(IN)    :: p_patch(:)
    TYPE(t_external_wave), ALLOCATABLE, INTENT(INOUT) :: wave_ext_data(:)
    TYPE(t_var_list_ptr),  ALLOCATABLE, INTENT(INOUT) :: wave_ext_data_list(:)

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

    INTEGER :: shape2d_c(2), shape2d_e(2)

    INTEGER :: ibits         !< "entropy" of horizontal slice

    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%alloc_cell_blocks
    nblks_e = p_patch%nblks_e

    ibits = DATATYPE_PACK16 ! packing accuracy of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)
    shape2d_e = (/ nproma, nblks_e /)

    CALL vlr_add(ext_data_wave_list, TRIM(listname), patch_id=p_patch%id, lrestart=.FALSE., &
      &          model_type=get_my_process_name())

    ! bathymetric height at cell center
    cf_desc    = t_cf_var('Model bathymetry at cell center', 'm', &
      &                   'Model bathymetry', datatype_flt)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( ext_data_wave_list, 'bathymetry_c', ext_data_wave%bathymetry_c,      &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! bathymetric height at cell edge
    !
    ! bathymetry_e  p_ext_oce%bathymetry_e(nproma,nblks_e)
    cf_desc    = t_cf_var('Model bathymetry at edge', 'm', &
                          'Model bathymetry', datatype_flt)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( ext_data_wave_list, 'bathymetry_e', ext_data_wave%bathymetry_e,      &
                  GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

  END SUBROUTINE new_ext_data_wave_list

END MODULE mo_wave_ext_data_state
