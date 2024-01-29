!
! Additional storage for Spectral Bin Microphysics (SBM)
!
! Stores temperature and humidity fields, which are specific to the
! Spectral Bin Microphysics scheme.
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

MODULE mo_sbm_storage

  USE mo_kind,                    ONLY: wp
  USE mo_impl_constants,          ONLY: success, max_char_length
  USE mo_exception,               ONLY: message, finish
  USE mo_model_domain,            ONLY: t_patch
  USE mo_var_list_register,       ONLY: vlr_add, vlr_del
  USE mo_var_list,                ONLY: add_var, t_var_list_ptr
  USE mo_zaxis_type,              ONLY: ZA_REFERENCE
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var, grib2_var
  USE mo_parallel_config,         ONLY: nproma
  USE mo_grid_config,             ONLY: n_dom
  USE mo_cdi,                     ONLY: DATATYPE_PACK16, GRID_UNSTRUCTURED, TSTEP_INSTANT, &
    &                                   DATATYPE_FLT32
  USE mo_cdi_constants,           ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL

  IMPLICIT NONE

  PRIVATE

  ! types
  PUBLIC :: t_sbm_storage

  ! subroutines/functions
  PUBLIC :: construct_sbm_storage
  PUBLIC :: destruct_sbm_storage
  PUBLIC :: get_sbm_storage


  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sbm_storage'


  ! SBM specific storage object
  TYPE t_sbm_storage
    REAL(wp), POINTER, CONTIGUOUS ::  &
     &  qv_before_satad   (:,:,:),    &    !< water vapour mass fraction after dynamics and transport
                                           !  right before the first satad call [kg/kg]
     &  temp_before_satad (:,:,:),    &    !< same for air temperature [K]
     &  qv_old            (:,:,:),    &    !< water vapour mass fraction at the end of a physics time step [kg/kg]
     &  temp_old          (:,:,:)          !< Temperature at the end of a physics time step [K]
  END TYPE t_sbm_storage


  ! SBM storage object
  TYPE (t_sbm_storage),  ALLOCATABLE, TARGET :: sbm_storage(:)

  ! variable list
  TYPE (t_var_list_ptr), ALLOCATABLE :: sbm_storage_list(:)

CONTAINS

  !>
  !! get pointer to domain-specific SBM storage
  !!
  FUNCTION get_sbm_storage(patch_id) RESULT(ptr_sbm_storage)
    INTEGER, INTENT(IN)          :: patch_id   ! domain ID
    TYPE(t_sbm_storage), POINTER :: ptr_sbm_storage

    ptr_sbm_storage => sbm_storage(patch_id)
  END FUNCTION get_sbm_storage


  !>
  !! Constructor for SBM storage
  !!
  SUBROUTINE construct_sbm_storage (p_patch)

    TYPE(t_patch), INTENT(in) :: p_patch(:)

    ! local variables
    INTEGER :: jg
    INTEGER :: ist                             !< error status
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(*), PARAMETER :: routine = modname//':construct_sbm_storage'


    ! Allocate pointer arrays sbm_storage, as well as the corresponding list arrays.
    !
    ALLOCATE(sbm_storage(n_dom), sbm_storage_list(n_dom),STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'allocation of sbm_storage array and list failed')
    ENDIF

    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'SBM_storage_of_domain_',jg
      CALL new_sbm_storage_list( p_patch(jg), listname, sbm_storage_list(jg), sbm_storage(jg))
    ENDDO

    CALL message(routine, 'construction of SBM storage finished')

  END SUBROUTINE construct_sbm_storage


  !>
  !! Destructor for SBM storage
  !!
  SUBROUTINE destruct_sbm_storage ()

    ! local variables
    INTEGER :: jg
    INTEGER :: ist                             !< error status
    CHARACTER(*), PARAMETER :: routine = modname//'destruct_sbm_storage'

    !--------------------------------------------------------------

    ! delete sbm_storage varlist
    DO jg = 1, n_dom
      CALL vlr_del(sbm_storage_list(jg))
    ENDDO

    DEALLOCATE(sbm_storage, sbm_storage_list, STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'deallocation of sbm_storage array and list failed')
    ENDIF

    CALL message(routine, 'destruction of sbm_storage state finished')

  END SUBROUTINE destruct_sbm_storage


  !>
  !! Constructor for SBM storage
  !!
  SUBROUTINE new_sbm_storage_list( p_patch, listname, sbm_storage_list, sbm_storage)

    TYPE(t_patch)         , INTENT(IN   ) :: p_patch
    CHARACTER(len=*)      , INTENT(IN   ) :: listname
    TYPE(t_var_list_ptr)  , INTENT(INOUT) :: sbm_storage_list
    TYPE(t_sbm_storage)   , INTENT(INOUT) :: sbm_storage

    ! local
    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt
    INTEGER :: shape3d_c(3)

!   CHARACTER(len=*), PARAMETER :: &
!     routine = modname//'new_sbm_storage_list'

    !
    ! Ensure that all pointers have a defined association status
    !
    NULLIFY( &
      &     sbm_storage%qv_before_satad,   &
      &     sbm_storage%temp_before_satad, &
      &     sbm_storage%qv_old,            &
      &     sbm_storage%temp_old           &
      &     )

    ibits        = DATATYPE_PACK16   ! "entropy" of horizontal slice
    datatype_flt = DATATYPE_FLT32

    shape3d_c = (/nproma, p_patch%nlev, p_patch%nblks_c /)

    ! Register a field list and apply default settings
    CALL vlr_add(sbm_storage_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE.)


    ! &      sbm_storage%qv_before_satad(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('qv_before_satad', 'kg kg-1', &
      &                   'qv before satad (after dynamics and transport)', datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sbm_storage_list, 'qv_before_satad', sbm_storage%qv_before_satad, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                & ldims=shape3d_c, loutput=.TRUE.,                           &
                & isteptype=TSTEP_INSTANT, lopenacc=.FALSE. )


    ! &      sbm_storage%temp_before_satad(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('temp_before_satad', 'K', &
      &                   'temperature before satad (after dynamics and transport)', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sbm_storage_list, 'temp_before_satad', sbm_storage%temp_before_satad, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                & ldims=shape3d_c, loutput=.TRUE.,                           &
                & isteptype=TSTEP_INSTANT, lopenacc=.FALSE. )


    ! &      sbm_storage%qv_old(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('qv_old', 'kg kg-1', &
      &                   'qv at the end of a physics time step', datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sbm_storage_list, 'qv_old', sbm_storage%qv_old, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                & ldims=shape3d_c, loutput=.TRUE.,                           &
                & isteptype=TSTEP_INSTANT, lopenacc=.FALSE. )


    ! &      sbm_storage%temp_old(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('temp_old', 'K', &
      &                   'temperature at the end of a physics time step', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sbm_storage_list, 'temp_old', sbm_storage%temp_old,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                & ldims=shape3d_c, loutput=.TRUE.,                           &
                & isteptype=TSTEP_INSTANT, lopenacc=.FALSE. )

  END SUBROUTINE new_sbm_storage_list

END MODULE mo_sbm_storage
