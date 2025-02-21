!
! Classes and functions for the turbulent mixing package (tmx)
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

MODULE mo_tmx_field_class

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_variable, ONLY: t_variable, allocate_variable
  USE mo_variable_list, ONLY: t_variable_list, t_variable_item, variable_list_id

  ! Todo: refactor so that t_patch is not needed
  USE mo_model_domain      ,ONLY: t_patch








  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_tmx_field, t_tmx_field_list, &! bind_tmx_field, &
    &       t_domain, isfc_oce, isfc_ice, isfc_lnd

  TYPE, EXTENDS(t_variable) :: t_tmx_field
    INTEGER :: id = -1
    INTEGER :: type = -1 ! Field type
    ! PROCEDURE(i_convert), POINTER :: convert_field => NULL()
  END TYPE t_tmx_field

  ! ABSTRACT INTERFACE
  !   SUBROUTINE i_convert(this, inputs)
  !     IMPORT t_tmx_field, t_variable_list
  !     CLASS(t_tmx_field), INTENT(in) :: this
  !     CLASS(t_variable_list), INTENT(in) :: inputs
  !   END SUBROUTINE
  ! END INTERFACE
  INTERFACE t_tmx_field
    MODULE PROCEDURE t_tmx_field_constructor
  END INTERFACE

  TYPE, EXTENDS(t_variable_list) :: t_tmx_field_list
  CONTAINS
    PROCEDURE :: search_field                  => t_tmx_field_list_search
  END TYPE t_tmx_field_list

  INTERFACE t_tmx_field_list
    MODULE PROCEDURE t_tmx_field_list_constructor
  END INTERFACE

  ! INTERFACE bind_tmx_field
  !   MODULE PROCEDURE bind_tmx_field_r2d
  !   MODULE PROCEDURE bind_tmx_field_r3d
  ! END INTERFACE

  TYPE t_domain
    INTEGER ::              &
      & nlev = 0,           &
      & ntiles = 0,         &
      & nproma = 1,         &
      & npromz = 1,         &
      & nblks_c = 1,        & ! Number of blocks (cells)
      & i_startblk_c = 1,   & ! Start block on cells
      & i_endblk_c = 1,     & ! End block on cells
      & nblks_e = 1,        & ! Number of blocks (edges)
      & i_startblk_e = 1,   & ! Start block on edges
      & i_endblk_e = 1,     & ! End block on edges
      & nblks_v = 1           ! Number of blocks (vertice)
    INTEGER, ALLOCATABLE :: &
      & sfc_types(:),       & ! Surface type for each tile if ntiles > 0
      & i_startidx_c(:),    & ! Start indices on cells (for each block)
      & i_endidx_c(:),      & ! End indices on cells (for each block)
      & i_startidx_e(:),    & ! Start indices on edges (for each block)
      & i_endidx_e(:)         ! End indices on edges (for each block)
    REAL(wp), ALLOCATABLE :: &
      & lon(:,:), lat(:,:)
    TYPE(t_patch), POINTER :: patch
  END TYPE

  INTERFACE t_domain
    PROCEDURE t_domain_constructor
  END INTERFACE

  ! Surface types
  ! Todo: currently, 1-3 need to be consistent with surface types in aes!
  ENUM, BIND(C)
    ENUMERATOR :: isfc_oce=1, isfc_ice, isfc_lnd
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_field_class'
  
CONTAINS

  FUNCTION t_tmx_field_constructor(name, dims, type) RESULT(field)

    CHARACTER(len=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: dims(:)
    INTEGER,          INTENT(in) :: type
    TYPE(t_tmx_field)            :: field

    TYPE(t_variable) :: var

    var = t_variable(name, dims, "", type_id="real" )
    field%name = var%name
    field%units = var%units
    field%type_id = var%type_id
    field%dim = var%dim
    field%dims = var%dims
    field%l_opt = var%l_opt

    field%type = type

  END FUNCTION t_tmx_field_constructor

  FUNCTION t_tmx_field_list_constructor(name) result(this_variable_list)

    TYPE(t_tmx_field_list) :: this_variable_list

    CHARACTER(len=*), INTENT(IN) :: name

    variable_list_id                      = variable_list_id + 1
    this_variable_list%variable_list_id   = variable_list_id
    this_variable_list%variable_list_name = name

    ALLOCATE(this_variable_list%variable_list)

  END FUNCTION t_tmx_field_list_constructor

  FUNCTION t_tmx_field_list_search(this, name) result(tv)
    CLASS (t_tmx_field_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    CLASS(t_tmx_field), POINTER :: tv
    
    TYPE(t_variable_item), POINTER :: item
    CLASS(*), POINTER :: variable

    tv => NULL()

    item => this%getFirstVariable()
    DO WHILE ( (.NOT. item%is_item_equal_to_key(name)) .AND. ASSOCIATED(item) )
      item => this%getNextVariable(item) 
    ENDDO

    variable => item%item_value

    SELECT TYPE (variable)
      CLASS is (t_tmx_field)
        tv => variable
      CLASS DEFAULT
        tv => NULL()
    END SELECT

  END FUNCTION t_tmx_field_list_search

  ! SUBROUTINE bind_tmx_field_r2d( tv, v )
  !   CLASS(t_tmx_field), POINTER :: tv
  !   REAL(wp), POINTER :: v(:,:)
  !   ! REAL(wp), TARGET, INTENT(in) :: v(:,:)
  !   IF (tv%bound) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
  !   ELSE
  !     tv%r2d => v
  !     tv%bound = .true.
  !   ENDIF
  !   IF ( tv%dim /= 2 ) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " is not two dimensional"
  !   ENDIF
  !   ! IF (.not. ASSOCIATED(v) ) THEN
  !   !     PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
  !   !   STOP
  !   ! ENDIF
  !   IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
  !   ENDIF
  ! END SUBROUTINE bind_tmx_field_r2d

  ! SUBROUTINE bind_tmx_field_r3d( tv, v )
  !   CLASS(t_tmx_field), POINTER :: tv
  !   REAL(wp), POINTER :: v(:,:,:)
  !   ! REAL(wp), TARGET, INTENT(in) :: v(:,:,:)
  !   IF (tv%bound) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
  !   ELSE
  !     tv%r3d => v
  !     tv%bound = .true.
  !   ENDIF
  !   IF ( tv%dim /= 3 ) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " is not three dimensional"
  !   ENDIF
  !   ! IF (.not. ASSOCIATED(v) ) THEN
  !   !     PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
  !   !   STOP
  !   ! ENDIF
  !   IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
  !   ENDIF
  ! END SUBROUTINE bind_tmx_field_r3d

  FUNCTION t_domain_constructor(patch, nproma, nlev, ntiles, sfc_types) RESULT(domain)

    USE mo_loopindices,        ONLY: get_indices_e, get_indices_c
    USE mo_impl_constants,     ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int
    USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
    USE mo_math_constants,     ONLY: rad2deg

    TYPE(t_domain), POINTER :: domain

    TYPE(t_patch), POINTER :: patch
    INTEGER, INTENT(in) :: &
      & nproma
    INTEGER, INTENT(in), OPTIONAL :: &
      & nlev, &
      & ntiles, &
      & sfc_types(:)

    INTEGER :: rl_start, rl_end, jb

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_domain_constructor'

    ALLOCATE(domain)

    domain%nproma = nproma

    domain%patch => patch
    

    domain%npromz = patch%npromz_c
    IF (PRESENT(nlev)) THEN
      domain%nlev = nlev
    END IF

    ALLOCATE(domain%lon(nproma,patch%nblks_c))
    ALLOCATE(domain%lat(nproma,patch%nblks_c))
    domain%lon(:,:) = rad2deg * patch%cells%center(:,:)%lon
    domain%lat(:,:) = rad2deg * patch%cells%center(:,:)%lat

    IF (PRESENT(ntiles)) THEN
      domain%ntiles = ntiles
    END IF
    IF (domain%ntiles > 0) THEN
      IF (.NOT. PRESENT(sfc_types)) CALL finish(routine, 'sfc_types required')
      ALLOCATE(domain%sfc_types, source=sfc_types)
    END IF

    domain%nblks_c = patch%nblks_c
    domain%nblks_e = patch%nblks_e
    domain%nblks_v = patch%nblks_v

    ALLOCATE(domain%i_startidx_c(patch%nblks_c))
    ALLOCATE(domain%i_endidx_c  (patch%nblks_c))
    ALLOCATE(domain%i_startidx_e(patch%nblks_e))
    ALLOCATE(domain%i_endidx_e  (patch%nblks_e))
    
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int
    domain%i_startblk_c = patch%cells%start_block(rl_start)
    domain%i_endblk_c   = patch%cells%end_block(rl_end)
    DO jb=domain%i_startblk_c,domain%i_endblk_c
      CALL get_indices_c(patch, jb, domain%i_startblk_c, domain%i_endblk_c, &
        &                           domain%i_startidx_c(jb), domain%i_endidx_c(jb), rl_start, rl_end)
    END DO

    rl_start = grf_bdywidth_e + 1
    rl_end   = min_rledge_int
    domain%i_startblk_e = patch%edges%start_block(rl_start)
    domain%i_endblk_e   = patch%edges%end_block(rl_end)
    DO jb=domain%i_startblk_e,domain%i_endblk_e
      CALL get_indices_e(patch, jb, domain%i_startblk_e, domain%i_endblk_e, &
        &                           domain%i_startidx_e(jb), domain%i_endidx_e(jb), rl_start, rl_end)
    END DO

    !$ACC ENTER DATA COPYIN(domain)
    !$ACC ENTER DATA COPYIN(domain%lon, domain%lat)
    !$ACC ENTER DATA COPYIN(domain%i_startidx_c, domain%i_endidx_c, domain%i_startidx_e, domain%i_endidx_e)
    !$ACC ENTER DATA COPYIN(domain%sfc_types)

  END FUNCTION t_domain_constructor

END MODULE mo_tmx_field_class
