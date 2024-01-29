! Types and procedures for save/restore of 5D fields
!
! This module contains types and procedures for saving and
! restoring arbitrary 5D fields to/from a key-value storage.
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

MODULE mo_save_restore

  USE mo_kind,                    ONLY: dp, sp
  USE mo_impl_constants,          ONLY: vname_len, REAL_T, SINGLE_T, INT_T, BOOL_T, SUCCESS, &
    &                                   TLEV_NNOW, TLEV_NNOW_RCF
  USE mo_cdi_constants,           ONLY: GRID_REGULAR_LONLAT
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_hash_table,              ONLY: t_HashTable
  USE mo_model_domain,            ONLY: t_patch
  USE mo_var_list_register,       ONLY: t_vl_register_iter
  USE mo_var_groups,              ONLY: var_groups_dyn
  USE mo_var,                     ONLY: t_var
  USE mo_var_metadata,            ONLY: get_var_timelevel, get_var_name
  USE mo_var_list_register_utils, ONLY: vlr_group
  USE mo_dynamics_config,         ONLY: nnow, nnow_rcf
  USE mo_run_config,              ONLY: msg_level
  USE mo_mpi,                     ONLY: my_process_is_stdio
  USE mo_util_string,             ONLY: pretty_print_string_list
  USE mo_fortran_tools,           ONLY: copy, init

#if defined(__PGI) || defined(__FLANG)
  USE mo_util_texthash,           ONLY: t_char_workaround
#endif

  IMPLICIT NONE

  PRIVATE

  !> Debug flag for generating additional debug output
  LOGICAL :: ldebug=.FALSE.


  ! types
  PUBLIC :: t_saved_field

  ! subroutines
  PUBLIC :: save_var_group_state 
  PUBLIC :: restore_var_group_state
  PUBLIC :: reinit_var_group_state

  ! An object of this type allows to save/restore a single arbitrary 5D field.
  !
  ! This type is targeted at the data stored in ICON's variable lists.
  ! It is similar to mo_var.f90::t_var, but unlike t_var it does NOT store
  ! variable meta-data (only the raw data). Moreover it makes use of
  ! allocatable arrays rather than pointers.
  !
  TYPE :: t_saved_field
    REAL(dp), ALLOCATABLE :: r(:,:,:,:,:)
    REAL(sp), ALLOCATABLE :: s(:,:,:,:,:)
    INTEGER,  ALLOCATABLE :: i(:,:,:,:,:)
    LOGICAL,  ALLOCATABLE :: l(:,:,:,:,:)
  CONTAINS
    PROCEDURE, PRIVATE :: t_saved_field_put_r, t_saved_field_put_s, t_saved_field_put_i, t_saved_field_put_l
    GENERIC :: put => t_saved_field_put_r, t_saved_field_put_s, t_saved_field_put_i, t_saved_field_put_l

    PROCEDURE, PRIVATE :: t_saved_field_get_r, t_saved_field_get_s, t_saved_field_get_i, t_saved_field_get_l
    GENERIC :: get => t_saved_field_get_r, t_saved_field_get_s, t_saved_field_get_i, t_saved_field_get_l

    PROCEDURE :: is_allocated => t_saved_field_is_allocated
  END TYPE t_saved_field



CONTAINS

  !>
  !! Helper function to extract a value from the hash table.
  FUNCTION get_val (tbl, key) RESULT(val)
    TYPE(t_HashTable), INTENT(IN) :: tbl !< Hash table.
    CHARACTER(*), TARGET, INTENT(IN) :: key !< Key string.
    CLASS(*), POINTER :: val !< Value associated with `key` or NULL

    CLASS(*), POINTER :: ptr

    ptr => key
    val => tbl%getEntry(ptr)
  END FUNCTION



  !>
  !! Save state of a variable group. This is e.g. used for iterative IAU (see mo_iau.f90).
  !!
  !! Fields that need to be stored are identified by their group id.
  !! The save/restore mechanism is based on a hash table and operates on
  !! character strings as keys. The character strings contain the field's
  !! name string, excluding its time level. Therefore this mechanism cannot
  !! be applied to groups which contain different variables of the same name
  !! (e.g. the same field defined on different grids)
  !!
  !! Note:
  !! * Variables defined on regular lat/lon grids are ignored.
  !!
  SUBROUTINE save_var_group_state (group_name, p_patch, fields)
    CHARACTER(len=*),           INTENT(IN) :: group_name  !< Name of the group to save.
    TYPE(t_patch),              INTENT(IN) :: p_patch     !< current patch
    TYPE(t_HashTable), POINTER, INTENT(IN) :: fields      !< storage

    CHARACTER(*), PARAMETER :: routine = 'save_var_group_state'

    TYPE(t_vl_register_iter) :: iter
    INTEGER :: group_id
    CHARACTER(len=vname_len) :: savedGroup(1000)
    INTEGER :: savedGroupSize
    INTEGER :: i, jg
#if defined(__PGI) || defined(__FLANG)
    TYPE(t_char_workaround), POINTER :: key_p
#endif


    ! skip patch, if inactive
    IF (.NOT. p_patch%ldom_active) RETURN

    ! get current patch ID
    jg = p_patch%id

    IF (msg_level >= 10 .AND. my_process_is_stdio()) THEN
      CALL vlr_group(group_name, savedGroup, savedGroupSize, loutputvars_only=.FALSE., &
        &            lremap_lonlat=.FALSE., opt_lskip_container=.FALSE.)
      !
      WRITE (0,*) " "
      WRITE (0,*) "---------------------------------------"
      WRITE (0,'(a,a,a,i2)') " save_var_group_state: Variable group info for group ", TRIM(group_name), " on DOM ", jg
      WRITE (0,*) "---------------------------------------"
      !
      CALL pretty_print_string_list(savedGroup(1:savedGroupSize))
      WRITE (0,*) " "
    ENDIF

    group_id = var_groups_dyn%group_id(TRIM(group_name))

    DO WHILE (iter%next())

      IF (iter%cur%p%patch_id /= jg) CYCLE    ! skip lists which do not belong to current patch

      DO i = 1, iter%cur%p%nvars
        CALL save(fields=fields, var=iter%cur%p%vl(i)%p)
      END DO
    END DO

  CONTAINS

    SUBROUTINE save(fields, var)
      TYPE(t_HashTable), INTENT(INOUT) :: fields
      TYPE(t_var),       INTENT(IN)    :: var

      CLASS(*), POINTER :: key
      CLASS(*), POINTER :: val
      TYPE(t_saved_field), POINTER :: field
      INTEGER :: tl

      !
      ! skip variables which do not belong to the desired group
      IF (.NOT. var%info%in_group(group_id)) RETURN
      !
      ! skip variables which are defined on a regular lat-lon grid
      IF (var%info%hgrid == GRID_REGULAR_LONLAT) RETURN
      !
      ! abort if a variable is a reference to a container slice and 
      ! if the corresponding container is in the same group.
      ! This avoids duplicate entries of the same thing.
      IF (ASSOCIATED(var%ref_to)) THEN  ! variable at hand is a reference
        IF (var%ref_to%info%in_group(group_id)) THEN
          WRITE (message_text,'(a,a,a,a,a)') 'container reference ', TRIM(get_var_name(var%info)), &
            &  ' and the container ', TRIM(get_var_name(var%ref_to%info)), ' must not be in the same group.'
          CALL finish(routine, message_text)
        ENDIF
      ENDIF


      tl = get_var_timelevel(var%info%name)

      IF (tl > 0) THEN
        IF (var%info%tlev_source == TLEV_NNOW .AND. tl /= nnow(jg)) RETURN
        IF (var%info%tlev_source == TLEV_NNOW_RCF .AND. tl /= nnow_rcf(jg)) RETURN
      END IF

      IF (ASSOCIATED(get_val(fields, TRIM(get_var_name(var%info))))) THEN
        CALL finish(routine, 'Variable ' // TRIM(var%info%name) // ' already saved!')
      END IF

      IF (ldebug .AND. my_process_is_stdio()) THEN
        IF (tl > 0) THEN
          WRITE (message_text,'(a,i2,a,a,a,i2)') 'DOM ', jg, ': save field ', TRIM(get_var_name(var%info)), ' TL ', tl
        ELSE
          WRITE (message_text,'(a,i2,a,a)')      'DOM ', jg, ': save field ', TRIM(get_var_name(var%info))
        ENDIF
        CALL message(routine,message_text)
      ENDIF

#if defined(__PGI) || defined(__FLANG)
      ALLOCATE(key_p)
      key_p%c = TRIM(get_var_name(var%info))
      key => key_p
#else
      ALLOCATE(key, SOURCE=TRIM(get_var_name(var%info)))
#endif
      ALLOCATE(field)

      SELECT CASE (var%info%data_type)
      CASE (REAL_T)
        CALL field%put(var%r_ptr)
      CASE (SINGLE_T)
        CALL field%put(var%s_ptr)
      CASE (INT_T)
        CALL field%put(var%i_ptr)
      CASE (BOOL_T)
        CALL field%put(var%l_ptr)
      END SELECT

      val => field
      CALL fields%setEntry(key, val)

    END SUBROUTINE save

  END SUBROUTINE save_var_group_state



  !>
  !! Restore state of a variable group. E.g. used for iterative IAU (see mo_iau.f90).
  !!
  SUBROUTINE restore_var_group_state (group_name, p_patch, fields)
    CHARACTER(len=*),           INTENT(IN) :: group_name  !< Name of the group to save.
    TYPE(t_patch),              INTENT(IN) :: p_patch     !< current patch
    TYPE(t_HashTable), POINTER, INTENT(IN) :: fields      !< storage

    CHARACTER(*), PARAMETER :: routine = 'restore_var_group_state'

    TYPE(t_vl_register_iter) :: iter
    INTEGER :: group_id

    INTEGER :: i, jg
    CHARACTER(len=vname_len) :: savedGroup(1000)
    INTEGER :: savedGroupSize

    IF (.NOT. p_patch%ldom_active) RETURN

    ! get current patch ID
    jg = p_patch%id

    IF (msg_level >= 10 .AND. my_process_is_stdio()) THEN
      CALL vlr_group(group_name, savedGroup, savedGroupSize, loutputvars_only=.FALSE., &
        &            lremap_lonlat=.FALSE., opt_lskip_container=.FALSE.)
      !
      WRITE (0,*) " "
      WRITE (0,*) "---------------------------------------"
      WRITE (0,'(a,a,a,i2)') " restore_var_group_state: Variable group info for group ", TRIM(group_name), " on DOM ", jg
      WRITE (0,*) "---------------------------------------"
      !
      CALL pretty_print_string_list(savedGroup(1:savedGroupSize))
      WRITE (0,*) " "
    ENDIF


    group_id = var_groups_dyn%group_id(TRIM(group_name))

    DO WHILE (iter%next())

      IF (iter%cur%p%patch_id /= jg) CYCLE    ! skip lists which do not belong to current patch

      ! The selection of variables is carried out by all threads redundantly. The actual work-
      ! sharing takes place inside `field%get`, which calls `copy`. The call to the nested
      ! subroutine is a workaround to some compilers (NEC) refusing CLASS(*) pointers in private
      ! clauses.
      !$OMP PARALLEL PRIVATE(i)
      DO i = 1, iter%cur%p%nvars
        CALL restore(fields, iter%cur%p%vl(i)%p)
      END DO
      !$OMP END PARALLEL
    END DO

  CONTAINS

    SUBROUTINE restore(fields, var)
      TYPE(t_HashTable), INTENT(IN) :: fields
      TYPE(t_var), INTENT(INOUT) :: var

      CLASS(*), POINTER :: val
      TYPE(t_saved_field), POINTER :: field
      INTEGER :: tl

      IF (.NOT. var%info%in_group(group_id)) RETURN

      tl = get_var_timelevel(var%info%name)

      IF (tl > 0) THEN
        IF (var%info%tlev_source == TLEV_NNOW .AND. tl /= nnow(jg)) RETURN
        IF (var%info%tlev_source == TLEV_NNOW_RCF .AND. tl /= nnow_rcf(jg)) RETURN
      END IF

      field => NULL()
      val => get_val(fields, TRIM(get_var_name(var%info)))

      IF (.NOT. ASSOCIATED(val)) THEN
        CALL finish(routine, 'Variable ' // TRIM(var%info%name) // ' has not been saved!')
      END IF

      IF (ldebug) THEN
        !$OMP SINGLE
        IF (tl > 0) THEN
          WRITE (message_text,'(a,i2,a,a,a,i2)') 'DOM ', jg, ': restore field ', TRIM(get_var_name(var%info)), ' TL ', tl
        ELSE
          WRITE (message_text,'(a,i2,a,a)')      'DOM ', jg, ': restore field ', TRIM(get_var_name(var%info))
        ENDIF
        CALL message(routine,message_text)
        !$OMP END SINGLE NOWAIT
      ENDIF


      SELECT TYPE (val)
      TYPE IS (t_saved_field)
        field => val
      END SELECT

      SELECT CASE (var%info%data_type)
      CASE (REAL_T)
        CALL field%get(var%r_ptr)
      CASE (SINGLE_T)
        CALL field%get(var%s_ptr)
      CASE (INT_T)
        CALL field%get(var%i_ptr)
      CASE (BOOL_T)
        CALL field%get(var%l_ptr)
      END SELECT

    END SUBROUTINE restore

  END SUBROUTINE restore_var_group_state


  !>
  !! Reinitialize state of a variable group. E.g. used for iterative IAU (see mo_iau.f90).
  !!
  SUBROUTINE reinit_var_group_state (group_name, p_patch)
    CHARACTER(len=*),    INTENT(IN) :: group_name  !< Name of the group to save.
    TYPE(t_patch),       INTENT(IN) :: p_patch     !< current patch

    CHARACTER(*), PARAMETER :: routine = 'reinit_var_group_state'

    TYPE(t_vl_register_iter) :: iter
    INTEGER :: group_id
    INTEGER :: i, jg
    CHARACTER(len=vname_len) :: savedGroup(1000)
    INTEGER :: savedGroupSize

    IF (.NOT. p_patch%ldom_active) RETURN

    ! get current patch ID
    jg = p_patch%id

    IF (msg_level >= 10 .AND. my_process_is_stdio()) THEN
      CALL vlr_group(group_name, savedGroup, savedGroupSize, loutputvars_only=.FALSE., &
        &            lremap_lonlat=.FALSE., opt_lskip_container=.FALSE.)
      !
      WRITE (0,*) " "
      WRITE (0,*) "---------------------------------------"
      WRITE (0,'(a,a,a,i2)') " reinit_var_group_state: Variable group info for group ", TRIM(group_name), " on DOM ", jg
      WRITE (0,*) "---------------------------------------"
      !
      CALL pretty_print_string_list(savedGroup(1:savedGroupSize))
      WRITE (0,*) " "
    ENDIF


    group_id = var_groups_dyn%group_id(TRIM(group_name))

    DO WHILE (iter%next())

      IF (iter%cur%p%patch_id /= jg) CYCLE    ! skip lists which do not belong to current patch

      ! The selection of variables is carried out by all threads redundantly. The actual work-
      ! sharing takes place inside the init routine.
      !
      !$OMP PARALLEL PRIVATE(i)
      DO i = 1, iter%cur%p%nvars
        CALL reinit(iter%cur%p%vl(i)%p)
      END DO
      !$OMP END PARALLEL
    END DO

  CONTAINS

    SUBROUTINE reinit(var)
      TYPE(t_var), INTENT(INOUT) :: var
      INTEGER :: tl

      IF (.NOT. var%info%in_group(group_id)) RETURN

      tl = get_var_timelevel(var%info%name)

      IF (tl > 0) THEN
        IF (var%info%tlev_source == TLEV_NNOW .AND. tl /= nnow(jg)) RETURN
        IF (var%info%tlev_source == TLEV_NNOW_RCF .AND. tl /= nnow_rcf(jg)) RETURN
      END IF

      IF (ldebug) THEN
        !$OMP SINGLE
        IF (tl > 0) THEN
          WRITE (message_text,'(a,i2,a,a,a,i2)') 'DOM ', jg, ': reinitialize field ', TRIM(get_var_name(var%info)), ' TL ', tl
        ELSE
          WRITE (message_text,'(a,i2,a,a)')      'DOM ', jg, ': reinitialize field ', TRIM(get_var_name(var%info))
        ENDIF
        CALL message(routine,message_text)
        !$OMP END SINGLE NOWAIT
      ENDIF

      SELECT CASE (var%info%data_type)
      CASE (REAL_T)
        CALL init(var%r_ptr, var%info%initval%rval)
      CASE (SINGLE_T)
        CALL init(var%s_ptr, var%info%initval%sval)
      CASE (INT_T)
        CALL init(var%i_ptr, var%info%initval%ival)
      CASE (BOOL_T)
        CALL init(var%l_ptr, var%info%initval%lval)
      END SELECT

    END SUBROUTINE reinit

  END SUBROUTINE reinit_var_group_state


  !> Check allocation status of a saved field.
  LOGICAL FUNCTION t_saved_field_is_allocated (self) RESULT(res)
    CLASS(t_saved_field), INTENT(IN) :: self

    res = ANY([ ALLOCATED(self%r), ALLOCATED(self%s), ALLOCATED(self%i), ALLOCATED(self%l) ])
  END FUNCTION t_saved_field_is_allocated


  !> Put a double-precision field into the saved field.
  SUBROUTINE t_saved_field_put_r (self, r)
    CLASS(t_saved_field), INTENT(INOUT) :: self
    REAL(dp), INTENT(IN) :: r(:,:,:,:,:)
    INTEGER :: ierrstat

    IF (self%is_allocated()) CALL finish('t_saved_field%put_r', 'Field is in use.')

    ALLOCATE(self%r, MOLD=r, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish('t_saved_field%put_r', 'Allocation failed for self%r')
    !
    !$OMP PARALLEL
    CALL copy(r,self%r)
    !$OMP END PARALLEL
  END SUBROUTINE t_saved_field_put_r

  !> Put a single-precision field into the saved field.
  SUBROUTINE t_saved_field_put_s (self, s)
    CLASS(t_saved_field), INTENT(INOUT) :: self
    REAL(sp), INTENT(IN) :: s(:,:,:,:,:)
    INTEGER :: ierrstat

    IF (self%is_allocated()) CALL finish('t_saved_field%put_s', 'Field is in use.')

    ALLOCATE(self%s, MOLD=s, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish('t_saved_field%put_s', 'Allocation failed for self%s')
    !
    !$OMP PARALLEL
    CALL copy(s,self%s)
    !$OMP END PARALLEL
  END SUBROUTINE t_saved_field_put_s

  !> Put an integer field into the saved field.
  SUBROUTINE t_saved_field_put_i (self, i)
    CLASS(t_saved_field), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: i(:,:,:,:,:)
    INTEGER :: ierrstat

    IF (self%is_allocated()) CALL finish('t_saved_field%put_i', 'Field is in use.')

    ALLOCATE(self%i, MOLD=i, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish('t_saved_field%put_i', 'Allocation failed for self%i')
    !
    !$OMP PARALLEL
    CALL copy(i,self%i)
    !$OMP END PARALLEL
  END SUBROUTINE t_saved_field_put_i

  !> Put a logical field into the saved field.
  SUBROUTINE t_saved_field_put_l (self, l)
    CLASS(t_saved_field), INTENT(INOUT) :: self
    LOGICAL, INTENT(IN) :: l(:,:,:,:,:)
    INTEGER :: ierrstat

    IF (self%is_allocated()) CALL finish('t_saved_field%put_l', 'Field is in use.')

    ALLOCATE(self%l, MOLD=l, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish('t_saved_field%put_l', 'Allocation failed for self%l')
    !
    !$OMP PARALLEL
    CALL copy(l,self%l)
    !$OMP END PARALLEL
  END SUBROUTINE t_saved_field_put_l

  !> Retrieve a double-precision field from the saved field.
  SUBROUTINE t_saved_field_get_r (self, r)
    CLASS(t_saved_field), INTENT(IN) :: self
    REAL(dp), INTENT(OUT) :: r(:,:,:,:,:)

    CHARACTER(len=*), PARAMETER :: procedure_name = 't_saved_field%get_r'

    IF (.NOT. self%is_allocated()) CALL finish(procedure_name, 'Field is not in use.')
    IF (.NOT. ALLOCATED(self%r)) CALL finish(procedure_name, 'Wrong data type.')

    CALL copy(self%r, r)
  END SUBROUTINE t_saved_field_get_r

  !> Retrieve a double-precision field from the saved field.
  SUBROUTINE t_saved_field_get_s (self, s)
    CLASS(t_saved_field), INTENT(IN) :: self
    REAL(sp), INTENT(OUT) :: s(:,:,:,:,:)

    CHARACTER(len=*), PARAMETER :: procedure_name = 't_saved_field%get_s'

    IF (.NOT. self%is_allocated()) CALL finish(procedure_name, 'Field is not in use.')
    IF (.NOT. ALLOCATED(self%s)) CALL finish(procedure_name, 'Wrong data type.')

    CALL copy(self%s, s)
  END SUBROUTINE t_saved_field_get_s

  !> Retrieve an integer field from the saved field.
  SUBROUTINE t_saved_field_get_i (self, i)
    CLASS(t_saved_field), INTENT(IN) :: self
    INTEGER, INTENT(OUT) :: i(:,:,:,:,:)

    CALL copy(self%i, i)
  END SUBROUTINE t_saved_field_get_i

  !> Retrieve a logical field from the saved field.
  SUBROUTINE t_saved_field_get_l (self, l)
    CLASS(t_saved_field), INTENT(IN) :: self
    LOGICAL, INTENT(OUT) :: l(:,:,:,:,:)

    CHARACTER(len=*), PARAMETER :: procedure_name = 't_saved_field%get_l'

    IF (.NOT. self%is_allocated()) CALL finish(procedure_name, 'Field is not in use.')
    IF (.NOT. ALLOCATED(self%l)) CALL finish(procedure_name, 'Wrong data type.')

    CALL copy(self%l, l)
  END SUBROUTINE t_saved_field_get_l

END MODULE mo_save_restore
