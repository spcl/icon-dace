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

MODULE mo_var_list_register
  USE mo_var_list,         ONLY: t_var_list_ptr
  USE mo_exception,        ONLY: message, finish
  USE mo_io_config,        ONLY: restart_file_type
  USE mo_key_value_store,  ONLY: t_key_value_store

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_vl_register_iter
  PUBLIC :: vlr_add, vlr_get, vlr_del
  PUBLIC :: get_nvl

  TYPE t_vl_register_iter
    INTEGER, PRIVATE :: cur_idx = -1
    TYPE(t_var_list_ptr), PUBLIC :: cur
    LOGICAL, PRIVATE :: anew = .TRUE.
  CONTAINS
    PROCEDURE, PUBLIC :: next => iter_next
  END TYPE t_vl_register_iter

  INTEGER, SAVE :: last_id = 0
  TYPE(t_key_value_store) :: map
  TYPE(t_var_list_ptr), ALLOCATABLE :: storage(:)
  CHARACTER(*), PARAMETER :: modname = 'mo_var_list_store'

CONTAINS
  !------------------------------------------------------------------------------------------------
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  SUBROUTINE vlr_add(list, vlname, model_type, output_type, restart_type,    &
    & post_suf, rest_suf, init_suf, loutput, lrestart, linitial, patch_id,   &
    & vlevel_type, filename, compression_type)
    TYPE(t_var_list_ptr), INTENT(OUT) :: list    ! anchor
    CHARACTER(*), INTENT(IN) :: vlname         ! name of output var_list
    CHARACTER(*), INTENT(IN) :: model_type
    INTEGER, INTENT(IN), OPTIONAL :: output_type, restart_type, patch_id, vlevel_type, compression_type
    CHARACTER(*), INTENT(IN), OPTIONAL :: post_suf, rest_suf, init_suf, filename
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, linitial  ! in standard output/restart/initial file
    INTEGER :: vln_len, ivl
    TYPE(t_var_list_ptr), ALLOCATABLE :: tmp_stor(:)

    vln_len = LEN_TRIM(vlname)
    CALL message('','')
    CALL message('','adding new var_list '//vlname(1:vln_len))
    IF (last_id .EQ. 0) THEN
      ALLOCATE(storage(8))
      CALL map%init(.false.)
    END IF
    CALL vlr_get(list, vlname)
    IF (ASSOCIATED(list%p)) CALL finish('new_var_list', ' >'//vlname(1:vln_len)//'< already in use.')
    ALLOCATE(list%p)
    last_id = last_id + 1
    CALL map%put(vlname, last_id)
    ! set default list characteristics
    list%p%vlname     = vlname(1:vln_len)
    list%p%model_type = model_type
    list%p%post_suf   = '_'//vlname(1:vln_len)
    list%p%rest_suf   = list%p%post_suf
    list%p%init_suf   = list%p%post_suf
    ! set non-default list characteristics
    list%p%restart_type = restart_file_type
    IF (PRESENT(output_type))      list%p%output_type      = output_type
    IF (PRESENT(restart_type))     list%p%restart_type     = restart_type
    IF (PRESENT(post_suf))         list%p%post_suf         = post_suf
    IF (PRESENT(rest_suf))         list%p%rest_suf         = rest_suf
    IF (PRESENT(init_suf))         list%p%init_suf         = init_suf
    IF (PRESENT(loutput))          list%p%loutput          = loutput
    IF (PRESENT(lrestart))         list%p%lrestart         = lrestart
    IF (PRESENT(linitial))         list%p%linitial         = linitial
    IF (PRESENT(patch_id))         list%p%patch_id         = patch_id
    IF (PRESENT(vlevel_type))      list%p%vlevel_type      = vlevel_type
    IF (PRESENT(filename))         list%p%filename         = filename
    IF (PRESENT(compression_type)) list%p%compression_type = compression_type
    IF (SIZE(storage) .LT. last_id) THEN
      CALL MOVE_ALLOC(storage, tmp_stor)
      ALLOCATE(storage(SIZE(tmp_stor) + 8))
      DO ivl = 1, SIZE(tmp_stor)
        IF (ASSOCIATED(tmp_stor(ivl)%p)) storage(ivl)%p => tmp_stor(ivl)%p
      END DO
    END IF
    storage(last_id)%p => list%p
  END SUBROUTINE vlr_add

  ! Get a reference to a memory buffer/output var_list
  SUBROUTINE vlr_get(list, vlname)
    TYPE(t_var_list_ptr), INTENT(OUT) :: list ! pointer
    CHARACTER(*), INTENT(IN) :: vlname
    INTEGER :: ivl, ierr

    CALL map%get(vlname, ivl, ierr)
    NULLIFY(list%p)
    IF (ierr .EQ. 0) THEN
      IF (ASSOCIATED(storage(ivl)%p)) list%p => storage(ivl)%p
    END IF
  END SUBROUTINE vlr_get

  LOGICAL FUNCTION iter_next(this) RESULT(valid)
    CLASS(t_vl_register_iter), INTENT(INOUT) :: this

    valid = .false.
    IF (last_id .GT. 0) THEN
      IF (this%anew) THEN
        this%cur_idx = 0
        this%anew = .false.
      END IF
      DO WHILE(this%cur_idx .LT. SIZE(storage) .AND. .NOT.valid)
        this%cur_idx = this%cur_idx + 1
        IF (ASSOCIATED(storage(this%cur_idx)%p)) THEN
          valid = .true.
          this%cur%p => storage(this%cur_idx)%p
        END IF
      END DO
      IF (.NOT.valid) THEN
        NULLIFY(this%cur%p)
        this%anew = .true.
      END IF
    END IF
  END FUNCTION iter_next

  ! Delete an output var_list, nullify the associated pointer
  SUBROUTINE vlr_del(list)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: list
    INTEGER :: ivl, ierr

    IF (ASSOCIATED(list%p)) THEN
      CALL map%get(list%p%vlname, ivl, ierr)
      IF (ierr .NE. 0) CALL finish(modname//":delete_var_list", &
        & "var_list <" // TRIM(list%p%vlname) // "> not registered")
      CALL list%delete()
      DEALLOCATE(storage(ivl)%p)
      NULLIFY(list%p)
    END IF
  END SUBROUTINE vlr_del

  INTEGER FUNCTION get_nvl()

    get_nvl = last_id
  END FUNCTION get_nvl

END MODULE mo_var_list_register
