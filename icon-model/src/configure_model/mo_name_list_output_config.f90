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

! @brief configuration setup for output name lists
!
! Configuration setup for output name lists
!
! @note This is only a preliminary implementation of a configure state
!       for output name lists; based on R. Johanni's name list
!       handling which was originally implemented in
!       "shared/mo_name_list_output"

MODULE mo_name_list_output_config

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish
  USE mo_impl_constants,        ONLY: max_var_ml, max_var_pl, &
    &                                 max_var_hl, max_var_il, &
    &                                 MAX_TIME_LEVELS, SUCCESS
  USE mo_grid_config,           ONLY: n_dom
  USE mo_cdi,                   ONLY: FILETYPE_GRB, FILETYPE_GRB2
  USE mo_var_metadata_types,    ONLY: t_var_metadata
  USE mo_util_string,           ONLY: toupper
  USE mo_name_list_output_types,ONLY: t_output_name_list
  USE mo_key_value_store,               ONLY: t_key_value_store

  IMPLICIT NONE

  PUBLIC :: is_grib_output,                                  &
    &       is_variable_in_output,                           &
    &       is_variable_in_output_dom,                       &
    &       is_variable_in_output_nml
  PUBLIC :: use_async_name_list_io
  PUBLIC :: first_output_name_list

  ! Pointer to a linked list of output name lists:
  TYPE(t_output_name_list), POINTER :: first_output_name_list => NULL()

  ! Flag whether async name_list I/O is used, it is set in the main program:
  LOGICAL :: use_async_name_list_io = .FALSE.

  TYPE :: t_storage_array
    TYPE(t_key_value_store), POINTER :: ptr
  END TYPE t_storage_array

  !> auxiliary data structures: (case insensitive) hash table containing
  !  all output variable names, either for all model domains or domain-specific
  TYPE(t_key_value_store), POINTER, PRIVATE       :: output_variables => NULL()
  TYPE(t_storage_array), POINTER, PRIVATE :: outputvar_dom(:) => NULL()

CONTAINS
  
  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if one of the output namelists has been specified with GRIB output.
  FUNCTION is_grib_output() RESULT(retval)
    LOGICAL                           :: retval
    TYPE(t_output_name_list), POINTER :: p_onl

    retval = .FALSE.
    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      retval = retval                           .OR.  &
        &      (p_onl%filetype == FILETYPE_GRB) .OR.  &
        &      (p_onl%filetype == FILETYPE_GRB2)
      p_onl => p_onl%next
    ENDDO
  END FUNCTION is_grib_output


  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output for a given variable is due for a given
  !!         namelist
  !!
  FUNCTION is_variable_in_output_nml(p_onl, var_name) RESULT(retval)
    LOGICAL :: retval

    TYPE(t_output_name_list), POINTER         :: p_onl      !< output name list
    CHARACTER(LEN=*), INTENT(IN) :: var_name   !< variable name
    ! local variables
    INTEGER :: ivar
    CHARACTER(len=LEN(var_name)) :: var_name_u

    ! if a specific variable name has been provided, loop over the
    ! variables for this output file
    retval = .FALSE.
    var_name_u = toupper(var_name)
    DO ivar=1,max_var_ml
      IF (p_onl%ml_varlist(ivar) /= ' ') THEN
        retval = toupper(p_onl%ml_varlist(ivar)) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
    DO ivar=1,max_var_pl
      IF (p_onl%pl_varlist(ivar) /= ' ') THEN
        retval = toupper(p_onl%pl_varlist(ivar)) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
    DO ivar=1,max_var_hl
      IF (p_onl%hl_varlist(ivar) /= ' ') THEN
        retval = toupper(p_onl%hl_varlist(ivar)) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
    DO ivar=1,max_var_il
      IF (p_onl%il_varlist(ivar) /= ' ') THEN
        retval = toupper(p_onl%il_varlist(ivar)) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
  END FUNCTION is_variable_in_output_nml


  !-------------------------------------------------------------------------------------------------
  !>
  !! Initialize hash table for efficient variable search
  !!
  SUBROUTINE init_hashtable(first_output_name_list,ptr_outvar,jg)

    !> head output namelist list
    TYPE(t_output_name_list), POINTER :: first_output_name_list
    TYPE(t_key_value_store), INTENT(INOUT)  :: ptr_outvar
    INTEGER, INTENT(IN), OPTIONAL :: jg

    ! local variables
    INTEGER :: ivar, ig
    TYPE (t_output_name_list), POINTER :: p_onl

    IF (PRESENT(jg)) THEN
      ig = jg
    ELSE
      ig = -1
    ENDIF

    CALL ptr_outvar%init(.FALSE.)
    p_onl => first_output_name_list
    DO WHILE (ASSOCIATED(p_onl))
      IF (p_onl%dom == ig .OR. ig == -1 .OR. p_onl%dom == -1) THEN
        DO ivar=1,max_var_ml
          IF (p_onl%ml_varlist(ivar) /= ' ')  CALL ptr_outvar%put(p_onl%ml_varlist(ivar), .TRUE.)
        END DO
        DO ivar=1,max_var_pl
          IF (p_onl%pl_varlist(ivar) /= ' ')  CALL ptr_outvar%put(p_onl%pl_varlist(ivar), .TRUE.)
        END DO
        DO ivar=1,max_var_hl
          IF (p_onl%hl_varlist(ivar) /= ' ')  CALL ptr_outvar%put(p_onl%hl_varlist(ivar), .TRUE.)
        END DO
        DO ivar=1,max_var_il
          IF (p_onl%il_varlist(ivar) /= ' ')  CALL ptr_outvar%put(p_onl%il_varlist(ivar), .TRUE.)
        END DO
      ENDIF
      p_onl => p_onl%next
    END DO

  END SUBROUTINE init_hashtable

  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output for a given variable is due in any
  !!         output namelist.
  !!
  !! @note This function will only return .TRUE. if the variable is given explicitely in an output
  !!       namelist, not if the variable is part of a group that's given in an output namelist!
  !!       In particular, if it is called in the memory initialization routines it will not catch
  !!       these cases. However, if you know that a variable is in a specific group "gname" this
  !!       function can be called with "group:gname" as var_name and it will return .TRUE. if
  !!       "group:gname" is in an output namelist. Same with "tiles:".
  !!
  FUNCTION is_variable_in_output(var_name) RESULT(retval)
    LOGICAL :: retval

    !> head output namelist list
    CHARACTER(LEN=*), INTENT(IN)  :: var_name   !< variable name

    ! local variables
    INTEGER :: ierror
    LOGICAL :: lval

    ! if called for the first time: set up (case insensitive) hash
    ! table containing all output variable names
    IF (.NOT. ASSOCIATED(output_variables)) THEN
      ALLOCATE(output_variables)
      CALL init_hashtable(first_output_name_list,output_variables)
    ENDIF

    ! hashtable look-up if "var_name" exists:
    CALL output_variables%get(var_name, lval, opt_err=ierror)
    retval = (ierror == SUCCESS)

  END FUNCTION is_variable_in_output

  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output for a given variable is due for a specific domain in any
  !!         output namelist.
  !!
  FUNCTION is_variable_in_output_dom(var_name, jg) RESULT(retval)
    LOGICAL :: retval

    CHARACTER(LEN=*), INTENT(IN)  :: var_name   !< variable name
    INTEGER, INTENT(in)           :: jg         !< domain index

    ! local variables
    INTEGER :: ierror, ig
    LOGICAL :: lval
    CHARACTER(1) :: jg_str

    ! if called for the first time: set up (case insensitive) hash
    ! table containing all output variable names
    IF (.NOT. ASSOCIATED(outputvar_dom)) THEN
      ALLOCATE(outputvar_dom(n_dom))
      DO ig = 1, n_dom
        ALLOCATE(outputvar_dom(ig)%ptr)
        CALL init_hashtable(first_output_name_list,outputvar_dom(ig)%ptr,ig)
      ENDDO
    ENDIF

    ! hashtable look-up if "var_name" exists:
    CALL outputvar_dom(jg)%ptr%get(var_name, lval, opt_err=ierror)
    retval = (ierror == SUCCESS)

  END FUNCTION is_variable_in_output_dom

END MODULE mo_name_list_output_config
