! Auxiliary module: Print a list of all output variables.
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

MODULE mo_name_list_output_printvars

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char
  USE mo_cdi,                               ONLY: cdi_max_name
  USE mo_gribout_config,                    ONLY: t_gribout_config
  USE mo_kind,                              ONLY: wp, dp
  USE mo_impl_constants,                    ONLY: SUCCESS, vname_len
  USE mo_cf_convention,                     ONLY: t_cf_var
  USE mo_exception,                         ONLY: finish, message_text
  USE mo_var_metadata_types,                ONLY: t_var_metadata
  USE mo_var_list_register,                 ONLY: t_vl_register_iter
  USE mo_var,                               ONLY: t_var, level_type_ml
  USE mo_dictionary,                        ONLY: t_dictionary
  USE mo_util_sort,                         ONLY: quicksort
  USE mo_util_string,                       ONLY: remove_duplicates, toupper, tolower, int2string

  IMPLICIT NONE
  PRIVATE

  ! subroutines
  PUBLIC :: print_var_list

  !------------------------------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_printvars'


CONTAINS

  !------------------------------------------------------------------------------------------------
  !> @return variable name without time level or tile suffix
  !
  CHARACTER(LEN=vname_len) FUNCTION get_var_basename(info)
    TYPE(t_var_metadata) :: info
    INTEGER :: endidx
    CHARACTER(LEN=1) :: suffix_str

    ! first cut off the time level suffix:
    endidx = INDEX(info%name,'.TL')
    IF (endidx .EQ. 0) THEN
      endidx = LEN_TRIM(info%name)
    ELSE
      endidx = endidx - 1
    END IF
    suffix_str = " "
    ! condense three-digit suffices "001, 002, 003, ..." into "*"
    ! for container variables
    IF (endidx > 3) THEN
      IF (info%lcontained .AND. &
          .NOT. is_number(info%name(endidx-3:endidx-3)) .AND. &
        &       is_number(info%name(endidx-2:endidx-2)) .AND. &
        &       is_number(info%name(endidx-1:endidx-1)) .AND. &
        &       is_number(info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 3
        suffix_str = "*"
      END IF
    END IF
    ! condense two-digit suffices "_01, _02, _03, ..." into "*"
    IF (endidx > 2) THEN
      IF ((info%name(endidx-2:endidx-2) == "_")    .AND. &
        &  is_number(info%name(endidx-1:endidx-1)) .AND. &
        &  is_number(info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 2
        suffix_str = "*"
      END IF
    END IF
    ! condense one-digit suffices "_1, _2, _3, ..." into "_*"
    IF ((endidx > 1) .AND. (suffix_str == " ")) THEN
      IF ((info%name(endidx-1:endidx-1) == "_") .AND. &
        &  is_number(info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 1
        suffix_str = "*"
      END IF
    END IF
    get_var_basename = info%name(1:endidx)//suffix_str

  CONTAINS
    LOGICAL FUNCTION is_number(chr)
      CHARACTER, INTENT(IN) :: chr
      is_number = (IACHAR(chr) - IACHAR('0')) <= 9
    END FUNCTION is_number
  END FUNCTION get_var_basename


  !------------------------------------------------------------------------------------------------
  !> Print list of all output variables (LaTeX table formatting).
  !
  SUBROUTINE print_var_list(out_varnames_dict, print_patch_id, gribout_config, &
    &                       i_lctype)
    TYPE(t_dictionary),     INTENT(IN) :: out_varnames_dict
    TYPE(t_gribout_config), INTENT(IN) :: gribout_config
    INTEGER,                INTENT(IN) :: i_lctype, print_patch_id

    CHARACTER(*), PARAMETER :: routine = modname//"::print_var_list"
    INTEGER,      PARAMETER :: max_str_len = cdi_max_name + 1 + 128 + vname_len + 99
    CHARACTER(*), PARAMETER :: varprefix = "\varname{", CR = " \\[0.5em]"
    CHARACTER(*), PARAMETER :: descrprefix = "\vardescr{"
    INTEGER,      PARAMETER :: PREF = LEN_TRIM(varprefix) + 1
    TYPE (t_var_metadata),    POINTER              :: info
    TYPE(t_var),     POINTER              :: elem
    TYPE(t_cf_var),           POINTER              :: this_cf
    INTEGER                                        :: i, iv, nout_vars, iout_var, ierrstat
    CHARACTER(kind=c_char, LEN = cdi_max_name + 1) :: vname
    CHARACTER(LEN=max_str_len), ALLOCATABLE        :: out_vars(:)
    CHARACTER(len=128)                             :: descr_string
    TYPE(t_vl_register_iter) :: vl_iter
    ! ---------------------------------------------------------------------------


    ! count the no. of output variables:
    nout_vars = 0
    DO WHILE(vl_iter%next())
      IF (vl_iter%cur%p%patch_id /= print_patch_id) CYCLE
      IF(.NOT.vl_iter%cur%p%loutput) CYCLE
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE
      DO iv = 1, vl_iter%cur%p%nvars
        IF (vl_iter%cur%p%vl(iv)%P%info%loutput) nout_vars = nout_vars + 1
      END DO
    END DO
    ! allocate sufficient space
    ALLOCATE(out_vars(nout_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    iout_var = 0
    DO WHILE(vl_iter%next())
      IF (vl_iter%cur%p%patch_id /= print_patch_id) CYCLE
      ! Inspect only model level variables
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE
      ! Do not inspect element if output is disabled
      IF(.NOT.vl_iter%cur%p%loutput) CYCLE
      DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)%p
        info => elem%info
        ! Do not inspect element if output is disabled
        IF (.NOT.info%loutput) CYCLE
        this_cf => info%cf
        IF (info%post_op%lnew_cf) this_cf => info%post_op%new_cf
        ! if no short is available and if the variable is a
        ! "reference" into another variable, then search for this
        ! source variable:
        IF ((LEN_TRIM(this_cf%long_name) == 0) .AND. info%lcontained) THEN
          this_cf => info%cf
          IF (ASSOCIATED(elem%ref_to)) THEN
            IF (elem%ref_to%info%ncontained > 0) this_cf => elem%ref_to%info%cf
          END IF
        END IF
        vname = info%name
        vname = tolower(vname)
        IF ((vname(1:3) == "var") .OR. (vname(1:5) == "param")) THEN
          vname = ""
        ELSE
          vname = varprefix//TRIM(vname)//"}"
        END IF

        descr_string = this_cf%long_name
        ! upcase first letter of description string:
        descr_string(1:1) = toupper(descr_string(1:1))

        iout_var = iout_var + 1
        WRITE (out_vars(iout_var),'(8a)') varprefix, &
          & tolower(get_var_basename(info)), "} & ", &
          & TRIM(vname), ' & ', descrprefix, TRIM(descr_string), "}"
      ENDDO
    ENDDO
    ! sort and remove duplicates
    CALL quicksort(out_vars(1:nout_vars))
    CALL remove_duplicates(out_vars(1:nout_vars), nout_vars)

    ! print table, but add a gap when new alphabetical letter starts:
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) "List of output variables"
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) " "
    DO i = 1, nout_vars-1
      ! check for the initial character of the variable name:
      WRITE(0,*) TRIM(out_vars(i)) // CR(1:MERGE(3,LEN(CR), &
        &  out_vars(i)(PREF:PREF) == out_vars(i+1)(PREF:PREF)))
    END DO
    WRITE (0,*) TRIM(out_vars(nout_vars)) // CR(1:3)
    WRITE (0,*) " "

    DEALLOCATE(out_vars, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE print_var_list

END MODULE mo_name_list_output_printvars
