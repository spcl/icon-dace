! Module handling the selection of vertical levels for the output
! module.
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
!
! --------------------------------
!    Details of the implementation
! --------------------------------
!
! Derived data type "t_level_selection":
!
! Vertical levels are selected via objects of the derived data type
! "t_level_selection".  If such a data object has been initialized
! for the output file, then the vertical axis definitions below
! create CDI axis objects for the selected levels only.
!
! The write routines in the module "mo_name_list_output" skip levels
! which are not part of the selection.
!
! Furthermore, if a "t_level_selection" object is present, then the
! memory windows for the asynchronous one-sided MPI communication are
! adjusted to the reduced level size
! "output_file%level_selection%n_selected".
!
! Creation of level selection data:
!
! During the setup phase, level selection objects are created in two
! different situations:
!
! a) Definition of namelist parameter "m_levels"
!
!    The user may specify levels and/or level ranges in the form of a
!    (string) namelist parameter. This string is parsed with the help
!    of the module "mo_util_string_parse" and converted into a
!    "t_level_selection" object.
!
! b) Definition of namelist parameters "p_levels", "h_levels", "i_levels"
!
!    Each output namelist "output_nml" may specify its own range of
!    pressure levels for output (or for height, isentropic levels as
!    well). However, the internal post-processing routines perform
!    vertical interpolation for the union set of all requested
!    pressure levels of a specific domain at once (merging the
!    information from several "output_nml" namelists. Afterwards, the
!    write routine for the output files merely copies the respective
!    levels. For this purpose, the described "t_level_selection"
!    mechanism is applied as well.

MODULE mo_level_selection

  USE mo_kind,                              ONLY: wp
  USE mo_impl_constants,                    ONLY: SUCCESS
  USE mo_var, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_exception,                         ONLY: finish
  USE mo_name_list_output_types,            ONLY: t_output_file
  USE mo_math_utilities,                    ONLY: t_value_set, find_values_in_set
  USE mo_run_config,                        ONLY: num_lev
  USE mo_util_string_parse,                 ONLY: util_do_parse_intlist
  USE mo_level_selection_types,             ONLY: t_level_selection

  USE mo_nh_pzlev_config,                   ONLY: nh_pzlev_config


  IMPLICIT NONE

  PRIVATE

  ! subroutines
  PUBLIC :: create_mipz_level_selections

  INTERFACE create_level_selection
    MODULE PROCEDURE create_level_selection_str
    MODULE PROCEDURE create_level_selection_set
  END INTERFACE

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_level_selection'


CONTAINS

  ! --------------------------------------------------------------------------------------
  !> Creates a level selection data object from a selection that is
  !  described as a character string, e.g. "1,5...10,15". See the
  !  module "mo_util_string_parse" for a detailed description of valid
  !  selection strings.
  !
  SUBROUTINE create_level_selection_str(selection_str, nlevs, level_selection, opt_nlev_value)
    CHARACTER(LEN=*),        INTENT(IN)    :: selection_str    !< selection described as string
    INTEGER,                 INTENT(IN)    :: nlevs            !< total no. of levels
    TYPE(t_level_selection), INTENT(INOUT), POINTER :: level_selection
    INTEGER, INTENT(IN), OPTIONAL :: opt_nlev_value            !< number to substitute for "N"/"nlev"
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::create_level_selection_str'
    INTEGER :: int_list(0:nlevs) ! note: 0-lower bound required
    INTEGER :: ierrstat, nlev_value

    IF (TRIM(selection_str) == "") RETURN ! do nothing

    ALLOCATE(level_selection)
    ALLOCATE(level_selection%s(nlevs), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! parse character string into a LOGICAL array:
    nlev_value = nlevs
    IF (PRESENT(opt_nlev_value)) nlev_value = opt_nlev_value
    CALL util_do_parse_intlist(selection_str, nlev_value, int_list, ierrstat)
    IF (ierrstat /= 0) CALL finish(routine, 'Parsing of level selection failed for "'//TRIM(selection_str)//'"')
    level_selection%s(1:nlevs) = (int_list(1:nlevs) == 1)
    ! count no. of selected levels
    level_selection%n_selected = COUNT(level_selection%s)
    ! get mapping: global level -> local level index
    ALLOCATE(level_selection%global_idx(level_selection%n_selected), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_level_indices_in_selection(level_selection%s, level_selection%n_selected, &
      &                                 level_selection%global_idx)
    ! get mapping: local level -> global level index
    ALLOCATE(level_selection%local_idx(nlevs), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_selection_index(level_selection%s, level_selection%local_idx)
  END SUBROUTINE create_level_selection_str


  ! --------------------------------------------------------------------------------------
  !  > Creates a level selection data object from a selection that is
  !    described as list of REAL(wp) that partly match the levels
  !    contained in a "t_value_set" object.
  !
  !  @note A list of level values starting with a negative value means
  !        that no level selection is to be created (this is the
  !        default).
  !
  SUBROUTINE create_level_selection_set(selection, value_set, level_selection)
    REAL(wp),                INTENT(IN)    :: selection(:)     !< selected levels
    TYPE (t_value_set),      INTENT(IN)    :: value_set        !< total set of available levels
    TYPE(t_level_selection), INTENT(INOUT), POINTER :: level_selection
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::create_level_selection_set'
    INTEGER :: ierrstat, nlevs

    IF (selection(1) < 0._wp)  RETURN ! do nothing
    ! count the no. of levels
    DO nlevs=1,SIZE(selection)
      IF (selection(nlevs) < 0._wp) EXIT
    END DO
    nlevs = nlevs - 1
    ! allocate the level selection object
    ALLOCATE(level_selection)
    ALLOCATE(level_selection%s(value_set%nvalues), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! determine the selected levels
    CALL find_values_in_set(nlevs, selection, value_set, level_selection%s)

    ! count no. of selected levels
    level_selection%n_selected = COUNT(level_selection%s)
    ! get mapping: global level -> local level index
    ALLOCATE(level_selection%global_idx(level_selection%n_selected), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_level_indices_in_selection(level_selection%s, level_selection%n_selected, &
      &                                 level_selection%global_idx)
    ! get mapping: local level -> global level index
    ALLOCATE(level_selection%local_idx(value_set%nvalues), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL get_selection_index(level_selection%s, level_selection%local_idx)
  END SUBROUTINE create_level_selection_set



  ! --------------------------------------------------------------------------------------
  !> Utility function: Assume that we have N vertical levels, out of
  !  which only a few levels are selected. This subroutine takes this
  !  selection as input in the form of a LOGICAL array s(1...N), where
  !  "s(i)=.TRUE." means that level "i" is selected. As an output, we
  !  get an integer list idx(1...n_selected) containing the selected
  !  level indices.
  !
  SUBROUTINE get_level_indices_in_selection(s, n_selected, idx)
    LOGICAL, INTENT(IN)  :: s(:)       !< level selection (LOGICAL array)
    INTEGER, INTENT(OUT) :: n_selected !< no. of selected levels
    INTEGER, INTENT(OUT) :: idx(:)     !< selected level indices
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_level_indices_in_selection'
    INTEGER :: nlevs, i

    nlevs = SIZE(s)
    n_selected = 0
    DO i=1,nlevs
      IF (s(i)) THEN
        n_selected = n_selected + 1
        IF (SIZE(idx) < n_selected)  CALL finish(routine, "Dimension mismatch!")
        idx(n_selected) = i
      END IF
    END DO
  END SUBROUTINE get_level_indices_in_selection


  ! --------------------------------------------------------------------------------------
  !> Utility function: Assume that we have N vertical levels, out of
  !  which only a few levels are selected. This subroutine takes this
  !  selection as input in the form of a LOGICAL array s(1...N), where
  !  "s(i)=.TRUE." means that level "i" is selected. As an output, we
  !  get an integer list idx(1...N) containing the local index in the
  !  list of selected level indices (i.e. an integer number in the
  !  range 1...n_selected).
  !
  SUBROUTINE get_selection_index(s, idx)
    LOGICAL, INTENT(IN)  :: s(:)       !< level selection (LOGICAL array)
    INTEGER, INTENT(OUT) :: idx(:)     !< selected level indices
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_selection_index'
    INTEGER :: nlevs, i, n_selected

    nlevs = SIZE(s)
    n_selected = 0
    DO i=1,nlevs
      IF (s(i)) THEN
        n_selected = n_selected + 1
        idx(i)     = n_selected
      ELSE
        idx(i)     = 0
      END IF
    END DO
  END SUBROUTINE get_selection_index


  !------------------------------------------------------------------------------------------------
  ! Loop over output file (p_of) and create the "selection" from the
  ! union set of vertical model/pressure/isentropic/height (MIPZ)
  ! levels:
  SUBROUTINE create_mipz_level_selections(output_file)
    TYPE(t_output_file), TARGET, INTENT(INOUT) :: output_file(:)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::create_mipz_level_selections'
    TYPE (t_output_file), POINTER   :: p_of
    INTEGER :: i, log_patch_id

    DO i=1,SIZE(output_file)
      p_of => output_file(i)
      log_patch_id = p_of%log_patch_id
      SELECT CASE(p_of%ilev_type)
      CASE (level_type_ml)
        ! note: "ml" comprises both full and half levels
        CALL create_level_selection(p_of%name_list%m_levels,              &
          &                         num_lev(log_patch_id)+1, p_of%level_selection,   &
          &                         opt_nlev_value = num_lev(log_patch_id))

      CASE (level_type_pl)
        CALL create_level_selection(p_of%name_list%p_levels, &
          &                         nh_pzlev_config(log_patch_id)%plevels, p_of%level_selection)
      CASE (level_type_hl)
        CALL create_level_selection(p_of%name_list%z_levels, &
          &                         nh_pzlev_config(log_patch_id)%zlevels, p_of%level_selection)
      CASE (level_type_il)
        CALL create_level_selection(p_of%name_list%i_levels, &
          &                         nh_pzlev_config(log_patch_id)%ilevels, p_of%level_selection)

      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT
    END DO

  END SUBROUTINE create_mipz_level_selections

END MODULE mo_level_selection

