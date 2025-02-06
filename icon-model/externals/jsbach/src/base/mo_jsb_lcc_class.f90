!> Types for lcc processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains types required for dealing with actively or passively relocated conserved quantities upon lcc
!>
MODULE mo_jsb_lcc_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_util,                ONLY: one_of
  USE mo_exception,           ONLY: finish, message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_math_constants,  ONLY: one_year

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: min_daily_cf_change, min_annual_cf_change, min_tolerated_fract_mismatch, min_tolerated_cf_mismatch,       &
    & t_jsb_lcc_proc, t_jsb_lcc_proc_p, t_jsb_lcc_var, t_jsb_lcc_var_p,                                               &
    & transfer_from_active_to_passive_var_onChunk, collect_matter_of_active_vars_onChunk,                             &
    & copy_from_active_to_passive_var_onChunk

  REAL(wp), PARAMETER :: min_daily_cf_change = 1.0E-16_wp
    !< min change in cover fractions for daily lcc, i.e. if the change is smaller than this threshold the lcc is not applied
  REAL(wp), PARAMETER :: min_annual_cf_change = min_daily_cf_change * one_year
    !< min change in cover fractions for annual lcc, i.e. if the change is smaller than this threshold the lcc is not applied
  REAL(wp), PARAMETER :: min_tolerated_cf_mismatch = 1.0E-11_wp
    !< min mismatch in cover fractions which is tolerated when comparing different cover fraction sums
  REAL(wp), PARAMETER :: min_tolerated_fract_mismatch = 1.0E-14_wp
    !< minimum tolearated mismatch in forest age fraction sums (only used when running with FAGE)

  !> Type used to enclose pointers and values for each conserved quantity with active or passive relocation
  TYPE :: t_jsb_lcc_var
    CHARACTER(len=VARNAME_LEN)    :: name             !< Name of the variable
!    INTEGER :: type_id = -1     !< one of the CQ_TYPEs ( mo_jsb_cqt_class ) -- JN-TODO: required -> also in encompassing type...
    TYPE(t_jsb_var_p), ALLOCATABLE :: var_on_tile(:)  !< Pointer to the var on each tile involved in the encompassing lcc process
    REAL(wp), ALLOCATABLE :: relocate_this(:,:,:)     !< One array (on domain) per tile
  END TYPE t_jsb_lcc_var
  TYPE :: t_jsb_lcc_var_p
    TYPE(t_jsb_lcc_var), POINTER :: p => NULL()
  END TYPE t_jsb_lcc_var_p

  !> Type used on lcc processes to collect conserved quantities with active and passive relocations
  TYPE :: t_jsb_lcc_proc
    CHARACTER(len=:), ALLOCATABLE :: name !< Name of the structure: process + _lcc
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: tile_names(:) !< Names of the tiles involved with this lcc process
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: unique_proc_names_active_vars(:)
        !< Collection of unique names of the processes carrying variables that are actively relocated by this lcc process
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: active_vars_names(:)
        !< List of active var names - same order as other lists for active vars
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: passive_vars_names(:)
        !< List of passive var names - same order as other lists for passive vars
    INTEGER, ALLOCATABLE :: active_vars_cqt(:)  ! JN-TODO: required?
    INTEGER, ALLOCATABLE :: passive_vars_cqt(:) ! JN-TODO: required?
    INTEGER, ALLOCATABLE :: active_vars_process_id(:)  !< id of the process carrying the corresponding active var
    INTEGER, ALLOCATABLE :: passive_vars_process_id(:) !< id of the process carrying the corresponding passive var
    INTEGER :: nr_of_tiles = 0          !< Number of tiles involved with this lcc process
    INTEGER :: nr_of_procs_with_active_vars = 0
        !< Number of processes carrying variables that are actively or passively relocated by this lcc process
    INTEGER :: nr_of_active_vars = 0    !< Number of variables with active relocation in the encompassing lcc process
    INTEGER :: nr_of_passive_vars = 0   !< Number of variables with passive relocation in the encompassing lcc process
    TYPE(t_jsb_lcc_var_p), ALLOCATABLE :: active_vars(:)
        !< Collection of lcc var types that are actively relocated for the encompassing lcc process
    TYPE(t_jsb_lcc_var_p), ALLOCATABLE :: passive_vars(:)
        !< Collection of lcc var types that are passively relocated for the encompassing lcc process
  END TYPE t_jsb_lcc_proc

  !> Type used to collect lcc process structures
  TYPE t_jsb_lcc_proc_p
    TYPE(t_jsb_lcc_proc), POINTER :: p
  END TYPE t_jsb_lcc_proc_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lcc_class'

CONTAINS


  ! ====================================================================================================== !
  !
  !> Relocates a certain fraction of the content of a specific given active to a specific given passive var
  !
  SUBROUTINE transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice,           &
      &                                                  active_var_name, passive_var_name, transfer_fract, &
      &                                                  conversion_fact)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing matter from an active to a passive var
    INTEGER,                    INTENT(in)    :: i_tile
        !< index of the tile in the lcc structure for which the relocation should be conducted
    INTEGER,                    INTENT(in)    :: iblk    !< Number of current block (chunk)
    INTEGER,                    INTENT(in)    :: ics     !< Index of chunk start
    INTEGER,                    INTENT(in)    :: ice     !< Index of chunk end
    CHARACTER(len=*), INTENT(in)    :: active_var_name   !< (active) source var
    CHARACTER(len=*), INTENT(in)    :: passive_var_name  !< (passive) sink var
    REAL(wp),  INTENT(in), OPTIONAL :: transfer_fract
        !< fraction of matter transferred from active to passive var - if not given: 1.0 is assumed
    REAL(wp),  INTENT(in), OPTIONAL :: conversion_fact
        !< conversion fraction if active var has different unit then passive var - if not given: 1.0 is assumed
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':transfer_from_active_to_passive_var_onChunk'

    REAL(wp) :: this_transfer_fract, this_conversion_fact
    INTEGER  :: i_active, i_passive
    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    i_active = one_of( active_var_name , lcc_relocations%active_vars_names)
    i_passive = one_of( passive_var_name , lcc_relocations%passive_vars_names)

    !>
    !> Assertion source variable needs to be active variable for this lcc process
    !>
    IF (i_active < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed active var ' // TRIM(active_var_name) &
        & // ' is not an active var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    !>
    !> Assertion sink variable needs to be passive variable for this lcc process
    !>
    IF (i_passive < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed passive var ' // TRIM(passive_var_name) &
        & // ' is not a passive var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    IF (PRESENT(transfer_fract)) THEN
      this_transfer_fract = transfer_fract
    ELSE
      this_transfer_fract = 1.0
    ENDIF

    IF (PRESENT(conversion_fact)) THEN
      this_conversion_fact = conversion_fact
    ELSE
      this_conversion_fact = 1.0
    ENDIF

    ! Add fraction of active content to be transferred to passive var,
    ! potentially converting to the unit of the passive var
    lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk)       &
      & = lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & + (lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)  &
      &    * this_transfer_fract * this_conversion_fact)

    ! And subtract from active - of course using the unit of the active var (i.e. no conversion)
    lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)        &
      & = lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)  &
      & - (lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk) &
      &    * this_transfer_fract)

  END SUBROUTINE transfer_from_active_to_passive_var_onChunk


  ! ====================================================================================================== !
  !
  !> Copies a fraction of the content of a given active to a given passive var -> for diagnostic vars!
  !
  !
  SUBROUTINE copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
      &                         active_var_name, passive_var_name, copy_fract, conversion_fact)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations !< lcc structure of the calling lcc process
    INTEGER,                    INTENT(in)    :: i_tile
        !< index of the tile in the lcc structure for which the action should be conducted
    INTEGER,                    INTENT(in)    :: iblk    !< Number of current block (chunk)
    INTEGER,                    INTENT(in)    :: ics     !< Index of chunk start
    INTEGER,                    INTENT(in)    :: ice     !< Index of chunk end
    CHARACTER(len=*), INTENT(in)    :: active_var_name   !< (active) source var
    CHARACTER(len=*), INTENT(in)    :: passive_var_name  !< (passive) sink var
    REAL(wp),  INTENT(in), OPTIONAL :: copy_fract
        !< fraction of matter copied from active to passive var - if not given: 1.0 is assumed
    REAL(wp),  INTENT(in), OPTIONAL :: conversion_fact
        !< conversion factor if active var has different unit then passive var - if not given: 1.0 is assumed
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':copy_from_active_to_passive_var_onChunk'

    REAL(wp) :: this_copy_fract, this_conversion_fact
    INTEGER  :: i_active, i_passive
    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    i_active = one_of( active_var_name , lcc_relocations%active_vars_names)
    i_passive = one_of( passive_var_name , lcc_relocations%passive_vars_names)

    !>
    !> Assertion source variable needs to be active variable for this lcc process
    !>
    IF (i_active < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed active var ' // TRIM(active_var_name) &
        & // ' is not an active var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    !>
    !> Assertion sink variable needs to be passive variable for this lcc process
    !>
    IF (i_passive < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed passive var ' // TRIM(passive_var_name) &
        & // ' is not a passive var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    IF (PRESENT(copy_fract)) THEN
      this_copy_fract = copy_fract
    ELSE
      this_copy_fract = 1.0
    ENDIF

    IF (PRESENT(conversion_fact)) THEN
      this_conversion_fact = conversion_fact
    ELSE
      this_conversion_fact = 1.0
    ENDIF

    ! Add active content to passive var
    lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & = lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & + (lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk) * this_copy_fract * this_conversion_fact)

  END SUBROUTINE copy_from_active_to_passive_var_onChunk


  SUBROUTINE collect_matter_of_active_vars_onChunk(collected_matter, lcc_relocations, &
      &                                       i_tile, iblk, ics, ice, active_var_names)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp),                   INTENT(INOUT) :: collected_matter(:)
        !< array on which the matter shall be collected
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations
        !< lcc structure of the calling lcc process from which the matter is collected
    INTEGER,                    INTENT(IN)    :: i_tile  !< index of the tile in lcc structure
    INTEGER,                    INTENT(in)    :: iblk    !< Number of current block (chunk)
    INTEGER,                    INTENT(in)    :: ics     !< Index of chunk start
    INTEGER,                    INTENT(in)    :: ice     !< Index of chunk end
    CHARACTER(VARNAME_LEN),     INTENT(IN)    :: active_var_names(:)
        !< active vars from which the matter shall be collected
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_matter_of_active_vars_onChunk'
    INTEGER :: i, i_active
    ! -------------------------------------------------------------------------------------------------- !

    DO i = 1,SIZE(active_var_names)
      i_active = one_of( TRIM(active_var_names(i)) , lcc_relocations%active_vars_names)
      IF(i_active > 0) THEN
        collected_matter(:) = collected_matter(:) + lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)
        lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk) = 0.0_wp
      ENDIF
    ENDDO

  END SUBROUTINE collect_matter_of_active_vars_onChunk

#endif
END MODULE mo_jsb_lcc_class
