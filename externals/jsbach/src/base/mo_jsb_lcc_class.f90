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
  USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_lcc_proc, t_jsb_lcc_var, t_jsb_lcc_var_p, &
    & transfer_from_active_to_passive_var_onChunk, collect_matter_of_active_vars_onChunk

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

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lcc_class'

CONTAINS


  ! ====================================================================================================== !
  !
  !> Relocates content of a specific given active to a specific given passive var
  !
  SUBROUTINE transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, options, &
      &                                                  active_var_name, passive_var_name)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing matter from an active to a passive var
    INTEGER,                    INTENT(in)    :: i_tile  
        !< index of the tile in the lcc structure for which the relocation should be conducted 
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< run-time parameters
    CHARACTER(len=*), INTENT(in)    :: active_var_name   !< (active) source var
    CHARACTER(len=*), INTENT(in)    :: passive_var_name  !< (passive) sink var
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':transfer_from_active_to_passive_var_onChunk'

    INTEGER :: i_active, i_passive, ics, ice, iblk
    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. options%iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

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

    ics = options%ics
    ice = options%ice
    iblk = options%iblk

    ! Add active content to passive var 
    lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & = lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & + lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)

    ! And set active to zero
    lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk) = 0.0_wp

  END SUBROUTINE transfer_from_active_to_passive_var_onChunk


  SUBROUTINE collect_matter_of_active_vars_onChunk(collected_matter, lcc_relocations, i_tile, options, active_var_names)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp),                   INTENT(INOUT) :: collected_matter(:) 
        !< array on which the matter shall be collected
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations
        !< lcc structure of the calling lcc process from which the matter is collected
    INTEGER,                    INTENT(IN)    :: i_tile  !< index of the tile in lcc structure
    TYPE(t_jsb_task_options),   INTENT(IN)    :: options !< run-time parameters
    CHARACTER(VARNAME_LEN),     INTENT(IN)    :: active_var_names(:) 
        !< active vars from which the matter shall be collected
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_matter_of_active_vars_onChunk'
    INTEGER :: i, i_active, ics, ice, iblk
    ! -------------------------------------------------------------------------------------------------- !

    ics = options%ics
    ice = options%ice
    iblk = options%iblk

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
