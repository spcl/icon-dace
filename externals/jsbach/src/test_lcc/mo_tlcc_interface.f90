!> interface to the tlcc process (test land cover change)
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
!>#### Contains the interfaces to the tlcc process
!>
!> Note: this test lcc process is a small test land cover change process moving just a bit of area
!>       and can be seen as example for lcc process implementations
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_tlcc_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_varlist,     ONLY: VARNAME_LEN

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes TLCC_, PPLCC_ !, TCQ_

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_tlcc_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tlcc_interface'
  CHARACTER(len=*), PARAMETER :: procname = 'tlcc'
  
  !> Type definition for tlcc pre task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_tlcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_tlcc     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_tlcc  !< Aggregates computed task variables
  END TYPE tsk_tlcc
  
  !> Constructor interface for tlcc pre task
  INTERFACE tsk_tlcc
    PROCEDURE Create_task_tlcc                        !< Constructor function for task
  END INTERFACE tsk_tlcc
  
CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for tlcc task
  !
  FUNCTION Create_task_tlcc(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "tlcc"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_tlcc::return_ptr)
    CALL return_ptr%Construct(name=procname, process_id=TLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_tlcc
    
  ! ====================================================================================================== !
  !
  !> Register tasks for tlcc process
  !
  SUBROUTINE Register_tlcc_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of tlcc process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !
    CALL this%Register_task(tsk_tlcc(model_id))
    
  END SUBROUTINE Register_tlcc_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "tlcc"
  !>
  !> Task "tlcc" simply shifts a bit of area
  !
  SUBROUTINE update_tlcc(tile, options)

    USE mo_jsb_time,          ONLY: is_newday
    USE mo_jsb_process_class, ONLY: Get_process_id
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc, transfer_from_active_to_passive_var_onChunk
    USE mo_jsb_lcc,           ONLY: init_lcc_reloc, start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_tlcc'
    LOGICAL, PARAMETER :: allAtOnce = .FALSE.

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile !, last_involved_tile
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations

    INTEGER  :: i, i_tile, iblk, ics, ice, nc, process_id, nr_of_tiles, nr_of_iterations
    REAL(wp) :: dtime
    REAL(wp), ALLOCATABLE :: lost_area(:,:), gained_area(:,:), initial_area(:,:)
    REAL(wp), DIMENSION(options%nc) :: moved_fract, current_fract

    CHARACTER(len=VARNAME_LEN) :: current_source, current_sink !,active_var_name, passive_var_name

    CHARACTER(len=VARNAME_LEN) :: source_tile(2)
    CHARACTER(len=VARNAME_LEN) :: sink_tile(2)

    LOGICAL  :: is_source, is_sink
    REAL(wp) :: shift_fract, remaining_fract
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    
    model => Get_model(tile%owner_model_id)
    process_id = Get_process_id(procname)

    !>
    !> 1. Get lcc structure
    !>
    lcc_relocations => model%processes(process_id)%p%lcc_relocations ! JN-TODO: DSL?
    nr_of_tiles = lcc_relocations%nr_of_tiles
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(TLCC_)) RETURN

    ! Assert, that PPLCC is active on the box tile
    CALL model%Get_top_tile(current_tile)
    IF (.NOT. current_tile%Is_process_active(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to run on the box tile')
    ENDIF

    ! If not newday, do nothing
    IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN
    
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert, that this tile is the VEG tile 
    ! - this lcc process works on the PFT tiles but is controlled on their parent tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: tlcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    ENDIF

    ! Allocate area vectors
    ALLOCATE(initial_area(nc, nr_of_tiles))
    ALLOCATE(lost_area(nc, nr_of_tiles))
    ALLOCATE(gained_area(nc, nr_of_tiles))
    initial_area(:,:) = 0.0_wp
    lost_area(:,:) = 0.0_wp
    gained_area(:,:) = 0.0_wp

    !>
    !> 2. Two possibilities:
    !>
    !>     2.a make all area movements at once (in that case there is no directed movement of matter
    !>         from one pft to another but only a movement of matter from all shrinking to all extending tiles)
    !>
    !>     2.b make an outer loop around init, start, end, dealloc and only move area from one pft to another
    !>         [TR: suggested to permute the order of tiles in the loop]
    !>
    source_tile = (/ 'pft04', 'pft03' /)
    sink_tile = (/ 'pft06', 'pft05' /) 

    IF (allAtOnce) THEN
      nr_of_iterations = 1
    ELSE ! not all at once but in a do loop
      nr_of_iterations = 2
    ENDIF

    DO i =1,nr_of_iterations
      IF (.NOT. allAtOnce) THEN
        current_source = source_tile(i)
        current_sink = sink_tile(i)
      ENDIF

      ! Collect inital area of all involved tiles
      CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)

      ! Move actual area and collect lost and gained area of involved tiles
      current_tile => tile%Get_first_child_tile()
      moved_fract(:) = 0.0_wp
      IF (.NOT. allAtOnce) THEN
        lost_area(:,:) = 0.0_wp
        gained_area(:,:) = 0.0_wp
      ENDIF

      i_tile = 0
      DO WHILE (ASSOCIATED(current_tile))

        ! Area movements are only possible among child tiles listed in the lcc structure!
        IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN !JN-TODO: does this work or TRIM and DO loop?
          i_tile = i_tile + 1

          is_source = .FALSE.
          is_sink = .FALSE.
          IF (allAtOnce .AND. ((current_tile%name .EQ. 'pft04') .OR. (current_tile%name .EQ. 'pft03'))) THEN
            is_source = .TRUE.
          ELSEIF (.NOT. allAtOnce .AND. (TRIM(current_tile%name) .EQ. TRIM(current_source))) THEN
            is_source = .TRUE.
          ELSEIF (allAtOnce .AND. (current_tile%name .EQ. 'pft05')) THEN
            is_sink = .TRUE.
            shift_fract = 0.3_wp
            remaining_fract = 0.7_wp 
          ELSEIF (allAtOnce .AND. (current_tile%name .EQ. 'pft06')) THEN
            is_sink = .TRUE.
            shift_fract = 1.0_wp
            remaining_fract = 0.0_wp 
          ELSEIF (.NOT. allAtOnce .AND. (TRIM(current_tile%name) .EQ. TRIM(current_sink))) THEN
            is_sink = .TRUE.
            shift_fract = 1.0_wp
            remaining_fract = 0.0_wp 
          ENDIF

          IF (is_source) THEN
            CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract(:))
            WHERE (current_fract(:) .GT. 0.00015_wp)
              lost_area(:, i_tile) = 0.00015_wp
              moved_fract(:) = moved_fract(:) + 0.00015_wp
            END WHERE 

            current_fract(:) = current_fract(:) - lost_area(:, i_tile)
            CALL current_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))

          ELSEIF(is_sink) THEN
            gained_area(:, i_tile) = moved_fract(:) * shift_fract

            CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract(:))
            current_fract(:) = current_fract(:) + gained_area(:, i_tile)
            CALL current_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))

            moved_fract(:) = moved_fract(:) * remaining_fract
          ENDIF

          !last_involved_tile => current_tile
        ENDIF

        current_tile => current_tile%Get_next_sibling_tile()
      END DO

      !>
      !> 3. Each lcc process needs to call start lcc to collect to be transferred matter 
      !>
      CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

      !>
      !> 4. Afterwards active content needs to be relocated to passive vars by manipulating the arrays in lcc_relocations
      !>
      !>     4.a This can be done by explicitly specifying how active variables should be treated
      !>         e.g. by calling transfer_from_active_to_passive_var_onChunk for named variables 
      !>
      ! 
      ! IF (last_involved_tile%Is_process_active(TCQ_)) THEN
      !   ! For the test cq a_veg is simply put in a_dead ...
      !   active_var_name = 'a_veg_c_tcq'
      !   passive_var_name = 'a_dead_c_tcq'
      !
      !   current_tile => tile%Get_first_child_tile()
      !   i_tile = 0
      !   DO WHILE (ASSOCIATED(current_tile))
      !     IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN
      !       i_tile = i_tile + 1
      !       CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, options, active_var_name, passive_var_name)
      !     ENDIF
      !   ENDIF
      !   current_tile => current_tile%Get_next_sibling_tile()
      ! ENDDO

      !>
      !>     4.b and finally by calling the more general transfer_active_to_passive_onChunk lcc routine where calls to
      !>         process specific relocation routines can be specified
      !>
      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      DO WHILE (ASSOCIATED(current_tile))
        IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN
          i_tile = i_tile + 1
          IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
            CALL transfer_active_to_passive_onChunk(lcc_relocations, current_tile, i_tile, options)
          ENDIF
          !last_involved_tile => current_tile
        ENDIF

        current_tile => current_tile%Get_next_sibling_tile()
      ENDDO

      !>
      !> 5. Each lcc process needs to call end lcc to make the passive matter transfer 
      !>
      CALL end_lcc_reloc(lcc_relocations, options, gained_area)
    ENDDO ! Number of iterations depending on if allAtOnce  

    !>
    !> 6. Deallocate area vectors
    !>
    DEALLOCATE(lost_area, gained_area, initial_area)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_tlcc
  
  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "tlcc"
  !
  SUBROUTINE aggregate_tlcc(tile, options)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_tlcc'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk = options%iblk
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    !> Currently nothing to do
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_tlcc

#endif
END MODULE mo_tlcc_interface
