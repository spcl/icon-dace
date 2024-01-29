!> interface to the alcc process (anthropogenic land cover change)
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
!>#### Contains the interfaces to the alcc process, currently annually moving land to be conserved matter 
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_alcc_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
!  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes ALCC_, PPLCC_

  ! Use of process configurations
!  dsl4jsb_Use_config(ALCC_)

  ! Use of process memories
  dsl4jsb_Use_memory(ALCC_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_alcc_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_interface'
  CHARACTER(len=*), PARAMETER :: procname = 'alcc'
  
  !> Type definition for alcc pre task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_alcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_alcc     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_alcc  !< Aggregates computed task variables
  END TYPE tsk_alcc
  
  !> Constructor interface for alcc pre task
  INTERFACE tsk_alcc
    PROCEDURE Create_task_alcc                        !< Constructor function for task
  END INTERFACE tsk_alcc
  
CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for alcc task
  !
  FUNCTION Create_task_alcc(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "alcc"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_alcc::return_ptr)
    CALL return_ptr%Construct(name=procname, process_id=ALCC_, owner_model_id=model_id)

  END FUNCTION Create_task_alcc
    
  ! ====================================================================================================== !
  !
  !> Register tasks for alcc process
  !
  SUBROUTINE Register_alcc_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of alcc process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !
    CALL this%Register_task(tsk_alcc(model_id))
    
  END SUBROUTINE Register_alcc_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "alcc"
  !>
  !> Task "alcc" currently annually moves land and to be conserved matter 
  !
  SUBROUTINE update_alcc(tile, options)

    USE mo_util,              ONLY: one_of
    USE mo_jsb_time,          ONLY: is_newday, is_newyear
    USE mo_jsb_process_class, ONLY: Get_process_id
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc
    USE mo_jsb_lcc,           ONLY: init_lcc_reloc, start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_alcc'
    LOGICAL, PARAMETER :: allAtOnce = .FALSE.

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations

    INTEGER  :: i_tile, iblk, ics, ice, nc, process_id, nr_of_tiles
    REAL(wp) :: dtime
    REAL(wp), ALLOCATABLE :: lost_area(:,:), gained_area(:,:), initial_area(:,:)
    REAL(wp), DIMENSION(options%nc) :: cf_diff

    dsl4jsb_Real3D_onChunk :: cf_current_year
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)
    process_id = Get_process_id(procname)

    dsl4jsb_Get_memory(ALCC_)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_current_year)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(alcc_)) RETURN

    !>
    !> Assert that PPLCC is active on the box tile
    !>
    CALL model%Get_top_tile(current_tile)
    IF (.NOT. current_tile%Is_process_active(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to run on the box tile')
    ENDIF

    !JN-TODO: as I understand in JSBACH3 cf changes according to luh were done on a daily timestep?!
!    ! If not newday, do nothing
!    IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN

    ! However, previous version of Rainer Schneck worked on an annual basis:
    IF( .NOT. is_newyear(options%current_datetime,dtime)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !>
    !> Assert, that this tile is the VEG tile 
    !>
    ! - this lcc process changes the PFT tiles but it runs and is controlled on their parent tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: alcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    ENDIF

    !>
    !> 1. Get lcc structure
    !>
    lcc_relocations => model%processes(process_id)%p%lcc_relocations ! JN-TODO: DSL?
    nr_of_tiles = lcc_relocations%nr_of_tiles

    ! Allocate area vectors
    ALLOCATE(initial_area(nc, nr_of_tiles))
    ALLOCATE(lost_area(nc, nr_of_tiles))
    ALLOCATE(gained_area(nc, nr_of_tiles))
    initial_area(:,:) = 0.0_wp
    lost_area(:,:) = 0.0_wp
    gained_area(:,:) = 0.0_wp

    ! Collect current area of all involved tiles
    CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)

    !>
    !> 2. Calculate area changes (according to current target cfs)
    !>
    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))

      ! Area movements are only possible among child tiles listed in the lcc structure!
      IF (one_of(current_tile%name, lcc_relocations%tile_names) > 0) THEN
        i_tile = i_tile + 1

        cf_diff(:) = cf_current_year(:, i_tile) - initial_area(:, i_tile)
        WHERE(cf_diff(:) > 0.0_wp)
          gained_area(:, i_tile) = cf_diff(:)
        ELSEWHERE
          lost_area(:, i_tile) = -1.0_wp * cf_diff(:)
        END WHERE
        CALL current_tile%Set_fraction(ics, ice, iblk, fract=cf_current_year(:, i_tile))
      ELSE
        ! In mo_alcc_init read_land_use_data assumes that all child tiles of the veg tile 
        ! are part of the lcc structure
        CALL finish(TRIM(routine), 'Violation of assertion: child tile of veg tile(' //TRIM(current_tile%name) &
          & //') is not part of the alcc lcc-structure. Current implementation assumes that all pfts are part of it.') 
      ENDIF

        current_tile => current_tile%Get_next_sibling_tile()
      END DO

      !>
      !>     Assert: sum of gained area should equal sum of lost area 
      !>
      cf_diff = SUM(gained_area,DIM=2) - SUM(lost_area,DIM=2) 
      IF (ANY( cf_diff > 1.E-11_wp)) THEN
        WRITE (message_text,*) 'Violation of assertion: gained area does not equal lost area! Please check. ' &
          & // 'Max difference of: ', MAXVAL(cf_diff)
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF
       
      !>
      !> 3. collect to be transferred matter 
      !>
      CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

      !>
      !> 4. transfer matter from active to passive vars 
      !>
      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      DO WHILE (ASSOCIATED(current_tile))
        IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN
          i_tile = i_tile + 1
          IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
            CALL transfer_active_to_passive_onChunk(lcc_relocations, current_tile, i_tile, options)
          ENDIF
        ENDIF

        current_tile => current_tile%Get_next_sibling_tile()
      ENDDO

      !>
      !> 5. make the passive matter transfer 
      !>
      CALL end_lcc_reloc(lcc_relocations, options, gained_area)

    DEALLOCATE(lost_area, gained_area, initial_area)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_alcc
  
  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "alcc"
  !
  SUBROUTINE aggregate_alcc(tile, options)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_alcc'
    ! -------------------------------------------------------------------------------------------------- !

    iblk = options%iblk
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    !> Currently nothing to do
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_alcc

#endif
END MODULE mo_alcc_interface
