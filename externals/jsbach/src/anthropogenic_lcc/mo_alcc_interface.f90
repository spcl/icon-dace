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
!>#### Contains the interfaces to the alcc process, currently annually moving area fractions and to-be-conserved matter
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
    CALL return_ptr%Construct(name='alcc', process_id=ALCC_, owner_model_id=model_id)

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
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc, min_cf_change, min_tolerated_cf_mismatch
    USE mo_jsb_lcc,           ONLY: init_lcc_reloc, start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_)

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile, age_class_tile, forest_pft_tile
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations

    LOGICAL  :: is_age_class
    INTEGER  :: i_tile, i_ac_index, i_cf_tile, iblk, ics, ice, nc, nr_of_tiles
    REAL(wp) :: dtime
    REAL(wp), ALLOCATABLE :: lost_area(:,:), gained_area(:,:), initial_area(:,:)
    REAL(wp), DIMENSION(options%nc) :: cf_diff, current_fract

    dsl4jsb_Real3D_onChunk :: cf_current_year

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_alcc'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(ALCC_)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_current_year)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(alcc_)) RETURN

    !>
    !> Assert that PPLCC process is active
    !>
    IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to be active')
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
    dsl4jsb_Get_lcc_relocations(ALCC_, lcc_relocations)
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
    i_cf_tile = 0

    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1
      i_cf_tile = i_cf_tile + 1

      ! if running with forest age classes (acs) the forest pfts themselves (i.e. the nodes above the acs)
      ! are not part of the lcc tiles
      IF ((TRIM(model%config%usecase) == 'jsbach_forest_age_classes') .AND. (current_tile%Has_children())) THEN
          ! Assertion: age-classes are subsequent parts of the lcc structure
          age_class_tile => current_tile%Get_first_child_tile()
          i_ac_index = i_tile
          DO WHILE (ASSOCIATED(age_class_tile))
            IF (.NOT. TRIM(age_class_tile%name) .EQ. TRIM(lcc_relocations%tile_names(i_ac_index))) THEN
              CALL finish(TRIM(routine), 'Violation of assertion: age class (' //TRIM(age_class_tile%name) &
                & //') is not listed in the expected place in the alcc lcc-structure. Found ' &
                & // TRIM(lcc_relocations%tile_names(i_ac_index)) // ' instead. Please check!')
            END IF
            age_class_tile => age_class_tile%Get_next_sibling_tile()
            i_ac_index = i_ac_index + 1
          END DO

          ! forest pfts are not part of the lcc structure and thus the initial area has not yet been automatically collected
          CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract)
          cf_diff(:) = cf_current_year(:, i_cf_tile) - current_fract(:)
          ! Prevent numerical issues
          WHERE (abs(cf_diff(:)) < min_cf_change)
            cf_diff(:) = 0.0_wp
            cf_current_year(:, i_cf_tile) = current_fract(:)
          END WHERE

          ! the losses / gains of forest pfts need to be redistributed to the age classes that are its children
          CALL redistribute_cover_fraction_changes_of_forest_pft(current_tile, options, cf_diff, i_tile, gained_area, lost_area)
          i_tile = i_tile + current_tile%Get_no_of_children() - 1

          ! Finally, update the fraction on the forest pft tile itself
          CALL current_tile%Set_fraction(ics, ice, iblk, fract=cf_current_year(:, i_cf_tile))

      ELSE
        ! Assert: pft without children
        IF (current_tile%Has_children()) THEN
          CALL finish(TRIM(routine), 'Violation of assertion: alcc only expects children of a lcc tile '  &
            & //'in case of the forest age class usecase, but usecase was '// TRIM(model%config%usecase)  &
            & // '. Tile with children: ' // trim(current_tile%name) // '. Please check.')
        END IF

        ! Assert: In mo_alcc_init read_land_use_data assumes that all child tiles of the veg tile
        ! are part of the lcc structure
        IF (one_of(current_tile%name, lcc_relocations%tile_names) <= 0) THEN
          CALL finish(TRIM(routine), 'Violation of assertion: child tile of veg tile(' //TRIM(current_tile%name) &
          & //') is not part of the alcc lcc-structure. Current implementation assumes that all pfts are part of it.')
        END IF

        cf_diff(:) = cf_current_year(:, i_cf_tile) - initial_area(:, i_tile)
        ! Prevent numerical issues
        WHERE(abs(cf_diff(:)) < min_cf_change)
          cf_diff(:) = 0.0_wp
          cf_current_year(:, i_cf_tile) = initial_area(:, i_tile)
        END WHERE
        WHERE(cf_diff(:) >= 0.0_wp)
          gained_area(:, i_tile) = cf_diff(:)
        ELSEWHERE
          lost_area(:, i_tile) = -1.0_wp * cf_diff(:)
        END WHERE

        CALL current_tile%Set_fraction(ics, ice, iblk, fract=cf_current_year(:, i_cf_tile))
      ENDIF ! If with jsbach_forest_age_classes usecase and forest pft or pft without leaves

      current_tile => current_tile%Get_next_sibling_tile()
    END DO

    !>
    !>     Assert: sum of gained area should equal sum of lost area
    !>
    cf_diff = SUM(gained_area,DIM=2) - SUM(lost_area,DIM=2)
    IF (ANY( cf_diff > min_tolerated_cf_mismatch)) THEN
      WRITE (message_text,*) 'Violation of assertion: gained area does not equal lost area! Please check input data. ' &
        & // 'Max difference of: ', MAXVAL(cf_diff)
      CALL finish(TRIM(routine), TRIM(message_text))
    END IF

    !>
    !> 3. collect to be transferred matter
    !>
    CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

    !>
    !> 4. transfer matter from active to passive vars
    !>
    current_tile => tile%Get_first_child_tile()
    is_age_class = .FALSE.
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))

      ! In case of forest age classes we have to decent one level further
      IF ((TRIM(model%config%usecase) == 'jsbach_forest_age_classes') .AND. (current_tile%Has_children())) THEN
        is_age_class = .TRUE.
        forest_pft_tile => current_tile
        current_tile => current_tile%Get_first_child_tile()
      END IF

      IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN
        i_tile = i_tile + 1
        IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
          CALL transfer_active_to_passive_onChunk(lcc_relocations, current_tile, i_tile, options)
        END IF
      END IF

      current_tile => current_tile%Get_next_sibling_tile()

      ! In case of forest age classes we have to go up one level and go to the next tile there
      IF (is_age_class .AND. .NOT. ASSOCIATED(current_tile)) THEN
        current_tile => forest_pft_tile
        current_tile => current_tile%Get_next_sibling_tile()
        is_age_class = .FALSE.
      END IF

    END DO

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

  ! ====================================================================================================== !
  !
  !> Redistribution of losses and gains of a forest pft to its age-classes if running with age classes
  !
  SUBROUTINE redistribute_cover_fraction_changes_of_forest_pft(tile, options, cf_diff, i_first_child, gained_area, lost_area)
    USE mo_fage_interface,         ONLY: apply_fract_per_age_change_for_alcc_forest_gain, &
      &                                  distribute_forest_loss_from_alcc
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile          !< forest pft for which the routine is called
    TYPE(t_jsb_task_options),   INTENT(in)    :: options       !< Additional run-time parameters
    REAL(wp),                   INTENT(in)    :: cf_diff(:)    !< to be distributed changes of cover fractions
    INTEGER,                    INTENT(in)    :: i_first_child
       !< index of first child tile in gained and lost area matrices
    REAL(wp),                   INTENT(inout) :: gained_area(:,:) !< gained_area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: lost_area(:,:)   !< lost_area matrix for lcc calculations
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':redistribute_cover_fraction_changes_of_forest_pft'

    CLASS(t_jsb_tile_abstract), POINTER :: age_class_tile
    TYPE(t_jsb_model), POINTER          :: model

    REAL(wp), DIMENSION(options%nc) :: this_loss, this_gain

    INTEGER  :: iblk, ics, ice, nc
    ! -------------------------------------------------------------------------------------------------- !

    nc = options%nc
    age_class_tile => tile%Get_first_child_tile()

    this_loss(:) = -1_wp * MIN(cf_diff(:), 0._wp)
    this_gain(:) = MAX(cf_diff(:), 0._wp)

    ! per definition, all gains are distributed to the first age class
    gained_area(:, i_first_child) = gained_area(:, i_first_child) + this_gain(:)
    CALL apply_fract_per_age_change_for_alcc_forest_gain(tile, options, this_gain(:))

    ! how losses are distributed depends on configurations of the FAGE process
    CALL distribute_forest_loss_from_alcc(tile, options, this_loss(:), i_first_child, lost_area(:, :))

  END SUBROUTINE redistribute_cover_fraction_changes_of_forest_pft

#endif
END MODULE mo_alcc_interface
