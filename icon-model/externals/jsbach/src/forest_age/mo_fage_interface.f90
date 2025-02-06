!> interface to the fage process (forest age)
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
!>#### Contains the interface to the fage process (forest age)
!>
MODULE mo_fage_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  USE mo_fage_util,          ONLY: get_mean_age
  USE mo_fage_process,       ONLY: recalc_fract_per_age_proportionally

  ! Use of processes in this module
  dsl4jsb_Use_processes FAGE_
  dsl4jsb_Use_memory(FAGE_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_fage_tasks, recalc_fract_per_age_upon_disturbance, apply_fract_per_age_change_for_alcc_forest_gain, &
    &       distribute_forest_loss_from_alcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_interface'
  CHARACTER(len=*), PARAMETER :: procname = 'fage'

  !> Type definition for fage pre task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_forest_age
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_forest_age     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_forest_age  !< Aggregates computed task variables
  END TYPE tsk_forest_age

  !> Constructor interface for fage pre task
  INTERFACE tsk_forest_age
    PROCEDURE Create_task_forest_age                        !< Constructor function for task
  END INTERFACE tsk_forest_age

CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for forest age task
  !
  FUNCTION Create_task_forest_age(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "fage"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_forest_age::return_ptr)
    CALL return_ptr%Construct(name=procname, process_id=FAGE_, owner_model_id=model_id)

  END FUNCTION Create_task_forest_age

  ! ====================================================================================================== !
  !
  !> Register tasks for fage process
  !
  SUBROUTINE Register_fage_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of fage process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !
    CALL this%Register_task(tsk_forest_age(model_id))

  END SUBROUTINE Register_fage_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "fage"
  !>
  !> Task "fage" shifts area from one age class to the next and cares for quantities conserved upon aging
  !>
  !> Note that ageing is conducted globally at the same time step 'newyear'
  !>  i.e. the update happens at different seasons, daytimes, etc
  !
  SUBROUTINE update_forest_age(tile, options)

    USE mo_util,               ONLY: int2string
    USE mo_jsb_time,           ONLY: is_newyear
    USE mo_jsb_lcc_class,      ONLY: t_jsb_lcc_proc
    USE mo_jsb_lcc,            ONLY: init_lcc_reloc, start_lcc_reloc, end_lcc_reloc

    USE mo_carbon_interface,   ONLY: merge_carbon_vars_upon_ageing
    USE mo_pheno_interface,    ONLY: merge_pheno_variables_upon_ageing
    USE mo_fage_process,       ONLY: recalc_fract_per_age_upon_ageing, synchronise_fracts

    dsl4jsb_Use_config(FAGE_)

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_forest_age'
    LOGICAL, PARAMETER :: allAtOnce = .FALSE.

    CLASS(t_jsb_tile_abstract), POINTER :: current_child, prev_child
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations

    INTEGER  :: i, iblk, ics, ice, nc
    INTEGER  :: nr_of_tiles, nacs, max_age
    REAL(wp) :: dtime
    REAL(wp), ALLOCATABLE :: lost_area(:,:), gained_area(:,:), initial_area(:,:)

    REAL(wp), DIMENSION(options%nc) :: current_fract
    REAL(wp), ALLOCATABLE :: ac_ubound(:)

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Real3D_onChunk :: fract_per_age
    dsl4jsb_Real2D_onChunk :: mean_age

    TYPE(t_jsb_vgrid), POINTER :: ac_grid
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)
    ac_grid => Get_vgrid('age_classes_grid')

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    dsl4jsb_Get_var3D_onChunk(FAGE_, fract_per_age)
    dsl4jsb_Get_var2D_onChunk(FAGE_, mean_age)

    !>
    !> 1. Get lcc structure
    !>
    dsl4jsb_Get_lcc_relocations(FAGE_, lcc_relocations)
    nr_of_tiles = lcc_relocations%nr_of_tiles
    nacs = dsl4jsb_Config(FAGE_)%nacs
    max_age = dsl4jsb_Config(FAGE_)%max_age

    ALLOCATE(ac_ubound(nacs))
    ac_ubound = ac_grid%ubounds

    ! Assert that the number of lcc tiles and of age classes are equal (might change at some point)
    IF (.NOT. nacs == nr_of_tiles) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: number of lcc tiles and of age classes ' &
        & // 'are expected to be equal. But number of lcc tiles is ' // int2string(nr_of_tiles)    &
        & // ' as opposed to the number of age classes: ' // int2string(nacs) // '. Please check!')
    END IF

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(FAGE_)) RETURN

    ! If not newyear, do nothing
    IF (.NOT. is_newyear(options%current_datetime,dtime)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert that this tile is a PFT tile
    ! - this lcc process changes the age class (pftXX_fcYY) tiles but is controlled on their parent tile
    IF (INDEX(tile%name, 'pft') == 0) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: fage processes is expected to run on pft tiles, instead' &
        & //' tried to run on '// trim(tile%name))
    END IF

    ! sync to avoid numerical drifts in the fractions on the forest pft tiles and their fract_per_age vector
    current_child => tile%Get_first_child_tile()
    DO i = 1,nacs
      CALL current_child%Get_fraction(ics, ice, iblk, fract=current_fract(:))

      CALL synchronise_fracts(nc, max_age, nacs, &
        & i, current_child%lcts(1)%lib_id, ac_grid%ubounds, current_fract(:), fract_per_age(:,:))

      CALL current_child%Set_fraction(ics, ice, iblk, fract=current_fract(:))

      current_child => current_child%Get_next_sibling_tile()
    END DO

    ! Allocate area vectors
    ALLOCATE(initial_area(nc, nr_of_tiles))
    ALLOCATE(lost_area(nc, nr_of_tiles))
    ALLOCATE(gained_area(nc, nr_of_tiles))
    initial_area(:,:) = 0.0_wp

    !>
    !> 2. loop including init, start, end, dealloc moving area only from one age class to the next
    !>

    ! get the last age class of the pft
    current_child => tile%Get_first_child_tile()

    DO i = 1,nacs-1
      prev_child => current_child
      current_child => prev_child%Get_next_sibling_tile()
    END DO
    ! Assert that currentChild is last child
    IF (.NOT. current_child%Is_last_child()) &
      & CALL finish(TRIM(routine), 'Violation of assertion: expected last child after iterating over nacs.')

    ! start with oldest age class and iterate down until second age class - the youngest age cohort cannot get area from below
    DO i = nacs,2,-1
      prev_child => current_child%Get_prev_sibling_tile()

      ! Collect initial area
      CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)

      ! init gained and lost area
      lost_area(:,:) = 0.0_wp
      gained_area(:,:) = 0.0_wp

      !----- care for shifts among age class boundaries
      ! area is moved from the younger to the older age class
      lost_area(:,i-1) = fract_per_age(:, INT(ac_ubound(i - 1)))
      gained_area(:,i) = lost_area(:,i-1)

      IF (ANY(lost_area(:, :) > 0.0_wp)) THEN
        !>     When acs are merged: also merge variables that are no cqs but used for bookkeeping
        CALL merge_pheno_variables_upon_ageing(current_child, prev_child, options, lost_area(:,i-1))
        CALL merge_carbon_vars_upon_ageing(current_child, prev_child, options, lost_area(:,i-1))
      END IF

      CALL prev_child%Get_fraction(ics, ice, iblk, fract=current_fract(:))
      current_fract(:) = current_fract(:) - lost_area(:,i-1)
      CALL prev_child%Set_fraction(ics, ice, iblk, fract=current_fract(:))

      CALL current_child%Get_fraction(ics, ice, iblk, fract=current_fract(:))
      current_fract(:) = current_fract(:) + lost_area(:,i-1)
      CALL current_child%Set_fraction(ics, ice, iblk, fract=current_fract(:))

      ! Only needs to be executed if there is area to be moved
      IF (ANY(lost_area(:, :) > 0.0_wp)) THEN
        !>
        !> 3. Each lcc process needs to call start lcc to collect to be transferred matter
        !>
        CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

        !>
        !> 4. Afterwards active content needs to be relocated to passive vars by manipulating the arrays in lcc_relocations
        !>
        !>     No need to care for active cqts upon ageing -> all passive

        !>
        !> 5. Each lcc process needs to call end lcc to make the passive matter transfer
        !>
        CALL end_lcc_reloc(lcc_relocations, options, gained_area)

      END IF ! ... only to be conducted if there is any area movement

      ! move further down in the age classes
      current_child => prev_child
    END DO ! i = nacs,2,-1 -- movements from one age-class to the text starting with older classes

    !>
    !> 6. Deallocate area vectors
    !>
    DEALLOCATE(lost_area, gained_area, initial_area, ac_ubound)

    !>
    !> 7. Further actions specific for ageing lcc process: shift fracts in fract_per_age vector
    !>
    CALL recalc_fract_per_age_upon_ageing(nc, dsl4jsb_Config(FAGE_)%max_age, fract_per_age(:,:), mean_age(:))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_forest_age

  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "fage"
  !
  SUBROUTINE aggregate_forest_age(tile, options)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk, ics, ice

    dsl4jsb_Def_memory(FAGE_)
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_forest_age'
    ! -------------------------------------------------------------------------------------------------- !

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(FAGE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_forest_age



  ! ====================================================================================================== !
  !
  !> recalculate the tracked age fractions [y] upon disturbance of given age class
  !
  SUBROUTINE recalc_fract_per_age_upon_disturbance(tile, options, i_tile, disturbed_area)

    dsl4jsb_Use_config(FAGE_)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)   :: tile           !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)      :: options        !< Additional run-time parameters
    INTEGER,  INTENT(in)                        ::  i_tile        !< index of the ac with disturbed area
    REAL(wp), INTENT(in), DIMENSION(options%nc) :: disturbed_area !< area fraction disturbed for this ac
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: ac_grid

    INTEGER :: iblk, ics, ice, nc

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Real3D_onChunk :: fract_per_age
    dsl4jsb_Real2D_onChunk :: mean_age

    CHARACTER(len=*), PARAMETER :: routine = modname//':recalc_fract_per_age_upon_disturbance'
    ! -------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    model => Get_model(tile%owner_model_id)
    ac_grid => Get_vgrid('age_classes_grid')

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    dsl4jsb_Get_var3D_onChunk(FAGE_, fract_per_age)
    dsl4jsb_Get_var2D_onChunk(FAGE_, mean_age)

    IF (ANY(disturbed_area(:) .GT. 0._wp)) THEN
      ! Put into age 0-1
      fract_per_age(:,1) = fract_per_age(:,1) + disturbed_area(:)

      ! And take from other according to defined scheme
      IF (dsl4jsb_Config(FAGE_)%disturbance_fracts_scheme == 'prop' ) THEN
        CALL recalc_fract_per_age_proportionally(nc, dsl4jsb_Config(FAGE_)%max_age, &
          & dsl4jsb_Config(FAGE_)%nacs, i_tile, disturbed_area(:), ac_grid%ubounds, fract_per_age(:,:))
      ELSE
        CALL finish(TRIM(routine), 'Currently only "prop" implemented as disturbance distribution scheme ' &
          & //'but tried '//TRIM(dsl4jsb_Config(FAGE_)%disturbance_fracts_scheme)//'. Please check!' )
      END IF
    END IF

    ! recalculate mean age
    mean_age (:) = get_mean_age(fract_per_age(:,:))

  END SUBROUTINE recalc_fract_per_age_upon_disturbance


  ! ====================================================================================================== !
  !
  !> apply the change due to alcc forest gain onto the first age class and the tract fract per age
  !
  SUBROUTINE apply_fract_per_age_change_for_alcc_forest_gain(tile, options, gain)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)   :: tile           !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)      :: options        !< Additional run-time parameters
    REAL(wp), INTENT(in), DIMENSION(:)          :: gain           !< change in cover fraction
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: ac_grid

    CLASS(t_jsb_tile_abstract), POINTER :: first_age_class_tile

    INTEGER :: iblk, ics, ice, nc

    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Real3D_onChunk :: fract_per_age
    dsl4jsb_Real2D_onChunk :: mean_age

    REAL(wp), DIMENSION(options%nc) :: current_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':apply_fract_per_age_change_for_alcc_forest_gain'
    ! -------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    model => Get_model(tile%owner_model_id)
    ac_grid => Get_vgrid('age_classes_grid')

    first_age_class_tile => tile%Get_first_child_tile()

    dsl4jsb_Get_memory(FAGE_)
    dsl4jsb_Get_var3D_onChunk(FAGE_, fract_per_age)
    dsl4jsb_Get_var2D_onChunk(FAGE_, mean_age)

    ! Update fractions in the age class tile
    CALL first_age_class_tile%Get_fraction(ics, ice, iblk, fract=current_fract)
    current_fract(:) = current_fract(:) + gain(:)
    CALL first_age_class_tile%Set_fraction(ics, ice, iblk, fract=current_fract)

    ! Since first class only covers age 0-1 it can be directly added
    fract_per_age(:, 1) = fract_per_age(:, 1) + gain(:)

    ! recalculate mean age
    mean_age(:) = get_mean_age(fract_per_age(:,:))

  END SUBROUTINE apply_fract_per_age_change_for_alcc_forest_gain


  ! ====================================================================================================== !
  !
  !> distribute changes due to alcc forest loss onto age classes and recalculate the tract fract per age
  !
  SUBROUTINE distribute_forest_loss_from_alcc(tile, options, cf_loss, ind_first_ac, lost_area)

    dsl4jsb_Use_config(FAGE_)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options      !< Additional run-time parameters
    REAL(wp), DIMENSION(:),     INTENT(in)    :: cf_loss      !< change in cover fraction
    INTEGER,                    INTENT(in)    :: ind_first_ac !< index of first ac in lost_area vector
    REAL(wp), DIMENSION(:,:),   INTENT(inout) :: lost_area    !< lost_area vector (dim2: nr of alcc tiles)
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: ac_grid

    CLASS(t_jsb_tile_abstract), POINTER :: age_class_tile

    INTEGER :: i_age_class, ind_age_class_tile, iblk, ics, ice, nc

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Real3D_onChunk :: fract_per_age
    dsl4jsb_Real2D_onChunk :: mean_age

    REAL(wp), DIMENSION(options%nc) :: current_fract, forest_pft_fract, this_cf_loss

    CHARACTER(len=*), PARAMETER :: routine = modname//':distribute_forest_loss_from_alcc'
    ! -------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    model => Get_model(tile%owner_model_id)
    ac_grid => Get_vgrid('age_classes_grid')

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    dsl4jsb_Get_var3D_onChunk(FAGE_, fract_per_age)
    dsl4jsb_Get_var2D_onChunk(FAGE_, mean_age)

    ! Get the initial area of the forest pft
    CALL tile%Get_fraction(ics, ice, iblk, fract=forest_pft_fract)
    IF (ANY(cf_loss(:) .GT. forest_pft_fract(:))) THEN
      WRITE (message_text,*) 'Violation of assertion: more area changed than available ' &
        & // 'for forest pft ', TRIM(tile%name)
      CALL finish (routine, message_text)
    END IF

    age_class_tile => tile%Get_first_child_tile()

    IF (ANY(cf_loss .GT. 0._wp)) THEN
      IF (dsl4jsb_Config(FAGE_)%alcc_age_fract_scheme == 'prop' ) THEN
        ! ---------------------------------------------------------------------------------------------- !
        ! iterate over each age class - with the prop (proportional scheme)
        ! they can be treated separately from each other with only the total
        ! area (i.e. forest_pft_fract) as weight
        i_age_class = 1
        DO WHILE(ASSOCIATED(age_class_tile))
          ! -------------------------------------------------------------------------------------------- !
          ! determine the share of the loss of this age class
          CALL age_class_tile%Get_fraction(ics, ice, iblk, fract=current_fract)
          this_cf_loss = 0._wp
          WHERE(forest_pft_fract(:) > 0._wp) !JN-TODO: replace WHERE!
            this_cf_loss(:) = cf_loss(:) * (current_fract(:) / forest_pft_fract(:))
          END WHERE
          ! -------------------------------------------------------------------------------------------- !
          ! adapt the fract_per_age tracing vector accordingly
          IF (.NOT. age_class_tile%Is_first_child()) THEN
            ! then call for this age-class:
            CALL recalc_fract_per_age_proportionally(nc, dsl4jsb_Config(FAGE_)%max_age, &
              & dsl4jsb_Config(FAGE_)%nacs, i_age_class, this_cf_loss(:), ac_grid%ubounds, fract_per_age(:,:))
          ELSE
            ! for the first age class it can directly be subtracted
            fract_per_age(:,1) = fract_per_age(:,1) - this_cf_loss(:)
          END IF
          ! -------------------------------------------------------------------------------------------- !
          ! Also adapt the cover fractions of this age class tile
          current_fract(:) = current_fract(:) - this_cf_loss(:)
          CALL age_class_tile%Set_fraction(ics, ice, iblk, fract=current_fract)
          ! -------------------------------------------------------------------------------------------- !
          ! Insert in lost_area vector
          ind_age_class_tile = ind_first_ac + i_age_class - 1
          lost_area(:,ind_age_class_tile) = lost_area(:,ind_age_class_tile) + this_cf_loss(:)
          ! -------------------------------------------------------------------------------------------- !
          ! Go to next age class (if any)
          i_age_class = i_age_class + 1
          age_class_tile => age_class_tile%Get_next_sibling_tile()
        END DO
      ELSE
        CALL finish(TRIM(routine), 'Currently only "prop" implemented as disturbance distribution scheme ' &
          & //'but tried '//TRIM(dsl4jsb_Config(FAGE_)%disturbance_fracts_scheme)//'. Please check!' )
      END IF

      ! recalculate mean age
       mean_age(:) = get_mean_age(fract_per_age(:,:))
    END IF

  END SUBROUTINE distribute_forest_loss_from_alcc

#endif
END MODULE mo_fage_interface
