!> Methods for lcc processes
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
!>#### Contains methods required to deal with actively or passively relocated conserved quantities upon lcc
!>
!> NOTE: This is currently only implemented to deal with leaves and down to two layers,
!>       if something more generic is required this needs to be re-engineered

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: finish, message
  USE mo_util,                ONLY: int2string, one_of
  USE mo_util_string,         ONLY: tolower
  USE mo_kind,                ONLY: wp

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN

  USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc, t_jsb_lcc_var_p
  USE mo_jsb_task_class,      ONLY: t_jsb_task_options
  USE mo_jsb_cqt_class,       ONLY: t_jsb_consQuan_p
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract

  USE mo_jsb_process_class,   ONLY: Get_process_id, Get_process_name

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: end_lcc_reloc, start_lcc_reloc, init_lcc_reloc, init_lcc, transfer_active_to_passive_onChunk
  PUBLIC :: count_descendants_for_up_to_2layers, collect_names_of_descendant_leaves_up_to_2layers

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lcc'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initializes/sets up the lcc instance for the given process on the given tile
  !
  SUBROUTINE init_lcc(process_name, tile, active_cqts, passive_cqts, involved_descendant_tiles)

    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid
    USE mo_jsb_grid,          ONLY: Get_grid
    USE mo_carbon_constants,  ONLY: carbon_potential_active_vars, carbon_required_passive_vars

    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*),  INTENT(in) :: process_name      !< name of the encompassing lcc process
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile  !< current tile on which the lcc process runs
    INTEGER,           INTENT(in) :: active_cqts(:)    !< cqts that are actively relocated by this lcc process
    INTEGER,           INTENT(in) :: passive_cqts(:)   !< cqts that are passively relocated by this lcc process
    CHARACTER(len=*),  INTENT(in) :: involved_descendant_tiles(:) !< descendants tiles involved in this lcc process
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_lcc'

    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_grid),  POINTER          :: hgrid
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations
    CLASS(t_jsb_tile_abstract), POINTER :: child, leave
    INTEGER :: nr_of_children, nr_of_leaves
    INTEGER :: nr_of_active_vars, nr_of_passive_vars, nr_of_tiles, nr_of_vars, nr_of_procs_with_active_vars
    INTEGER :: nr_of_active_vars_this_proc, nr_of_passive_vars_this_proc
    INTEGER :: lcc_process_id, this_process_id
    INTEGER :: ind, i_child, i_leave, i_cqt, i_var, i_active, i_active_this_tile, i_passive, i_passive_this_tile, i_tile, i_proc
    INTEGER :: nproma, nblks
    LOGICAL :: is_involved_tile
    LOGICAL, ALLOCATABLE :: is_descendant(:)
    INTEGER, ALLOCATABLE :: active_var_process_ids(:)
    CHARACTER(len=:), ALLOCATABLE :: this_proc_name
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: active_vars_this_proc(:), passive_vars_this_proc(:)
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    hgrid => Get_grid(model%grid_id)
    lcc_process_id = Get_process_id(process_name)

    ALLOCATE(t_jsb_lcc_proc::lcc_relocations)
    tile%lcc_processes(lcc_process_id)%p => lcc_relocations

    IF (debug_on()) CALL message(TRIM(routine), 'For proc '//TRIM(process_name)//' ...')

    nproma = hgrid%nproma
    nblks  = hgrid%nblks

    nr_of_tiles = SIZE(involved_descendant_tiles(:))

    lcc_relocations%nr_of_tiles = nr_of_tiles
    ALLOCATE(is_descendant(nr_of_tiles))
    is_descendant = .FALSE.

    lcc_relocations%name = TRIM(process_name) // "_lcc"
    lcc_relocations%tile_names = involved_descendant_tiles(:)

    ! for flcc a CQT can be active and passive at the same time: jsbach4: ag woody litter is burning upon fire, but then at the same
    ! time filled with new woody litter which then needs to be passively allocated, therefore no assertion on disjunct cqt groups here

    nr_of_active_vars = 0
    nr_of_passive_vars = 0
    nr_of_procs_with_active_vars = 0

    ! Iterate over the descendants of the tile to collect the CQ variables of
    ! the actively and passively relocated CQTs
    IF (tile%Has_children()) THEN
      nr_of_children = tile%Get_no_of_children()
      child => tile%Get_first_child_tile()
    ELSE
      nr_of_children = 1
      child => tile%self_tile
    END IF

    i_tile = 0
    DO i_child = 1,nr_of_children !@TODO: no do loop with counter required: do while!

      ! check if this is a leave (currently lcc can only run on leaves!) and if not go to the children of this node (i.e. the leaves)
      IF (child%Has_children()) THEN
        leave => child%Get_first_child_tile()
        nr_of_leaves = child%Get_no_of_children()
      ELSE
        leave => child
        nr_of_leaves = 1
      END IF

      DO i_leave = 1,nr_of_leaves !@TODO: no do loop with counter required: do while!

        ! only descendant tiles listed as involved are used
        is_involved_tile = .FALSE.
        ind = one_of(leave%name, involved_descendant_tiles)
        IF (ind > 0) THEN
          i_tile = i_tile + 1
          is_involved_tile = .TRUE.
          is_descendant(ind) = .TRUE.
        END IF

        ! The lcc structure is created according to active and passive cqs on the first involved tile
        IF (is_involved_tile .AND. (i_tile == 1)) THEN

          ! In the first step only collect the number of variables
          DO i_cqt = 1, leave%nr_of_cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == active_cqts(:) )) THEN
              nr_of_active_vars = nr_of_active_vars &
                & + leave%conserved_quantities(i_cqt)%p%no_of_vars
            END IF
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == passive_cqts(:) )) THEN
              nr_of_passive_vars = nr_of_passive_vars &
                & + leave%conserved_quantities(i_cqt)%p%no_of_vars
            END IF
          END DO

          lcc_relocations%nr_of_active_vars = nr_of_active_vars
          lcc_relocations%nr_of_passive_vars = nr_of_passive_vars

          ! allocate
          ALLOCATE(lcc_relocations%active_vars(nr_of_active_vars))
          ALLOCATE(lcc_relocations%active_vars_names(nr_of_active_vars))
          ALLOCATE(lcc_relocations%active_vars_cqt(nr_of_active_vars))
          ALLOCATE(lcc_relocations%active_vars_process_id(nr_of_active_vars))
          ALLOCATE(lcc_relocations%passive_vars(nr_of_passive_vars))
          ALLOCATE(lcc_relocations%passive_vars_names(nr_of_passive_vars))
          ALLOCATE(lcc_relocations%passive_vars_cqt(nr_of_passive_vars))
          ALLOCATE(lcc_relocations%passive_vars_process_id(nr_of_passive_vars))

          ALLOCATE(active_var_process_ids(nr_of_active_vars))
          active_var_process_ids = 0

          ! and then collect vars and pointers for the first involved tile
          i_active = 0
          i_passive = 0
          DO i_cqt = 1, leave%nr_of_cqts
            ! for the active cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == active_cqts(:) )) THEN
              nr_of_vars = leave%conserved_quantities(i_cqt)%p%no_of_vars
              DO i_var = 1,nr_of_vars
                i_active = i_active + 1

                this_process_id = leave%conserved_quantities(i_cqt)%p%associated_process(i_var)
                IF (.NOT. ANY(active_var_process_ids == this_process_id)) THEN
                  nr_of_procs_with_active_vars = nr_of_procs_with_active_vars + 1
                  active_var_process_ids(nr_of_procs_with_active_vars) = this_process_id
                END IF

                CALL add_cq_to_lcc_structure(nr_of_tiles, nproma, nblks, i_tile, i_cqt, i_var, i_active,                 &
                  & leave%conserved_quantities, lcc_relocations%active_vars_cqt, lcc_relocations%active_vars_process_id, &
                  & lcc_relocations%active_vars, lcc_relocations%active_vars_names)

              END DO
            END IF
            ! and separately for the passive cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == passive_cqts(:) )) THEN
              nr_of_vars = leave%conserved_quantities(i_cqt)%p%no_of_vars
              DO i_var = 1,nr_of_vars
                i_passive = i_passive + 1

                CALL add_cq_to_lcc_structure(nr_of_tiles, nproma, nblks, i_tile, i_cqt, i_var, i_passive,                   &
                  & leave%conserved_quantities, lcc_relocations%passive_vars_cqt,  lcc_relocations%passive_vars_process_id, &
                  & lcc_relocations%passive_vars, lcc_relocations%passive_vars_names )

              END DO
            ENDIF !... for active or passive CQTs

          ENDDO !... for all CQTs

          ! Collect the names of the processes carrying variables that are actively relocated by this lcc process
          lcc_relocations%nr_of_procs_with_active_vars = nr_of_procs_with_active_vars
          ALLOCATE(lcc_relocations%unique_proc_names_active_vars(nr_of_procs_with_active_vars))
          i_proc = 0
          DO ind = 1,nr_of_active_vars
            IF (active_var_process_ids(ind) /= 0) THEN
              i_proc = i_proc + 1
              lcc_relocations%unique_proc_names_active_vars(i_proc) &
                & = Get_process_name(active_var_process_ids(ind))
            ELSE
              EXIT
            END IF
          END DO
        END IF !... for the first involved tile

        ! For all other involved tiles (i.e. other than the first) only the pointers need to be collected
        ! + assert that the demanded active and passive variables are part of all involved tiles
        IF (is_involved_tile .AND. (i_tile > 1)) THEN

          i_active_this_tile = 0
          i_passive_this_tile = 0

          DO i_cqt = 1, leave%nr_of_cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == active_cqts(:) )) THEN
              CALL collect_pointers_for_cqs(i_cqt, i_tile, i_active_this_tile, leave%name, &
                & TRIM(involved_descendant_tiles(1)), leave%conserved_quantities,          &
                & lcc_relocations%active_vars, lcc_relocations%nr_of_active_vars, 'active')
            END IF
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == passive_cqts(:) )) THEN
              CALL collect_pointers_for_cqs(i_cqt, i_tile, i_passive_this_tile, leave%name, &
                & TRIM(involved_descendant_tiles(1)), leave%conserved_quantities,           &
                & lcc_relocations%passive_vars, lcc_relocations%nr_of_passive_vars, 'passive')
            END IF
          END DO

          ! Assertion: all variables found on the first tile need to also be cqs for all other tiles
          IF (.NOT. i_active_this_tile == nr_of_active_vars ) THEN
            CALL finish(TRIM(routine), &
              & 'Violation of assertion: Not all active variables found on first tile '//TRIM(involved_descendant_tiles(1)) &
              & //' were also active cqs of this tile '//TRIM(involved_descendant_tiles(i_tile))//'. Found '                &
              & //int2string(nr_of_active_vars)//' on first and '//int2string(i_active_this_tile)//' on this. Please check.')
          ENDIF
          IF (.NOT. i_passive_this_tile == nr_of_passive_vars ) THEN
            CALL finish(TRIM(routine), &
              & 'Violation of assertion: Not all passive variables found on first tile '//TRIM(involved_descendant_tiles(1)) &
              & //' were also passive cqs of this tile '//TRIM(involved_descendant_tiles(i_tile))//'. Found ' &
              & //int2string(nr_of_passive_vars)//' on first and '//int2string(i_passive_this_tile)//' on this. Please check.')
          ENDIF
        ENDIF

        leave => leave%Get_next_sibling_tile()
      END DO ! i_leave = 1,nr_of_leaves

      child => child%Get_next_sibling_tile()
    ENDDO ! i_child = 1,nr_of_children

    !>
    !> Assertion: all involved tiles are actually descendant tiles
    !>
    IF (.NOT. ALL(is_descendant)) THEN
      DO i_tile = 1,nr_of_tiles
        IF (.NOT. is_descendant(i_tile)) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: Conceptually each involved tile needs to be a leave and a descendant tile '&
            & // 'of the tile on which the lcc proc is running -- did not found ' &
            & // TRIM(involved_descendant_tiles(i_tile))//'. Please check.' )
        ENDIF
      ENDDO
    ENDIF

    DEALLOCATE(is_descendant, active_var_process_ids)

    ! JN-TODO: following code is only a check! Could e.g. be surrounded by debug flag or similar (to else disable check)
    !>
    !> Assert: check against potential active and required passive variables for each process with active vars
    !>
    !
    ! ... to enable active relocation, a process needs to 'expect' the variables assumed as active by a lcc process
    !     in order to provide a method to relocate them to passive variables. Passive variables used in an active
    !     relocation are therefore required to be passive on the calling lcc process.
    DO i_proc = 1, lcc_relocations%nr_of_procs_with_active_vars
      this_proc_name = TRIM(lcc_relocations%unique_proc_names_active_vars(i_proc))
      this_process_id = Get_process_id(this_proc_name)

      ! Count and collect active vars for this process
      nr_of_active_vars_this_proc = COUNT( lcc_relocations%active_vars_process_id == this_process_id )
      ALLOCATE(active_vars_this_proc(nr_of_active_vars_this_proc))
      ind = 0
      DO i_var = 1, lcc_relocations%nr_of_active_vars
        IF (lcc_relocations%active_vars_process_id(i_var) == this_process_id) THEN
          ind = ind + 1
          active_vars_this_proc(ind) = TRIM(lcc_relocations%active_vars(i_var)%p%name)
        ENDIF
      ENDDO

      ! and of passive vars of this proc
      nr_of_passive_vars_this_proc = COUNT(lcc_relocations%passive_vars_process_id == this_process_id)
      IF(nr_of_passive_vars_this_proc == 0) THEN
        ! Assertion: if there is an active var, at least one passive var is required, too.
        CALL finish(TRIM(routine), TRIM(lcc_relocations%name) &
          & // ' contains actively relocated variables for ' // TRIM(this_proc_name) &
          & // ' but no passively relocated sink variable for this proc - please check!')
      ENDIF
      ALLOCATE(passive_vars_this_proc(nr_of_passive_vars_this_proc))
      ind = 0
      DO i_var = 1, lcc_relocations%nr_of_passive_vars
        IF (lcc_relocations%passive_vars_process_id(i_var) == this_process_id) THEN
          ind = ind + 1
          passive_vars_this_proc(ind) = TRIM(lcc_relocations%passive_vars(i_var)%p%name)
        ENDIF
      ENDDO

      ! and check them
      SELECT CASE(tolower(this_proc_name))

      CASE ('carbon')
        CALL check_active_lcc_vars(lcc_relocations, this_proc_name, active_vars_this_proc, carbon_potential_active_vars)
        CALL check_passive_lcc_vars(lcc_relocations, this_proc_name, passive_vars_this_proc, carbon_required_passive_vars)

      CASE DEFAULT
        !>
        !> Assertion: Each process that has actively relocated conserved quantities needs
        !>            to provide lists with potential active and required passive vars
        !>
        CALL finish(TRIM(routine), this_proc_name &
          & // ' contains variables that are actively relocated by ' // TRIM(lcc_relocations%name) &
          & // ' but did not specify lists with potential active and required passive vars - please check!')
      END SELECT

      DEALLOCATE(active_vars_this_proc, passive_vars_this_proc)
    ENDDO ! i_proc

  END SUBROUTINE init_lcc

  ! ====================================================================================================== !
  !
  !> Init the relocation upon lcc process task update (currently: collect initial areas and reset relocation arrays)
  !
  SUBROUTINE init_lcc_reloc(lcc_relocations, options, tile, initial_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations  !< lcc structure of the calling lcc process
    TYPE(t_jsb_task_options),   INTENT(in) :: options             !< run-time parameters
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
        !< tile on which the lcc process is running for which the lcc structure is initialized here
    REAL(wp),                   INTENT(inout) :: initial_area(:,:)  !< current area of involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_lcc_reloc'

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile, parent
    LOGICAL :: on_second_level
    REAL(wp), DIMENSION(options%nc) :: current_fract
    INTEGER :: i_tile, i_var, iblk, ics, ice, nc, nr_of_tiles
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles

    !>
    !> Assertion: initial area size needs to match number of tiles and vector length
    !>
    IF ((SIZE(initial_area,DIM=2) /= nr_of_tiles) .OR. (SIZE(initial_area,DIM=1) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area array has unexpected size. Please check!')
    ENDIF

    ! Collect initial areas of all involved tiles
    on_second_level = .FALSE.
    IF (tile%Has_children()) THEN
      current_tile => tile%Get_first_child_tile()
    ELSE
      current_tile => tile%self_tile
    END IF

    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile) .AND. i_tile < lcc_relocations%nr_of_tiles)

      IF (current_tile%Has_children()) THEN
        on_second_level = .TRUE.
        parent => current_tile
        current_tile => parent%Get_first_child_tile()
      END IF

      IF (one_of(current_tile%name, lcc_relocations%tile_names) > 0) THEN
        i_tile = i_tile + 1

        CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract(:))

        initial_area(:, i_tile) = current_fract(:)
      ENDIF

      IF (.NOT. ASSOCIATED(current_tile%next_sibling_tile) .AND. on_second_level) THEN
        current_tile => parent
        on_second_level = .FALSE.
      END IF

      current_tile => current_tile%Get_next_sibling_tile()
    ENDDO

    ! reset relocation arrays
    DO i_var = 1,lcc_relocations%nr_of_active_vars
      lcc_relocations%active_vars(i_var)%p%relocate_this(ics:ice,:,iblk) = 0.0_wp
    ENDDO
    DO i_var = 1,lcc_relocations%nr_of_passive_vars
      lcc_relocations%passive_vars(i_var)%p%relocate_this(ics:ice,:,iblk) = 0.0_wp
    ENDDO

  END SUBROUTINE init_lcc_reloc

  ! ====================================================================================================== !
  !
  !> Start the relocation upon lcc (collect to be transferred matter from all active and passive vars)
  !
  SUBROUTINE start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for collecting to be transferred matter
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< run-time parameters
    REAL(wp),                   INTENT(in)    :: lost_area(:,:)    !< area moved from involved tiles
    REAL(wp),                   INTENT(in)    :: initial_area(:,:) !< initial area of involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':start_lcc_reloc'

    INTEGER :: iblk, nc, nr_of_tiles
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles
    !>
    !> Assertion: dimensions of lost_area and initial_area need to match
    !>
    IF ((SIZE(lost_area,DIM=2) /= nr_of_tiles) .OR. (SIZE(initial_area,DIM=2) /= nr_of_tiles) &
      & .OR. (SIZE(lost_area,DIM=1) /= nc) .OR. (SIZE(initial_area,DIM=1) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area arrays have unexpected size. Please check!')
    ENDIF
    !>
    !> Assertion: lost area is >= 0.0_wp and <= 1.0_wp
    !>
    IF (ANY(lost_area < 0.0_wp) .OR. ANY(lost_area > 1.0_wp) &
      & .OR. ANY(initial_area < 0.0_wp) .OR. ANY(initial_area > 1.0_wp)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area arrays have values out of range ([0.0, 1.0]). Please check!')
    ENDIF
    !>
    !> Assertion: lost area <= initial_area
    !>
    ! -> JN-TODO: we might want/need to discuss the concept of the e-10 areas here
    IF (ANY((initial_area - lost_area) < -EPSILON(1._wp) )) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Lost area exceeds available area ' &
        & //'for lcc process '// TRIM(lcc_relocations%name) //'! Please check!')
    ENDIF

    !>
    !> NOTE: in the current implementation it is crucial that the active vars are filled first,
    !>       because implementation of flcc made it necessary that variables can be active and passive
    !>       and currently active vars are relocated first and then passive vars are relocated
    !> HOWEVER: as soon as lcc is realised as passive -> active -> passive relocation via bare tile
    !>          it'll of course be the other way round
    !>
    !@TODO: switch order as soon as lcc is realised as passive -> active -> passive relocation via bare tile
    CALL collect_to_be_transferred_matter(options, lost_area, initial_area, &
      & nr_of_tiles, lcc_relocations%nr_of_active_vars, lcc_relocations%active_vars)

    CALL collect_to_be_transferred_matter(options, lost_area, initial_area, &
      & nr_of_tiles, lcc_relocations%nr_of_passive_vars, lcc_relocations%passive_vars)

  END SUBROUTINE start_lcc_reloc

  ! ====================================================================================================== !
  !
  !> End the relocation upon lcc (distribute all to be transferred matter - passive relocation only!)
  !>
  !> Note: the concept here is that the to be distributed area is collected of ALL shrinking tiles
  !>       i.e. if there are several shrinking tiles in the call
  !>            then there is no directed movement from tile_x -> tile_y
  !>            but sum (shrinking tiles) -> per extending tiles prop. to sum (extending tiles)
  !
  SUBROUTINE end_lcc_reloc(lcc_relocations, options, gained_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing to be relocated matter
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< run-time parameters
    REAL(wp),                   INTENT(in)    :: gained_area(:,:)  !< area moved to involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':end_lcc_reloc'

    INTEGER :: i_var, i_tile, iblk, ics, ice, nc, nr_of_tiles
    REAL(wp) :: passive_var_sum(options%nc), gained_area_sum(options%nc)
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles
    gained_area_sum = SUM(gained_area(:, :), DIM=2)

    !>
    !> Assertion: dimensions of gained_area need to be tiles x nc
    !>
    IF ((SIZE(gained_area,DIM=2) /= nr_of_tiles) .OR. (SIZE(gained_area,DIM=1) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area array has unexpected size. Please check!')
    ENDIF
    !>
    !> Assertion: reloc matter of active variables expected to be zero (i.e. relocated before routine is called)
    !>
    DO i_var = 1,lcc_relocations%nr_of_active_vars
      DO i_tile = 1,nr_of_tiles
        IF(ANY(lcc_relocations%active_vars(i_var)%p%relocate_this(ics:ice,i_tile,iblk) /= 0.0_wp)) THEN
          CALL finish(TRIM(routine), 'Violation of assertion: All active variables need to already have been ' &
            & // 'relocated upon call of this routine. However, found != 0.0 values for active var '           &
            & // TRIM(lcc_relocations%active_vars(i_var)%p%name) // ' on tile '                                &
            & // TRIM(lcc_relocations%tile_names(i_tile)) // ' for lcc process '                               &
            & // TRIM(lcc_relocations%name) // '. Please check!')
        ENDIF
      ENDDO
    ENDDO

    DO i_var = 1,lcc_relocations%nr_of_passive_vars
      ! Matter from passive variables is proportionally distributed to these passive variables on extending tiles
      passive_var_sum = SUM(lcc_relocations%passive_vars(i_var)%p%relocate_this(ics:ice,:,iblk), DIM=2)
      DO i_tile = 1,nr_of_tiles
        WHERE(gained_area(:, i_tile) > 0.0_wp)
          lcc_relocations%passive_vars(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics:ice,iblk) &
            & = lcc_relocations%passive_vars(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics:ice,iblk) &
            & + passive_var_sum(:) * (gained_area(:, i_tile) / gained_area_sum(:))
        END WHERE
      ENDDO
      lcc_relocations%passive_vars(i_var)%p%relocate_this(ics:ice,:,iblk) = 0.0_wp
    ENDDO

  END SUBROUTINE end_lcc_reloc


  ! ====================================================================================================== !
  !
  !> Attempts to relocate content of active to passive vars by calling functions on the encompassing processes
  !
  SUBROUTINE transfer_active_to_passive_onChunk(lcc_relocations, tile, i_tile, options)
    USE mo_carbon_interface,    ONLY: carbon_transfer_from_active_to_passive_vars_onChunk

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing matter from an active to a passive var
    CLASS(t_jsb_tile_abstract), INTENT(in)    :: tile    !< tile on which relocation should be conducted
    INTEGER,                    INTENT(in)    :: i_tile  !< index of the tile in lcc structure
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':transfer_active_to_passive_onChunk'

    CHARACTER(len=:), ALLOCATABLE :: this_proc_name
    INTEGER :: i_proc, i_var, proc_id
    INTEGER :: ics, ice, iblk
    LOGICAL :: has_non_zero_active_var
    ! -------------------------------------------------------------------------------------------------- !
    ics = options%ics
    ice = options%ice
    iblk = options%iblk

    IF (debug_on() .AND. options%iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    DO i_proc = 1, lcc_relocations%nr_of_procs_with_active_vars

      this_proc_name = TRIM(lcc_relocations%unique_proc_names_active_vars(i_proc))
      proc_id = Get_process_id(this_proc_name)

      !JN-TODO: required?
      ! Count number of active vars with non-zero relocations for this process
      has_non_zero_active_var = .FALSE.
      DO i_var = 1, lcc_relocations%nr_of_active_vars
        IF (lcc_relocations%active_vars_process_id(i_var) == proc_id) THEN
          IF (ANY(lcc_relocations%active_vars(i_var)%p%relocate_this(ics:ice,i_tile,iblk) .NE. 0.0_wp ) ) THEN
            has_non_zero_active_var = .TRUE.
          ENDIF
        ENDIF
      ENDDO

      IF(has_non_zero_active_var) THEN

        SELECT CASE(tolower(this_proc_name))

        CASE ('carbon')
          CALL carbon_transfer_from_active_to_passive_vars_onChunk(lcc_relocations, tile, i_tile, options)

        CASE DEFAULT
          !>
          !> Assertion: Each process that has conserved quantities used as actively relocated CQs needs
          !>            to provide a specific routine for the active relocation
          !>
          CALL finish(TRIM(routine), this_proc_name &
            & // ' contains variables that are actively relocated by ' // TRIM(lcc_relocations%name) &
            & // ' but no specific routine is defined - please check!')
        END SELECT

      ENDIF ! has_non_zero_active_var
    ENDDO ! i_proc

  END SUBROUTINE transfer_active_to_passive_onChunk

  ! ====================================================================================================== !
  !
  !> Helper: count descendants for up to two layers
  !
  SUBROUTINE count_descendants_for_up_to_2layers(tile, nr_of_descendants)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(INOUT) :: tile
    INTEGER,                    INTENT(OUT)   :: nr_of_descendants
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER   :: current_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':count_descendants_for_up_to_2layers'
    ! -------------------------------------------------------------------------------------------------- !

    ! count leaves
    nr_of_descendants = 0
    current_tile => tile%Get_first_child_tile()
    DO WHILE (ASSOCIATED(current_tile))
      IF(current_tile%Has_children()) THEN
        nr_of_descendants = nr_of_descendants + current_tile%Get_no_of_children()
      ELSE
        nr_of_descendants = nr_of_descendants + 1
      ENDIF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

  END SUBROUTINE count_descendants_for_up_to_2layers

  ! ====================================================================================================== !
  !
  !> Helper: collect names of descendant leaves up to two layers
  !
  SUBROUTINE collect_names_of_descendant_leaves_up_to_2layers(tile, descendant_leave_tiles)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(INOUT)  :: tile
    CHARACTER(len=SHORT_NAME_LEN), INTENT(OUT) :: descendant_leave_tiles(:)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: ind_list, ind_end
    CLASS(t_jsb_tile_abstract), POINTER   :: current_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_names_of_descendant_leaves_up_to_2layers'
    ! -------------------------------------------------------------------------------------------------- !

    current_tile => tile%Get_first_child_tile()
    ind_list = 1
    DO WHILE (ASSOCIATED(current_tile))
      IF (current_tile%Has_children()) THEN
        !>
        !> Only implemented for 1 or two hierarchical layers!
        !>
        IF (current_tile%first_child%Has_children()) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: number of involved layers exceeds two below given tile '&
            & //TRIM(tile%name)//'. Which is not implemented. Please check.')
        END IF

        ind_end = ind_list+current_tile%Get_no_of_children()-1
        CALL current_tile%Get_children_names(descendant_leave_tiles(ind_list:ind_end))
        ind_list = ind_end + 1
      ELSE
        ! Add name of the tile itself
        descendant_leave_tiles(ind_list) = current_tile%name
        ind_list = ind_list + 1
      END IF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

  END SUBROUTINE collect_names_of_descendant_leaves_up_to_2layers

  ! ====================================================================================================== !
  !
  !> Helper: adds cq to lcc structure upon lcc init
  !
  SUBROUTINE add_cq_to_lcc_structure(nr_of_tiles, nproma, nblks, i_tile, i_cqt, i_var, i_of_kind, &
    & conserved_quantities, cqt_of_vars, vars_process_id, var_list, var_names )

    USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_var

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                INTENT(in)    :: nr_of_tiles, nproma, nblks, i_tile, i_cqt, i_var, i_of_kind
    TYPE(t_jsb_consQuan_p), INTENT(inout) :: conserved_quantities(:)
    INTEGER,                INTENT(inout) :: cqt_of_vars(:), vars_process_id(:)
    TYPE(t_jsb_lcc_var_p),  INTENT(inout) :: var_list(:)
    CHARACTER(len=*),       INTENT(inout) :: var_names(:)
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_var), POINTER  :: a_lcc_var
    ! -------------------------------------------------------------------------------------------------- !

    cqt_of_vars(i_of_kind) = conserved_quantities(i_cqt)%p%type_id
    vars_process_id(i_of_kind) = conserved_quantities(i_cqt)%p%associated_process(i_var)

    ALLOCATE(a_lcc_var)
    var_list(i_of_kind)%p => a_lcc_var
    var_list(i_of_kind)%p%name = conserved_quantities(i_cqt)%p%cq_vars_2D(i_var)%p%name
    var_names(i_of_kind) = var_list(i_of_kind)%p%name

    ALLOCATE(var_list(i_of_kind)%p%var_on_tile(nr_of_tiles))
    var_list(i_of_kind)%p%var_on_tile(i_tile)%p => conserved_quantities(i_cqt)%p%cq_vars_2D(i_var)%p

    ALLOCATE(var_list(i_of_kind)%p%relocate_this(nproma,nr_of_tiles,nblks))
    var_list(i_of_kind)%p%relocate_this(:,:,:) = 0.0_wp

  END SUBROUTINE add_cq_to_lcc_structure


  ! ====================================================================================================== !
  !
  !> Helper: collect pointers for cqts upon lcc init
  !
  SUBROUTINE collect_pointers_for_cqs(i_cqt, i_tile, i_var_this_tile, this_tile_name, &
    & first_involved_child_name, conserved_quantities, var_list, nr_of_vars, this_case)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                INTENT(in)    :: i_cqt, i_tile
    INTEGER,                INTENT(inout) :: i_var_this_tile
    CHARACTER(len=*),       INTENT(in)    :: this_tile_name
    CHARACTER(len=*),       INTENT(in)    :: first_involved_child_name
    TYPE(t_jsb_consQuan_p), INTENT(inout) :: conserved_quantities(:)
    TYPE(t_jsb_lcc_var_p),  INTENT(inout) :: var_list(:)
    INTEGER,                INTENT(in)    :: nr_of_vars
    CHARACTER(len=*),       INTENT(in)    :: this_case
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: nr_of_vars_on_cq, i_var, i
    LOGICAL :: has_searched_var
    CHARACTER(len=VARNAME_LEN) :: thisName

    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_pointers_for_cqs'
    ! -------------------------------------------------------------------------------------------------- !

    nr_of_vars_on_cq = conserved_quantities(i_cqt)%p%no_of_vars
    DO i_var = 1,nr_of_vars_on_cq
      i_var_this_tile = i_var_this_tile + 1

      thisName = conserved_quantities(i_cqt)%p%cq_vars_2D(i_var)%p%name
      has_searched_var = .FALSE.
      DO i = 1,nr_of_vars
        IF(TRIM(var_list(i)%p%name) .EQ. TRIM(thisName)) THEN
          has_searched_var = .TRUE.
          var_list(i)%p%var_on_tile(i_tile)%p => conserved_quantities(i_cqt)%p%cq_vars_2D(i_var)%p
          EXIT
        ENDIF
      ENDDO
      ! Assertion: all passive and active variables need to exist on all tiles
      IF (.NOT. has_searched_var) THEN
        CALL finish(TRIM(routine), &
          & 'Violation of assertion: Conceptually each tile needs to contain the same active and passive vars. ' &
          & //this_case //' var ' // TRIM(thisName) // ' of tile ' // TRIM(this_tile_name) // &
          & ' does not exist on tile ' // first_involved_child_name //'. Please check.' )
      ENDIF
    ENDDO

  END SUBROUTINE collect_pointers_for_cqs


  ! ====================================================================================================== !
  !
  !> Helper: collect to be transferred matter upon start of lcc
  !
  SUBROUTINE collect_to_be_transferred_matter(options, lost_area, initial_area, &
    & nr_of_tiles, nr_of_vars, var_list)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< run-time parameters
    REAL(wp),                   INTENT(in)    :: lost_area(:,:)    !< area moved from involved tiles
    REAL(wp),                   INTENT(in)    :: initial_area(:,:) !< initial area of involved tiles
    INTEGER,                    INTENT(in)    :: nr_of_tiles, nr_of_vars
    TYPE(t_jsb_lcc_var_p),      INTENT(inout) :: var_list(:)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: ics, ice, iblk, i_var, i_tile
    REAL(wp) :: to_be_transferred(options%nc)
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    DO i_var = 1,nr_of_vars
      DO i_tile = 1,nr_of_tiles
        to_be_transferred(:) = 0.0_wp
        WHERE(lost_area(:, i_tile) > 0.0_wp)
          ! collect to be transferred matter
          to_be_transferred(:) = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics:ice,iblk) &
            &   * (lost_area(:, i_tile) / initial_area(:, i_tile))
          var_list(i_var)%p%relocate_this(ics:ice,i_tile,iblk) = to_be_transferred(:)
          ! and subtract it
          var_list(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics:ice,iblk) &
            & = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics:ice,iblk) - to_be_transferred(:)
        END WHERE
      ENDDO
    ENDDO

  END SUBROUTINE collect_to_be_transferred_matter


  ! ====================================================================================================== !
  !
  !> Assert that all vars that a lcc process assumes as active for another process are possible active vars
  !
  SUBROUTINE check_active_lcc_vars(lcc_relocations, this_proc_name, active_vars_on_this_lcc_proc, potential_active_vars)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations
        !< lcc structure of the calling lcc process
    CHARACTER(len=*),           INTENT(IN)    :: this_proc_name
        !< name of the process whose active variables are checked
    CHARACTER(len=*),           INTENT(IN)    :: active_vars_on_this_lcc_proc(:)
        !< variables active for the calling lcc process that live on the tested process
    CHARACTER(len=*),           INTENT(IN)    :: potential_active_vars(:)
        !< variables potentially active for the tested process
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_active_lcc_vars'
    INTEGER                     :: i_var
    ! -------------------------------------------------------------------------------------------------- !

    ! Assertion: calling lcc process should not call this routine with active vars
    ! that are not contained in the specified (and handled) potential active vars
    DO i_var = 1,SIZE(active_vars_on_this_lcc_proc)
      IF (.NOT. ANY(active_vars_on_this_lcc_proc(i_var) .EQ. potential_active_vars)) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: ' // TRIM(lcc_relocations%name)       &
          & // ' called with unexpected active var ' // TRIM(active_vars_on_this_lcc_proc(i_var)) &
          & // ' of process ' // this_proc_name)
      ENDIF
    ENDDO

  END SUBROUTINE check_active_lcc_vars

  ! ====================================================================================================== !
  !
  !> Assert that passive vars required for active relocation on a tested proc are passive in the calling lcc proc
  !
  SUBROUTINE check_passive_lcc_vars(lcc_relocations, this_proc_name, passive_vars_on_this_lcc_proc , required_passive_vars)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations
        !< lcc structure of the calling lcc process
    CHARACTER(len=*),           INTENT(IN)    :: this_proc_name
        !< name of the process whose active variables are checked
    CHARACTER(len=*),           INTENT(IN)    :: passive_vars_on_this_lcc_proc(:)
        !< variables passive for the calling lcc process that live on the tested process
    CHARACTER(len=*),           INTENT(IN)    :: required_passive_vars(:)
        !< passive variables required for active relocation on the tested process
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_passive_lcc_vars'
    INTEGER                     :: i_var
    ! -------------------------------------------------------------------------------------------------- !

    ! Assertion: calling lcc process needs to call this routine with all required passive vars listed
    DO i_var = 1,SIZE(required_passive_vars)
      IF (.NOT. ANY(required_passive_vars(i_var) .EQ. passive_vars_on_this_lcc_proc)) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: ' // TRIM(lcc_relocations%name) &
          & // ' did not list required passive var ' // TRIM(required_passive_vars(i_var))  &
          & // ' of process ' // this_proc_name)
      ENDIF
    ENDDO

  END SUBROUTINE check_passive_lcc_vars

#endif
END MODULE mo_jsb_lcc
