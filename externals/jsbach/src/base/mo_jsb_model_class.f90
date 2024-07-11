!> Contains main JSBACH model class.
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_var_class.pp-jsb.f90"

MODULE mo_jsb_model_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_kind,              ONLY: wp
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_parallel,      ONLY: Get_omp_no_of_threads, Get_omp_thread, Is_omp_inside_serial

  USE mo_jsb_control,       ONLY: debug_on
  USE mo_hsm_class,         ONLY: t_Hsm, t_State, t_Message, t_Action
  USE mo_jsb_config_class,  ONLY: t_jsb_model_config
  USE mo_jsb_process_class, ONLY: t_jsb_process, t_jsb_process_p, Get_process_name
  USE mo_jsb_task_class,    ONLY: t_jsb_process_task_p, t_Task_queue, t_jsb_task_options
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_lctlib_class,  ONLY: t_lctlib_element
  USE mo_jsb_impl_constants,ONLY: SHORT_NAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_model_m, t_jsb_model
  PUBLIC :: new_model, print_model, delete_model
  PUBLIC :: MODEL_JSBACH, MODEL_QUINCY

  !>
  !! Type for model structure
  !!
  !! This derived type defines one instance of the JSBACH model
  !!
  TYPE, EXTENDS(t_Hsm) :: t_jsb_model
    ! From base type:
    ! CHARACTER(len=30)           :: name
    ! CLASS(t_State), POINTER     :: top     => NULL()    !< Topmost state
    ! CLASS(t_State), POINTER     :: current => NULL()    !< Current state
    ! CLASS(t_State), POINTER     :: next    => NULL()    !< Next state (non_null if transition taken)
    ! CLASS(t_State), POINTER     :: source  => NULL()    !< Source state during last transition
    ! TYPE(t_Action)              :: actions(MAX_ACTIONS) !< List of possible actions (events)
    ! Extensions:
    INTEGER                            :: id                !< Model id
    CHARACTER(len=SHORT_NAME_LEN)      :: shortname         !< Model short name
    CHARACTER(len=132)                 :: description       !< Model description
    CHARACTER(len=filename_max)        :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model_config),  POINTER :: config            !< Configuration
    TYPE(t_jsb_task_options),  ALLOCATABLE :: options(:)    !< Options needed by tasks
    INTEGER                            :: grid_id
    TYPE(t_jsb_process_p),     POINTER :: processes(:)
    TYPE(t_Task_queue),        POINTER :: queue        => NULL()
    TYPE(t_jsb_process_task_p), ALLOCATABLE :: current_task(:)
    TYPE(t_lctlib_element),    POINTER :: lctlib(:)
  CONTAINS
    ! From base type:
    ! PROCEDURE :: Get_current
    ! PROCEDURE :: Get_levels_to_LCA
    ! PROCEDURE :: Next_state
    ! PROCEDURE :: Get_state
    ! PROCEDURE :: Reset_state
    ! PROCEDURE :: Register_action
    ! PROCEDURE :: Print_actions
    ! PROCEDURE :: Process_message
    ! Extensions:
    PROCEDURE :: Configure_processes => Configure_model_processes
    PROCEDURE :: Do_fractions_change => Any_process_changes_fractions
    PROCEDURE :: Get_shortname       => Get_model_shortname
    PROCEDURE :: Get_top_tile
    PROCEDURE :: Goto_next_tile
    PROCEDURE :: Get_next_tile
    PROCEDURE :: Reset
    PROCEDURE :: Reset_tiles
    PROCEDURE :: Get_tile
    PROCEDURE :: Get_process
    PROCEDURE :: Is_process_enabled
    PROCEDURE :: Queue_task
!!$    PROCEDURE :: Get_process_configs
!!$    PROCEDURE :: Process_task
    PROCEDURE :: Run_tasks
    PROCEDURE :: Set_options
    PROCEDURE :: Set_subset
    PROCEDURE :: Associate_var_pointers
#ifndef __NO_QUINCY__
    PROCEDURE :: Write_back_to_bgcms
#endif
#ifdef _OPENACC
    PROCEDURE :: gpu_to_cpu => Sync_vars_from_gpu_to_cpu
    PROCEDURE :: cpu_to_gpu => Sync_vars_from_cpu_to_gpu
#endif
  END TYPE t_jsb_model

  !>
  !! Container CLASS for JSBACH model
  !!
  !! Contains pointer to an instance of type t_jsb_model in order to allow arrays
  !! of pointers to models.
  !!
  TYPE t_jsb_model_m
    TYPE(t_jsb_model), POINTER :: m
  END TYPE t_jsb_model_m


  !>model schemes
  !!  model_scheme namlist option in mo_jsb_config_class:t_jsb_model_config
  !!  NOTE in mo_jsb_config_class:new_model_config this ENUM is not used, but replaced by INTEGER, to avoid a dependency cycle
  ENUM, BIND(C)
    ENUMERATOR :: MODEL_JSBACH=1, MODEL_QUINCY
  END ENUM


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_class'

CONTAINS

  !>
  !! Constructor method for model instance.
  !!
  !! Subroutine to construct new instance of basic model structure.
  !!
  FUNCTION new_model(id, name, shortname, description, namelist_filename) RESULT(return_ptr)

!!$    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process

    INTEGER,                    INTENT(in) :: id
    CHARACTER(len=*),           INTENT(in) :: name
    CHARACTER(len=*),           INTENT(in) :: shortname
    CHARACTER(len=*),           INTENT(in) :: description
    CHARACTER(len=*),           INTENT(in) :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model), POINTER             :: return_ptr

    INTEGER :: nthreads

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:new_model'

    CALL message(TRIM(routine), 'starting construction of new JSBACH model instance.')

    ALLOCATE(return_ptr)

    CALL return_ptr%Init_hsm(name, debug=.FALSE.)

    CALL return_ptr%Register_action(t_Action("INTEGRATE"), debug=.FALSE.)
    CALL return_ptr%Register_action(t_Action("AGGREGATE"), debug=.FALSE.)

    return_ptr%id                = id
    return_ptr%shortname         = TRIM(shortname)
    return_ptr%description       = TRIM(description)
    return_ptr%namelist_filename = TRIM(namelist_filename)
    return_ptr%config => NULL()

    nthreads = Get_omp_no_of_threads()
    ALLOCATE(return_ptr%options(nthreads))
    ALLOCATE(return_ptr%current_task(nthreads))

!!$    return_ptr%Hsm => return_ptr%t_Hsm

!!$    ALLOCATE(return_ptr%processes(max_no_of_processes))
!!$    DO iproc=1,max_no_of_processes
!!$      return_ptr%processes(iproc)%p => Create_process(iproc, return_ptr%id, TRIM(namelist_filename))
!!$    END DO

    ALLOCATE(return_ptr%queue)

    CALL message(TRIM(routine), 'construction of JSBACH model completed.')

  END FUNCTION new_model

  FUNCTION Get_process(this, iproc) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(in) :: this
    INTEGER, INTENT(in) :: iproc
    CLASS(t_jsb_process), POINTER :: return_ptr

    return_ptr => this%processes(iproc)%p

  END FUNCTION Get_process

  LOGICAL FUNCTION Is_process_enabled(this, iproc)

    CLASS(t_jsb_model), INTENT(in) :: this
    INTEGER, INTENT(in) :: iproc

    Is_process_enabled = .FALSE.
    IF (ASSOCIATED(this%processes(iproc)%p)) THEN
      IF (ASSOCIATED(this%processes(iproc)%p%config)) THEN
        Is_process_enabled = this%processes(iproc)%p%config%active
      ELSE
        ! If the process is associated but the configuration is not associated
        ! we assume here that the process is active on the model -- currently this is the case for l2a and a2l
        Is_process_enabled = .TRUE.
      ENDIF
    ENDIF

  END FUNCTION Is_process_enabled

  SUBROUTINE Queue_task(this, iproc, name)

    CLASS(t_jsb_model),        INTENT(inout) :: this
    INTEGER,                   INTENT(in)    :: iproc
    CHARACTER(LEN=*),          INTENT(in)    :: name

    CALL this%queue%Append(this%processes(iproc)%p%Get_task(TRIM(name)), Get_process_name(iproc))

  END SUBROUTINE Queue_task

  SUBROUTINE Configure_model_processes(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    INTEGER :: iproc

    DO iproc=1,SIZE(this%processes)
      IF (ASSOCIATED(this%processes(iproc)%p)) &
        & CALL this%processes(iproc)%p%Configure()
    END DO

  END SUBROUTINE Configure_model_processes

  ! ======================================================================================================= !
  !>
  !> Checks if there is an active process that changes fractions on this tile
  !>
  FUNCTION Any_process_changes_fractions(this) RESULT(l_fractions_change)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_model), INTENT(inout) :: this                 !< tile to check
    LOGICAL                           :: l_fractions_change
      !< true if there is any active land cover change (lcc) process
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: iproc
    ! -------------------------------------------------------------------------------------------------- !
    l_fractions_change = .FALSE.
    DO iproc=1,SIZE(this%processes)
      IF (ASSOCIATED(this%processes(iproc)%p)) THEN
        IF (ASSOCIATED(this%processes(iproc)%p%config)) THEN
          IF (this%processes(iproc)%p%config%active .AND. this%processes(iproc)%p%l_changes_fractions) THEN
            l_fractions_change = .TRUE.
            EXIT
          END IF
        END IF
      END IF
    END DO

  END FUNCTION Any_process_changes_fractions


  SUBROUTINE Run_tasks(this, debug)

    USE mo_jsb_control, ONLY: timer_on, timer_process_msg
    USE mo_timer,       ONLY: timer_start, timer_stop

    CLASS(t_jsb_model), INTENT(inout) :: this
    LOGICAL, OPTIONAL,  INTENT(in)    :: debug

    CLASS(t_Task_queue),  POINTER :: current_queue_element
    CLASS(t_Message), POINTER :: msg
    LOGICAL :: l_debug
    INTEGER :: no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Run_tasks'

    no_omp_thread = Get_omp_thread()

    l_debug = .FALSE.
    IF (PRESENT(debug)) l_debug = debug

    msg => t_Message('', this%Get_action('INTEGRATE'))

    ! loop over all tasks in the task queue
    current_queue_element => this%queue
    DO WHILE (ASSOCIATED(current_queue_element))

      ASSOCIATE(task => current_queue_element%task)
        IF (this%processes(task%process_id)%p%config%active) THEN  ! Only if process this task belongs to is active
          this%current_task(no_omp_thread)%p => task
          CALL this%Start() ! Start state machine
          !IF (ASSOCIATED(msg%action)) DEALLOCATE(msg%action)
          !msg%action = this%Get_action('INTEGRATE')
          IF (timer_on('all') .OR. (timer_on() .AND. debug)) CALL timer_start(timer_process_msg(this%id))
          CALL this%Process_message(msg, l_debug)
          IF (timer_on('all') .OR. (timer_on() .AND. debug)) CALL timer_stop(timer_process_msg(this%id))
        END IF
      END ASSOCIATE

      ! go to the next task
      current_queue_element => current_queue_element%next

    END DO

    DEALLOCATE(msg)

  END SUBROUTINE Run_tasks

  FUNCTION Get_model_name(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%name

  END FUNCTION Get_model_name

  FUNCTION Get_model_shortname(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%shortname

  END FUNCTION Get_model_shortname

  SUBROUTINE Get_top_tile(this, return_ptr)

    CLASS(t_jsb_model), INTENT(in) :: this
    CLASS(t_jsb_tile_abstract), POINTER, INTENT(inout) :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Get_model_top_tile'

    return_ptr => NULL()

    SELECT TYPE (top => this%top)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => top
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'tile is of type t_State, should be t_jsb_tile')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

  END SUBROUTINE Get_top_tile

  ! Set current tile to top on current OMP thread and set %visited to .FALSE.
  SUBROUTINE Reset(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Reset'

    CALL this%Reset_state()

  END SUBROUTINE Reset

  ! Reset %visited(:) to .FALSE. for all tiles and OMP threads, and current tile to top for all threads
  SUBROUTINE Reset_tiles(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    ! CLASS(t_jsb_tile_abstract), POINTER :: top

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Reset_tiles'

    CALL this%Reset_all()

  END SUBROUTINE Reset_tiles

  SUBROUTINE Goto_next_tile(this, return_ptr)

    CLASS(t_jsb_model), INTENT(inout)  :: this
    CLASS(t_jsb_tile_abstract), POINTER, INTENT(inout) :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Goto_next_tile'

    next => this%Goto_next()
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

  END SUBROUTINE Goto_next_tile

  FUNCTION Get_next_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(inout)  :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Get_next_tile'

    next => this%Get_next()
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find next tile')

  END FUNCTION Get_next_tile

  FUNCTION Get_tile(this, tile_id) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(in) :: this
    INTEGER,            INTENT(in) :: tile_id(:)
    CLASS(t_jsb_tile_abstract),  POINTER    :: return_ptr

    CLASS(t_State), POINTER :: state

    state => this%Get_state(tile_id)

    SELECT TYPE (state)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => state
    END SELECT

    NULLIFY(state)

  END FUNCTION Get_tile

  SUBROUTINE Set_options(this, current_datetime, dtime, steplen, alpha, iblk, ics, ice, nc, &
    &                    nsoil_e, nsoil_w, nsnow_e)

    USE mo_jsb_utils_iface, ONLY: assign_if_present
    USE mo_jsb_time,        ONLY: t_datetime

    CLASS(t_jsb_model),                  INTENT(inout) :: this
    TYPE(t_datetime), OPTIONAL, POINTER, INTENT(in)    :: current_datetime
    REAL(wp),         OPTIONAL,          INTENT(in)    :: dtime
    REAL(wp),         OPTIONAL,          INTENT(in)    :: steplen
    REAL(wp),         OPTIONAL,          INTENT(in)    :: alpha
    INTEGER,          OPTIONAL,          INTENT(in)    :: iblk
    INTEGER,          OPTIONAL,          INTENT(in)    :: ics
    INTEGER,          OPTIONAL,          INTENT(in)    :: ice
    INTEGER,          OPTIONAL,          INTENT(in)    :: nc
    INTEGER,          OPTIONAL,          INTENT(in)    :: nsoil_e
    INTEGER,          OPTIONAL,          INTENT(in)    :: nsoil_w
    INTEGER,          OPTIONAL,          INTENT(in)    :: nsnow_e

    INTEGER :: no_omp_thread

    no_omp_thread = Get_omp_thread()

    CALL assign_if_present(this%options(no_omp_thread)%dtime,   dtime)
    CALL assign_if_present(this%options(no_omp_thread)%steplen, steplen)
    CALL assign_if_present(this%options(no_omp_thread)%alpha,   alpha)
    CALL assign_if_present(this%options(no_omp_thread)%iblk,    iblk)
    CALL assign_if_present(this%options(no_omp_thread)%ics,     ics)
    CALL assign_if_present(this%options(no_omp_thread)%ice,     ice)
    CALL assign_if_present(this%options(no_omp_thread)%nc,      nc)
    IF (PRESENT(current_datetime)) this%options(no_omp_thread)%current_datetime => current_datetime
    ! Set number of soil and snow layers to options on all threads
    ! Note: The assign_if_present interface does not support assigning a value to a vector.
    IF (PRESENT(nsoil_e)) this%options(:)%nsoil_e = nsoil_e
    IF (PRESENT(nsoil_w)) this%options(:)%nsoil_w = nsoil_w
    IF (PRESENT(nsnow_e)) this%options(:)%nsnow_e = nsnow_e

  END SUBROUTINE Set_options

  SUBROUTINE Set_subset(this, type, iblk, ics, ice)

    USE mo_jsb_subset,     ONLY: jsbach_subsets, ON_DOMAIN, ON_CHUNK
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    CLASS(t_jsb_model), INTENT(inout) :: this
    INTEGER,            INTENT(in)    :: type
    INTEGER, OPTIONAL,  INTENT(in)    :: iblk
    INTEGER, OPTIONAL,  INTENT(in)    :: ics
    INTEGER, OPTIONAL,  INTENT(in)    :: ice

    TYPE(t_jsb_grid), POINTER :: grid

    INTEGER :: no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_subset'

    no_omp_thread = Get_omp_thread()

    SELECT CASE (type)
    CASE (ON_CHUNK)
      IF (.NOT. (PRESENT(iblk) .AND. PRESENT(ics) .AND. PRESENT(ice))) &
        & CALL finish(routine, 'Parameters missing for subset type ON_CHUNK')
      jsbach_subsets(this%id)%sub(no_omp_thread)%iblk =  iblk
      jsbach_subsets(this%id)%sub(no_omp_thread)%ics  =  ics
      jsbach_subsets(this%id)%sub(no_omp_thread)%ice  =  ice
      jsbach_subsets(this%id)%sub(no_omp_thread)%nb   =  1
      jsbach_subsets(this%id)%sub(no_omp_thread)%nc   =  ice - ics + 1
    CASE (ON_DOMAIN)
      grid => Get_grid(this%grid_id)
      jsbach_subsets(this%id)%sub(no_omp_thread)%iblk =  0
      jsbach_subsets(this%id)%sub(no_omp_thread)%ics  =  0
      jsbach_subsets(this%id)%sub(no_omp_thread)%ice  =  0
      jsbach_subsets(this%id)%sub(no_omp_thread)%nb   =  grid%nblks
      jsbach_subsets(this%id)%sub(no_omp_thread)%nc   =  grid%nproma
    END SELECT
    jsbach_subsets(this%id)%sub(no_omp_thread)%type = type

  END SUBROUTINE Set_subset

  SUBROUTINE Associate_var_pointers(this, ic_start, ic_end, iblk_start, iblk_end)

    USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d
    USE mo_jsb_subset,        ONLY: ON_DOMAIN, ON_CHUNK
    USE mo_jsb_process_class, ONLY: INHERIT_

    CLASS(t_jsb_model), INTENT(inout) :: this
    INTEGER, OPTIONAL,  INTENT(in)    :: ic_start
    INTEGER, OPTIONAL,  INTENT(in)    :: ic_end
    INTEGER, OPTIONAL,  INTENT(in)    :: iblk_start
    INTEGER, OPTIONAL,  INTENT(in)    :: iblk_end

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    INTEGER :: iproc, ivar, itype, no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':Associate_var_pointers'

    no_omp_thread = Get_omp_thread()

    IF (iblk_start == iblk_end) THEN
      itype = ON_CHUNK
    ELSE
      itype = ON_DOMAIN
    END IF

    CALL this%Reset()
    CALL this%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL this%Goto_next_tile(tile)
        CYCLE
      END IF

      DO iproc=1,SIZE(tile%mem)
        IF (.NOT. ASSOCIATED(tile%mem(iproc)%p)) CYCLE
        IF (tile%process_action(iproc) == INHERIT_) CYCLE
        DO ivar=1,tile%mem(iproc)%p%no_of_vars
          SELECT TYPE (var => tile%mem(iproc)%p%vars(ivar)%p)
          CLASS IS (t_jsb_var_real2d)
            CALL var%Associate_pointers(ic_start, ic_end, iblk_start, iblk_end)
            var%subset_type = itype
          CLASS IS (t_jsb_var_real3d)
            CALL var%Associate_pointers(ic_start, ic_end, iblk_start, iblk_end)
            var%subset_type = itype
          END SELECT
        END DO

        ! mem%pools
        IF (ASSOCIATED(tile%mem(iproc)%p%pools)) &
          & CALL tile%mem(iproc)%p%pools%Associate_var_pointers(ic_start, ic_end, iblk_start, iblk_end)
        ! mem%bgc_material
        IF (ASSOCIATED(tile%mem(iproc)%p%bgc_material)) &
          & CALL tile%mem(iproc)%p%bgc_material%Associate_var_pointers(ic_start, ic_end, iblk_start, iblk_end)

      END DO

      CALL this%Goto_next_tile(tile)

    END DO

  END SUBROUTINE Associate_var_pointers


#ifndef __NO_QUINCY__
  ! ====================================================================================================== !
  !>
  !> Write all values from matrices back to the bgcms for all tiles that have a bgcm store
  !>
  SUBROUTINE Write_back_to_bgcms(this)
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid
    USE mo_jsb_grid,          ONLY: Get_grid
    USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d
    USE mo_jsb_subset,        ONLY: ON_DOMAIN, ON_CHUNK
    USE mo_jsb_process_class, ONLY: INHERIT_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_model), INTENT(inout) :: this
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_jsb_grid),  POINTER          :: grid       ! Horizontal grid
    INTEGER                             :: startblk, endblk, iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':Store_bgcms'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    CALL this%Reset_tiles()
    CALL this%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    grid  => get_grid(this%grid_id)
    startblk = grid%Get_blk_start()
    endblk   = grid%Get_blk_end()

    DO WHILE (ASSOCIATED(tile))
      IF (ASSOCIATED(tile%bgcm_store)) THEN
        DO iblk = startblk, endblk
          ! index of gridcell at start of the chunk
          ics = grid%Get_col_start(iblk)
          ! index of gridcell at end of the chunk
          ice = grid%Get_col_end(iblk)
          CALL tile%bgcm_store%Write_stored_matrices_to_bgcms(tile%name, ics, ice, iblk)
        ENDDO
      ENDIF

      CALL this%Goto_next_tile(tile)
    END DO

  END SUBROUTINE Write_back_to_bgcms
#endif

  ! ====================================================================================================== !
  !>
  !! Destructor method for model instance.
  !!
  !! Subroutine to delete model instance from memory.
  !!
  SUBROUTINE delete_model(self)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:delete_model'

    CALL message(TRIM(routine), 'starting destruction of JSBACH model instance')

    DEALLOCATE(self)

    CALL message(TRIM(routine), 'destruction of JSBACH model instance completed.')

  END SUBROUTINE delete_model

  !>
  !! Print method for model instance.
  !!
  !! Subroutine to print description of model instance.
  !!
  SUBROUTINE print_model(self)

    USE mo_jsb_control,    ONLY: model_base_dir
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    TYPE(t_jsb_model), INTENT(in) :: self

    TYPE(t_jsb_grid), POINTER :: grid
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:print_model'

    CALL message(TRIM(routine), 'print JSBACH model description')

    WRITE(message_text,'(A,I2,A,A,A,A,A)') 'Model id: ', self%id, ', model name: ', TRIM(self%name), &
                                     ' (', TRIM(self%shortname), ')'
    CALL message('', message_text)
    CALL message('', '... '//TRIM(self%description))

    CALL message('', 'Namelist file: '//TRIM(model_base_dir)//'/'//TRIM(self%namelist_filename))

    grid => Get_grid(self%grid_id)
    IF (ASSOCIATED(grid)) CALL grid%Print

  END SUBROUTINE print_model

#ifdef _OPENACC
  SUBROUTINE Sync_vars_from_gpu_to_cpu(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_process_class,   ONLY: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_

    CLASS(t_jsb_model), INTENT(inout) :: this

    CLASS(t_jsb_memory), POINTER  :: mem
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                       :: iproc, i, no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':Sync_vars_from_gpu_to_cpu'

    IF (debug_on()) CALL message(TRIM(routine), 'Syncing variables from GPU to CPU '//TRIM(this%name))

    no_omp_thread = Get_omp_thread()
    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    CALL this%Reset_tiles()
    CALL this%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL this%Goto_next_tile(tile)
        CYCLE
      END IF

      ! Loop over all processes and initialize memory for this tile
      DO iproc=1,SIZE(tile%mem)

        ! if the current process is not on leafs of current tile exit DO loop
        SELECT CASE (tile%process_action(iproc))
        CASE (SKIP_, INHERIT_)
          CYCLE
        END SELECT

        IF (.NOT. tile%Has_process_memory(iproc)) CYCLE

        mem => tile%mem(iproc)%p
        DO i=1,mem%no_of_vars
          IF (ASSOCIATED(mem%vars(i)%p%ptr2d)) THEN
            !$ACC UPDATE HOST(mem%vars(i)%p%ptr2d) ASYNC(1) IF(mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart)
            ! IF (mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart) CALL message(routine//': gpu2cpu', TRIM(mem%vars(i)%p%full_name))
            CONTINUE
          ELSE IF (ASSOCIATED(mem%vars(i)%p%ptr3d)) THEN
            !$ACC UPDATE HOST(mem%vars(i)%p%ptr3d) ASYNC(1) IF(mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart)
            ! IF (mem%vars(i)%p%is_in_output .OR. mem%vars(i)%p%is_in_restart) CALL message(routine//': gpu2cpu', TRIM(mem%vars(i)%p%full_name))
            CONTINUE
          END IF
        END DO

      END DO

      CALL this%Goto_next_tile(tile)

    END DO

    !$ACC WAIT(1)

  END SUBROUTINE Sync_vars_from_gpu_to_cpu

  SUBROUTINE Sync_vars_from_cpu_to_gpu(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_process_class,   ONLY: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_

    CLASS(t_jsb_model), INTENT(inout) :: this

    CLASS(t_jsb_memory), POINTER  :: mem
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                       :: iproc, i, no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':Sync_vars_from_cpu_to_gpu'

    IF (debug_on()) CALL message(TRIM(routine), 'Syncing variables from CPU to GPU '//TRIM(this%name))

    no_omp_thread = Get_omp_thread()
    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    CALL this%Reset_tiles()
    CALL this%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL this%Goto_next_tile(tile)
        CYCLE
      END IF

      ! Loop over all processes and initialize memory for this tile
      DO iproc=1,SIZE(tile%mem)

        ! if the current process is not on leafs of current tile exit DO loop
        SELECT CASE (tile%process_action(iproc))
        CASE (SKIP_, INHERIT_)
          CYCLE
        END SELECT

        IF (.NOT. tile%Has_process_memory(iproc)) CYCLE

        mem => tile%mem(iproc)%p
        DO i=1,mem%no_of_vars
          ! IF (mem%vars(i)%p%is_in_restart) CALL message(routine//': cpu2gpu', TRIM(mem%vars(i)%p%full_name))
          IF (ASSOCIATED(mem%vars(i)%p%ptr2d)) THEN
            !$ACC UPDATE DEVICE(mem%vars(i)%p%ptr2d) ASYNC(1) IF(mem%vars(i)%p%is_in_restart)
            CONTINUE
          ELSE IF (ASSOCIATED(mem%vars(i)%p%ptr3d)) THEN
            !$ACC UPDATE DEVICE(mem%vars(i)%p%ptr3d) ASYNC(1) IF(mem%vars(i)%p%is_in_restart)
            CONTINUE
          END IF
        END DO

      END DO

      CALL this%Goto_next_tile(tile)

    END DO

  END SUBROUTINE Sync_vars_from_cpu_to_gpu
#endif

#endif
END MODULE mo_jsb_model_class
