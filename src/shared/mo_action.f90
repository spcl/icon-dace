! Routines for defining, initializing and executing action events
!
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
!
! *****************************************************************
!            RECIPE FOR CREATING A NEW ACTION EVENT
! *****************************************************************
! 1. Define a new actionTyp such as ACTION_RESET in mo_action_types
! 2. Assign new actionTyp to variables of your choice.
!    E.g. add the following code snippet to an add_var/add_ref of your choice:
!    action_list=actions(new_action(ACTION_XXX,'PTXXH'), new_action(...), ...)
!    ACTION_XXX is the actionTyp defined in step 1, and PTXXH is the
!    interval at which the action should be triggered.
! 3. Create an extension of the abstract type t_action_obj and overwrite
!    the deferred procedure 'kernel' with your action-specific kernel-routine
!    (to be defined in step 5).
! 4. Create a variable (object) of the type defined in step 3.
! 5. Write your own action-Routine (action-kernel). This is the routine which actually does
!    the work. (see e.g. routine 'reset_kernel' for actionTyp=ACTION_RESET)
! 6. Initialize the new action object by invoking the type-bound procedure 'initialize'.
!    (CALL act_obj%initialize(actionTyp)). The actiontyp defines the specific action to be
!    initialized. By this, you assign all matching fields to your particular action.
!    I.e. this is the reverse operation of assigning actions to fields as done in step 2.
! 7. Execute your newly defined action object at a suitable place by invoking the
!    type-bound procedure 'execute' (CALL act_obj%execute(slack)). 'Slack' is the user-defined
!    maximum allowed time mismatch for executing the action.

MODULE mo_action

  USE mo_kind,               ONLY: wp, i8
  USE mo_mpi,                ONLY: my_process_is_stdio, p_pe, p_io, p_comm_work, p_bcast
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: vname_len
  USE mtime,                 ONLY: event, newEvent, datetime, isCurrentEventActive,  &
    &                              MAX_DATETIME_STR_LEN, datetimetostring,           &
    &                              MAX_EVENTNAME_STR_LEN, timedelta, newTimedelta,   &
    &                              deallocateTimedelta, getPTStringFromMS,           &
    &                              getTriggeredPreviousEventAtDateTime,              &
    &                              getPTStringFromMS, datetimetostring, OPERATOR(+)
  USE mo_util_mtime,         ONLY: is_event_active, getElapsedSimTimeInSeconds
  USE mo_util_string,        ONLY: remove_duplicates
  USE mo_util_table,         ONLY: initialize_table, finalize_table, add_table_column, &
    &                              set_table_entry, print_table, t_table
  USE mo_action_types,       ONLY: t_var_action, action_names, ACTION_RESET, new_action, actions
  USE mo_grid_config,        ONLY: n_dom
  USE mo_run_config,         ONLY: msg_level
  USE mo_parallel_config,    ONLY: proc0_offloading
  USE mo_var_list_register,  ONLY: t_vl_register_iter
  USE mo_var,                ONLY: t_var
  USE mo_var_list,           ONLY: t_var_list_ptr, find_list_element
  USE mo_fortran_tools,      ONLY: init, assign_if_present
  USE mo_timer,              ONLY: ltimer, timers_level, timer_start, timer_stop, timer_action

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_action"

  PUBLIC :: reset_act
  PUBLIC :: get_prev_trigger_time

  ! FIXME
  ! Removal of this line requires a change in the ART external
  ! USE mo_action,        ONLY: action_reset, new_action, actions
  ! -- >
  ! USE mo_action_types,  ONLY: action_reset, new_action, actions
  !
  PUBLIC :: action_reset, new_action, actions

  INTEGER, PARAMETER :: NMAX_VARS = 220 ! maximum number of fields that can be
                                        ! assigned to a single action (should be
                                        ! sufficient for 4 domains in one nested run)

  ! type for generating an array of pointers of type t_var_list_element
  TYPE t_var_element_ptr
    TYPE(t_var), POINTER :: p
    TYPE(event), POINTER :: mevent     ! event from mtime library
    INTEGER              :: patch_id  ! patch on which field lives
  END TYPE t_var_element_ptr

  ! base type for action objects
  TYPE, ABSTRACT :: t_action_obj
    INTEGER                    :: actionTyp                   ! Type of action
    TYPE(t_var_element_ptr)    :: var_element_ptr(NMAX_VARS)  ! assigned variables
    INTEGER                    :: var_action_index(NMAX_VARS) ! index in var_element_ptr(10)%action
    INTEGER                    :: nvars                 ! number of variables for which
                                                        ! this action is to be performed
  CONTAINS
    PROCEDURE :: initialize => action__collect_vars  ! initialize action object
#if defined (__SX__) || defined (__NEC_VH__)
    PROCEDURE :: execute    => action__execute_SX    ! execute action object
#else
    PROCEDURE :: execute    => action__execute       ! execute action object
#endif
    PROCEDURE :: print_setup=> action__print_setup   ! Screen print out of action object setup
    ! deferred routine for action specific kernel (to be defined in extended type)
    PROCEDURE(kernel), deferred :: kernel
  END TYPE t_action_obj
  !
  ! kernel interface
  ABSTRACT INTERFACE
    SUBROUTINE kernel(act_obj, ivar)
      IMPORT                    :: t_action_obj
      CLASS(t_action_obj)       :: act_obj
      INTEGER, INTENT(IN)       :: ivar
    END SUBROUTINE kernel
  END INTERFACE

  ! extension of the action base type for the purpose of creating objects of that type.
  ! create specific type for reset-action
  TYPE, EXTENDS(t_action_obj) :: t_reset_obj
  CONTAINS
    PROCEDURE :: kernel => reset_kernel     ! type-specific action kernel (to be defined by user)
  END TYPE t_reset_obj

  ! create action object
  TYPE(t_reset_obj) :: reset_act  ! action which resets field to resetval%rval

CONTAINS

  !>
  !! Initialize action object
  !!
  !! Assign variables to specific actions/initialize action-object.
  !! When generating a new field via add_var, it is possible to assign various
  !! actions to this field. Here, we go the other way around. For a specific
  !! action, we loop over all fields and check, whether this action must
  !! be performed for the particular field. If this is the case, this field
  !! is assigned to the action.
  !! The action to be initialized is identified via the actionTyp.
  !!
  !! Loop over all variables and collect the variables names
  !! corresponding to the action @p action%actionTyp
  !!
  SUBROUTINE action__collect_vars(act_obj, actionTyp)
    CLASS(t_action_obj) :: act_obj
    INTEGER, INTENT(IN) :: actionTyp
    INTEGER :: iact, iv, nv, slen, tlen, vlen
    TYPE(t_var), POINTER :: elem
    TYPE(t_var_action), POINTER :: action_list
    TYPE(t_vl_register_iter) :: vl_iter
    CHARACTER(LEN=MAX_EVENTNAME_STR_LEN):: event_name
    CHARACTER(LEN=vname_len) :: varlist(NMAX_VARS)
    CHARACTER(*), PARAMETER :: sep = ', '
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_action:action_collect_vars")

    nv = 0
    act_obj%actionTyp = actionTyp
    ! loop over all variable lists and variables
    DO WHILE(vl_iter%next())
      DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)%p
        ! point to variable specific action list
        action_list => elem%info%action_list
        ! Loop over all variable-specific actions
        DO iact = 1, action_list%n_actions
          ! If the variable specific action fits, assign variable
          ! to the corresponding action.
          IF (action_list%action(iact)%actionTyp == actionTyp) THEN
            ! Add field to action object
            nv = nv + 1
            IF (nv > NMAX_VARS) THEN
              CALL finish(routine, "Max. number of vars for action"// &
                   &       TRIM(ACTION_NAMES(act_obj%actionTyp))// &
                   &       " exceeded! Increase NMAX_VARS in mo_action.f90!")
            END IF
            act_obj%var_element_ptr(nv)%p => elem
            act_obj%var_action_index(nv) = iact
            act_obj%var_element_ptr(nv)%patch_id = vl_iter%cur%p%patch_id
            ! Create event for this specific field
            WRITE(event_name,'(a,i2,2a)') 'act_TYP', actionTyp, &
                 '_', TRIM(action_list%action(iact)%intvl)
            act_obj%var_element_ptr(nv)%mevent => newEvent(event_name, &
              & action_list%action(iact)%ref, action_list%action(iact)%start, &
              & action_list%action(iact)%end, action_list%action(iact)%intvl)
          END IF
        ENDDO ! loop over variable-specific actions
      ENDDO ! loop over vlist "i"
    ENDDO ! i = 1, SIZE(var_lists)
    ! set nvars
    act_obj%nvars = nv
    IF (msg_level >= 11) THEN
      ! remove duplicate variable names
      DO iv = 1, nv
        varlist(iv) = act_obj%var_element_ptr(iv)%p%info%name
      ENDDO
      CALL remove_duplicates(varlist, nv)
      WRITE(message_text,'(a)') 'Variables assigned to action '//TRIM(ACTION_NAMES(act_obj%actionTyp))//':'
      tlen = LEN_TRIM(message_text)
      DO iv = 1, nv
        slen = MERGE(1, 2, iv .EQ. 1)
        vlen = LEN_TRIM(varlist(iv))

        IF (tlen+slen+vlen > LEN(message_text)) THEN
          CALL message('', message_text)
          message_text = '... '
          tlen = 4
        END IF

        WRITE(message_text(tlen+1:tlen+slen+vlen),'(2a)') sep(3-slen:2), varlist(iv)(1:vlen)
        tlen = tlen + slen + vlen
      ENDDO
      CALL message('',message_text)
      IF(my_process_is_stdio()) CALL act_obj%print_setup()
    ENDIF
  END SUBROUTINE action__collect_vars

  !>
  !! Screen print out of action event setup.
  SUBROUTINE action__print_setup (act_obj)
    CLASS(t_action_obj)  :: act_obj  !< action for which setup will be printed
    TYPE(t_table)   :: table
    INTEGER         :: ivar, irow, jg, var_action_idx
    CHARACTER(LEN=2):: str_patch_id

    ! table-based output
    CALL initialize_table(table)
    ! the latter is no longer mandatory
    CALL add_table_column(table, "VarName")
    CALL add_table_column(table, "PID")
    CALL add_table_column(table, "Ref date")
    CALL add_table_column(table, "Start date")
    CALL add_table_column(table, "End date")
    CALL add_table_column(table, "Interval")
    irow = 0
    ! print event info sorted by patch ID in ascending order
    DO jg = 1, n_dom
      DO ivar=1,act_obj%nvars
        IF (act_obj%var_element_ptr(ivar)%patch_id /= jg) CYCLE
        var_action_idx = act_obj%var_action_index(ivar)
        irow = irow + 1
        CALL set_table_entry(table,irow,"VarName", TRIM(act_obj%var_element_ptr(ivar)%p%info%name))
        write(str_patch_id,'(i2)')  act_obj%var_element_ptr(ivar)%patch_id
        CALL set_table_entry(table,irow,"PID", TRIM(str_patch_id))
        CALL set_table_entry(table,irow,"Ref date", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%ref))
        CALL set_table_entry(table,irow,"Start date", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%start))
        CALL set_table_entry(table,irow,"End date", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%end))
        CALL set_table_entry(table,irow,"Interval", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%intvl))
      ENDDO
    ENDDO  ! jg
    CALL print_table(table, opt_delimiter=' | ')
    CALL finalize_table(table)
    WRITE (0,*) " " ! newline
  END SUBROUTINE action__print_setup

  !>
  !! Execute action
  !!
  !! For each variable attached to this action it is checked, whether the action should
  !! be executed at the datetime given. This routine does not make any assumption
  !! about the details of the action to be executed. The action itself is encapsulated
  !! in the kernel-routine.
  !!
  SUBROUTINE action__execute(act_obj, slack, mtime_date)
    CLASS(t_action_obj)       :: act_obj
    REAL(wp), INTENT(IN)      :: slack     !< allowed slack for event triggering  [s]
    INTEGER :: iv, ia
    TYPE(t_var), POINTER :: field
    TYPE(event), POINTER :: mevent
    TYPE(datetime), POINTER :: mtime_date
    LOGICAL :: isactive
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: cur_dtime_str, slack_str
    TYPE(timedelta), POINTER :: p_slack                    ! slack in 'timedelta'-Format
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_action:action_execute")

  !-------------------------------------------------------------------------
    IF (ltimer .AND. timers_level > 4) CALL timer_start(timer_action)

    ! compute allowed slack in PT-Format
    ! Use factor 999 instead of 1000, since no open interval is available
    ! needed [trigger_date, trigger_date + slack[
    ! used   [trigger_date, trigger_date + slack]
    CALL getPTStringFromMS(INT(999.0_wp*slack,i8), slack_str)
    ! get slack in 'timedelta'-format appropriate for isCurrentEventActive
    p_slack => newTimedelta(slack_str)
    ! Loop over all fields attached to this action
    DO iv = 1, act_obj%nvars
      ia = act_obj%var_action_index(iv)
      field => act_obj%var_element_ptr(iv)%p
      mevent => act_obj%var_element_ptr(iv)%mevent
      ! Check whether event-pointer is associated.
      IF (.NOT.ASSOCIATED(mevent)) THEN
        WRITE (message_text,'(a,i2,a)') 'WARNING: action event ', ia, &
          & ' of field ' // TRIM(field%info%name) // ': no event-ptr associated!'
        CALL message(routine,message_text)
      ENDIF
      ! Note that a second call to isCurrentEventActive will lead to
      ! a different result! Is this a bug or a feature?
      ! triggers in interval [trigger_date + slack]
      isactive = isCurrentEventActive(mevent, mtime_date, plus_slack=p_slack)
      ! Check whether the action should be triggered for variable
      ! under consideration
      IF (isactive) THEN
        ! store latest true triggering date
        field%info%action_list%action(ia)%prevActive = mtime_date
        ! store latest intended triggering date
        CALL getTriggeredPreviousEventAtDateTime(mevent, &
          field%info%action_list%action(ia)%EventPrevTriggerDate)
        IF (msg_level >= 12) THEN
          CALL datetimeToString(mtime_date, cur_dtime_str)
          WRITE(message_text,'(5a,i2,a,a)') 'action ',TRIM(ACTION_NAMES(act_obj%actionTyp)),&
            &  ' triggered for ', TRIM(field%info%name),' (PID ',               &
            &  act_obj%var_element_ptr(iv)%patch_id, ') at ', TRIM(cur_dtime_str)
          CALL message(TRIM(routine),message_text)
        ENDIF
        ! perform action on element number 'ivar'
        CALL act_obj%kernel(iv)
      ENDIF
    ENDDO
    ! cleanup
    CALL deallocateTimedelta(p_slack)

    IF (ltimer .AND. timers_level > 4) CALL timer_stop(timer_action)
  END SUBROUTINE action__execute

  !! Execute action
  !! !!! 'optimized' variant for SX-architecture                                   !!!
  !! !!! isCurrentEventActive is only executed by PE0 and broadcasted to the other !!!
  !! !!! workers on the vector engine.                                             !!!
  !!
  !! For each variable attached to this action it is checked, whether the action should
  !! be executed at the datetime given. This routine does not make any assumption
  !! about the details of the action to be executed. The action itself is encapsulated
  !! in the kernel-routine.
  !!
  SUBROUTINE action__execute_SX(act_obj, slack, mtime_date)
    CLASS(t_action_obj)       :: act_obj
    REAL(wp), INTENT(IN)      :: slack     !< allowed slack for event triggering  [s]
    ! local variables
    INTEGER :: iv, ia
    TYPE(t_var), POINTER :: field
    TYPE(event), POINTER :: mevent
    TYPE(datetime), POINTER :: mtime_date
    LOGICAL :: isactive(NMAX_VARS)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: cur_dtime_str, slack_str
    TYPE(timedelta), POINTER :: p_slack                    ! slack in 'timedelta'-Format
    CHARACTER(*), PARAMETER :: routine = "mo_action:action_execute_SX"

  !-------------------------------------------------------------------------
    IF (ltimer .AND. timers_level > 4) CALL timer_start(timer_action)

    isactive(:) =.FALSE.
    ! compute allowed slack in PT-Format
    ! Use factor 999 instead of 1000, since no open interval is available
    ! needed [trigger_date, trigger_date + slack[
    ! used   [trigger_date, trigger_date + slack]
    IF (.NOT. proc0_offloading .OR. p_pe == p_io) THEN
      CALL getPTStringFromMS(INT(999.0_wp*slack,i8),slack_str)
      ! get slack in 'timedelta'-format appropriate for isCurrentEventActive
      p_slack => newTimedelta(slack_str)
      ! Check for all action events whether they are due at the datetime given.
      DO iv = 1, act_obj%nvars
        ia = act_obj%var_action_index(iv)
        mevent => act_obj%var_element_ptr(iv)%mevent
        ! Check whether event-pointer is associated.
        IF (.NOT.ASSOCIATED(mevent)) THEN
          field => act_obj%var_element_ptr(iv)%p
          WRITE (message_text,'(a,i2,a,a,a)')                           &
               'WARNING: action event ', ia, ' of field ',  &
               TRIM(field%info%name),': Event-Ptr is disassociated!'
          CALL message(routine,message_text)
        ENDIF
        ! Note that a second call to isCurrentEventActive will lead to
        ! a different result! Is this a bug or a feature?
        ! triggers in interval [trigger_date + slack]
        isactive(iv) = is_event_active(mevent, mtime_date, proc0_offloading, plus_slack=p_slack, opt_lasync=.TRUE.)
      ENDDO
      ! cleanup
      CALL deallocateTimedelta(p_slack)
    ENDIF
    ! broadcast
    IF (proc0_offloading) CALL p_bcast(isactive, p_io, p_comm_work)
    ! Loop over all fields attached to this action
    DO iv = 1, act_obj%nvars
      ! Check whether the action should be triggered for var under consideration
      IF (isactive(iv)) THEN
        ia = act_obj%var_action_index(iv)
        field => act_obj%var_element_ptr(iv)%p
        mevent => act_obj%var_element_ptr(iv)%mevent
        ! Check whether event-pointer is associated.
        IF (.NOT. ASSOCIATED(mevent)) THEN
          WRITE (message_text,'(a,i2,a,a,a)')                           &
               'WARNING: action event ', ia, ' of field ',  &
                TRIM(field%info%name),': Event-Ptr is disassociated!'
          CALL message(routine,message_text)
        ENDIF
        ! store latest true triggering date
        field%info%action_list%action(ia)%prevActive = mtime_date
        ! store latest intended triggering date
        CALL getTriggeredPreviousEventAtDateTime(mevent, &
          field%info%action_list%action(ia)%EventPrevTriggerDate)
        IF (msg_level >= 12) THEN
          CALL datetimeToString(mtime_date, cur_dtime_str)
          WRITE(message_text,'(5a,i2,a,a)') 'action ',TRIM(ACTION_NAMES(act_obj%actionTyp)),&
            &  ' triggered for ', TRIM(field%info%name),' (PID ',               &
            &  act_obj%var_element_ptr(iv)%patch_id,') at ', TRIM(cur_dtime_str)
          CALL message(routine, message_text)
        ENDIF
        ! perform action on element number 'ivar'
        CALL act_obj%kernel(iv)
      ENDIF
    ENDDO

    IF (ltimer .AND. timers_level > 4) CALL timer_stop(timer_action)
  END SUBROUTINE action__execute_SX

  !>
  !! Reset-action kernel
  !!
  SUBROUTINE reset_kernel(act_obj, ivar)
    CLASS (t_reset_obj)  :: act_obj
    INTEGER, INTENT(IN) :: ivar    ! element number
    CHARACTER(*), PARAMETER :: routine = modname//'::reset_kernel'

    ! re-set field to its pre-defined reset-value
    IF (ASSOCIATED(act_obj%var_element_ptr(ivar)%p%r_ptr)) THEN
!$OMP PARALLEL
      CALL init(act_obj%var_element_ptr(ivar)%p%r_ptr, &
           act_obj%var_element_ptr(ivar)%p%info%resetval%rval, lacc=.TRUE.)
!$OMP END PARALLEL
    ELSE IF (ASSOCIATED(act_obj%var_element_ptr(ivar)%p%i_ptr)) THEN
!$OMP PARALLEL
      CALL init(act_obj%var_element_ptr(ivar)%p%i_ptr, &
           act_obj%var_element_ptr(ivar)%p%info%resetval%ival, lacc=.TRUE.)
!$OMP END PARALLEL
    ELSE
      CALL finish (routine, 'Field not allocated for '//TRIM(act_obj%var_element_ptr(ivar)%p%info%name))
    ENDIF
  END SUBROUTINE reset_kernel


  !>
  !! Returns the forecast lead time of the last previous execution
  !! of a specific action
  !!
  !! Limitations:
  !! Application of this routine to variables which have multiple action events
  !! of the same type can lead to erroneous results, as the routine is not
  !! able to distinguish them. Therefore the routine aborts in this case.
  !!
  SUBROUTINE get_prev_trigger_time(var_list, var_name, res_time, opt_act_typ)

    TYPE(t_var_list_ptr), INTENT(IN) :: var_list(:)
    CHARACTER(*),         INTENT(IN) :: var_name
    REAL(wp),             INTENT(OUT):: res_time(:)
    INTEGER, OPTIONAL,    INTENT(IN) :: opt_act_typ   ! type of action

    ! Local
    TYPE(t_var), POINTER        :: var      ! pointer to specific var_list element
    INTEGER                     :: jg, ia   ! loop variables
    INTEGER                     :: act_typ  ! type of action
    INTEGER                     :: nact_typ ! number of actions of desired type
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//'::get_prev_trigger_time'

    IF (.NOT. proc0_offloading .OR. my_process_is_stdio()) THEN
      act_typ = ACTION_RESET
      CALL assign_if_present(act_typ, opt_act_typ)

      dom_loop: DO jg = 1, n_dom
        var => find_list_element (var_list(jg), TRIM(var_name))

        ASSOCIATE(action_list => var%info%action_list)
          ! get number of available actions of requested type
          nact_typ = action_list%getNumActions(act_typ)

          ! abort, if the variable has multiple actions of requested type attached to it
          IF (nact_typ > 1) THEN
            WRITE(message_text,'(a,a,a,i2)')                              &
              &  'Subroutine not applicable to variable', TRIM(var_name), &
              &  'since it has more than one action event of type ', act_typ
            CALL finish(routine, message_text)
          ENDIF

          ! abort, if the variable has no action of requested type attached to it
          IF (nact_typ == 0) THEN
            WRITE(message_text,'(a,a,a,i2)')                              &
              &  'Subroutine not applicable to variable', TRIM(var_name), &
              &  'since it has no action event of type ', act_typ
            CALL finish(routine, message_text)
          ENDIF

          act_loop: DO ia=1,action_list%n_actions
            IF (action_list%action(ia)%actionTyp == ACTION_RESET) THEN
              res_time(jg) = getElapsedSimTimeInSeconds       &
              &            (                                  &
              &             action_list%action(ia)%prevActive &
              &            )
              EXIT act_loop
            ENDIF
          END DO act_loop
        END ASSOCIATE
      ENDDO dom_loop
    ENDIF

    IF (proc0_offloading) CALL p_bcast(res_time(1:n_dom), p_io, p_comm_work)

  END SUBROUTINE get_prev_trigger_time


END MODULE mo_action

