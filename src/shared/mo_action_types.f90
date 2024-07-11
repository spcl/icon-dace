! Type definition for action events
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

MODULE mo_action_types

  USE mo_exception,          ONLY: finish
  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN, datetime, newDatetime, &
    &                              timedelta, newTimedelta,                     &
    &                              deallocateDatetime, deallocateTimeDelta,     &
    &                              datetimetostring, OPERATOR(>=), OPERATOR(<=), OPERATOR(+)
  USE mo_time_config,        ONLY: time_config

  IMPLICIT NONE
  PRIVATE


  ! types
  PUBLIC :: t_var_action_element
  PUBLIC :: t_var_action

  ! subroutines/functions
  PUBLIC :: new_action
  PUBLIC :: actions

  ! maximum number of actions for a single variable
  INTEGER, PARAMETER :: NMAX_ACTION=5

  ! defines single variable specific action
  TYPE t_var_action_element
    INTEGER                             :: actionTyp             ! Type of action
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: start                 ! action start time (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: end                   ! action end time   (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: ref                   ! action reference time (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: intvl                 ! action interval (duration)
                                                                 ! [PnYnMnDTnHnMn]
    TYPE(datetime)                      :: prevActive            ! actual previous triggering date
    TYPE(datetime)                      :: EventPrevTriggerDate  ! intended previous trigger date
                                                               ! Differs from prevActive in the sense
                                                               ! that it is the intended trigger date, whereas
                                                               ! prevActive is the TRUE trigger date. These two
                                                               ! can differ by the allowed 'slack'.
  END TYPE t_var_action_element

  ! defines list of variable specific actions
  !
  TYPE t_var_action
    TYPE(t_var_action_element) :: action(NMAX_ACTION)
    INTEGER                    :: n_actions = 0       ! number of variable specific actions
  CONTAINS
    ! returns number of actions of requested type
    PROCEDURE :: getNumActions   => var_action__getNumActions
    !
    ! Returns index of potentially active action-event
    PROCEDURE :: getActiveAction => var_action__getActiveAction
  END TYPE t_var_action

  ! List of available action types
  INTEGER, PARAMETER, PUBLIC :: ACTION_RESET = 1   ! re-set field to 0
  ! corresponding array of action names
  CHARACTER(LEN=10), PARAMETER, PUBLIC :: ACTION_NAMES(1) =(/"RESET     "/)

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! Creation of action events
  !------------------------------------------------------------------------------------------------
  !>
  !! Initialize single variable specific action
  !!
  !! Initialize single variable specific action. A variable named 'var_action'
  !! of type t_var_action_element is initialized.
  !!
  FUNCTION new_action(actionTyp, intvl, opt_start, opt_end, opt_ref) RESULT(var_action)
    INTEGER, INTENT(IN)                :: actionTyp ! type of action
    CHARACTER(*), INTENT(IN)           :: intvl     ! action interval
    CHARACTER(*), OPTIONAL, INTENT(IN) :: opt_start, opt_end, opt_ref ! action times [ISO_8601]
    CHARACTER(*), PARAMETER             :: routine = 'mo_action:new_action'
    TYPE(datetime), POINTER             :: dummy_ptr
    TYPE(t_var_action_element)          :: var_action
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: start0, end0, ref0
#ifdef _MTIME_DEBUG
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: start1, end1, ref1, itime
    TYPE(datetime), POINTER             :: itime_dt
#endif

    start0 = get_time_str(time_config%tc_startdate, &
      &                        time_config%tc_startdate,     opt_start)
    end0   = get_time_str(time_config%tc_stopdate, &
      &                        time_config%tc_startdate,     opt_end)
    ref0   = get_time_str(time_config%tc_exp_startdate, &
      &                        time_config%tc_startdate,     opt_ref)
#ifdef _MTIME_DEBUG
    ! CONSISTENCY CHECK:
    CALL dateTimeToString(time_config%tc_startdate, itime)
    itime_dt => newDatetime(TRIM(itime))
    start1 = get_time_str(time_config%tc_startdate,     itime_dt, opt_start)
    end1   = get_time_str(time_config%tc_stopdate,      itime_dt, opt_end)
    ref1   = get_time_str(time_config%tc_exp_startdate, itime_dt, opt_ref)
    CALL deallocateDatetime(itime_dt)

    IF ((TRIM(start0) /= TRIM(start1)) .OR.   &
      & (TRIM(end0)   /= TRIM(end1))   .OR.   &
      & (TRIM(ref0)   /= TRIM(ref1))) THEN
      CALL finish(routine, "Error in mtime consistency check!")
    END IF
#endif
    ! define var_action
    var_action%actionTyp  = actionTyp
    var_action%intvl      = TRIM(intvl)               ! interval
    var_action%start      = TRIM(start0)               ! start
    var_action%end        = TRIM(end0)                 ! end
    var_action%ref        = TRIM(ref0)                 ! ref date
    ! convert start datetime from ISO_8601 format to type datetime
    dummy_ptr => newDatetime(TRIM(start0))
    IF (.NOT. ASSOCIATED(dummy_ptr)) &
      & CALL finish(routine, "date/time conversion error: "//TRIM(start0))
    var_action%prevActive           = dummy_ptr    ! arbitrary init
    var_Action%EventPrevTriggerDate = dummy_ptr    ! arbitrary init
    CALL deallocateDatetime(dummy_ptr)
  CONTAINS

    FUNCTION get_time_str(a_time, init_time, opt_offset) RESULT(time_str)
      TYPE(datetime), POINTER, INTENT(IN) :: a_time, init_time
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: time_str
      CHARACTER(*), OPTIONAL, INTENT(IN)  :: opt_offset
      TYPE(timedelta), POINTER            :: offset_td
      TYPE(datetime), TARGET              :: time_dt

      IF (PRESENT(opt_offset)) THEN
        offset_td => newTimedelta(TRIM(opt_offset))
        time_dt = init_time + offset_td
        dummy_ptr => time_dt
        CALL dateTimeToString(dummy_ptr, time_str)
        CALL deallocateTimeDelta(offset_td)
      ELSE
        CALL dateTimeToString(a_time, time_str)
      END IF
    END FUNCTION get_time_str
  END FUNCTION new_action

  !>
  !! Generate list (array) of variable specific actions
  !!
  !! Generate list (array) of variable specific actions.
  !! Creates array 'action_list' of type t_var_action
  !!
  FUNCTION actions(a01, a02, a03, a04, a05)  RESULT(action_list)
    TYPE(t_var_action_element), INTENT(IN), OPTIONAL :: a01, a02, a03, a04, a05
    TYPE(t_var_action)             :: action_list
    INTEGER :: n_act             ! action counter
    ! create action list
    n_act = 0
    CALL add_action_item(a01)
    CALL add_action_item(a02)
    CALL add_action_item(a03)
    CALL add_action_item(a04)
    CALL add_action_item(a05)
    action_list%n_actions = n_act
  CONTAINS

    SUBROUTINE add_action_item(aX)
      TYPE(t_var_action_element), INTENT(IN), OPTIONAL :: aX

      IF (PRESENT(aX)) THEN
        n_act = n_act + 1
        action_list%action(n_act) = aX
      END IF
    END SUBROUTINE add_action_item
  END FUNCTION actions


  !>
  !! Returns the number of actions of type act_typ for a given action list
  !!
  INTEGER FUNCTION var_action__getNumActions(var_action, act_typ) RESULT(n_act)
    CLASS(t_var_action) :: var_action
    INTEGER, INTENT(IN) :: act_typ    ! type of action
    ! local
    INTEGER :: ia
 
    n_act = 0
    DO ia=1,var_action%n_actions
      IF (var_action%action(ia)%actionTyp == act_typ) THEN
        n_act = n_act + 1
      ENDIF
    ENDDO

  END FUNCTION var_action__getNumActions


  !>
  !! Get index of potentially active action-event
  !!
  !! For a specific variable,
  !! get index of potentially active action-event of selected action-type.
  !!
  !! The variable's info state and the action-type must be given.
  !! The function returns the active action index within the variable's array
  !! of actions. If no matching action is found, the function returns
  !! the result -1.
  !!
  INTEGER FUNCTION var_action__getActiveAction(var_action, actionTyp, cur_date) RESULT(actionId)
    CLASS(t_var_action)             :: var_action
    INTEGER           , INTENT(IN)  :: actionTyp    ! type of action to be searched for
    TYPE(datetime)    , INTENT(IN)  :: cur_date     ! current datetime (mtime format)
    ! local
    INTEGER :: iact                                 ! loop counter
    TYPE(datetime), POINTER :: start_date, end_date ! action-event start/end datetime

    actionId = -1
    ! loop over all variable-specific actions
    ! We unconditionally take the first active one found, even if there are more active ones.
    ! (which however would normally make little sense)
    DO iact = 1,var_action%n_actions
      IF (var_action%action(iact)%actionTyp /= actionTyp ) CYCLE  ! skip all non-matching action types
      start_date => newDatetime(var_action%action(iact)%start)
      end_date   => newDatetime(var_action%action(iact)%end)
      IF ((cur_date >= start_date) .AND. (cur_date <= end_date)) THEN
        actionId = iact   ! found active action
        CALL deallocateDatetime(start_date)
        CALL deallocateDatetime(end_date)
        EXIT      ! exit loop
      ENDIF
      CALL deallocateDatetime(start_date)
      CALL deallocateDatetime(end_date)
    ENDDO
  END FUNCTION var_action__getActiveAction

END MODULE mo_action_types

