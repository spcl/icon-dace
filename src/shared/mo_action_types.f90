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

  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN, datetime, &
    &                              newDatetime, deallocateDatetime,&
    &                              OPERATOR(>=), OPERATOR(<=)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_var_action_element, t_var_action
  PUBLIC :: getActiveAction

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
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: lastActive            ! date of last triggering (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    TYPE(datetime)                      :: EventLastTriggerDate  ! last intended trigger date.
                                                               ! Differs from lastActive in the sense 
                                                               ! that it is the intended trigger date, whereas 
                                                               ! lastActive is the TRUE trigger date. These two 
                                                               ! can differ by the allowed 'slack'.  
  END TYPE t_var_action_element

  ! defines list of variable specific actions
  !
  TYPE t_var_action
    TYPE(t_var_action_element) :: action(NMAX_ACTION)
    INTEGER                    :: n_actions = 0       ! number of variable specific actions  
  END TYPE t_var_action

  ! List of available action types
  INTEGER, PARAMETER, PUBLIC :: ACTION_RESET = 1   ! re-set field to 0
  ! corresponding array of action names
  CHARACTER(LEN=10), PARAMETER, PUBLIC :: ACTION_NAMES(1) =(/"RESET     "/)

CONTAINS
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
  INTEGER FUNCTION getActiveAction(actions, actionTyp, cur_date) RESULT(actionId)
    TYPE(t_var_action), INTENT(IN)  :: actions      ! var metadata
    INTEGER           , INTENT(IN)  :: actionTyp    ! type of action to be searched for
    TYPE(datetime)    , INTENT(IN)  :: cur_date     ! current datetime (mtime format)
    INTEGER :: iact             ! loop counter
    TYPE(datetime), POINTER :: start_date, end_date ! action-event start/end datetime

    actionId = -1
    ! loop over all variable-specific actions
    ! We unconditionally take the first active one found, even if there are more active ones.
    ! (which however would normally make little sense)
    DO iact = 1,actions%n_actions
      IF (actions%action(iact)%actionTyp /= actionTyp ) CYCLE  ! skip all non-matching action types
      start_date => newDatetime(actions%action(iact)%start)
      end_date   => newDatetime(actions%action(iact)%end)
      IF ((cur_date >= start_date) .AND. (cur_date <= end_date)) THEN
        actionId = iact   ! found active action
        CALL deallocateDatetime(start_date)
        CALL deallocateDatetime(end_date)
        EXIT      ! exit loop
      ENDIF
      CALL deallocateDatetime(start_date)
      CALL deallocateDatetime(end_date)
    ENDDO
  END FUNCTION getActiveAction

END MODULE mo_action_types

