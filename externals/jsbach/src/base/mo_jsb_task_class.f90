!> Contains process task class
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_task_class
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_exception,   ONLY: message, finish
  USE mo_timer,       ONLY: new_timer, timer_start, timer_stop

  USE mo_jsb_control, ONLY: timer_aggregate, timer_integrate, timer_on, debug_on

  ! USE mo_hsm_class,       ONLY: t_jsb_task_msg => t_Message
  USE mo_hsm_class,       ONLY: t_Message
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile_abstract
  USE mo_jsb_time,        ONLY: t_datetime

  IMPLICIT NONE
  PRIVATE

  ! PUBLIC :: t_jsb_process_task, t_jsb_process_task_p, t_Task_queue, t_jsb_task_msg, t_jsb_task_options
  PUBLIC :: t_jsb_process_task, t_jsb_process_task_p, t_Task_queue, t_jsb_task_options

  TYPE t_jsb_task_options
    REAL(wp)                  :: dtime   = -1._wp
    REAL(wp)                  :: steplen = -1._wp
    TYPE(t_datetime), POINTER :: current_datetime
    INTEGER                   :: iblk    = -1  ! Number of current block (chunk)
    INTEGER                   :: ics     = -1  ! Index of chunk start
    INTEGER                   :: ice     = -1  ! Index of chunk end
    INTEGER                   :: nc      = -1  ! Length of chunk
  END TYPE t_jsb_task_options

!!$  TYPE, EXTENDS(t_Message) :: t_jsb_task_msg
!!$    ! From base type:
!!$    ! CHARACTER(len=50) :: name
!!$    ! TYPE(t_Action)    :: action
!!$  END TYPE t_jsb_task_msg

  TYPE, ABSTRACT :: t_jsb_process_task
    CHARACTER(len=50)              :: name
    INTEGER                        :: process_id     = 0
    INTEGER                        :: owner_model_id = 0
    INTEGER                        :: timer_integrate = 0
    INTEGER                        :: timer_aggregate = 0
    !CLASS(t_jsb_task_msg), POINTER :: msg => NULL()
!!$    PROCEDURE(Integrate_task), POINTER, NOPASS :: Integrate
!!$    PROCEDURE(Aggregate_task), POINTER, NOPASS :: Aggregate
  CONTAINS
    PROCEDURE(Interface_task), DEFERRED, NOPASS :: Integrate
    PROCEDURE(Interface_task), DEFERRED, NOPASS :: Aggregate
    PROCEDURE                                   :: Do_it
    PROCEDURE                                   :: Construct => Construct_task
  END TYPE t_jsb_process_task

  ABSTRACT INTERFACE
    SUBROUTINE Interface_task(tile, options)
      IMPORT :: t_jsb_tile_abstract, t_jsb_task_options
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_task_options),   INTENT(in) :: options
    END SUBROUTINE Interface_task
  END INTERFACE

  TYPE t_jsb_process_task_p
    CLASS(t_jsb_process_task), POINTER :: p
  END TYPE t_jsb_process_task_p

  TYPE t_Task_queue
    CLASS(t_jsb_process_task), POINTER :: task    => NULL()
    TYPE(t_Task_queue),        POINTER :: next    => NULL()
  CONTAINS
    PROCEDURE :: Append => Append_task
  END TYPE t_Task_queue

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_task_class'

CONTAINS

  SUBROUTINE Construct_task(this, name, process_id, owner_model_id)

    CLASS(t_jsb_process_task), INTENT(inout) :: this
    CHARACTER(len=*),          INTENT(in)    :: name
    INTEGER,                   INTENT(in)    :: process_id
    INTEGER,                   INTENT(in)    :: owner_model_id

    ! ALLOCATE(t_jsb_process_task::task)
    this%name = name
    this%process_id = process_id
    this%owner_model_id = owner_model_id

    IF (timer_on('all')) THEN
      this%timer_integrate = new_timer('jsb:update_'//TRIM(name))
      this%timer_aggregate = new_timer('jsb:aggregate_'//TRIM(name))
    END IF

  END SUBROUTINE Construct_task

  SUBROUTINE Do_it(this, msg, tile, options)

    CLASS(t_jsb_process_task),  INTENT(in)    :: this
    ! CLASS(t_jsb_task_msg),      INTENT(in)    :: msg
    CLASS(t_Message),           INTENT(in)    :: msg
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':Do_it'

    ! If process the task belongs to is not active on this tile, do nothing.
    IF (.NOT. tile%Is_process_active(this%process_id)) RETURN

    SELECT TYPE (this)
    CLASS IS (t_jsb_process_task)

      IF (msg%action%name == '') &
        & CALL finish(TRIM(routine), 'Action not set for task '//TRIM(this%name))

!!$      IF (.NOT. ASSOCIATED(task%Integrate)) &
!!$        & CALL finish(TRIM(routine), 'Procedure pointer for Integrate not set')
!!$
!!$      IF (.NOT. ASSOCIATED(task%Aggregate)) &
!!$        & CALL finish(TRIM(routine), 'Procedure pointer for Aggregate not set')
!!$
      SELECT CASE(TRIM(msg%action%name))
      CASE ('INTEGRATE')
        IF (timer_on('detail')) THEN
          CALL timer_start(timer_integrate(this%owner_model_id))
          IF (this%timer_integrate > 0) CALL timer_start(this%timer_integrate)
        END IF
        CALL this%Integrate(tile, options)
        IF (timer_on('detail')) THEN
          IF (this%timer_integrate > 0) CALL timer_stop(this%timer_integrate)
          CALL timer_stop(timer_integrate(this%owner_model_id))
        END IF
      CASE ('AGGREGATE')
        IF (timer_on('detail')) THEN
          CALL timer_start(timer_aggregate(this%owner_model_id))
          IF (this%timer_aggregate > 0) CALL timer_start(this%timer_aggregate)
        END IF
        CALL this%Aggregate(tile, options)
        IF (timer_on('detail')) THEN
          IF (this%timer_aggregate > 0) CALL timer_stop(this%timer_aggregate)
          CALL timer_stop(timer_aggregate(this%owner_model_id))
        END IF
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unknown action for task'//TRIM(this%name))
      END SELECT

    END SELECT

  END SUBROUTINE Do_it

  SUBROUTINE Append_task(this, task, process_name)

    CLASS(t_Task_queue), TARGET, INTENT(inout)  :: this
    CLASS(t_jsb_process_task), TARGET           :: task
    CHARACTER(len=*),            INTENT(in)     :: process_name

    CLASS(t_Task_queue), POINTER :: current, next

    CHARACTER(len=*), PARAMETER :: routine = modname//':Append_task'

    IF (.NOT. ASSOCIATED(this%task)) THEN
      ! If queue is empty just add task as first element
      this%task => task
    ELSE
      ! Else find last task in queue
      current => this
      DO WHILE (ASSOCIATED(current%next))
        current => current%next
      END DO

      ! Append new element
      ALLOCATE(next)
      next%task => task
      current%next => next
    END IF

    IF (debug_on()) CALL message(TRIM(routine), &
      & 'Appended task "'//TRIM(process_name)//':'//TRIM(task%name)//'" to model queue')

  END SUBROUTINE Append_task

#endif
END MODULE mo_jsb_task_class
