! A small timer specialization geared towards easy use as a throw-away local instrumentation timer.
!
! Usage
!      Simply add a save variable to a function/subroutine like this:
!          type(t_localTimer), save :: myTimer = localTimerDefault
!
!      Then you can simply start and stop the timer using
!          call startTimer(myTimer, "description for the timer output table")
!          ! do stuff
!          call stopTimer(myTimer)
!
!      Note, however, that you must call startTimer() either with all processes or with none due to the way mo_real_timer builds the timer output table.
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

MODULE mo_local_timer
    USE mo_real_timer, ONLY: new_timer, timer_start, timer_stop
    
    IMPLICIT NONE

PUBLIC

    TYPE t_localTimer
        LOGICAL :: initialized
        INTEGER :: timerId
    END TYPE
    TYPE(t_localTimer), PARAMETER :: localTimerDefault = t_localTimer(.false., 0)

CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------
    !> start the timer
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE startTimer(me, description)
        TYPE(t_localTimer), INTENT(INOUT) :: me
        CHARACTER(len=*), INTENT(IN) :: description

        IF(.not. me%initialized) THEN
            me%timerId = new_timer(description)
            me%initialized = .true.
        ENDIF
        CALL timer_start(me%timerId)
    END SUBROUTINE startTimer

    !-------------------------------------------------------------------------------------------------------------------------------
    !> stop the timer
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE stopTimer(me)
        TYPE(t_localTimer), INTENT(INOUT) :: me

        CALL timer_stop(me%timerId)
    END SUBROUTINE stopTimer

END MODULE
