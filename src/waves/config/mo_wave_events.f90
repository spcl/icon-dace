! Creation and destruction of mtime events for the wave model
!
! Creation and destruction of mtime events for the wave model
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

MODULE mo_wave_events

  USE mtime,                       ONLY: datetime, timedelta, newTimedelta, &
    &                                    deallocateTimedelta,               &
    &                                    event, eventGroup, newEvent, addEventToEventGroup
  USE mo_event_manager,            ONLY: addEventGroup, getEventGroup, printEventGroup
  USE mo_time_config,              ONLY: t_time_config

  IMPLICIT NONE

  PRIVATE

  ! subroutine
  PUBLIC :: create_wave_events

  ! events
  PUBLIC :: dummyWaveEvent
  PUBLIC :: waveEventGroup

  TYPE(event),      POINTER :: dummyWaveEvent  => NULL()
  TYPE(eventGroup), POINTER :: waveEventGroup  => NULL()

CONTAINS

  !>
  !! Create mtime events for wave model
  !!
  !! This routine creates mtime events for the wave model and
  !! puts them into suitable event groups.
  !!
  SUBROUTINE create_wave_events (time_config)

    TYPE(t_time_config), INTENT(IN) :: time_config  !< information for time control

    TYPE(datetime), POINTER         :: eventStartDate    => NULL(), &
      &                                eventEndDate      => NULL(), &
      &                                eventRefDate      => NULL()
    TYPE(timedelta), POINTER        :: eventInterval     => NULL()
    INTEGER                         :: waveEventsGroupID   !< might be changed to PUBLIC module variable,
                                                           !< if needed outside of this module
    INTEGER                         :: ierr
    LOGICAL                         :: lret


    ! create an event group for the wave model
    waveEventsGroupID = addEventGroup('waveEventGroup')
    waveEventGroup    => getEventGroup(waveEventsGroupID)


    ! create an example dummy event
    eventRefDate   => time_config%tc_exp_startdate
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate
    eventInterval  => newTimedelta("PT01h")

    dummyWaveEvent => newEvent('dummyWave', eventRefDate, eventStartDate, &
      &                           eventEndDate, eventInterval, errno=ierr)

    lret = addEventToEventGroup(dummyWaveEvent, waveEventGroup)

    CALL printEventGroup(waveEventsGroupID)
    ! cleanup
    CALL deallocateTimedelta(eventInterval)

  END SUBROUTINE create_wave_events

END MODULE mo_wave_events

