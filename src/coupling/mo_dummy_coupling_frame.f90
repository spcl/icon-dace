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

! @brief Routines for the dummy initialisation of YAC
!
! The purpose of routines construct_dummy_coupling and destruct_dummy_coupling is
! to initialise YAC on asynchronous processes that do not take part in any
! coupling exchanges. The component definition and enddef operations are collective
! operations in the MPI sense. Thus all MPI processes must take part in the
! component initialisation and enddef operation.

MODULE mo_dummy_coupling_frame

  USE mo_coupling_config, ONLY: is_coupled_run
  USE mo_exception,       ONLY: message, message_text
  USE mo_timer,           ONLY: ltimer, timer_start, timer_stop, &
    &                           timer_coupling_init
  USE mo_coupling_utils,  ONLY: cpl_def_main_dummy, cpl_enddef

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_dummy_coupling_frame'

  PUBLIC :: construct_dummy_coupling
  PUBLIC :: destruct_dummy_coupling
 
CONTAINS

  !>
  !! Initializes the coupling
  SUBROUTINE construct_dummy_coupling ( comp_name )

    CHARACTER(LEN=*), INTENT(IN) :: comp_name
    CHARACTER(*), PARAMETER :: &
      routine = str_module // ":construct_dummy_coupling"

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start(timer_coupling_init)

    WRITE(message_text,*) &
      "YAC dummy initialisation for processes of type ",  TRIM(comp_name)
    CALL message(routine, message_text)

    ! Inform YAC about what we are
    CALL cpl_def_main_dummy(routine, TRIM(comp_name))

    ! All processes have to participate in the enddef operation
    CALL cpl_enddef(routine)

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_dummy_coupling

  !>
  !! Finalizes the coupling
  SUBROUTINE destruct_dummy_coupling ( comp_name )

    CHARACTER(LEN=*), INTENT(IN) :: comp_name
    CHARACTER(*), PARAMETER :: &
      routine = str_module // ":destruct_dummy_coupling"

    IF ( is_coupled_run() ) THEN

      WRITE(message_text,*) &
        "YAC termination of processes of type ", TRIM(comp_name)
      CALL message(routine, message_text)

    ENDIF

  END SUBROUTINE destruct_dummy_coupling

END MODULE mo_dummy_coupling_frame
