!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
MODULE mo_exception

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message

CONTAINS

  SUBROUTINE message(leading_text, message_text)
    CHARACTER(len=*), INTENT(in) :: leading_text
    CHARACTER(len=*), INTENT(in) :: message_text

    WRITE (0, *) TRIM(leading_text)//': '//TRIM(message_text)

  END SUBROUTINE message

END MODULE mo_exception

PROGRAM callback_test

  USE mtime

  USE mo_exception

  CALL setCalendar(proleptic_gregorian)

!LK  call register_print_mtime_procedure(message)
!LK
!LK  current_date => newDatetime('1999-01-01T00:00:00')
!LK  call print_mtime('Message','Works as expected!', current_date)
!LK
!LK  model_time_step => newTimedelta('PT2H')
!LK  call print_mtime('Message','Works as expected!', model_time_step)

END PROGRAM callback_test
