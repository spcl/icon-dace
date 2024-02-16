!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
module mo_exception

  implicit none
  
  private

  public :: message
  
contains

  subroutine message(leading_text, message_text)
    character(len=*), intent(in) :: leading_text
    character(len=*), intent(in) :: message_text    

    write (0,*) trim(leading_text)//': '//trim(message_text)
    
  end subroutine message

end module mo_exception

program callback_test

  use mtime
  
  use mo_exception

  type(datetime), pointer :: current_date => null() 
  type(timedelta), pointer :: model_time_step => null()

  call setCalendar(proleptic_gregorian)
  
!LK  call register_print_mtime_procedure(message)
!LK
!LK  current_date => newDatetime('1999-01-01T00:00:00')
!LK  call print_mtime('Message','Works as expected!', current_date)
!LK
!LK  model_time_step => newTimedelta('PT2H')
!LK  call print_mtime('Message','Works as expected!', model_time_step)  
  
end program callback_test
