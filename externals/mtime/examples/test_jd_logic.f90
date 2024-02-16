!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program test_jd_logic

  use mtime
  
  implicit none

  type(julianday) :: jd1, jd2
  
  jd1 = julianday( 2459904, 0 )
  jd2 = julianday( 2459914, 0 )
  call jd1_less_than_jd2(jd1, jd2)

  jd1 = julianday( 2459914, 0 )
  jd2 = julianday( 2459904, 0 )
  call jd1_greater_than_jd2(jd1, jd2)

  jd1 = julianday( 2459924, 0 )
  jd2 = julianday( 2459924, 0 )
  call jd1_equals_jd2(jd1, jd2)

  
contains

  subroutine jd1_less_than_jd2(jd1, jd2)
    type(julianday), intent(in) :: jd1, jd2
    write(0,*) 'result ', jd1 < jd2,  ' (expect T) ' , jd1%day, ' <  ',  jd2%day
    write(0,*) 'result ', jd1 > jd2,  ' (expect F) ' , jd1%day, ' >  ',  jd2%day
    write(0,*) 'result ', jd1 <= jd2, ' (expect T) ' , jd1%day, ' <= ', jd2%day
    write(0,*) 'result ', jd1 >= jd2, ' (expect F) ' , jd1%day, ' >= ', jd2%day
    write(0,*) 'result ', jd1 == jd2, ' (expect F) ' , jd1%day, ' == ', jd2%day
    write(0,*) 'result ', jd1 /= jd2, ' (expect T) ' , jd1%day, ' /= ', jd2%day            
  end subroutine jd1_less_than_jd2

  subroutine jd1_greater_than_jd2(jd1, jd2)
    type(julianday), intent(in) :: jd1, jd2
    write(0,*) 'result ', jd1 < jd2,  ' (expect F) ' , jd1%day, ' <  ', jd2%day
    write(0,*) 'result ', jd1 > jd2,  ' (expect T) ' , jd1%day, ' >  ', jd2%day
    write(0,*) 'result ', jd1 <= jd2, ' (expect F) ' , jd1%day, ' <= ', jd2%day
    write(0,*) 'result ', jd1 >= jd2, ' (expect T) ' , jd1%day, ' >= ', jd2%day
    write(0,*) 'result ', jd1 == jd2, ' (expect F) ' , jd1%day, ' == ', jd2%day
    write(0,*) 'result ', jd1 /= jd2, ' (expect T) ' , jd1%day, ' /= ', jd2%day        
  end subroutine jd1_greater_than_jd2

  subroutine jd1_equals_jd2(jd1, jd2)
    type(julianday), intent(in) :: jd1, jd2
    write(0,*) 'result ', jd1 < jd2,  ' (expect F) ' , jd1%day, ' <  ', jd2%day
    write(0,*) 'result ', jd1 > jd2,  ' (expect F) ' , jd1%day, ' >  ', jd2%day
    write(0,*) 'result ', jd1 <= jd2, ' (expect T) ' , jd1%day, ' <= ', jd2%day
    write(0,*) 'result ', jd1 >= jd2, ' (expect T) ' , jd1%day, ' >= ', jd2%day
    write(0,*) 'result ', jd1 == jd2, ' (expect T) ' , jd1%day, ' == ', jd2%day
    write(0,*) 'result ', jd1 /= jd2, ' (expect F) ' , jd1%day, ' /= ', jd2%day        
  end subroutine jd1_equals_jd2
  
end program test_jd_logic
