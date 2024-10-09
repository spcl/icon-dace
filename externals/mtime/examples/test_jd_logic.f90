!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM test_jd_logic

  USE mtime

  IMPLICIT NONE

  TYPE(julianday) :: jd1, jd2

  jd1 = julianday(2459904, 0)
  jd2 = julianday(2459914, 0)
  CALL jd1_less_than_jd2(jd1, jd2)

  jd1 = julianday(2459914, 0)
  jd2 = julianday(2459904, 0)
  CALL jd1_greater_than_jd2(jd1, jd2)

  jd1 = julianday(2459924, 0)
  jd2 = julianday(2459924, 0)
  CALL jd1_equals_jd2(jd1, jd2)

CONTAINS

  SUBROUTINE jd1_less_than_jd2(jd1, jd2)
    TYPE(julianday), INTENT(in) :: jd1, jd2
    WRITE (0, *) 'result ', jd1 < jd2, ' (expect T) ', jd1%day, ' <  ', jd2%day
    WRITE (0, *) 'result ', jd1 > jd2, ' (expect F) ', jd1%day, ' >  ', jd2%day
    WRITE (0, *) 'result ', jd1 <= jd2, ' (expect T) ', jd1%day, ' <= ', jd2%day
    WRITE (0, *) 'result ', jd1 >= jd2, ' (expect F) ', jd1%day, ' >= ', jd2%day
    WRITE (0, *) 'result ', jd1 == jd2, ' (expect F) ', jd1%day, ' == ', jd2%day
    WRITE (0, *) 'result ', jd1 /= jd2, ' (expect T) ', jd1%day, ' /= ', jd2%day
  END SUBROUTINE jd1_less_than_jd2

  SUBROUTINE jd1_greater_than_jd2(jd1, jd2)
    TYPE(julianday), INTENT(in) :: jd1, jd2
    WRITE (0, *) 'result ', jd1 < jd2, ' (expect F) ', jd1%day, ' <  ', jd2%day
    WRITE (0, *) 'result ', jd1 > jd2, ' (expect T) ', jd1%day, ' >  ', jd2%day
    WRITE (0, *) 'result ', jd1 <= jd2, ' (expect F) ', jd1%day, ' <= ', jd2%day
    WRITE (0, *) 'result ', jd1 >= jd2, ' (expect T) ', jd1%day, ' >= ', jd2%day
    WRITE (0, *) 'result ', jd1 == jd2, ' (expect F) ', jd1%day, ' == ', jd2%day
    WRITE (0, *) 'result ', jd1 /= jd2, ' (expect T) ', jd1%day, ' /= ', jd2%day
  END SUBROUTINE jd1_greater_than_jd2

  SUBROUTINE jd1_equals_jd2(jd1, jd2)
    TYPE(julianday), INTENT(in) :: jd1, jd2
    WRITE (0, *) 'result ', jd1 < jd2, ' (expect F) ', jd1%day, ' <  ', jd2%day
    WRITE (0, *) 'result ', jd1 > jd2, ' (expect F) ', jd1%day, ' >  ', jd2%day
    WRITE (0, *) 'result ', jd1 <= jd2, ' (expect T) ', jd1%day, ' <= ', jd2%day
    WRITE (0, *) 'result ', jd1 >= jd2, ' (expect T) ', jd1%day, ' >= ', jd2%day
    WRITE (0, *) 'result ', jd1 == jd2, ' (expect T) ', jd1%day, ' == ', jd2%day
    WRITE (0, *) 'result ', jd1 /= jd2, ' (expect F) ', jd1%day, ' /= ', jd2%day
  END SUBROUTINE jd1_equals_jd2

END PROGRAM test_jd_logic
