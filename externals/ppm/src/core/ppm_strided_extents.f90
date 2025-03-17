!>
!! @file ppm_strided_extents.f90
!! @brief extend extent type to handle strided access
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords:
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

MODULE ppm_strided_extents



  USE ppm_extents, ONLY: extent, extent_size, extent_start, &
       extent_end, OPERATOR(==)



  IMPLICIT NONE
  PRIVATE



  ! succinct description of subarrays
  TYPE strided_extent
    TYPE(extent) :: ext
    INTEGER :: stride
  END TYPE strided_extent
  !> string representation of extent size/position takes
  !! 11 decimal places (10 + sign)
  INTEGER, PARAMETER :: sext_i2s_len=11

  PUBLIC :: strided_extent, extent_size, char, extent_start, extent_end, &
       OPERATOR(==)

  INTERFACE extent_size
    MODULE PROCEDURE strided_extent_size_1d
    MODULE PROCEDURE strided_extent_size_nd
  END INTERFACE

  INTERFACE extent_start
    MODULE PROCEDURE strided_extent_start_1d
  END INTERFACE

  INTERFACE extent_end
    MODULE PROCEDURE strided_extent_end_1d
  END INTERFACE

  INTERFACE char
    MODULE PROCEDURE char_auto
  END INTERFACE

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE strided_extent_equality
  END INTERFACE




  CHARACTER(len=*), PARAMETER :: filename = 'ppm_strided_extents.f90'
CONTAINS

  FUNCTION strided_extent_size_1d(sext) RESULT(sext_size)
    INTEGER :: sext_size
    TYPE(strided_extent), INTENT(in) :: sext

    sext_size = (extent_size(sext%ext) + sext%stride - SIGN(1, sext%stride)) &
         / sext%stride
  END FUNCTION strided_extent_size_1d

  FUNCTION strided_extent_size_nd(sext) RESULT(sext_size)
    INTEGER :: sext_size
    TYPE(strided_extent), INTENT(in) :: sext(:)

    INTEGER :: i
    sext_size = 1
    DO i = 1, SIZE(sext)
      sext_size = sext_size * extent_size(sext(i)%ext)/sext(i)%stride
    END DO
  END FUNCTION strided_extent_size_nd

  ELEMENTAL FUNCTION strided_extent_start_1d(sext) RESULT(sext_start)
    INTEGER :: sext_start
    TYPE(strided_extent), INTENT(in) :: sext
    sext_start = extent_start(sext%ext)
  END FUNCTION strided_extent_start_1d

  ELEMENTAL FUNCTION strided_extent_end_1d(sext) RESULT(sext_end)
    INTEGER :: sext_end
    TYPE(strided_extent), INTENT(in) :: sext
    sext_end = extent_end(sext%ext) - MOD(extent_end(sext%ext), sext%stride)
  END FUNCTION strided_extent_end_1d

  ELEMENTAL FUNCTION char_auto(sext) RESULT(str)
    CHARACTER(len=3*sext_i2s_len+5) :: str
    TYPE(strided_extent), INTENT(in) :: sext
    IF (sext%ext%size /= 0) THEN
      WRITE (str, '(3(a,i0),a)') '[', sext%ext%first, ',', extent_end(sext), &
           ',', sext%stride, ']'
    ELSE
      str = '{}'
    END IF
  END FUNCTION char_auto

  ELEMENTAL FUNCTION strided_extent_equality(a, b) RESULT(l)
    TYPE(strided_extent), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%ext == b%ext .AND. a%stride == b%stride
  END FUNCTION strided_extent_equality


END MODULE ppm_strided_extents
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
