!>
!! @file ppm_statistics.f90
!! @brief compute statistical properties
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: ScalES PPM statistics mean
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
! Commentary:
!
!
!
! Code:
!

MODULE ppm_statistics
  IMPLICIT NONE
  PRIVATE
  INTERFACE arithmetic_mean
    MODULE PROCEDURE arithmetic_mean_1d_r_p
    MODULE PROCEDURE arithmetic_mean_2d_r
  END INTERFACE
  PUBLIC :: arithmetic_mean
CONTAINS
  ! FIXME: replace internal computation with ddp or exact computation
  PURE FUNCTION arithmetic_mean_1d_r_p(a, mask) RESULT(mean)
    REAL, INTENT(in) :: a(:)
    LOGICAL, INTENT(in), OPTIONAL :: mask(SIZE(a))
    REAL :: mean
    IF (PRESENT(mask)) THEN
      mean = SUM(a, mask)/REAL(MAX(1, COUNT(mask)))
    ELSE
      mean = SUM(a)/REAL(MAX(1, SIZE(a)))
    END IF
  END FUNCTION arithmetic_mean_1d_r_p

  PURE FUNCTION arithmetic_mean_2d_r(a, mask) RESULT(mean)
    REAL, INTENT(in) :: a(:,:)
    LOGICAL, INTENT(in), OPTIONAL :: mask(SIZE(a, 1),SIZE(a, 2))
    REAL :: mean
    IF (PRESENT(mask)) THEN
      mean = SUM(a, mask)/REAL(MAX(1, COUNT(mask)))
    ELSE
      mean = SUM(a)/REAL(MAX(1, SIZE(a)))
    END IF
  END FUNCTION arithmetic_mean_2d_r
END MODULE ppm_statistics
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
