!>
!! @file ppm_math_extensions_internal.f90
!! @brief wrapper for precision control functions needed in ppm_math_extensions
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
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

MODULE ppm_math_extensions_internal
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER, PUBLIC :: fpu_precision_sp = 1, &
       fpu_precision_dp = 2, fpu_precision_ep = 3
  INTERFACE
    ! these routines are not actually pure, but we want to be able to
    ! use them in some that are, and if used correctly, i.e. when
    ! always restoring the control registers before returning, they
    ! are without side-effects
    PURE SUBROUTINE ppm_fpu_save_cw(fpu_cw)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(out) :: fpu_cw
    END SUBROUTINE ppm_fpu_save_cw

    PURE SUBROUTINE ppm_fpu_restore_cw(fpu_cw)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(in) :: fpu_cw
    END SUBROUTINE ppm_fpu_restore_cw

    PURE SUBROUTINE ppm_fpu_save_mxcsr(old_mxcsr)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(out) :: old_mxcsr
    END SUBROUTINE ppm_fpu_save_mxcsr

    PURE SUBROUTINE ppm_fpu_restore_mxcsr(old_mxcsr)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(in) :: old_mxcsr
    END SUBROUTINE ppm_fpu_restore_mxcsr

  END INTERFACE
  PUBLIC :: ppm_fpu_save_cw, ppm_fpu_restore_cw, fpu_set_precision
  PUBLIC :: ppm_fpu_save_mxcsr, ppm_fpu_restore_mxcsr, fpu_set_abrupt_underflow
CONTAINS
  PURE SUBROUTINE fpu_set_precision(fpu_precision, old_fpu_cw)
    INTEGER, INTENT(in) :: fpu_precision
    INTEGER(i4), OPTIONAL, INTENT(out) :: old_fpu_cw

    INTERFACE
      PURE SUBROUTINE ppm_fpu_set_precision_c(fpu_precision, old_fpu_cw)
        IMPORT :: i4
        INTEGER, INTENT(in) :: fpu_precision
        INTEGER(i4), INTENT(out) :: old_fpu_cw
      END SUBROUTINE ppm_fpu_set_precision_c
    END INTERFACE

    INTEGER(i4) :: temp

    CALL ppm_fpu_set_precision_c(fpu_precision, temp)
    IF (PRESENT(old_fpu_cw)) old_fpu_cw = temp
  END SUBROUTINE fpu_set_precision

  PURE SUBROUTINE fpu_set_abrupt_underflow(abrupt_underflow, old_fpu_mxcsr)
    LOGICAL, INTENT(in) :: abrupt_underflow
    INTEGER(i4), OPTIONAL, INTENT(out) :: old_fpu_mxcsr

    INTERFACE
      PURE SUBROUTINE ppm_fpu_set_abrupt_underflow_c(old_fpu_mxcsr, &
           abrupt_underflow)
        IMPORT :: i4
        INTEGER(i4), INTENT(out) :: old_fpu_mxcsr
        LOGICAL, INTENT(in) :: abrupt_underflow
      END SUBROUTINE ppm_fpu_set_abrupt_underflow_c
    END INTERFACE
    INTEGER(i4) :: temp

    CALL ppm_fpu_set_abrupt_underflow_c(temp, abrupt_underflow)
    IF (PRESENT(old_fpu_mxcsr)) old_fpu_mxcsr = temp
  END SUBROUTINE fpu_set_abrupt_underflow
END MODULE ppm_math_extensions_internal
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
