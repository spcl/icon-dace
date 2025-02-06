!> @file ppm_real_sp_dp_edit_descriptor.f90
!! @brief get g edit descriptors for non-truncating I/O of
!! REAL(sp) and REAL(dp) variables
!!
!! @copyright Copyright  (C)  2018  Thomas Jahns <jahns@dkrz.de>
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
#include "fc_feature_defs.inc"
MODULE ppm_real_sp_dp_edit_descriptor
  USE ppm_std_type_kinds, ONLY: sp, dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: get_edit_descriptor_sp, get_edit_descriptor_dp

CONTAINS

  SUBROUTINE get_edit_descriptor_sp(de_g_sp, de_g_sp_width)
    CHARACTER(len=*), INTENT(out) :: de_g_sp
    INTEGER, INTENT(out) :: de_g_sp_width
    INTEGER :: field_w, field_d, field_e, min_e, max_e
    REAL(sp) :: f
    ! width of fraction
    field_d = MERGE(DIGITS(f), &
         CEILING(1.0 + REAL(DIGITS(f)) * LOG10(REAL(RADIX(f)))), &
         RADIX(f) == 10)
    max_e = MAXEXPONENT(f)
    min_e = MINEXPONENT(f)
    ! width of exponent field
    field_e = MAX(&
         CEILING(LOG10(REAL(max_e)*LOG10(REAL(RADIX(f))))), &
         CEILING(LOG10(ABS(REAL(min_e))*LOG10(REAL(RADIX(f))))))
    ! extra space: sign, decimal separator
    ! and literal 'e' plus optional exponent sign
    field_w = 2 + field_d + 2 + field_e
    ! certain compilers versions fail with the above computed minimum width
    ! because of weak implementations, therefore we increase the field width
    ! slightly
    field_w = field_w + 1
    de_g_sp_width = field_w
    IF (field_e > 2) THEN
      WRITE(de_g_sp, '(3(a,i0))') 'g', field_w, '.', field_d, 'e', field_e
    ELSE
      WRITE(de_g_sp, '(2(a,i0))') 'g', field_w, '.', field_d
    END IF
  END SUBROUTINE get_edit_descriptor_sp

  SUBROUTINE get_edit_descriptor_dp(de_g_dp, de_g_dp_width)
    CHARACTER(len=*), INTENT(out) :: de_g_dp
    INTEGER, INTENT(out) :: de_g_dp_width
    INTEGER :: field_w, field_d, field_e, min_e, max_e
    REAL(dp) :: d

    ! width of fraction
    field_d = MERGE(DIGITS(d), &
         CEILING(1.0 + REAL(DIGITS(d)) * LOG10(REAL(RADIX(d)))), &
         RADIX(d) == 10)
    max_e = MAXEXPONENT(d)
    min_e = MINEXPONENT(d)
    ! width of exponent field
    field_e = MAX(&
         CEILING(LOG10(REAL(max_e)*LOG10(REAL(RADIX(d))))), &
         CEILING(LOG10(ABS(REAL(min_e))*LOG10(REAL(RADIX(d))))))
    ! extra space: sign, decimal separator
    ! and literal 'e' plus optional exponent sign
    field_w = 2 + field_d + 2 + field_e
    ! certain compiler versions fail for some edge cases,
    ! e.g. ifort cannot print -6.0191114046301362d+145 for above field width
    field_w = field_w + 1
    de_g_dp_width = field_w
    IF (field_e > 2) THEN
      WRITE(de_g_dp, '(3(a,i0))') 'g', field_w, '.', field_d, 'e', field_e
    ELSE
      WRITE(de_g_dp, '(2(a,i0))') 'g', field_w, '.', field_d
    END IF
  END SUBROUTINE get_edit_descriptor_dp

END MODULE ppm_real_sp_dp_edit_descriptor
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
