!>
!! @file test_ieee_emulation.f90
!! @brief test ieee_arithmetic module implementation
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: IEEE 754 floating-point arithmetics
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
PROGRAM test_ieee_emulation
  USE ppm_base, ONLY: assertion
  USE ppm_std_type_kinds, ONLY: sp, dp
  USE ieee_arithmetic, ONLY: ieee_value, &
       ieee_quiet_nan, ieee_copy_sign, &
       ieee_negative_normal, ieee_negative_denormal, ieee_negative_zero, &
       ieee_positive_zero, ieee_positive_denormal, ieee_positive_normal, &
       ieee_positive_inf, ieee_negative_inf, ieee_is_normal, &
       ieee_support_denormal, ieee_support_nan, ieee_support_inf
  IMPLICIT NONE
  REAL(sp) :: vs, vs_sign
  REAL(dp) :: vd, vd_sign
  ! PGI 16.x and later might have this problem still, but I don't have a way to test
#if defined __PGI && __PGIC__ <= 17
  REAL(sp), PARAMETER :: min_normal_sp = 2.0_sp**MINEXPONENT(1.0_sp)
  REAL(dp), PARAMETER :: min_normal_dp = 2.0_dp**MINEXPONENT(1.0_dp)
#else
  REAL(sp), PARAMETER :: min_normal_sp = SCALE(1.0_sp, MINEXPONENT(1.0_sp))
  REAL(dp), PARAMETER :: min_normal_dp = SCALE(1.0_dp, MINEXPONENT(1.0_dp))
#endif
  CHARACTER(len=*), PARAMETER :: filename = 'test_ieee_emulation.f90'
  IF (ieee_support_nan(1.0_sp)) THEN
    CALL assertion(.NOT. ieee_is_normal(ieee_value(1.0_sp, ieee_quiet_nan)), &
         filename, __LINE__, 'quiet nan is normal?')
  END IF
  IF (ieee_support_inf(1.0_sp)) THEN
    vs = ieee_value(1.0_sp, ieee_negative_inf)
    CALL assertion(vs < 0.0_sp .AND. .NOT. ieee_is_normal(vs), &
         filename, __LINE__, 'negative infinity is normal or >= 0.0 ?')
  END IF
  vs = ieee_value(1.0_sp, ieee_negative_normal)
  CALL assertion(vs < 0.0_sp .AND. ieee_is_normal(vs) &
       .AND. vs <= min_normal_sp, &
       filename, __LINE__, 'negative normal is not normal or > 0.0 ?')
  IF (ieee_support_denormal(1.0_sp)) THEN
    vs = ieee_value(1.0_sp, ieee_negative_denormal)
    vs_sign = ieee_copy_sign(1.0_sp, vs)
    CALL assertion(vs > -min_normal_sp .AND. vs <= 0.0_sp &
         .AND. vs_sign <= -1.0_sp .AND. vs_sign >= -1.0_sp &
         .AND. .NOT. ieee_is_normal(vs), &
         filename, __LINE__, 'negative denormal is normal or > 0.0 ?')
  END IF
  vs = ieee_value(1.0_sp, ieee_negative_zero)
  CALL assertion(vs >= 0.0_sp .AND. vs <= 0.0_sp .AND. ieee_is_normal(vs), &
       filename, __LINE__, 'negative zero is not normal or /= 0.0 ?')
  vs = ieee_value(1.0_sp, ieee_positive_zero)
  CALL assertion(vs <= 0.0_sp .AND. vs >= 0.0_sp .AND. ieee_is_normal(vs), &
       filename, __LINE__, 'negative zero is not normal or /= 0.0 ?')
  IF (ieee_support_denormal(1.0_sp)) THEN
    vs = ieee_value(1.0_sp, ieee_positive_denormal)
    vs_sign = ieee_copy_sign(1.0_sp, vs)
    CALL assertion(vs < min_normal_sp .AND. vs >= 0.0_sp &
         .AND. vs_sign <= 1.0_sp .AND. vs_sign >= 1.0_sp &
         .AND. .NOT. ieee_is_normal(vs), &
         filename, __LINE__, 'positive denormal is normal or <= 0.0 ?')
  END IF
  vs = ieee_value(1.0_sp, ieee_positive_normal)
  CALL assertion(vs > 0.0_sp .AND. ieee_is_normal(vs), &
       filename, __LINE__, 'positive normal is not normal or <= 0.0 ?')
  IF (ieee_support_inf(1.0_sp)) THEN
    vs = ieee_value(1.0_sp, ieee_positive_inf)
    CALL assertion(vs > 0.0_sp .AND. .NOT. ieee_is_normal(vs), &
         filename, __LINE__, 'positive infinity is normal or <= 0.0 ?')
  END IF
  IF (ieee_support_nan(1.0_dp)) THEN
    CALL assertion(.NOT. ieee_is_normal(ieee_value(1.0_dp, ieee_quiet_nan)), &
         filename, __LINE__, 'quiet nan is normal?')
  END IF
  IF (ieee_support_inf(1.0_dp)) THEN
    vd = ieee_value(1.0_dp, ieee_negative_inf)
    CALL assertion(vd < 0.0_dp .AND. .NOT. ieee_is_normal(vd), &
         filename, __LINE__, 'negative infinity is normal or >= 0.0 ?')
  END IF
  vd = ieee_value(1.0_dp, ieee_negative_normal)
  CALL assertion(vd < 0.0_dp .AND. ieee_is_normal(vd), &
       filename, __LINE__, 'negative normal is not normal or >= 0.0 ?')
  IF (ieee_support_denormal(1.0_dp)) THEN
    vd = ieee_value(1.0_dp, ieee_negative_denormal)
    vd_sign = ieee_copy_sign(1.0_dp, vd)
    CALL assertion(vd > -min_normal_dp .AND. vd <= 0.0_dp &
         .AND. vd_sign <= -1.0_dp .AND. vd_sign >= -1.0_dp &
         .AND. .NOT. ieee_is_normal(vd), &
         filename, __LINE__, 'negative denormal is normal or > 0.0 ?')
  END IF
  vd = ieee_value(1.0_dp, ieee_negative_zero)
  CALL assertion(vd >= 0.0_dp .AND. vd <= 0.0_dp .AND. ieee_is_normal(vd), &
       filename, __LINE__, 'negative zero is not normal or /= 0.0 ?')
  vd = ieee_value(1.0_dp, ieee_positive_zero)
  CALL assertion(vd >= 0.0_dp .AND. vd <= 0.0_dp .AND. ieee_is_normal(vd), &
       filename, __LINE__, 'negative zero is not normal or /= 0.0 ?')
  IF (ieee_support_denormal(1.0_dp)) THEN
    vd = ieee_value(1.0_dp, ieee_positive_denormal)
    vd_sign = ieee_copy_sign(1.0_dp, vd)
    CALL assertion(vd < min_normal_dp .AND. vd >= 0.0_dp &
         .AND. vd_sign <= 1.0_dp .AND. vd_sign >= 1.0_dp &
         .AND. .NOT. ieee_is_normal(vd), &
         filename, __LINE__, 'positive denormal is normal or <= 0.0 ?')
  END IF
  vd = ieee_value(1.0_dp, ieee_positive_normal)
  CALL assertion(vd > 0.0_dp .AND. ieee_is_normal(vd), &
       filename, __LINE__, 'positive normal is not normal or <= 0.0 ?')
  IF (ieee_support_inf(1.0_dp)) THEN
    vd = ieee_value(1.0_dp, ieee_positive_inf)
    CALL assertion(vd > 0.0_dp .AND. .NOT. ieee_is_normal(vd), &
         filename, __LINE__, 'positive infinity is normal or <= 0.0 ?')
  END IF
END PROGRAM test_ieee_emulation
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
