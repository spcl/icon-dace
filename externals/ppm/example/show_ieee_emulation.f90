!> @file show_ieee_emulation.f90
!! @brief print attributes of IEEE_ARITHMETIC implementation
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
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
  USE ppm_std_type_kinds, ONLY: sp, dp
  USE ieee_arithmetic, ONLY: ieee_value, ieee_negative_inf, &
#ifdef HAVE_IEEE_SIGNALING_NAN
       ieee_signaling_nan, &
#endif
       ieee_quiet_nan, &
       ieee_negative_normal, ieee_negative_denormal, ieee_negative_zero, &
       ieee_positive_zero, ieee_positive_denormal, ieee_positive_normal, &
       ieee_positive_inf
  IMPLICIT NONE
#ifdef HAVE_IEEE_SIGNALING_NAN
  PRINT *, 'signaling NaN:     ', ieee_value(1.0, ieee_signaling_nan)
#endif
  PRINT *, 'quiet NaN:         ', ieee_value(1.0_sp, ieee_quiet_nan)
  PRINT *, 'negative infinity: ', ieee_value(1.0_sp, ieee_negative_inf)
  PRINT *, 'negative normal:   ', ieee_value(1.0_sp, ieee_negative_normal)
  PRINT *, 'negative denormal: ', ieee_value(1.0_sp, ieee_negative_denormal)
  PRINT *, 'negative zero:     ', ieee_value(1.0_sp, ieee_negative_zero)
  PRINT *, 'positive zero:     ', ieee_value(1.0_sp, ieee_positive_zero)
  PRINT *, 'positive denormal: ', ieee_value(1.0_sp, ieee_positive_denormal)
  PRINT *, 'positive normal:   ', ieee_value(1.0_sp, ieee_positive_normal)
  PRINT *, 'positive infinity: ', ieee_value(1.0_sp, ieee_positive_inf)
#ifdef HAVE_IEEE_SIGNALING_NAN
  PRINT *, 'signaling NaN:     ', ieee_value(1.0_dp, ieee_signaling_nan)
#endif
  PRINT *, 'quiet NaN:         ', ieee_value(1.0_dp, ieee_quiet_nan)
  PRINT *, 'negative infinity: ', ieee_value(1.0_dp, ieee_negative_inf)
  PRINT *, 'negative normal:   ', ieee_value(1.0_dp, ieee_negative_normal)
  PRINT *, 'negative denormal: ', ieee_value(1.0_dp, ieee_negative_denormal)
  PRINT *, 'negative zero:     ', ieee_value(1.0_dp, ieee_negative_zero)
  PRINT *, 'positive zero:     ', ieee_value(1.0_dp, ieee_positive_zero)
  PRINT *, 'positive denormal: ', ieee_value(1.0_dp, ieee_positive_denormal)
  PRINT *, 'positive normal:   ', ieee_value(1.0_dp, ieee_positive_normal)
  PRINT *, 'positive infinity: ', ieee_value(1.0_dp, ieee_positive_inf)
END PROGRAM test_ieee_emulation
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
