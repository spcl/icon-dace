!>
!! @file ieee_arithmetic.f90
!! @brief emulation of standard functions on older compilers
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: IEEE arithmetic emulation module
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
!> emulation of standard functions on older compilers

MODULE ieee_arithmetic
  USE ppm_std_type_kinds, ONLY: dp, sp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ieee_is_normal
  INTERFACE ieee_is_normal
    ELEMENTAL FUNCTION ppm_ieee_is_normal_sp(x)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(in) :: x
      LOGICAL :: ppm_ieee_is_normal_sp
    END FUNCTION ppm_ieee_is_normal_sp
    ELEMENTAL FUNCTION ppm_ieee_is_normal_dp(x)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(in) :: x
      LOGICAL :: ppm_ieee_is_normal_dp
    END FUNCTION ppm_ieee_is_normal_dp
  END INTERFACE
  PUBLIC :: ieee_is_nan
  INTERFACE ieee_is_nan
    ELEMENTAL FUNCTION ppm_ieee_is_nan_sp(x)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(in) :: x
      LOGICAL :: ppm_ieee_is_nan_sp
    END FUNCTION ppm_ieee_is_nan_sp
    ELEMENTAL FUNCTION ppm_ieee_is_nan_dp(x)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(in) :: x
      LOGICAL :: ppm_ieee_is_nan_dp
    END FUNCTION ppm_ieee_is_nan_dp
  END INTERFACE ieee_is_nan
  PUBLIC :: ieee_class_type
  TYPE ieee_class_type
    SEQUENCE
    INTEGER :: enum_val
  END TYPE ieee_class_type
  TYPE(ieee_class_type), PUBLIC, PARAMETER :: &

       ieee_signaling_nan = ieee_class_type(1), &

       ieee_quiet_nan = ieee_class_type(2), &
       ieee_negative_inf = ieee_class_type(3), &
       ieee_negative_normal = ieee_class_type(4), &
       ieee_negative_denormal = ieee_class_type(5), &
       ieee_negative_zero = ieee_class_type(6), &
       ieee_positive_zero = ieee_class_type(7), &
       ieee_positive_denormal = ieee_class_type(8), &
       ieee_positive_normal = ieee_class_type(9), &
       ieee_positive_inf = ieee_class_type(10)
  PUBLIC :: ieee_value
  INTERFACE ieee_value
    ELEMENTAL FUNCTION ppm_ieee_value_sp(x, v)
      USE ppm_std_type_kinds, ONLY: sp
      IMPORT :: ieee_class_type
      REAL(sp), INTENT(in) :: x
      TYPE(ieee_class_type), INTENT(in) :: v
      REAL(sp) :: ppm_ieee_value_sp
    END FUNCTION ppm_ieee_value_sp
    ELEMENTAL FUNCTION ppm_ieee_value_dp(x, v)
      USE ppm_std_type_kinds, ONLY: dp
      IMPORT :: ieee_class_type
      REAL(dp), INTENT(in) :: x
      TYPE(ieee_class_type), INTENT(in) :: v
      REAL(dp) :: ppm_ieee_value_dp
    END FUNCTION ppm_ieee_value_dp
  END INTERFACE ieee_value
  PUBLIC :: ieee_copy_sign
  INTERFACE ieee_copy_sign
    ELEMENTAL FUNCTION ppm_copy_sign_sp(v, s)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(in) :: v, s
      REAL(sp) :: ppm_copy_sign_sp
    END FUNCTION ppm_copy_sign_sp
    ELEMENTAL FUNCTION ppm_copy_sign_dp(v, s)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(in) :: v, s
      REAL(dp) :: ppm_copy_sign_dp
    END FUNCTION ppm_copy_sign_dp
  END INTERFACE ieee_copy_sign
  PUBLIC :: ieee_support_denormal
  INTERFACE ieee_support_denormal
    ELEMENTAL FUNCTION ppm_support_denormal_sp(x)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(in) :: x
      LOGICAL :: ppm_support_denormal_sp
    END FUNCTION ppm_support_denormal_sp
    ELEMENTAL FUNCTION ppm_support_denormal_dp(x)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(in) :: x
      LOGICAL :: ppm_support_denormal_dp
    END FUNCTION ppm_support_denormal_dp
  END INTERFACE ieee_support_denormal
  PUBLIC :: ieee_support_nan
  INTERFACE ieee_support_nan
    ELEMENTAL FUNCTION ppm_support_nan_sp(x)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(in) :: x
      LOGICAL :: ppm_support_nan_sp
    END FUNCTION ppm_support_nan_sp
    ELEMENTAL FUNCTION ppm_support_nan_dp(x)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(in) :: x
      LOGICAL :: ppm_support_nan_dp
    END FUNCTION ppm_support_nan_dp
  END INTERFACE ieee_support_nan
  PUBLIC :: ieee_support_inf
  INTERFACE ieee_support_inf
    ELEMENTAL FUNCTION ppm_support_inf_sp(x)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(in) :: x
      LOGICAL :: ppm_support_inf_sp
    END FUNCTION ppm_support_inf_sp
    ELEMENTAL FUNCTION ppm_support_inf_dp(x)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(in) :: x
      LOGICAL :: ppm_support_inf_dp
    END FUNCTION ppm_support_inf_dp
  END INTERFACE ieee_support_inf
END MODULE ieee_arithmetic
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
