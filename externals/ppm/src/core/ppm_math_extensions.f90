!>
!! @file ppm_math_extensions.f90
!! @brief utility routines for floating-point math
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
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




MODULE ppm_math_extensions





  USE ppm_std_type_kinds, ONLY: dp, sp
  USE ppm_math_extensions_internal, ONLY: fpu_save_cw => ppm_fpu_save_cw, &
       fpu_restore_cw => ppm_fpu_restore_cw, &
       fpu_set_precision, fpu_precision_sp, fpu_precision_dp, fpu_precision_ep,&
       fpu_restore_mxcsr => ppm_fpu_restore_mxcsr, &
       fpu_save_mxcsr => ppm_fpu_save_mxcsr, &
       fpu_set_abrupt_underflow



  IMPLICIT NONE
  PRIVATE



  REAL(dp), PARAMETER, PUBLIC :: &
       m_pi_dp=3.14159265358979323846_dp
  REAL(sp), PARAMETER, PUBLIC :: &
       m_pi_sp=3.14159265358979323846_sp

! ppm_real_sp_dp_edit_descriptor.inc.in --- define edit descriptors for REAL kinds sp and dp
!
! Copyright  (C)  2018  Thomas Jahns <jahns@dkrz.de>
!
! Version: 1.0
! Author: Thomas Jahns <jahns@dkrz.de>
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
!
! Local Variables:
! mode: f90
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
  !> data edit descriptors for real kinds sp and dp
  CHARACTER(len=*), PARAMETER :: &
       de_g_sp = "g16.9", &
       de_g_dp = "g25.17e3"
  !> character width needed for the corresponding data edit descriptor
  INTEGER, PARAMETER :: de_g_sp_width = 16, &
       de_g_dp_width = 25

  !   The ddp functions compute complex results that store a Kahan
  !   correction term in the imaginary part
  INTERFACE ddp_add
    MODULE PROCEDURE ddp_add_ddp_ddp
    MODULE PROCEDURE ddp_add_dp_dp
  END INTERFACE ddp_add
  INTERFACE ddp_sum
    MODULE PROCEDURE ddp_sum_dp_1d
    MODULE PROCEDURE ddp_sum_dp_2d
    MODULE PROCEDURE ddp_sum_dp_3d
  END INTERFACE ddp_sum
  PUBLIC :: ddp_add, ddp_sum, ddp_abs
  PUBLIC :: fpu_save_cw, fpu_restore_cw, fpu_set_precision, &
       fpu_precision_sp, fpu_precision_dp, &
       fpu_precision_ep, fpu_set_abrupt_underflow, &
       fpu_save_mxcsr, fpu_restore_mxcsr
  PUBLIC :: initialize_math_extensions, finalize_math_extensions
  PUBLIC :: de_g_sp, de_g_dp, de_g_sp_width, de_g_dp_width

  INTERFACE
    PURE SUBROUTINE ppm_ddp_sum_dp(n, a, s)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(in) :: a(n)
      COMPLEX(dp), INTENT(out) :: s
    END SUBROUTINE ppm_ddp_sum_dp

    ELEMENTAL SUBROUTINE ppm_ddp_add_dp_dp(a, b, s)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(in) :: a, b
      COMPLEX(dp), INTENT(out) :: s
    END SUBROUTINE ppm_ddp_add_dp_dp

    ELEMENTAL SUBROUTINE ppm_ddp_add_ddp_ddp(a, b, s)
      USE ppm_std_type_kinds, ONLY: dp
      COMPLEX(dp), INTENT(in) :: a, b
      COMPLEX(dp), INTENT(out) :: s
    END SUBROUTINE ppm_ddp_add_ddp_ddp
  END INTERFACE

  INTERFACE assign_nan
    PURE SUBROUTINE ppm_assign_nan_dp(v)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(out) :: v
    END SUBROUTINE ppm_assign_nan_dp
    PURE SUBROUTINE ppm_assign_nan_sp(v)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(out) :: v
    END SUBROUTINE ppm_assign_nan_sp
  END INTERFACE assign_nan
  PUBLIC :: assign_nan
  CHARACTER(len=*), PARAMETER :: filename = 'ppm_math_extensions.f90'
CONTAINS
  !   Modification of original codes written by David H. Bailey
  !>  This function computes c = a(i)+b(i) but gives as result
  !!  a complex that stores a Kahan correction term in the imaginary part
  ELEMENTAL FUNCTION ddp_add_ddp_ddp(a, b) RESULT(c)
    COMPLEX(dp) :: c
    COMPLEX(dp), INTENT(in) :: a, b
    CALL ppm_ddp_add_ddp_ddp(a, b, c)
  END FUNCTION ddp_add_ddp_ddp

  ELEMENTAL FUNCTION ddp_add_dp_dp(a, b) RESULT(c)
    COMPLEX(dp) :: c
    REAL(dp), INTENT(in) :: a, b
    CALL ppm_ddp_add_dp_dp(a, b, c)
  END FUNCTION ddp_add_dp_dp

  ELEMENTAL FUNCTION ddp_abs(a) RESULT(c)
    COMPLEX(dp) :: c
    COMPLEX(dp), INTENT(in) :: a
    c = MERGE(a, -a, REAL(a, dp) >= 0.0_dp)
  END FUNCTION ddp_abs

  !> compute double-double-precision corrected sum of 1d-array, should
  !! give better results than sum(a) for arrays a where cancellation
  !! occurs.
  !! @param a array a(1)..a(n) to sum up
  !! @return \f$\sum^n_{i=1} a_i\f$
  PURE FUNCTION ddp_sum_dp_1d(a) RESULT(c)
    COMPLEX(dp) :: c
    REAL(dp), INTENT(in) :: a(:)

    CALL ppm_ddp_sum_dp(SIZE(a), a, c)
  END FUNCTION ddp_sum_dp_1d

  !> compute double-double-precision corrected sum of 2d-array
  !> @see ddp_sum_dp_1d
  !> @param a array a(1,1)..a(m,n) to sum up
  !> @return \f$\sum^m_{i=1}\sum^n_{j=1} a_{i,j}\f$
  PURE FUNCTION ddp_sum_dp_2d(a) RESULT(c)
    COMPLEX(dp) :: c
    REAL(dp), INTENT(in) :: a(:, :)
    CALL ppm_ddp_sum_dp(SIZE(a), a, c)
  END FUNCTION ddp_sum_dp_2d

  !> compute double-double-precision corrected sum of 3d-array
  !> @see ddp_sum_dp_1d
  !> @param a array a(1,1,1)..a(m,n,o) to sum up
  !> @return \f$\sum^m_{i=1}\sum^n_{j=1}\sum^o_{l=1} a_{i,j,l}\f$
  PURE FUNCTION ddp_sum_dp_3d(a) RESULT(c)
    COMPLEX(dp) :: c
    REAL(dp), INTENT(in) :: a(:, :, :)
    CALL ppm_ddp_sum_dp(SIZE(a), a, c)
  END FUNCTION ddp_sum_dp_3d


  SUBROUTINE initialize_math_extensions

  END SUBROUTINE initialize_math_extensions

  SUBROUTINE finalize_math_extensions
  END SUBROUTINE finalize_math_extensions

END MODULE ppm_math_extensions
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
