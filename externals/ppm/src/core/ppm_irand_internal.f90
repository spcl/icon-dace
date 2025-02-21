!>
!! @file ppm_irand_internal.f90
!! @brief internal definitions of random number generator
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

MODULE ppm_irand_internal
  USE ppm_std_type_kinds, ONLY: i4, i8, dp, sp
  USE ppm_extents, ONLY: iinterval, iinterval_sp, iinterval_dp



  IMPLICIT NONE
  PRIVATE
  !> irand returns integers in the range irand_min..irand_max
  INTEGER(i4), PARAMETER :: irand_min=-2147483647, irand_max=2147483647
  !> irand8 returns integers in the range irand8_min..irand8_max
  INTEGER(i8), PARAMETER :: irand8_min=-9223372036854775807_i8, &
       irand8_max=9223372036854775807_i8
  !> this function is implemented in an OpenMP-thread-safe means,
  !! returns integers within range [irand_min,irand_max]
  INTERFACE
    FUNCTION ppm_irand()
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4) :: ppm_irand
    END FUNCTION ppm_irand
  END INTERFACE
  !> this function is implemented in an OpenMP-thread-safe means,
  !! returns integers within range [irand8_min,irand8_max]
  INTERFACE
    FUNCTION ppm_irand8()
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER(i8) :: ppm_irand8
    END FUNCTION ppm_irand8
  END INTERFACE
  !> this function is implemented in an OpenMP-thread-safe means,
  !! returns integers within range [0,irand_max]
  INTERFACE
    FUNCTION ppm_irandp()
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4) :: ppm_irandp
    END FUNCTION ppm_irandp
  END INTERFACE

  !> this function is implemented in an OpenMP-thread-safe means,
  !! returns integers within range [0,irand8_max]
  INTERFACE
    FUNCTION ppm_irandp8()
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER(i8) :: ppm_irandp8
    END FUNCTION ppm_irandp8
  END INTERFACE

  !> these functions are implemented in an OpenMP-thread-safe means,
  !! return integers within the specified range, which must be non-empty
  INTERFACE ppm_irandr
    FUNCTION ppm_irandr(range)
      USE ppm_extents, ONLY: iinterval
      USE ppm_std_type_kinds, ONLY: i4
      TYPE(iinterval), INTENT(in) :: range
      INTEGER(i4) :: ppm_irandr
    END FUNCTION ppm_irandr

    FUNCTION ppm_irandr8(range)
      USE ppm_extents, ONLY: iinterval_i8
      USE ppm_std_type_kinds, ONLY: i8
      TYPE(iinterval_i8), INTENT(in) :: range
      INTEGER(i8) :: ppm_irandr8
    END FUNCTION ppm_irandr8
  END INTERFACE

  !> these functions are implemented in an OpenMP-thread-safe means,
  !! return REALs within range [0.0_dp,1.0_dp) or [0.0_sp,1.0_sp) respectively
  INTERFACE
    FUNCTION ppm_drandp()
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp) :: ppm_drandp
    END FUNCTION ppm_drandp
    FUNCTION ppm_frandp()
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp) :: ppm_frandp
    END FUNCTION ppm_frandp
  END INTERFACE

  !> these functions are implemented in an OpenMP-thread-safe means,
  !! return reals within range (-1.0_dp,1.0_dp) or (-1.0_sp,1.0_sp) respectively
  INTERFACE
    FUNCTION ppm_drand()
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp) :: ppm_drand
    END FUNCTION ppm_drand
    FUNCTION ppm_frand()
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp) :: ppm_frand
    END FUNCTION ppm_frand
  END INTERFACE

  !> these functions are implemented in an OpenMP-thread-safe means,
  !! return REALs within specified range
  INTERFACE
    FUNCTION ppm_drandr(range)
      USE ppm_std_type_kinds, ONLY: dp
      USE ppm_extents, ONLY: iinterval_dp
      TYPE(iinterval_dp), INTENT(in) :: range
      REAL(dp) :: ppm_drandr
    END FUNCTION ppm_drandr
    FUNCTION ppm_frandr(range)
      USE ppm_std_type_kinds, ONLY: sp
      USE ppm_extents, ONLY: iinterval_sp
      TYPE(iinterval_sp), INTENT(in) :: range
      REAL(sp) :: ppm_frandr
    END FUNCTION ppm_frandr
  END INTERFACE


  !> unfortunately, Fortrans random number generator is only prepared
  !! to produce REAL-type results, this add similar capabilities for
  !! INTEGER results in range [irand_min,irand_max] and
  !! [irand8_min,irand8_max] respectively and REALs in the range
  !! (-1.0,1.0)
  INTERFACE a_rand
    SUBROUTINE ppm_irand_a(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_irand_a
    SUBROUTINE ppm_irand_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_irand_a_2d
    SUBROUTINE ppm_irand_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_irand_a_3d
    SUBROUTINE ppm_irand8_a(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER(i8), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_irand8_a
    SUBROUTINE ppm_irand8_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER, INTENT(in) :: n
      INTEGER(i8), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_irand8_a_2d
    SUBROUTINE ppm_irand8_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER, INTENT(in) :: n
      INTEGER(i8), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_irand8_a_3d
    SUBROUTINE ppm_drand_a(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_drand_a
    SUBROUTINE ppm_drand_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_drand_a_2d
    SUBROUTINE ppm_drand_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_drand_a_3d
    SUBROUTINE ppm_frand_a(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_frand_a
    SUBROUTINE ppm_frand_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_frand_a_2d
    SUBROUTINE ppm_frand_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_frand_a_3d
    MODULE PROCEDURE irand_a_1d
    MODULE PROCEDURE irand_a_2d
    MODULE PROCEDURE irand_a_3d
    MODULE PROCEDURE irand8_a_1d
    MODULE PROCEDURE irand8_a_2d
    MODULE PROCEDURE irand8_a_3d
    MODULE PROCEDURE drand_a_1d
    MODULE PROCEDURE drand_a_2d
    MODULE PROCEDURE drand_a_3d
    MODULE PROCEDURE frand_a_1d
    MODULE PROCEDURE frand_a_2d
    MODULE PROCEDURE frand_a_3d
  END INTERFACE a_rand

  !> unfortunately, Fortrans random number generator is only prepared
  !! to produce REAL-type results, this add similar capabilities for
  !! INTEGER results in range [0,irand_max]
  !! and REALs in the range [0.0,1.0)
  INTERFACE a_randp
    SUBROUTINE ppm_irandp_a(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_irandp_a
    SUBROUTINE ppm_irandp_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_irandp_a_2d
    SUBROUTINE ppm_irandp_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_irandp_a_3d
    SUBROUTINE ppm_irandp8_a(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER(i8), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_irandp8_a
    SUBROUTINE ppm_irandp8_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER, INTENT(in) :: n
      INTEGER(i8), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_irandp8_a_2d
    SUBROUTINE ppm_irandp8_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER, INTENT(in) :: n
      INTEGER(i8), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_irandp8_a_3d
    SUBROUTINE ppm_drandp_a(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_drandp_a
    SUBROUTINE ppm_drandp_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_drandp_a_2d
    SUBROUTINE ppm_drandp_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_drandp_a_3d
    SUBROUTINE ppm_frandp_a(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_frandp_a
    SUBROUTINE ppm_frandp_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_frandp_a_2d
    SUBROUTINE ppm_frandp_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_frandp_a_3d
    MODULE PROCEDURE irandp_a_1d
    MODULE PROCEDURE irandp_a_2d
    MODULE PROCEDURE irandp_a_3d
    MODULE PROCEDURE drandp_a_1d
    MODULE PROCEDURE drandp_a_2d
    MODULE PROCEDURE drandp_a_3d
    MODULE PROCEDURE frandp_a_1d
    MODULE PROCEDURE frandp_a_2d
    MODULE PROCEDURE frandp_a_3d
  END INTERFACE a_randp

  !> these functions must be called by openmp teams (if compiled with OpenMP)
  !! a is filled with INTEGERs in range [irand_min,irand_max]
  !! or REALs of range (-1.0,1.0)
  INTERFACE a_rand_mt
    SUBROUTINE ppm_irand_mt_a(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_irand_mt_a
    SUBROUTINE ppm_irand_mt_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_irand_mt_a_2d
    SUBROUTINE ppm_irand_mt_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_irand_mt_a_3d
    SUBROUTINE ppm_irand8_mt_a(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER(i8), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_irand8_mt_a
    SUBROUTINE ppm_irand8_mt_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER, INTENT(in) :: n
      INTEGER(i8), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_irand8_mt_a_2d
    SUBROUTINE ppm_irand8_mt_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: i8
      INTEGER, INTENT(in) :: n
      INTEGER(i8), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_irand8_mt_a_3d
    SUBROUTINE ppm_drand_mt_a(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_drand_mt_a
    SUBROUTINE ppm_drand_mt_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_drand_mt_a_2d
    SUBROUTINE ppm_drand_mt_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_drand_mt_a_3d
    SUBROUTINE ppm_frand_mt_a(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_frand_mt_a
    SUBROUTINE ppm_frand_mt_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_frand_mt_a_2d
    SUBROUTINE ppm_frand_mt_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_frand_mt_a_3d
    MODULE PROCEDURE irand_mt_a_1d
    MODULE PROCEDURE irand_mt_a_2d
    MODULE PROCEDURE irand_mt_a_3d
    MODULE PROCEDURE irand8_mt_a_1d
    MODULE PROCEDURE irand8_mt_a_2d
    MODULE PROCEDURE irand8_mt_a_3d
    MODULE PROCEDURE drand_mt_a_1d
    MODULE PROCEDURE drand_mt_a_2d
    MODULE PROCEDURE drand_mt_a_3d
    MODULE PROCEDURE frand_mt_a_1d
    MODULE PROCEDURE frand_mt_a_2d
    MODULE PROCEDURE frand_mt_a_3d
  END INTERFACE a_rand_mt

  !> generate arrays filled with positive random numbers
  !! these functions must be called by all threads in an OpenMP team!
  !! a is filled with INTEGERs in range [0,irand_max]
  !! or REALs of range [0.0,1.0)
  INTERFACE a_randp_mt
    SUBROUTINE ppm_irandp_mt_a(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_irandp_mt_a
    SUBROUTINE ppm_irandp_mt_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_irandp_mt_a_2d
    SUBROUTINE ppm_irandp_mt_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_irandp_mt_a_3d
    SUBROUTINE ppm_drandp_mt_a(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_drandp_mt_a
    SUBROUTINE ppm_drandp_mt_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_drandp_mt_a_2d
    SUBROUTINE ppm_drandp_mt_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_drandp_mt_a_3d
    SUBROUTINE ppm_frandp_mt_a(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
    END SUBROUTINE ppm_frandp_mt_a
    SUBROUTINE ppm_frandp_mt_a_2d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,*)
    END SUBROUTINE ppm_frandp_mt_a_2d
    SUBROUTINE ppm_frandp_mt_a_3d(a, n)
      USE ppm_std_type_kinds, ONLY: sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,1,*)
    END SUBROUTINE ppm_frandp_mt_a_3d
    MODULE PROCEDURE irandp_mt_a_1d
    MODULE PROCEDURE irandp_mt_a_2d
    MODULE PROCEDURE irandp_mt_a_3d
    MODULE PROCEDURE drandp_mt_a_1d
    MODULE PROCEDURE drandp_mt_a_2d
    MODULE PROCEDURE drandp_mt_a_3d
    MODULE PROCEDURE frandp_mt_a_1d
    MODULE PROCEDURE frandp_mt_a_2d
    MODULE PROCEDURE frandp_mt_a_3d
  END INTERFACE a_randp_mt

  !> generate arrays filled with random numbers in given range
  INTERFACE a_randr
    SUBROUTINE ppm_irandr_a(a, n, range)
      USE ppm_std_type_kinds, ONLY: i4
      USE ppm_extents, ONLY: iinterval
      INTEGER(i4), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
      TYPE(iinterval), INTENT(in) :: range
    END SUBROUTINE ppm_irandr_a
    SUBROUTINE ppm_irandr_a_2d(a, n, range)
      USE ppm_std_type_kinds, ONLY: i4
      USE ppm_extents, ONLY: iinterval
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,*)
      TYPE(iinterval), INTENT(in) :: range
    END SUBROUTINE ppm_irandr_a_2d
    SUBROUTINE ppm_irandr_a_3d(a, n, range)
      USE ppm_std_type_kinds, ONLY: i4
      USE ppm_extents, ONLY: iinterval
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n,1,*)
      TYPE(iinterval), INTENT(in) :: range
    END SUBROUTINE ppm_irandr_a_3d
    SUBROUTINE ppm_drandr_a(a, n, range)
      USE ppm_std_type_kinds, ONLY: dp
      USE ppm_extents, ONLY: iinterval_dp
      REAL(dp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
      TYPE(iinterval_dp), INTENT(in) :: range
    END SUBROUTINE ppm_drandr_a
    SUBROUTINE ppm_drandr_a_2d(a, n, range)
      USE ppm_std_type_kinds, ONLY: dp
      USE ppm_extents, ONLY: iinterval_dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,*)
      TYPE(iinterval_dp), INTENT(in) :: range
    END SUBROUTINE ppm_drandr_a_2d
    SUBROUTINE ppm_drandr_a_3d(a, n, range)
      USE ppm_std_type_kinds, ONLY: dp
      USE ppm_extents, ONLY: iinterval_dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n,1,*)
      TYPE(iinterval_dp), INTENT(in) :: range
    END SUBROUTINE ppm_drandr_a_3d
    SUBROUTINE ppm_frandr_a(a, n, range)
      USE ppm_std_type_kinds, ONLY: sp
      USE ppm_extents, ONLY: iinterval_sp
      REAL(sp), INTENT(out) :: a(*)
      INTEGER, INTENT(in) :: n
      TYPE(iinterval_sp), INTENT(in) :: range
    END SUBROUTINE ppm_frandr_a
    SUBROUTINE ppm_frandr_a_2d(a, n, range)
      USE ppm_std_type_kinds, ONLY: sp
      USE ppm_extents, ONLY: iinterval_sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,*)
      TYPE(iinterval_sp), INTENT(in) :: range
    END SUBROUTINE ppm_frandr_a_2d
    SUBROUTINE ppm_frandr_a_3d(a, n, range)
      USE ppm_std_type_kinds, ONLY: sp
      USE ppm_extents, ONLY: iinterval_sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n,1,*)
      TYPE(iinterval_sp), INTENT(in) :: range
    END SUBROUTINE ppm_frandr_a_3d
    MODULE PROCEDURE irandr_a_1d
    MODULE PROCEDURE irandr_a_2d
    MODULE PROCEDURE irandr_a_3d
    MODULE PROCEDURE drandr_a_1d
    MODULE PROCEDURE drandr_a_2d
    MODULE PROCEDURE drandr_a_3d
    MODULE PROCEDURE frandr_a_1d
    MODULE PROCEDURE frandr_a_2d
    MODULE PROCEDURE frandr_a_3d
  END INTERFACE a_randr

  !> generate arrays filled with random numbers in given range
  !! these functions must be called by all threads in an OpenMP team!
  INTERFACE a_randr_mt
    SUBROUTINE ppm_irandr_mt_a(a, n, range)
      USE ppm_std_type_kinds, ONLY: i4
      USE ppm_extents, ONLY: iinterval
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(n)
      TYPE(iinterval), INTENT(in) :: range
    END SUBROUTINE ppm_irandr_mt_a
    SUBROUTINE ppm_irandr_mt_a_2d(a, n, range)
      USE ppm_std_type_kinds, ONLY: i4
      USE ppm_extents, ONLY: iinterval
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(1,n)
      TYPE(iinterval), INTENT(in) :: range
    END SUBROUTINE ppm_irandr_mt_a_2d
    SUBROUTINE ppm_irandr_mt_a_3d(a, n, range)
      USE ppm_std_type_kinds, ONLY: i4
      USE ppm_extents, ONLY: iinterval
      INTEGER, INTENT(in) :: n
      INTEGER(i4), INTENT(out) :: a(1,1,n)
      TYPE(iinterval), INTENT(in) :: range
    END SUBROUTINE ppm_irandr_mt_a_3d
    SUBROUTINE ppm_drandr_mt_a(a, n, range)
      USE ppm_std_type_kinds, ONLY: dp
      USE ppm_extents, ONLY: iinterval_dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(n)
      TYPE(iinterval_dp), INTENT(in) :: range
    END SUBROUTINE ppm_drandr_mt_a
    SUBROUTINE ppm_drandr_mt_a_2d(a, n, range)
      USE ppm_std_type_kinds, ONLY: dp
      USE ppm_extents, ONLY: iinterval_dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(1,n)
      TYPE(iinterval_dp), INTENT(in) :: range
    END SUBROUTINE ppm_drandr_mt_a_2d
    SUBROUTINE ppm_drandr_mt_a_3d(a, n, range)
      USE ppm_std_type_kinds, ONLY: dp
      USE ppm_extents, ONLY: iinterval_dp
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(out) :: a(1,1,n)
      TYPE(iinterval_dp), INTENT(in) :: range
    END SUBROUTINE ppm_drandr_mt_a_3d
    SUBROUTINE ppm_frandr_mt_a(a, n, range)
      USE ppm_std_type_kinds, ONLY: sp
      USE ppm_extents, ONLY: iinterval_sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(n)
      TYPE(iinterval_sp), INTENT(in) :: range
    END SUBROUTINE ppm_frandr_mt_a
    SUBROUTINE ppm_frandr_mt_a_2d(a, n, range)
      USE ppm_std_type_kinds, ONLY: sp
      USE ppm_extents, ONLY: iinterval_sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(1,n)
      TYPE(iinterval_sp), INTENT(in) :: range
    END SUBROUTINE ppm_frandr_mt_a_2d
    SUBROUTINE ppm_frandr_mt_a_3d(a, n, range)
      USE ppm_std_type_kinds, ONLY: sp
      USE ppm_extents, ONLY: iinterval_sp
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(out) :: a(1,1,n)
      TYPE(iinterval_sp), INTENT(in) :: range
    END SUBROUTINE ppm_frandr_mt_a_3d
    MODULE PROCEDURE irandr_mt_a_1d
    MODULE PROCEDURE irandr_mt_a_2d
    MODULE PROCEDURE irandr_mt_a_3d
    MODULE PROCEDURE drandr_mt_a_1d
    MODULE PROCEDURE drandr_mt_a_2d
    MODULE PROCEDURE drandr_mt_a_3d
    MODULE PROCEDURE frandr_mt_a_1d
    MODULE PROCEDURE frandr_mt_a_2d
    MODULE PROCEDURE frandr_mt_a_3d
  END INTERFACE a_randr_mt

  INTERFACE
    SUBROUTINE ppm_initialize_irand(comm, random_seed)
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(inout) :: random_seed
    END SUBROUTINE ppm_initialize_irand

    SUBROUTINE ppm_finalize_irand()
    END SUBROUTINE ppm_finalize_irand
  END INTERFACE
  PUBLIC :: ppm_irand, ppm_irandp, ppm_drand, ppm_drandp, ppm_frand, ppm_frandp
  PUBLIC :: ppm_irandr, ppm_drandr, ppm_frandr
  PUBLIC :: initialize_irand, finalize_irand
  PUBLIC :: irand_min, irand_max
  PUBLIC :: ppm_irand8, irand8_min, irand8_max, ppm_irandp8, ppm_irandr8
  PUBLIC :: a_rand, a_randp, a_randr, a_rand_mt, a_randp_mt, a_randr_mt
CONTAINS
  SUBROUTINE irand_a_1d(a)
    INTEGER(i4), INTENT(out) :: a(:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE irand_a_1d

  SUBROUTINE irand_a_2d(a)
    INTEGER(i4), INTENT(out) :: a(:,:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE irand_a_2d

  SUBROUTINE irand_a_3d(a)
    INTEGER(i4), INTENT(out) :: a(:,:,:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE irand_a_3d

  SUBROUTINE irand8_a_1d(a)
    INTEGER(i8), INTENT(out) :: a(:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE irand8_a_1d

  SUBROUTINE irand8_a_2d(a)
    INTEGER(i8), INTENT(out) :: a(:,:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE irand8_a_2d

  SUBROUTINE irand8_a_3d(a)
    INTEGER(i8), INTENT(out) :: a(:,:,:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE irand8_a_3d

  SUBROUTINE drand_a_1d(a)
    REAL(dp), INTENT(out) :: a(:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE drand_a_1d

  SUBROUTINE drand_a_2d(a)
    REAL(dp), INTENT(out) :: a(:, :)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE drand_a_2d

  SUBROUTINE drand_a_3d(a)
    REAL(dp), INTENT(out) :: a(:, :, :)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE drand_a_3d

  SUBROUTINE frand_a_1d(a)
    REAL(sp), INTENT(out) :: a(:)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE frand_a_1d

  SUBROUTINE frand_a_2d(a)
    REAL(sp), INTENT(out) :: a(:, :)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE frand_a_2d

  SUBROUTINE frand_a_3d(a)
    REAL(sp), INTENT(out) :: a(:, :, :)
    CALL a_rand(a, SIZE(a))
  END SUBROUTINE frand_a_3d

  SUBROUTINE irandp_a_1d(a)
    INTEGER(i4), INTENT(out) :: a(:)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE irandp_a_1d

  SUBROUTINE irandp_a_2d(a)
    INTEGER(i4), INTENT(out) :: a(:, :)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE irandp_a_2d

  SUBROUTINE irandp_a_3d(a)
    INTEGER(i4), INTENT(out) :: a(:, :, :)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE irandp_a_3d

  SUBROUTINE drandp_a_1d(a)
    REAL(dp), INTENT(out) :: a(:)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE drandp_a_1d

  SUBROUTINE drandp_a_2d(a)
    REAL(dp), INTENT(out) :: a(:, :)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE drandp_a_2d

  SUBROUTINE drandp_a_3d(a)
    REAL(dp), INTENT(out) :: a(:, :, :)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE drandp_a_3d

  SUBROUTINE frandp_a_1d(a)
    REAL(sp), INTENT(out) :: a(:)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE frandp_a_1d

  SUBROUTINE frandp_a_2d(a)
    REAL(sp), INTENT(out) :: a(:, :)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE frandp_a_2d

  SUBROUTINE frandp_a_3d(a)
    REAL(sp), INTENT(out) :: a(:, :, :)
    CALL a_randp(a, SIZE(a))
  END SUBROUTINE frandp_a_3d

  SUBROUTINE irand_mt_a_1d(a)
    INTEGER(i4), INTENT(out) :: a(:)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE irand_mt_a_1d

  SUBROUTINE irand_mt_a_2d(a)
    INTEGER(i4), INTENT(out) :: a(:, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE irand_mt_a_2d

  SUBROUTINE irand_mt_a_3d(a)
    INTEGER(i4), INTENT(out) :: a(:, :, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE irand_mt_a_3d

  SUBROUTINE irand8_mt_a_1d(a)
    INTEGER(i8), INTENT(out) :: a(:)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE irand8_mt_a_1d

  SUBROUTINE irand8_mt_a_2d(a)
    INTEGER(i8), INTENT(out) :: a(:, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE irand8_mt_a_2d

  SUBROUTINE irand8_mt_a_3d(a)
    INTEGER(i8), INTENT(out) :: a(:, :, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE irand8_mt_a_3d

  SUBROUTINE drand_mt_a_1d(a)
    REAL(dp), INTENT(out) :: a(:)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE drand_mt_a_1d

  SUBROUTINE drand_mt_a_2d(a)
    REAL(dp), INTENT(out) :: a(:, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE drand_mt_a_2d

  SUBROUTINE drand_mt_a_3d(a)
    REAL(dp), INTENT(out) :: a(:, :, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE drand_mt_a_3d

  SUBROUTINE frand_mt_a_1d(a)
    REAL(sp), INTENT(out) :: a(:)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE frand_mt_a_1d

  SUBROUTINE frand_mt_a_2d(a)
    REAL(sp), INTENT(out) :: a(:, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE frand_mt_a_2d

  SUBROUTINE frand_mt_a_3d(a)
    REAL(sp), INTENT(out) :: a(:, :, :)
    CALL a_rand_mt(a, SIZE(a))
  END SUBROUTINE frand_mt_a_3d

  SUBROUTINE irandp_mt_a_1d(a)
    INTEGER(i4), INTENT(out) :: a(:)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE irandp_mt_a_1d

  SUBROUTINE irandp_mt_a_2d(a)
    INTEGER(i4), INTENT(out) :: a(:, :)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE irandp_mt_a_2d

  SUBROUTINE irandp_mt_a_3d(a)
    INTEGER(i4), INTENT(out) :: a(:, :, :)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE irandp_mt_a_3d

  SUBROUTINE drandp_mt_a_1d(a)
    REAL(dp), INTENT(out) :: a(:)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE drandp_mt_a_1d

  SUBROUTINE drandp_mt_a_2d(a)
    REAL(dp), INTENT(out) :: a(:, :)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE drandp_mt_a_2d

  SUBROUTINE drandp_mt_a_3d(a)
    REAL(dp), INTENT(out) :: a(:, :, :)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE drandp_mt_a_3d

  SUBROUTINE frandp_mt_a_1d(a)
    REAL(sp), INTENT(out) :: a(:)
    CALL a_randp_mt(a(1:SIZE(a)), SIZE(a))
  END SUBROUTINE frandp_mt_a_1d

  SUBROUTINE frandp_mt_a_2d(a)
    REAL(sp), INTENT(out) :: a(:, :)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE frandp_mt_a_2d

  SUBROUTINE frandp_mt_a_3d(a)
    REAL(sp), INTENT(out) :: a(:, :, :)
    CALL a_randp_mt(a, SIZE(a))
  END SUBROUTINE frandp_mt_a_3d

  SUBROUTINE irandr_a_1d(a, range)
    INTEGER(i4), INTENT(out) :: a(:)
    TYPE(iinterval), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE irandr_a_1d

  SUBROUTINE irandr_a_2d(a, range)
    INTEGER(i4), INTENT(out) :: a(:, :)
    TYPE(iinterval), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE irandr_a_2d

  SUBROUTINE irandr_a_3d(a, range)
    INTEGER(i4), INTENT(out) :: a(:, :, :)
    TYPE(iinterval), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE irandr_a_3d

  SUBROUTINE drandr_a_1d(a, range)
    REAL(dp), INTENT(out) :: a(:)
    TYPE(iinterval_dp), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE drandr_a_1d

  SUBROUTINE drandr_a_2d(a, range)
    REAL(dp), INTENT(out) :: a(:, :)
    TYPE(iinterval_dp), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE drandr_a_2d

  SUBROUTINE drandr_a_3d(a, range)
    REAL(dp), INTENT(out) :: a(:, :, :)
    TYPE(iinterval_dp), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE drandr_a_3d

  SUBROUTINE frandr_a_1d(a, range)
    REAL(sp), INTENT(out) :: a(:)
    TYPE(iinterval_sp), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE frandr_a_1d

  SUBROUTINE frandr_a_2d(a, range)
    REAL(sp), INTENT(out) :: a(:, :)
    TYPE(iinterval_sp), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE frandr_a_2d

  SUBROUTINE frandr_a_3d(a, range)
    REAL(sp), INTENT(out) :: a(:, :, :)
    TYPE(iinterval_sp), INTENT(in) :: range
    CALL a_randr(a, SIZE(a), range)
  END SUBROUTINE frandr_a_3d

  SUBROUTINE irandr_mt_a_1d(a, range)
    INTEGER(i4), INTENT(out) :: a(:)
    TYPE(iinterval), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE irandr_mt_a_1d

  SUBROUTINE irandr_mt_a_2d(a, range)
    INTEGER(i4), INTENT(out) :: a(:, :)
    TYPE(iinterval), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE irandr_mt_a_2d

  SUBROUTINE irandr_mt_a_3d(a, range)
    INTEGER(i4), INTENT(out) :: a(:, :, :)
    TYPE(iinterval), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE irandr_mt_a_3d

  SUBROUTINE drandr_mt_a_1d(a, range)
    REAL(dp), INTENT(out) :: a(:)
    TYPE(iinterval_dp), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE drandr_mt_a_1d

  SUBROUTINE drandr_mt_a_2d(a, range)
    REAL(dp), INTENT(out) :: a(:, :)
    TYPE(iinterval_dp), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE drandr_mt_a_2d

  SUBROUTINE drandr_mt_a_3d(a, range)
    REAL(dp), INTENT(out) :: a(:, :, :)
    TYPE(iinterval_dp), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE drandr_mt_a_3d

  SUBROUTINE frandr_mt_a_1d(a, range)
    REAL(sp), INTENT(out) :: a(:)
    TYPE(iinterval_sp), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE frandr_mt_a_1d

  SUBROUTINE frandr_mt_a_2d(a, range)
    REAL(sp), INTENT(out) :: a(:, :)
    TYPE(iinterval_sp), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE frandr_mt_a_2d

  SUBROUTINE frandr_mt_a_3d(a, range)
    REAL(sp), INTENT(out) :: a(:, :, :)
    TYPE(iinterval_sp), INTENT(in) :: range
    CALL a_randr_mt(a, SIZE(a), range)
  END SUBROUTINE frandr_mt_a_3d

  SUBROUTINE initialize_irand(comm, random_seed, seed_output)
    INTEGER, INTENT(inout) :: comm
    INTEGER, OPTIONAL, INTENT(in) :: random_seed
    INTEGER, OPTIONAL, INTENT(out) :: seed_output

    INTEGER :: temp

    IF (PRESENT(random_seed)) THEN; temp = random_seed; ELSE



      temp = 9

    END IF

    CALL ppm_initialize_irand(comm, temp)
    IF (PRESENT(seed_output)) seed_output = temp
  END SUBROUTINE initialize_irand

  SUBROUTINE finalize_irand
    CALL ppm_finalize_irand
  END SUBROUTINE finalize_irand
END MODULE ppm_irand_internal
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
