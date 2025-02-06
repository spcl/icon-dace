!>
!! @file ppm_std_type_kinds.f90
!! @brief determine fortran 90 type kinds corresponding to
!!                            single/double precision real and 4/8 byte integer
!!
!! @copyright Copyright  (C)  2001  Luis Kornblueh <luis.kornblueh@we dont want spamzmaw.de>
!!
!! @version 1.0
!! @author Luis Kornblueh <luis.kornblueh@we dont want spamzmaw.de>
!
! Keywords:
! Maintainer: Thomas Jahns <jahns@dkrz.de>
!    that is: for the version in scales-ppm, it's originally from ECHAM
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
MODULE ppm_std_type_kinds
  IMPLICIT NONE
  PUBLIC

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15
  !                   exponent  = 37   exponent  =  307
  !
  ! Most likely this are the only possible models.

  ! Floating point section:
  INCLUDE 'ppm_real_sp_dp.inc'

  ! Floating point working precision

  INTEGER, PARAMETER :: wp = dp

  ! Integer section

  INTEGER, PARAMETER :: pi4 = 9
  INTEGER, PARAMETER :: pi8 = 14

  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)

  ! Working precision for index variables
  !
  ! predefined preprocessor macros:
  !
  ! xlf         __64BIT__   checked with P6 and AIX
  ! gfortran    __LP64__    checked with Darwin and Linux
  ! Intel, PGI  __x86_64__  checked with Linux
  ! Sun         __x86_64    checked with Linux

#if defined (__64BIT__) || defined (__LP64__) || defined (__x86_64__) || defined (__x86_64)
  INTEGER, PARAMETER :: widx = i8
#else
  INTEGER, PARAMETER :: widx = i4
#endif
  INCLUDE 'ftype_size.inc'
END MODULE ppm_std_type_kinds
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
