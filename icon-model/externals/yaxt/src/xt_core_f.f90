!>
!! @file xt_core_f.f90
!! @brief Fortran interface to yaxt core declarations
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
MODULE xt_core
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int
#ifdef XT_INT_FC_KIND_IN_ISO_C_BINDING
  USE, INTRINSIC :: iso_c_binding, ONLY: XT_INT_FC_KIND
#endif
  USE xt_mpi, ONLY: XT_INT_FC_MPIDT, xt_mpi_fint_kind
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC, PARAMETER :: xt_int_kind   = XT_INT_FC_KIND
  INTEGER, PUBLIC, PARAMETER :: pi2 = 4
  INTEGER, PUBLIC, PARAMETER :: pi4 = 9
  INTEGER, PUBLIC, PARAMETER :: pi8 = 14
  INTEGER, PUBLIC, PARAMETER :: i2 = SELECTED_INT_KIND(pi2)
  INTEGER, PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)
  PUBLIC :: xt_initialize, xt_finalize, xt_abort, xt_get_default_comm, char
  PUBLIC :: xt_initialized, xt_finalized
  PUBLIC :: xt_slice_c_loc
  PUBLIC :: OPERATOR(==), OPERATOR(/=)

  PUBLIC :: xt_mpi_fint_kind
  INTEGER, PUBLIC, PARAMETER :: xt_int_mpidt = XT_INT_FC_MPIDT
  INTEGER(xt_int_kind), PARAMETER :: dummy = 0_xt_int_kind
  !> number of decimal places needed to print any variable of type
  !! INTEGER(xt_int_kind)
  INTEGER, PUBLIC, PARAMETER :: xt_int_dec_len &
       = CEILING(1.0 + REAL(DIGITS(dummy)) * LOG10(REAL(RADIX(dummy))))
  CHARACTER(9), PARAMETER :: xt_stripe_tag = 'xt_stripe'
  !> maximal length of string xt_stripe(a, b, c)
  INTEGER, PUBLIC, PARAMETER :: xt_stripe2s_len &
       = LEN(xt_stripe_tag) + 2 + 4 + 3 * xt_int_dec_len

  TYPE, BIND(C), PUBLIC :: xt_stripe
    INTEGER(xt_int_kind) :: start
    INTEGER(xt_int_kind) :: stride
    INTEGER(c_int) :: nstrides
  END TYPE xt_stripe

  TYPE, BIND(C), PUBLIC :: xt_bounds
    INTEGER(xt_int_kind) :: start, size
  END TYPE xt_bounds

  !> describes range of positions starting with start up to start + size - 1
  !! i.e. [start,start+size) if size is positive and down to start + size + 1
  !! i.e. (start+size,start] if size is negative
  TYPE, BIND(c), PUBLIC :: xt_pos_ext
    INTEGER(c_int) :: start, size
  END TYPE xt_pos_ext

  INTERFACE

    FUNCTION xt_get_default_comm() RESULT(comm) &
         BIND(c, name='xt_get_default_comm_f')
      IMPORT :: xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind) :: comm
    END FUNCTION xt_get_default_comm

    SUBROUTINE xt_initialize(default_comm) BIND(C, name='xt_initialize_f')
      IMPORT:: xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind), INTENT(in) :: default_comm
    END SUBROUTINE xt_initialize

    SUBROUTINE xt_finalize() BIND(C, name='xt_finalize')
    END SUBROUTINE xt_finalize

    SUBROUTINE xt_restore_default_abort_hndl
    END SUBROUTINE xt_restore_default_abort_hndl

  END INTERFACE

  INTERFACE xt_abort
    MODULE PROCEDURE xt_abort4
    MODULE PROCEDURE xt_abort3
  END INTERFACE xt_abort

  INTERFACE
    SUBROUTINE xt_abort_c(comm, msg, source, line) BIND(c, name='xt_abort_f')
      IMPORT :: c_char, xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in):: comm
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: msg
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: source
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: line
    END SUBROUTINE xt_abort_c
  END INTERFACE

  INTERFACE char
    MODULE PROCEDURE xt_stripe2char
  END INTERFACE char

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE xt_pos_ext_eq
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE xt_pos_ext_ne
  END INTERFACE OPERATOR(/=)

  EXTERNAL :: xt_slice_c_loc

  PUBLIC :: set_abort_handler, xt_restore_default_abort_hndl

  ENUM, BIND( C )
    ENUMERATOR :: xt_lib_pre_init, &
         xt_lib_initialized, &
         xt_lib_finalized
  END ENUM
  INTEGER(c_int), PUBLIC, BIND(c, name='xt_lib_state') :: xt_lib_state

CONTAINS

  SUBROUTINE xt_abort4(comm, msg, source, line)
    INTEGER, INTENT(in) :: comm
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_c(comm, TRIM(msg)//c_null_char, &
         TRIM(source)//c_null_char, line)
  END SUBROUTINE xt_abort4

  SUBROUTINE xt_abort3(msg, source, line)
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_c(xt_get_default_comm(), TRIM(msg)//c_null_char, &
         TRIM(source)//c_null_char, line)
  END SUBROUTINE xt_abort3

  ELEMENTAL FUNCTION xt_stripe2char(stripe) RESULT(str)
    CHARACTER(len=xt_stripe2s_len) :: str
    TYPE(xt_stripe), INTENT(in) :: stripe
    WRITE (str, '(2a,3(i0,a))') xt_stripe_tag, '(', stripe%start, ', ', &
         stripe%stride, ', ', stripe%nstrides, ')'
  END FUNCTION xt_stripe2char

  PURE FUNCTION xt_initialized() RESULT(is_initialized)
    LOGICAL :: is_initialized
    is_initialized = xt_lib_state > xt_lib_pre_init
  END FUNCTION xt_initialized

  PURE FUNCTION xt_finalized() RESULT(is_finalized)
    LOGICAL :: is_finalized
    is_finalized = xt_lib_state == xt_lib_finalized
  END FUNCTION xt_finalized

  ELEMENTAL FUNCTION xt_pos_ext_eq(a, b) RESULT(p)
    TYPE(xt_pos_ext), INTENT(in) :: a, b
    LOGICAL :: p
    p = a%start == b%start .AND. (a%size == b%size &
         .OR. (ABS(a%size) == 1 .AND. ABS(a%size) == ABS(b%size)))
  END FUNCTION xt_pos_ext_eq

  ELEMENTAL FUNCTION xt_pos_ext_ne(a, b) RESULT(p)
    TYPE(xt_pos_ext), INTENT(in) :: a, b
    LOGICAL :: p
    p = a%start /= b%start .OR. (a%size /= b%size &
         .AND. .NOT. (ABS(a%size) == 1 .AND. ABS(a%size) == ABS(b%size)))
  END FUNCTION xt_pos_ext_ne

  !> set routine f to use as abort function which is called on xt_abort
  SUBROUTINE set_abort_handler(f)
    INTERFACE
      SUBROUTINE f(comm, msg, source, line)
        INTEGER, INTENT(in) :: comm, line
        CHARACTER(len=*), INTENT(in) :: msg, source
      END SUBROUTINE f
      SUBROUTINE xt_set_abort_handler(f)
        INTERFACE
          SUBROUTINE f(comm, msg, source, line)
            INTEGER, INTENT(in) :: comm, line
            CHARACTER(len=*), INTENT(in) :: msg, source
          END SUBROUTINE f
        END INTERFACE
      END SUBROUTINE xt_set_abort_handler
    END INTERFACE
    CALL xt_set_abort_handler(f)
  END SUBROUTINE set_abort_handler

END MODULE xt_core
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
! license-default: "bsd"
! End:
!
