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
#ifdef __xlC__
@PROCESS STRICT(precision)
#endif
#include "fc_feature_defs.inc"
MODULE ppm_math_extensions
#ifdef USE_MPI
  USE ppm_base, ONLY: abort_ppm, ppm_default_comm, &
       calls_to_mpi_are_allowed
  USE ppm_std_type_kinds_mp, ONLY: mp_dp
#endif
  USE ppm_std_type_kinds, ONLY: dp, sp
  USE ppm_math_extensions_internal, ONLY: fpu_save_cw => ppm_fpu_save_cw, &
       fpu_restore_cw => ppm_fpu_restore_cw, &
       fpu_set_precision, fpu_precision_sp, fpu_precision_dp, fpu_precision_ep,&
       fpu_restore_mxcsr => ppm_fpu_restore_mxcsr, &
       fpu_save_mxcsr => ppm_fpu_save_mxcsr, &
       fpu_set_abrupt_underflow
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  REAL(dp), PARAMETER, PUBLIC :: &
       m_pi_dp=3.14159265358979323846_dp
  REAL(sp), PARAMETER, PUBLIC :: &
       m_pi_sp=3.14159265358979323846_sp

#include "ppm_real_sp_dp_edit_descriptor.inc"
  !> data edit descriptors for real kinds sp and dp
  CHARACTER(len=*), PARAMETER :: &
       de_g_sp = PPM_de_g_sp, &
       de_g_dp = PPM_de_g_dp
  !> character width needed for the corresponding data edit descriptor
  INTEGER, PARAMETER :: de_g_sp_width = PPM_de_g_sp_width, &
       de_g_dp_width = PPM_de_g_dp_width

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
#ifdef USE_MPI
  INTEGER :: ddpdd_op
  INTERFACE ddp_sum_mp
    MODULE PROCEDURE ddp_sum_mp_dp_1d_es
    MODULE PROCEDURE ddp_sum_mp_dp_1d
    MODULE PROCEDURE ddp_sum_mp_dp_2d
    MODULE PROCEDURE ddp_sum_mp_dp_3d
  END INTERFACE ddp_sum_mp
  PUBLIC :: ddp_sum_mp
#endif
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

#ifdef USE_MPI
  SUBROUTINE ddpdd(a, b, len, itype)
    INTEGER, INTENT(in) :: len, itype
    COMPLEX(dp), INTENT(in) :: a(len)
    COMPLEX(dp), INTENT(inout) :: b(len)
    INTEGER :: i
    SELECT CASE (itype)
    CASE(mpi_double_complex)
      DO i = 1, len
        b(i) = ddp_add(a(i), b(i))
      END DO
    CASE default
      CALL abort_ppm("invalid MPI data type in double-double summation", &
           filename, __LINE__)
    END SELECT
  END SUBROUTINE ddpdd

  !> compute double-double-precision corrected sum of distributed 1d-array
  !!
  !!
  !! explicit shape function to consolidate common code of multi-dimensional
  !! versions below.
  !!
  !! This is a collective MPI-Call.
  !! @param n size of array a
  !! @param a array a(1)..a(n) to sum up
  !! @return \f$\sum^n_{i=1} a_i\f$
  FUNCTION ddp_sum_mp_dp_1d_es(a, n, comm) RESULT(c)
    COMPLEX(dp) :: c
    INTEGER, INTENT(in) :: comm, n
    REAL(dp), INTENT(in) :: a(n)
    INTEGER :: ierror
    c = ddp_sum(a)
    CALL mpi_allreduce(mpi_in_place, c, 1, mpi_double_complex, &
         ddpdd_op, comm, ierror)
    CALL handle_mpi_error(ierror, comm, filename, __LINE__)
  END FUNCTION ddp_sum_mp_dp_1d_es

  !> compute double-double-precision corrected sum of distributed 1d-array
  !!
  !! This is a collective MPI-Call.
  !! @param a array a(1)..a(n) to sum up
  !! @return \f$\sum^n_{i=1} a_i\f$
  FUNCTION ddp_sum_mp_dp_1d(a, comm) RESULT(c)
    COMPLEX(dp) :: c
    INTEGER, INTENT(in) :: comm
    REAL(dp), INTENT(in) :: a(:)
    INTEGER :: n
    n = SIZE(a)
    c = ddp_sum_mp_dp_1d_es(a, n, comm)
  END FUNCTION ddp_sum_mp_dp_1d

  !> compute double-double-precision corrected sum of distributed 2d-array
  !> @see ddp_sum_mp_dp_1d
  !> @param a array a(1,1)..a(m,n) to sum up
  !> @return \f$\sum^m_{i=1}\sum^n_{j=1} a_{i,j}\f$
  FUNCTION ddp_sum_mp_dp_2d(a, comm) RESULT(c)
    COMPLEX(dp) :: c
    REAL(dp), INTENT(in) :: a(:, :)
    INTEGER, INTENT(in) :: comm
    INTEGER :: n
    n = SIZE(a)
    c = ddp_sum_mp_dp_1d_es(a, n, comm)
  END FUNCTION ddp_sum_mp_dp_2d

  !> compute double-double-precision corrected sum of distributed 3d-array
  !> @see ddp_sum_mp_dp_1d
  !> @param a array a(1,1,1)..a(m,n,o) to sum up
  !> @return \f$\sum^m_{i=1}\sum^n_{j=1}\sum^o_{l=1} a_{i,j,l}\f$
  FUNCTION ddp_sum_mp_dp_3d(a, comm) RESULT(c)
    COMPLEX(dp) :: c
    REAL(dp), INTENT(in) :: a(:, :, :)
    INTEGER, INTENT(in) :: comm
    INTEGER :: n
    n = SIZE(a)
    c = ddp_sum_mp_dp_1d_es(a, n, comm)
  END FUNCTION ddp_sum_mp_dp_3d
#endif

  SUBROUTINE initialize_math_extensions
#ifdef USE_MPI
    INTEGER :: ierror
#endif
#ifdef USE_MPI
    IF (calls_to_mpi_are_allowed()) THEN
      CALL mpi_op_create(ddpdd, .TRUE., ddpdd_op, ierror)
      CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    END IF
#endif

  END SUBROUTINE initialize_math_extensions

  SUBROUTINE finalize_math_extensions
#ifdef USE_MPI
    INTEGER :: ierror
    CALL mpi_op_free(ddpdd_op, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
#endif
  END SUBROUTINE finalize_math_extensions

#ifdef USE_MPI
  SUBROUTINE handle_mpi_error(ierror, comm, source, line)
    INTEGER, INTENT(in) :: ierror, comm, line
    CHARACTER, INTENT(in) :: source

    INTEGER :: msg_len, ierror_
    CHARACTER(len=mpi_max_error_string) :: msg

    IF (ierror /= MPI_SUCCESS) THEN
      CALL mpi_error_string(ierror, msg, msg_len, ierror_)
      CALL abort_ppm(msg(1:msg_len), source, line, comm)
    END IF
  END SUBROUTINE handle_mpi_error
#endif
END MODULE ppm_math_extensions
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
