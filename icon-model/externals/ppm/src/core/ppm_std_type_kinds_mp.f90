!>
!! @file ppm_std_type_kinds_mp.f90
!! @brief map type kinds supported by library to corresponding MPI datatypes
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
#include "fc_feature_defs.inc"
MODULE ppm_std_type_kinds_mp
#ifdef USE_MPI
#include "mpi_fc_conf.inc"
#endif
  USE ppm_base, ONLY: abort_ppm
#ifndef HAVE_MPI_INTEGER4
  USE ppm_std_type_kinds, ONLY: pi4
#endif
#ifndef HAVE_MPI_INTEGER8
  USE ppm_std_type_kinds, ONLY: pi8
#endif
#if ! defined(HAVE_MPI_INTEGER4) || ! defined(HAVE_MPI_INTEGER8)
  USE ppm_std_type_kinds, ONLY: i4, i8
#endif
#ifndef HAVE_MPI_REAL4
  USE ppm_std_type_kinds, ONLY: ps, rs
#endif
#ifndef HAVE_MPI_REAL4
  USE ppm_std_type_kinds, ONLY: pd, rd
#endif
#if ! defined(HAVE_MPI_REAL4) || ! defined(HAVE_MPI_REAL8)
  USE ppm_std_type_kinds, ONLY: sp, dp
#endif
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  !> MPI datatype corresponding to 4 byte integer
#ifdef HAVE_MPI_INTEGER4
  INTEGER, PARAMETER :: mp_i4 = mpi_integer4
#else
  INTEGER, SAVE :: mp_i4
#endif
  !> MPI datatype corresponding to 8 byte integer
#ifdef HAVE_MPI_INTEGER8
  INTEGER, PARAMETER :: mp_i8 = mpi_integer8
#else
  INTEGER, SAVE :: mp_i8
#endif
  !> MPI datatype corresponding to 4 byte float
#ifdef HAVE_MPI_REAL4
  INTEGER, PARAMETER :: mp_sp = mpi_real4
#else
  INTEGER, SAVE :: mp_sp
#endif
  !> MPI datatype corresponding to 8 byte float
#ifdef HAVE_MPI_REAL8
  INTEGER, PARAMETER :: mp_dp = mpi_real8
#else
  INTEGER, SAVE :: mp_dp
#endif
  INTEGER(mpi_address_kind), SAVE :: &
       mp_i4_extent, mp_i8_extent, mp_sp_extent, mp_dp_extent, mp_l_extent
  PUBLIC :: mp_i4, mp_i8, mp_sp, mp_dp
  PUBLIC :: mp_i4_extent, mp_i8_extent, mp_sp_extent, mp_dp_extent, mp_l_extent
  PUBLIC :: create_types_mp, destroy_types_mp

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_std_type_kinds_mp.f90'
CONTAINS
  SUBROUTINE create_types_mp
    INTEGER :: ierror, default_errh, comm_world
    INTEGER(mpi_address_kind) :: lb
    ! save default error handler
    ! note: this does not use ppm_default_comm, because the
    ! mpi_type_create_f90_* routines do not use a specific communicator
    CALL mpi_comm_get_errhandler(mpi_comm_world, default_errh, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_get_errhandler failed", filename, __LINE__)
    ! the use of comm_world is a work-around for old OpenMPI versions
    comm_world = mpi_comm_world
    CALL mpi_comm_set_errhandler(comm_world, mpi_errors_return, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_set_errhandler failed", filename, __LINE__)
    ! initialize C part
    CALL ppm_initialize_std_type_kinds_mp

#ifndef HAVE_MPI_INTEGER4
    CALL mpi_type_create_f90_integer(pi4, mp_i4, ierror)
    IF (ierror /= mpi_success) &
         CALL type_create_fallback_int(pi4, mp_i4, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_create_f90_integer failed", filename, &
         __LINE__)
    CALL mpi_type_commit(mp_i4, ierror)
#endif

#ifndef HAVE_MPI_INTEGER8
    CALL mpi_type_create_f90_integer(pi8, mp_i8, ierror)
    IF (ierror /= mpi_success) &
         CALL type_create_fallback_int(pi8, mp_i8, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_create_f90_integer failed", filename, &
         __LINE__)
    CALL mpi_type_commit(mp_i8, ierror)
#endif

#ifndef HAVE_MPI_REAL4
    CALL mpi_type_create_f90_real(ps, rs, mp_sp, ierror)
    IF (ierror /= mpi_success) &
         CALL type_create_fallback_real(ps, rs, mp_sp, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_create_f90_real failed", filename, __LINE__)
    CALL mpi_type_commit(mp_sp, ierror)
#endif

#ifndef HAVE_MPI_REAL8
    CALL mpi_type_create_f90_real(pd, rd, mp_dp, ierror)
    IF (ierror /= mpi_success) &
         CALL type_create_fallback_real(pd, rd, mp_dp, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_create_f90_real failed", filename, __LINE__)
    CALL mpi_type_commit(mp_dp, ierror)
#endif
    ! query extent of base datatypes
    CALL mpi_type_get_extent(mp_i4, lb, mp_i4_extent, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_get_extent failed", &
         filename, __LINE__)
    CALL mpi_type_get_extent(mp_i8, lb, mp_i8_extent, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_get_extent failed", filename, __LINE__)
    CALL mpi_type_get_extent(mp_sp, lb, mp_sp_extent, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_get_extent failed", filename, __LINE__)
    CALL mpi_type_get_extent(mp_dp, lb, mp_dp_extent, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_get_extent failed", filename, __LINE__)
    CALL mpi_type_get_extent(mpi_logical, lb, mp_l_extent, ierror)
    IF (ierror /= mpi_success) CALL abort_ppm("mpi_type_get_extent failed", &
         filename, __LINE__)

    ! restore default error handler
    CALL mpi_comm_set_errhandler(comm_world, default_errh, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_set_errhandler failed", filename, __LINE__)
#if ! defined(HAVE_MPI_INTEGER4) || ! defined(HAVE_MPI_INTEGER8) || ! defined(HAVE_MPI_REAL4) || ! defined(HAVE_MPI_REAL8)
  CONTAINS
#endif
#if ! defined(HAVE_MPI_INTEGER4) || ! defined(HAVE_MPI_INTEGER8)
    SUBROUTINE type_create_fallback_int(r, dtype, ierror)
      INTEGER, INTENT(in) :: r
      INTEGER, INTENT(out) :: dtype, ierror

      INTEGER, PARAMETER :: dummy=1
      INTEGER, PARAMETER :: integer_byte_size = BIT_SIZE(dummy)/8
      INTEGER :: io_size_integer, io_size, io_bytes

      INTEGER(i4) :: ii4
      INTEGER(i8) :: ii8

      INQUIRE (iolength=io_size_integer) dummy

      SELECT CASE(r)
      CASE(1:9)
        INQUIRE (iolength=io_size) ii4
        io_bytes = io_size / io_size_integer * integer_byte_size
      CASE(10:18)
        INQUIRE (iolength=io_size) ii8
        io_bytes = io_size / io_size_integer * integer_byte_size
      CASE default
        ierror = mpi_err_type
        RETURN
      END SELECT

      CALL mpi_type_match_size(mpi_typeclass_integer, io_bytes, dtype, ierror)
      IF (ierror /= mpi_success) RETURN
    END SUBROUTINE type_create_fallback_int
#endif

#if ! defined(HAVE_MPI_REAL4) || ! defined(HAVE_MPI_REAL8)
    SUBROUTINE type_create_fallback_real(p, r, dtype, ierror)
      INTEGER, INTENT(in) :: r, p
      INTEGER, INTENT(out) :: dtype, ierror

      INTEGER, PARAMETER :: dummy=1
      INTEGER, PARAMETER :: integer_byte_size = BIT_SIZE(dummy)/8
      INTEGER :: io_size_integer, io_size, io_bytes
      INTEGER :: ref_type

      REAL(sp) :: r_sp
      REAL(dp) :: r_dp

      INQUIRE (iolength=io_size_integer) dummy

      SELECT CASE(r)
      CASE(1:6)
        INQUIRE (iolength=io_size) r_sp
        io_bytes = io_size / io_size_integer * integer_byte_size
      CASE(7:15)
        INQUIRE (iolength=io_size) r_dp
        io_bytes = io_size / io_size_integer * integer_byte_size
      CASE default
        ierror = mpi_err_type
        RETURN
      END SELECT

      CALL mpi_type_match_size(mpi_typeclass_real, io_bytes, dtype, ierror)
      IF (ierror /= mpi_success) RETURN
    END SUBROUTINE type_create_fallback_real
#endif
  END SUBROUTINE create_types_mp

  SUBROUTINE destroy_types_mp
#if ! defined(HAVE_MPI_INTEGER4) || ! defined(HAVE_MPI_INTEGER8) || ! defined(HAVE_MPI_REAL4) || ! defined(HAVE_MPI_REAL8)
    INTEGER :: ierror
#endif
#ifndef HAVE_MPI_INTEGER4
    CALL mpi_type_free(mp_i4, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_free failed", filename, __LINE__)
#endif
#ifndef HAVE_MPI_INTEGER8
    CALL mpi_type_free(mp_i8, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_free failed", filename, __LINE__)
#endif
#ifndef HAVE_MPI_REAL4
    CALL mpi_type_free(mp_sp, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_free failed", filename, __LINE__)
#endif
#ifndef HAVE_MPI_REAL8
    CALL mpi_type_free(mp_dp, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_free failed", filename, __LINE__)
#endif
    ! finalize C part
    CALL ppm_finalize_std_type_kinds_mp

  END SUBROUTINE destroy_types_mp

END MODULE ppm_std_type_kinds_mp
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
