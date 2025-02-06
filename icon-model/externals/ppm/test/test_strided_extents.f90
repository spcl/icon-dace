!>
!! @file test_strided_extents.f90
!! @brief test transfer of array slice
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
#include "fc_feature_defs.inc"
PROGRAM test_strided_extents
  USE ppm_std_type_kinds, ONLY: dp
  USE ppm_extents, ONLY: extent
  USE ppm_strided_extents, ONLY: strided_extent, subarray_mpi_datatype
  USE ppm_base, ONLY: abort_ppm
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  INTEGER, PARAMETER :: dim_max(3) = (/ 10, 12, 8 /), islice=2, jslice=3
  REAL(dp) :: a(dim_max(1), dim_max(2), dim_max(3))
  INTEGER :: i, j, k, slice_mdt, ierror
  TYPE(extent) :: a_shape(3)
  TYPE(strided_extent) :: slice_shape(2)
  INTEGER(mpi_address_kind) :: lb, slice_extent
  INTEGER :: req, status(mpi_status_size, 2), comm_rank, comm_size, dst, src
  CHARACTER(len=*), PARAMETER :: filename = 'test_strided_extents.f90'
  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_init failed", filename, __LINE__)
  DO k = 1, dim_max(3)
    DO j = 1, dim_max(2)
      DO i = 1, dim_max(1)
        a(i, j, k) = REAL(i * 10000 + j * 100 + k, dp)
      END DO
    END DO
  END DO
  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_rank failed", filename, __LINE__)
  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_size failed", filename, __LINE__)
  dst = MOD(comm_rank + 1, comm_size)
  src = MOD(comm_rank - 1 + comm_size, comm_size)
  DO i = 1, 3
    a_shape(i) = extent(1, dim_max(i))
  END DO
  slice_shape(1)%ext = extent(1, islice)
  slice_shape(1)%stride = 1
  slice_shape(2)%ext = extent(1, jslice)
  slice_shape(2)%stride = 1
  slice_mdt = subarray_mpi_datatype(mpi_double_precision, a_shape, &
       slice_shape, (/ 0, 0, 1 /))
  CALL mpi_type_get_extent(slice_mdt, lb, slice_extent, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_type_get_extent failed", filename, __LINE__)
  ! WRITE(diag_out, '(2(a,i10))') 'extent of slice: ', slice_extent, ', lb: ', lb
  CALL mpi_type_get_true_extent(slice_mdt, lb, slice_extent, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_type_get_extent failed", filename, __LINE__)
  ! WRITE(20, '(2(a,i10))') 'true extent of slice: ', slice_extent, ', lb: ', lb
  ! FLUSH(20)
  CALL mpi_isend(a(4, 5, 1), 3, slice_mdt, dst, 0, mpi_comm_world, req, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_isend failed", filename, __LINE__)
  CALL mpi_recv(a(6, 2, 6), 3, slice_mdt, src, 0, mpi_comm_world, &
       status(:,1), ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_recv failed", filename, __LINE__)
  CALL mpi_wait(req, status(:, 2), ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_wait failed", filename, __LINE__)
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_finalize failed", filename, __LINE__)
END PROGRAM test_strided_extents
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
