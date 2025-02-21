!>
!! @file ppm_extents_mp.f90
!! @brief helper module to create mpi datatypes for extent type
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

MODULE ppm_extents_mp
  USE ppm_base, ONLY: abort_ppm, assertion
  USE ppm_extents, ONLY: extent, OPERATOR(==)
  USE ppm_std_type_kinds, ONLY: i4
  USE ppm_std_type_kinds_mp, ONLY: mp_i4



  IMPLICIT NONE
  PRIVATE




  INTEGER, SAVE :: extent_mp

  PUBLIC :: extent_mp, create_extents_mp, destroy_extents_mp

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_extents_mp.f90'
CONTAINS
  SUBROUTINE create_extents_mp
    INTEGER :: ierror, comm_self_clone, i, request
    INTEGER, PARAMETER :: msg_count = 5
    TYPE(extent) :: a(msg_count), b(msg_count)
    CALL mpi_type_create_struct(1, (/ 2 /), (/ 0_mpi_address_kind /), &
         (/ mp_i4 /), &
         extent_mp, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_create_struct failed", filename, 69)
    CALL mpi_type_commit(extent_mp, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_commit failed", filename, 72)
    CALL mpi_comm_dup(mpi_comm_self, comm_self_clone, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_dup failed", filename, 75)
    a(1) = extent(123456_i4, 78901_i4)
    DO i = 2, SIZE(a)
      a(i)%first = a(i - 1)%first + 333
      a(i)%size  = a(i - 1)%size + 555
    END DO
    CALL mpi_isend(a, msg_count, extent_mp, 0, 1, comm_self_clone, &
         request, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_isend failed", filename, 84)
    CALL mpi_recv(b, msg_count, extent_mp, 0, 1, comm_self_clone, &
         mpi_status_ignore, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_recv failed", filename, 88)
    CALL mpi_wait(request, mpi_status_ignore, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_wait failed", filename, 91)
    CALL assertion(request == mpi_request_null .AND. ALL(a == b), &
         filename, 93, 'error in transfer')
    CALL mpi_comm_free(comm_self_clone, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_free failed", filename, 96)
    CALL ppm_create_extents_mp
  END SUBROUTINE create_extents_mp

  SUBROUTINE destroy_extents_mp
    INTEGER :: ierror
    CALL mpi_type_free(extent_mp, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_type_free failed", filename, 104)
  END SUBROUTINE destroy_extents_mp
END MODULE ppm_extents_mp
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
