!>
!! @file test_redist_single_array_base_f.f90
!!
!! @copyright Copyright  (C)  2017 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!
!  Keywords:
!  Maintainer: Jörg Behrens <behrens@dkrz.de>
!              Moritz Hanke <hanke@dkrz.de>
!              Thomas Jahns <jahns@dkrz.de>
!  URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are  permitted provided that the following conditions are
!  met:
!
!  Redistributions of source code must retain the above copyright notice,
!  this list of conditions and the following disclaimer.
!
!  Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!
!  Neither the name of the DKRZ GmbH nor the names of its contributors
!  may be used to endorse or promote products derived from this software
!  without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
#include "fc_feature_defs.inc"
PROGRAM test_redist_single_array_base_f
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_redist_msg, &
       xt_config, xt_config_delete

  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_redist_common, ONLY: &
       test_redist_single_array_base, redist_exchanger_option
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_redist_single_array_base_f.f90'
  TYPE(xt_config) :: config

  ! init mpi
  CALL init_mpi

  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  ! single double
  CALL test_single_double(mpi_comm_world, config)
  ! reverse order of some doubles
  CALL test_reverse_doubles(mpi_comm_world, config)

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)

  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE test_single_double(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config


    TYPE(xt_redist_msg) :: send_msgs(1)
    TYPE(xt_redist_msg) :: recv_msgs(1)

    INTEGER, PARAMETER :: num_elem = 1
    DOUBLE PRECISION, PARAMETER :: src_data(num_elem) &
         = (/ 0.0d0 /)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(num_elem) &
         = (/ 0.0d0 /)

    send_msgs(1)%rank = 0
    send_msgs(1)%datatype = MPI_DOUBLE_PRECISION
    recv_msgs(1)%rank = 0
    recv_msgs(1)%datatype = MPI_DOUBLE_PRECISION

    CALL test_redist_single_array_base(send_msgs, recv_msgs, src_data, &
         ref_dst_data, comm, config)

  END SUBROUTINE test_single_double

  SUBROUTINE test_reverse_doubles(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config


    TYPE(xt_redist_msg) :: send_msgs(1)
    TYPE(xt_redist_msg) :: recv_msgs(1)

    INTEGER :: i, ierror
    INTEGER, PARAMETER :: num_elem = 10
    INTEGER, PARAMETER :: displ(num_elem) &
         = (/ (i, i = num_elem - 1, 0, -1) /)
#ifndef __PGI
    DOUBLE PRECISION, PARAMETER :: src_data(num_elem) &
         = (/ (DBLE(i), i = 1, num_elem) /)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(num_elem) &
         = (/ (DBLE(i), i = num_elem, 1, -1) /)
#else
    DOUBLE PRECISION :: src_data(num_elem), ref_dst_data(num_elem)
#endif

#ifdef __PGI
    DO i = 1, num_elem
      src_data(i) = DBLE(i)
      ref_dst_data(i) = DBLE(num_elem - i + 1)
    END DO
#endif
    send_msgs(1)%rank = 0
    CALL MPI_Type_contiguous( &
      num_elem, MPI_DOUBLE_PRECISION, send_msgs(1)%datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_type_contiguous", &
         filename, __LINE__)
    CALL MPI_Type_commit(send_msgs(1)%datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_type_commit", &
         filename, __LINE__)
    recv_msgs(1)%rank = 0
    CALL MPI_Type_create_indexed_block(num_elem, 1, displ, &
         MPI_DOUBLE_PRECISION, recv_msgs(1)%datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_type_create_indexed_block", &
         filename, __LINE__)
    CALL MPI_Type_commit(recv_msgs(1)%datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_type_commit", &
         filename, __LINE__)

    CALL test_redist_single_array_base(send_msgs, recv_msgs, src_data, &
         ref_dst_data, comm, config)

    CALL MPI_Type_free(recv_msgs(1)%datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_type_free", filename, __LINE__)
    CALL MPI_Type_free(send_msgs(1)%datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_type_free", filename, __LINE__)

  END SUBROUTINE test_reverse_doubles

END PROGRAM test_redist_single_array_base_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
