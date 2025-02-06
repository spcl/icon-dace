!>
!! @file test_redist_single_array_base_parallel_f.f90
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
PROGRAM test_redist_single_array_base_parallel_f
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_redist_msg, &
       xt_config, xt_config_delete

  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_redist_common, ONLY: &
       test_redist_single_array_base, redist_exchanger_option
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_redist_single_array_base_parallel_f.f90'
  TYPE(xt_config) :: config

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL test_round_robin(mpi_comm_world, config)
  CALL test_allgather(mpi_comm_world, config)
  CALL test_scatter(mpi_comm_world, config)

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)

  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi


CONTAINS

  SUBROUTINE test_round_robin(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    TYPE(xt_redist_msg) :: send_msgs(1), recv_msgs(1)

    INTEGER, PARAMETER :: num_elem = 1
    DOUBLE PRECISION :: src_data(num_elem)
    DOUBLE PRECISION :: ref_dst_data(num_elem)
    INTEGER :: comm_rank, comm_size, ierror

    CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)
    CALL mpi_comm_size(comm, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)

    send_msgs(1)%rank = MOD(comm_rank + 1, comm_size)
    send_msgs(1)%datatype = MPI_DOUBLE_PRECISION
    recv_msgs(1)%rank = MOD(comm_rank + comm_size - 1, comm_size)
    recv_msgs(1)%datatype = MPI_DOUBLE_PRECISION

    src_data(1) = DBLE(comm_rank)
    ref_dst_data(1) = DBLE(MOD(comm_rank + comm_size - 1, comm_size))

    CALL test_redist_single_array_base(send_msgs, recv_msgs, src_data, &
         ref_dst_data, comm, config)

  END SUBROUTINE test_round_robin

  SUBROUTINE test_allgather(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    TYPE(xt_redist_msg), ALLOCATABLE :: send_msgs(:), recv_msgs(:)

    DOUBLE PRECISION :: src_data(1)
    DOUBLE PRECISION, ALLOCATABLE :: ref_dst_data(:)

    INTEGER :: comm_rank, comm_size, i, ierror

    CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)

    CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)

    ALLOCATE(send_msgs(comm_size), recv_msgs(comm_size), &
         ref_dst_data(comm_size))
    DO i = 1, comm_size
      send_msgs(i)%rank = i - 1
      send_msgs(i)%datatype = MPI_DOUBLE_PRECISION
      recv_msgs(i)%rank = i - 1
      CALL MPI_Type_create_indexed_block( &
           1, 1, (/i - 1/), MPI_DOUBLE_PRECISION, recv_msgs(i)%datatype, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_type_create_indexed_block", &
           filename, __LINE__)
      CALL MPI_Type_commit(recv_msgs(i)%datatype, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_type_commit", filename, __LINE__)
    END DO

    src_data(1) = DBLE(comm_rank)
    DO i = 1, comm_size
      ref_dst_data(i) = DBLE(i-1)
    END DO

    CALL test_redist_single_array_base(send_msgs, recv_msgs, src_data, &
         ref_dst_data, comm, config)

    DO i = 1, comm_size
      CALL MPI_Type_free(recv_msgs(i)%datatype, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_type_free", filename, __LINE__)
    END DO

  END SUBROUTINE test_allgather

  SUBROUTINE test_scatter(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    TYPE(xt_redist_msg), ALLOCATABLE :: send_msgs(:)
    TYPE(xt_redist_msg) :: recv_msgs(1)

    DOUBLE PRECISION, ALLOCATABLE :: src_data(:)
    DOUBLE PRECISION :: ref_dst_data(1)

    INTEGER :: comm_size, comm_rank, i, ierror, nsend, rank, displ(1)

    CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)
    ref_dst_data(1) = DBLE(comm_rank)

    CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)

    nsend = MERGE(comm_size, 0, comm_rank == 0)
    ALLOCATE(send_msgs(nsend))
    DO i = 1, nsend
      rank = i - 1
      send_msgs(i)%rank = rank
      displ(1) = rank
      CALL MPI_Type_create_indexed_block( &
           1, 1, displ, MPI_DOUBLE_PRECISION, send_msgs(i)%datatype, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_type_create_indexed_block", &
           filename, __LINE__)
      CALL MPI_Type_commit(send_msgs(i)%datatype, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_type_commit", filename, __LINE__)
    END DO
    recv_msgs(1)%rank = 0
    recv_msgs(1)%datatype = MPI_DOUBLE_PRECISION

    ALLOCATE(src_data(nsend))
    DO i = 1, nsend
      src_data(i) = DBLE(i-1)
    END DO

    CALL test_redist_single_array_base(send_msgs, recv_msgs, src_data, &
         ref_dst_data, comm, config)

    DO i = 1, nsend
      CALL MPI_Type_free(send_msgs(i)%datatype, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_type_free", filename, __LINE__)
    END DO

  END SUBROUTINE test_scatter

END PROGRAM test_redist_single_array_base_parallel_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
