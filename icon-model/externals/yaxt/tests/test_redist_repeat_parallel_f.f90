!>
!! @file test_redist_repeat_parallel_f.f90
!! @brief Parallelized Fortran test of redist_repeat class
!!
!! @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
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
PROGRAM test_redist_repeat_parallel
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xi => xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_stripe, &
       xt_idxfsection_new, xt_idxlist_collection_new, xt_idxstripes_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_repeat_new, &
       xt_redist_delete, xt_redist_s_exchange, &
       xt_idxlist_get_indices, xt_int_mpidt, &
       xt_request, xt_redist_a_exchange, xt_config, xt_config_delete
    ! older PGI compilers do not handle generic interface correctly
#if defined __PGI
  USE xt_redist_int_i2, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_int_i4, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_int_i8, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
#endif
  USE test_redist_common, ONLY: check_wait_request, redist_exchanger_option
  USE iso_c_binding, ONLY: c_int
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: filename = 'test_redist_repeat_parallel_f.f90'
  CHARACTER(len=*), PARAMETER :: err_msg(2) = &
       (/ "error on xt_redist_s_exchange", "error on xt_redist_a_exchange" /)
  TYPE(xt_config) :: config
  INTEGER :: comm_size, ierror
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_size failed', filename, __LINE__)

  IF (comm_size > 1) THEN
    CALL test_4redist(mpi_comm_world, config, 2*comm_size**2)
  END IF

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  ! create index lists to exchange data sections from a global 4D
  ! array.
  !
  ! For the source side, the global array is of size W x X x Y x Z,
  ! where W=X=Y=comm_size and Z=2. The source data is decomposed
  ! into two shards of size comm_size*comm_size per process, where one
  ! shard is positioned at (1, comm_rank, 1, 1), the other at (1,
  ! comm_size-comm_rank, 1, 2), i.e. the decomposition of the last
  ! dimension is decomposed anti-symmetrically.
  !
  ! The destination decomposition is a contiguous subset of the
  ! interval [0,W*X*Y*Z-1], the stripe [S, S+2*comm_size^2] with S =
  ! comm_rank * 2 * comm_size**2. This corresponds to a section of a
  ! 4D array reshaped to [W,X,Z,Y], decomposed along the Y-axis
  ! according to comm_rank (but enumerated differently from source
  ! array).
  SUBROUTINE build_idxlists(indices_a, indices_b, comm_size, comm_rank)
    TYPE(xt_idxlist), INTENT(out) :: indices_a, indices_b
    INTEGER, INTENT(in) :: comm_size, comm_rank

    INTEGER, PARAMETER :: glob_rank = 4
    TYPE(xt_idxlist) :: indices_a_(2)
    INTEGER :: i
    INTEGER(xt_int_kind), PARAMETER :: start = 0
    INTEGER(xt_int_kind) :: global_size(glob_rank), local_start(glob_rank, 2)
    INTEGER :: local_size(glob_rank)

    TYPE(xt_stripe) :: stripe

    global_size(1) = INT(comm_size, xi)
    global_size(2) = INT(comm_size, xi)
    global_size(3) = INT(comm_size, xi)
    global_size(4) = 2_xi
    local_size(1) = comm_size
    local_size(2) = 1
    local_size(3) = comm_size
    local_size(4) = 1
    local_start(1, 1) = 1_xi
    local_start(2, 1) = INT(comm_rank + 1, xi)
    local_start(3, 1) = 1_xi
    local_start(4, 1) = 1_xi
    !
    local_start(1, 2) = 1_xi
    local_start(2, 2) = INT(comm_size-comm_rank, xi)
    local_start(3, 2) = 1_xi
    local_start(4, 2) = 2_xi

    DO i = 1, 2
      indices_a_(i) = xt_idxfsection_new(start, global_size, local_size, &
           local_start(:, i))
    END DO
    indices_a = xt_idxlist_collection_new(indices_a_)

    CALL xt_idxlist_delete(indices_a_(1))
    CALL xt_idxlist_delete(indices_a_(2))

    stripe = xt_stripe(start = INT(comm_rank * 2 * comm_size**2, xi), &
         &             stride = 1_xi, &
         &             nstrides = INT(2*comm_size**2, c_int))
    indices_b = xt_idxstripes_new(stripe)
  END SUBROUTINE build_idxlists

  ! redist test for 4 level repetition of redist (i.e. 3D extension of 2D
  ! redist)
  SUBROUTINE test_4redist(comm, config, dim1)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    INTEGER, INTENT(in) :: dim1
    TYPE(xt_idxlist) :: indices_a, indices_b
    INTEGER(xt_int_kind) :: index_vector_a(dim1), &
                            index_vector_b(dim1)
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist_repeat, redist_repeat_2, redist_p2p
    INTEGER, PARAMETER :: dim2a = 9, rpt_cnt = 4
    INTEGER(xt_int_kind) :: results_1(dim1,rpt_cnt), &
                            results_2(dim1,dim2a), dim1_xi
    INTEGER(xt_int_kind) :: input_data(dim1,dim2a)
    INTEGER(xt_int_kind) :: ref_results_1(dim1,rpt_cnt), &
                            ref_results_2(dim1,dim2a)
    INTEGER(mpi_address_kind) :: extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int), PARAMETER :: &
         displacements(rpt_cnt, 2) &
         = RESHAPE((/ 0_c_int, 1_c_int, 2_c_int, 3_c_int, &
         &            1_c_int, 2_c_int, 4_c_int, 8_c_int /), (/ rpt_cnt, 2 /))
    ! skip_lev_2 must correspond to the levels skipped via displacements_2
    LOGICAL, PARAMETER :: skip_lev_2(9) &
         = (/  .TRUE., .FALSE., .FALSE., &
         &     .TRUE., .FALSE.,  .TRUE., &
         &     .TRUE.,  .TRUE., .FALSE. /)
    INTEGER :: i, j, ierror
    TYPE(xt_request) :: request1, request2
    INTEGER :: iexch
    INTEGER :: comm_rank, comm_size

    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort('mpi_comm_rank failed', filename, __LINE__)
    CALL mpi_comm_size(comm, comm_size, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort('mpi_comm_size failed', filename, __LINE__)

    CALL build_idxlists(indices_a, indices_b, comm_size, comm_rank)

    CALL xt_idxlist_get_indices(indices_a, index_vector_a)
    CALL xt_idxlist_get_indices(indices_b, index_vector_b)

    xmap = xt_xmap_all2all_new(indices_a, indices_b, comm)

    CALL xt_idxlist_delete(indices_a)
    CALL xt_idxlist_delete(indices_b)

    redist_p2p = xt_redist_p2p_new(xmap, xt_int_mpidt)
    CALL xt_xmap_delete(xmap)

    CALL mpi_get_address(input_data(1,1), base_address, ierror)
    CALL mpi_get_address(input_data(1,2), temp_address, ierror)
    extent = temp_address - base_address

    redist_repeat = xt_redist_repeat_new(redist_p2p, extent, extent, &
         rpt_cnt, displacements(:, 1), config)
    redist_repeat_2 = xt_redist_repeat_new(redist_p2p, extent, extent, &
         rpt_cnt, displacements(:, 2), config)

    CALL xt_redist_delete(redist_p2p)

    dim1_xi = INT(dim1, xi)
    DO j = 1, dim2a
      DO i = 1, dim1
        input_data(i, j) = index_vector_a(i) + INT(j-1, xi) * 2_xi * dim1_xi
      END DO
    END DO

    DO j = 1, rpt_cnt
      DO i = 1, dim1
        ref_results_1(i, j) = index_vector_b(i) + INT(j-1, xi) * 2_xi * dim1_xi
      END DO
    END DO
    DO j = 1, dim2a
      IF (skip_lev_2(j)) THEN
        ref_results_2(:, j) = -1_xi
      ELSE
        DO i = 1, dim1
          ref_results_2(i, j) &
               = index_vector_b(i) + INT(j-1, xi) * 2_xi * dim1_xi
        END DO
      END IF
    END DO

    DO iexch = 1, 2
      results_1 = -1
      results_2 = -1

      IF (iexch == 1) THEN
        CALL xt_redist_s_exchange(redist_repeat, input_data, results_1)
        CALL xt_redist_s_exchange(redist_repeat_2, input_data, results_2)
      ELSE
        CALL xt_redist_a_exchange(redist_repeat, input_data, results_1, &
             request1)
        CALL xt_redist_a_exchange(redist_repeat_2, input_data, results_2, &
             request2)
        CALL check_wait_request(request1, filename, __LINE__)
        CALL check_wait_request(request2, filename, __LINE__)
      ENDIF

      ! check results
      IF (ANY(results_1 /= ref_results_1)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)
      IF (ANY(results_2 /= ref_results_2)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)
    ENDDO
    ! clean up

    CALL xt_redist_delete(redist_repeat)
    CALL xt_redist_delete(redist_repeat_2)
  END SUBROUTINE test_4redist

END PROGRAM test_redist_repeat_parallel
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
