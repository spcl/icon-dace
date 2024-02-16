!>
!! @file test_redist_collection_static_parallel_f.f90
!! @brief Fortran test of redist_collection_static class
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
PROGRAM test_redist_collection_static_parallel
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xi => xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_stripe, xt_idxvec_new, &
       xt_idxsection_new, xt_idxlist_collection_new, xt_idxstripes_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_static_new, &
       xt_redist_copy, xt_redist_delete, xt_redist_s_exchange, &
       xt_idxlist_get_indices, xt_int_mpidt, &
       xt_request, xt_redist_a_exchange, xt_config, xt_config_delete
  USE test_redist_common, ONLY: check_redist_xi, check_wait_request, &
       redist_exchanger_option
  USE iso_c_binding, ONLY: c_loc, c_ptr
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && (__PGIC__ < 12 || (__PGIC__ ==  12 && __PGIC_MINOR__ <= 10))
  USE xt_redist_base, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
#endif
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_redist_collection_static_parallel_f.f90'
  CHARACTER(len=*), PARAMETER :: err_msg(2) = &
       (/ "xt_redist_s_exchange", "xt_redist_a_exchange" /)
  TYPE(xt_config) :: config
  INTEGER :: rank, comm_size, ierror
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_rank failed', filename, __LINE__)
  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_size failed', filename, __LINE__)

  IF (comm_size > 1) THEN
    CALL test_4redist(mpi_comm_world, config)
    CALL test_rr_exchange(mpi_comm_world, config)
  END IF

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE build_idxlists(indices_a, indices_b, indices_all)
    ! redist test with four different redists
    TYPE(xt_idxlist), INTENT(out) :: indices_a, indices_b, indices_all

    TYPE(xt_idxlist) :: indices_a_(2)
    INTEGER :: i
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi
    INTEGER(xt_int_kind) :: global_size(2), local_start(2, 2)
    INTEGER :: local_size(2)

    TYPE(xt_stripe) :: stripe

    global_size(1) = INT(2 * comm_size, xi)
    global_size(2) = INT(comm_size**2, xi)
    local_size = comm_size
    local_start = RESHAPE((/ 0_xi, INT(rank*comm_size, xi), &
         INT(comm_size, xi), INT(comm_size**2-(rank+1)*comm_size, xi) /), &
         (/ 2, 2 /))

    DO i = 1, 2
      indices_a_(i) = xt_idxsection_new(start, global_size, local_size, &
           local_start(:, i))
    END DO
    indices_a = xt_idxlist_collection_new(indices_a_)

    CALL xt_idxlist_delete(indices_a_(1))
    CALL xt_idxlist_delete(indices_a_(2))

    stripe = xt_stripe(INT(rank * 2 * comm_size**2, xi), 1_xi, 2*comm_size**2)
    indices_b = xt_idxstripes_new(stripe)

    stripe = xt_stripe(0_xi, 1_xi, 2*comm_size**3)
    indices_all = xt_idxstripes_new(stripe)
  END SUBROUTINE build_idxlists

  SUBROUTINE test_4redist(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    INTEGER, PARAMETER :: num_tx = 4
    TYPE(xt_idxlist) :: indices_a, indices_b, indices_all
    INTEGER(xt_int_kind), ALLOCATABLE, TARGET :: src(:), dst(:)
    INTEGER(xt_int_kind), POINTER :: index_vector_a(:), &
         index_vector_b(:), index_vector_all(:)
    TYPE(xt_xmap) :: xmaps(num_tx)
    TYPE(xt_redist) :: redists(num_tx), redist, redist_copy
    INTEGER(mpi_address_kind) :: src_displacements(num_tx), &
         dst_displacements(num_tx)
    INTEGER :: i, ierror, size_a, size_b, size_all
    INTEGER(xt_int_kind), POINTER :: results_1(:), &
         results_2(:), results_3(:), results_4(:)

    size_a = 2 * comm_size**2
    size_b = 2 * comm_size**2
    size_all = 2 * comm_size**3

    ALLOCATE(src(size_a + size_b + size_all), dst(size_b + size_a + 2*size_all))

    index_vector_a => src(1:size_a)
    index_vector_b => src(size_a+1:size_a+size_b)
    index_vector_all => src(size_a+size_b+1:)

    results_1 => dst(1:size_b)
    results_2 => dst(size_b+1:size_b+size_a)
    results_3 => dst(size_b+size_a+1:size_b+size_a+size_all)
    results_4 => dst(size_b+size_a+size_all+1:size_b+size_a+2*size_all)

    CALL build_idxlists(indices_a, indices_b, indices_all)

    CALL xt_idxlist_get_indices(indices_a, index_vector_a)
    CALL xt_idxlist_get_indices(indices_b, index_vector_b)
    CALL xt_idxlist_get_indices(indices_all, index_vector_all)

    xmaps(1) = xt_xmap_all2all_new(indices_a, indices_b, comm)
    xmaps(2) = xt_xmap_all2all_new(indices_b, indices_a, comm)
    xmaps(3) = xt_xmap_all2all_new(indices_a, indices_all, comm)
    xmaps(4) = xt_xmap_all2all_new(indices_b, indices_all, comm)

    CALL xt_idxlist_delete(indices_a)
    CALL xt_idxlist_delete(indices_b)
    CALL xt_idxlist_delete(indices_all)

    DO i = 1, num_tx
      redists(i) = xt_redist_p2p_new(xmaps(i), xt_int_mpidt)
      CALL xt_xmap_delete(xmaps(i))
    END DO

    CALL mpi_get_address(index_vector_a, src_displacements(1), ierror)
    CALL mpi_get_address(index_vector_b, src_displacements(2), ierror)
    CALL mpi_get_address(index_vector_a, src_displacements(3), ierror)
    CALL mpi_get_address(index_vector_b, src_displacements(4), ierror)

    src_displacements = src_displacements - src_displacements(1)

    CALL mpi_get_address(results_1, dst_displacements(1), ierror)
    CALL mpi_get_address(results_2, dst_displacements(2), ierror)
    CALL mpi_get_address(results_3, dst_displacements(3), ierror)
    CALL mpi_get_address(results_4, dst_displacements(4), ierror)

    dst_displacements = dst_displacements - dst_displacements(1)

    redist = xt_redist_collection_static_new(redists, num_tx, &
         src_displacements, dst_displacements, comm, config)

    ! test communicator of redist
    ! if (!test_communicator(xt_redist_get_MPI_Comm(redist), COMM))
    !   PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    CALL xt_redist_delete(redists)

    CALL test_transpose_gather(redist, dst, size_a, size_b, size_all, &
         index_vector_a, index_vector_b, index_vector_all)
    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    CALL test_transpose_gather(redist_copy, dst, size_a, size_b, size_all, &
         index_vector_a, index_vector_b, index_vector_all)

    ! clean up
    CALL xt_redist_delete(redist_copy)
  END SUBROUTINE test_4redist

  SUBROUTINE test_transpose_gather(redist, dst, size_a, size_b, &
       size_all, index_vector_a, index_vector_b, index_vector_all)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: size_a, size_b, size_all
    INTEGER(xi), TARGET, INTENT(inout) :: dst(size_b+size_a+2*size_all)
    INTEGER(xi), TARGET, INTENT(in) :: index_vector_a(size_a)
    INTEGER(xi), INTENT(in) :: index_vector_b(size_b), &
         index_vector_all(size_all)

    INTEGER(xi), POINTER :: results_1(:), &
         results_2(:), results_3(:), results_4(:)
    TYPE(c_ptr) :: results(1), input(1)
    INTEGER :: iexch
    TYPE(xt_request) :: request

    results_1 => dst(1:size_b)
    results_2 => dst(size_b+1:size_b+size_a)
    results_3 => dst(size_b+size_a+1:size_b+size_a+size_all)
    results_4 => dst(size_b+size_a+size_all+1:size_b+size_a+2*size_all)

    input(1) = C_LOC(index_vector_a(1))
    results(1) = C_LOC(results_1(1))

    DO iexch = 1, 2
      dst = 0_xi

      IF (iexch == 1) THEN
        CALL xt_redist_s_exchange(redist, 1, input, results)
      ELSE
        CALL xt_redist_a_exchange(redist, 1, input, results, request)
        CALL check_wait_request(request, filename, __LINE__)
      ENDIF
      ! check results
      IF (ANY(results_1(:) /= index_vector_b)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)

      IF (ANY(results_2(:) /= index_vector_a)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)

      IF (ANY(results_3(:) /= index_vector_all)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)

      IF (ANY(results_4(:) /= index_vector_all)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)
    ENDDO
  END SUBROUTINE test_transpose_gather

  ! redist test with two redists that do a round robin exchange in
  ! different directions
  SUBROUTINE test_rr_exchange(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    INTEGER, PARAMETER :: num_local_indices = 5
    INTEGER(xi) :: src_indices(num_local_indices)
    INTEGER(xi) :: i_xi, temp, dst_indices(num_local_indices, 2)
    INTEGER(xi) :: num_indices_global
    INTEGER :: i
    TYPE(xt_xmap) :: xmaps(2)
    TYPE(xt_redist) :: redists(2), redist
    INTEGER(xi) :: results(num_local_indices, 2)
    INTEGER(mpi_address_kind) :: src_displacements(2), dst_displacements(2), &
         addr_temp
    INTEGER :: ierror

    num_indices_global = INT(comm_size, xi) * INT(num_local_indices, xi)
    DO i = 1, num_local_indices
      i_xi = INT(i, xi)
      src_indices(i) &
           = INT(rank, xi) * INT(num_local_indices, xi) + (i_xi - 1_xi)
      dst_indices(i, 1) = MOD(src_indices(i) + 1_xi, num_indices_global)
      temp = src_indices(i) - 1_xi
      dst_indices(i, 2) = MERGE(num_indices_global - 1_xi, temp, temp < 0_xi)
    END DO

    src_idxlist = xt_idxvec_new(src_indices, num_local_indices)
    DO i = 1, 2
      dst_idxlist = xt_idxvec_new(dst_indices(:, i))
      xmaps(i) = xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm)
      CALL xt_idxlist_delete(dst_idxlist)
      redists(i) = xt_redist_p2p_new(xmaps(i), xt_int_mpidt)
      CALL xt_xmap_delete(xmaps(i))
    END DO

    CALL xt_idxlist_delete(src_idxlist)

    src_displacements = 0_mpi_address_kind
    dst_displacements(1) = 0_mpi_address_kind
    CALL mpi_get_address(results(1, 2), dst_displacements(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error in mpi_get_address", filename, __LINE__)
    CALL mpi_get_address(results(1, 1), addr_temp, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error in mpi_get_address", filename, __LINE__)
    dst_displacements(2) = dst_displacements(2) - addr_temp

    redist = xt_redist_collection_static_new(redists, 2, src_displacements, &
         dst_displacements, comm, config)

    ! test communicator of redist
    ! IF (!test_communicator(xt_redist_get_MPI_Comm(redist), COMM))
    !     PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    CALL xt_redist_delete(redists)

    CALL check_redist_xi(redist, num_local_indices, src_indices, &
         SIZE(results), results, dst_indices)

    ! clean up
    CALL xt_redist_delete(redist)
  END SUBROUTINE test_rr_exchange

END PROGRAM test_redist_collection_static_parallel
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
