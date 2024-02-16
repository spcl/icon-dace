!>
!! @file row2col_f.f90
!! @brief Fortran example of row to column redistribution
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
PROGRAM test_redist_collection
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, &
       xt_idxlist, xt_idxfsection_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_delete, xt_redist_s_exchange
  USE iso_c_binding, ONLY: c_int
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && ( __PGIC__ == 15 || __PGIC__ == 14 )
  USE xt_redist_real_sp, ONLY: xt_redist_s_exchange
  USE xt_redist_real_dp, ONLY: xt_redist_s_exchange
#endif
  IMPLICIT NONE
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

  INTEGER :: rank, comm_size
  INTEGER(i4), PARAMETER :: global_shape(2) = (/ 1024, 512 /), &
       global_range(2, 2) = RESHAPE((/ 1_i4, global_shape(1), &
       &                               1_i4, global_shape(2) /), (/ 2, 2 /))
  INTEGER :: row_part_range(2, 2), col_part_range(2, 2)
  INTEGER :: i, j
  INTEGER(i4) :: row_part_shape(2), col_part_shape(2)
  INTEGER :: ierror
  TYPE(xt_idxlist) :: src_idxlist, tgt_idxlist
  TYPE(xt_xmap) :: xmap
  TYPE(xt_redist) :: redist
  DOUBLE PRECISION, ALLOCATABLE :: src_array(:, :), tgt_array(:, :)

  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) STOP 1
  CALL xt_initialize(mpi_comm_world)
  CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
  IF (ierror /= mpi_success) STOP 1
  CALL mpi_comm_size (mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) STOP 1

  ! row-wise partitioning
  row_part_range(1, 1) = &
       INT(uniform_partition_start(global_range(:, 1), comm_size, rank + 1))
  row_part_range(2, 1) = &
       INT(uniform_partition_start(global_range(:, 1), comm_size, rank + 2)) - 1
  row_part_range(:, 2) = INT(global_range(:, 2))
  row_part_shape(:) = row_part_range(2, :) - row_part_range(1, :) + 1

  ! column-wise partitioning
  col_part_range(:, 1) = INT(global_range(:, 1))
  col_part_range(1, 2) = &
       INT(uniform_partition_start(global_range(:, 2), comm_size, rank + 1))
  col_part_range(2, 2) = &
       INT(uniform_partition_start(global_range(:, 2), comm_size, rank + 2)) - 1
  col_part_shape(:) = col_part_range(2, :) - col_part_range(1, :) + 1

  ! create decomposition descriptors
  src_idxlist = xt_idxfsection_new(0_xt_int_kind, &
       INT(global_shape, xt_int_kind), INT(row_part_shape, c_int), &
       INT(row_part_range(1, :), xt_int_kind))
  tgt_idxlist = xt_idxfsection_new(0_xt_int_kind, &
       INT(global_shape, xt_int_kind), INT(col_part_shape, c_int), &
       INT(col_part_range(1, :), xt_int_kind))

  ! generate exchange map
  xmap = xt_xmap_all2all_new(src_idxlist, tgt_idxlist, mpi_comm_world)

  ! generate redistribution object
  redist = xt_redist_p2p_new(xmap, mpi_double_precision)

  ! prepare arrays
  ALLOCATE(src_array(row_part_range(1, 1):row_part_range(2, 1), &
       &             row_part_range(1, 2):row_part_range(2, 2)), &
       &   tgt_array(col_part_range(1, 1):col_part_range(2, 1), &
       &             col_part_range(1, 2):col_part_range(2, 2)))

  DO j = row_part_range(1, 2), row_part_range(2, 2)
    DO i = row_part_range(1, 1), row_part_range(2, 1)
      src_array(i, j) = DBLE(i * j)
    END DO
  END DO

  ! do the exchange
  CALL xt_redist_s_exchange(redist, src_array, tgt_array)

  ! clean up
  DEALLOCATE(tgt_array, src_array)
  CALL xt_redist_delete(redist)
  CALL xt_xmap_delete(xmap)
  CALL xt_idxlist_delete(tgt_idxlist)
  CALL xt_idxlist_delete(src_idxlist)

  ! finalise
  CALL xt_finalize()
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) STOP 1
CONTAINS
  FUNCTION uniform_partition_start(set_interval, nparts, part_idx) &
       RESULT(start)
    INTEGER(i4), INTENT(in) :: nparts
    INTEGER(i4), INTENT(in) :: set_interval(2)
    INTEGER(i4), INTENT(in) :: part_idx

    INTEGER(i4) :: start, part_offset

    part_offset = INT((INT(set_interval(2) - set_interval(1) + 1, i8) &
         &             * INT(part_idx - 1, i8)) / INT(nparts, i8))
    start = set_interval(1) + part_offset
  END FUNCTION uniform_partition_start

END PROGRAM test_redist_collection
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
