!>
!! @file test_set_repartition_mp.f90
!! @brief test whether multi-process repartitioning produced better balance
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
PROGRAM test_set_repartition_mp
  USE ppm_std_type_kinds, ONLY: dp, i4
  USE ppm_std_type_kinds_mp, ONLY: mp_dp, mp_i4
  USE ppm_base, ONLY: abort_ppm
  USE ppm_random, ONLY: irandp, drand_normal, lrand
  USE ppm_set_partition_base, ONLY: balance_of_max_mp
  USE ppm_set_repartition, ONLY: repartition_swap_mp
  USE ppm_extents, ONLY: extent, extent_start, extent_end, ASSIGNMENT(=), &
       extent_size, iinterval
  USE ppm_combinatorics, ONLY: is_permutation
  USE scales_ppm, ONLY: initialize_scales_ppm, finalize_scales_ppm
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif

  INTEGER, ALLOCATABLE :: part_sizes(:), offsets(:)

  INTEGER :: rand_seed
  INTEGER :: comm, comm_size, comm_rank, ierror
  TYPE(extent) :: my_range
  REAL :: balance
  CHARACTER(len=*), PARAMETER :: filename = 'test_set_repartition_mp.f90'

  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_init failed", filename, __LINE__)

  comm = mpi_comm_world
  CALL initialize_scales_ppm(comm, 0, rand_seed)

  CALL mpi_comm_rank(comm, comm_rank, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_comm_rank failed", filename, __LINE__)

  PRINT '(2(a,i0))', 'rank: ', comm_rank, ', random seed=', rand_seed

  CALL mpi_comm_size(comm, comm_size, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_comm_size failed", filename, __LINE__)

  ALLOCATE(offsets(0:comm_size), part_sizes(comm_size))
  !===========================================================================
  ! test for integer weights
  CALL test_repartition_i4

  !===========================================================================
  ! test for double precision weights
  CALL test_repartition_dp

  CALL finalize_scales_ppm
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_finalize failed", filename, __LINE__)

  DEALLOCATE(offsets, part_sizes)

CONTAINS
  SUBROUTINE test_repartition_i4
    INTEGER(i4), ALLOCATABLE :: weight(:), all_weights(:)
    INTEGER, ALLOCATABLE :: part(:)
    INTEGER :: num_elems, i, ierror, rs, re, rsize

    num_elems = my_part_size()

    ALLOCATE(weight(num_elems))

    DO i = 1, num_elems-2, 3
      weight(i) = MAX(INT(drand_normal(1000.0_dp, 200.0_dp)), 800)
      weight(i+1) = MAX(INT(drand_normal(1000.0_dp, 200.0_dp)), 800)
      weight(i+2) = MAX(INT(drand_normal(2350.0_dp, 500.0_dp)), 800)
    END DO
    weight(num_elems - 1) = MAX(INT(drand_normal(1000.0_dp, 200.0_dp)), 800)
    weight(num_elems) = MAX(INT(drand_normal(2350.0_dp, 500.0_dp)), 800)

    CALL mpi_allgather(num_elems, 1, mpi_integer, part_sizes, 1, mpi_integer, &
         comm, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_allgather failed", filename, __LINE__)
    offsets(0) = 0
    DO i = 1, comm_size
      offsets(i) = offsets(i - 1) + part_sizes(i)
    END DO
    my_range = extent(offsets(comm_rank)+1, num_elems)
    rsize = extent_size(my_range)
    ALLOCATE(part(rsize))
    rs = extent_start(my_range)
    re = extent_end(my_range)
    DO i = 1, rsize
      part(i) = i-1+rs
    END DO
    balance = balance_of_max_mp(weight)
    IF (comm_rank == 0) &
         PRINT '(a, f8.5)', 'imbalance of uniform partition: ', balance

    CALL repartition_swap_mp(part, weight)

    ALLOCATE(all_weights(offsets(comm_size)))
    CALL mpi_allgatherv(weight, num_elems, mp_i4, &
         all_weights, part_sizes, offsets, mp_i4, &
         comm, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_allgatherv failed", filename, __LINE__)

    DO i = 1, rsize
      weight(i) = all_weights(part(i))
    END DO
    DEALLOCATE(all_weights)
    balance = balance_of_max_mp(weight)
    DEALLOCATE(weight)

    IF (comm_rank == 0) &
         PRINT '(a,f8.5)', 'imbalance after repartitioning: ', balance
    CALL mpi_barrier(comm, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_barrier failed", filename, __LINE__)
    CALL check_partition_property(part, iinterval(1, offsets(comm_size)))


    DEALLOCATE(part)

  END SUBROUTINE test_repartition_i4

  SUBROUTINE test_repartition_dp
    REAL(dp), ALLOCATABLE :: weight(:), all_weights(:)
    INTEGER, ALLOCATABLE :: part(:)
    INTEGER :: num_elems, i, ierror, rs, re, rsize

    num_elems = my_part_size()

    ALLOCATE(weight(num_elems))

    DO i = 1, num_elems-2, 3
      weight(i) = MAX(drand_normal(1000.0_dp, 200.0_dp), 800.0_dp)
      weight(i+1) = MAX(drand_normal(1000.0_dp, 200.0_dp), 800.0_dp)
      weight(i+2) = MAX(drand_normal(2350.0_dp, 500.0_dp), 800.0_dp)
    END DO
    weight(num_elems - 1) = MAX(drand_normal(1000.0_dp, 200.0_dp), 800.0_dp)
    weight(num_elems) = MAX(drand_normal(2350.0_dp, 500.0_dp), 800.0_dp)

    CALL mpi_allgather(num_elems, 1, mpi_integer, part_sizes, 1, mpi_integer, &
         comm, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_allgather failed", filename, __LINE__)
    offsets(0) = 0
    DO i = 1, comm_size
      offsets(i) = offsets(i - 1) + part_sizes(i)
    END DO
    my_range = extent(offsets(comm_rank)+1, num_elems)
    rsize = extent_size(my_range)
    ALLOCATE(part(rsize))
    rs = extent_start(my_range)
    re = extent_end(my_range)
    DO i = 1, rsize
      part(i) = i-1+rs
    END DO

    balance = balance_of_max_mp(weight)
    IF (comm_rank == 0) &
         PRINT '(a, f8.5)', 'imbalance of uniform partition: ', balance

    CALL repartition_swap_mp(part, weight)

    ALLOCATE(all_weights(offsets(comm_size)))
    CALL mpi_allgatherv(weight, num_elems, mp_dp, &
         all_weights, part_sizes, offsets, &
         mp_dp, comm, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_allgatherv failed", filename, __LINE__)

    DO i = 1, rsize
      weight(i) = all_weights(part(i))
    END DO
    DEALLOCATE(all_weights)
    balance = balance_of_max_mp(weight)
    DEALLOCATE(weight, part)

    IF (comm_rank == 0) &
         PRINT '(a,f8.5)', 'imbalance after repartitioning: ', balance
  END SUBROUTINE test_repartition_dp

  FUNCTION my_part_size()
    INTEGER :: my_part_size

    INTEGER, PARAMETER :: max_max_num_elems = 30
    INTEGER, PARAMETER :: min_min_num_elems = 5, max_min_num_elems = 20
    INTEGER :: max_num_elems, min_num_elems, buf(2)

    IF (comm_rank == 0) THEN
      ! decide whether to use uniform or randomly sized partitions
      min_num_elems = MOD(irandp(), max_min_num_elems - min_min_num_elems + 1) &
           + min_min_num_elems
      max_num_elems = MOD(irandp(), max_max_num_elems - min_num_elems + 1) &
           + min_num_elems
      IF (lrand()) THEN
        PRINT '(a,2(i0,a))', 'using random size parts between ', &
             min_num_elems, ' and ', max_num_elems, ' elements in size'
      ELSE
        min_num_elems = max_num_elems
        PRINT '(a,i0,a)', 'using uniform parts ', &
             max_num_elems, ' elements in size'
      END IF
    END IF
    buf(1) = min_num_elems; buf(2) = max_num_elems
    CALL mpi_bcast(buf, 2, mpi_integer, 0, comm, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_bcast failed", filename, __LINE__)
    min_num_elems = buf(1); max_num_elems = buf(2)
    my_part_size &
         = MOD(irandp(), max_num_elems - min_num_elems + 1) + min_num_elems
  END FUNCTION my_part_size

  SUBROUTINE check_partition_property(part, global_range)
    INTEGER, ALLOCATABLE, INTENT(in) :: part(:)
    TYPE(iinterval), INTENT(in) :: global_range

    INTEGER(i4), ALLOCATABLE :: all_parts(:)
    INTEGER :: i

    IF (comm_rank == 0) THEN
      ALLOCATE(all_parts(offsets(comm_size)))
    ELSE
      ALLOCATE(all_parts(1))
    END IF
    CALL mpi_gatherv(part, SIZE(part), mpi_integer, &
         all_parts, part_sizes, offsets, mpi_integer, &
         0, comm, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_gatherv failed", filename, __LINE__)
    IF (comm_rank == 0) THEN
      IF (.NOT. is_permutation(all_parts, global_range)) THEN
        DO i = 0, comm_size-1
          PRINT '(a,i0,a,100(i0,", "))', 'rank: ', i, &
               ', part: ', all_parts(offsets(i)+1:offsets(i+1))
        END DO
        CALL abort_ppm(source=filename, line=__LINE__, &
             msg="ppm_repartition_swap_mp failed: result is not &
             &partition")
      END IF
    END IF
  END SUBROUTINE check_partition_property

END PROGRAM test_set_repartition_mp
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
