!>
!! @file test_set_repartition.f90
!! @brief test whether repartitioning actually produced more balanced partition
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
PROGRAM test_set_repartition
  USE ppm_std_type_kinds, ONLY: dp
  USE ppm_base, ONLY: ppm_default_comm
  USE ppm_random, ONLY: irandp, initialize_irand, drand_normal
  USE ppm_set_partition_base, ONLY: partition_vec, balance_of_max
  USE ppm_set_repartition, ONLY: repartition_swap
  USE ppm_extents, ONLY: extent, extent_start, extent_end
  USE ppm_uniform_partition, ONLY: uniform_decomposition
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_num_elems = 30000, min_num_elems = 2, &
       max_num_parts = 8000
  INTEGER(i4), ALLOCATABLE :: weight_i(:)
  REAL(dp), ALLOCATABLE :: weight_dp(:)
  INTEGER :: rand_seed, num_elems, num_parts, i, j, lb, ub
  TYPE(partition_vec) :: partitioning
  TYPE(extent), ALLOCATABLE :: uniform_partitions(:)

  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  num_elems = MOD(irandp(), max_num_elems - min_num_elems + 1) + min_num_elems
  num_parts = MOD(irandp(), max_num_parts + 1) + 1

  ALLOCATE(partitioning%start(num_parts+1), partitioning%elements(num_elems), &
       uniform_partitions(num_parts), weight_i(num_elems))

  CALL uniform_decomposition(extent(1, num_elems), num_parts, &
       uniform_partitions)
  DO i = 1, num_parts
    partitioning%start(i) = extent_start(uniform_partitions(i))
    lb = extent_start(uniform_partitions(i))
    ub = extent_end(uniform_partitions(i))
    DO j = lb, ub
      partitioning%elements(j) = j
    END DO
  END DO
  partitioning%start(num_parts + 1) &
       = extent_end(uniform_partitions(num_parts)) + 1

  DO i = 1, num_elems-2, 3
    weight_i(i) = MAX(INT(drand_normal(1000.0_dp, 200.0_dp)), 800)
    weight_i(i+1) = MAX(INT(drand_normal(1000.0_dp, 200.0_dp)), 800)
    weight_i(i+2) = MAX(INT(drand_normal(2350.0_dp, 500.0_dp)), 800)
  END DO
  weight_i(num_elems - 1) = MAX(INT(drand_normal(1000.0_dp, 200.0_dp)), 800)
  weight_i(num_elems) = MAX(INT(drand_normal(2350.0_dp, 500.0_dp)), 800)

  PRINT *, 'number of partitions: ', num_parts, 'number of elements: ', &
       num_elems
  PRINT *, 'imbalance of uniform partition:', &
       balance_of_max(partitioning, weight_i)

  CALL repartition_swap(partitioning, weight_i)

  PRINT *, 'imbalance after repartitioning:', &
       balance_of_max(partitioning, weight_i)

  DEALLOCATE(partitioning%start, partitioning%elements, uniform_partitions, &
       weight_i)

  num_elems = MOD(irandp(), max_num_elems - min_num_elems + 1) + min_num_elems
  num_parts = MOD(irandp(), max_num_parts + 1) + 1

  ALLOCATE(partitioning%start(num_parts+1), partitioning%elements(num_elems), &
       uniform_partitions(num_parts), weight_dp(num_elems))

  CALL uniform_decomposition(extent(1, num_elems), num_parts, &
       uniform_partitions)
  DO i = 1, num_parts
    partitioning%start(i) = extent_start(uniform_partitions(i))
    lb = extent_start(uniform_partitions(i))
    ub = extent_end(uniform_partitions(i))
    DO j = lb, ub
      partitioning%elements(j) = j
    END DO
  END DO
  partitioning%start(num_parts + 1) &
       = extent_end(uniform_partitions(num_parts)) + 1

  DO i = 1, num_elems-2, 3
    weight_dp(i) = MAX(drand_normal(1000.0_dp, 200.0_dp), 800.0_dp)
    weight_dp(i+1) = MAX(drand_normal(1000.0_dp, 200.0_dp), 800.0_dp)
    weight_dp(i+2) = MAX(drand_normal(2350.0_dp, 500.0_dp), 800.0_dp)
  END DO
  weight_dp(num_elems - 1) = MAX(drand_normal(1000.0_dp, 200.0_dp), 800.0_dp)
  weight_dp(num_elems) = MAX(drand_normal(2350.0_dp, 500.0_dp), 800.0_dp)

  PRINT *, 'number of partitions: ', num_parts, 'number of elements: ', &
       num_elems
  PRINT *, 'imbalance of uniform partition:', &
       balance_of_max(partitioning, weight_dp)

  CALL repartition_swap(partitioning, weight_dp)

  PRINT *, 'imbalance after repartitioning:', &
       balance_of_max(partitioning, weight_dp)

END PROGRAM test_set_repartition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
