!>
!! @file test_set_partition.f90
!! @brief test balance of partition produced by partition routine
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
PROGRAM test_set_partition
  USE ppm_base, ONLY: abort_ppm, assertion, ppm_default_comm
  USE ppm_extents, ONLY: iinterval
  USE ppm_random, ONLY: a_randp, initialize_irand, irandr
  USE ppm_set_partition, ONLY: greedy_partitioning
  USE ppm_set_partition_base, ONLY: set_i4, partition_vec, OPERATOR(==)
  USE ppm_std_type_kinds, ONLY: i4, i8
  IMPLICIT NONE
  INTEGER, PARAMETER :: maxasize = 30000, minasize = 2, max_nparts=5000
  TYPE(iinterval), PARAMETER :: &
       set_size_range = iinterval(minasize, maxasize), &
       num_parts_range = iinterval(1, max_nparts)
  INTEGER(i4), ALLOCATABLE :: a(:)
  INTEGER(i8), ALLOCATABLE :: weight_sum(:)
  INTEGER :: asize, i, j, m, rand_seed, nparts
  TYPE(set_i4), ALLOCATABLE :: partition(:)
  TYPE(partition_vec) :: partitioning
  LOGICAL, ALLOCATABLE :: seen(:)
  INTEGER :: psize_sum
  CHARACTER(len=*), PARAMETER :: filename = 'test_set_partition.f90'

  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  asize = irandr(set_size_range)
  ALLOCATE(a(asize))
  CALL a_randp(a)

  nparts = irandr(num_parts_range)
  ALLOCATE(partition(nparts))

  CALL greedy_partitioning(partition, a)
  ALLOCATE(partitioning%start(nparts+1), partitioning%elements(asize))
  CALL greedy_partitioning(partitioning, a)

  CALL assertion(partition == partitioning, filename, __LINE__, &
       "same algorithm gave different results")
  psize_sum = 0
  DO i = 1, nparts
    psize_sum = psize_sum + SIZE(partition(i)%elem)
  END DO
  CALL assertion(psize_sum == asize, line=__LINE__, source=filename, &
       msg='not a partitioning: number of elements in total incorrect')

  ALLOCATE(seen(asize))
  ALLOCATE(weight_sum(nparts))
  seen = .FALSE.
  DO i = 1, nparts
    weight_sum(i) = 0_i8
    DO j = 1, SIZE(partition(i)%elem)
      m  = partition(i)%elem(j)
      IF (seen(m)) CALL abort_ppm(msg='duplicate element in partitioning', &
           line=__LINE__, source=filename)
      seen(m) = .TRUE.
      weight_sum(i) = weight_sum(i) + INT(a(m), i8)
    END DO
  END DO
  CALL assertion(ALL(seen), line=__LINE__, source=filename, &
       msg='not a partition: not surjective')
  PRINT '(3(a,i0))', 'maximum weight=', MAXVAL(weight_sum), &
       ', mean weight=', SUM(weight_sum)/INT(nparts, i8), &
       ', min weight=', MINVAL(weight_sum)
END PROGRAM test_set_partition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
