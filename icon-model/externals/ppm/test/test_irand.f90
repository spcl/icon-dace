!>
!! @file test_irand.f90
!! @brief random number generator unit-tests
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Maintainer: Thomas Jahns <jahns@dkrz.de>
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
PROGRAM test_irand
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_num_threads, omp_get_thread_num
#endif
  USE ppm_base, ONLY: ppm_default_comm, abort_ppm
  USE ppm_std_type_kinds, ONLY: i4, sp, dp
  USE ppm_extents, ONLY: iinterval, extent, ASSIGNMENT(=), &
       iinterval_sp, iinterval_dp
  USE ppm_random, ONLY: irand, initialize_irand, a_rand, a_randr, irand_max, &
       irandp, frand, frandp, drand, drandp
  IMPLICIT NONE
  INTEGER(i4), ALLOCATABLE :: irand_res(:, :)
  REAL(sp), ALLOCATABLE :: frand_res(:, :)
  REAL(dp), ALLOCATABLE :: drand_res(:, :)
  INTEGER :: num_rands, num_threads, tid, i, ref, rand_seed
  LOGICAL :: p, p_all
  INTEGER :: range_start, range_size
  TYPE(iinterval) :: rand_range
  TYPE(iinterval_sp) :: frand_range
  TYPE(iinterval_dp) :: drand_range
  REAL(sp) :: frange_start, frange_size
  REAL(dp) :: drange_start, drange_size
  CHARACTER(len=*), PARAMETER :: filename = 'test_irand.f90'
!$omp parallel private(tid, num_threads, i, ref, p, rand_seed, &
!$omp rand_range, range_size, range_start, frand_range, frange_start, frange_size,&
!$omp drand_range, drange_start, drange_size) &
!$omp shared(p_all, num_rands, irand_res, frand_res, drand_res)
#ifdef _OPENMP
  num_threads = omp_get_num_threads()
  tid = omp_get_thread_num()
#else
  num_threads = 1
  tid = 0
#endif
  ! test for reproducability in multi-thread setup
  CALL initialize_irand(ppm_default_comm, 9)
  p = .FALSE.
!$omp master
  p_all = .FALSE.
  num_rands = MOD(IAND(irand(), 2147483647), 10000)
  ALLOCATE(irand_res(num_rands, 0:num_threads-1), &
       frand_res(num_rands, 0:num_threads-1), &
       drand_res(num_rands, 0:num_threads-1))
!$omp end master
!$omp barrier
  CALL initialize_irand(ppm_default_comm, 9)
  DO i = 1, num_rands
    irand_res(i, tid) = irand()
  END DO
!$omp barrier
  ref = MOD(tid + 1, num_threads)
  p = p &
       .OR. ANY(irand_res(:, tid) /= irand_res(:, ref))
!$omp barrier
  CALL a_rand(irand_res(:, tid))
!$omp barrier
  p = p &
       .OR. ANY(irand_res(:, tid) /= irand_res(:, ref))
  CALL reduce_p(p, p_all, "Inequality found in random numbers!", &
       filename, __LINE__)
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
!$omp critical
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed
!$omp end critical
  range_start = irand()
  range_size = MOD(irandp(), MERGE(irand_max, irand_max - range_start, &
       range_start <= 0)) + 1
  rand_range = extent(range_start, range_size)
  CALL a_randr(irand_res(:, tid), rand_range)
  p = ANY(irand_res(:, tid) > rand_range%last &
       .OR. irand_res(:, tid) < rand_range%first)
  CALL reduce_p(p, p_all, "Range violation in generated random integers!", &
       filename, __LINE__)

  frange_start = frand() * 2500.0_sp
  frange_size = frandp() * 2500.0_sp
  frand_range = iinterval_sp(frange_start, frange_start + frange_size)
  CALL a_randr(frand_res(:, tid), frand_range)
  p = ANY(frand_res(:, tid) > frand_range%last &
       .OR. frand_res(:, tid) < frand_range%first)
  CALL reduce_p(p, p_all, &
       "Range violation in generated random single precision reals!", &
       filename, __LINE__)

  drange_start = drand() * 2500.0_dp
  drange_size = drandp() * 2500.0_dp
  drand_range = iinterval_dp(drange_start, drange_start + drange_size)
  CALL a_randr(drand_res(:, tid), drand_range)
  p = ANY(drand_res(:, tid) > drand_range%last &
       .OR. drand_res(:, tid) < drand_range%first)
  CALL reduce_p(p, p_all, &
       "Range violation in generated random double precision reals!", &
       filename, __LINE__)

!$omp end parallel
CONTAINS
  SUBROUTINE reduce_p(p, p_all, msg, file, line)
    LOGICAL, INTENT(in) :: p
    LOGICAL, INTENT(inout) :: p_all
    CHARACTER(len=*), INTENT(in) :: msg, file
    INTEGER, INTENT(in) :: line
    IF (p) THEN
!$omp critical
      p_all = .TRUE.
!$omp end critical
    END IF
!$omp barrier
!$omp master
    IF (p_all) THEN
      CALL abort_ppm(msg, file, line, ppm_default_comm)
    END IF
!$omp end master
  END SUBROUTINE reduce_p
END PROGRAM test_irand
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
