!>
!! @file test_ddp_math.f90
!! @brief Verify that the DDP corrected summation yields satisfactory results
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
PROGRAM test_ddp_math
  USE ppm_base, ONLY: assertion, ppm_default_comm
  USE ppm_std_type_kinds, ONLY: dp
  USE ppm_combinatorics, ONLY: permute
  USE ppm_random, ONLY: initialize_irand, drandp, irandp
  USE ppm_math_extensions, ONLY: ddp_sum
  IMPLICIT NONE
  INTEGER, PARAMETER :: max_num_elems = 30000, min_num_elems = 4
  INTEGER :: rand_seed, i, j, n, step, pos(4)
  REAL(dp), ALLOCATABLE :: a(:)
  REAL(dp) :: t1, t2, sum_ddp, sum_naive
  CHARACTER(len=*), PARAMETER :: filename = 'test_ddp_math.f90'


  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  CALL assertion(RADIX(a(1)) == 2, filename, __LINE__, 'binary math assumption not true!')

  n = MOD(irandp(), max_num_elems - min_num_elems + 1) + min_num_elems
  n = IAND(n + 3, NOT(3))
  ALLOCATE(a(n+5))
  ! first fill part of array with reals cancelling in groups of 4
  DO i = 0, n-4, 4
    t2 = drandp() * 0.5_dp + 0.5_dp
    a(i + 1) = t2
    t1 = drandp()
    t1 = MERGE(t1, t1 + SPACING(t1), MODULO(t1, SPACING(t1)) >= SPACING(t1))
    ! FIXME: how do we get binary digits in mantissa in Fortran?
    !        it's 53 for double precision, but where do we find this
    t1 = SCALE(t1, -50)
    a(i + 2) = t1
    a(i + 3) = -t2
    a(i + 4) = -t1
  END DO
  ! then append cancelling group of 5 reals,
  ! to prevent power-of-2 vector registers from generating precise sums
  t2 = drandp() * 0.5_dp + 0.5_dp
  a(n+1) = t2
  t1 = drandp()
  t1 = MERGE(t1, t1 + SPACING(t1), MODULO(t1, SPACING(t1)) >= SPACING(t1))
  ! FIXME: how do we get binary digits in mantissa in Fortran?
  !        it's 53 for double precision, but where do we find this
  t1 = SCALE(t1, -50)
  a(n+2) = 0.0_dp
  a(n+3) = t1
  a(n+4) = -t2
  a(n+5) = -t1
  sum_naive = SUM(a)
  sum_ddp = REAL(ddp_sum(a), dp)
  PRINT *, 'simple test, naive sum:          ', sum_naive
  PRINT *, 'simple test, DDP   sum:          ', sum_ddp

  CALL assertion(sum_naive > 0.0_dp .OR. sum_naive < 0.0_dp, &
       filename, __LINE__, &
       'Naive summation gave correct result, did you configure the FPU correctly?')
  CALL assertion(ABS(sum_ddp) <= 0.0_dp, filename, __LINE__, &
       'Kahan summation failed, did you configure the FPU correctly?')


  ! next distribute reals from random sample,
  ! but form groups of similar magnitude
  step = n / 4
  DO i = 0, step-1
    DO j = 1, 4
      pos(j) = MOD(i + step * (j-1), n) + 1
    END DO
    t2 = drandp()
    a(pos(1)) = t2
    t1 = drandp() * 0.5_dp + 0.5_dp
    t1 = SCALE(t1, -40)
    a(pos(2)) = t1
    a(pos(3)) = -t2
    a(pos(4)) = -t1
  END DO
  sum_naive = SUM(a)
  sum_ddp = REAL(ddp_sum(a), dp)
  PRINT *, 'staggered test, naive sum:       ', sum_naive
  PRINT *, 'staggered test, DDP   sum:       ', sum_ddp

  CALL assertion(ABS(sum_ddp) <= 0.0_dp, filename, __LINE__, &
       'Kahan summation failed, have you' &
       // ' configured the FPU correctly?')

  ! next distribute reals from random sample,
  ! do no longer group but permute randomly
  step = n / 4
  DO i = 0, step-1
    DO j = 1, 4
      pos(j) = MOD(i + step * (j-1), n) + 1
    END DO
    CALL permute(pos)
    t2 = drandp()
    a(pos(1)) = t2
    t1 = drandp() * 0.5_dp + 0.5_dp
    t1 = SCALE(t1, -40)
    a(pos(2)) = t1
    a(pos(3)) = -t2
    a(pos(4)) = -t1
  END DO
  sum_naive = SUM(a)
  sum_ddp = REAL(ddp_sum(a), dp)
  PRINT *, 'random sequence test, naive sum: ', sum_naive
  PRINT *, 'random sequence test, DDP   sum: ', sum_ddp

  CALL assertion(ABS(sum_ddp) <= 0.0_dp, &
       filename, &
       __LINE__, 'Kahan summation failed, have you' &
       // ' configured the FPU correctly?')

END PROGRAM test_ddp_math
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
