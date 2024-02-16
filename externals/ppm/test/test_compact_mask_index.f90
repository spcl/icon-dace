!>
!! @file test_compact_mask_index.f90
!! @brief build and test use of mask index data structure
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
MODULE compact_mask_index_tests
  USE ppm_base, ONLY: assertion
  USE ppm_compact_mask_index, ONLY: range_compact_2d, range_compact_3d
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_cmi_equiv_2d, test_cmi_equiv_3d
  CHARACTER(len=*), PARAMETER :: filename = 'test_compact_mask_index.f90'
CONTAINS
  SUBROUTINE test_cmi_equiv_2d(num_ranges, idx2d, strt, ends, sseq, &
       mask_2d, accum_2d, accum_ref_2d)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: num_ranges, strt(2), ends(2), sseq(2)
    LOGICAL, INTENT(in) :: mask_2d(strt(1):ends(1), strt(2):ends(2))
    REAL, DIMENSION(strt(1):ends(1), strt(2):ends(2)), INTENT(inout) :: &
         accum_2d, accum_ref_2d
    TYPE(range_compact_2d), INTENT(in) :: idx2d(num_ranges)

    INTEGER :: i, j, n, r1
    LOGICAL :: p
#if defined _CRAYFTN && ! defined _OPENMP
#  if _RELEASE_MAJOR >= 10 && _RELEASE_MAJOR <= 12
    ! in non-OpenMP mode, crayftn 10 to 12 botch this function (incorrect
    ! index computations when strt has large absolute value
!DIR$ OPTIMIZE (-O0)
#  endif
#endif

!$omp do
    DO j = strt(2), ends(2)
      DO i = strt(1), ends(1)
        accum_ref_2d(i, j) = accum_ref_2d(i, j) + MERGE(1.0, 0.0, mask_2d(i, j))
      END DO
    END DO
!$omp end do

    n = num_ranges
    IF (sseq(1) == 1) THEN
!$omp do
      DO r1 = 1, n
        j = idx2d(r1)%tl_ss
        DO i = idx2d(r1)%ll_range%first, idx2d(r1)%ll_range%last
          accum_2d(i, j) = accum_2d(i, j) + 1.0
        END DO
      END DO
!$omp end do
    ELSE ! sseq(1) == 2
!$omp do
      DO r1 = 1, n
        i = idx2d(r1)%tl_ss
        DO j = idx2d(r1)%ll_range%first, idx2d(r1)%ll_range%last
          accum_2d(i, j) = accum_2d(i, j) + 1.0
        END DO
      END DO
!$omp end do
    END IF
    p = .TRUE.
    n = ends(2)
!$omp do
    DO j = strt(2), n
      p = p .AND. ALL(ABS(accum_2d(:, j) - accum_ref_2d(:, j)) <= 0.0)
    END DO
!$omp end do nowait
    CALL assertion(p, filename, __LINE__, &
         'result of WHERE(mask) and index iteration differs')
  END SUBROUTINE test_cmi_equiv_2d

  SUBROUTINE test_cmi_equiv_3d(num_ranges, idx3d, strt, ends, sseq, &
       mask, accum, accum_ref)
    INTEGER, INTENT(in) :: num_ranges, strt(3), ends(3), sseq(3)
    TYPE(range_compact_3d), INTENT(in) :: idx3d(num_ranges)
    LOGICAL, INTENT(in) :: &
         mask(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3))
    REAL, DIMENSION(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3)), &
         INTENT(inout) :: accum, accum_ref
    INTEGER :: i, j, k, r1, indices(3)
    LOGICAL :: p
#if defined _CRAYFTN
#  if _RELEASE_MAJOR == 10 \
    || _RELEASE_MAJOR == 11 && _RELEASE_MINOR == 0 && _RELEASE_PATCHLEVEL < 4
    ! crayftn 10 and 11 botch this function (incorrect
    ! index computations when strt has large absolute value
!DIR$ OPTIMIZE (-O0)
#  endif
#endif
    ! since older intel compilers cannot parallelize WHERE correctly,
    ! use an explicit DO
!$omp do
    DO k = strt(3), ends(3)
      DO j = strt(2), ends(2)
        DO i = strt(1), ends(1)
          accum_ref(i, j, k) = accum_ref(i, j, k) &
               + MERGE(1.0, 0.0, mask(i, j, k))
        END DO
      END DO
    END DO
!$omp end do

!$omp do
    DO r1 = 1, num_ranges
      j = idx3d(r1)%tl_ss(1)
      k = idx3d(r1)%tl_ss(2)
      indices(sseq(2)) = j
      indices(sseq(3)) = k
      DO i = idx3d(r1)%ll_range%first, idx3d(r1)%ll_range%last
        ! do not allocate temporary
        indices(sseq(1)) = i
        accum(indices(1), indices(2), indices(3)) &
             = accum(indices(1), indices(2), indices(3)) + 1.0
      END DO
    END DO
!$omp end do
  p = .TRUE.
!$omp do
  DO k = strt(3), ends(3)
    p = p .AND. ALL(ABS(accum(:, :, k) - accum_ref(:, :, k)) <= 0.0)
  END DO
!$omp end do nowait
  CALL assertion(p, filename, __LINE__, &
       'result of WHERE(mask) and index iteration differs')

  END SUBROUTINE test_cmi_equiv_3d
END MODULE compact_mask_index_tests

PROGRAM test_compact_mask_index
  USE ppm_base, ONLY: ppm_default_comm, abort_ppm
  USE ppm_random, ONLY: irandr, initialize_irand, a_randp_mt, &
       frandp, lrand
  USE ppm_compact_mask_index, ONLY: range_compact_2d, range_compact_3d, &
       index_from_mask, index_from_mask_mt
  USE ppm_combinatorics, ONLY: permute
  USE ppm_extents, ONLY: extent, extent_start, extent_end, iinterval
  USE random_data, ONLY: random_rectilinear
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  USE compact_mask_index_tests, ONLY: test_cmi_equiv_2d, &
       test_cmi_equiv_3d
  IMPLICIT NONE
  TYPE(range_compact_2d), ALLOCATABLE :: idx2d(:)
  TYPE(range_compact_3d), ALLOCATABLE :: idx3d(:)

  INTEGER :: n, sseq(3), rand_seed, tid, strt(3), ends(3)
  INTEGER, PARAMETER :: max_iterations = 1000
  TYPE(iinterval), PARAMETER :: mask_size_max=iinterval(0, 40000000)
  TYPE(extent) :: mask_dims(3)
  INTEGER :: allocation_error, mask_size_ub, i, j, k
  LOGICAL, ALLOCATABLE :: mask(:, :, :)
  REAL, ALLOCATABLE :: accum(:, :, :), accum_ref(:, :, :)
#if defined __INTEL_COMPILER && __INTEL_COMPILER >= 1400 && __INTEL_COMPILER <= 1600 && defined _OPENMP && ! defined __OPTIMIZE__
#define NEED_LB_WORKAROUND
  TARGET :: accum
  REAL, POINTER :: accum_p(:, :, :)
#endif
  REAL :: mask_fraction
  CHARACTER(len=132) :: msg
  CHARACTER(len=*), PARAMETER :: filename = 'test_compact_mask_index.f90'

!$omp parallel private(tid, n, rand_seed, mask_size_ub, &
!$omp & allocation_error, strt, ends) &
#ifdef NEED_LB_WORKAROUND
!$omp & private(accum_p) &
#endif
!$omp shared(accum, accum_ref, idx2d, idx3d, mask, mask_fraction, &
!$omp        mask_dims, ppm_default_comm, sseq)
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed
!$omp master
  n = 0
  allocation_loop: DO WHILE(.NOT. ALLOCATED(mask) .AND. n <max_iterations)
    mask_size_ub = irandr(mask_size_max)
    CALL random_rectilinear(mask_dims, mask_size_ub)
#ifdef __G95__
    ! g95 cannot handle arbitrary lower array bounds
    mask_dims%first = 1
#endif
    strt = extent_start(mask_dims)
    ends = extent_end(mask_dims)
    ALLOCATE(mask(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3)), &
         accum(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3)), &
         accum_ref(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3)), &
         stat=allocation_error)
    n = n + 1
  END DO allocation_loop
  IF (.NOT. ALLOCATED(mask)) THEN
    WRITE (msg, '(a,i0)') &
         'failed to generate allocatable mask, allocation_error=', &
         allocation_error
    CALL abort_ppm(msg, filename, __LINE__)
  END IF
  ! the following is redundant, but needed for compiler analysis of allocation
  IF (.NOT. ALLOCATED(accum_ref)) STOP 1
!$omp end master
!$omp barrier
  strt = LBOUND(mask)
  ends = UBOUND(mask)
#ifdef NEED_LB_WORKAROUND
  ! this work-around is needed because ifort 14.x-16.x with flags -openmp -O0 -g
  ! might loose the lbounds of an array when calling a_randp_mt
  accum_p => accum
  CALL a_randp_mt(accum_p)
#else
  CALL a_randp_mt(accum)
#endif
!$omp master
  ! ensure mask is densely-populated
  mask_fraction = 0.25 * frandp()
!$omp end master
!$omp barrier
!$omp do
  DO k = strt(3), ends(3)
    DO j = strt(2), ends(2)
      DO i = strt(1), ends(1)
        accum_ref(i, j, k) = accum(i, j, k)
        mask(i, j, k) = accum(i, j, k) >= mask_fraction
      END DO
    END DO
  END DO
!$omp end do
  ! 2D
!$omp single
  ! g95 incorrectly evaluates array merges when the mask argument is a scalar
  sseq(1) = MERGE(2, 1, lrand())
  sseq(2) = 3 - sseq(1)
  CALL index_from_mask(idx2d, mask(:, :, strt(3)), offsets=strt(1:2), &
       sseq=sseq(1:2))
!$omp end single
  CALL test_cmi_equiv_2d(SIZE(idx2d), idx2d, strt, ends, sseq, &
       mask, accum, accum_ref)
!$omp end parallel
  ! 3D
  sseq = (/ 1, 2, 3 /)
  CALL permute(sseq)
  strt = LBOUND(mask)
  ends = UBOUND(mask)
  CALL index_from_mask_mt(idx3d, mask(:,:,:), offsets=strt, sseq=sseq)
!$omp parallel firstprivate(sseq, strt, ends) &
!$omp shared(accum, accum_ref, idx2d, idx3d, mask, mask_fraction, &
!$omp        mask_dims, ppm_default_comm)

  CALL test_cmi_equiv_3d(SIZE(idx3d), idx3d, strt, ends, sseq, &
       mask, accum, accum_ref)
!$omp end parallel
END PROGRAM test_compact_mask_index
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
