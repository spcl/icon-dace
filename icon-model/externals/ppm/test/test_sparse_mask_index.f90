!>
!! @file test_sparse_mask_index.f90
!! @brief test sparse mask index
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: sparse mask index test
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
MODULE sparse_mask_index_tests
  USE ppm_base, ONLY: abort_ppm
  USE ppm_sparse_mask_index, ONLY: index_sparse_nd
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_smi_equiv_1d, test_smi_equiv_2d, test_smi_equiv_3d
  CHARACTER(len=*), PARAMETER :: filename = 'test_sparse_mask_index.f90'
CONTAINS
  SUBROUTINE test_smi_equiv_1d(idx, strt, ends, mask_1d, accum_1d, accum_ref_1d)
    INTEGER, INTENT(in) :: strt(1), ends(1)
    LOGICAL, INTENT(in) :: mask_1d(strt(1):ends(1))
    REAL, DIMENSION(strt(1):ends(1)), INTENT(inout) :: accum_1d, accum_ref_1d
    TYPE(index_sparse_nd), INTENT(in) :: idx
    INTEGER :: i, n, r1
    LOGICAL :: p

    ! since older intel compilers cannot parallelize the above where correctly,
    ! use an explicit DO
!$omp do
    DO i = strt(1), ends(1)
      accum_ref_1d(i) = accum_ref_1d(i) + MERGE(1.0, 0.0, mask_1d(i))
    END DO
!$omp end do

    n = idx%tl_ranges
!$omp do
    DO r1 = 1, n
      DO i = idx%ranges(r1)%first, idx%ranges(r1)%last
        accum_1d(i) = accum_1d(i) + 1.0
      END DO
    END DO
!$omp end do
    p = .TRUE.
    n = ends(1)
!$omp do
    DO i = strt(1), n
      p = p .AND. ABS(accum_1d(i) - accum_ref_1d(i)) <= 0.0
    END DO
!$omp end do
    IF (.NOT. p) CALL abort_ppm(source=filename, line=__LINE__, &
         msg='result of WHERE(mask) and index iteration differs')
  END SUBROUTINE test_smi_equiv_1d

  SUBROUTINE where_add_2d(strt, ends, mask, a, y)
    INTEGER, INTENT(in) :: strt(2), ends(2)
    LOGICAL, INTENT(in) :: mask(strt(1):ends(1), strt(2):ends(2))
    REAL, INTENT(inout) :: a(strt(1):ends(1), strt(2):ends(2))
    REAL, INTENT(in) :: y
    INTEGER :: i, j
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
        IF (mask(i, j)) a(i, j) = a(i, j) + y
      END DO
    END DO
!$omp end do nowait

  END SUBROUTINE where_add_2d

  SUBROUTINE test_smi_equiv_2d(idx, strt, ends, sseq, &
       mask_2d, accum_2d, accum_ref_2d)
    INTEGER, INTENT(in) :: strt(2), ends(2), sseq(2)
    LOGICAL, INTENT(in) :: mask_2d(strt(1):ends(1), strt(2):ends(2))
    REAL, DIMENSION(strt(1):ends(1), strt(2):ends(2)), INTENT(inout) :: &
         accum_2d, accum_ref_2d
    TYPE(index_sparse_nd), INTENT(in) :: idx

    INTEGER :: i, j, n, r1, r2, iblock
    LOGICAL :: p

    CALL where_add_2d(strt, ends, mask_2d, accum_ref_2d, 1.0)

    n = 2 * idx%tl_ranges
    IF (sseq(1) == 1) THEN
!$omp do
      DO r1 = 1, n, 2
        DO j = idx%ranges(r1)%first, idx%ranges(r1)%last
          r2 = idx%ranges(r1+1)%first + j - idx%ranges(r1)%first
          DO iblock = idx%ranges(r2)%first, idx%ranges(r2)%last
            DO i = idx%ranges(iblock)%first, idx%ranges(iblock)%last
              accum_2d(i, j) = accum_2d(i, j) + 1.0
            END DO
          END DO
        END DO
      END DO
!$omp end do
    ELSE ! sseq(1) == 2
!$omp do
      DO r1 = 1, n, 2
        DO i = idx%ranges(r1)%first, idx%ranges(r1)%last
          r2 = idx%ranges(r1+1)%first + i - idx%ranges(r1)%first
          DO iblock = idx%ranges(r2)%first, idx%ranges(r2)%last
            DO j = idx%ranges(iblock)%first, idx%ranges(iblock)%last
              accum_2d(i, j) = accum_2d(i, j) + 1.0
            END DO
          END DO
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
!$omp end do
    IF (.NOT. p) CALL abort_ppm(source=filename, line=__LINE__, &
         msg='result of mask and index iteration differs')
  END SUBROUTINE test_smi_equiv_2d

  SUBROUTINE test_smi_equiv_3d(idx, strt, ends, sseq, mask, accum, accum_ref)
    INTEGER, INTENT(in) :: strt(3), ends(3), sseq(3)
    LOGICAL, INTENT(in) :: mask(strt(1):ends(1),strt(2):ends(2),strt(3):ends(3))
    REAL, DIMENSION(strt(1):ends(1), strt(2):ends(2),strt(3):ends(3)), &
         INTENT(inout) :: accum, accum_ref
    TYPE(index_sparse_nd), INTENT(in) :: idx

    INTEGER :: i, j, k, r1, r2, r3, n, jblock, iblock, indices(3)
    LOGICAL :: p
#if defined _CRAYFTN
#  if _RELEASE_MAJOR == 10 \
    || _RELEASE_MAJOR == 11 && _RELEASE_MINOR == 0 && _RELEASE_PATCHLEVEL < 4 \
    || _RELEASE_MAJOR == 12 && _RELEASE_MINOR == 0 && _RELEASE_PATCHLEVEL == 3 \
    || _RELEASE_MAJOR == 13 && _RELEASE_MINOR == 0 && _RELEASE_PATCHLEVEL == 1
    ! In non-OpenMP mode, crayftn 10 and 11 botch this function (incorrect
    ! index computations when strt has large absolute value).
    ! The bug was fixed in 11.0.4, but was found again in 12.0.3
    ! and 13.0.1. crayftn 14.0.0 seems to be fixed, but so far
    ! versions in between those mentioned were not tested
!DIR$ OPTIMIZE (-O0)
#  elif _RELEASE_MAJOR == 12 || _RELEASE_MINOR == 13
#    warning "suspect, but not yet tested version of Cray ftn detected."
#  endif
#endif

    DO k = strt(3), ends(3)
      CALL where_add_2d(strt, ends, mask(:, :, k), accum_ref(:, :, k), 1.0)
    END DO

    n = 2 * idx%tl_ranges
!$omp do
    DO r1 = 1, n, 2
      DO k = idx%ranges(r1)%first, idx%ranges(r1)%last
        indices(sseq(3)) = k
        r2 = idx%ranges(r1+1)%first + k - idx%ranges(r1)%first
        DO jblock = idx%ranges(r2)%first, idx%ranges(r2)%last, 2
          DO j = idx%ranges(jblock)%first, idx%ranges(jblock)%last
            indices(sseq(2)) = j
            r3 = idx%ranges(jblock+1)%first + j - idx%ranges(jblock)%first
            DO iblock = idx%ranges(r3)%first, idx%ranges(r3)%last
              DO i = idx%ranges(iblock)%first, idx%ranges(iblock)%last
                indices(sseq(1)) = i
                accum(indices(1), indices(2), indices(3)) &
                     = accum(indices(1), indices(2), indices(3)) + 1.0
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do
    p = .TRUE.
    n = ends(3)
!$omp do
    DO k = strt(3), n
      p = p .AND. ALL(ABS(accum(:, :, k) - accum_ref(:, :, k)) <= 0.0)
    END DO
!$omp end do nowait
    IF (.NOT. p) CALL abort_ppm(source=filename, line=__LINE__, &
         msg='result of mask and index iteration differs')
  END SUBROUTINE test_smi_equiv_3d
END MODULE sparse_mask_index_tests

PROGRAM test_sparse_mask_index
  USE ppm_base, ONLY: ppm_default_comm, assertion, abort_ppm
  USE ppm_combinatorics, ONLY: permute
  USE ppm_extents, ONLY: extent, extent_start, extent_end, iinterval
  USE ppm_random, ONLY: irandr, lrand, initialize_irand, a_randp_mt, frandp
  USE ppm_sparse_mask_index, ONLY: index_sparse_nd, index_from_mask
  USE random_data, ONLY: random_rectilinear
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  USE sparse_mask_index_tests, ONLY: test_smi_equiv_1d, test_smi_equiv_2d, &
       test_smi_equiv_3d
  IMPLICIT NONE
  TYPE(index_sparse_nd) :: idx
  INTEGER :: n, sseq(3), strt(3), ends(3)
  INTEGER, PARAMETER :: max_iterations = 1000
  TYPE(iinterval), PARAMETER :: mask_size_range=iinterval(0, 400000)
  TYPE(extent) :: mask_dims(3)
  INTEGER :: mask_size_ub
  INTEGER :: allocation_error, rand_seed, tid, i, j, k
  LOGICAL, ALLOCATABLE :: mask(:, :,:)
  REAL, ALLOCATABLE :: accum(:, :, :), accum_ref(:, :, :)
#if defined __INTEL_COMPILER && __INTEL_COMPILER >= 1400 && __INTEL_COMPILER <= 1600 && defined _OPENMP && ! defined __OPTIMIZE__
#define NEED_LB_WORKAROUND
  TARGET :: accum
  REAL, POINTER :: accum_p(:, :, :)
#endif
  REAL :: mask_fraction
  CHARACTER(len=132) :: msg
  CHARACTER(len=*), PARAMETER :: filename = 'test_sparse_mask_index.f90'

!$omp parallel &
!$omp private(tid, rand_seed, n, mask_size_ub, allocation_error, &
!$omp & strt, ends) &
#ifdef NEED_LB_WORKAROUND
!$omp & private(accum_p) &
#endif
!$omp & shared(accum, accum_ref, mask, ppm_default_comm, mask_dims, &
!$omp & mask_fraction, idx, sseq)
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed
!$omp master
  n = 0
  allocation_loop: DO WHILE(.NOT. ALLOCATED(mask) .AND. n < max_iterations)
    mask_size_ub = irandr(mask_size_range)
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
  ! status
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
  ! ensure mask is sparse
  mask_fraction = 0.75 + 0.25 * frandp()
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
  ! 1D-case
!$omp single
  CALL index_from_mask(idx, mask(:, strt(2), strt(3)), offsets=strt(1:1))
!$omp end single
  CALL test_smi_equiv_1d(idx, strt, ends, mask, accum, accum_ref)

  ! 2D
!$omp single
  ! g95 incorrectly evaluates array merges when the mask argument is a scalar
  sseq(1) = MERGE(1, 2, lrand())
  sseq(2) = 3 - sseq(1)
  CALL index_from_mask(idx, mask(:, :, strt(3)), sseq=sseq(1:2), &
       offsets=strt(1:2))
!$omp end single
  CALL test_smi_equiv_2d(idx, strt, ends, sseq, mask, accum, accum_ref)
  ! 3D
!$omp single
  sseq = (/ 1, 2, 3 /)
  CALL permute(sseq)
  CALL index_from_mask(idx, mask, sseq=sseq, offsets=strt)
!$omp end single
  CALL test_smi_equiv_3d(idx, strt, ends, sseq, mask, accum, accum_ref)
!$omp end parallel
END PROGRAM test_sparse_mask_index
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
