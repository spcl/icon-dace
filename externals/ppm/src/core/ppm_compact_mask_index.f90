!> @file ppm_compact_mask_index.f90
!! @brief compute index from densely-populated mask
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords:
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
!> ppm_compact_mask_index.f90
!! supposing you have a mask array and use it in something like a
!! WHERE, then it can be beneficial to only visit those indices of the
!! array expression with the condition evaluating to .TRUE., in case
!! this is relatively sparse or the ensuing conditional is expensive
!! in comparison to the WHERE/FORALL body
!!
!! This index variation addresses the case where each range in the
!! subscript varied most rapidly (ll) corresponds to a change in at least
!! one slower (hl) subscript.

MODULE ppm_compact_mask_index
  USE ppm_extents, ONLY: iinterval
  USE ppm_sparse_mask_index, ONLY: build_ranges_1d, count_ranges_1d

  IMPLICIT NONE
  PRIVATE
  !> describe contiguous range in 2d mask
  TYPE range_compact_2d
    INTEGER :: tl_ss              !< top-level subscript
    TYPE(iinterval) :: ll_range   !< range of lower-level subscripts
  END TYPE range_compact_2d

  TYPE range_compact_3d
    INTEGER :: tl_ss(2)           !< top-level subscript
    TYPE(iinterval) :: ll_range   !< range of lower-level subscripts
  END TYPE range_compact_3d

  !> compute index from mask array
  INTERFACE index_from_mask
    MODULE PROCEDURE index_from_mask_2d
    MODULE PROCEDURE index_from_mask_3d
  END INTERFACE index_from_mask

  !> compute index from mask array, multi-threaded variant
  INTERFACE index_from_mask_mt
    ! MODULE PROCEDURE index_from_mask_mt_2d
    MODULE PROCEDURE index_from_mask_mt_3d
  END INTERFACE index_from_mask_mt

  PUBLIC :: index_from_mask, index_from_mask_mt, &
       range_compact_2d, range_compact_3d

CONTAINS
  PURE SUBROUTINE build_ranges_2d(mask, offsets, sseq, idx)
    LOGICAL, INTENT(in) :: mask(:, :)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(2)
    INTEGER, INTENT(in) :: sseq(2)
    TYPE(range_compact_2d), INTENT(out) :: idx(:)

    INTEGER :: i, j, m, n, p, q
    INTEGER :: ofs(2)

    ofs = 1
    IF (PRESENT(offsets)) ofs = offsets

    p = 1
    IF (sseq(2) == 2) THEN
      n = SIZE(mask, 2)
      DO j = 1, n
        CALL build_ranges_1d(mask(:, j), ofs(1), idx%ll_range, p, rcount=q)
        idx(p:p+q-1)%tl_ss = j + ofs(2) - 1
        p = p + q
      END DO
    ELSE ! sseq(2) == 1
      m = SIZE(mask, 1)
      DO i = 1, m
        CALL build_ranges_1d(mask(i, :), ofs(2), idx%ll_range, p, rcount=q)
        idx(p:p+q-1)%tl_ss = i + ofs(1) - 1
        p = p + q
      END DO
    END IF
  END SUBROUTINE build_ranges_2d

  PURE FUNCTION count_ranges_2d(mask, sseq) RESULT(num_ranges)
    LOGICAL, INTENT(in) :: mask(:, :)
    INTEGER, INTENT(in) :: sseq(2)

    INTEGER :: i, j, m, n, num_ranges

    num_ranges = 0
    IF (sseq(2) == 2) THEN
      n = SIZE(mask, 2)
      DO j = 1, n
        num_ranges = num_ranges + count_ranges_1d(mask(:, j))
      END DO
    ELSE ! sseq(2) == 1
      m = SIZE(mask, 1)
      DO i = 1, m
        num_ranges = num_ranges + count_ranges_1d(mask(i, :))
      END DO
    END IF
  END FUNCTION count_ranges_2d

  !> compute index from mask array
  !! @param idx index to construct
  !! @param mask to be indexed
  !! @param offsets low bounds of mask (if not present assumed to be 1)
  !! @param sseq sequence in which to travel indices,
  !! e.g. if mask is transposed, one would pass (/ 2, 1 /) to the routine
  PURE SUBROUTINE index_from_mask_2d(idx, mask, offsets, sseq)
    TYPE(range_compact_2d), ALLOCATABLE, INTENT(out) :: idx(:)
    LOGICAL, INTENT(in) :: mask(:, :)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(2), sseq(2)

    INTEGER :: num_ranges, subscript_sequence(2)

    subscript_sequence = (/ 1, 2 /)
    IF (PRESENT(sseq)) subscript_sequence = sseq

    num_ranges = count_ranges_2d(mask, subscript_sequence)
    ALLOCATE(idx(num_ranges))
    CALL build_ranges_2d(mask, offsets, subscript_sequence, idx)
  END SUBROUTINE index_from_mask_2d

  FUNCTION build_ranges_3d(mask, sizes, strides, sseq, offsets, idx) &
       RESULT(num_ranges)
    LOGICAL, INTENT(in) :: mask(*)
    INTEGER, INTENT(in) :: sizes(3), strides(3), sseq(3)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(3)
    TYPE(range_compact_3d), OPTIONAL, TARGET, INTENT(out) :: idx(:)

    TYPE(iinterval), POINTER :: d_ranges(:)
    INTEGER :: i, j, k, m, n, o, p, q, num_ranges
    INTEGER :: ofs(3), base
    INTEGER, ALLOCATABLE :: rcounts(:,:), d_ranges_ofs(:,:)

    ofs = 1
    IF (PRESENT(offsets)) ofs = offsets

    m = sizes(1)
    n = sizes(2)
    o = sizes(3)

    IF (PRESENT(idx)) THEN
      ALLOCATE(d_ranges_ofs(n, o), rcounts(n, o))
      DO k = 1, o
        DO j = 1, n
          base = (k-1)*strides(3) + (j - 1) * strides(2)
          rcounts(j, k) &
               = count_ranges_1d(mask(base+1:base+m*strides(1):strides(1)))
        END DO
      END DO
      p = 1
      DO k = 1, o
        DO j = 1, n
          d_ranges_ofs(j, k) = p
          p = p + rcounts(j, k)
        END DO
      END DO
      DO k = 1, o
        DO j = 1, n
          IF (rcounts(j, k) > 0) THEN
            p = d_ranges_ofs(j, k)
            q = p + rcounts(j, k) - 1
            d_ranges => idx(p:q)%ll_range
            base = (k-1)*strides(3) + (j - 1) * strides(2)
            CALL build_ranges_1d(ofs=ofs(sseq(1)), &
                 mask=mask(base+1:base+m*strides(1):strides(1)), &
                 d_ranges=d_ranges, rcount=q)
            DO i = p, p + q - 1
              idx(i)%tl_ss(1) = j + ofs(sseq(2)) - 1
              idx(i)%tl_ss(2) = k + ofs(sseq(3)) - 1
            END DO
          END IF
        END DO
      END DO
      num_ranges = d_ranges_ofs(n, o) + rcounts(n, o) - 1
      DEALLOCATE(d_ranges_ofs, rcounts)
    ELSE
      num_ranges = 0
      DO k = 1, o
        DO j = 1, n
          base = (k-1)*strides(3) + (j - 1) * strides(2)
          num_ranges = num_ranges &
               + count_ranges_1d(mask(base+1:base+m*strides(1):strides(1)))
        END DO
      END DO
    END IF
  END FUNCTION build_ranges_3d

  !> construct compact index from mask
  !! @param idx index to construct
  !! @param mask to be indexed
  !! @param offsets low bounds of mask (if not present assumed to be 1)
  !! @param sseq subscript sequence i.e. sequence by which to store
  !! ranges in index
  SUBROUTINE index_from_mask_3d(idx, mask, offsets, sseq)
    TYPE(range_compact_3d), ALLOCATABLE, INTENT(out) :: idx(:)
    LOGICAL, INTENT(in) :: mask(:, :, :)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(3), sseq(3)

    INTEGER :: num_ranges, sseq_(3), strides(3), sizes(3), temp(3), &
         i, accum, sseqi(3)

    sizes = SHAPE(mask)

    IF (PRESENT(sseq)) THEN
      DO i = 1, 3
        sseqi(sseq(i)) = i
      END DO
      temp = sizes
      accum = 1
      DO i = 1, 3
        sizes(i) = temp(sseq(i))
        strides(sseqi(i)) = accum
        accum = accum * temp(i)
      END DO
      sseq_ = sseq
    ELSE
      sseq_ = (/ 1, 2, 3 /)
      strides(1) = 1
      strides(2) = sizes(1)
      strides(3) = sizes(1) * sizes(2)
    END IF

    num_ranges = build_ranges_3d(mask, sizes, strides, sseq_, offsets)
    ALLOCATE(idx(num_ranges))
    num_ranges = build_ranges_3d(mask, sizes, strides, sseq_, offsets, idx)
  END SUBROUTINE index_from_mask_3d

  FUNCTION build_ranges_mt_3d(mask, sizes, strides, sseq, offsets, idx) &
       RESULT(num_ranges)
    LOGICAL, INTENT(in) :: mask(*)
    INTEGER, INTENT(in) :: sizes(3), strides(3), sseq(3)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(3)
    TYPE(range_compact_3d), OPTIONAL, TARGET, INTENT(out) :: idx(:)

    TYPE(iinterval), POINTER :: d_ranges(:)
    INTEGER :: i, j, k, m, n, o, p, q, num_ranges
    INTEGER :: ofs(3), base
    INTEGER, ALLOCATABLE :: rcounts(:,:), d_ranges_ofs(:,:)

    ofs = 1
    IF (PRESENT(offsets)) ofs = offsets

    num_ranges = 0

!$omp parallel shared(d_ranges_ofs, idx, mask, rcounts) &
!$omp private(d_ranges, i, j, k, m, n, o, p, q, base) &
!$omp firstprivate(ofs, sseq, sizes, strides) reduction(+: num_ranges)
    m = sizes(1)
    n = sizes(2)
    o = sizes(3)
    IF (PRESENT(idx)) THEN
!$omp single
      ALLOCATE(rcounts(n, o))
!$omp end single
!$omp do
      DO k = 1, o
        DO j = 1, n
          base = (k-1)*strides(3) + (j - 1) * strides(2)
          rcounts(j, k) &
               = count_ranges_1d(mask(base+1:base+m*strides(1):strides(1)))
        END DO
      END DO
!$omp end do
!$omp single
      ALLOCATE(d_ranges_ofs(n, o))
      p = 1
      DO k = 1, o
        DO j = 1, n
          d_ranges_ofs(j, k) = p
          p = p + rcounts(j, k)
        END DO
      END DO
!$omp end single
!$omp do
      DO k = 1, o
        DO j = 1, n
          IF (rcounts(j, k) > 0) THEN
            p = d_ranges_ofs(j, k)
            q = p + rcounts(j, k) - 1
            d_ranges => idx(p:q)%ll_range
            base = (k-1)*strides(3) + (j - 1) * strides(2)
            CALL build_ranges_1d(ofs=ofs(sseq(1)), &
                 mask=mask(base+1:base+m*strides(1):strides(1)), &
                 d_ranges = d_ranges, rcount=q)
            DO i = p, p + q - 1
              idx(i)%tl_ss(1) = j + ofs(sseq(2)) - 1
              idx(i)%tl_ss(2) = k + ofs(sseq(3)) - 1
            END DO
          END IF
        END DO
      END DO
!$omp end do
!$omp master
      num_ranges = num_ranges + (d_ranges_ofs(n, o) + rcounts(n, o) - 1)
      DEALLOCATE(d_ranges_ofs, rcounts)
!$omp end master
    ELSE
!$omp do
      DO k = 1, o
        DO j = 1, n
          base = (k-1)*strides(3) + (j - 1) * strides(2)
          num_ranges = num_ranges &
               + count_ranges_1d(mask(base+1:base+m*strides(1):strides(1)))
        END DO
      END DO
!$omp end do
    END IF
!$omp end parallel
  END FUNCTION build_ranges_mt_3d

  !> construct compact index from mask
  !!
  !! this version must be called from outside an OpenMP parallel
  !! region because it will open its own region
  !!
  !! @param idx index to construct
  !! @param mask to be indexed
  !! @param offsets low bounds of mask (if not present assumed to be 1)
  !! @param sseq subscript sequence i.e. sequence by which to store
  !! ranges in index
  SUBROUTINE index_from_mask_mt_3d(idx, mask, offsets, sseq)
    TYPE(range_compact_3d), ALLOCATABLE, INTENT(out) :: idx(:)
    LOGICAL, INTENT(in) :: mask(:, :, :)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(3), sseq(3)

    INTEGER :: num_ranges, sseq_(3), strides(3), sizes(3), temp(3)

    sizes = SHAPE(mask)
    strides(1) = 1
    strides(2) = sizes(1)
    strides(3) = sizes(1) * sizes(2)

    IF (PRESENT(sseq)) THEN
      sseq_ = sseq
      ! temp is needed to account for a bug of NAG Fortran
      temp = sizes(sseq_) ; sizes = temp
      temp = strides(sseq_) ; strides = temp
    ELSE
      sseq_ = (/ 1, 2, 3 /)
    END IF

    num_ranges = build_ranges_mt_3d(mask, sizes, strides, sseq_, offsets)
    ALLOCATE(idx(num_ranges))
    num_ranges = build_ranges_mt_3d(mask, sizes, strides, sseq_, offsets, idx)
  END SUBROUTINE index_from_mask_mt_3d

END MODULE ppm_compact_mask_index
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
