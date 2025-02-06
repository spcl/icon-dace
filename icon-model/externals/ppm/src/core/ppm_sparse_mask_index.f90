!>
!! @file ppm_sparse_mask_index.f90
!! @brief sparsely populated mask array indices
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: mask index
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
!> ppm_sparse_mask_index
!! supposing you have a mask array and use it in something like a
!! WHERE, then it can be beneficial to only visit those indices of the
!! array expression with the condition evaluating to .TRUE., in case
!! this is relatively sparse or the ensuing conditional is expensive
!! in comparison to the WHERE/FORALL body
#include "fc_feature_defs.inc"
MODULE ppm_sparse_mask_index
  USE ppm_extents, ONLY: iinterval
  IMPLICIT NONE
  PRIVATE

  TYPE index_sparse_nd
    !> rank of original mask
    INTEGER :: rank
    !> number of extents describing the highest stride (top level) index
    INTEGER :: tl_ranges
    !> linear list of ranges
    TYPE(iinterval), ALLOCATABLE :: ranges(:)
  END TYPE index_sparse_nd

  INTERFACE index_from_mask
    MODULE PROCEDURE index_from_mask_1d
    MODULE PROCEDURE index_from_mask_2d
    MODULE PROCEDURE index_from_mask_3d
  END INTERFACE index_from_mask

  INTERFACE fold_mask
    MODULE PROCEDURE fold_mask_2d
    MODULE PROCEDURE fold_mask_3d
  END INTERFACE fold_mask

  PUBLIC :: index_sparse_nd, index_from_mask, build_ranges_1d, count_ranges_1d

CONTAINS

  !> construct ranges to represent true elements of mask
  PURE SUBROUTINE build_ranges_1d(mask, ofs, d_ranges, d_ranges_base, rcount)
    LOGICAL, INTENT(in) :: mask(:)
    INTEGER, INTENT(in) :: ofs
    TYPE(iinterval), INTENT(inout) :: d_ranges(:)
    INTEGER, OPTIONAL, INTENT(in) :: d_ranges_base
    INTEGER, OPTIONAL, INTENT(out) :: rcount

    INTEGER :: i, m, p, p_ofs

    m = SIZE(mask)
    p_ofs = 1
    IF (PRESENT(d_ranges_base)) p_ofs = d_ranges_base
    p = p_ofs
    i = 0
    DO WHILE(i < m)
      IF (mask(i + 1)) THEN
        d_ranges(p)%first = i + ofs
        DO WHILE (i + 1 < m)
          IF (.NOT. mask(i + 2)) EXIT
          i = i + 1
        END DO
        d_ranges(p)%last = i + ofs
        p = p + 1
      END IF
      i = i + 1
    END DO
    IF (PRESENT(rcount)) rcount = p - p_ofs
  END SUBROUTINE build_ranges_1d

  !> @result number of ranges required to represent true values in mask
  PURE FUNCTION count_ranges_1d(mask) RESULT(num_ranges)
    LOGICAL, INTENT(in) :: mask(:)

    INTEGER :: i, m, num_ranges

    m = SIZE(mask)
    num_ranges = 0
    i = 0
    DO WHILE(i < m)
      IF (mask(i + 1)) THEN
        DO WHILE (i + 1 < m)
          IF (.NOT. mask(i + 2)) EXIT
          i = i + 1
        END DO
        num_ranges = num_ranges + 1
      END IF
      i = i + 1
    END DO
  END FUNCTION count_ranges_1d

  PURE SUBROUTINE index_from_mask_1d(idx, mask, offsets)
    TYPE(index_sparse_nd), INTENT(out) :: idx
    LOGICAL, INTENT(in) :: mask(:)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(1)

    INTEGER :: offsets_dummy(1)

    idx%rank = 1

    offsets_dummy = 1
    IF (PRESENT(offsets)) offsets_dummy = offsets

    idx%tl_ranges = count_ranges_1d(mask)
    ALLOCATE(idx%ranges(idx%tl_ranges))

    CALL build_ranges_1d(mask, offsets_dummy(1), idx%ranges)
  END SUBROUTINE index_from_mask_1d

  ! this is equivalent to reduction = ANY(mask, dim=3-dim) but does not
  ! require a temporary array, and most versions of the PGI compiler
  ! do not support ANY with a dim argument anyway
  PURE SUBROUTINE fold_mask_2d(mask, dim, reduction)
    LOGICAL, INTENT(in) :: mask(:,:)
    INTEGER, INTENT(in) :: dim
    LOGICAL, INTENT(out) :: reduction(:)
    LOGICAL :: p
    INTEGER :: i, j, rsize, asize
    asize = SIZE(mask, 3-dim)
    rsize = SIZE(mask, dim)
    SELECT CASE(dim)
    CASE(2)
      DO i = 1, rsize
        p = .FALSE.
        DO j = 1, asize
          p = p .OR. mask(j, i)
        END DO
        reduction(i) = p
      END DO
    CASE(1)
      DO i = 1, rsize
        p = .FALSE.
        DO j = 1, asize
          p = p .OR. mask(i, j)
        END DO
        reduction(i) = p
      END DO
    END SELECT
  END SUBROUTINE fold_mask_2d

  ! this is equivalent to
  ! reduction = ANY(ANY(mask, dim=mod(dim+1,3)+1), dim=mod(dim+2,2)+1)
  ! but does not
  ! require a temporary array, and most versions of the PGI compiler
  ! do not support ANY with a dim argument anyway
  PURE SUBROUTINE fold_mask_3d(mask, dim, reduction)
    LOGICAL, INTENT(in) :: mask(:,:,:)
    INTEGER, INTENT(in) :: dim
    LOGICAL, INTENT(out) :: reduction(:)
    LOGICAL :: p
    INTEGER :: i, j, k, m, n, rsize
    rsize = SIZE(mask, dim)
    SELECT CASE(dim)
    CASE(3)
      m = SIZE(mask, 1)
      n = SIZE(mask, 2)
      DO k = 1, rsize
        p = .FALSE.
        DO j = 1, n
          DO i = 1, m
            p = p .OR. mask(i, j, k)
          END DO
        END DO
        reduction(k) = p
      END DO
    CASE(2)
      m = SIZE(mask, 1)
      n = SIZE(mask, 3)
      DO k = 1, rsize
        p = .FALSE.
        DO j = 1, n
          DO i = 1, m
            p = p .OR. mask(i, k, j)
          END DO
        END DO
        reduction(k) = p
      END DO
    CASE(1)
      m = SIZE(mask, 2)
      n = SIZE(mask, 3)
      DO k = 1, rsize
        p = .FALSE.
        DO j = 1, n
          DO i = 1, m
            p = p .OR. mask(k, i, j)
          END DO
        END DO
        reduction(k) = p
      END DO
    END SELECT
  END SUBROUTINE fold_mask_3d

  PURE SUBROUTINE build_ranges_2d(mask, subscript_sequence, offsets, &
       d_ranges, num_tl_ranges, d_ranges_base, rcount)
    LOGICAL, INTENT(in) :: mask(:, :)
    INTEGER, INTENT(in) :: subscript_sequence(2)
    INTEGER, INTENT(in) :: offsets(2)
    TYPE(iinterval), INTENT(inout) :: d_ranges(:)
    INTEGER, OPTIONAL, INTENT(out) :: num_tl_ranges, rcount
    INTEGER, OPTIONAL, INTENT(in) :: d_ranges_base

    LOGICAL :: mask_reduce(SIZE(mask, subscript_sequence(2)))
    INTEGER :: p, iblock, tl_ssdim, ll_ssdim, i, ntlr
    INTEGER :: p_ofs, q, subrange_count

    tl_ssdim = subscript_sequence(2)
    ll_ssdim = subscript_sequence(1)

    CALL fold_mask(mask, dim=tl_ssdim, reduction=mask_reduce)

    p_ofs = 1
    IF (PRESENT(d_ranges_base)) p_ofs = d_ranges_base
    p = p_ofs

    CALL build_ranges_1d(mask_reduce(:), offsets(tl_ssdim), &
         d_ranges=d_ranges(p_ofs::2), rcount=ntlr)
    IF (PRESENT(num_tl_ranges)) num_tl_ranges = ntlr
    p = p + 2 * ntlr
    SELECT CASE(tl_ssdim)
    CASE(2)
      DO iblock = 1, 2*ntlr, 2
        q = p_ofs + iblock - 1
        d_ranges(q + 1)%first = p
        p = p + d_ranges(q)%last - d_ranges(q)%first + 1
        d_ranges(q + 1)%last = p - 1
        DO i = d_ranges(q)%first, d_ranges(q)%last
          d_ranges(d_ranges(q + 1)%first &
               &   + i - d_ranges(q)%first)%first = p
          CALL build_ranges_1d(mask(:, i - offsets(2) + 1), &
               offsets(ll_ssdim), d_ranges(p:), rcount=subrange_count)
          p = p + subrange_count
          d_ranges(d_ranges(q + 1)%first &
               &   + i - d_ranges(q)%first)%last = p - 1
        END DO
      END DO
    CASE(1)
      DO iblock = 1, 2*ntlr, 2
        q = p_ofs + iblock - 1
        d_ranges(q + 1)%first = p
        p = p + d_ranges(q)%last - d_ranges(q)%first + 1
        d_ranges(q + 1)%last = p - 1
        DO i = d_ranges(q)%first, d_ranges(q)%last
          d_ranges(d_ranges(q + 1)%first &
               &   + i - d_ranges(q)%first)%first = p
          CALL build_ranges_1d(mask(i - offsets(1) + 1, :), &
               offsets(ll_ssdim), d_ranges(p:), rcount=subrange_count)
          p = p + subrange_count
          d_ranges(d_ranges(q + 1)%first &
               &   + i - d_ranges(q)%first)%last = p - 1
        END DO
      END DO
    END SELECT
    IF (PRESENT(rcount)) rcount = p - p_ofs
  END SUBROUTINE build_ranges_2d

  PURE FUNCTION count_ranges_2d(mask, subscript_sequence) RESULT(p)
    LOGICAL, INTENT(in) :: mask(:, :)
    INTEGER, INTENT(in) :: subscript_sequence(2)

    LOGICAL :: mask_reduce(SIZE(mask, subscript_sequence(2)))
    INTEGER :: p, tl_size, tl_ssdim, i, ntlr

    tl_ssdim = subscript_sequence(2)
    tl_size = SIZE(mask, tl_ssdim)

    CALL fold_mask(mask, dim=tl_ssdim, reduction=mask_reduce)

    p = 0

    ntlr = count_ranges_1d(mask_reduce)
    p = p + 2 * ntlr + COUNT(mask_reduce)
    SELECT CASE(tl_ssdim)
    CASE(2)
      DO i = 1, tl_size
        p = p + count_ranges_1d(mask(:, i))
      END DO
    CASE(1)
      DO i = 1, tl_size
        p = p + count_ranges_1d(mask(i, :))
      END DO
    END SELECT
  END FUNCTION count_ranges_2d

  PURE SUBROUTINE index_from_mask_2d(idx, mask, offsets, sseq)
    TYPE(index_sparse_nd), INTENT(out) :: idx
    LOGICAL, INTENT(in) :: mask(:, :)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(2), sseq(2)

    INTEGER :: offsets_dummy(2), subscript_sequence(2)
    INTEGER :: num_ranges

    offsets_dummy = 1
    IF (PRESENT(offsets)) offsets_dummy = offsets

    subscript_sequence(1) = 1
    subscript_sequence(2) = 2
    IF (PRESENT(sseq)) subscript_sequence = sseq

    idx%rank = 2
    num_ranges = count_ranges_2d(mask, subscript_sequence)
    ALLOCATE(idx%ranges(num_ranges))
    CALL build_ranges_2d(mask, subscript_sequence, offsets_dummy, &
         idx%ranges, idx%tl_ranges)
  END SUBROUTINE index_from_mask_2d

  PURE SUBROUTINE build_ranges_3d(mask, subscript_sequence, offsets, &
       d_ranges, num_tl_ranges, rcount)
    LOGICAL, INTENT(in) :: mask(:, :, :)
    INTEGER, INTENT(in) :: subscript_sequence(3)
    INTEGER, INTENT(in) :: offsets(3)
    TYPE(iinterval), INTENT(inout) :: d_ranges(:)
    INTEGER, OPTIONAL, INTENT(out) :: num_tl_ranges, rcount

    LOGICAL :: mask_reduce(SIZE(mask, subscript_sequence(3)))
    INTEGER :: ll_offsets(2), ll_subscripts(2)
    INTEGER :: p, iblock, tl_ssdim, i, j, ntlr, nllr, r, subrange_count

    tl_ssdim = subscript_sequence(3)
    j = 1
    DO i = 1, 3
      IF (i /= tl_ssdim) THEN
        ll_offsets(j) = offsets(i)
        j = j + 1
      END IF
    END DO
    DO i = 1, 2
      ll_subscripts(i) = subscript_sequence(i) &
           - MERGE(1, 0, subscript_sequence(i) > tl_ssdim)
    END DO

    CALL fold_mask(mask, dim=tl_ssdim, reduction=mask_reduce)

    CALL build_ranges_1d(mask_reduce(:), offsets(tl_ssdim), &
         d_ranges(1::2), rcount=ntlr)
    IF (PRESENT(num_tl_ranges)) num_tl_ranges = ntlr
    p = 2 * ntlr + 1
    SELECT CASE(tl_ssdim)
    CASE(3)
      DO iblock = 1, 2*ntlr, 2
        d_ranges(iblock + 1)%first = p
        p = p + d_ranges(iblock)%last - d_ranges(iblock)%first + 1
        d_ranges(iblock + 1)%last = p - 1
        DO i = d_ranges(iblock)%first, d_ranges(iblock)%last
          r = d_ranges(iblock + 1)%first + i - d_ranges(iblock)%first
          d_ranges(r)%first = p
          CALL build_ranges_2d(mask(:, :, i - offsets(3) + 1), &
               ll_subscripts, ll_offsets, d_ranges, num_tl_ranges=nllr, &
               d_ranges_base=p, rcount=subrange_count)
          p = p + subrange_count
          d_ranges(r)%last = d_ranges(r)%first + 2 * (nllr - 1)
        END DO
      END DO
    CASE(2)
      DO iblock = 1, 2*ntlr, 2
        d_ranges(iblock + 1)%first = p
        p = p + d_ranges(iblock)%last - d_ranges(iblock)%first + 1
        d_ranges(iblock + 1)%last = p - 1
        DO i = d_ranges(iblock)%first, d_ranges(iblock)%last
          r = d_ranges(iblock + 1)%first + i - d_ranges(iblock)%first
          d_ranges(r)%first = p
          CALL build_ranges_2d(mask(:, i - offsets(2) + 1, :), &
               ll_subscripts, &
               ll_offsets, d_ranges, num_tl_ranges=nllr, d_ranges_base=p, &
               rcount=subrange_count)
          p = p + subrange_count
          d_ranges(r)%last = d_ranges(r)%first + 2 * (nllr - 1)
        END DO
      END DO
    CASE(1)
      DO iblock = 1, 2*ntlr, 2
        d_ranges(iblock + 1)%first = p
        p = p + d_ranges(iblock)%last - d_ranges(iblock)%first + 1
        d_ranges(iblock + 1)%last = p - 1
        DO i = d_ranges(iblock)%first, d_ranges(iblock)%last
          r = d_ranges(iblock + 1)%first + i - d_ranges(iblock)%first
          d_ranges(r)%first = p
          CALL build_ranges_2d(mask(i - offsets(1) + 1, :, :), &
               ll_subscripts, ll_offsets, d_ranges, num_tl_ranges=nllr, &
               d_ranges_base=p, rcount=subrange_count)
          p = p + subrange_count
          d_ranges(r)%last = d_ranges(r)%first + 2 * (nllr - 1)
        END DO
      END DO
    END SELECT
    IF (PRESENT(rcount)) rcount = p - 1
  END SUBROUTINE build_ranges_3d

  PURE FUNCTION count_ranges_3d(mask, subscript_sequence) RESULT(num_ranges)
    LOGICAL, INTENT(in) :: mask(:, :, :)
    INTEGER, INTENT(in) :: subscript_sequence(3)

    LOGICAL :: mask_reduce(SIZE(mask, subscript_sequence(3)))
    INTEGER :: ll_subscripts(2), num_ranges, tl_size, tl_ssdim, i, ntlr

    tl_ssdim = subscript_sequence(3)
    DO i = 1, 2
      ll_subscripts(i) = subscript_sequence(i) &
           - MERGE(1, 0, subscript_sequence(i) > tl_ssdim)
    END DO
    tl_size = SIZE(mask, tl_ssdim)

    CALL fold_mask(mask, tl_ssdim, mask_reduce)

    ntlr = count_ranges_1d(mask_reduce)
    num_ranges = 2 * ntlr + COUNT(mask_reduce)
    SELECT CASE(tl_ssdim)
    CASE(3)
      DO i = 1, tl_size
        num_ranges = num_ranges + count_ranges_2d(mask(:, :, i), ll_subscripts)
      END DO
    CASE(2)
      DO i = 1, tl_size
        num_ranges = num_ranges + count_ranges_2d(mask(:, i, :), ll_subscripts)
      END DO
    CASE(1)
      DO i = 1, tl_size
        num_ranges = num_ranges + count_ranges_2d(mask(i, :, :), ll_subscripts)
      END DO
    END SELECT
  END FUNCTION count_ranges_3d

  PURE SUBROUTINE index_from_mask_3d(idx, mask, offsets, sseq)
    TYPE(index_sparse_nd), INTENT(out) :: idx
    LOGICAL, INTENT(in) :: mask(:, :, :)
    INTEGER, OPTIONAL, INTENT(in) :: offsets(3), sseq(3)

    INTEGER :: offsets_dummy(3), subscript_sequence(3)
    INTEGER :: num_ranges

    offsets_dummy = 1
    IF (PRESENT(offsets)) offsets_dummy = offsets

    subscript_sequence = (/ 1, 2, 3 /)
    IF (PRESENT(sseq)) subscript_sequence = sseq

    idx%rank = 3
    num_ranges = count_ranges_3d(mask, subscript_sequence)
    ALLOCATE(idx%ranges(num_ranges))
    CALL build_ranges_3d(mask, subscript_sequence, offsets_dummy, &
         idx%ranges, idx%tl_ranges)
  END SUBROUTINE index_from_mask_3d

END MODULE ppm_sparse_mask_index
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
