!> @file ppm_set_partition_base.f90
!! @brief basic routines and data structures for handling partitions
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
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
!> basic routines and data structures for handling partitions

MODULE ppm_set_partition_base
  USE ppm_std_type_kinds, ONLY: dp, i4, i8




  USE ppm_base, ONLY: assertion, abort_ppm
  USE ppm_extents, ONLY: extent, extent_start, extent_end, empty_extent, &
       iinterval, ASSIGNMENT(=)
  USE ppm_f90_io_lun, ONLY: next_free_unit



  IMPLICIT NONE
  PRIVATE



  !> succinct representation of partitioning, where
  !! <tt>elements(start(p):start(p+1)-1)</tt> contains
  !! the indices \c i of partition \c p
  TYPE partition_vec
    INTEGER, ALLOCATABLE :: start(:), elements(:)
  END TYPE partition_vec

  !> easily changeable partition descriptor, to be used as array of
  !! size n for description of partitioning into n partitions
  TYPE set_i4
    INTEGER(i4), ALLOCATABLE :: elem(:)
  END TYPE set_i4

  !> denotes partitioning by tabulating for each element i partition p_i
  !! as p_i == assigned(i), also part_range must denote the legal range for p_i
  TYPE partition_assignment
    TYPE(extent) :: part_range
    INTEGER(i4), ALLOCATABLE :: assigned(:)
  END TYPE partition_assignment

  !> describe range decomposed into non-overlapping contiguous ranges
  TYPE block_decomposition
    TYPE(extent), ALLOCATABLE :: partition(:)
  END TYPE block_decomposition


  INTERFACE OPERATOR(==)
    MODULE PROCEDURE eq_set_i4_partition_vec
    MODULE PROCEDURE eq_partition_vec_set_i4
    MODULE PROCEDURE eq_pa_set_i4
    MODULE PROCEDURE eq_set_i4_pa
    MODULE PROCEDURE eq_pv_pa
    MODULE PROCEDURE eq_pa_pv
    MODULE PROCEDURE eq_pa_pa
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE neq_set_i4_partition_vec
    MODULE PROCEDURE neq_partition_vec_set_i4
    MODULE PROCEDURE neq_pa_set_i4
    MODULE PROCEDURE neq_set_i4_pa
    MODULE PROCEDURE neq_pv_pa
    MODULE PROCEDURE neq_pa_pv
    MODULE PROCEDURE neq_pa_pa
  END INTERFACE OPERATOR(/=)

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE assign_pa2pvec
    MODULE PROCEDURE assign_pvec2pa
    MODULE PROCEDURE assign_pa2set_i4
  END INTERFACE ASSIGNMENT(=)

  INTERFACE balance_of_max
    MODULE PROCEDURE balance_of_max_partitioning_i4
    MODULE PROCEDURE balance_of_max_partitioning_dp
    MODULE PROCEDURE balance_of_max_weight_sums_i8
    MODULE PROCEDURE balance_of_max_weight_sums_dp
  END INTERFACE balance_of_max








  INTERFACE partition_weight_sums
    MODULE PROCEDURE partition_weight_sums_i8i4
    MODULE PROCEDURE partition_weight_sums_dpdp
  END INTERFACE partition_weight_sums

  INTERFACE part_size
    MODULE PROCEDURE part_size_pv
  END INTERFACE part_size

  PUBLIC :: partition_vec, set_i4, partition_assignment, block_decomposition
  PUBLIC :: balance_of_max, OPERATOR(==), OPERATOR(/=), ASSIGNMENT(=)
  PUBLIC :: part_size



  PUBLIC :: partition_weight_sums
  ! these cannot be part of an generic assignment interface
  PUBLIC :: assign_set_i4_2_pv, assign_pv_2_set_i4

  INTERFACE
    SUBROUTINE ppm_read_int_array(filename, destcount, dest, icnt)
      USE ppm_std_type_kinds, ONLY: i4
      CHARACTER(len=*), INTENT(in) :: filename
      INTEGER, INTENT(in) :: destcount
      INTEGER(i4), INTENT(out) :: dest(*)
      INTEGER, INTENT(out) :: icnt
    END SUBROUTINE ppm_read_int_array
  END INTERFACE

  INTERFACE read_partitioning
    MODULE PROCEDURE read_partitioning_pa
    MODULE PROCEDURE read_partitioning_pv
    MODULE PROCEDURE read_partitioning_set_i4
  END INTERFACE read_partitioning

  INTERFACE write_partition
    MODULE PROCEDURE write_partition_pa
  END INTERFACE write_partition

  PUBLIC :: read_partitioning
  PUBLIC :: write_partition


  CHARACTER(len=*), PARAMETER :: filename = 'ppm_set_partition_base.f90'
CONTAINS
  FUNCTION eq_set_i4_partition_vec(a, b) RESULT(p)
    TYPE(partition_vec), INTENT(IN) :: a
    TYPE(set_i4), INTENT(in) :: b(LBOUND(a%start,1):)
    LOGICAL :: p

    INTEGER :: nparts, i, p_lb, p_ub

    CALL assertion(ALLOCATED(a%start) .AND. ALLOCATED(a%elements), filename, &
         177, "partition_vec must be initialized")
    nparts = SIZE(b)
    p = nparts == SIZE(a%start) - 1

    p_lb = LBOUND(a%start, 1)
    p_ub = UBOUND(a%start, 1) - 1
    DO i = p_lb, p_ub
      p = p .AND. SIZE(b(i)%elem) == a%start(i+1) - a%start(i)
      IF (.NOT. p) EXIT
      p = p .AND. ALL(a%elements(a%start(i):a%start(i+1)-1) == b(i)%elem)
    END DO
  END FUNCTION eq_set_i4_partition_vec

  FUNCTION eq_partition_vec_set_i4(a, b) RESULT(p)
    TYPE(set_i4), INTENT(in) :: a(:)
    TYPE(partition_vec), INTENT(IN) :: b
    LOGICAL :: p
    p = b == a
  END FUNCTION eq_partition_vec_set_i4

  FUNCTION eq_pa_set_i4(a, b) RESULT(p)
    TYPE(partition_assignment), INTENT(in) :: a
    TYPE(set_i4), INTENT(in) :: b(:)
    LOGICAL :: p

    INTEGER :: part_id, e, i, lb_e, ub_e, lb_p, ub_p, lb_s, ub_s, size_s
    LOGICAL, ALLOCATABLE :: seen(:)

    lb_s = extent_start(a%part_range)
    ub_s = extent_end(a%part_range)
    size_s = SIZE(b)
    p = size_s == ub_s - lb_s + 1
    IF (p) THEN
      lb_e = LBOUND(a%assigned, 1)
      ub_e = UBOUND(a%assigned, 1)
      ALLOCATE(seen(lb_e:ub_e))
      seen = .FALSE.

      compare_loop: DO part_id = 1, size_s
        lb_p = LBOUND(b(part_id)%elem, 1)
        ub_p = UBOUND(b(part_id)%elem, 1)
        DO i = lb_p, ub_p
          e = b(part_id)%elem(i)
          IF (seen(e) .OR. a%assigned(e) /= part_id) THEN
            p = .FALSE.
            EXIT compare_loop
          END IF
          seen(e) = .TRUE.
        END DO
      END DO compare_loop
      p = p .AND. ALL(seen)
      DEALLOCATE(seen)
    END IF
  END FUNCTION eq_pa_set_i4

  FUNCTION eq_set_i4_pa(a, b) RESULT(p)
    TYPE(set_i4), INTENT(in) :: a(:)
    TYPE(partition_assignment), INTENT(in) :: b
    LOGICAL :: p
    p = b == a
  END FUNCTION eq_set_i4_pa

  FUNCTION eq_pv_pa(a, b) RESULT(p)
    TYPE(partition_vec), INTENT(in) :: a
    TYPE(partition_assignment), INTENT(in) :: b
    LOGICAL :: p

    INTEGER :: part_id, e, i, lb_e, ub_e, lb_p, ub_p, lb_s, ub_s, size_s
    LOGICAL, ALLOCATABLE :: seen(:)

    lb_s = extent_start(b%part_range)
    ub_s = extent_end(b%part_range)
    size_s = SIZE(a%start) - 1
    p = lb_s == LBOUND(a%start, 1) .AND. size_s == ub_s - lb_s + 1
    IF (p) THEN
      lb_e = LBOUND(b%assigned, 1)
      ub_e = UBOUND(b%assigned, 1)
      ALLOCATE(seen(lb_e:ub_e))
      seen = .FALSE.

      compare_loop: DO part_id = lb_s, ub_s
        lb_p = a%start(part_id)
        ub_p = a%start(part_id+1)-1
        DO i = lb_p, ub_p
          e = a%elements(i)
          IF (seen(e) .OR. b%assigned(e) /= part_id) THEN
            p = .FALSE.
            EXIT compare_loop
          END IF
          seen(e) = .TRUE.
        END DO
      END DO compare_loop
      p = p .AND. ALL(seen)
      DEALLOCATE(seen)
    END IF
  END FUNCTION eq_pv_pa

  FUNCTION eq_pa_pv(a, b) RESULT(p)
    TYPE(partition_assignment), INTENT(in) :: a
    TYPE(partition_vec), INTENT(in) :: b
    LOGICAL :: p
    p = b == a
  END FUNCTION eq_pa_pv

  FUNCTION eq_pa_pa(a, b) RESULT(p)
    TYPE(partition_assignment), INTENT(in) :: a, b
    LOGICAL :: p
    IF (ALLOCATED(a%assigned) .AND. ALLOCATED(b%assigned)) THEN
      p = SIZE(a%assigned) == SIZE(b%assigned) &
           .AND. LBOUND(a%assigned, 1) == LBOUND(b%assigned, 1)
      IF (.NOT. p) RETURN
      p = ALL(a%assigned == b%assigned)
    ELSE IF (.NOT. ALLOCATED(a%assigned) .AND. .NOT. ALLOCATED(b%assigned)) THEN
      p = .TRUE.
    ELSE
      p = .FALSE.
    END IF
  END FUNCTION eq_pa_pa

  FUNCTION neq_set_i4_partition_vec(a, b) RESULT(p)
    TYPE(partition_vec), INTENT(IN) :: a
    TYPE(set_i4), INTENT(in) :: b(LBOUND(a%start,1):)
    LOGICAL :: p
    p = .NOT. a == b
  END FUNCTION neq_set_i4_partition_vec

  FUNCTION neq_partition_vec_set_i4(a, b) RESULT(p)
    TYPE(set_i4), INTENT(in) :: a(:)
    TYPE(partition_vec), INTENT(IN) :: b
    LOGICAL :: p
    p = .NOT. a == b
  END FUNCTION neq_partition_vec_set_i4

  FUNCTION neq_pa_set_i4(a, b) RESULT(p)
    TYPE(partition_assignment), INTENT(in) :: a
    TYPE(set_i4), INTENT(in) :: b(:)
    LOGICAL :: p
    p = .NOT. a == b
  END FUNCTION neq_pa_set_i4

  FUNCTION neq_set_i4_pa(a, b) RESULT(p)
    TYPE(set_i4), INTENT(in) :: a(:)
    TYPE(partition_assignment), INTENT(in) :: b
    LOGICAL :: p
    p = .NOT. a == b
  END FUNCTION neq_set_i4_pa

  FUNCTION neq_pv_pa(a, b) RESULT(p)
    TYPE(partition_vec), INTENT(in) :: a
    TYPE(partition_assignment), INTENT(in) :: b
    LOGICAL :: p
    p = .NOT. a == b
  END FUNCTION neq_pv_pa

  FUNCTION neq_pa_pv(a, b) RESULT(p)
    TYPE(partition_assignment), INTENT(in) :: a
    TYPE(partition_vec), INTENT(in) :: b
    LOGICAL :: p
    p = .NOT. a == b
  END FUNCTION neq_pa_pv

  FUNCTION neq_pa_pa(a, b) RESULT(p)
    TYPE(partition_assignment), INTENT(in) :: a, b
    LOGICAL :: p
    p = .NOT. a == b
  END FUNCTION neq_pa_pa

  FUNCTION balance_of_max_partitioning_i4(partitioning, weight) &
       RESULT(imbalance)
    TYPE(partition_vec), INTENT(in) :: partitioning
    INTEGER(i4), INTENT(in) :: weight(:)
    REAL :: imbalance

    INTEGER(i8) :: mean_weight_sum, max_weight_sum
    INTEGER(i8), ALLOCATABLE :: weight_sums(:)
    INTEGER :: num_parts

    num_parts = SIZE(partitioning%start) - 1
    ALLOCATE(weight_sums(num_parts))
    CALL partition_weight_sums(weight_sums, partitioning, weight, &
         mean_weight_sum)
    max_weight_sum = MAXVAL(weight_sums)

    DEALLOCATE(weight_sums)
    imbalance = REAL(max_weight_sum)/REAL(mean_weight_sum)

  END FUNCTION balance_of_max_partitioning_i4

  FUNCTION balance_of_max_partitioning_dp(partitioning, weight) &
       RESULT(imbalance)
    TYPE(partition_vec), INTENT(in) :: partitioning
    REAL(dp), INTENT(in) :: weight(:)
    REAL :: imbalance

    REAL(dp) :: mean_weight_sum, max_weight_sum
    REAL(dp), ALLOCATABLE :: weight_sums(:)
    INTEGER :: num_parts

    num_parts = SIZE(partitioning%start) - 1
    ALLOCATE(weight_sums(num_parts))
    CALL partition_weight_sums(weight_sums, partitioning, weight, &
         mean_weight_sum)
    max_weight_sum = MAXVAL(weight_sums)

    DEALLOCATE(weight_sums)
    imbalance = REAL(max_weight_sum/mean_weight_sum)

  END FUNCTION balance_of_max_partitioning_dp

  FUNCTION balance_of_max_weight_sums_i8(weight_sums) RESULT(imbalance)
    INTEGER(i8), INTENT(in) :: weight_sums(:)
    REAL :: imbalance

    INTEGER(i8) :: mean_weight_sum, max_weight_sum

    max_weight_sum = MAXVAL(weight_sums)
    mean_weight_sum = SUM(weight_sums)/INT(SIZE(weight_sums), i8)

    imbalance = REAL(max_weight_sum)/REAL(mean_weight_sum)

  END FUNCTION balance_of_max_weight_sums_i8

  FUNCTION balance_of_max_weight_sums_dp(weight_sums) RESULT(imbalance)
    REAL(dp), INTENT(in) :: weight_sums(:)
    REAL :: imbalance

    REAL(dp) :: mean_weight_sum, max_weight_sum

    max_weight_sum = MAXVAL(weight_sums)
    mean_weight_sum = SUM(weight_sums)/REAL(SIZE(weight_sums), dp)

    imbalance = REAL(max_weight_sum/mean_weight_sum)

  END FUNCTION balance_of_max_weight_sums_dp


  SUBROUTINE partition_weight_sums_i8i4(weight_sums, partitioning, weight, &
       mean_weight_sum)
    TYPE(partition_vec), INTENT(in) :: partitioning
    INTEGER(i4), INTENT(in) :: weight(:)
    INTEGER(i8), INTENT(out) :: weight_sums(:)
    INTEGER(i8), OPTIONAL, INTENT(out) :: mean_weight_sum

    INTEGER(i8) :: accum
    INTEGER :: i, j, lb, ub, num_parts

    num_parts = SIZE(partitioning%start) - 1
    CALL assertion(num_parts == SIZE(weight_sums), filename, 424, &
         'number of partitions must match number of weight sums')
    CALL assertion(SIZE(partitioning%elements) == SIZE(weight), filename, &
         427, 'number of elements must match number of weights')

    DO j = 1, num_parts
      lb = partitioning%start(j)
      ub = partitioning%start(j+1) - 1
      accum = 0_i8
      DO i = lb, ub
        accum = accum + INT(weight(partitioning%elements(i)), i8)
      END DO
      weight_sums(j) = accum
    END DO
    IF (PRESENT(mean_weight_sum)) mean_weight_sum &
         = SUM(weight_sums)/INT(num_parts, i8)
  END SUBROUTINE partition_weight_sums_i8i4


  SUBROUTINE partition_weight_sums_dpdp(weight_sums, partitioning, weight, &
       mean_weight_sum)
    TYPE(partition_vec), INTENT(in) :: partitioning
    REAL(dp), INTENT(in) :: weight(:)
    REAL(dp), INTENT(out) :: weight_sums(:)
    REAL(dp), OPTIONAL, INTENT(out) :: mean_weight_sum

    REAL(dp) :: accum
    INTEGER :: i, j, lb, ub, num_parts

    num_parts = SIZE(partitioning%start) - 1
    CALL assertion(num_parts == SIZE(weight_sums), filename, 454, &
         'number of partitions must match number of weight sums')
    CALL assertion(SIZE(partitioning%elements) == SIZE(weight), filename, &
         457, 'number of elements must match number of weights')

    DO j = 1, num_parts
      lb = partitioning%start(j)
      ub = partitioning%start(j+1) - 1
      accum = 0.0_dp
      DO i = lb, ub
        accum = accum + weight(partitioning%elements(i))
      END DO
      weight_sums(j) = accum
    END DO
    IF (PRESENT(mean_weight_sum)) mean_weight_sum &
         = SUM(weight_sums)/REAL(num_parts, dp)
  END SUBROUTINE partition_weight_sums_dpdp


  ELEMENTAL SUBROUTINE assign_pa2pvec(pvec, pa)
    TYPE(partition_vec), INTENT(out) :: pvec
    TYPE(partition_assignment), INTENT(in) :: pa

    INTEGER :: lb_p, ub_p, partid, i, j
    INTEGER, ALLOCATABLE :: part_size(:)
    lb_p = extent_start(pa%part_range)
    ub_p = extent_end(pa%part_range)
    ALLOCATE(pvec%start(lb_p:ub_p+1), &
         pvec%elements(LBOUND(pa%assigned, 1):UBOUND(pa%assigned, 1)))
    ALLOCATE(part_size(lb_p:ub_p))
    part_size = 0
    DO i = LBOUND(pa%assigned, 1), UBOUND(pa%assigned, 1)
      partid = pa%assigned(i)
      part_size(partid) = part_size(partid) + 1
    END DO
    pvec%start(lb_p) = 1
    DO partid = lb_p, ub_p
      pvec%start(partid+1) = pvec%start(partid) + part_size(partid)
    END DO
    part_size = 0
    DO i = LBOUND(pa%assigned, 1), UBOUND(pa%assigned, 1)
      partid = pa%assigned(i)
      j = pvec%start(partid) + part_size(partid)
      pvec%elements(j) = i
      part_size(partid) = part_size(partid) + 1
    END DO
  END SUBROUTINE assign_pa2pvec

  ELEMENTAL SUBROUTINE assign_pvec2pa(pa, pvec)
    TYPE(partition_assignment), INTENT(out) :: pa
    TYPE(partition_vec), INTENT(in) :: pvec

    INTEGER :: partid, i, j, lb_s, ub_s

    ALLOCATE(pa%assigned(LBOUND(pvec%elements, 1):UBOUND(pvec%elements, 1)))
    lb_s = LBOUND(pvec%start, 1)
    ub_s = UBOUND(pvec%start, 1) - 1
    pa%part_range = iinterval(lb_s, ub_s)

    DO partid = lb_s, ub_s
      DO i = pvec%start(partid), pvec%start(partid + 1) - 1
        j = pvec%elements(i)
        pa%assigned(j) = partid
      END DO
    END DO

  END SUBROUTINE assign_pvec2pa

  SUBROUTINE assign_pa2set_i4(sets, pa)
    TYPE(set_i4), ALLOCATABLE, INTENT(inout) :: sets(:)
    TYPE(partition_assignment), INTENT(in) :: pa

    INTEGER :: num_parts, part_id, i, lb_s, ub_s, lb_e, ub_e
    INTEGER, ALLOCATABLE :: part_size(:)

    num_parts = SIZE(sets)

    lb_s = extent_start(pa%part_range)
    ub_s = extent_end(pa%part_range)
    IF (LBOUND(sets, 1) /= lb_s .OR. num_parts /= ub_s - lb_s + 1) THEN
      DEALLOCATE(sets)
      ALLOCATE(sets(lb_s:ub_s))
      num_parts = ub_s - lb_s + 1
    END IF

    lb_e = LBOUND(pa%assigned, 1)
    ub_e = UBOUND(pa%assigned, 1)
    ALLOCATE(part_size(lb_s:ub_s))
    part_size = 0

    DO i = lb_e, ub_e
      part_id = pa%assigned(i)
      part_size(part_id) = part_size(part_id) + 1
    END DO

    DO part_id = lb_s, ub_s
      IF (ALLOCATED(sets(part_id)%elem)) DEALLOCATE(sets(part_id)%elem)
      ALLOCATE(sets(part_id)%elem(part_size(part_id)))
    END DO

    part_size = 0

    DO i = lb_e, ub_e
      part_id = pa%assigned(i)
      part_size(part_id) = part_size(part_id) + 1
      sets(part_id)%elem(part_size(part_id)) = i
    END DO
    DEALLOCATE(part_size)
  END SUBROUTINE assign_pa2set_i4

  SUBROUTINE assign_set_i4_2_pv(pv, sets)
    TYPE(partition_vec), INTENT(out) :: pv
    TYPE(set_i4), ALLOCATABLE, INTENT(in) :: sets(:)

    INTEGER :: lb_s, ub_s, sum_set_sizes, part_id

    lb_s = LBOUND(sets, 1)
    ub_s = UBOUND(sets, 1)

    IF (ALLOCATED(pv%start)) DEALLOCATE(pv%start)
    ALLOCATE(pv%start(lb_s:ub_s+1))
    sum_set_sizes = 0
    DO part_id = lb_s, ub_s
      sum_set_sizes = sum_set_sizes + SIZE(sets(part_id)%elem)
    END DO
    IF (ALLOCATED(pv%elements)) DEALLOCATE(pv%elements)
    ALLOCATE(pv%elements(sum_set_sizes))

    pv%start(lb_s) = 1
    DO part_id = lb_s, ub_s
      pv%start(part_id + 1) = pv%start(part_id) + SIZE(sets(part_id)%elem)
      pv%elements(pv%start(part_id):pv%start(part_id + 1)-1) &
           = sets(part_id)%elem
    END DO
  END SUBROUTINE assign_set_i4_2_pv

  SUBROUTINE assign_pv_2_set_i4(sets, pv)
    TYPE(set_i4), ALLOCATABLE, INTENT(out) :: sets(:)
    TYPE(partition_vec), INTENT(in) :: pv

    INTEGER :: part_id, lb_s, ub_s, size_p

    lb_s = LBOUND(pv%start, 1)
    ub_s = UBOUND(pv%start, 1) - 1
    IF (ALLOCATED(sets)) THEN
      IF (LBOUND(sets, 1) /= lb_s .OR. UBOUND(sets, 1) /= ub_s) THEN
        DEALLOCATE(sets)
        ALLOCATE(sets(lb_s:ub_s))
      END IF
    ELSE
      ALLOCATE(sets(lb_s:ub_s))
    END IF
    DO part_id = lb_s, ub_s
      size_p = pv%start(part_id+1) - pv%start(part_id)
      IF (ALLOCATED(sets(part_id)%elem)) THEN
        IF (SIZE(sets(part_id)%elem) /= size_p) THEN
          DEALLOCATE(sets(part_id)%elem)
          ALLOCATE(sets(part_id)%elem(size_p))
        END IF
      ELSE
        ALLOCATE(sets(part_id)%elem(size_p))
      END IF
      sets(part_id)%elem = pv%elements(pv%start(part_id):pv%start(part_id+1)-1)
    END DO
  END SUBROUTINE assign_pv_2_set_i4

  SUBROUTINE read_partitioning_pa(filename, partitioning, ierror)
    CHARACTER(len=*), INTENT(in) :: filename
    TYPE(partition_assignment), INTENT(out) :: partitioning
    INTEGER, INTENT(out) :: ierror

    INTEGER :: node_count, lb_p, ub_p
    INTEGER(i4) :: temp(1)

    CALL ppm_read_int_array(filename, 0, temp, ierror)
    IF (ierror > 0) THEN
      node_count = ierror
      ALLOCATE(partitioning%assigned(node_count))
      CALL ppm_read_int_array(filename, node_count, &
           partitioning%assigned, ierror)
      lb_p = MINVAL(partitioning%assigned)
      ub_p = MAXVAL(partitioning%assigned)
      partitioning%part_range = iinterval(lb_p, ub_p)
    ELSE IF (ierror == 0) THEN
      partitioning%part_range = empty_extent
    END IF
  END SUBROUTINE read_partitioning_pa

  SUBROUTINE read_partitioning_pv(filename, partitioning, ierror)
    CHARACTER(len=*), INTENT(in) :: filename
    TYPE(partition_vec), INTENT(out) :: partitioning
    INTEGER, INTENT(out) :: ierror

    TYPE(partition_assignment) :: element_assignment

    CALL read_partitioning(filename, element_assignment, ierror)
    IF (ierror >= 0) THEN
      partitioning = element_assignment
      IF (ALLOCATED(element_assignment%assigned)) &
           DEALLOCATE(element_assignment%assigned)
    END IF
  END SUBROUTINE read_partitioning_pv

  SUBROUTINE read_partitioning_set_i4(filename, partitioning, ierror)
    CHARACTER(len=*), INTENT(in) :: filename
    TYPE(set_i4), ALLOCATABLE, INTENT(out) :: partitioning(:)
    INTEGER, INTENT(out) :: ierror

    TYPE(partition_assignment) :: element_assignment

    CALL read_partitioning(filename, element_assignment, ierror)
    IF (ierror >= 0) THEN
      partitioning = element_assignment
      IF (ALLOCATED(element_assignment%assigned)) &
           DEALLOCATE(element_assignment%assigned)
    END IF
  END SUBROUTINE read_partitioning_set_i4

  SUBROUTINE write_partition_pa(filename, partition, ierror)
    CHARACTER(len=*), INTENT(in) :: filename
    TYPE(partition_assignment), INTENT(in) :: partition
    INTEGER, OPTIONAL, INTENT(out) :: ierror

    INTEGER :: i, unit, ierr
    LOGICAL :: p

    ierr = -1
    CALL assertion(ALLOCATED(partition%assigned), filename, 742, &
         'invalid partition argument')
    DO i = 1, 1
      unit = next_free_unit()
      OPEN(unit=unit, file=filename, action='write', iostat=ierr)
      IF (ierr /= 0) EXIT
      WRITE (unit=unit, fmt='(i0)', iostat=ierr) partition%assigned
      IF (ierr /= 0) EXIT
      CLOSE(unit, iostat=ierr)
    END DO
    IF (PRESENT(ierror)) THEN
      ierror = ierr
    ELSE IF (ierr /= 0) THEN
      INQUIRE(unit=unit, opened=p)
      IF (p) CLOSE(unit)
      CALL abort_ppm('file error when writing partition', filename, 757)
    END IF
  END SUBROUTINE write_partition_pa

  FUNCTION part_size_pv(pv, part) RESULT(psize)
    INTEGER :: psize
    TYPE(partition_vec), INTENT(in) :: pv
    INTEGER, INTENT(in) :: part
    psize = pv%start(part+1) - pv%start(part)
  END FUNCTION part_size_pv

END MODULE ppm_set_partition_base
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
