!> @file ppm_set_repartition.f90
!! @brief routines for repartitioning
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords: repartition
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
MODULE ppm_set_repartition
  USE ppm_base, ONLY: assertion
  USE ppm_std_type_kinds, ONLY: i4, i8, dp, sp, sizeof_integer
  USE ppm_compare, ONLY: cmp_i4_indirect_i8, cmp_i4_indirect_dp, &
       rcmp_i4_indirect_dp, rcmp_i4_indirect_i8
  USE ppm_set_partition_base, ONLY: partition_vec, partition_weight_sums, &
       part_size
  USE ppm_sort, ONLY: qsort_r
  USE ppm_heap, ONLY: heap_elem_increase_sort, heapify, heap_leaf_minimize
#ifdef USE_MPI
  USE ppm_base, ONLY: abort_ppm, ppm_default_comm
  USE ppm_std_type_kinds_mp, ONLY: mp_dp, mp_i8, mp_i4, mp_dp_extent
  USE ppm_compare, ONLY: rcmp_i4, rcmp_dp
  USE ppm_sort, ONLY: insertion_sort_once
#endif
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif

  PUBLIC :: initialize_set_repartition, finalize_set_repartition

  INTERFACE repartition_swap
    MODULE PROCEDURE repartition_swap_i4
    MODULE PROCEDURE repartition_swap_dp
  END INTERFACE repartition_swap

  PUBLIC :: repartition_swap

  TYPE swap_pair
    INTEGER :: a, b
  END TYPE swap_pair

#ifdef USE_MPI

  TYPE part_desc_i4
    SEQUENCE
    INTEGER(i8) :: weight_sum
    INTEGER(i4) :: part_size, part_rank
  END TYPE part_desc_i4

  TYPE part_desc_dp
    SEQUENCE
    REAL(dp) :: weight_sum
    INTEGER(i4) :: part_size, part_rank
  END TYPE part_desc_dp

  TYPE elem_desc_i4
    SEQUENCE
    INTEGER(i4) :: elem, weight
  END TYPE elem_desc_i4

  TYPE elem_desc_dp
    SEQUENCE
    INTEGER(i4) :: elem, fill
    REAL(dp) :: weight
  END TYPE elem_desc_dp

  !> MPI datatype used to describe type part_desc
  INTEGER, SAVE :: mp_part_desc_i4, sizeof_part_desc_i4, &
       mp_elem_desc_i4, sizeof_elem_desc_i4, &
       mp_part_desc_dp, sizeof_part_desc_dp, &
       mp_elem_desc_dp, sizeof_elem_desc_dp

  INTERFACE repartition_swap_mp
    MODULE PROCEDURE repartition_swap_mp_i4
    MODULE PROCEDURE repartition_swap_mp_dp
  END INTERFACE repartition_swap_mp
  PUBLIC :: repartition_swap_mp
#endif

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_set_repartition.f90'
CONTAINS
#ifdef USE_MPI
  !> assumptions:
  !!  -# number of partitions equals number of callers
  !!  -# every caller passes one partition
  !! @param partition local members of global (implicit) set
  !! @param weight element weight, must be same size as partition
  SUBROUTINE repartition_swap_mp_i4(partition, weight, comm, &
       efficiency_threshold)
    INTEGER, INTENT(inout) :: partition(:)
    INTEGER(i4), INTENT(in) :: weight(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(sp), OPTIONAL, INTENT(in) :: efficiency_threshold

    INTEGER :: part_comm, num_parts, local_part_size, ierror, comm_rank
    TYPE(elem_desc_i4), ALLOCATABLE :: elems(:)
    INTEGER(i4), ALLOCATABLE :: rsort_part(:)
    INTEGER :: dummy, max_part_size, partner_part_size
    INTEGER :: i, swap_partner
    TYPE(swap_pair) :: flip
    TYPE(part_desc_i4), ALLOCATABLE :: part_desc(:)
    INTEGER(i8) :: weight_sum_high, weight_sum_low
    INTEGER(i4), ALLOCATABLE :: weight_high(:), weight_low(:)
    INTEGER :: num_low_weights, num_high_weights, ofs_low_weights, &
         ofs_high_weights
    TYPE(elem_desc_i4) :: t_elem
    TYPE(part_desc_i4) :: local_part_desc
    INTEGER(i4) :: t_weight
    LOGICAL :: sufficient_improval
    REAL(dp) :: mean_weight_sum_i
    REAL(sp) :: ethresholdi

    ethresholdi = 1.0_sp ! unrealistic optimum
    IF (PRESENT(efficiency_threshold)) &
         ethresholdi = 1.0_sp/efficiency_threshold

    local_part_size = SIZE(partition)
    CALL assertion(local_part_size == SIZE(weight), filename, __LINE__, &
         'arguments sizes must match')
    part_comm = ppm_default_comm
    IF (PRESENT(comm)) part_comm = comm
    ! simple serial implementation at this time
    CALL mpi_comm_size(part_comm, num_parts, ierror)
    CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    CALL mpi_comm_rank(part_comm, comm_rank, ierror)
    CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    ALLOCATE(part_desc(0:num_parts-1), rsort_part(0:num_parts-1))
    local_part_desc%weight_sum = SUM(INT(weight, i8))
    local_part_desc%part_size = local_part_size
    local_part_desc%part_rank = comm_rank
    sufficient_improval = .TRUE.
    CALL mpi_allgather(local_part_desc, 1, mp_part_desc_i4, &
         part_desc, 1, mp_part_desc_i4, part_comm, ierror)
    CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    max_part_size = 0
    mean_weight_sum_i = 0.0_dp
    DO i = 0, num_parts-1
      IF (part_desc(i)%part_size > max_part_size) &
           max_part_size = part_desc(i)%part_size
      mean_weight_sum_i = mean_weight_sum_i + REAL(part_desc(i)%weight_sum, dp)
    END DO
    mean_weight_sum_i = REAL(num_parts, dp)/mean_weight_sum_i
    ALLOCATE(elems(local_part_size + max_part_size), &
         weight_high(max_part_size), weight_low(max_part_size))
    DO i = 1, local_part_size
      elems(i)%weight = weight(i)
      elems(i)%elem = partition(i)
    END DO
    ! sort elem description by weight to ease swap lookup
    CALL qsort_r(elems(1:local_part_size), local_part_size, &
         sizeof_elem_desc_i4, dummy, rcmp_elem_desc_i4)
    swap_partner = 0
    DO WHILE (sufficient_improval)
      CALL qsort_r(part_desc, num_parts, sizeof_part_desc_i4, dummy, &
           rcmp_part_desc_i4)
      DO i = 0, num_parts-1
        rsort_part(part_desc(i)%part_rank) = i
      END DO
      sufficient_improval = .FALSE.
      IF (REAL(REAL(part_desc(0)%weight_sum, dp) * mean_weight_sum_i, sp) &
           <= ethresholdi) EXIT
      swap_partner = part_desc(num_parts - rsort_part(comm_rank) - 1)%part_rank
      IF (swap_partner /= comm_rank) THEN
        partner_part_size = part_desc(rsort_part(swap_partner))%part_size
        CALL mpi_sendrecv(elems(1:local_part_size), local_part_size, &
             mp_elem_desc_i4, swap_partner, &
             0, elems(local_part_size+1:local_part_size+partner_part_size), &
             partner_part_size, mp_elem_desc_i4, swap_partner, 0, &
             part_comm, mpi_status_ignore, ierror)
        CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
        IF (rsort_part(comm_rank) >= num_parts/2) THEN
          ! our partition is low weight
          weight_sum_high = part_desc(rsort_part(swap_partner))%weight_sum
          weight_sum_low = part_desc(rsort_part(comm_rank))%weight_sum
          weight_high(1:partner_part_size) = &
               elems(local_part_size+1:local_part_size+partner_part_size)%weight
          weight_low(1:local_part_size) &
               = elems(1:local_part_size)%weight
          num_low_weights = local_part_size
          num_high_weights = partner_part_size
          ofs_low_weights = 0
          ofs_high_weights = local_part_size
        ELSE
          ! our partition is high weight
          weight_sum_high = part_desc(rsort_part(comm_rank))%weight_sum
          weight_sum_low = part_desc(rsort_part(swap_partner))%weight_sum
          weight_high(1:local_part_size) &
               = elems(1:local_part_size)%weight
          num_high_weights = local_part_size
          weight_low(1:partner_part_size) = &
               elems(local_part_size+1:local_part_size+partner_part_size)%weight
          num_low_weights = partner_part_size
          num_high_weights = local_part_size
          ofs_low_weights = local_part_size
          ofs_high_weights = 0
        END IF
        DO WHILE(weight_sum_high /= weight_sum_low)
          flip = find_swap_pair_sorted_i4(weight_sum_high, weight_sum_low, &
               weight_high(1:num_high_weights), weight_low(1:num_low_weights))
          IF (flip%a /= 0) THEN
            weight_sum_low = weight_sum_low &
                 - INT(weight_low(flip%b), i8) &
                 + INT(weight_high(flip%a), i8)
            weight_sum_high = weight_sum_high &
                 - INT(weight_high(flip%a), i8) &
                 + INT(weight_low(flip%b), i8)
            t_elem = elems(ofs_high_weights+flip%a)
            elems(ofs_high_weights+flip%a) = &
                 elems(ofs_low_weights+flip%b)
            elems(ofs_low_weights+flip%b) = t_elem
            t_weight = weight_high(flip%a)
            weight_high(flip%a) = weight_low(flip%b)
            weight_low(flip%b) = t_weight
            sufficient_improval = .TRUE.
            CALL insertion_sort_once(elems(1:local_part_size), &
                 local_part_size, sizeof_elem_desc_i4, &
                 MERGE(flip%a, flip%b, ofs_high_weights == 0), &
                 dummy, rcmp_elem_desc_i4)
            CALL insertion_sort_once(elems(local_part_size+1:), &
                 partner_part_size, sizeof_elem_desc_i4, &
                 MERGE(flip%b, flip%a, ofs_high_weights == 0), &
                 dummy, rcmp_elem_desc_i4)
            CALL insertion_sort_once(weight_high, num_high_weights, &
                 4, flip%a, dummy, rcmp_i4)
            CALL insertion_sort_once(weight_low, num_low_weights, &
                 4, flip%b, dummy, rcmp_i4)
          ELSE
            EXIT
          END IF
        END DO
        IF (rsort_part(comm_rank) >= num_parts/2) THEN
          ! our partition is low weight
          local_part_desc%weight_sum = weight_sum_low
        ELSE
          ! our partition is high weight
          local_part_desc%weight_sum = weight_sum_high
        END IF
      END IF
      CALL mpi_bcast(sufficient_improval, 1, mpi_logical, &
           part_desc(num_parts-1)%part_rank, part_comm, ierror)
      CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
      CALL mpi_allgather(local_part_desc, 1, mp_part_desc_i4, part_desc, &
           1, mp_part_desc_i4, part_comm, ierror)
      CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    END DO
    partition(:) = elems(1:local_part_size)%elem

  END SUBROUTINE repartition_swap_mp_i4

  !> assumptions:
  !!  -# number of partitions equals number of callers
  !!  -# every caller passes one partition
  !! @param partition local members of global (implicit) set
  !! @param weight element weight, must be same size as partition
  SUBROUTINE repartition_swap_mp_dp(partition, weight, comm, &
       efficiency_threshold)
    INTEGER, INTENT(inout) :: partition(:)
    REAL(dp), INTENT(in) :: weight(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(sp), OPTIONAL, INTENT(in) :: efficiency_threshold

    INTEGER :: part_comm, num_parts, local_part_size, ierror, comm_rank
    TYPE(elem_desc_dp), ALLOCATABLE :: elems(:)
    INTEGER(i4), ALLOCATABLE :: rsort_part(:)
    INTEGER :: dummy, max_part_size, partner_part_size
    INTEGER :: i, swap_partner
    TYPE(swap_pair) :: flip
    TYPE(part_desc_dp), ALLOCATABLE :: part_desc(:)
    REAL(dp) :: weight_sum_high, weight_sum_low
    REAL(dp), ALLOCATABLE :: weight_high(:), weight_low(:)
    INTEGER :: num_low_weights, num_high_weights, ofs_low_weights, &
         ofs_high_weights
    TYPE(elem_desc_dp) :: t_elem
    TYPE(part_desc_dp) :: local_part_desc
    REAL(dp) :: t_weight
    LOGICAL :: sufficient_improval
    REAL(dp) :: mean_weight_sum_i
    REAL(sp) :: ethresholdi
    REAL(dp) :: epsilon

    ethresholdi = 1.0_sp ! unrealistic optimum
    IF (PRESENT(efficiency_threshold)) &
         ethresholdi = 1.0_sp/efficiency_threshold

    local_part_size = SIZE(partition)
    CALL assertion(local_part_size == SIZE(weight), filename, __LINE__, &
         'arguments sizes must match')
    part_comm = ppm_default_comm
    IF (PRESENT(comm)) part_comm = comm
    ! simple serial implementation at this time
    CALL mpi_comm_size(part_comm, num_parts, ierror)
    CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    CALL mpi_comm_rank(part_comm, comm_rank, ierror)
    CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    ALLOCATE(part_desc(0:num_parts-1), rsort_part(0:num_parts-1))
    local_part_desc%weight_sum = SUM(weight)
    local_part_desc%part_size = local_part_size
    local_part_desc%part_rank = comm_rank
    sufficient_improval = .TRUE.
    CALL mpi_allgather(local_part_desc, 1, mp_part_desc_dp, &
         part_desc, 1, mp_part_desc_dp, part_comm, ierror)
    CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    epsilon = MAX(MAXVAL(SPACING(2.0_dp * weight)), &
         &        SPACING(local_part_desc%weight_sum))
    CALL mpi_allreduce(mpi_in_place, epsilon, 1, mp_dp, mpi_max, part_comm, &
         ierror)
    CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    max_part_size = 0
    mean_weight_sum_i = 0.0_dp
    DO i = 0, num_parts-1
      IF (part_desc(i)%part_size > max_part_size) &
           max_part_size = part_desc(i)%part_size
      mean_weight_sum_i = mean_weight_sum_i + part_desc(i)%weight_sum
    END DO
    mean_weight_sum_i = REAL(num_parts, dp)/mean_weight_sum_i
    ALLOCATE(elems(local_part_size + max_part_size), &
         weight_high(max_part_size), weight_low(max_part_size))
    DO i = 1, local_part_size
      elems(i)%weight = weight(i)
      elems(i)%elem = partition(i)
    END DO
    ! sort elem description by weight to ease swap lookup
    CALL qsort_r(elems(1:local_part_size), local_part_size, &
         sizeof_elem_desc_dp, dummy, rcmp_elem_desc_dp)
    swap_partner = 0
    DO WHILE (sufficient_improval)
      CALL qsort_r(part_desc, num_parts, sizeof_part_desc_dp, dummy, &
           rcmp_part_desc_dp)
      DO i = 0, num_parts-1
        rsort_part(part_desc(i)%part_rank) = i
      END DO
      sufficient_improval = .FALSE.
      IF (REAL(part_desc(0)%weight_sum * mean_weight_sum_i, sp) &
           <= ethresholdi) EXIT
      swap_partner = part_desc(num_parts - rsort_part(comm_rank) - 1)%part_rank
      IF (swap_partner /= comm_rank) THEN
        partner_part_size = part_desc(rsort_part(swap_partner))%part_size
        CALL mpi_sendrecv(elems(1:local_part_size), local_part_size, &
             mp_elem_desc_dp, swap_partner, &
             0, elems(local_part_size+1:local_part_size+partner_part_size), &
             partner_part_size, mp_elem_desc_dp, swap_partner, 0, &
             part_comm, mpi_status_ignore, ierror)
        CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
        IF (rsort_part(comm_rank) >= num_parts/2) THEN
          ! our partition is low weight
          weight_sum_high = part_desc(rsort_part(swap_partner))%weight_sum
          weight_sum_low = part_desc(rsort_part(comm_rank))%weight_sum
          weight_high(1:partner_part_size) = &
               elems(local_part_size+1:local_part_size+partner_part_size)%weight
          weight_low(1:local_part_size) &
               = elems(1:local_part_size)%weight
          num_low_weights = local_part_size
          num_high_weights = partner_part_size
          ofs_low_weights = 0
          ofs_high_weights = local_part_size
        ELSE
          ! our partition is high weight
          weight_sum_high = part_desc(rsort_part(comm_rank))%weight_sum
          weight_sum_low = part_desc(rsort_part(swap_partner))%weight_sum
          weight_high(1:local_part_size) &
               = elems(1:local_part_size)%weight
          num_high_weights = local_part_size
          weight_low(1:partner_part_size) = &
               elems(local_part_size+1:local_part_size+partner_part_size)%weight
          num_low_weights = partner_part_size
          num_high_weights = local_part_size
          ofs_low_weights = local_part_size
          ofs_high_weights = 0
        END IF
        DO WHILE(ABS(weight_sum_high - weight_sum_low) > 0.0_dp)
          flip = find_swap_pair_sorted_dp(weight_sum_high, weight_sum_low, &
               weight_high(1:num_high_weights), weight_low(1:num_low_weights), &
               epsilon)
          IF (flip%a /= 0) THEN
            weight_sum_low = weight_sum_low &
                 - weight_low(flip%b) + weight_high(flip%a)
            weight_sum_high = weight_sum_high &
                 - weight_high(flip%a) + weight_low(flip%b)
            t_elem = elems(ofs_high_weights+flip%a)
            elems(ofs_high_weights+flip%a) = &
                 elems(ofs_low_weights+flip%b)
            elems(ofs_low_weights+flip%b) = t_elem
            t_weight = weight_high(flip%a)
            weight_high(flip%a) = weight_low(flip%b)
            weight_low(flip%b) = t_weight
            sufficient_improval = .TRUE.
            CALL insertion_sort_once(elems(1:local_part_size), &
                 local_part_size, sizeof_elem_desc_dp, &
                 MERGE(flip%a, flip%b, ofs_high_weights == 0), &
                 dummy, rcmp_elem_desc_dp)
            CALL insertion_sort_once(elems(local_part_size+1:), &
                 partner_part_size, sizeof_elem_desc_dp, &
                 MERGE(flip%b, flip%a, ofs_high_weights == 0), &
                 dummy, rcmp_elem_desc_dp)
            CALL insertion_sort_once(weight_high, num_high_weights, &
                 mp_dp_extent, flip%a, dummy, rcmp_dp)
            CALL insertion_sort_once(weight_low, num_low_weights, &
                 mp_dp_extent, flip%b, dummy, rcmp_dp)
          ELSE
            EXIT
          END IF
        END DO
        IF (rsort_part(comm_rank) >= num_parts/2) THEN
          ! our partition is low weight
          local_part_desc%weight_sum = weight_sum_low
        ELSE
          ! our partition is high weight
          local_part_desc%weight_sum = weight_sum_high
        END IF
      END IF
      CALL mpi_bcast(sufficient_improval, 1, mpi_logical, &
           part_desc(num_parts-1)%part_rank, part_comm, ierror)
      CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
      CALL mpi_allgather(local_part_desc, 1, mp_part_desc_dp, part_desc, &
           1, mp_part_desc_dp, part_comm, ierror)
      CALL handle_mpi_error(ierror, part_comm, filename, __LINE__)
    END DO
    partition(:) = elems(1:local_part_size)%elem

  END SUBROUTINE repartition_swap_mp_dp

  SUBROUTINE handle_mpi_error(ierror, comm, source, line)
    INTEGER, INTENT(in) :: ierror, comm, line
    CHARACTER, INTENT(in) :: source

    INTEGER :: msg_len, ierror_
    CHARACTER(len=mpi_max_error_string) :: msg

    IF (ierror /= MPI_SUCCESS) THEN
      CALL mpi_error_string(ierror, msg, msg_len, ierror_)
      CALL abort_ppm(msg(1:msg_len), source, line, comm)
    END IF
  END SUBROUTINE handle_mpi_error

#endif

  SUBROUTINE repartition_swap_i4(partitioning, weight, efficiency_threshold)
    TYPE(partition_vec), INTENT(inout) :: partitioning
    INTEGER(i4), INTENT(in) :: weight(:)
    REAL(sp), OPTIONAL, INTENT(in) :: efficiency_threshold

    INTEGER(i8), ALLOCATABLE :: weight_sums(:)
    INTEGER :: i, j, num_elems, lb_low, ub_low, num_parts, lb_high, ub_high, &
         swap_low, swap_high, swap_temp, swap_part
    INTEGER(i8) :: mean_weight_sum
    INTEGER, ALLOCATABLE :: sort_part(:)
    INTEGER(i4), ALLOCATABLE :: weight_s(:)
    INTEGER(i4) :: weight_temp
    INTEGER :: max_part_size, psize
    LOGICAL :: exch
    TYPE(swap_pair) :: flip
    REAL(dp) :: mean_weight_sum_i
    REAL(sp) :: ethresholdi, efficiencyi


    ethresholdi = 1.0_sp ! unrealistic optimum
    IF (PRESENT(efficiency_threshold)) &
         ethresholdi = 1.0_sp/efficiency_threshold
    num_elems = SIZE(partitioning%elements)
    num_parts = SIZE(partitioning%start) - 1
    CALL assertion(num_elems == SIZE(weight), filename, __LINE__, &
         'number of elements must match number of weights')

    max_part_size = 0
    DO j = LBOUND(partitioning%start, 1), UBOUND(partitioning%start, 1) - 1
      psize = part_size(partitioning, j)
      IF (psize > max_part_size) &
           max_part_size = psize
    END DO
    ALLOCATE(weight_s(num_elems), weight_sums(num_parts), sort_part(num_parts))

    CALL partition_weight_sums(weight_sums, partitioning, weight, &
         mean_weight_sum)
    mean_weight_sum_i = 1.0_dp/REAL(mean_weight_sum, dp)
    DO i = 1, num_parts
      sort_part(i) = i
    END DO

    ! re-order weights to match partitions
    DO i = 1, num_elems
      weight_s(i) = weight(partitioning%elements(i))
    END DO

    exch = .TRUE.
    CALL qsort_r(sort_part, num_parts, sizeof_integer, &
         weight_sums, rcmp_i4_indirect_i8)
    DO WHILE(exch)
      ! improve to bipolar heap later
      exch = .FALSE.
      efficiencyi = REAL(REAL(weight_sums(sort_part(1)), dp) &
           * mean_weight_sum_i, sp)
      IF (efficiencyi <= ethresholdi) EXIT

      candidate_search: DO j = num_parts,num_parts,-1
        lb_high = partitioning%start(sort_part(1))
        ub_high = partitioning%start(sort_part(1)+1) - 1
        lb_low = partitioning%start(sort_part(j))
        ub_low = partitioning%start(sort_part(j)+1) - 1

        flip = find_swap_pair_i4(weight_sums(sort_part(1)), &
             weight_sums(sort_part(j)), &
             weight_s(lb_high:ub_high), weight_s(lb_low:ub_low))
        IF (flip%a /= 0) THEN
          exch = .TRUE.
          swap_part = j
          swap_high = partitioning%start(sort_part(1)) + flip%a - 1
          swap_low = partitioning%start(sort_part(j)) + flip%b - 1
          EXIT candidate_search
        END IF
      END DO candidate_search
      IF (exch) THEN
        j = swap_part
        swap_part = sort_part(j)
        weight_sums(swap_part) = weight_sums(swap_part) &
             - INT(weight(partitioning%elements(swap_low)), i8) &
             + INT(weight(partitioning%elements(swap_high)), i8)
        CALL heap_elem_increase_sort(sort_part, num_parts, sizeof_integer, &
             j, weight_sums, cmp_i4_indirect_i8)
        weight_sums(sort_part(1)) = weight_sums(sort_part(1)) &
             - INT(weight(partitioning%elements(swap_high)), i8) &
             + INT(weight(partitioning%elements(swap_low)), i8)
        CALL heapify(sort_part, num_parts, sizeof_integer, &
             1, weight_sums, cmp_i4_indirect_i8)
        swap_temp = partitioning%elements(swap_high)
        partitioning%elements(swap_high) = partitioning%elements(swap_low)
        partitioning%elements(swap_low) = swap_temp
        weight_temp = weight_s(swap_high)
        weight_s(swap_high) = weight_s(swap_low)
        weight_s(swap_low) = weight_temp

        CALL heap_leaf_minimize(sort_part, num_parts, sizeof_integer, &
             weight_sums, cmp_i4_indirect_i8)
      END IF
    END DO

    DEALLOCATE(weight_sums, sort_part, weight_s)

  END SUBROUTINE repartition_swap_i4

  !> Finds optimal pair of weights to swap to balance two parts.
  !! This routine accepts unsorted parts and requires
  !! size(weights_a) * size(weights_b) comparisons. It is thus
  !! best used on small parts only.
  PURE FUNCTION find_swap_pair_i4(weight_sum_a, weight_sum_b, &
       weights_a, weights_b) &
       RESULT(flip)
    TYPE(swap_pair) :: flip
    INTEGER(i8), INTENT(in) :: weight_sum_a, weight_sum_b
    INTEGER(i4), INTENT(in) :: weights_a(:), weights_b(:)

    INTEGER(i8) :: balance_improve, t
    INTEGER :: i, j, size_a, size_b

    flip = swap_pair(0, 0)
    size_a = SIZE(weights_a)
    size_b = SIZE(weights_b)
    balance_improve = 0_i8
    DO j = 1, size_b
      DO i = 1, size_a
        t = ABS(weight_sum_a &
             &  - INT(weights_a(i), i8) &
             &  + INT(weights_b(j), i8) &
             &  - (weight_sum_b &
             &     - INT(weights_b(j), i8) &
             &     + INT(weights_a(i), i8))) &
             - ABS(weight_sum_a - weight_sum_b)

        IF (t < balance_improve) THEN
          flip%a = i
          flip%b = j
          balance_improve = t
        END IF
      END DO
    END DO
  END FUNCTION find_swap_pair_i4

  !> differs from \link find_swap_pair_i4\endlink in that the
  !! weights_[ab] arrays have to be sorted (from high to low) so
  !! that candidate search can be terminated early and carried out
  !! efficiently with bisection.
  RECURSIVE PURE FUNCTION find_swap_pair_sorted_i4(weight_sum_a, weight_sum_b, &
       weights_a, weights_b) &
       RESULT(flip)
    TYPE(swap_pair) :: flip
    INTEGER(i8), INTENT(in) :: weight_sum_a, weight_sum_b
    INTEGER(i4), INTENT(in) :: weights_a(:), weights_b(:)

    INTEGER(i8) :: balance_improve, t, optimum_improve, target_weight
    INTEGER :: i, j, size_a, size_b, il, ir, im

    size_a = SIZE(weights_a)
    size_b = SIZE(weights_b)
    IF (size_a == 0 .OR. size_b == 0) RETURN
    IF (size_a < size_b) THEN
      flip = find_swap_pair_sorted_i4(weight_sum_b, weight_sum_a, &
           weights_b, weights_a)
      i = flip%a
      flip%a = flip%b
      flip%b = i
    ELSE
      flip = swap_pair(0, 0)
      balance_improve = 0_i8
      ! this amount should be added to a (subtracted from b) for
      ! optimal rebalancing
      optimum_improve = weight_sum_b - weight_sum_a
      IF (optimum_improve /= 0_i8) THEN
        swap_search: DO j = 1, size_b
          target_weight = 2_i8 * INT(weights_b(j), i8) - optimum_improve
          il = 1
          ir = size_a + 1
          im = (size_a + 1) / 2
          DO WHILE (.TRUE.)
            ! compute change of a
            IF (2_i8 * INT(weights_a(im), i8) > target_weight) THEN
              il = im
              im = (im + ir)/2
              IF (im == ir - 1) EXIT
            ELSE IF (2_i8 * INT(weights_a(im), i8) < target_weight) THEN
              ir = im
              im = (im + il)/2
              IF (im == il) EXIT
            ELSE ! 2 * weights_a(im) == target_weight => optimal partner
              flip%a = im
              flip%b = j
              EXIT swap_search
            END IF
          END DO
          ! the best candidate must be one of im-1,im,im+1
          DO i = MAX(im-1, 1), MIN(im+1, size_a)
            t = 2_i8 * (-INT(weights_a(i), i8) + INT(weights_b(j), i8))
            IF (ABS(t - optimum_improve) &
                 < ABS(balance_improve - optimum_improve)) THEN
              flip%a = i
              flip%b = j
              balance_improve = t
            END IF
          END DO
        END DO swap_search
      END IF
    END IF
  END FUNCTION find_swap_pair_sorted_i4

  !> differs from \link find_swap_pair_dp\endlink in that the
  !! weights_[ab] arrays have to be sorted (from high to low) so
  !! that candidate search can be terminated early and carried out
  !! efficiently with bisection.
  RECURSIVE PURE FUNCTION find_swap_pair_sorted_dp(weight_sum_a, weight_sum_b, &
       weights_a, weights_b, epsilon) &
       RESULT(flip)
    TYPE(swap_pair) :: flip
    REAL(dp), INTENT(in) :: weight_sum_a, weight_sum_b, epsilon
    REAL(dp), INTENT(in) :: weights_a(:), weights_b(:)

    REAL(dp) :: balance_improve, t, optimum_improve, target_weight
    INTEGER :: i, j, size_a, size_b, il, ir, im

    size_a = SIZE(weights_a)
    size_b = SIZE(weights_b)
    IF (size_a == 0 .OR. size_b == 0) RETURN
    IF (size_a < size_b) THEN
      flip = find_swap_pair_sorted_dp(weight_sum_b, weight_sum_a, &
           weights_b, weights_a, epsilon)
      i = flip%a
      flip%a = flip%b
      flip%b = i
    ELSE
      flip = swap_pair(0, 0)
      balance_improve = 0.0_dp
      ! this amount should be added to a (subtracted from b) for
      ! optimal rebalancing
      optimum_improve = weight_sum_b - weight_sum_a
      IF (ABS(optimum_improve) > 0.0_dp) THEN
        swap_search: DO j = 1, size_b
          target_weight = 2.0_dp * weights_b(j) - optimum_improve
          il = 1
          ir = size_a + 1
          im = (size_a + 1) / 2
          DO WHILE (.TRUE.)
            ! compute change of a
            IF (2.0_dp * weights_a(im) > target_weight) THEN
              il = im
              im = (im + ir)/2
              IF (im == ir - 1) EXIT
            ELSE IF (2.0_dp * weights_a(im) < target_weight) THEN
              ir = im
              im = (im + il)/2
              IF (im == il) EXIT
            ELSE ! 2 * weights_a(im) == target_weight => optimal partner
              flip%a = im
              flip%b = j
              EXIT swap_search
            END IF
          END DO
          ! the best candidate must be one of im-1,im,im+1
          DO i = MAX(im-1, 1), MIN(im+1, size_a)
            t = 2.0_dp * (-weights_a(i) + weights_b(j))
            IF (ABS(t - optimum_improve) &
                 < ABS(balance_improve - optimum_improve) &
                 .AND. ABS(t - balance_improve) > epsilon) THEN
              flip%a = i
              flip%b = j
              balance_improve = t
            END IF
          END DO
        END DO swap_search
      END IF
    END IF
  END FUNCTION find_swap_pair_sorted_dp

  SUBROUTINE repartition_swap_dp(partitioning, weight, efficiency_threshold)
    TYPE(partition_vec), INTENT(inout) :: partitioning
    REAL(dp), INTENT(in) :: weight(:)
    REAL(sp), OPTIONAL, INTENT(in) :: efficiency_threshold

    REAL(dp), ALLOCATABLE :: weight_sums(:), weights_low(:), weights_high(:)
    INTEGER :: i, j, num_elems, lb_low, ub_low, num_parts, lb_high, ub_high, &
         swap_low, swap_high, swap_temp, swap_part, part_high, part_low, &
         psize, max_part_size
    REAL(dp) :: mean_weight_sum, mean_weight_sum_i
    REAL(dp) :: epsilon
    INTEGER, ALLOCATABLE :: sort_part(:)
    LOGICAL :: exch
    TYPE(swap_pair) :: flip
    REAL(sp) :: ethresholdi, efficiencyi

    ethresholdi = 1.0_sp ! unrealistic optimum
    IF (PRESENT(efficiency_threshold)) &
         ethresholdi = 1.0_sp/efficiency_threshold

    num_elems = SIZE(partitioning%elements)
    num_parts = SIZE(partitioning%start) - 1
    CALL assertion(num_elems == SIZE(weight), filename, __LINE__, &
         'number of elements must match number of weights')

    max_part_size = 0
    DO j = LBOUND(partitioning%start, 1), UBOUND(partitioning%start, 1) - 1
      psize = part_size(partitioning, j)
      IF (psize > max_part_size) &
           max_part_size = psize
    END DO
    ALLOCATE(weights_low(0:max_part_size-1), weights_high(0:max_part_size-1))
    ALLOCATE(weight_sums(num_parts), sort_part(num_parts))

    CALL partition_weight_sums(weight_sums, partitioning, weight, &
         mean_weight_sum)
    mean_weight_sum_i = 1.0_dp/mean_weight_sum
    DO i = 1, num_parts
      sort_part(i) = i
    END DO
    exch = .TRUE.
    epsilon = MAX(MAXVAL(SPACING(2.0_dp * weight)), &
         &        MAXVAL(SPACING(weight_sums)))
    CALL qsort_r(sort_part, num_parts, sizeof_integer, &
         weight_sums, rcmp_i4_indirect_dp)
    DO WHILE(exch)
      exch = .FALSE.
      part_high = sort_part(1)
      efficiencyi = REAL(weight_sums(part_high) * mean_weight_sum_i, sp)
      IF (efficiencyi <= ethresholdi) EXIT
      lb_high = partitioning%start(part_high)
      ub_high = partitioning%start(part_high+1) - 1
      candidate_search: DO j = num_parts,num_parts,-1
        part_low = sort_part(j)
        lb_high = partitioning%start(part_high)
        ub_high = partitioning%start(part_high+1) - 1
        lb_low = partitioning%start(part_low)
        ub_low = partitioning%start(part_low+1) - 1
        DO i = 0, ub_low-lb_low
          weights_low(i) = weight(partitioning%elements(lb_low+i))
        END DO
        DO i = 0, ub_high-lb_high
          weights_high(i) = weight(partitioning%elements(lb_high+i))
        END DO
        flip = find_swap_pair_dp(weight_sums(part_high), &
             weight_sums(part_low), weights_high(0:ub_high-lb_high), &
             weights_low(0:ub_low-lb_low), epsilon)
        IF (flip%a /= 0) THEN
          exch = .TRUE.
          swap_part = j
          swap_high = partitioning%start(part_high) + flip%a - 1
          swap_low = partitioning%start(part_low) + flip%b - 1
          EXIT candidate_search
        END IF
      END DO candidate_search
      IF (exch) THEN
        j = swap_part
        swap_part = sort_part(j)
        weight_sums(swap_part) = weight_sums(swap_part) &
             - weight(partitioning%elements(swap_low)) &
             + weight(partitioning%elements(swap_high))
        CALL heap_elem_increase_sort(sort_part, num_parts, sizeof_integer, &
             j, weight_sums, cmp_i4_indirect_dp)
        weight_sums(part_high) = weight_sums(part_high) &
             - weight(partitioning%elements(swap_high)) &
             + weight(partitioning%elements(swap_low))
        CALL heapify(sort_part, num_parts, sizeof_integer, &
             1, weight_sums, cmp_i4_indirect_dp)
        swap_temp = partitioning%elements(swap_high)
        partitioning%elements(swap_high) = partitioning%elements(swap_low)
        partitioning%elements(swap_low) = swap_temp
        CALL heap_leaf_minimize(sort_part, num_parts, sizeof_integer, &
             weight_sums, cmp_i4_indirect_dp)
      END IF
    END DO

    DEALLOCATE(weight_sums, sort_part, weights_low, weights_high)

  END SUBROUTINE repartition_swap_dp

  PURE FUNCTION find_swap_pair_dp(weight_sum_a, weight_sum_b, &
       weights_a, weights_b, epsilon) &
       RESULT(flip)
    TYPE(swap_pair) :: flip
    REAL(dp), INTENT(in) :: weight_sum_a, weight_sum_b, epsilon
    REAL(dp), INTENT(in) :: weights_a(:), weights_b(:)

    REAL(dp) :: balance_improve, t
    INTEGER :: i, j, size_a, size_b

    flip = swap_pair(0, 0)
    size_a = SIZE(weights_a)
    size_b = SIZE(weights_b)
    balance_improve = 0.0_dp
    DO j = 1, size_b
      DO i = 1, size_a
        t = ABS((weight_sum_a - weight_sum_b) &
             &  - 2.0_dp * weights_a(i) + 2.0_dp * weights_b(j)) &
             - ABS(weight_sum_a - weight_sum_b)
        IF (t < balance_improve .AND. ABS(t) > epsilon) THEN
          flip%a = i
          flip%b = j
          balance_improve = t
        END IF
      END DO
    END DO
  END FUNCTION find_swap_pair_dp

  SUBROUTINE initialize_set_repartition
#ifdef USE_MPI
    INTEGER :: ierror, block_lens(0:1)
    INTEGER(mpi_address_kind) :: offsets(0:1), lb, extent

    ! build data types for repartition_swap_mp_i4
    offsets(0) = 0_mpi_address_kind
    CALL mpi_type_get_extent(mp_i8, lb, offsets(1), ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    block_lens = (/ 1, 2 /)
    CALL mpi_type_create_struct(2, block_lens, offsets, (/ mp_i8, mp_i4 /), &
         mp_part_desc_i4, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    CALL mpi_type_get_extent(mp_part_desc_i4, lb, extent, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    sizeof_part_desc_i4 = INT(extent)
    CALL mpi_type_commit(mp_part_desc_i4, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)

    block_lens(0) = 2
    CALL mpi_type_create_struct(1, block_lens, offsets, (/ mp_i4 /), &
         mp_elem_desc_i4, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    CALL mpi_type_get_extent(mp_elem_desc_i4, lb, extent, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    sizeof_elem_desc_i4 = INT(extent)
    CALL mpi_type_commit(mp_elem_desc_i4, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)

    ! build data types for repartition_swap_mp_dp
    offsets(0) = 0_mpi_address_kind
    CALL mpi_type_get_extent(mp_dp, lb, offsets(1), ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    block_lens(0) = 1
    block_lens(1) = 2
    CALL mpi_type_create_struct(2, block_lens, offsets, (/ mp_dp, mp_i4 /), &
         mp_part_desc_dp, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    CALL mpi_type_get_extent(mp_part_desc_dp, lb, extent, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    sizeof_part_desc_dp = INT(extent)
    CALL mpi_type_commit(mp_part_desc_dp, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)

    block_lens(0) = 2
    block_lens(1) = 1
    offsets(0) = 0_mpi_address_kind
    CALL mpi_type_get_extent(mp_i4, lb, offsets(1), ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    offsets(1) = 2_mpi_address_kind * offsets(1)
    CALL mpi_type_create_struct(2, block_lens, offsets, (/ mp_i4, mp_dp /), &
         mp_elem_desc_dp, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    CALL mpi_type_get_extent(mp_elem_desc_dp, lb, extent, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
    sizeof_elem_desc_dp = INT(extent)
    CALL mpi_type_commit(mp_elem_desc_dp, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
#endif
  END SUBROUTINE initialize_set_repartition

  SUBROUTINE finalize_set_repartition
#ifdef USE_MPI
    INTEGER :: ierror

    ! deallocate data types for repartition_swap_mp_i4
    CALL mpi_type_free(mp_part_desc_i4, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)

    CALL mpi_type_free(mp_elem_desc_i4, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)

    ! deallocate data types for repartition_swap_mp_dp
    CALL mpi_type_free(mp_part_desc_dp, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)

    CALL mpi_type_free(mp_elem_desc_dp, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, filename, __LINE__)
#endif
  END SUBROUTINE finalize_set_repartition

#ifdef USE_MPI
  FUNCTION rcmp_part_desc_i4(a, b) RESULT(cmp)
    TYPE(part_desc_i4), INTENT(in) :: a, b

    INTEGER :: cmp

    cmp = MERGE(-1, MERGE(1, 0, a%weight_sum < b%weight_sum), &
         &      a%weight_sum > b%weight_sum)
  END FUNCTION rcmp_part_desc_i4

  FUNCTION rcmp_elem_desc_i4(a, b) RESULT(cmp)
    TYPE(elem_desc_i4), INTENT(in) :: a, b

    INTEGER :: cmp

    cmp = MERGE(-1, MERGE(1, 0, a%weight < b%weight), &
         &      a%weight > b%weight)
  END FUNCTION rcmp_elem_desc_i4

  FUNCTION rcmp_part_desc_dp(a, b) RESULT(cmp)
    TYPE(part_desc_dp), INTENT(in) :: a, b

    INTEGER :: cmp

    cmp = MERGE(-1, MERGE(1, 0, a%weight_sum < b%weight_sum), &
         &      a%weight_sum > b%weight_sum)
  END FUNCTION rcmp_part_desc_dp

  FUNCTION rcmp_elem_desc_dp(a, b) RESULT(cmp)
    TYPE(elem_desc_dp), INTENT(in) :: a, b

    INTEGER :: cmp

    cmp = MERGE(-1, MERGE(1, 0, a%weight < b%weight), &
         &      a%weight > b%weight)
  END FUNCTION rcmp_elem_desc_dp
#endif

END MODULE ppm_set_repartition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
