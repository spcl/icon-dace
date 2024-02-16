!>
!! @file scales_ppm.f90
!! @brief gathers facilities from all the other modules to
!!        provide an installable module interface file
!!
!! @copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
! Version: 1.0
! Keywords: library Fortran 90 interface
! Author: Thomas Jahns <jahns@dkrz.de>
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
MODULE scales_ppm
  ! we need iso_c_binding without only because configure will
  ! substitute ppm_metis_int according to the actual library type
  USE iso_c_binding
  USE ppm_base, ONLY: abort_ppm, ppm_default_comm, set_default_comm, &
       ppm_default_comm, assertion
  USE ppm_std_type_kinds, ONLY: ps, rs, pd, rd, pi4, pi8, sp, dp, &
       i4, i8, sizeof_integer
  USE ppm_compare, ONLY: cmp_i4, cmp_i4_indirect_i4, cmp_i4_indirect_i8, &
       cmp_i4_indirect_dp, &
       rcmp_i4, rcmp_i4_indirect_dp, rcmp_i4_indirect_i8, &
       cmp_i8, &
       cmp_dp, rcmp_dp
  USE ppm_random, ONLY: a_rand, a_randp, irand, irandp, &
       drand, drandp, frand, frandp, lrand
  USE ppm_sort, ONLY: qsort_r, qsort_r_mt, insertion_sort, sorted_insertion
  USE ppm_search, ONLY: bsearch_r, bsearch_el_r
  USE ppm_heap, ONLY: build_heap, heapify, is_heap, &
       heap_elem_increase_sort, heap_remove_top, &
       heap_leaf_minimize
  USE ppm_extents, ONLY: extent, extent_size, extent_end, &
       extent_start, rebased_extent, &
       extent_set_iinterval, extent_from_iinterval, char, iinterval, &
       ASSIGNMENT(=), OPERATOR(==), iinterval_from_extent, &
       extent_intersect, extents_do_intersect
#ifdef USE_MPI
  USE ppm_extents_mp, ONLY: extent_mp
#endif
  USE ppm_strided_extents, ONLY: char, OPERATOR(==), &
       strided_extent, extent_size, extent_start, extent_end
  USE ppm_set_partition_base, ONLY: partition_vec, set_i4, &
       partition_assignment, balance_of_max, OPERATOR(==), &
       OPERATOR(/=), ASSIGNMENT(=), part_size, partition_weight_sums, &
       assign_set_i4_2_pv, assign_pv_2_set_i4, block_decomposition
  USE ppm_set_repartition, ONLY: repartition_swap
#ifdef USE_MPI
  USE ppm_strided_extents, ONLY: subarray_mpi_datatype
  USE ppm_set_partition_base, ONLY: balance_of_max_mp
  USE ppm_set_repartition, ONLY: repartition_swap_mp
#endif /* USE_MPI */
  USE ppm_uniform_partition, ONLY: uniform_decomposition, &
       uniform_partition
#ifdef USE_PARMETIS
  USE ppm_graph_partition_mpi, ONLY: graph_partition_parmetis
#endif
#ifdef USE_METIS
  USE ppm_graph_partition_serial, ONLY: graph_partition_metis
#endif
  USE ppm_checksum, ONLY: init_digests, hex_checksum, ppm_md5, &
       ppm_sha1, hashes
  IMPLICIT NONE
#if defined(USE_PARMETIS) || defined(USE_METIS)
#include <ppm.inc>
#endif
#ifdef USE_PARMETIS
  PUBLIC :: graph_partition_parmetis
#endif
#ifdef USE_METIS
  PUBLIC :: graph_partition_metis
#endif
  PUBLIC :: abort_ppm, initialize_scales_ppm, finalize_scales_ppm
  PUBLIC :: char
CONTAINS

  SUBROUTINE initialize_scales_ppm(default_comm, random_seed, seed_output)
#ifdef USE_MPI
    USE ppm_std_type_kinds_mp, ONLY: create_types_mp
    USE ppm_extents_mp, ONLY: create_extents_mp
#endif
    USE ppm_math_extensions, ONLY: initialize_math_extensions
    USE ppm_f90_io_lun, ONLY: setup_lun_table
    USE ppm_set_repartition, ONLY: initialize_set_repartition
    USE ppm_random, ONLY: initialize_irand
    USE ppm_extents, ONLY: extents_init_io_formats => init_io_formats
    INTEGER, OPTIONAL, INTENT(in) :: default_comm, random_seed
    INTEGER, OPTIONAL, INTENT(out) :: seed_output

!$omp single
    IF (PRESENT(default_comm)) THEN
      ppm_default_comm = default_comm
      CALL set_default_comm(default_comm)
    END IF

    CALL setup_lun_table
    CALL initialize_math_extensions
    CALL extents_init_io_formats
#ifdef USE_MPI
    CALL create_types_mp
    CALL create_extents_mp
#endif
!$omp end single
    CALL initialize_irand(ppm_default_comm, random_seed, seed_output)
!$omp single
    CALL initialize_set_repartition
    CALL init_digests
!$omp end single
  END SUBROUTINE initialize_scales_ppm

  SUBROUTINE finalize_scales_ppm
#ifdef USE_MPI
    USE ppm_std_type_kinds_mp, ONLY: destroy_types_mp
    USE ppm_extents_mp, ONLY: destroy_extents_mp
#endif
    USE ppm_set_repartition, ONLY: finalize_set_repartition
    USE ppm_math_extensions, ONLY: finalize_math_extensions
    USE ppm_f90_io_lun, ONLY: take_down_lun_table
    USE ppm_random, ONLY: finalize_irand
!$omp single
    CALL finalize_set_repartition
!$omp end single
    CALL finalize_irand
!$omp single
#ifdef USE_MPI
    CALL destroy_extents_mp
    CALL destroy_types_mp
#endif
    CALL finalize_math_extensions
    CALL take_down_lun_table
!$omp end single
  END SUBROUTINE finalize_scales_ppm

END MODULE scales_ppm
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
