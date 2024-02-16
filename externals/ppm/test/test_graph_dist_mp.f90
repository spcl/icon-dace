!>
!! @file test_graph_dist_mp.f90
!! @brief test distributed graph data structure
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
PROGRAM test_graph_dist_mp
  USE ppm_base, ONLY: abort_ppm, assertion
  USE ppm_std_type_kinds, ONLY: i4
  USE ppm_combinatorics, ONLY: prime_factorization
  USE ppm_distributed, ONLY: graph_csr_dist_i4, build_graph, &
       graph_gather, num_edges, num_nodes
  USE ppm_extents, ONLY: extent, iinterval
  USE ppm_extents_mp, ONLY: extent_mp
  USE ppm_graph_csr, ONLY: graph_csr, build_graph, OPERATOR(==), &
       num_edges, num_nodes
  USE ppm_random, ONLY: a_randr, irandr
  USE ppm_rectilinear, ONLY: lidx2rlcoord
  USE random_data, ONLY: random_rectilinear
  USE ppm_uniform_partition, ONLY: uniform_partition
  USE scales_ppm, ONLY: initialize_scales_ppm, finalize_scales_ppm
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  TYPE(graph_csr) :: graph_g_ref, graph_g_gathered
  TYPE(graph_csr_dist_i4) :: graph_g

  INTEGER :: ierror, comm_rank, comm_size
  INTEGER :: ndims, i
  INTEGER(i4), PARAMETER :: rect_size_lim = 10
  INTEGER(i4) :: random_seed
  INTEGER(i4), ALLOCATABLE :: factors(:), deco_coord(:)
  INTEGER, ALLOCATABLE :: factor_assign(:)
  TYPE(extent), ALLOCATABLE :: rect_g_shape(:), rect_l_shape(:), &
       rect_g_decompose(:)
  CHARACTER(len=*), PARAMETER :: filename = 'test_graph_dist_mp.f90'
  ! init mpi
  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_init failed", filename, __LINE__)
  ! init scales ppm
  CALL initialize_scales_ppm(random_seed=0, seed_output=random_seed)
  ! find rank and number of tasks
  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_comm_rank failed", filename, __LINE__)
  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_comm_size failed", filename, __LINE__)
  PRINT '(i0,a,i0)', comm_rank, ': random_seed=', random_seed
  IF (comm_rank == 0) THEN
    CALL prime_factorization(comm_size, factors)
    ndims = irandr(iinterval(1, MAX(1, MIN(SIZE(factors), 7))))
  END IF
  CALL mpi_bcast(ndims, 1, mpi_integer, 0, mpi_comm_world, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_bcast failed", filename, __LINE__)
  ALLOCATE(rect_g_decompose(ndims), rect_g_shape(ndims), rect_l_shape(ndims))
  IF (comm_rank == 0) THEN
    ALLOCATE(factor_assign(SIZE(factors)))
    CALL a_randr(factor_assign, iinterval(1, ndims))
    rect_g_decompose = extent(1, 1)
    DO i = 1, SIZE(factors)
      rect_g_decompose(factor_assign(i))%size &
           = rect_g_decompose(factor_assign(i))%size * factors(i)
    END DO
    CALL random_rectilinear(rect_g_shape, rect_size_lim)
    DEALLOCATE(factor_assign, factors)
  END IF
  CALL mpi_bcast(rect_g_decompose, ndims, extent_mp, 0, mpi_comm_world, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_bcast failed", filename, __LINE__)
  CALL mpi_bcast(rect_g_shape, ndims, extent_mp, 0, mpi_comm_world, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_bcast failed", filename, __LINE__)
  ! compute each slice of uniform decomposition
  ALLOCATE(deco_coord(ndims))
  deco_coord = lidx2rlcoord(rect_g_decompose, comm_rank + 1)
  rect_l_shape = uniform_partition(rect_g_shape, rect_g_decompose%size, &
       deco_coord)
  DEALLOCATE(deco_coord)
  ! build distributed graph parts
  CALL build_graph(graph_g, rect_g_shape, rect_l_shape)
  ! TODO: gather distributed graph
  CALL graph_gather(graph_g, 0, mpi_comm_world, graph_g_gathered)
  ! build global graph on rank 0
  CALL build_graph(graph_g_ref, rect_g_shape)
  ! compare gathered to global graph,
  ! test global properties of graph_g vs. graph_g_ref
  CALL assertion(num_nodes(graph_g, mpi_comm_world) == num_nodes(graph_g_ref), &
       filename, __LINE__, 'number of nodes is inequal')
  CALL assertion(num_edges(graph_g, mpi_comm_world) == num_edges(graph_g_ref), &
       filename, __LINE__, 'number of edges is inequal')
  ! test for equality of graph_g_gathered vs. graph_g_ref
  IF (comm_rank == 0) THEN
    CALL assertion(graph_g_gathered == graph_g_ref, filename, __LINE__, &
         'equality of gathered distributed graph and local full graph failed')
  END IF
  CALL finalize_scales_ppm
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_finalize failed", filename, __LINE__)
  DEALLOCATE(rect_g_decompose, rect_g_shape, rect_l_shape)
END PROGRAM test_graph_dist_mp
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
