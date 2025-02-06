!>
!! @file ppm_graph_partition_serial.f90
!! @brief wrapper for single task graph partitioner
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!
! Version: 1.0
! Keywords: graph partitioning
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
!> perform partitioning of graph from serial code
#include "fc_feature_defs.inc"
MODULE ppm_graph_partition_serial
  USE iso_c_binding
  USE ppm_base, ONLY: assertion
  USE ppm_extents, ONLY: extent
  USE ppm_graph_csr, ONLY: graph_csr, num_nodes, num_edges
  USE ppm_set_partition_base, ONLY: partition_assignment
  USE ppm_std_type_kinds, ONLY: i4, i8
  IMPLICIT NONE
  PRIVATE
#include <ppm.inc>
#ifdef HAVE_METIS_V4
  INTEGER, PARAMETER :: ppm_metis_int = c_int
#else
  INTEGER, PARAMETER :: ppm_metis_int = ppm_metis_idx
#endif
#ifdef HAVE_METIS_V4
  INTERFACE
    SUBROUTINE metis_mcpartgraphkway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, &
         wgtflag, numflag, nparts, rubvec, options, edgecut, part) BIND(C)
      USE iso_c_binding, ONLY: c_int, c_ptr
      IMPORT :: ppm_metis_idx
      IMPORT :: ppm_metis_real
      INTEGER(ppm_metis_idx), INTENT(in) :: xadj(*), adjncy(*)
      TYPE(c_ptr), VALUE, INTENT(in) :: vwgt, adjwgt
      INTEGER(c_int), INTENT(in) :: nvtxs, wgtflag, numflag, ncon, nparts, options(*)
      REAL(ppm_metis_real), INTENT(in) :: rubvec(ncon)
      INTEGER(c_int), INTENT(out) :: edgecut
      INTEGER(ppm_metis_idx), INTENT(out) :: part(*)
    END SUBROUTINE metis_mcpartgraphkway
  END INTERFACE
#else
  INTERFACE
    SUBROUTINE metis_setdefaultoptions(options) BIND(C)
      IMPORT :: ppm_metis_idx
      INTEGER(ppm_metis_idx), INTENT(out) :: options(*)
    END SUBROUTINE metis_setdefaultoptions
  END INTERFACE
#endif
  INTERFACE
#ifdef HAVE_METIS_V4
    SUBROUTINE metis_partgraphkway(nvtxs, xadj, adjncy, vwgt, adjwgt, &
         wgtflag, numflag, nparts, options, edgecut, part) BIND(C)
      USE iso_c_binding, ONLY: c_int
#else
    SUBROUTINE metis_partgraphkway(nvtxs, ncon, xadj, adjncy, vwgt, vsize, &
         adjwgt, nparts, tpwgts, ubvec, options, edgecut, part) BIND(C)
#endif
      USE iso_c_binding, ONLY: c_ptr
      IMPORT :: ppm_metis_idx, ppm_metis_int
      INTEGER(ppm_metis_idx), INTENT(in) :: xadj(*), adjncy(*)
      TYPE(c_ptr), VALUE, INTENT(in) :: vwgt, adjwgt
      INTEGER(ppm_metis_int), INTENT(in) :: nvtxs, nparts, options(*)
#ifdef HAVE_METIS_V4
      INTEGER(c_int), INTENT(in) :: wgtflag, numflag
#else
      INTEGER(ppm_metis_idx), INTENT(in) :: ncon
      TYPE(c_ptr), VALUE, INTENT(in) :: vsize
      TYPE(c_ptr), VALUE, INTENT(in) :: tpwgts, ubvec
#endif
      INTEGER(ppm_metis_int), INTENT(out) :: edgecut
      INTEGER(ppm_metis_idx), INTENT(out) :: part(*)
    END SUBROUTINE metis_partgraphkway
  END INTERFACE
  PUBLIC :: graph_partition_metis
  INTERFACE graph_partition_metis
    MODULE PROCEDURE graph_partition_metis_base
    MODULE PROCEDURE graph_partition_metis_csr_nowgt
    MODULE PROCEDURE graph_partition_metis_csr_i4
    MODULE PROCEDURE graph_partition_metis_csr_i8
    MODULE PROCEDURE graph_partition_metis_csr_mweights_i4
    MODULE PROCEDURE graph_partition_metis_csr_mweights_i8
  END INTERFACE graph_partition_metis
  PUBLIC :: graph_partition_metis_csr

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_graph_partition_serial.f90'
CONTAINS
#define UNUSED(x) IF (SIZE( (/(x)/) ) < 0) CONTINUE
  SUBROUTINE graph_partition_metis_base(num_vertices, edge_list_lens, &
       edge_lists, partition_out, num_partitions, &
       imbalance_tolerance, vertex_weights, edge_weights)
    INTEGER(ppm_metis_int), INTENT(in) :: num_vertices
    INTEGER(ppm_metis_idx), INTENT(in) :: edge_list_lens(*)
    INTEGER(ppm_metis_idx), INTENT(in) :: edge_lists(*)
    INTEGER(ppm_metis_idx), INTENT(out) :: partition_out(*)
    INTEGER(ppm_metis_int), INTENT(in) :: num_partitions
    REAL(ppm_metis_real), OPTIONAL, TARGET, INTENT(in) :: imbalance_tolerance(:)
    INTEGER(ppm_metis_idx), OPTIONAL, TARGET, INTENT(in) :: vertex_weights(*)
    INTEGER(ppm_metis_idx), OPTIONAL, TARGET, INTENT(in) :: edge_weights(*)
    TYPE(c_ptr) :: vwgt, adjwgt
#ifdef HAVE_METIS_V4
    INTEGER(c_int) :: wgtflag
    INTEGER(c_int) :: metis_options(0:4), edge_cut

    metis_options = 0
#else
    TYPE(c_ptr) :: vsize = c_null_ptr
    TYPE(c_ptr) :: tpwgts = c_null_ptr, ubvec
    INTEGER(ppm_metis_idx) :: metis_options(0:39), edge_cut

    CALL metis_setdefaultoptions(metis_options)
    ! METIS_OPTION_NUMBERING : use Fortran-style
    metis_options(17) = 1
#endif

    IF (PRESENT(vertex_weights)) THEN
      vwgt = c_loc(vertex_weights(1))
    ELSE
      vwgt = c_null_ptr
    END IF
    IF (PRESENT(edge_weights)) THEN
      adjwgt = c_loc(edge_weights(1))
    ELSE
      adjwgt = c_null_ptr
    END IF
    IF (PRESENT(imbalance_tolerance)) THEN
      CALL assertion(PRESENT(vertex_weights), line=__LINE__, source=filename, &
           msg="when imbalance_tolerance is provided, vertex weights&
           & are also required")
#ifdef HAVE_METIS_V4
      wgtflag = MERGE(1, 0, PRESENT(edge_weights))
      CALL metis_mcpartgraphkway(num_vertices, &
           INT(SIZE(imbalance_tolerance), c_int), edge_list_lens, edge_lists, &
           vwgt, adjwgt, wgtflag, 1, num_partitions, imbalance_tolerance, &
           metis_options, edge_cut, partition_out)
#else
      ubvec = c_loc(imbalance_tolerance(1))
#endif
    ELSE
#ifdef HAVE_METIS_V4
      wgtflag = MERGE(2, 0, PRESENT(vertex_weights))
      wgtflag = IOR(wgtflag, MERGE(1, 0, PRESENT(edge_weights)))
      CALL metis_partgraphkway(num_vertices, edge_list_lens, edge_lists, vwgt, &
           adjwgt, wgtflag, 1, num_partitions, metis_options, edge_cut, &
           partition_out)
#else
      ubvec = c_null_ptr
#endif
    END IF
#ifndef HAVE_METIS_V4
    CALL metis_partgraphkway(num_vertices, 1_ppm_metis_idx, edge_list_lens, &
         edge_lists, vwgt, vsize, adjwgt, num_partitions, tpwgts, ubvec, &
         metis_options, edge_cut, partition_out)
#endif
  END SUBROUTINE graph_partition_metis_base

  !> This function is meant for binary compatibility with previous releases
  !! that only supported a single integer kind for weights
  SUBROUTINE graph_partition_metis_csr(partition, graph, num_partitions, &
       vertex_weights, edge_weights)
    TYPE(partition_assignment), INTENT(out) :: partition
    TYPE(graph_csr), INTENT(in) :: graph
    INTEGER, INTENT(in) :: num_partitions
    INTEGER(ppm_metis_idx), OPTIONAL, INTENT(in) :: vertex_weights(*)
    INTEGER(ppm_metis_idx), OPTIONAL, INTENT(in) :: edge_weights(*)
    INTEGER(ppm_metis_idx) :: kind_select

    IF (PRESENT(vertex_weights) .OR. PRESENT(edge_weights)) THEN
      CALL graph_partition_metis(partition, graph, num_partitions, &
           kind_select, vertex_weights, edge_weights)
    ELSE
      CALL graph_partition_metis(partition, graph, num_partitions)
    END IF
  END SUBROUTINE graph_partition_metis_csr

  SUBROUTINE graph_partition_metis_csr_nowgt(partition, graph, num_partitions)
    TYPE(partition_assignment), INTENT(out) :: partition
    TYPE(graph_csr), INTENT(in) :: graph
    INTEGER, INTENT(in) :: num_partitions

    INTEGER :: n
    INTEGER(ppm_metis_int) :: n_, num_partitions_
#if METIS_IDX_KIND_MATCH_I8
    INTEGER :: i
    INTEGER(ppm_metis_idx), ALLOCATABLE :: xadj(:), adjncy(:), part_assign(:)
#endif

    n = num_nodes(graph)
    n_ = INT(n, ppm_metis_int)
    num_partitions_ = INT(num_partitions, ppm_metis_int)
    partition%part_range = extent(1, num_partitions)
    ALLOCATE(partition%assigned(1:n))
#if METIS_IDX_KIND_MATCH_I4
    CALL graph_partition_metis_base(n_, graph%edges_of_vtx, &
         graph%edges, partition%assigned, num_partitions_)
#else
    ALLOCATE(xadj(SIZE(graph%edges_of_vtx)), adjncy(SIZE(graph%edges)), &
         part_assign(n))
    xadj = INT(graph%edges_of_vtx, ppm_metis_idx)
    adjncy = INT(graph%edges, ppm_metis_idx)
    CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
         num_partitions_)
    DO i = 1, n
      partition%assigned(i) = INT(part_assign(i), i4)
    END DO
#endif
  END SUBROUTINE graph_partition_metis_csr_nowgt

  SUBROUTINE graph_partition_metis_csr_i4(partition, graph, num_partitions, &
       wgt_kind, vertex_weights, edge_weights)
    TYPE(partition_assignment), INTENT(out) :: partition
    TYPE(graph_csr), INTENT(in) :: graph
    INTEGER, INTENT(in) :: num_partitions
    !> intent is none, because the contents are not referencec but rather used
    !! for interface disambiguation
    INTEGER(i4) :: wgt_kind
    INTEGER(i4), OPTIONAL, INTENT(in) :: vertex_weights(*)
    INTEGER(i4), OPTIONAL, INTENT(in) :: edge_weights(*)

    INTEGER :: n, e, wgtflag
    INTEGER(ppm_metis_int) :: n_, num_partitions_
#if METIS_IDX_KIND_MATCH_I8
    INTEGER :: i
    INTEGER(ppm_metis_idx), ALLOCATABLE :: xadj(:), adjncy(:), &
         part_assign(:), vwgt(:), ewgt(:)
#endif

    UNUSED(wgt_kind)
    n = num_nodes(graph)
    e = num_edges(graph)
    n_ = INT(n, ppm_metis_int)
    num_partitions_ = INT(num_partitions, ppm_metis_int)
    wgtflag = MERGE(2, 0, PRESENT(vertex_weights))
    wgtflag = IOR(wgtflag, MERGE(1, 0, PRESENT(edge_weights)))
    partition%part_range = extent(1, num_partitions)
#if METIS_IDX_KIND_MATCH_I4
    ALLOCATE(partition%assigned(n))
    CALL graph_partition_metis_base(n_, graph%edges_of_vtx, &
         graph%edges, partition%assigned, num_partitions_, &
         vertex_weights=vertex_weights, edge_weights=edge_weights)
#else
    ALLOCATE(xadj(n), adjncy(e), part_assign(n))
    IF (PRESENT(vertex_weights)) THEN
      ALLOCATE(vwgt(n))
      DO i = 1, n
        vwgt(i) = INT(vertex_weights(i), ppm_metis_idx)
      END DO
    END IF
    IF (PRESENT(edge_weights)) THEN
      ALLOCATE(ewgt(e))
      DO i = 1, e
        ewgt(i) = INT(edge_weights(i), ppm_metis_idx)
      END DO
    END IF
    DO i = 1, n
      xadj(i) = INT(graph%edges_of_vtx(i), ppm_metis_idx)
    END DO
    DO i = 1, e
      adjncy(i) = INT(graph%edges(i), ppm_metis_idx)
    END DO
    SELECT CASE (wgtflag)
    CASE(0)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_)
    CASE(1)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, edge_weights=ewgt)
    CASE(2)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, vertex_weights=vwgt)
    CASE(3)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, vertex_weights=vwgt, edge_weights=ewgt)
    END SELECT
    DEALLOCATE(xadj, adjncy)
    ALLOCATE(partition%assigned(n))
    DO i = 1, n
      partition%assigned(i) = INT(part_assign(i), i4)
    END DO
#endif
  END SUBROUTINE graph_partition_metis_csr_i4

  SUBROUTINE graph_partition_metis_csr_i8(partition, graph, num_partitions, &
       wgt_kind, vertex_weights, edge_weights)
    TYPE(partition_assignment), INTENT(out) :: partition
    TYPE(graph_csr), INTENT(in) :: graph
    INTEGER, INTENT(in) :: num_partitions
    !> intent is none, because the contents are not referencec but rather used
    !! for interface disambiguation
    INTEGER(i8) :: wgt_kind
    INTEGER(i8), OPTIONAL, INTENT(in) :: vertex_weights(*)
    INTEGER(i8), OPTIONAL, INTENT(in) :: edge_weights(*)

    INTEGER :: i, n, e, wgtflag
    INTEGER(ppm_metis_idx), ALLOCATABLE :: xadj(:), adjncy(:), part_assign(:)
    INTEGER(ppm_metis_int) :: n_, num_partitions_
#if METIS_IDX_KIND_MATCH_I4
    INTEGER(ppm_metis_idx), ALLOCATABLE :: vwgt(:), ewgt(:)
#endif

    UNUSED(wgt_kind)
    n = num_nodes(graph)
    e = num_edges(graph)
    n_ = INT(n, ppm_metis_int)
    num_partitions_ = INT(num_partitions, ppm_metis_int)
    wgtflag = MERGE(2, 0, PRESENT(vertex_weights))
    wgtflag = IOR(wgtflag, MERGE(1, 0, PRESENT(edge_weights)))
    partition%part_range = extent(1, num_partitions)
    ALLOCATE(xadj(n), adjncy(e), part_assign(n))
    DO i = 1, n
      xadj(i) = INT(graph%edges_of_vtx(i), ppm_metis_idx)
    END DO
    DO i = 1, e
      adjncy(i) = INT(graph%edges(i), ppm_metis_idx)
    END DO
#if METIS_IDX_KIND_MATCH_I8
    CALL graph_partition_metis_base(n_, xadj, adjncy, &
         part_assign, num_partitions_, &
         vertex_weights=vertex_weights, edge_weights=edge_weights)
#else
    IF (PRESENT(vertex_weights)) THEN
      ALLOCATE(vwgt(n))
      DO i = 1, n
        vwgt(i) = INT(vertex_weights(i), ppm_metis_idx)
      END DO
    END IF
    IF (PRESENT(edge_weights)) THEN
      ALLOCATE(ewgt(e))
      DO i = 1, e
        ewgt(i) = INT(edge_weights(i), ppm_metis_idx)
      END DO
    END IF
    SELECT CASE (wgtflag)
    CASE(0)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_)
    CASE(1)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, edge_weights=ewgt)
    CASE(2)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, vertex_weights=vwgt)
    CASE(3)
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, vertex_weights=vwgt, edge_weights=ewgt)
    END SELECT
#endif
    DEALLOCATE(xadj, adjncy)
    ALLOCATE(partition%assigned(n))
    DO i = 1, n
      partition%assigned(i) = INT(part_assign(i), i4)
    END DO
  END SUBROUTINE graph_partition_metis_csr_i8

  SUBROUTINE graph_partition_metis_csr_mweights_i4(partition, graph, &
       num_partitions, imbalance_tolerance, vertex_weights, edge_weights)
    TYPE(partition_assignment), INTENT(out) :: partition
    TYPE(graph_csr), TARGET, INTENT(in) :: graph
    INTEGER, INTENT(in) :: num_partitions
    REAL(ppm_metis_real), INTENT(in) :: imbalance_tolerance(:)
    INTEGER(i4), INTENT(in) :: vertex_weights(:,:)
    INTEGER(i4), OPTIONAL, INTENT(in) :: edge_weights(*)

    INTEGER :: n
    INTEGER(ppm_metis_int) :: n_, num_partitions_
#if METIS_IDX_KIND_MATCH_I8
    INTEGER :: i, j, e, m
    INTEGER(ppm_metis_idx), ALLOCATABLE :: xadj(:), adjncy(:), part_assign(:), &
         vwgt(:,:), ewgt(:)
#endif

    n = num_nodes(graph)
    n_ = INT(n, ppm_metis_int)
    num_partitions_ = INT(num_partitions, ppm_metis_idx)
    partition%part_range = extent(1, num_partitions)
#if METIS_IDX_KIND_MATCH_I4
    ALLOCATE(partition%assigned(n))
    IF (PRESENT(edge_weights)) THEN
      CALL graph_partition_metis_base(n_, graph%edges_of_vtx, &
           graph%edges, partition%assigned, &
           num_partitions_, imbalance_tolerance, vertex_weights, edge_weights)
    ELSE
      CALL graph_partition_metis_base(n_, graph%edges_of_vtx, &
           graph%edges, partition%assigned, &
           num_partitions_, imbalance_tolerance, vertex_weights)
    END IF
#elif METIS_IDX_KIND_MATCH_I8
    e = num_edges(graph)
    m = SIZE(vertex_weights, 1)
    ALLOCATE(xadj(n), adjncy(e), part_assign(n), vwgt(m, n))
    DO i = 1, n
      xadj(i) = INT(graph%edges_of_vtx(i), ppm_metis_idx)
    END DO
    adjncy = INT(graph%edges, ppm_metis_idx)
    DO j = 1, n
      DO i = 1, m
        vwgt(i, j) = INT(vertex_weights(i, j), ppm_metis_idx)
      END DO
    END DO
    IF (PRESENT(edge_weights)) THEN
      e = num_edges(graph)
      ALLOCATE(ewgt(e))
      DO i = 1, e
        ewgt(i) = INT(edge_weights(i), ppm_metis_idx)
      END DO
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, imbalance_tolerance, vwgt, ewgt)
      DEALLOCATE(ewgt)
    ELSE
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, imbalance_tolerance, vwgt)
    END IF
    DEALLOCATE(xadj, adjncy, vwgt)
    ALLOCATE(partition%assigned(n))
    DO i = 1, n
      partition%assigned(i) = INT(part_assign(i), i4)
    END DO
#endif
  END SUBROUTINE graph_partition_metis_csr_mweights_i4

  SUBROUTINE graph_partition_metis_csr_mweights_i8(partition, graph, &
       num_partitions, imbalance_tolerance, vertex_weights, edge_weights)
    TYPE(partition_assignment), INTENT(out) :: partition
    TYPE(graph_csr), TARGET, INTENT(in) :: graph
    INTEGER, INTENT(in) :: num_partitions
    REAL(ppm_metis_real), INTENT(in) :: imbalance_tolerance(:)
    INTEGER(i8), INTENT(in) :: vertex_weights(:,:)
    INTEGER(i8), OPTIONAL, INTENT(in) :: edge_weights(*)

    INTEGER :: e, i, n
    INTEGER(ppm_metis_int) :: n_, num_partitions_

#if METIS_IDX_KIND_MATCH_I4
    INTEGER :: j, m
    INTEGER(ppm_metis_idx), ALLOCATABLE :: vwgt(:,:), ewgt(:)
#elif METIS_IDX_KIND_MATCH_I8
    INTEGER(ppm_metis_idx), ALLOCATABLE :: xadj(:), adjncy(:), part_assign(:)
#endif

    n = num_nodes(graph)
    n_ = INT(n, ppm_metis_int)
    num_partitions_ = INT(num_partitions, ppm_metis_int)
    partition%part_range = extent(1, num_partitions)
#if METIS_IDX_KIND_MATCH_I4
    m = SIZE(vertex_weights, 1)
    ALLOCATE(partition%assigned(n), vwgt(m, n))
    DO j = 1, n
      DO i = 1, m
        vwgt(i, j) = INT(vertex_weights(i, j), ppm_metis_idx)
      END DO
    END DO
    IF (PRESENT(edge_weights)) THEN
      e = num_edges(graph)
      ALLOCATE(ewgt(e))
      DO i = 1, e
        ewgt(i) = INT(edge_weights(i), ppm_metis_idx)
      END DO
      CALL graph_partition_metis_base(n_, graph%edges_of_vtx, &
           graph%edges, partition%assigned, &
           num_partitions_, imbalance_tolerance, vwgt, ewgt)
    ELSE
      CALL graph_partition_metis_base(n_, graph%edges_of_vtx, &
           graph%edges, partition%assigned, &
           num_partitions_, imbalance_tolerance, vwgt)
    END IF
#elif METIS_IDX_KIND_MATCH_I8
    e = num_edges(graph)
    ALLOCATE(xadj(n), adjncy(e), part_assign(n))
    DO i = 1, n
      xadj(i) = INT(graph%edges_of_vtx(i), ppm_metis_idx)
    END DO
    DO i = 1, e
      adjncy(i) = INT(graph%edges(i), ppm_metis_idx)
    END DO
    IF (PRESENT(edge_weights)) THEN
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, imbalance_tolerance, vertex_weights, edge_weights)
    ELSE
      CALL graph_partition_metis_base(n_, xadj, adjncy, part_assign, &
           num_partitions_, imbalance_tolerance, vertex_weights)
    END IF
    DEALLOCATE(xadj, adjncy)
    ALLOCATE(partition%assigned(n))
    DO i = 1, SIZE(partition%assigned)
      partition%assigned(i) = INT(part_assign(i), i4)
    END DO
#endif
  END SUBROUTINE graph_partition_metis_csr_mweights_i8

END MODULE ppm_graph_partition_serial
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
