!> @file ppm_distributed.f90
!! @brief distributed graph data structure
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords: mpi distributed data structure
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
!> distributed data structures base module
#include "fc_feature_defs.inc"
MODULE ppm_distributed
  USE ppm_base, ONLY: assertion, abort_ppm
  USE ppm_std_type_kinds, ONLY: i4
  USE ppm_std_type_kinds_mp, ONLY: mp_i4
  USE ppm_graph_csr, ONLY: graph_csr, num_nodes, num_edges
  USE ppm_extents, ONLY: extent, extent_size
  USE ppm_rectilinear, ONLY: rlcoord2lidx, lidx_nb_indices, &
       num_neighbours_of_rect_elem, lidx2rlcoord
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  !> proxy object for distributed graph-structured data
  TYPE graph_csr_dist_i4
    !> number of nodes in global graph (or -1 if not yet computed)
    INTEGER(i4) :: num_vertices_g
    !> number of edges in global graph (or -1 if not yet computed)
    INTEGER(i4) :: num_edges_g
    !> locally present nodes are kept here with edge connectivity
    !! referring to global node indices
    TYPE(graph_csr) :: edges_l
    !> mapping of local node indices to global node index
    INTEGER(i4), ALLOCATABLE :: idxl2idxg(:)
  END TYPE graph_csr_dist_i4

  !> construct graph from rectilinear slices
  INTERFACE build_graph
    MODULE PROCEDURE build_graph_from_rect_dist_i4
  END INTERFACE build_graph

  !> number of edges
  INTERFACE num_edges
    MODULE PROCEDURE num_edges_csr_dist_i4_mp
    MODULE PROCEDURE num_edges_csr_dist_i4
  END INTERFACE num_edges

  !> number of vertices
  INTERFACE num_nodes
    MODULE PROCEDURE num_nodes_csr_dist_i4_mp
    MODULE PROCEDURE num_nodes_csr_dist_i4
  END INTERFACE num_nodes

  PUBLIC :: graph_csr_dist_i4, build_graph, num_nodes, num_edges

  !> gather distributed data structure on one process in MPI program
  INTERFACE graph_gather
    MODULE PROCEDURE gather_graph_dist_i4_mp
  END INTERFACE graph_gather
  PUBLIC :: graph_gather

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_distributed.f90'
CONTAINS
  ! TODO: missing functionality: holes and extra edges
  ! TODO: multiple rect_l
  !> Construct graph of distributed global rectilinear where
  !! the calling process has data for one rectilinear section.
  !!
  !! @param graph_g proxy object for local representation of
  !! shared global graph
  !! @param rect_g distributed rectilinear extents
  !! @param rect_l local rectilinear part
  SUBROUTINE build_graph_from_rect_dist_i4(graph_g, rect_g, rect_l)
    TYPE(graph_csr_dist_i4), INTENT(out) :: graph_g
    TYPE(extent), INTENT(in) :: rect_g(:), rect_l(:)

    INTEGER(i4) :: num_vertices_l, rlcoord(SIZE(rect_g)), nnb
    INTEGER :: dims, i

    dims = SIZE(rect_g)
    CALL assertion(dims == SIZE(rect_l), filename, __LINE__, &
         "distributed rectilinear descriptors must have same size")
    num_vertices_l = extent_size(rect_l)
    ! fixme: assert local rect is contained in global rect
    ALLOCATE(graph_g%edges_l%edges_of_vtx(num_vertices_l+1), &
         graph_g%idxl2idxg(num_vertices_l))
    graph_g%edges_l%edges_of_vtx(1) = 1
    DO i = 1, num_vertices_l
      rlcoord = lidx2rlcoord(rect_l, i)
      nnb = num_neighbours_of_rect_elem(rect_g, rlcoord)
      graph_g%idxl2idxg(i) = rlcoord2lidx(rect_g, rlcoord)
      graph_g%edges_l%edges_of_vtx(i+1) = nnb
    END DO
    DO i = 1, num_vertices_l
      graph_g%edges_l%edges_of_vtx(i+1) = &
           graph_g%edges_l%edges_of_vtx(i) &
           + graph_g%edges_l%edges_of_vtx(i+1)
    END DO
    ALLOCATE(graph_g%edges_l%edges(&
         graph_g%edges_l%edges_of_vtx(num_vertices_l + 1) - 1))
    DO i = 1, num_vertices_l
      CALL lidx_nb_indices(rect_g, graph_g%idxl2idxg(i), &
           graph_g%edges_l%edges(&
           graph_g%edges_l%edges_of_vtx(i):&
           graph_g%edges_l%edges_of_vtx(i+1)-1))
    END DO
    ! evaluate these lazily if requested
    graph_g%num_vertices_g = -1
    graph_g%num_edges_g = -1
  END SUBROUTINE build_graph_from_rect_dist_i4

  !> the communicator can be omitted if num_nodes has already been
  !! called collectively (i.e. with a comm object) and the cached
  !! result is already present
  !!
  !! @param csr_dist proxy object for distributed graph data structure
  !! @return number of nodes of graph in total
  FUNCTION num_nodes_csr_dist_i4(csr_dist) RESULT(n)
    TYPE(graph_csr_dist_i4), INTENT(in) :: csr_dist
    INTEGER(i4) :: n

    CALL assertion(csr_dist%num_vertices_g /= -1, filename, __LINE__, &
         'num_vertices_g must have previously been evaluated')

    n = csr_dist%num_vertices_g
  END FUNCTION num_nodes_csr_dist_i4

  !> use this if num_nodes has not yet been called
  !!
  !! @param csr_dist proxy object for distributed graph data structure
  !! @param comm communicator handle for processes sharing the data structure
  !! @return number of nodes of graph in total
  FUNCTION num_nodes_csr_dist_i4_mp(csr_dist, comm) RESULT(n)
    TYPE(graph_csr_dist_i4), INTENT(inout) :: csr_dist
    INTEGER, INTENT(in) :: comm
    INTEGER(i4) :: n

    INTEGER :: ierror

    IF (csr_dist%num_vertices_g == -1) THEN
      CALL mpi_allreduce(num_nodes(csr_dist%edges_l), n, 1, mp_i4, mpi_sum, &
           comm, ierror)
      IF (ierror /= mpi_success) &
           CALL abort_ppm("mpi_allreduce failed", filename, __LINE__)
      csr_dist%num_vertices_g = n
    ELSE
      n = csr_dist%num_vertices_g
    END IF
  END FUNCTION num_nodes_csr_dist_i4_mp

  !> use this if num_edges has not yet been called collectively
  !!
  !! @param csr_dist proxy object for distributed graph data structure
  !! @param comm communicator handle for processes sharing the data structure
  !! @return number of edges of graph in total
  FUNCTION num_edges_csr_dist_i4_mp(csr_dist, comm) RESULT(n)
    TYPE(graph_csr_dist_i4), INTENT(inout) :: csr_dist
    INTEGER, INTENT(in) :: comm
    INTEGER(i4) :: n

    INTEGER :: ierror

    IF (csr_dist%num_edges_g == -1) THEN
      CALL mpi_allreduce(num_edges(csr_dist%edges_l), n, 1, mp_i4, mpi_sum, &
           comm, ierror)
      IF (ierror /= mpi_success) &
           CALL abort_ppm("mpi_allreduce failed", filename, __LINE__)
      csr_dist%num_edges_g = n
    ELSE
      n = csr_dist%num_edges_g
    END IF
  END FUNCTION num_edges_csr_dist_i4_mp

  !> query number of edges in distributed graph
  !! the communicator can be omitted IF num_edges has already been
  !! called collectively (i.e. with a comm object) and the cached
  !! result is already present
  !!
  !! @param csr_dist proxy object for distributed graph data structure
  !! @return number of edges of graph in total
  FUNCTION num_edges_csr_dist_i4(csr_dist) RESULT(n)
    TYPE(graph_csr_dist_i4), INTENT(in) :: csr_dist
    INTEGER(i4) :: n

    CALL assertion(csr_dist%num_edges_g /= -1, filename, __LINE__, &
         'num_edges_g must have previously been evaluated')

    n = csr_dist%num_edges_g
  END FUNCTION num_edges_csr_dist_i4

  !> gather distributed graph object into local representation
  !!
  !! @param graph_g_gather local destination object on \c dest
  !! @param graph_g_dist proxy object of graph distributed over
  !! processes in \c comm
  !! @param dest process to gather graph on
  !! @param comm communicator handle for processes sharing \c graph_g_dist
  SUBROUTINE gather_graph_dist_i4_mp(graph_g_dist, dest, comm, graph_g_gather)
    TYPE(graph_csr_dist_i4), INTENT(in) :: graph_g_dist
    INTEGER, INTENT(in) :: dest, comm
    TYPE(graph_csr), OPTIONAL, INTENT(out) :: graph_g_gather

    INTEGER :: comm_gather, comm_size, comm_rank, ierror, part_desc_size
    INTEGER(i4) :: num_vertices_l, num_edges_l, num_other_vertices_l, &
         num_vertices_g, num_edges_g, node_edges, node_idx_g
    INTEGER :: i, j, p, q
    INTEGER(i4), ALLOCATABLE :: part_sizes(:,:), buf(:), edge_buf(:)

    CALL mpi_comm_dup(comm, comm_gather, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_dup failed", filename, __LINE__)
    CALL mpi_comm_rank(comm_gather, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_rank failed", filename, __LINE__)
    CALL mpi_comm_size(comm_gather, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_comm_size failed", filename, __LINE__)
    ALLOCATE(part_sizes(2, 0:MERGE(0, comm_size - 1, comm_rank /= dest)), &
         buf(2))
    num_vertices_l = SIZE(graph_g_dist%edges_l%edges_of_vtx) - 1
    num_edges_l = graph_g_dist%edges_l%edges_of_vtx(num_vertices_l+1) - 1
    buf(1) = num_vertices_l
    buf(2) = num_edges_l
    CALL mpi_gather(buf, 2, mp_i4, &
         part_sizes, 2, mp_i4, dest, comm_gather, ierror)
    IF (ierror /= mpi_success) &
         CALL abort_ppm("mpi_gather failed", filename, __LINE__)
    IF (comm_rank == dest) THEN
      CALL assertion(PRESENT(graph_g_gather), filename, __LINE__, &
           'output object must be present on dest process')
      num_vertices_g = SUM(part_sizes(1,:))
      num_edges_g = SUM(part_sizes(2,:))
      ALLOCATE(graph_g_gather%edges_of_vtx(num_vertices_g+1), &
           graph_g_gather%edges(num_edges_g))
      ALLOCATE(edge_buf(2*num_vertices_g+num_edges_g))
      p = 1
      DO j = 0, comm_size - 1
        IF (j /= dest) THEN
          part_desc_size = 2 * part_sizes(1, j) + part_sizes(2, j)
          DEALLOCATE(buf)
          ALLOCATE(buf(part_desc_size))
          CALL mpi_recv(buf, part_desc_size, mp_i4, j, 1, comm_gather, &
               mpi_status_ignore, ierror)
          IF (ierror /= mpi_success) &
               CALL abort_ppm("mpi_recv failed", filename, __LINE__)
          num_other_vertices_l = part_sizes(1, j)
          q = 1
          DO i = 1, num_other_vertices_l
            node_edges = buf(num_other_vertices_l+i)
            graph_g_gather%edges_of_vtx(buf(i)) = node_edges
            edge_buf(p) = buf(i)
            edge_buf(p+1) = node_edges
            IF (node_edges > 0) THEN
              edge_buf(p+2:p+1+node_edges) &
                   = buf(2 * num_other_vertices_l + q &
                   &     :2 * num_other_vertices_l + q - 1 + node_edges)
            END IF
            q = q + node_edges
            p = p + 2 + node_edges
          END DO
        ELSE
          q = 1
          DO i = 1, num_vertices_l
            node_idx_g = graph_g_dist%idxl2idxg(i)
            node_edges = graph_g_dist%edges_l%edges_of_vtx(i+1) &
                 & - graph_g_dist%edges_l%edges_of_vtx(i)
            graph_g_gather%edges_of_vtx(node_idx_g) = node_edges
            edge_buf(p) = node_idx_g
            edge_buf(p+1) = node_edges
            IF (node_edges > 0) THEN
              edge_buf(p+2:p+1+node_edges) &
                   = graph_g_dist%edges_l%edges(q:q+node_edges-1)
            END IF
            p = p + 2 + node_edges
            q = q + node_edges
          END DO
        END IF
      END DO
      p = graph_g_gather%edges_of_vtx(1)
      graph_g_gather%edges_of_vtx(1) = 1
      DO i = 2, num_vertices_g+1
        q = graph_g_gather%edges_of_vtx(i)
        graph_g_gather%edges_of_vtx(i) = p + graph_g_gather%edges_of_vtx(i - 1)
        p = q
      END DO
      p = 1
      DO j = 1, num_vertices_g
        i = edge_buf(p)
        q = graph_g_gather%edges_of_vtx(i)
        node_edges = graph_g_gather%edges_of_vtx(i + 1) &
             &       - graph_g_gather%edges_of_vtx(i)
        CALL assertion(node_edges == edge_buf(p + 1))
        IF (node_edges > 0) THEN
          graph_g_gather%edges(q:q+node_edges-1) = &
               edge_buf(p + 2:p + 1 + node_edges)
        END IF
        p = p + 2 + node_edges
      END DO
    ELSE
      DEALLOCATE(buf)
      part_desc_size = 2 * num_vertices_l + num_edges_l
      ALLOCATE(buf(part_desc_size))
      buf(1:num_vertices_l) = graph_g_dist%idxl2idxg
      buf(num_vertices_l+1:2*num_vertices_l) &
           = graph_g_dist%edges_l%edges_of_vtx(2:num_vertices_l+1) &
           - graph_g_dist%edges_l%edges_of_vtx(1:num_vertices_l)
      buf(2*num_vertices_l+1:part_desc_size) = graph_g_dist%edges_l%edges
      CALL mpi_send(buf, part_desc_size, mp_i4, dest, 1, comm_gather, ierror)
      IF (ierror /= mpi_success) &
           CALL abort_ppm("mpi_send failed", filename, __LINE__)
    ENDIF
    DEALLOCATE(part_sizes, buf)
    CALL mpi_comm_free(comm_gather, ierror)
  END SUBROUTINE gather_graph_dist_i4_mp

END MODULE ppm_distributed
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
