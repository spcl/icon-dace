!> @file ppm_graph_csr.f90 --- CSR data structure
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <Thomas.Jahns@gmx.net>
!!
!! @version 1.0
!! @author Thomas Jahns <Thomas.Jahns@gmx.net>
! Keywords: CSR compressed sparse row graph
! Maintainer: Thomas Jahns <Thomas.Jahns@gmx.net>
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
!> data structure for representation of graph in csr format

MODULE ppm_graph_csr
  USE ppm_base, ONLY: assertion, abort_ppm
  USE ppm_combinatorics, ONLY: is_permutation
  USE ppm_extents, ONLY: extent, extent_size, iinterval, ASSIGNMENT(=)
  USE ppm_strio, ONLY: sscana, get_address, fmt_elem, arg_i4
  USE ppm_rectilinear, ONLY: lidx2rlcoord, num_neighbours_of_rect_elem, &
       lidx_nb_indices
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE
  PRIVATE

  TYPE graph_csr
    !> Vertex i has edges(edges_of_vtx(i):edges_of_vtx(i+1)-1
    !! also referred to as xadj in MeTiS
    INTEGER, ALLOCATABLE :: edges_of_vtx(:)
    !> End vertices of graph edges
    !! also referred to as adjncy in MeTiS
    INTEGER, ALLOCATABLE :: edges(:)
  END TYPE graph_csr

  !> number of edges
  INTERFACE num_edges
    MODULE PROCEDURE num_edges_csr
  END INTERFACE num_edges

  !> number of vertices
  INTERFACE num_nodes
    MODULE PROCEDURE num_nodes_csr
  END INTERFACE num_nodes

  !> serialize graph to file
  INTERFACE write_graph
    MODULE PROCEDURE write_graph_csr_nwmo_ewo
  END INTERFACE write_graph

  !> construct graph from rectilinear or adjacency matrix
  INTERFACE build_graph
    MODULE PROCEDURE graph_csr_from_rect_i4
    MODULE PROCEDURE graph_csr_from_adj_matrix_i4
  END INTERFACE build_graph

  !> construct graph from rectilinear in multi-threaded program
  INTERFACE build_graph_mt
    MODULE PROCEDURE graph_csr_from_rect_mt_i4
    MODULE PROCEDURE graph_csr_from_irect_mt_i4
  END INTERFACE build_graph_mt

  !> compare whether two graphs are equal
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE graph_csr_equal_i4
  END INTERFACE OPERATOR(==)

  !> @return whether graph is symmetric
  INTERFACE graph_is_symmetric
    MODULE PROCEDURE graph_is_symmetric_csr
    MODULE PROCEDURE graph_is_symmetric_csr_ew
  END INTERFACE graph_is_symmetric

  PUBLIC :: num_edges, num_nodes, OPERATOR(==)
  PUBLIC :: graph_csr, build_graph, build_graph_mt
  PUBLIC :: read_graph, write_graph
  PUBLIC :: assign_edge_weight, assign_symmetric_edge_weight
  PUBLIC :: graph_is_symmetric

  INTEGER, PUBLIC, PARAMETER :: graph_io_undirected=1

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_graph_csr.f90'
CONTAINS

  PURE FUNCTION num_edges_csr(csr)
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER :: num_edges_csr

    num_edges_csr = SIZE(csr%edges)
  END FUNCTION num_edges_csr

  PURE FUNCTION num_nodes_csr(csr)
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER :: num_nodes_csr

    num_nodes_csr = SIZE(csr%edges_of_vtx) - 1
  END FUNCTION num_nodes_csr

  !> Construct CSR graph object from adjacency matrix
  !! @param csr graph representation will be written to this argument
  !! @param adj for each edge from node u to node v, adj(v, u) must be true
  !! @param assert_undirected if given and .true., check symmetry of adj
  !! @param node_offset enter all node references with a bias of
  !! node_offset, defaults to 1
  SUBROUTINE graph_csr_from_adj_matrix_i4(csr, adj, assert_undirected, &
       node_offset)
    TYPE(graph_csr), INTENT(inout) :: csr
    LOGICAL, INTENT(in) :: adj(:, :)
    LOGICAL, OPTIONAL, INTENT(in) :: assert_undirected
    INTEGER, OPTIONAL, INTENT(in) :: node_offset

    INTEGER :: i, j, k, m, ofs

    m = SIZE(adj, 1)
    ofs = 1
    IF (PRESENT(node_offset)) ofs = node_offset
    CALL assertion(m == SIZE(adj, 2), filename, 141, &
         'adjacency matrix must be rectangular')
!$omp single
    IF (ALLOCATED(csr%edges_of_vtx)) DEALLOCATE(csr%edges_of_vtx)
    ALLOCATE(csr%edges_of_vtx(ofs:ofs+m))
    csr%edges_of_vtx(ofs) = 1
    IF (ALLOCATED(csr%edges)) DEALLOCATE(csr%edges)
!$omp end single
!$omp do
    DO j = 1, m
      csr%edges_of_vtx(j + ofs) = COUNT(adj(:, j))
    END DO
!$omp end do
!$omp single
    DO j = 1, m
      csr%edges_of_vtx(ofs + j) = csr%edges_of_vtx(ofs + j) &
           + csr%edges_of_vtx(ofs + j - 1)
    END DO
    ALLOCATE(csr%edges(csr%edges_of_vtx(m + ofs) - 1))
!$omp end single
!$omp do
    DO j = 1, m
      k = csr%edges_of_vtx(ofs + j - 1)
      DO i = 1, m
        IF (adj(i, j)) THEN
          csr%edges(k) = i + ofs - 1
          k = k + 1
        END IF
      END DO
      CALL assertion(csr%edges_of_vtx(ofs + j) == k, filename, 170, &
           'error in building adjacency lists')
    END DO
!$omp end do nowait
    IF (PRESENT(assert_undirected)) THEN
!$omp do
      DO j = 1, m
        DO i = j + 1, m
          IF (adj(i, j) .NEQV. adj(j, i)) THEN
            CALL abort_ppm('input matrix must be symmetric', filename, &
                 180)
          END IF
        END DO
      END DO
!$omp end do
    END IF
  END SUBROUTINE graph_csr_from_adj_matrix_i4

  SUBROUTINE graph_csr_from_rect_i4(csr, rect)
    TYPE(graph_csr), INTENT(inout) :: csr
    TYPE(extent), INTENT(in) :: rect(:)

    INTEGER :: i, node_ub, node_lb

    node_lb = 1
    node_ub = node_lb + extent_size(rect) - 1
    CALL csr_rect_alloc(csr, rect)
    DO i = node_lb, node_ub
      CALL lidx_nb_indices(rect, i, &
           csr%edges(csr%edges_of_vtx(i):csr%edges_of_vtx(i+1)-1))
    END DO
  END SUBROUTINE graph_csr_from_rect_i4

  SUBROUTINE graph_csr_from_rect_mt_i4(csr, rect)
    TYPE(graph_csr), INTENT(inout) :: csr
    TYPE(extent), INTENT(in) :: rect(:)

    INTEGER :: i, node_ub, node_lb

    node_lb = 1
    node_ub = node_lb + extent_size(rect) - 1
!$omp single
    CALL csr_rect_alloc(csr, rect)
!$omp end single
!$omp do
    DO i = node_lb, node_ub
      CALL lidx_nb_indices(rect, i, &
           csr%edges(csr%edges_of_vtx(i):csr%edges_of_vtx(i+1)-1))
    END DO
!$omp end do
  END SUBROUTINE graph_csr_from_rect_mt_i4

  SUBROUTINE graph_csr_from_irect_mt_i4(csr, rect)
    TYPE(graph_csr), INTENT(inout) :: csr
    TYPE(iinterval), INTENT(in) :: rect(:)

    TYPE(extent) :: rect_e(SIZE(rect))

    rect_e = rect
    CALL graph_csr_from_rect_mt_i4(csr, rect_e)
  END SUBROUTINE graph_csr_from_irect_mt_i4

  SUBROUTINE csr_rect_alloc(csr, rect)
    TYPE(graph_csr), INTENT(inout) :: csr
    TYPE(extent), INTENT(in) :: rect(:)

    INTEGER :: i, node_ub, node_lb, accum
    INTEGER :: rlcoord(SIZE(rect))

    node_lb = 1
    node_ub = node_lb + extent_size(rect) - 1

    IF (ALLOCATED(csr%edges_of_vtx)) DEALLOCATE(csr%edges_of_vtx)
    ALLOCATE(csr%edges_of_vtx(node_lb:node_ub+1))
    csr%edges_of_vtx(node_lb) = 1
    IF (ALLOCATED(csr%edges)) DEALLOCATE(csr%edges)

    accum = 1
    DO i = node_lb, node_ub
      rlcoord = lidx2rlcoord(rect, i)
      accum = accum + num_neighbours_of_rect_elem(rect, rlcoord)
      csr%edges_of_vtx(i + 1) = accum
    END DO
    ALLOCATE(csr%edges(csr%edges_of_vtx(node_ub + 1) - 1))
  END SUBROUTINE csr_rect_alloc

  SUBROUTINE write_graph_csr_plain(lun, csr, node_bias, flags)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER, OPTIONAL, INTENT(in) :: node_bias, flags

    INTEGER :: i, last, num, ofs, n_bias, eflags
    CHARACTER(len=11+1+9) :: fmt

    eflags = 0
    IF (PRESENT(flags)) eflags = flags

    ofs = LBOUND(csr%edges_of_vtx, 1)
    n_bias = -ofs + 1
    IF (PRESENT(node_bias)) n_bias = node_bias
    last = UBOUND(csr%edges_of_vtx, 1) - 1

    IF (IAND(eflags, graph_io_undirected) /= 0) THEN
      num = SIZE(csr%edges)/2
    ELSE
      num = SIZE(csr%edges)
    END IF
    WRITE (unit=lun, fmt='(i0,1x,i0)') SIZE(csr%edges_of_vtx) - 1, num
    DO i = ofs, last
      num = csr%edges_of_vtx(i + 1) - csr%edges_of_vtx(i)
      IF (num > 1) THEN
        WRITE(fmt, '(a,i0,a)') '(', num-1, '(i0,1x)i0)'
      ELSE
        fmt = '(i0)'
      END IF
      WRITE(unit=lun, fmt=fmt) &
           csr%edges(csr%edges_of_vtx(i):csr%edges_of_vtx(i+1)-1) + n_bias
    END DO
  END SUBROUTINE write_graph_csr_plain

  SUBROUTINE write_graph_csr_nw(lun, csr, node_weights, node_bias, flags)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER(i4), INTENT(in) :: node_weights(:)
    INTEGER, OPTIONAL, INTENT(in) :: node_bias, flags

    INTEGER :: i, last, num, ofs, n_bias, eflags, weight_fmt
    CHARACTER(len=11+1+9) :: fmt

    eflags = 0
    IF (PRESENT(flags)) eflags = flags

    ofs = LBOUND(csr%edges_of_vtx, 1)
    n_bias = -ofs + 1
    IF (PRESENT(node_bias)) n_bias = node_bias
    last = UBOUND(csr%edges_of_vtx, 1) - 1

    IF (IAND(eflags, graph_io_undirected) /= 0) THEN
      num = SIZE(csr%edges)/2
    ELSE
      num = SIZE(csr%edges)
    END IF
    weight_fmt = 10

    CALL assertion(SIZE(node_weights) == num_nodes(csr), filename, &
           315, 'number of node weights must equal number of nodes')

    WRITE (unit=lun, fmt='(i0,2(1x,i0))') SIZE(csr%edges_of_vtx) - 1, &
         num, weight_fmt
    DO i = ofs, last
      num = csr%edges_of_vtx(i + 1) - csr%edges_of_vtx(i)
      IF (num > 1) THEN
        WRITE(fmt, '(a,i0,a)') '(', num, '(i0,1x)i0)'
      ELSE
        fmt = '(i0,1x,i0)'
      END IF
      WRITE(unit=lun, fmt=fmt) &
           node_weights(i), &
           csr%edges(csr%edges_of_vtx(i):csr%edges_of_vtx(i+1)-1) + n_bias
    END DO
  END SUBROUTINE write_graph_csr_nw

  SUBROUTINE write_graph_csr_ew(lun, csr, edge_weights, node_bias, flags)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER(i4), INTENT(in) :: edge_weights(:)
    INTEGER, OPTIONAL, INTENT(in) :: node_bias, flags

    INTEGER :: i, last, lb, ub, num, ofs, n_bias, eflags, weight_fmt
    CHARACTER(len=11+1+9) :: fmt
    TYPE edge
      INTEGER :: node
      INTEGER(i4) :: weight
    END TYPE edge
    TYPE(edge), ALLOCATABLE :: edge_weight_vec(:)

    eflags = 0
    IF (PRESENT(flags)) eflags = flags

    ofs = LBOUND(csr%edges_of_vtx, 1)
    n_bias = -ofs + 1
    IF (PRESENT(node_bias)) n_bias = node_bias
    last = UBOUND(csr%edges_of_vtx, 1) - 1

    IF (IAND(eflags, graph_io_undirected) /= 0) THEN
      num = SIZE(csr%edges)/2
    ELSE
      num = SIZE(csr%edges)
    END IF
    weight_fmt = 1

    CALL assertion(SIZE(edge_weights) == num_edges(csr), filename, &
           362, 'number of edge weights must equal number of edges')

    WRITE (unit=lun, fmt='(i0,2(1x,i0))') SIZE(csr%edges_of_vtx) - 1, &
         num, weight_fmt
    ALLOCATE(edge_weight_vec(10))
    DO i = ofs, last
      lb = csr%edges_of_vtx(i)
      ub = csr%edges_of_vtx(i+1)-1
      num = csr%edges_of_vtx(i + 1) - csr%edges_of_vtx(i)
      IF (num > 0) THEN
        WRITE(fmt, '(a,i0,a)') '(', 2 * num - 1, '(i0,1x)i0)'
      END IF
      IF (num > SIZE(edge_weight_vec)) THEN
        DEALLOCATE(edge_weight_vec)
        ALLOCATE(edge_weight_vec(num))
      END IF
      edge_weight_vec(1:num)%node = csr%edges(lb:ub) + n_bias
      edge_weight_vec(1:num)%weight = edge_weights(lb:ub)
      WRITE(unit=lun, fmt=fmt) edge_weight_vec(1:num)
    END DO
    DEALLOCATE(edge_weight_vec)
  END SUBROUTINE write_graph_csr_ew

  SUBROUTINE write_graph_csr_nwm(lun, csr, node_weights, &
       node_bias, flags)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER(i4), INTENT(in) :: node_weights(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: node_bias, flags

    INTEGER :: i, last, edge_count, ofs, n_bias, eflags, weight_fmt, nnw, &
         adj_count, node_count
    CHARACTER(len=11+1+9) :: fmt

    eflags = 0
    IF (PRESENT(flags)) eflags = flags

    ofs = LBOUND(csr%edges_of_vtx, 1)
    n_bias = -ofs + 1
    IF (PRESENT(node_bias)) n_bias = node_bias
    last = UBOUND(csr%edges_of_vtx, 1) - 1

    IF (IAND(eflags, graph_io_undirected) /= 0) THEN
      edge_count = SIZE(csr%edges)/2
    ELSE
      edge_count = SIZE(csr%edges)
    END IF
    weight_fmt = 10
    node_count = num_nodes(csr)

    CALL assertion(SIZE(node_weights, 2) == node_count, filename, 412, &
         'number of node weights must match number of nodes')

    nnw = SIZE(node_weights, 1)
    CALL assertion(nnw > 0, filename, 416, &
           'number of node weights must be greater than zero if given')

    WRITE (unit=lun, fmt='(i0,3(1x,i0))') SIZE(csr%edges_of_vtx) - 1, &
         edge_count, weight_fmt, nnw
    DO i = ofs, last
      adj_count = csr%edges_of_vtx(i + 1) - csr%edges_of_vtx(i)
      IF (adj_count + nnw > 1) THEN
        WRITE(fmt, '(a,i0,a)') '(', adj_count + nnw - 1, '(i0,1x)i0)'
      ELSE
        fmt = '(i0)'
      END IF
      WRITE(unit=lun, fmt=fmt) &
           node_weights(:,i), &
           csr%edges(csr%edges_of_vtx(i):csr%edges_of_vtx(i+1)-1) + n_bias
    END DO
  END SUBROUTINE write_graph_csr_nwm

  SUBROUTINE write_graph_csr_nwm_ew(lun, csr, node_weights, edge_weights, &
       node_bias, flags)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER(i4), INTENT(in) :: node_weights(:,:), edge_weights(:)
    INTEGER, OPTIONAL, INTENT(in) :: node_bias, flags

    INTEGER :: i, last, lb, ub, num, ofs, n_bias, eflags, weight_fmt, nnw
    CHARACTER(len=11+1+9) :: fmt
    TYPE edge
      INTEGER :: node
      INTEGER(i4) :: weight
    END TYPE edge
    TYPE(edge), ALLOCATABLE :: edge_weight_vec(:)

    eflags = 0
    IF (PRESENT(flags)) eflags = flags

    ofs = LBOUND(csr%edges_of_vtx, 1)
    n_bias = -ofs + 1
    IF (PRESENT(node_bias)) n_bias = node_bias
    last = UBOUND(csr%edges_of_vtx, 1) - 1

    IF (IAND(eflags, graph_io_undirected) /= 0) THEN
      num = SIZE(csr%edges)/2
    ELSE
      num = SIZE(csr%edges)
    END IF
    weight_fmt = 11

    CALL assertion(SIZE(node_weights, 2) == num_nodes(csr), &
           filename, 465, &
           'number of node weights must match number of nodes')

    nnw = SIZE(node_weights, 1)
    CALL assertion(nnw > 0, filename, 469, &
           'number of node weights must be greater than zero if given')

    WRITE (unit=lun, fmt='(i0,3(1x,i0))') SIZE(csr%edges_of_vtx) - 1, &
         num, weight_fmt, nnw

    ALLOCATE(edge_weight_vec(10))
    DO i = ofs, last
      num = csr%edges_of_vtx(i + 1) - csr%edges_of_vtx(i)
      IF (2*num + nnw > 1) THEN
        WRITE(fmt, '(a,i0,a)') '(', 2*num - 1 + nnw, '(i0,1x)i0)'
      ELSE
        fmt = '(i0)'
      END IF
      IF (num > SIZE(edge_weight_vec)) THEN
        DEALLOCATE(edge_weight_vec)
        ALLOCATE(edge_weight_vec(num))
      END IF
      lb = csr%edges_of_vtx(i)
      ub = csr%edges_of_vtx(i + 1) - 1
      edge_weight_vec(1:num)%node = csr%edges(lb:ub) + n_bias
      edge_weight_vec(1:num)%weight = edge_weights(lb:ub)
      WRITE(unit=lun, fmt=fmt) node_weights(:,i), edge_weight_vec(1:num)
    END DO
    DEALLOCATE(edge_weight_vec)
  END SUBROUTINE write_graph_csr_nwm_ew

  SUBROUTINE write_graph_csr_nwmo_ewo(lun, csr, node_weights, edge_weights, &
       node_bias, flags)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_csr), INTENT(in) :: csr
    INTEGER(i4), OPTIONAL, INTENT(in) :: node_weights(:,:), edge_weights(:)
    INTEGER, OPTIONAL, INTENT(in) :: node_bias, flags

    IF (PRESENT(node_weights) .AND. PRESENT(edge_weights)) THEN
      CALL write_graph_csr_nwm_ew(lun, csr, node_weights, edge_weights, &
           node_bias=node_bias, flags=flags)
    ELSE IF (PRESENT(node_weights)) THEN
      CALL write_graph_csr_nwm(lun, csr, node_weights, &
           node_bias=node_bias, flags=flags)
    ELSE IF (PRESENT(edge_weights)) THEN
      CALL write_graph_csr_ew(lun, csr, edge_weights=edge_weights, &
           node_bias=node_bias, flags=flags)
    ELSE ! neither edge nor node weights given
      CALL write_graph_csr_plain(lun, csr, node_bias=node_bias, flags=flags)
    END IF

  END SUBROUTINE write_graph_csr_nwmo_ewo


  ELEMENTAL &

       FUNCTION graph_csr_equal_i4(a, b) RESULT(p)
    TYPE(graph_csr), INTENT(in) :: a, b
    LOGICAL :: p

    INTEGER(i4) :: i, nv, ne

    nv = num_nodes(a)
    p = nv == num_nodes(b)
    IF (.NOT. p) RETURN

    ne = num_edges(a)
    p = ne == num_edges(b)
    IF (.NOT. p) RETURN

    DO i = 1, nv
      p = is_permutation(a%edges(a%edges_of_vtx(i):a%edges_of_vtx(i+1)-1), &
           a%edges(b%edges_of_vtx(i):b%edges_of_vtx(i+1)-1))
      IF (.NOT. p) RETURN
    END DO
  END FUNCTION graph_csr_equal_i4

  SUBROUTINE read_graph(lun, csr, node_weights, edge_weights, ierror)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_csr), INTENT(out) :: csr
    INTEGER(i4), OPTIONAL, ALLOCATABLE, INTENT(out) :: node_weights(:,:), &
         edge_weights(:)
    INTEGER(i4), OPTIONAL, INTENT(out) :: ierror

    INTEGER(i4) :: inum_nodes, inum_edges, weight_fmt, nnw, count, dummy
    TYPE(fmt_elem), ALLOCATABLE :: var_fmt(:)
    INTEGER(i4), ALLOCATABLE :: node_desc(:)
    CHARACTER(len=1024) :: s
    INTEGER :: i, n, p, q, ew_stride, adj_len

    READ(lun, '(a)') s
    ALLOCATE(var_fmt(5))
    var_fmt(1)%argtype = arg_i4
    CALL get_address(inum_nodes, var_fmt(1)%addr)
    var_fmt(2)%argtype = arg_i4
    CALL get_address(inum_edges, var_fmt(2)%addr)
    var_fmt(3)%argtype = arg_i4
    CALL get_address(weight_fmt, var_fmt(3)%addr)
    var_fmt(4)%argtype = arg_i4
    CALL get_address(nnw, var_fmt(4)%addr)
    var_fmt(5)%argtype = arg_i4
    CALL get_address(dummy, var_fmt(5)%addr)
    CALL sscana(s, var_fmt, count, ierror)
    IF (PRESENT(ierror)) THEN; IF (ierror /= 0) RETURN; END IF
    SELECT CASE(count)
    CASE(0:1)
      IF (PRESENT(ierror)) THEN
        ierror = 1
        RETURN
      END IF
      CALL abort_ppm("insufficient graph description", filename, 575)
    CASE(2)
      weight_fmt = 0
      nnw = 0
    CASE(3)
      nnw = MERGE(1, 0, weight_fmt >= 10)
    CASE(5)
      IF (PRESENT(ierror)) THEN
        ierror = 1
        RETURN
      END IF
      CALL abort_ppm("erroneous graph description", filename, 586)
    END SELECT
    ew_stride = 1
    SELECT CASE(weight_fmt)
    CASE(0)
    CASE(1,11)
      ew_stride = 2
    CASE(10)
      IF (nnw < 1) THEN
        IF (PRESENT(ierror)) THEN
          ierror = 1
          RETURN
        END IF
        CALL abort_ppm("invalid weight format specification", filename, &
             600)
      END IF
    CASE default
      IF (PRESENT(ierror)) THEN
        ierror = 1
        RETURN
      END IF
      CALL abort_ppm("invalid weight format specification", filename, &
           608)
    END SELECT
    DEALLOCATE(var_fmt)
    n = nnw + 2 * MIN(inum_edges, inum_nodes)
    ALLOCATE(csr%edges_of_vtx(inum_nodes + 1), csr%edges(inum_edges * 2))
    ALLOCATE(var_fmt(n), node_desc(n))
    IF (PRESENT(node_weights) .AND. nnw > 0) &
         ALLOCATE(node_weights(nnw, inum_nodes))
    IF (PRESENT(edge_weights) .AND. ew_stride == 2) &
         ALLOCATE(edge_weights(inum_edges * 2))
    DO i = 1, n
      var_fmt(i)%argtype = arg_i4
      CALL get_address(node_desc(i), var_fmt(i)%addr)
    END DO
    i = 1
    csr%edges_of_vtx(1) = 1
    q = 1
    DO WHILE (i <= inum_nodes)
      READ(lun, '(a)') s
      ! skip comments
      IF (s(1:1) == '%') CYCLE
      CALL sscana(s, var_fmt, count)
      IF (count < nnw .OR. MOD(count - nnw, ew_stride) == 1) THEN
        IF (PRESENT(ierror)) THEN
          ierror = 1
          RETURN
        END IF
        CALL abort_ppm("malformed line in graph description", filename, &
             636)
      END IF
      p = q
      adj_len = (count - nnw)/ew_stride
      q = p + adj_len
      IF (PRESENT(node_weights) .AND. nnw > 0) &
           node_weights(1:nnw, i) = node_desc(1:nnw)
      csr%edges(p:q-1) = node_desc(nnw+1:count:ew_stride)
      IF (PRESENT(edge_weights) .AND. ew_stride == 2) &
           edge_weights(p:q-1) = node_desc(nnw+2:count:ew_stride)
      i = i + 1
      csr%edges_of_vtx(i) = q
    END DO
    DEALLOCATE(var_fmt, node_desc)
    IF (q < SIZE(csr%edges)) THEN
      ALLOCATE(node_desc(q))
      node_desc = csr%edges(1:q)
      DEALLOCATE(csr%edges)
      ALLOCATE(csr%edges(1:q))
      csr%edges = node_desc
      node_desc = edge_weights(1:q)
      DEALLOCATE(edge_weights)
      ALLOCATE(edge_weights(1:q))
      edge_weights = node_desc
      DEALLOCATE(node_desc)
    END IF
  END SUBROUTINE read_graph

  FUNCTION edge_index(graph, edge) RESULT(idx)
    TYPE(graph_csr), INTENT(in) :: graph
    INTEGER(i4), INTENT(in) :: edge(2)
    INTEGER(i4) :: idx

    INTEGER(i4) :: i, p, q

    p = graph%edges_of_vtx(edge(1))
    q = graph%edges_of_vtx(edge(1) + 1) - 1
    idx = UBOUND(graph%edges, 1) + 1
    DO i = p, q
      IF (graph%edges(i) == edge(2)) THEN
        idx = i
        EXIT
      END IF
    END DO
  END FUNCTION edge_index

  SUBROUTINE assign_edge_weight(graph, edge, weights, weight)
    TYPE(graph_csr), INTENT(inout) :: graph
    INTEGER(i4), INTENT(in) :: weight, edge(2)
    INTEGER(i4), INTENT(inout) :: weights(:)

    INTEGER(i4) :: i

    CALL assertion(edge(1) >= LBOUND(graph%edges_of_vtx, 1) &
         .AND. edge(1) < UBOUND(graph%edges_of_vtx, 1), &
         filename, 691, 'edge violates node constraints')
    i = edge_index(graph, edge)
    CALL assertion(i <= UBOUND(graph%edges, 1), filename, 693, &
         'invalid edge specified')
    weights(i) = weight
  END SUBROUTINE assign_edge_weight

  SUBROUTINE assign_symmetric_edge_weight(graph, edge, weights, weight)
    TYPE(graph_csr), INTENT(inout) :: graph
    INTEGER(i4), INTENT(in) :: weight, edge(2)
    INTEGER(i4), INTENT(inout) :: weights(:)
    INTEGER(i4) :: iedge(2)
    CALL assign_edge_weight(graph, edge, weights, weight)
    iedge(1) = edge(2)
    iedge(2) = edge(1)
    CALL assign_edge_weight(graph, iedge, weights, weight)
  END SUBROUTINE assign_symmetric_edge_weight

  FUNCTION graph_is_symmetric_csr(graph) RESULT(symmetry)
    TYPE(graph_csr), INTENT(in) :: graph
    LOGICAL :: symmetry

    INTEGER(i4) :: edge_inv, i, j, n, p, q, edge(2)

    symmetry = .TRUE.
    n = num_nodes(graph)
    edge_loop: DO i = 1_i4, n
      p = graph%edges_of_vtx(i)
      q = graph%edges_of_vtx(i + 1) - 1
      DO j = p, q
        IF (graph%edges(j) > i) THEN
          edge(1) = i
          edge(2) = graph%edges(j)
          edge_inv = edge_index(graph, edge)
          symmetry = edge_inv <= UBOUND(graph%edges, 1)
          IF (.NOT. symmetry) EXIT edge_loop
        END IF
      END DO
    END DO edge_loop
  END FUNCTION graph_is_symmetric_csr

  FUNCTION graph_is_symmetric_csr_ew(graph, edge_weights) RESULT(symmetry)
    TYPE(graph_csr), INTENT(in) :: graph
    INTEGER(i4), INTENT(in) :: edge_weights(:)
    LOGICAL :: symmetry

    INTEGER(i4) :: edge_inv, i, j, n, p, q, edge(2)

    symmetry = .TRUE.
    n = num_nodes(graph)
    CALL assertion(SIZE(edge_weights) == SIZE(graph%edges), &
         filename, 742, 'non-matching edge weights argument')
    edge_loop: DO i = 1_i4, n
      p = graph%edges_of_vtx(i)
      q = graph%edges_of_vtx(i + 1) - 1
      DO j = p, q
        IF (graph%edges(j) > i) THEN
          edge(1) = graph%edges(j)
          edge(2) = i
          edge_inv = edge_index(graph, edge)
          symmetry = edge_inv <= UBOUND(graph%edges, 1)
          IF (.NOT. symmetry) EXIT edge_loop
          symmetry = edge_weights(j) == edge_weights(edge_inv)
          IF (.NOT. symmetry) EXIT edge_loop
        END IF
      END DO
    END DO edge_loop
  END FUNCTION graph_is_symmetric_csr_ew

END MODULE ppm_graph_csr
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
