!>
!! @file ppm_graph_alist.f90
!! @brief describe graph in terms of adjacency lists
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: graph adjacency list
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
!

MODULE ppm_graph_alist
  USE ppm_base, ONLY: assertion, abort_ppm
  IMPLICIT NONE
  PRIVATE
  TYPE adjl
    INTEGER, ALLOCATABLE :: node(:)
  END TYPE adjl
  TYPE graph_adjl
    TYPE(adjl), ALLOCATABLE :: adj(:)
    INTEGER :: edge_count
  END TYPE graph_adjl
  CHARACTER(len=*), PARAMETER :: filename = 'ppm_graph_alist.f90'
CONTAINS
  SUBROUTINE build_graph_alist_from_adj_matrix(list, adj, &
       assert_undirected, node_offset)
    TYPE(graph_adjl), INTENT(inout) :: list
    LOGICAL, INTENT(in) :: adj(:, :)
    LOGICAL, OPTIONAL, INTENT(in) :: assert_undirected
    INTEGER, OPTIONAL, INTENT(in) :: node_offset

    INTEGER :: i, j, k, m, nedges, ofs, edge_count

    m = SIZE(adj, 1)
    ofs = 1
    IF (PRESENT(node_offset)) ofs = node_offset
    CALL assertion(m == SIZE(adj, 2), filename, 68, &
         'adjacency matrix must be rectangular')
!$omp single
    IF (ALLOCATED(list%adj)) DEALLOCATE(list%adj)
    ALLOCATE(list%adj(ofs:ofs+m))
    list%edge_count = 0
!$omp end single
    edge_count = 0
!$omp do
    DO j = ofs, ofs+m
      k = 0
      nedges = COUNT(adj(:, j))
      ALLOCATE(list%adj(j)%node(nedges))
      DO i = 1, m
        IF (adj(i, j - ofs + 1)) THEN
          k = k + 1
          list%adj(j)%node(k) = i + ofs - 1
          IF (k == nedges) EXIT
        END IF
      END DO
      edge_count = edge_count + nedges
    END DO
!$omp end do
!$omp atomic
    list%edge_count = list%edge_count + edge_count

    IF (PRESENT(assert_undirected)) THEN
!$omp do
      DO j = 1, m
        DO i = j + 1, m
          IF (adj(i, j) .NEQV. adj(j, i)) THEN
            CALL abort_ppm('input matrix must be symmetric', filename, &
                 100)
          END IF
        END DO
      END DO
!$omp end do
    END IF
  END SUBROUTINE build_graph_alist_from_adj_matrix

  SUBROUTINE print_graph_alist(lun, list)
    INTEGER, INTENT(in) :: lun
    TYPE(graph_adjl), INTENT(in) :: list

    INTEGER :: i, ofs, last
    CHARACTER(len=9+4) :: fmt

    ofs = LBOUND(list%adj, 1)
    last = UBOUND(list%adj, 1)
    WRITE (unit=lun, fmt='(i0,1x,i0)') SIZE(list%adj), list%edge_count
    DO i = ofs, last
      WRITE(fmt, '(a,i0,a)') '(', SIZE(list%adj(i)%node), 'i0)'
      WRITE(unit=lun, fmt=fmt) list%adj(i)%node - ofs + 1
    END DO
  END SUBROUTINE print_graph_alist
END MODULE ppm_graph_alist
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
