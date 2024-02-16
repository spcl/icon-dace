!> @file random_data.f90
!! @brief generate data structures randomly for test purposes
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
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
#include "fc_feature_defs.inc"
MODULE random_data
  USE ppm_extents, ONLY: extent, iinterval
  USE ppm_graph_csr, ONLY: graph_csr, assign_symmetric_edge_weight, &
       num_nodes, num_edges, build_graph
  USE ppm_random, ONLY: a_randr, a_randp, irandr, lrand
  USE ppm_std_type_kinds, ONLY: dp, i4
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: random_rect_graph
  PUBLIC :: random_rectilinear
CONTAINS
  !> construct rectilinear within given size limit
  SUBROUTINE random_rectilinear(rect, size_lim)
    TYPE(extent), INTENT(out) :: rect(:)
    INTEGER, INTENT(in) :: size_lim

    REAL(dp) :: size_dist(SIZE(rect))
    REAL(dp), PARAMETER :: p = 0.0_dp
    INTEGER :: i, n
    INTEGER(i4), PARAMETER :: j = 0_i4

    n = SIZE(rect)
    ! compute relative distribution of sizes
    CALL a_randp(size_dist)
    size_dist = size_dist + EPSILON(p)
    ! scale to produce full size
    size_dist = size_dist/PRODUCT(size_dist)**(1.0_dp/REAL(n, dp)) &
         * REAL(size_lim, dp)**(1.0_dp/REAL(n, dp))
    rect%size = MAX(1_i4, INT(size_dist, i4))
    DO i = 1, n
      rect(i)%first = irandr(iinterval(-HUGE(j), HUGE(j) - rect(i)%size))
    END DO
  END SUBROUTINE random_rectilinear

  SUBROUTINE random_rect_graph(graph, node_weights, edge_weights, rect, &
       num_node_weights, with_edge_weights)
    TYPE(graph_csr), INTENT(out) :: graph
    INTEGER(i4), ALLOCATABLE, OPTIONAL, INTENT(out) :: node_weights(:,:), &
         edge_weights(:)
    TYPE(extent), OPTIONAL, INTENT(in) :: rect(:)
    INTEGER, OPTIONAL, INTENT(in) :: num_node_weights
    LOGICAL, OPTIONAL, INTENT(in) :: with_edge_weights

    INTEGER, PARAMETER :: max_rank=3
    INTEGER(i4) :: i, j, m, n, p, q, nnw, edge(2)
    TYPE(extent), ALLOCATABLE :: rect_(:)
    LOGICAL :: use_edge_weights
    TYPE(iinterval), PARAMETER :: weight_range = iinterval(1, 1000)

    IF (PRESENT(rect)) THEN
      m = SIZE(rect)
      ALLOCATE(rect_(m))
      rect_ = rect
    ELSE
      m = irandr(iinterval(1, max_rank))
      ALLOCATE(rect_(m))
      DO i = 1, m
        rect_(i) = extent(1, irandr(iinterval(1, 10)))
      END DO
    END IF
    CALL build_graph(graph, rect_)
    n = num_nodes(graph)
    IF (PRESENT(node_weights)) THEN
      IF (PRESENT(num_node_weights)) THEN
        nnw = num_node_weights
      ELSE
        nnw = irandr(iinterval(0, 5))
      END IF
      IF (nnw > 0) THEN
        ALLOCATE(node_weights(nnw, n))
        CALL a_randr(node_weights, weight_range)
      END IF
    ELSE
      nnw = 0
    END IF
    IF (PRESENT(edge_weights)) THEN
      IF (PRESENT(with_edge_weights)) THEN
        use_edge_weights = with_edge_weights
      ELSE
        use_edge_weights = lrand()
      END IF
    ELSE
      use_edge_weights = .FALSE.
    END IF
    IF (use_edge_weights) THEN
      ALLOCATE(edge_weights(num_edges(graph)))
      DO i = 1, n
        p = graph%edges_of_vtx(i)
        q = graph%edges_of_vtx(i + 1) - 1
        edge(1) = i
        DO j = p, q
          IF (graph%edges(j) >= i) THEN
            edge(2) = graph%edges(j)
            CALL assign_symmetric_edge_weight(graph, &
                 edge, edge_weights, irandr(weight_range))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE random_rect_graph

END MODULE random_data
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
