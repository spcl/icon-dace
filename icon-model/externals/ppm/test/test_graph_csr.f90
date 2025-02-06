!>
!! @file test_graph_csr.f90
!! @brief black box test for csr data structure
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: CSR compressed sparse row
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
PROGRAM test_graph_csr
  USE ppm_base, ONLY: ppm_default_comm, abort_ppm
  USE ppm_extents, ONLY: iinterval
  USE ppm_graph_csr, ONLY: graph_csr, build_graph, build_graph_mt
  USE ppm_random, ONLY: initialize_irand, irand, irandr, lrand
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nodes = 4000
  TYPE(graph_csr) :: graph
  LOGICAL, ALLOCATABLE :: adj(:,:)
#if defined __INTEL_COMPILER && __INTEL_COMPILER >= 1400 && __INTEL_COMPILER <= 1600 && defined _OPENMP && ! defined __OPTIMIZE__
#define NEED_LB_WORKAROUND
  TARGET :: adj
  LOGICAL, POINTER :: adj_p(:, :)
#endif
  INTEGER :: e, i, j, m, n, ofs, rand_seed
  INTEGER :: tid
  TYPE(iinterval) :: bnds(2)
  LOGICAL :: p
  CHARACTER(len=*), PARAMETER :: filename = 'test_graph_csr.f90'

!$omp parallel private(e, i, j, m, p, rand_seed, tid, bnds) &
#ifdef NEED_LB_WORKAROUND
!$omp & private(adj_p) &
#endif
!$omp & shared(adj)
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
!$omp critical
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed
!$omp end critical

!$omp master
#ifdef __G95__
  ! g95 cannot handle arbitrary lower array bounds
  ofs = 1
#else
  ofs = HUGE(ofs)
  DO WHILE (ofs > HUGE(ofs) - max_nodes)
    ofs = irand()
  END DO
#endif
  n = irandr(iinterval(0, max_nodes-1)) + ofs
  ALLOCATE(adj(ofs:n, ofs:n))
!$omp end master
!$omp barrier
  ! first test undirected graph
!$omp do firstprivate(n)
  DO j = ofs, n
    adj(j, j) = lrand()
    DO i = j + 1, n
      p = lrand()
      adj(i, j) = p
      adj(j, i) = p
    END DO
  END DO
!$omp end do nowait
#ifdef NEED_LB_WORKAROUND
  ! this work-around is needed because ifort 14.x-16.x with flags -openmp -O0 -g
  ! might loose the lbounds of an array when calling build_graph
  adj_p => adj
  CALL build_graph(graph, adj_p, assert_undirected = .TRUE., &
       node_offset=ofs)
#else
  CALL build_graph(graph, adj, assert_undirected = .TRUE., &
       node_offset=ofs)
#endif

  p = .FALSE.
!$omp do firstprivate(n)
  DO j = ofs, n
    m = graph%edges_of_vtx(j + 1) - 1
    DO e = graph%edges_of_vtx(j), m
      p = p .OR. .NOT. adj(graph%edges(e), j)
    END DO
  END DO
!$omp end do
  IF (p) THEN
    CALL abort_ppm('graph CSR construction failed', filename, __LINE__)
  END IF

  bnds = iinterval(ofs,n)
  CALL build_graph_mt(graph, bnds)
!$omp end parallel

END PROGRAM test_graph_csr
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
