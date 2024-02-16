!> @file graph_build.f90 --- construct graph data structure from rectilinear
!!                           specification
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords: rect graph construction
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
PROGRAM graph_build
  USE ppm_base, ONLY: ppm_default_comm
  USE ppm_extents, ONLY: iinterval
  USE ppm_f90_io_lun, ONLY: setup_lun_table, next_free_unit
  USE ppm_extents, ONLY: extent, extent_shape, empty_extent, iinterval
  USE ppm_graph_csr, ONLY: graph_csr, write_graph, graph_io_undirected
  USE ppm_random, ONLY: initialize_irand, irandr, lrand
  USE ppm_std_type_kinds, ONLY: i4
  USE random_data, ONLY: random_rect_graph
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE

  TYPE(graph_csr) :: graph
  INTEGER :: graph_out_lun, i, m, rand_seed, tid, param_in, ierror
  CHARACTER(len=255) :: fname, tempfmt
  CHARACTER(len=*), PARAMETER :: config_file="graph_build.conf"
  INTEGER, PARAMETER :: max_rank=3
  TYPE(extent) :: rect(max_rank)
  INTEGER :: num_node_weights
  INTEGER(i4), ALLOCATABLE :: node_weight(:,:), edge_weight(:)
  LOGICAL :: use_edge_weights
  NAMELIST /rect_dim/ m, rect
  NAMELIST /weight_conf/ num_node_weights, use_edge_weights
!$omp parallel shared(graph, m, rect, num_node_weights, use_edge_weights) &
!$omp private(param_in, tid, graph_out_lun, rand_seed, &
!$omp tempfmt, fname)
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed
!$omp master
  CALL setup_lun_table
  param_in = next_free_unit()
  OPEN(unit=param_in, file=config_file, iostat=ierror, status='old', &
       access='sequential', form='formatted', action='read')
  num_node_weights = irandr(iinterval(0, 4))
  use_edge_weights = lrand()
  IF (ierror /= 0) THEN
    PRINT '(a)', 'using random dimensions instead of ' // config_file &
         // ' graph size determination'
    m = irandr(iinterval(1, max_rank))
    DO i = 1, m
      rect(i) = extent(1, irandr(iinterval(1, 10)))
    END DO
  ELSE
    m = 0
    rect = empty_extent
    READ(unit=param_in, nml=rect_dim)
    ! errors on weight constitution are non-fatal, because this namelist is
    ! optional
    READ(unit=param_in, nml=weight_conf, iostat=ierror)
  END IF
  graph_out_lun = next_free_unit()
  IF (m > 1) THEN
    WRITE (tempfmt, '(a,i0,a)') '(a,i0,', m-1, '("x",i0),a)'
    WRITE (fname, tempfmt) 'rect_graph_', extent_shape(rect(1:m)), &
         '.metis'
  ELSE
    WRITE (fname, '(a,i0,a)') 'rect_graph_', extent_shape(rect(1)), '.metis'
  END IF
  OPEN(graph_out_lun, file=fname, access='sequential', &
       form='formatted', action='write')
!$omp end master
!$omp barrier
!$omp master
  CALL random_rect_graph(graph, node_weight, edge_weight, rect(1:m), &
       num_node_weights, use_edge_weights)
  IF (use_edge_weights) THEN
    IF (num_node_weights > 0) THEN
      CALL write_graph(graph_out_lun, graph, node_weights=node_weight, &
           edge_weights=edge_weight, flags=graph_io_undirected)
    ELSE
      CALL write_graph(graph_out_lun, graph, edge_weights=edge_weight, &
           flags=graph_io_undirected)
    END IF
  ELSE
    IF (num_node_weights > 0) THEN
      CALL write_graph(graph_out_lun, graph, node_weights=node_weight, &
           flags=graph_io_undirected)
    ELSE
      CALL write_graph(graph_out_lun, graph, flags=graph_io_undirected)
    END IF
  END IF
!$omp end master
!$omp end parallel
END PROGRAM graph_build
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! license-markup: "doxygen"
! End:
