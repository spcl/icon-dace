!> @file graph_partition.f90
!! @brief short description
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! ENDDOXYGENPART
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
PROGRAM graph_partition
  USE iso_c_binding, ONLY: c_float
  USE ppm_base, ONLY: ppm_default_comm, assertion
  USE ppm_extents, ONLY: iinterval
  USE ppm_f90_io_lun, ONLY: setup_lun_table, next_free_unit
  USE ppm_graph_csr, ONLY: graph_csr, read_graph, graph_is_symmetric
  USE ppm_graph_partition_serial, ONLY: graph_partition_metis
  USE ppm_posix, ONLY: path_max
  USE ppm_random, ONLY: initialize_irand, irandr
  USE ppm_set_partition_base, ONLY: partition_assignment
  USE ppm_std_type_kinds, ONLY: i4
  USE random_data, ONLY: random_rect_graph
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE
  INTEGER :: tid, ierror, param_in, graph_in
  INTEGER(i4) :: rand_seed
  CHARACTER(len=*), PARAMETER :: config_file="graph_partition.conf"
  CHARACTER(len=path_max) :: graph_file
  TYPE(graph_csr) :: graph
  INTEGER(i4), ALLOCATABLE :: node_weights(:,:), edge_weights(:)
  TYPE(partition_assignment) :: partition
  INTEGER(i4) :: num_parts
  INTEGER :: nnw
  REAL(c_float), ALLOCATABLE :: imbalance_tolerance(:)
  CHARACTER(len=*), PARAMETER :: filename = 'graph_partition.f90'
  NAMELIST /graph_conf/ graph_file, num_parts
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  CALL initialize_irand(ppm_default_comm, 0_i4, rand_seed)
  WRITE (0, '(2(a,i0))') 'thread id=', tid, ', random seed=', rand_seed
  CALL setup_lun_table

  num_parts = -1
  graph_file = ' '

  param_in = next_free_unit()
  OPEN(unit=param_in, file=config_file, iostat=ierror, status='old', &
       action='read')
  IF (ierror /= 0) THEN
    PRINT '(a)', 'using random graph instead of ' // config_file &
         // ' graph file determination'
    CALL random_rect_graph(graph, node_weights, edge_weights)
  ELSE
    ! errors on weight constitution are non-fatal, because this namelist is
    ! optional
    READ(unit=param_in, nml=graph_conf, iostat=ierror)
    IF (ierror /= 0) STOP 1
    CLOSE(param_in)
    graph_in = next_free_unit()
    OPEN(unit=graph_in, file=graph_file, status='old')
    CALL read_graph(graph_in, graph, node_weights, edge_weights)
  END IF

  IF (num_parts < 1) num_parts = irandr(iinterval(1, 10))

  IF (ALLOCATED(edge_weights)) THEN
    CALL assertion(graph_is_symmetric(graph, edge_weights), &
         filename, __LINE__, 'graph is asymmetric')
  ELSE
    CALL assertion(graph_is_symmetric(graph), &
         filename, __LINE__, 'graph is asymmetric')
  END IF

  IF (ALLOCATED(node_weights)) THEN
    nnw = SIZE(node_weights, 1)
    IF (nnw > 1) THEN
      ALLOCATE(imbalance_tolerance(nnw))
      imbalance_tolerance = 1.05
      IF (ALLOCATED(edge_weights)) THEN
        CALL graph_partition_metis(partition, graph, &
             num_parts, imbalance_tolerance, node_weights, edge_weights)
      ELSE
        CALL graph_partition_metis(partition, graph, &
             num_parts, imbalance_tolerance, node_weights)
      END IF
    ELSE
      IF (ALLOCATED(edge_weights)) THEN
        CALL graph_partition_metis(partition, graph, &
             num_parts, 1_i4, vertex_weights = node_weights(1, :), &
             edge_weights = edge_weights)
      ELSE
        CALL graph_partition_metis(partition, graph, &
             num_parts, 1_i4, vertex_weights = node_weights(1, :))
      END IF
    END IF
  ELSE
    IF (ALLOCATED(edge_weights)) THEN
      CALL graph_partition_metis(partition, graph, &
           num_parts, 1_i4, edge_weights = edge_weights)
    ELSE
      CALL graph_partition_metis(partition, graph, num_parts)
    END IF
  END IF
  PRINT '(i0)', partition%assigned
END PROGRAM graph_partition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
