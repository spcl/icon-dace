!> @file ppm_graph_partition_mpi.f90
!! @brief generic graph partition method interface
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords: graph partitioning
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
!> This is currently only a convenient wrapper of ParMeTis, other
!! heuristics are to follow later.
#include "fc_feature_defs.inc"
MODULE ppm_graph_partition_mpi
  ! we need iso_c_binding without only because configure will
  ! substitute ppm_metis_int according to the actual library type
  USE iso_c_binding
  USE ppm_base, ONLY: abort_ppm
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
#include <ppm.inc>
#ifdef HAVE_PARMETIS_V3
  INTEGER, PARAMETER :: ppm_metis_int = c_int
#else
  INTEGER, PARAMETER :: ppm_metis_int = ppm_metis_idx
#endif
  INTERFACE
    SUBROUTINE parmetis_v3_partkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
         wgtflag, numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, &
         part, comm) BIND(C)
      IMPORT :: ppm_metis_idx, ppm_metis_int, ppm_metis_real, c_ptr, &
           ppm_mpi_fint_fc_kind
      INTEGER(ppm_metis_idx), INTENT(in) :: vtxdist(*), xadj(*), adjncy(*)
      TYPE(c_ptr), VALUE, INTENT(in) :: vwgt, adjwgt
      INTEGER(ppm_metis_int), INTENT(in) :: wgtflag, numflag, ncon, nparts, options(*)
      TYPE(c_ptr), VALUE, INTENT(in) :: tpwgts
      REAL(ppm_metis_real), INTENT(in) :: ubvec(ncon)
      INTEGER(ppm_metis_int), INTENT(out) :: edgecut
      INTEGER(ppm_metis_idx), INTENT(out) :: part(*)
      INTEGER(ppm_mpi_fint_fc_kind), INTENT(in) :: comm
    END SUBROUTINE parmetis_v3_partkway
  END INTERFACE
  PUBLIC :: graph_partition_parmetis

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_graph_partition_mpi.f90'
CONTAINS
  !> call parmetis to partition graph in distributed CSR format
  !!
  !! @param num_vertices number of vertices
  !! @param edge_list_lens cumulative number of edges per vertex
  !! @param edge_lists concatenated list of vertices connected to each vertex
  !! @param partition_out for each vertex passed write partition assignment here
  !! @param comm optional communicator object identifying participating
  !! processes, defaults to ppm_default_comm
  !! @param num_partitions optional number of desired partitions,
  !! defaults to size of comm
  !! @param balance specify desired distribution of loads
  !! @param num_vertex_weights number of weighting factors per vertex
  !! @param vertex_weights
  SUBROUTINE graph_partition_parmetis(num_vertices, edge_list_lens, &
       edge_lists, partition_out, comm, num_partitions, &
       balance, num_vertex_weights, vertex_weights, edge_weights)
    INTEGER(ppm_metis_idx), INTENT(in) :: num_vertices
    INTEGER(ppm_metis_idx), INTENT(in) :: edge_list_lens(*)
    INTEGER(ppm_metis_idx), INTENT(in) :: edge_lists(*)
    INTEGER(ppm_metis_idx), INTENT(out) :: partition_out(*)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER(ppm_metis_int), OPTIONAL, INTENT(in) :: num_partitions
    REAL(ppm_metis_real), OPTIONAL, TARGET, INTENT(in) :: balance(1, *)
    INTEGER(ppm_metis_int), OPTIONAL, INTENT(in) :: num_vertex_weights
    INTEGER(ppm_metis_idx), OPTIONAL, TARGET, INTENT(in) :: vertex_weights(*)
    INTEGER(ppm_metis_idx), OPTIONAL, TARGET, INTENT(in) :: edge_weights(*)
    INTEGER(ppm_metis_int) :: wgtflag
    INTEGER :: part_comm, comm_size, comm_rank, ierror, i, ierror_
    INTEGER(ppm_metis_idx), ALLOCATABLE :: vtxdist(:)
    INTEGER(ppm_metis_int) :: metis_options(0:2), edge_cut, ncon, num_parts, &
         accum
    INTEGER :: msg_len
    CHARACTER(len=mpi_max_error_string) :: msg
    TYPE(c_ptr) :: vwgt, adjwgt
    TYPE(c_ptr) :: tpwgts
#ifndef HAVE_PARMETIS_V3
    REAL(ppm_metis_real), ALLOCATABLE, TARGET :: tpwgts_balance(:, :)
#endif
    INTEGER :: ppm_metis_idx_mpidt = mpi_datatype_null

    IF (PRESENT(comm)) THEN; part_comm = comm; ELSE; part_comm = mpi_comm_world
    END IF
    CALL mpi_comm_size(part_comm, comm_size, ierror)
    IF (ierror /= MPI_SUCCESS) THEN
      CALL mpi_error_string(ierror, msg, msg_len, ierror_)
      CALL abort_ppm(msg(1:msg_len), filename, __LINE__, comm)
    END IF
    CALL mpi_comm_rank(part_comm, comm_rank, ierror)
    IF (ierror /= MPI_SUCCESS) THEN
      CALL mpi_error_string(ierror, msg, msg_len, ierror_)
      CALL abort_ppm(msg(1:msg_len), filename, __LINE__, comm)
    END IF

    ! build table of node distribution
    ALLOCATE(vtxdist(0:comm_size))

    IF (ppm_metis_idx_mpidt == mpi_datatype_null) THEN
      CALL mpi_type_create_f90_integer(RANGE(num_vertices), &
           &                           ppm_metis_idx_mpidt, ierror)
      IF (ierror /= MPI_SUCCESS) THEN
        CALL mpi_error_string(ierror, msg, msg_len, ierror_)
        CALL abort_ppm(msg(1:msg_len), filename, __LINE__, comm)
      END IF
    END IF
    CALL mpi_allgather(num_vertices, 1, ppm_metis_idx_mpidt, &
         vtxdist(1), 1, ppm_metis_idx_mpidt, part_comm, ierror)
    IF (ierror /= MPI_SUCCESS) then
      CALL mpi_error_string(ierror, msg, msg_len, ierror_)
      CALL abort_ppm(msg(1:msg_len), filename, __LINE__, comm)
    END IF

    ! compute partial sums of vertices over ranks
    accum = 1_ppm_metis_idx
    DO i = 0, comm_size-1
      vtxdist(i) = accum
      accum = accum + vtxdist(i+1)
    END DO
    vtxdist(comm_size) = accum

    wgtflag = 0_ppm_metis_int
    IF (PRESENT(vertex_weights)) wgtflag = 2_ppm_metis_int
    IF (PRESENT(edge_weights)) wgtflag = IOR(wgtflag, 1_ppm_metis_int)
    IF (PRESENT(num_vertex_weights)) THEN
      ncon = num_vertex_weights
    ELSE
      ncon = 1
    END IF
    IF (PRESENT(num_partitions)) THEN
      num_parts = num_partitions
    ELSE
      num_parts = comm_size
    END IF
    metis_options(0) = 0
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
    IF (PRESENT(balance)) THEN
      tpwgts = c_LOC(balance(1,1))
    ELSE
#ifdef HAVE_PARMETIS_V3
      tpwgts = c_null_ptr
#else
      ALLOCATE(tpwgts_balance(ncon, num_parts))
      tpwgts_balance = 1.0_ppm_metis_real / REAL(num_parts, ppm_metis_real)
      tpwgts = c_loc(tpwgts_balance(1, 1))
#endif
    END IF
    CALL parmetis_v3_partkway(vtxdist, edge_list_lens, edge_lists, &
         vwgt, adjwgt, wgtflag, 1_ppm_metis_int, ncon, num_parts, tpwgts, &
         (/ 1.05_ppm_metis_real /), metis_options, edge_cut, &
         partition_out, part_comm)
  END SUBROUTINE graph_partition_parmetis
END MODULE ppm_graph_partition_mpi
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
