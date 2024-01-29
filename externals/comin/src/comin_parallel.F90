!> @file comin_parallel.F90
!! @brief Definitions for parallel communication.
!
!  @authors 09/2022 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_parallel

  USE mpi
  USE comin_state,          ONLY: state
  USE comin_errhandler,     ONLY: comin_plugin_finish
  USE iso_c_binding,        ONLY: c_int
  USE comin_parallel_types, ONLY: t_comin_parallel_info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_parallel_get_host_mpi_comm
  PUBLIC :: comin_parallel_get_plugin_mpi_comm
  PUBLIC :: comin_parallel_get_host_mpi_rank
  PUBLIC :: comin_parallel_free_mpi_comms
  PUBLIC :: comin_parallel_mpi_handshake
  PUBLIC :: comin_parallel_handle_mpi_errcode


#include "comin_global.inc"

  ! max. character string length (e.g. MPI group name)
  INTEGER, PARAMETER :: MAX_GRPNAMELEN = 256

  INTEGER, PARAMETER :: MAX_GROUPS = 64

  ! level for debug output verbosity (0: quiet)
  INTEGER, PARAMETER :: debug_level = 0


CONTAINS

  !> Procedure for the communicator splitting ("MPI handshake") that
  !> has been harmonized with the respective algorithm of the YAC
  !> coupler software
  !!
  !! @ingroup host_interface
  SUBROUTINE comin_parallel_mpi_handshake(comm, group_names, component_name)
    USE mo_mpi_handshake, ONLY: mpi_handshake

    INTEGER, INTENT(IN)                       :: comm
    CHARACTER(LEN=MAX_GRPNAMELEN), INTENT(IN) :: group_names(:)
    CHARACTER(LEN=*), INTENT(IN)              :: component_name
    !
    INTEGER :: num_plugin_comms
    CHARACTER(LEN=MAX_GRPNAMELEN), ALLOCATABLE :: all_group_names(:)
    INTEGER, ALLOCATABLE                       :: all_comms(:)

    num_plugin_comms = COUNT(group_names /= "")
    ALLOCATE(all_group_names(num_plugin_comms+1))
    IF(num_plugin_comms > 0) &
         all_group_names(1:num_plugin_comms) = PACK(group_names, group_names /= "")
    all_group_names(num_plugin_comms+1) = "comin_host_"//component_name

    ALLOCATE(all_comms(num_plugin_comms+1))
    CALL mpi_handshake(comm, all_group_names, all_comms)

    ALLOCATE(state%parallel_info%component_name_bilateral(num_plugin_comms))
    state%parallel_info%component_name_bilateral = all_group_names(1:num_plugin_comms)
    state%parallel_info%mpi_comm_bilateral = all_comms(1:num_plugin_comms)
    state%parallel_info%host_comm = all_comms(num_plugin_comms+1)

    DEALLOCATE(all_group_names)
    DEALLOCATE(all_comms)
  END SUBROUTINE comin_parallel_mpi_handshake


  SUBROUTINE comin_parallel_free_mpi_comms()
    INTEGER :: i, ierr, nmpi_comm

    IF (.NOT. ALLOCATED(state%parallel_info%mpi_comm_bilateral))  RETURN
    nmpi_comm = SIZE(state%parallel_info%mpi_comm_bilateral)
    DO i=1,nmpi_comm
       CALL MPI_Comm_free(state%parallel_info%mpi_comm_bilateral(i), ierr)
       CALL comin_parallel_handle_mpi_errcode(ierr)
       state%parallel_info%mpi_comm_bilateral(i) = MPI_COMM_NULL
    END DO
    DEALLOCATE(state%parallel_info%mpi_comm_bilateral)
    DEALLOCATE(state%parallel_info%component_name_bilateral)
  END SUBROUTINE comin_parallel_free_mpi_comms


  !> Returns the communicator with all ICON processes.
  !! @ingroup common
  !!
  INTEGER FUNCTION comin_parallel_get_host_mpi_comm() RESULT(mpi_comm)
    mpi_comm = state%parallel_info%host_comm
  END FUNCTION comin_parallel_get_host_mpi_comm

  ! Wrapper function that converts the communicator to C_INT
  INTEGER(KIND=C_INT) FUNCTION comin_parallel_get_host_mpi_comm_c() &
       RESULT(mpi_comm) BIND(C, name="comin_parallel_get_host_mpi_comm")
    mpi_comm = comin_parallel_get_host_mpi_comm()
  END FUNCTION comin_parallel_get_host_mpi_comm_c


  !> Called within a plugin's callback function: get MPI communicator
  !> which contains all MPI tasks of the host model together with the
  !> plugin's external MPI partners (if any).
  !! @ingroup plugin_interface
  !!
  INTEGER FUNCTION comin_parallel_get_plugin_mpi_comm() RESULT(mpi_comm)
    mpi_comm = comin_parallel_find_mpi_comm(state%parallel_info, state%current_plugin%comm)
    IF(mpi_comm == MPI_COMM_NULL) THEN
       CALL comin_plugin_finish("comin_parallel_get_plugin_mpi_comm", &
            &                   "Error: To use the plugin_comm, 'comm' must be set in the plugin namelist")
    END IF
  END FUNCTION comin_parallel_get_plugin_mpi_comm

  ! Wrapper function that converts the communicator to C_INT
  INTEGER(KIND=C_INT) FUNCTION comin_parallel_get_plugin_mpi_comm_c() &
       RESULT(mpi_comm) BIND(C, name="comin_parallel_get_plugin_mpi_comm")
    mpi_comm = comin_parallel_get_plugin_mpi_comm()
  END FUNCTION comin_parallel_get_plugin_mpi_comm_c


  !> Called within a plugin's callback function: get MPI rank with
  !> respect to the "host" MPI communicator.
  !! @ingroup plugin_interface
  !!
  INTEGER(KIND=C_INT) FUNCTION comin_parallel_get_host_mpi_rank() RESULT(mpi_rank) BIND(C)
    INTEGER :: mpi_comm, ierr, rank_f
    mpi_comm = comin_parallel_get_host_mpi_comm()
    CALL MPI_COMM_RANK(mpi_comm, rank_f, ierr)
    CALL comin_parallel_handle_mpi_errcode(ierr)
    mpi_rank = rank_f
  END FUNCTION comin_parallel_get_host_mpi_rank


  !> Auxiliary routine for finding an MPI intra-communicator, formed
  !  by a PE of each group plus a "root group".
  !
  INTEGER FUNCTION comin_parallel_find_mpi_comm(parallel_info, component_name)
    TYPE(t_comin_parallel_info), INTENT(IN) :: parallel_info
    CHARACTER(LEN=*), INTENT(IN) :: component_name
    !
    INTEGER :: i, nmpi_comm

    comin_parallel_find_mpi_comm = MPI_COMM_NULL
    IF (.NOT. ALLOCATED(parallel_info%mpi_comm_bilateral))  RETURN
    nmpi_comm = SIZE(parallel_info%mpi_comm_bilateral)
    LOOP : DO i=1,nmpi_comm
      IF (TRIM(parallel_info%component_name_bilateral(i)) == component_name) THEN
        comin_parallel_find_mpi_comm = parallel_info%mpi_comm_bilateral(i)
        EXIT LOOP
      END IF
    END DO LOOP
  END FUNCTION comin_parallel_find_mpi_comm


  !> Utility function.
  SUBROUTINE comin_parallel_handle_mpi_errcode(errcode)
    INTEGER, INTENT(IN) :: errcode
    INTEGER :: ierr
    IF (errcode .NE. MPI_SUCCESS) THEN
      WRITE (0,*) "Error in MPI program. Terminating."
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    END IF
  END SUBROUTINE comin_parallel_handle_mpi_errcode

END MODULE comin_parallel
