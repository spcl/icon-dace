! Test plugin for the ICON Community Interface (ComIn)
! test the parallel descriptive data of ComIn
!
!  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

MODULE test_parallel_plugin

  USE comin_plugin_interface, ONLY: comin_parallel_get_host_mpi_comm, &
       & t_comin_plugin_info, t_comin_setup_version_info, &
       & comin_setup_get_version, comin_plugin_finish, &
       & t_comin_descrdata_domain, comin_descrdata_get_domain
  USE mpi,                    ONLY: MPI_IN_PLACE, MPI_Allreduce, MPI_Abort, MPI_INTEGER, MPI_SUM, &
       & MPI_COMM_WORLD, MPI_MAX

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: pluginname = "test_parallel_plugin"

  !> working precision (will be compared to ComIn's and ICON's)
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  TYPE(t_comin_setup_version_info) :: version

CONTAINS
  SUBROUTINE test_parallel_plugin_setup()  BIND(C)
    USE iso_c_binding, ONLY: C_BOOL, C_INT

    INTEGER :: ierr, no_owned_cells, host_comm, max_glb_index
    TYPE(t_comin_descrdata_domain), POINTER :: p_patch

    version = comin_setup_get_version()
    IF (version%version_no_major > 1) CALL comin_plugin_finish("parallel test", "version mismatch")

    p_patch => comin_descrdata_get_domain(1)
    IF (.NOT. ASSOCIATED(p_patch)) CALL comin_plugin_finish("parallel test ", "parallel data missing")

    ! check if the sum of all owned cells matches the global size
    no_owned_cells = COUNT(p_patch%cells%decomp_domain .EQ. 0)
    host_comm = comin_parallel_get_host_mpi_comm()
    CALL MPI_Allreduce(MPI_IN_PLACE, no_owned_cells, 1, MPI_INTEGER, MPI_SUM, host_comm, ierr)
    IF (no_owned_cells .NE. p_patch%cells%ncells_global) &
      & CALL comin_plugin_finish("parallel test", "There are more owned cells as global_size allows")

    ! check if the maximal global index matches the global size
    max_glb_index = MAXVAL(p_patch%cells%glb_index)
    CALL MPI_Allreduce(MPI_IN_PLACE, max_glb_index, 1, MPI_INTEGER, MPI_MAX, host_comm, ierr)
    IF(max_glb_index .NE. p_patch%cells%ncells_global) &
      & CALL comin_plugin_finish("parallel test", "The maximal global index is larger than global_size")

  END SUBROUTINE test_parallel_plugin_setup

END MODULE test_parallel_plugin
