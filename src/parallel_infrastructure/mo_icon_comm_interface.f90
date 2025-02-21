! A collection of MPI communication tools
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_icon_comm_interface

  USE mo_parallel_config, ONLY: p_test_run, use_icon_comm
  USE mo_model_domain,    ONLY: t_patch
  USE mo_mpi,             ONLY: my_process_is_mpi_parallel, &
    & p_barrier, p_comm_work_test, p_comm_work
  USE mo_icon_comm_lib





  IMPLICIT NONE

  PRIVATE

  ! public constants
  PUBLIC :: construct_icon_communication
  PUBLIC :: destruct_icon_communication

  PUBLIC :: icon_comm_barrier



  !-------------------------------------------------------------------------
CONTAINS

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE construct_icon_communication(patch, n_dom)
    TYPE(t_patch), TARGET ::  patch(:)
    INTEGER :: n_dom

    INTEGER :: grid_id

    CHARACTER(*), PARAMETER :: method_name = "construct_icon_communication"

    !-------------------------------------------------------------------------------------
    ! constrcut patch communicators
    ! this is not the right place, this should be aware of the local decomposition features
    ! as well the coupled setup, parallel io, etc.
    ! Will be moved in the future
    patch(:)%compute_is_parallel = my_process_is_mpi_parallel()
    patch(:)%is_in_parallel_test = p_test_run

    patch(:)%work_communicator = p_comm_work
    patch(:)%parallel_test_communicator = p_comm_work_test
    ! WRITE(0,*) "work_communicator, test=", p_comm_work, p_comm_work_test
    !-------------------------------------------------------------------------------------
    IF (.NOT. use_icon_comm) RETURN

    CALL construct_icon_comm_lib()

    DO grid_id = 1, n_dom
      ! create the communication patterns
      CALL init_icon_std_comm_patterns(patch(grid_id))
    ENDDO

!     CALL work_mpi_barrier()
!     CALL finish("barrier returns", "")

    RETURN

  END SUBROUTINE construct_icon_communication
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE destruct_icon_communication()

    IF (.NOT. use_icon_comm) RETURN

    CALL destruct_icon_comm_lib()

  END SUBROUTINE destruct_icon_communication
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE icon_comm_barrier(for_patch)
    TYPE(t_patch), TARGET ::  for_patch

!    write(0,*) get_my_mpi_all_id(), ": enter icon_comm_barrier, for_patch%compute_is_parallel:",for_patch%compute_is_parallel, &
!      & for_patch%work_communicator
    ! CALL  flush(0)
    RETURN
    ! CALL  flush(0)
  END SUBROUTINE icon_comm_barrier
  !-----------------------------------------------------------------------

END MODULE mo_icon_comm_interface
