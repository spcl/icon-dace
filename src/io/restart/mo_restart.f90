! This is the high-level module that is used to write a restart.
! Basically it reexports t_RestartDescriptor and provides a factory function to create restart descriptors.
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

MODULE mo_restart



  USE mo_exception, ONLY: finish, message
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_io_config, ONLY: restartWritingParameters, kSyncRestartModule, kAsyncRestartModule, kMultifileRestartModule
  USE mo_mpi, ONLY: stop_mpi, my_process_is_restart, process_mpi_restart_size
  USE mo_multifile_restart, ONLY: t_MultifileRestartDescriptor



  USE mo_restart_descriptor, ONLY: t_RestartDescriptor
  USE mo_sync_restart, ONLY: t_SyncRestartDescriptor




  USE mo_timer, ONLY: print_timer, timer_stop, timer_model_init, ltimer

  IMPLICIT NONE
  PRIVATE
  ! documentation for t_RestartDescriptor IS found IN mo_restart_descriptor
  PUBLIC :: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor, detachRestartProcs

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart"

CONTAINS

  ! Factory FUNCTION to create the appropriate restart descriptor.
  FUNCTION createRestartDescriptor(modelType) RESULT(resultVar)
    CHARACTER(*), INTENT(IN) :: modelType
    CLASS(t_RestartDescriptor), POINTER :: resultVar
    INTEGER :: error, restartModule
    LOGICAL :: lDedicatedProcMode
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":createRestartDescriptor"

    CALL restartWritingParameters(opt_restartModule      = restartModule, &
      &                           opt_lDedicatedProcMode = lDedicatedProcMode)
    SELECT CASE(restartModule)
    CASE(kSyncRestartModule)
      CALL message('','synchronous restart writing selected.')
      ALLOCATE(t_SyncRestartDescriptor :: resultVar, STAT = error)
    CASE(kAsyncRestartModule)
      CALL message('','asynchronous restart writing selected.')

      CALL finish(routine, "this executable was compiled without MPI support, hence async restart writing is not &
                           &available")



    CASE(kMultifileRestartModule)
      IF (lDedicatedProcMode) THEN
        CALL message('','multifile restart writing selected, with dedicated procs.')
      ELSE
        CALL message('','multifile restart writing selected, joint proc mode.')
      END IF
      ALLOCATE(t_MultifileRestartDescriptor :: resultVar, STAT = error)
    END SELECT
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    CALL resultVar%construct(modelType)
  END FUNCTION createRestartDescriptor

  ! Convenience FUNCTION for destroying a restart descriptor.
  SUBROUTINE deleteRestartDescriptor(descriptor)
    CLASS(t_RestartDescriptor), POINTER, INTENT(INOUT) :: descriptor

    CALL descriptor%destruct()
    DEALLOCATE(descriptor)
  END SUBROUTINE deleteRestartDescriptor

  ! Enter the restart main proc IF this IS a pure restart PE AND set use_async_restart_output accordingly.
  !
  ! This routine does NOT RETURN on dedicated restart processes.
  SUBROUTINE detachRestartProcs(timer_started)
    LOGICAL, INTENT(IN) :: timer_started
    INTEGER :: restartModule
    CHARACTER(*), PARAMETER :: routine = modname//":detachRestartProcs"

    IF(process_mpi_restart_size <= 0) RETURN    ! no dedicated restart processes configured -> noop
    IF(.NOT.my_process_is_restart()) RETURN ! this IS NOT a dedicated restart process -> noop

    CALL finish(routine, 'no MPI-support -> no restart procs -> go away!')
  END SUBROUTINE detachRestartProcs

END MODULE mo_restart
