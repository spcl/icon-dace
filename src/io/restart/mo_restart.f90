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
#ifndef NOMPI
  USE mo_async_restart, ONLY: t_AsyncRestartDescriptor, asyncRestart_mainLoop
#endif
  USE mo_exception, ONLY: finish, message
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_io_config, ONLY: restartWritingParameters, kSyncRestartModule, kAsyncRestartModule, kMultifileRestartModule
  USE mo_mpi, ONLY: stop_mpi, my_process_is_restart, process_mpi_restart_size
  USE mo_multifile_restart, ONLY: t_MultifileRestartDescriptor
#ifndef NOMPI
  USE mo_multifile_restart, ONLY: multifileRestart_mainLoop
#endif
  USE mo_restart_descriptor, ONLY: t_RestartDescriptor
  USE mo_sync_restart, ONLY: t_SyncRestartDescriptor
#ifdef YAC_coupling
  USE mo_coupling_config, ONLY: is_coupled_run
  USE mo_io_coupling_frame, ONLY: construct_io_coupling, destruct_io_coupling
#endif
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
#ifdef NOMPI
      CALL finish(routine, "this executable was compiled without MPI support, hence async restart writing is not &
                           &available")
#else
      ALLOCATE(t_AsyncRestartDescriptor :: resultVar, STAT = error)
#endif
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
#ifdef NOMPI
    CALL finish(routine, 'no MPI-support -> no restart procs -> go away!')
#else
    ! Actually detach the restart processes.
#ifdef YAC_coupling
    ! The initialisation of YAC needs to be called by all (!) MPI processes
    ! in MPI_COMM_WORLD. Thus we do it here for the restart processes.
    IF ( is_coupled_run() ) CALL construct_io_coupling ( "dummy" )
#endif
    IF (timer_started) CALL timer_stop(timer_model_init)
    CALL restartWritingParameters(opt_restartModule = restartModule)
    SELECT CASE(restartModule)
    CASE(kSyncRestartModule)
      CALL finish(routine, "sync mode, but on restart proc !?!")
    CASE(kAsyncRestartModule)
      CALL asyncRestart_mainLoop()
    CASE(kMultifileRestartModule)
      CALL multifileRestart_mainLoop()
    END SELECT
#ifdef YAC_coupling
    IF ( is_coupled_run() ) CALL destruct_io_coupling ( "dummy" )
#endif
    ! This is the end of all things!
    IF(ltimer) CALL print_timer
    CALL stop_mpi
#endif
  END SUBROUTINE detachRestartProcs

END MODULE mo_restart
