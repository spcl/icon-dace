! Module for writing restart files (synchronously)
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
!
!
! Note: The asynchronous implementation of the restart output can be
!       found in the module "mo_async_restart"
!
! --------------------------------------------------------------------------------
!
! Generated files:
! ================
!
! 1. For each domain 1, ..., n_dom, and for each restart output time step:
!
!    Restart data file
!    -----------------
!      "<gridfile>_restart_<modeltype>_<timestamp>.nc",     e.g.
!      "iconR2B06_DOM01_restart_atm_20110101T001200Z.nc"    (NetCDF format)
!
!    This filename can be customized using the namelist parameter
!    "mo_run_nml/restart_filename".
!
!    This file contains
!    -  data
!    -  namelists
!    -  several attributes
!
!    Note:
!    -  We read the namelists only once and assume that these
!       are identical for all domains.
!    -  Since we do not know about the total number of domains at startup,
!       we have to ask the current restart file for the attribute "n_dom"
!
! 2. For each domain 1, ..., n_dom, and for the LAST restart output time step:
!
!    Symbolic link to data file
!    --------------------------
!      "restart_<modeltype>_DOMxx.nc"
!
!    Note:
!    -  The domain-dependent suffix "...DOMxx" is also required for non-nested setups.
!
! --------------------------------------------------------------------------------
!
!OPTION! -pvctl conflict!

MODULE mo_sync_restart

  USE mtime,                        ONLY: datetime
  USE mo_exception,                 ONLY: finish
  USE mo_grid_config,               ONLY: n_dom
  USE mo_impl_constants,            ONLY: SUCCESS
  USE mo_restart_descriptor,        ONLY: t_RestartDescriptor
  USE mo_restart_util,              ONLY: t_restart_args
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_restart, timer_write_restart_io, &
                                        & timer_write_restart_communication, timers_level

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_SyncRestartDescriptor

  TYPE, EXTENDS(t_RestartDescriptor) :: t_SyncRestartDescriptor
  CONTAINS
    PROCEDURE :: construct => syncRestartDescriptor_construct   ! override
    PROCEDURE :: writeRestart => syncRestartDescriptor_writeRestart ! override
    PROCEDURE :: destruct => syncRestartDescriptor_destruct ! override
  END TYPE t_SyncRestartDescriptor

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sync_restart'

  !------------------------------------------------------------------------------------------------
CONTAINS

  SUBROUTINE syncRestartDescriptor_construct(me, modelType)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT), TARGET :: me
    INTEGER :: jg, error
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartDescriptor_construct'

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    IF(timers_level >= 7) THEN
      CALL timer_start(timer_write_restart_io)
      CALL timer_stop(timer_write_restart_io)
      CALL timer_start(timer_write_restart_communication)
      CALL timer_stop(timer_write_restart_communication)
    END IF
    me%modelType = modelType
    ALLOCATE(me%sPatchData(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    me%patchData => me%sPatchData
    DO jg = 1, n_dom
      CALL me%sPatchData(jg)%construct(modelType, jg)
    END DO
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE syncRestartDescriptor_construct

  SUBROUTINE syncRestartDescriptor_writeRestart(me, this_datetime, jstep, opt_output_jfile)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT), TARGET :: me
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    INTEGER :: jg
    TYPE(t_restart_args) :: restartArgs

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    DO jg = 1, n_dom
      CALL me%sPatchData(jg)%description%setTimeLevels() !update the time levels
    END DO
    CALL restartArgs%construct(this_datetime, jstep, me%modelType, opt_output_jfile)
    CALL me%writeFiles(restartArgs, isSync=.true.)
    CALL restartArgs%destruct()
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE syncRestartDescriptor_writeRestart

  SUBROUTINE syncRestartDescriptor_destruct(me)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    INTEGER :: i

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    IF (ALLOCATED(me%sPatchData)) THEN
      DO i = 1, SIZE(me%sPatchData, 1)
        CALL me%sPatchData(i)%destruct()
      END DO
    END IF
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE syncRestartDescriptor_destruct

END MODULE mo_sync_restart
