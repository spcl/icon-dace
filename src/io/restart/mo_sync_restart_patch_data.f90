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

MODULE mo_sync_restart_patch_data

  USE mo_communication,             ONLY: t_comm_gather_pattern, exchange_data
  USE mo_exception,                 ONLY: finish, message
  USE mo_dynamics_config,           ONLY: nnow, nnow_rcf
  USE mo_impl_constants,            ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_kind,                      ONLY: dp, sp
  USE mo_mpi,                       ONLY: my_process_is_mpi_workroot, my_process_is_mpi_test
  USE mo_restart_patch_data,        ONLY: t_RestartPatchData
  USE mo_restart_var_data,          ONLY: get_var_3d_ptr, has_valid_time_level
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_restart_io, &
                                        & timer_write_restart_communication, timers_level
  USE mo_var_metadata_types,        ONLY: t_var_metadata
  USE mo_var_list_register_utils,   ONLY: vlr_select_restart_vars
  USE mo_netcdf_errhandler,         ONLY: nf
  USE mo_netcdf

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_SyncPatchData

  TYPE, EXTENDS(t_RestartPatchData) :: t_SyncPatchData
  CONTAINS
    PROCEDURE :: construct => syncPatchData_construct
    PROCEDURE :: destruct => syncPatchData_destruct
    PROCEDURE :: writeData => syncPatchData_writeData  ! override
  END TYPE t_SyncPatchData

  CHARACTER(*), PARAMETER :: modname = 'mo_sync_restart_patch_data'

CONTAINS

  SUBROUTINE syncPatchData_destruct(me)
    CLASS(t_SyncPatchData), INTENT(INOUT) :: me

    IF (ALLOCATED(me%varData)) DEALLOCATE(me%varData)
  END SUBROUTINE syncPatchData_destruct

  SUBROUTINE syncPatchData_construct(me, modelType, jg)
    CLASS(t_syncPatchData), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, INTENT(IN) :: jg

    CALL me%description%init(jg)
    CALL vlr_select_restart_vars(me%varData, jg, modelType, me%restartType)
  END SUBROUTINE syncPatchData_construct

  ! loop over all var_lists for restart
  SUBROUTINE syncPatchData_writeData(me, ncid)
    CLASS(t_SyncPatchData), INTENT(INOUT), TARGET :: me
    INTEGER, INTENT(IN) :: ncid
    INTEGER :: domain, i, gridSize, error, level, st(3), ct(3), nd
    TYPE(t_var_metadata), POINTER        :: info
    TYPE(t_comm_gather_pattern), POINTER :: gatherPattern
    REAL(dp), ALLOCATABLE                :: gatherBuffer_dp(:)
    REAL(sp), ALLOCATABLE                :: gatherBuffer_sp(:)
    INTEGER, ALLOCATABLE                 :: gatherBuffer_int(:)
    REAL(dp), POINTER                    :: r_ptr_3d(:,:,:)
    REAL(sp), POINTER                    :: s_ptr_3d(:,:,:)
    INTEGER, POINTER                     :: i_ptr_3d(:,:,:)
    CHARACTER(*), PARAMETER              :: routine = modname//":syncPatchData_writeData"
    LOGICAL :: is_mpi_workroot

    is_mpi_workroot = my_process_is_mpi_workroot()
    IF(my_process_is_mpi_test()) RETURN
    domain = me%description%id
    DO i = 1, SIZE(me%varData)
      info => me%varData(i)%p%info
      IF(.NOT.has_valid_time_level(info, domain, nnow(domain), nnow_rcf(domain))) CYCLE
#ifdef DEBUG
      IF(is_mpi_workroot) write (0,*)' ... write '//TRIM(info%name)
#endif
      gridSize = me%description%n_patch_elem_g(me%description%hmap(info%hgrid))
      gatherPattern => me%description%gpat(me%description%hmap(info%hgrid))%p
      ! get pointers to the local DATA
      nd = info%ndims
      ct(:) = [gridSize,1,1]
      SELECT CASE(info%data_type)
      CASE(REAL_T)
        ALLOCATE(gatherBuffer_dp(MERGE(gridSize, 0, is_mpi_workroot)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        CALL get_var_3d_ptr(me%varData(i)%p, r_ptr_3d)
        ! gather the data in the master process and write it to disk
        DO level = 1, SIZE(r_ptr_3d, 2)
          IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
          CALL exchange_data(in_array = r_ptr_3d(:,level,:), &
            &                out_array = gatherBuffer_dp, &
            &                gather_pattern = gatherPattern)
          IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
          IF(is_mpi_workroot) THEN
            st(:) = [1,level,1]
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
            CALL nf(nf_put_vara_double(ncid, info%cdiVarID, st(:nd), ct(:nd), gatherBuffer_dp), routine)
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
          END IF
        END DO
        DEALLOCATE(gatherBuffer_dp)
      CASE(SINGLE_T)
        ALLOCATE(gatherBuffer_sp(MERGE(gridSize, 0, is_mpi_workroot)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        CALL get_var_3d_ptr(me%varData(i)%p, s_ptr_3d)
        ! gather the data in the master process and write it to disk
        DO level = 1, SIZE(s_ptr_3d, 2)
          IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
          CALL exchange_data(in_array = s_ptr_3d(:,level,:), &
            &                out_array = gatherBuffer_sp, &
            &                gather_pattern = gatherPattern)
          IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
          IF(is_mpi_workroot) THEN
            st(:) = [1,level,1]
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
            CALL nf(nf_put_vara_real(ncid, info%cdiVarID, st(:nd), ct(:nd), gatherBuffer_sp), routine)
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
          END IF
        END DO
        DEALLOCATE(gatherBuffer_sp)
      CASE(INT_T)
        ALLOCATE(gatherBuffer_int(MERGE(gridSize, 0, is_mpi_workroot)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        CALL get_var_3d_ptr(me%varData(i)%p, i_ptr_3d)
        ! gather the data in the master process and write it to disk
        DO level = 1, SIZE(i_ptr_3d, 2)
          IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
          CALL exchange_data(in_array = i_ptr_3d(:,level,:), &
            &                out_array = gatherBuffer_int, &
            &                gather_pattern = gatherPattern)
          IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
          IF(is_mpi_workroot) THEN
            st(:) = [1,level,1]
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
            CALL nf(nf_put_vara_int(ncid, info%cdiVarID, st(:nd), ct(:nd), gatherBuffer_int), routine)
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
          END IF
        END DO
        DEALLOCATE(gatherBuffer_int)
      CASE DEFAULT
        CALL finish(routine, "Internal error! Variable "//TRIM(info%name))
      END SELECT
    END DO
    CALL message('syncPatchData_writeData','Finished writing restart')
  END SUBROUTINE syncPatchData_writeData

END MODULE mo_sync_restart_patch_data
