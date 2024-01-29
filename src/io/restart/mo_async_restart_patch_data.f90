! Contains routines for asynchronous restart output
!
! Note: The synchronous implementation of the restart output can be
!       found in the module "mo_sync_restart". See module header for
!       more details on generated files.
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

MODULE mo_async_restart_patch_data
#ifndef NOMPI
  USE ISO_C_BINDING,                ONLY: C_PTR, C_F_POINTER, C_LOC
  USE mo_async_restart_comm_data,   ONLY: t_AsyncRestartCommData
  USE mo_exception,                 ONLY: finish
  USE mo_kind,                      ONLY: wp, dp, sp
  USE mo_impl_constants,            ONLY: SINGLE_T, REAL_T, INT_T
  USE mo_parallel_config,           ONLY: crcs => restart_chunk_size
  USE mo_run_config,                ONLY: msg_level
  USE mo_restart_var_data,          ONLY: has_valid_time_level
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_restart_io, timers_level
  USE mo_mpi,                       ONLY: num_work_procs, my_process_is_restart, p_mpi_wtime,  &
    &                                     p_real_sp_byte, p_real_dp_byte, p_int_byte
  USE mpi,                          ONLY: MPI_ADDRESS_KIND
  USE mo_restart_patch_data,        ONLY: t_RestartPatchData
  USE mo_var_list_register_utils,   ONLY: vlr_select_restart_vars
  USE mo_netcdf_errhandler,         ONLY: nf
  USE mo_netcdf
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_var_metadata_types,        ONLY: t_var_metadata

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_AsyncPatchData

  TYPE, EXTENDS(t_RestartPatchData) :: t_AsyncPatchData
    TYPE(t_AsyncRestartCommData) :: commData
  CONTAINS
    PROCEDURE :: construct => asyncPatchData_construct
    PROCEDURE :: writeData => asyncPatchData_writeData
    PROCEDURE :: destruct  => asyncPatchData_destruct
  END TYPE t_AsyncPatchData

  CHARACTER(*), PARAMETER :: modname = 'mo_async_restart_patch_data'

CONTAINS

  ! collective across restart AND worker PEs
  SUBROUTINE asyncPatchData_construct(me, modelType, jg)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, INTENT(IN) :: jg

    CALL me%description%init(jg)
    CALL vlr_select_restart_vars(me%varData, jg, modelType, me%restartType)
    CALL me%commData%construct(jg, me%varData)
    CALL me%description%transferToRestart()
  END SUBROUTINE asyncPatchData_construct

  SUBROUTINE asyncPatchData_destruct(me)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me

    CALL me%commData%destruct()
    IF(ALLOCATED(me%varData)) DEALLOCATE(me%varData)
  END SUBROUTINE asyncPatchData_destruct

  !------------------------------------------------------------------------------------------------
  !
  ! Write restart variable list for a restart PE.
  !
  SUBROUTINE asyncPatchData_writeData(me, ncid)
    CLASS(t_AsyncPatchData), INTENT(INOUT), TARGET :: me
    INTEGER, INTENT(IN) :: ncid
    INTEGER(MPI_ADDRESS_KIND) :: ioff(0:num_work_procs-1), bGet, bWrite
    REAL(dp), ALLOCATABLE, TARGET :: buf(:)
    REAL(dp), CONTIGUOUS, POINTER :: p_dp(:,:)
    REAL(sp), CONTIGUOUS, POINTER :: p_sp(:,:)
    INTEGER, CONTIGUOUS, POINTER :: p_i(:,:)
    INTEGER :: ichunk, rcs, maxl, iv, nlevs, nd, st(3), ct(3)
    REAL(dp) :: t_get, t_write
    CHARACTER(*), PARAMETER :: routine = modname//':asyncPatchData_writeData'
    TYPE(c_ptr) :: cptr
    TYPE(t_restart_patch_description), POINTER :: desc
    TYPE(t_var_metadata), POINTER :: ci

    IF (.NOT. my_process_is_restart()) CALL finish(routine, 'Must be called on a restart PE!')
    t_get = 0.d0; t_write = 0.d0
    bGet = 0_mpi_address_kind; bWrite = 0_mpi_address_kind
    maxl = 0
    desc => me%description
    DO iv = 1, SIZE(me%varData)
      ci => me%varData(iv)%p%info
      maxl = MAX(maxl, MERGE(1, ci%used_dimensions(2), ci%ndims .EQ. 2))
    END DO
    rcs = MERGE(MIN(crcs, maxl), maxl, crcs > 0)
    ALLOCATE(buf(me%commData%maxLevelSize*rcs))
    cptr = C_LOC(buf(1))
    ioff(:) = 0
    st(:) = 1
    ct(:) = 1
    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(me%varData)
      ci => me%varData(iv)%p%info
      IF (.NOT.has_valid_time_level(ci, desc%id, desc%nnow, desc%nnow_rcf)) CYCLE
      nd = ci%ndims
      nlevs = MERGE(1, ci%used_dimensions(2), nd .EQ.2)
      ct(1) = desc%n_patch_elem_g(desc%hmap(ci%hgrid))
      ! loop over all chunks (of levels)
      DO ichunk = 1, (nlevs-1)/rcs + 1
        st(2) = (ichunk-1)*rcs + 1
        ct(2) = MIN(st(2)+rcs-1, nlevs) - st(2) + 1
        SELECT CASE(ci%data_type)
        CASE(REAL_T)
          CALL C_F_POINTER(cptr, p_dp, ct(1:2))
          CALL me%commData%collectData(desc%hmap(ci%hgrid), ct(2), p_dp, ioff, t_get, bGet)
          t_write = t_write - p_mpi_wtime()
          IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
          CALL nf(nf_put_vara_double(ncid, ci%cdiVarID, st(:nd), ct(:nd), p_dp), routine)
          bWrite = bWrite + ct(1)*ct(2)*p_real_dp_byte
        CASE(SINGLE_T)
          CALL C_F_POINTER(cptr, p_sp, ct(1:2))
          CALL me%commData%collectData(desc%hmap(ci%hgrid), ct(2), p_sp, ioff, t_get, bGet)
          t_write = t_write - p_mpi_wtime()
          IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
          CALL nf(nf_put_vara_real(ncid, ci%cdiVarID, st(:nd), ct(:nd), p_sp), routine)
          bWrite = bWrite + ct(1)*ct(2)*p_real_sp_byte
        CASE(INT_T)
          CALL C_F_POINTER(cptr, p_i, ct(1:2))
          CALL me%commData%collectData(desc%hmap(ci%hgrid), ct(2), p_i, ioff, t_get, bGet)
          t_write = t_write - p_mpi_wtime()
          IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
          CALL nf(nf_put_vara_int(ncid, ci%cdiVarID, st(:nd), ct(:nd), p_i), routine)
          bWrite = bWrite + ct(1)*ct(2)*p_int_byte
        CASE DEFAULT
          CALL finish(routine, "Internal error! Variable "//TRIM(ci%name))
        END SELECT
        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
        t_write = t_write + p_mpi_wtime()
      END DO
    END DO
    IF (msg_level >= 7) &
      & WRITE (0,'(10(a,f10.3))') ' Checkpointing: Got ', REAL(bGet, dp)*1.d-6, ' MB, time get: ', t_get, ' s [', &
           & REAL(bGet, dp)*1.d-6/MAX(1.e-6_wp, t_get), ' MB/s], time write: ', t_write, ' s [', &
           & REAL(bWrite, dp)*1.d-6/MAX(1.e-6_wp,t_write), ' MB/s]'
  END SUBROUTINE asyncPatchData_writeData
#endif
END MODULE mo_async_restart_patch_data
