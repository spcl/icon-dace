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

! Module for reading single file restart files
!
! Note: The asynchronous implementation of the restart output can be
!       found in the module "mo_async_restart"

#include "omp_definitions.inc"
MODULE mo_load_singlefile_restart
  USE mo_exception,          ONLY: message, warning, finish
  USE mo_impl_constants,     ONLY: SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
  USE mo_kind,               ONLY: dp, sp
  USE mo_model_domain,       ONLY: t_patch
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot, p_bcast, p_comm_work
  USE mo_restart_nml_and_att,ONLY: getAttributesForRestarting, ocean_initFromRestart_OVERRIDE
  USE mo_key_value_store,    ONLY: t_key_value_store
  USE mo_restart_util,       ONLY: restartSymlinkName
  USE mo_restart_var_data,   ONLY: get_var_3d_ptr, has_valid_time_level
  USE mo_var,                ONLY: t_var_ptr
  USE mo_timer, ONLY: timer_start, timer_stop, timer_load_restart_io, timers_level
  USE mo_util_string,        ONLY: int2string
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data, distrib_nf_open, &
    & distrib_read, distrib_nf_close, idx_lvl_blk
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf
  USE mo_dynamics_config,    ONLY: nnow, nnow_rcf
  USE mo_fortran_tools,      ONLY: t_ptr_3d, t_ptr_3d_int, t_ptr_3d_sp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: singlefileCheckRestartFiles, singlefileReadPatch

  CHARACTER(*), PARAMETER :: modname = "mo_load_singlefile_restart"

CONTAINS

  ! This IS just an existence check for the different patch files, a
  ! warning IS generated each domain for which no corresponding
  ! restart file IS found.
  !
  ! XXX: The code that I found READ the restart attributes from all
  ! the files, *appending* them to the list of attributes. Since all
  ! restart files contain the same attributes, this ONLY served to
  ! duplicate them. The current code just ignores the attributes IN
  ! all but the first file, just like the original code ignored the
  ! namelists IN all but the first file.  However, it might be a good
  ! idea to add some consistency checking on the attributes/namelists
  ! of the other files.
  SUBROUTINE singlefileCheckRestartFiles(modelType)
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER :: jg, n_dom
    LOGICAL :: lexists
    CHARACTER(:), ALLOCATABLE :: filename
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = modname//":singlefileCheckRestartFiles"

    CALL getAttributesForRestarting(restartAttributes)
    CALL restartAttributes%get('n_dom', n_dom)
    ! check whether we have all the restart files we expect
    DO jg = 1, n_dom
        CALL restartSymlinkName(modelType, jg, filename, n_dom)
        INQUIRE(file = filename, exist = lexists)
        IF (lexists) THEN
            CALL message(routine, "found restart file at '"//filename//"'")
        ELSE
            CALL warning(routine, 'domain '//TRIM(int2string(jg))//' is not active at restart time')
        END IF
    END DO
  END SUBROUTINE singlefileCheckRestartFiles

  SUBROUTINE singlefileReadPatch(vDat, modelType, ptc, ndom)
    TYPE(t_var_ptr), INTENT(in) :: vDat(:)
    CHARACTER(*), INTENT(IN) :: modelType
    TYPE(t_patch), TARGET, INTENT(IN) :: ptc
    INTEGER, INTENT(IN) :: ndom
    LOGICAL :: int_is_int
    TYPE(t_ptr_3d) :: r(1)
    TYPE(t_ptr_3d_sp) :: s(1)
    TYPE(t_ptr_3d_int) :: i(1)
    CHARACTER(:), ALLOCATABLE :: restart_filename
    INTEGER :: fID, nread, iV, dt, dummy
    TYPE(t_distrib_read_data), POINTER :: dio
    TYPE(t_key_value_store), POINTER :: restartAttributes
    LOGICAL :: skip(SIZE(vDat))
    CHARACTER(*), PARAMETER :: routine = modname // ":singlefileReadPatch"

    CALL restartSymlinkName(modelType, ptc%id, restart_filename, ndom)
    CALL getAttributesForRestarting(restartAttributes)
    CALL restartAttributes%get('int_is_int', int_is_int, opt_err=dummy)
    IF (dummy .NE. 0) int_is_int = .FALSE.
    IF(my_process_is_mpi_workroot()) THEN
      WRITE(0, "(a)") "opening " // restart_filename
      CALL nf(nf90_open(restart_filename, NF90_NOWRITE, fID), routine)
      DO iV = 1, SIZE(vDat)
        skip(iV) = .NOT.has_valid_time_level(vDat(iV)%p%info, ptc%id, nnow(ptc%id), nnow_rcf(ptc%id))
        IF (.NOT.skip(iV)) THEN
          skip(iV) = nf90_inq_varid(fID, TRIM(vDat(iV)%p%info%name), dummy) .NE. NF90_NOERR
          IF ((ocean_initFromRestart_OVERRIDE .OR. vDat(iV)%p%info%lrestart_cont) .AND. skip(iV)) THEN
            CALL warning(routine, "variable not found: '"//TRIM(vDat(iV)%p%info%name))
          ELSE IF (skip(iV)) THEN
            CALL finish(routine, "variable not found: "//TRIM(vDat(iV)%p%info%name))
          END IF
        END IF
      END DO
      CALL nf(nf90_close(fID), routine)
    END IF
    CALL p_bcast(skip, 0, comm=p_comm_work)
    fID = distrib_nf_open(restart_filename) 
    nread = 0
    DO iV = 1, SIZE(vDat)
      IF (skip(iV)) CYCLE
      SELECT CASE(vDat(iV)%p%info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        dio => ptc%cells%dist_io_data
      CASE (GRID_UNSTRUCTURED_VERT)
        dio => ptc%verts%dist_io_data
      CASE (GRID_UNSTRUCTURED_EDGE)
        dio => ptc%edges%dist_io_data
      CASE default
        CALL finish(routine, 'unknown grid type')
      END SELECT
      dt = vDat(iV)%p%info%data_type
      IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
      SELECT CASE(dt)
      CASE(REAL_T)
        CALL get_var_3d_ptr(vDat(iV)%p, r(1)%p)
        IF (SIZE(r(1)%p, 3) .GT. 0) r(1)%p(:,:,SIZE(r(1)%p, 3)) = 0._dp
        CALL distrib_read(fID, vDat(iV)%p%info%name, r, (/dio/), &
          & edim=(/SIZE(r(1)%p, 2)/), dimo=idx_lvl_blk)
      CASE(SINGLE_T)
        CALL get_var_3d_ptr(vDat(iV)%p, s(1)%p)
        IF (SIZE(s(1)%p, 3) .GT. 0) s(1)%p(:,:,SIZE(s(1)%p, 3)) = 0._sp
        CALL distrib_read(fID, vDat(iV)%p%info%name, s, (/dio/), &
          & edim=(/SIZE(s(1)%p, 2)/), dimo=idx_lvl_blk)
      CASE(INT_T)
        CALL get_var_3d_ptr(vDat(iV)%p, i(1)%p)
        IF (int_is_int) THEN
          IF (SIZE(i(1)%p, 3) .GT. 0) i(1)%p(:,:,SIZE(i(1)%p, 3)) = 0
          CALL distrib_read(fID, vDat(iV)%p%info%name, i, (/dio/), &
            & edim=(/SIZE(i(1)%p, 2)/), dimo=idx_lvl_blk)
        ELSE
          ALLOCATE(r(1)%p(SIZE(i(1)%p, 1), SIZE(i(1)%p, 2), SIZE(i(1)%p, 3)))
          IF (SIZE(r(1)%p, 3) .GT. 0) r(1)%p(:,:,SIZE(r(1)%p, 3)) = 0
          CALL distrib_read(fID, vDat(iV)%p%info%name, r, (/dio/), &
            & edim=(/SIZE(r(1)%p, 2)/), dimo=idx_lvl_blk)
          IF (SIZE(r(1)%p) .GT. 0 .AND. dt .EQ. INT_T) THEN
            !ICON_OMP PARALLEL WORKSHARE
            i(1)%p(:,:,:) = INT(r(1)%p(:,:,:))
            !ICON_OMP END PARALLEL WORKSHARE
          END IF
          IF (dt .EQ. INT_T) DEALLOCATE(r(1)%p)
        END IF
      CASE DEFAULT
        CALL finish(routine, "Internal error! Variable " // TRIM(vDat(iV)%p%info%name))
      END SELECT
      IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
      nread = nread + 1
#ifdef DEBUG
      IF (my_process_is_mpi_workroot()) WRITE (0,*) ' ... read ' // TRIM(vDat(iV)%p%info%name)
#endif
    END DO
    IF (my_process_is_mpi_workroot()) WRITE (0,'(a,i5,a)') ' ... read ', nread, ' variables'
    CALL distrib_nf_close(fID)
  END SUBROUTINE singlefileReadPatch

END MODULE mo_load_singlefile_restart
