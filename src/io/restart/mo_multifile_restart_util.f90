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

! Utilities for multifile restart writing.
! These are mostly convenience functions for accessing the restart configuration.

#include "handle_mpi_error.inc"

MODULE mo_multifile_restart_util
  USE mo_kind,           ONLY: dp, sp
  USE mo_exception,      ONLY: finish,message_text
  USE mo_impl_constants, ONLY: SUCCESS, REAL_T, SINGLE_T, INT_T
  USE mo_util_libc,      ONLY: strerror
  USE mo_util_file,      ONLY: createSymlink
  USE mo_util_string,    ONLY: int2string
  USE mo_io_config,      ONLY: restartWritingParameters, ALL_WORKERS_INVOLVED
  USE mo_mpi,            ONLY: num_work_procs, p_comm_work_restart, p_comm_rank, &
    &                          p_real_dp, p_real_sp, p_int
#ifndef NOMPI
  HANDLE_MPI_ERROR_USE
  USE mpi,               ONLY: addr => MPI_ADDRESS_KIND
#else
  USE mo_kind,           ONLY: addr => i8
#endif


  IMPLICIT NONE

  PUBLIC :: multifileRestartLinkName, createMultifileRestartLink
  PUBLIC :: multifilePayloadPath
  PUBLIC :: isAsync, restartProcCount, rBuddy, rGroup, commonBuf_t
  PUBLIC :: iAmRestartMaster, iAmRestartWriter, initUtil

  PRIVATE

  TYPE :: commonBuf_t
    REAL(dp), POINTER :: d(:)
    REAL(sp), POINTER :: s(:)
    INTEGER,  POINTER :: i(:)
  END TYPE commonBuf_t

  CHARACTER(*), PARAMETER, PUBLIC :: vNames_glbIdx(3) = &
    ["global_cell_indices", "global_vert_indices", "global_edge_indices"]

  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_util"

  INTEGER, PUBLIC, PARAMETER :: typeID(3) = (/ REAL_T, SINGLE_T, INT_T /)
  INTEGER, PUBLIC, PARAMETER :: typeMax = MAX(REAL_T, SINGLE_T, INT_T)
  INTEGER, PUBLIC, SAVE :: typeMap(typeMax) = -1, mpiDtype(3) = -1
  INTEGER(addr), PUBLIC, SAVE :: typeByte(3) = -1, facTtoSp(3) = -1
  LOGICAL, SAVE :: before_init = .true.

CONTAINS

  SUBROUTINE initUtil()
    INTEGER :: i, ierr
    INTEGER(addr) :: tLB
    CHARACTER(*), PARAMETER :: routine = "multifileRestart_initUtil"

    IF (before_init) THEN
      mpiDtype(:) = [p_real_dp, p_real_sp, p_int]
      DO i = 1, 3
        typeMap(typeID(i)) = i
#ifndef NOMPI
        CALL MPI_Type_get_extent(mpiDtype(i), tLB, typeByte(i), ierr)
        HANDLE_MPI_ERROR(ierr, 'MPI_Type_get_extent')
#endif
      END DO
#ifdef NOMPI
      typeByte(:) = (/ 8, 4, 4 /) ! set some defaults (e.g. p_int_byte would be zero, if NOMPI)
#endif
      facTtoSp(:) = typeByte(:) / typeByte(2)
      IF (facTtoSp(1) .NE. (typeByte(1)+typeByte(2)-1_addr) / typeByte(2)) &
        CALL finish(routine, "sizeof(DP) not an integer multiple of sizeof(SP)!")
      IF (facTtoSp(3) .NE. (typeByte(2)+typeByte(3)-1_addr) / typeByte(3)) &
        CALL finish(routine, "sizeof(INT) not an integer multiple of sizeof(SP)!")
    END IF
  END SUBROUTINE initUtil

  SUBROUTINE multifileRestartLinkName(modelType, resultVar)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: modelType

    resultVar = "multifile_restart_"//modelType//".mfr"
  END SUBROUTINE multifileRestartLinkName

  SUBROUTINE createMultifileRestartLink(filename, modelType)
    CHARACTER(*), INTENT(IN) :: filename, modelType
    CHARACTER(:), ALLOCATABLE :: linkname
    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":createMultifileRestartLink"

    CALL multifileRestartLinkName(modelType, linkname)
    error = createSymlink(filename, linkname)
    IF(error /= SUCCESS) THEN
      WRITE(message_text, '(a,a)') "error creating symlink to restart file: ", strerror(error)
      CALL finish(routine,message_text)
    ENDIF
  END SUBROUTINE createMultifileRestartLink

  !XXX: this IS NOT the ONLY place where this path IS defined, it IS also generated/recognized IN c_restart_util.c
  SUBROUTINE multifilePayloadPath(multifilePath, jg, procId, resultVar)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: multifilePath
    INTEGER, INTENT(in) :: jg, procId

    resultVar = multifilePath//"/patch"//TRIM(int2string(jg))//"_"//TRIM(int2string(procId))//".nc"
  END SUBROUTINE multifilePayloadPath

  INTEGER FUNCTION restartProcCount() RESULT(resultVar)
    CALL restartWritingParameters(opt_restartProcCount = resultVar)
    IF (resultVar == ALL_WORKERS_INVOLVED)  resultVar = num_work_procs
  END FUNCTION restartProcCount

  LOGICAL FUNCTION isAsync() RESULT(resultVar)
    CALL restartWritingParameters(opt_lDedicatedProcMode = resultVar)
  END FUNCTION isAsync

  INTEGER FUNCTION rBuddy(pe_in, nr_in, nw_in) RESULT(res)
    INTEGER, OPTIONAL,INTENT(IN) :: pe_in, nr_in, nw_in
    INTEGER :: pe, nr, nw

    nw = num_work_procs
    nr = restartProcCount()
    pe = p_comm_rank(p_comm_work_restart)
    IF (PRESENT(nr_in)) nr = nr_in
    IF (PRESENT(pe_in)) pe = pe_in
    IF (PRESENT(nw_in)) nw = nw_in
    res = rGroup(pe_in=pe, nr_in=nr, nw_in=nw)
    IF (isAsync() .AND. .NOT.PRESENT(nr_in)) THEN
      res = nw + res
    ELSE
      res = CEILING(REAL(res * nw, dp)/REAL(nr, dp))
    END IF
  END FUNCTION rBuddy

  INTEGER FUNCTION rGroup(pe_in, nr_in, nw_in) RESULT(res)
    INTEGER, OPTIONAL,INTENT(IN) :: pe_in, nr_in, nw_in
    INTEGER :: pe, nr, nw

    nw = num_work_procs
    nr = restartProcCount()
    pe = p_comm_rank(p_comm_work_restart)
    IF (PRESENT(nr_in)) nr = nr_in
    IF (PRESENT(pe_in)) pe = pe_in
    IF (PRESENT(nw_in)) nw = nw_in
    res = MERGE(pe - nw, &
      &         FLOOR(REAL(nr * pe, dp)/REAL(nw, dp)), &
      &         isAsync() .AND. pe .GE. nw)
  END FUNCTION rGroup

  LOGICAL FUNCTION iAmRestartMaster() RESULT(resultVar)

    resultVar = rGroup() == rGroup(pe_in=0) .AND. &
      &         iAmRestartWriter()
  END FUNCTION iAmRestartMaster

  LOGICAL FUNCTION iAmRestartWriter() RESULT(resultVar)
    INTEGER :: the_pe

    the_pe = p_comm_rank(p_comm_work_restart)
    resultVar = rBuddy() == the_pe
  END FUNCTION iAmRestartWriter

END MODULE mo_multifile_restart_util
