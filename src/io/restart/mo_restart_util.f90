! Contains common helper routines for(a)synchronous restart
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

MODULE mo_restart_util
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_fortran_tools,      ONLY: assign_if_present_allocatable
  USE mo_impl_constants,     ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T, MAX_CHAR_LENGTH
  USE mo_io_config,          ONLY: restartWritingParameters, kMultifileRestartModule
  USE mo_kind,               ONLY: i8, dp
  USE mo_packed_message,     ONLY: t_PackedMessage, kPackOp
  USE mo_run_config,         ONLY: restart_filename
  USE mo_util_libc,          ONLY: strerror
  USE mo_util_file,          ONLY: get_filename_noext, createSymlink
  USE mo_util_string,        ONLY: int2string, associate_keyword, with_keywords, t_keyword_list
  USE mtime,                 ONLY: datetime, newDatetime, deallocateDatetime
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_comm_work, p_bcast
#ifndef NOMPI
  USE mo_mpi, ONLY: p_pe_work, my_process_is_restart
  USE mpi, ONLY: MPI_PROC_NULL, MPI_ROOT
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_restart_args, t_rfids
  PUBLIC :: getRestartFilename, restartSymlinkName, create_restart_file_link
  PUBLIC :: restartBcastRoot, check_for_checkpoint
  ! patch independent restart arguments
  TYPE t_restart_args
    TYPE(datetime), POINTER :: restart_datetime => NULL()
    INTEGER :: jstep
    CHARACTER(LEN=32) :: modelType
    INTEGER, ALLOCATABLE :: output_jfile(:)
  CONTAINS
    PROCEDURE :: construct => restartArgs_construct
    PROCEDURE :: packer => restartArgs_packer   ! unpacking IS considered construction
    PROCEDURE :: destruct => restartArgs_destruct
  END TYPE t_restart_args

  TYPE t_rfids
    LOGICAL :: isinit = .FALSE.
    INTEGER :: ncid, nlayids, ftid, gids(3)
    INTEGER, ALLOCATABLE :: layids(:), nlays(:)
  CONTAINS
    PROCEDURE :: init => rfids_init
    PROCEDURE :: def_ncdfvar => rfids_def_ncdfvar
  END TYPE t_rfids

  CHARACTER(*), PARAMETER :: modname = "mo_restart_util"

CONTAINS

  ! Broadcast root for intercommunicator broadcasts from compute PEs to restart
  ! PEs using p_comm_work_2_restart.
  INTEGER FUNCTION restartBcastRoot() RESULT(resultVar)
    resultVar = 0
#ifndef NOMPI
    ! Special root setting for intercommunicators:
    ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL.
    IF (.NOT.my_process_is_restart()) &
      & resultVar = MERGE(MPI_ROOT, MPI_PROC_NULL, p_pe_work == 0)
#endif
  END FUNCTION restartBcastRoot

  SUBROUTINE getRestartFilename(baseName, jg, restartArgs, resultVar, date_dayas)
    CHARACTER(*), INTENT(IN) :: baseName
    INTEGER, INTENT(IN) :: jg
    TYPE(t_restart_args), INTENT(IN) :: restartArgs
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    REAL(dp), INTENT(OUT) :: date_dayas
    CHARACTER(LEN=32) :: datetimeString
    INTEGER :: restartModule
    TYPE(t_keyword_list), POINTER :: keywords
    TYPE(datetime), POINTER :: dt
    INTEGER :: date_int
    REAL(dp) :: date_frac
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: restart_fname 

    dt => restartArgs%restart_datetime
    WRITE (datetimeString,'(i4.4,2(i2.2),a,3(i2.2),a)')    &
       & dt%date%year, dt%date%month, dt%date%day , 'T', &
       & dt%time%hour, dt%time%minute, dt%time%second, 'Z'
    date_int = dt%date%year*10000 + dt%date%month*100 + dt%date%day
    date_frac = (dt%time%hour*3600.0 + dt%time%minute*60.0 + dt%time%second)/86400.0
    date_dayas = date_int+date_frac
    NULLIFY(keywords)
    ! build the keyword list
    CALL associate_keyword("<gridfile>", TRIM(get_filename_noext(baseName)), keywords)
    CALL associate_keyword("<idom>", TRIM(int2string(jg, "(i2.2)")), keywords)
    CALL associate_keyword("<rsttime>", TRIM(datetimeString), keywords)
    CALL associate_keyword("<mtype>", TRIM(restartArgs%modelType), keywords)
    CALL restartWritingParameters(opt_restartModule = restartModule)
    IF(restartModule == kMultifileRestartModule) THEN
      CALL associate_keyword("<extension>", "mfr", keywords)
    ELSE
      CALL associate_keyword("<extension>", "nc", keywords)
    END IF
    ! replace keywords in file name
    resultVar = TRIM(with_keywords(keywords, TRIM(restart_filename)))
  END SUBROUTINE getRestartFilename

  SUBROUTINE restartSymlinkName(modelType, jg, resultVar, opt_ndom)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, VALUE :: jg
    INTEGER, INTENT(IN) :: opt_ndom
    CHARACTER(2) :: sg

    IF (opt_ndom == 1) jg = 1
    WRITE (sg, '(i2.2)') jg
    resultVar = 'restart_'//modelType//"_DOM"//sg//'.nc'
  END SUBROUTINE restartSymlinkName

  SUBROUTINE create_restart_file_link(filename, modelType, jg, opt_ndom)
    CHARACTER(*), INTENT(IN) :: filename, modelType
    INTEGER, INTENT(in) :: jg
    INTEGER, INTENT(IN) :: opt_ndom
    INTEGER :: ierr
    CHARACTER(:), ALLOCATABLE :: linkname
    CHARACTER(*), PARAMETER :: routine = modname//':create_restart_file_link'

    CALL restartSymlinkName(modelType, jg, linkname, opt_ndom)
    ierr = createSymlink(filename, linkname)
    IF(ierr /= SUCCESS) THEN
      WRITE(message_text,'(a,a,a,a)') "error creating symlink at ",linkname,":", strerror(ierr)
      CALL finish(routine,message_text)
    ENDIF
  END SUBROUTINE create_restart_file_link

  SUBROUTINE restartArgs_construct(me, this_datetime, jstep, modelType, opt_output_jfile)
    CLASS(t_restart_args), INTENT(INOUT) :: me
    TYPE(datetime), INTENT(IN) :: this_datetime
    INTEGER, INTENT(in) :: jstep
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    integer :: ierr

    me%restart_datetime => newDatetime(this_datetime%date%year, this_datetime%date%month,    &
      &                                this_datetime%date%day, this_datetime%time%hour,      &
      &                                this_datetime%time%minute, this_datetime%time%second, &
      &                                this_datetime%time%ms, ierr)
    me%jstep = jstep
    me%modelType = modelType
    CALL assign_if_present_allocatable(me%output_jfile, opt_output_jfile)
  END SUBROUTINE restartArgs_construct

  SUBROUTINE restartArgs_packer(me, operation, packedMessage)
    CLASS(t_restart_args), INTENT(INOUT) :: me
    INTEGER, INTENT(in) :: operation
    TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage
    CHARACTER(len=*), PARAMETER ::  routine = modname//':restartArgs_packer'
    
    IF (operation == kPackOp .AND. .NOT. ASSOCIATED(me%restart_datetime)) THEN
      CALL finish(routine, 'Assertion failed: cannot pack unconstructed object.')
    ENDIF
    IF (.NOT. ASSOCIATED(me%restart_datetime)) THEN
      me%restart_datetime => newDatetime(1878_i8,1,1,0,0,0,0)
    ENDIF
    CALL packedMessage%packer(operation, me%restart_datetime%date%year)
    CALL packedMessage%packer(operation, me%restart_datetime%date%month)
    CALL packedMessage%packer(operation, me%restart_datetime%date%day)
    CALL packedMessage%packer(operation, me%restart_datetime%time%hour)
    CALL packedMessage%packer(operation, me%restart_datetime%time%minute)
    CALL packedMessage%packer(operation, me%restart_datetime%time%second)
    CALL packedMessage%packer(operation, me%jstep)
    CALL packedMessage%packer(operation, me%modelType)
    CALL packedMessage%packer(operation, me%output_jfile)
  END SUBROUTINE restartArgs_packer

  SUBROUTINE restartArgs_destruct(me)
    CLASS(t_restart_args), INTENT(INOUT) :: me

    IF(ASSOCIATED(me%restart_datetime)) CALL deallocateDatetime(me%restart_datetime)
    IF(ALLOCATED(me%output_jfile)) DEALLOCATE(me%output_jfile)
  END SUBROUTINE restartArgs_destruct

  SUBROUTINE rfids_init(rfids, ncid_in, nelem, tvid)
    CLASS(t_rfids), INTENT(OUT) :: rfids
    INTEGER, INTENT(IN) :: ncid_in, nelem(3)
    INTEGER, INTENT(OUT) :: tvid
    INTEGER :: dummy
    CHARACTER(*), PARAMETER :: routine = modname//":rfids_init"

    rfids%ncid = ncid_in
    IF (rfids%isinit) DEALLOCATE(rfids%layids, rfids%nlays)
    rfids%nlayids = 0
    ALLOCATE(rfids%layids(32), rfids%nlays(32))
    CALL nf(nf90_def_dim(rfids%ncid, "time", 1, rfids%ftid), routine)
    CALL nf(nf90_def_var(rfids%ncid, "time", NF90_DOUBLE, rfids%ftid, tvid), routine)
    CALL nf(nf90_put_att(rfids%ncid, tvid, "axis", "T"), routine)
    CALL nf(nf90_put_att(rfids%ncid, tvid, "units", "day as %Y%m%d.%f"), routine)
    CALL nf(nf90_def_dim(rfids%ncid, "cells", nelem(1), rfids%gids(1)), routine)
    CALL nf(nf90_def_var(rfids%ncid, "cells", NF90_INT, rfids%gids(1), dummy), routine)
    CALL nf(nf90_def_dim(rfids%ncid, "verts", nelem(2), rfids%gids(2)), routine)
    CALL nf(nf90_def_var(rfids%ncid, "verts", NF90_INT, rfids%gids(2), dummy), routine)
    CALL nf(nf90_def_dim(rfids%ncid, "edges", nelem(3), rfids%gids(3)), routine)
    CALL nf(nf90_def_var(rfids%ncid, "edges", NF90_INT, rfids%gids(3), dummy), routine)
    rfids%isinit = .TRUE.
  END SUBROUTINE rfids_init

  SUBROUTINE rfids_def_ncdfvar(rfids, ci, gid)
    CLASS(t_rfids), INTENT(INOUT) :: rfids
    INTEGER, INTENT(IN) :: gid
    TYPE(t_var_metadata), POINTER, INTENT(INOUT) :: ci
    INTEGER :: dtid, ivid, iivid, dids(3), nd, la, dummy
    INTEGER, ALLOCATABLE :: tmpa(:), tmpb(:)
    CHARACTER(LEN=10) :: layname
    CHARACTER(*), PARAMETER :: routine = modname//":rfids_def_ncdfvar"

    IF (.NOT.rfids%isinit) CALL finish(routine, "not initialized")
    nd = ci%ndims
    dids(1) = rfids%gids(gid)
    SELECT CASE(ci%data_type)
    CASE(REAL_T)
      dtid = NF90_DOUBLE
    CASE(SINGLE_T)
      dtid = NF90_REAL
    CASE(INT_T)
      dtid = NF90_INT
    END SELECT
    SELECT CASE(nd)
    CASE(2)
      dids(2) = rfids%ftid
    CASE(3)
      ivid = 0
      DO iivid = 1, rfids%nlayids
        IF (rfids%nlays(iivid) .EQ. ci%used_dimensions(2)) THEN
          ivid = iivid
          EXIT
        END IF
      END DO
      IF (ivid .EQ. 0) THEN
        IF (rfids%nlayids .EQ. SIZE(rfids%nlays)) THEN
          CALL MOVE_ALLOC(rfids%nlays, tmpa)
          CALL MOVE_ALLOC(rfids%layids, tmpb)
          ALLOCATE(rfids%nlays(rfids%nlayids+32), rfids%layids(rfids%nlayids+32))
          rfids%nlays(:rfids%nlayids) = tmpa
          rfids%layids(:rfids%nlayids) = tmpb
        END IF
        rfids%nlayids = rfids%nlayids + 1
        rfids%nlays(rfids%nlayids) = ci%used_dimensions(2)
        WRITE(layname, "(a,i0)") "layers_", rfids%nlays(rfids%nlayids)
        CALL nf(nf90_def_dim(rfids%ncid, TRIM(layname), rfids%nlays(rfids%nlayids), &
          & rfids%layids(rfids%nlayids)), routine)
        CALL nf(nf90_def_var(rfids%ncid, TRIM(layname), NF90_INT, &
          &                  rfids%layids(rfids%nlayids:rfids%nlayids), dummy), routine)
        CALL nf(nf90_put_att(rfids%ncid, dummy, "axis", "z"), routine)
        ivid = rfids%nlayids
      END IF
      dids(2) = rfids%layids(ivid)
      dids(3) = rfids%ftid
    END SELECT
    CALL nf(nf90_def_var(rfids%ncid, TRIM(ci%name), dtid, dids(:nd), ci%cdiVarId), routine)
    la = LEN_TRIM(ci%cf%units)
    CALL nf(nf90_put_att(rfids%ncid, ci%cdiVarId, "units", ci%cf%units(1:la)), routine)
    la = LEN_TRIM(ci%cf%long_name)
    CALL nf(nf90_put_att(rfids%ncid, ci%cdiVarId, "long_name", ci%cf%long_name(1:la)), routine)
  END SUBROUTINE rfids_def_ncdfvar

  SUBROUTINE check_for_checkpoint(lready_for_checkpoint, lchkp_allowed, lstop_on_demand)

    LOGICAL, INTENT(IN)    :: lready_for_checkpoint
    LOGICAL, INTENT(INOUT) :: lchkp_allowed
    LOGICAL, INTENT(OUT)   :: lstop_on_demand

    IF (my_process_is_stdio()) THEN
      ! The purpose of this file is to inform other processes that ICON is ready for checkpointing
      IF (lready_for_checkpoint .AND. .NOT. lchkp_allowed) THEN
        OPEN(123,file='ready_for_checkpoint')
        lchkp_allowed = .TRUE.
        CLOSE(123)
      ENDIF
      IF (lready_for_checkpoint) THEN
        INQUIRE(file='stop_icon',exist=lstop_on_demand)
      ELSE
        lstop_on_demand = .FALSE.
      ENDIF
    ENDIF
    CALL p_bcast(lstop_on_demand,p_io,p_comm_work)

  END SUBROUTINE check_for_checkpoint

END MODULE mo_restart_util
