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

! The multifile version of t_RestartPatchData.

#include "omp_definitions.inc"

MODULE mo_multifile_restart_patch_data
  USE mo_exception,                   ONLY: finish
  USE mo_impl_constants,              ONLY: SINGLE_T, REAL_T, INT_T
  USE mo_kind,                        ONLY: dp, sp, i8
  USE mo_multifile_restart_collector, ONLY: t_MultifileRestartCollector
  USE mo_multifile_restart_util,      ONLY: multifilePayloadPath, isAsync, vNames_glbIdx, commonBuf_t, &
    & typeMax, typeID, typeMap, typeByte
  USE mo_restart_patch_data,          ONLY: t_RestartPatchData
  USE mo_restart_patch_description,   ONLY: t_restart_patch_description
  USE mo_restart_var_data,            ONLY: has_valid_time_level, get_var_3d_ptr
  USE mo_timer,                       ONLY: timer_start, timer_stop, timer_write_restart_io,                   &
    &                                       timer_write_restart_communication, timers_level,                   &
    &                                       timer_write_restart_setup, timer_write_restart_wait
  USE mo_var_metadata_types,          ONLY: t_var_metadata
  USE mo_var,                         ONLY: t_var_ptr
  USE mo_parallel_config,             ONLY: restart_chunk_size
  USE mo_var_list_register_utils,     ONLY: vlr_select_restart_vars
  USE mo_netcdf_errhandler,           ONLY: nf
  USE mo_netcdf
  USE mo_restart_util,                ONLY: t_rfids
  USE mo_util_sort,                   ONLY: quicksort

  IMPLICIT NONE
  PRIVATE
  
  TYPE, EXTENDS(t_RestartPatchData), PUBLIC :: t_MultifilePatchData
    INTEGER :: cnkLvs
    LOGICAL :: shortcut
    TYPE(t_MultifileRestartCollector) :: coll
  CONTAINS
    PROCEDURE :: construct => multifilePatchData_construct
    PROCEDURE :: createCollectors => multifilePatchData_createCollectors
    PROCEDURE :: start_local_access  => multifilePatchData_local_access
    PROCEDURE :: start_remote_access => multifilePatchData_remote_access
    PROCEDURE :: exposeData  => multifilePatchData_exposeData    
    PROCEDURE :: destruct => multifilePatchData_destruct    ! override
    PROCEDURE :: writeData => multifilePatchData_writeData
    PROCEDURE :: fileStuff => multifilePatchData_fileStuff
  END TYPE t_MultifilePatchData
  
  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_patch_data"

CONTAINS

  SUBROUTINE multifilePatchData_writeData(me, ncid)
    CLASS(t_MultifilePatchData), INTENT(INOUT), TARGET :: me
    INTEGER,                     INTENT(IN)    :: ncid

    CALL finish(modname//"writeData", "not implemented!")
  END SUBROUTINE multifilePatchData_writeData

  SUBROUTINE multifilePatchData_construct(me, modelType, jg)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    CHARACTER(*),                INTENT(IN)    :: modelType
    INTEGER,                     INTENT(IN)    :: jg

    CALL me%description%init(jg)
    CALL vlr_select_restart_vars(me%varData, jg, modelType, me%restartType)
    IF (isAsync()) CALL me%description%transferToRestart()
  END SUBROUTINE multifilePatchData_construct

  SUBROUTINE multifilePatchData_createCollectors(me, wRnk, srcRnks, lactive)
    CLASS(t_MultifilePatchData), INTENT(INOUT), TARGET :: me
    INTEGER,                     INTENT(IN) :: wRnk, srcRnks(:)
    LOGICAL,                     INTENT(IN) :: lactive
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_createCollectors"
    INTEGER :: iV, jg, iT, nVar, i, j, tf(typeMax), tl(typeMax), nl(SIZE(me%varData)), perm(SIZE(me%varData))
    INTEGER(i8)                     :: iOffset(typeMax)
    TYPE(t_var_metadata), POINTER   :: ci
    TYPE(t_var_ptr), ALLOCATABLE :: varReordered(:)

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_setup)
    jg = me%description%id
    nVar = SIZE(me%varData)
    me%shortcut = .FALSE.
    IF (SIZE(srcRnks) .EQ. 1) &
      me%shortcut = .NOT.isAsync() .AND. wRnk .EQ. srcRnks(1)
    CALL me%coll%construct(jg, nVar, wRnk, srcRnks, lactive)
    ioffset(:) = 0_i8
    ALLOCATE(varReordered(nVar))
    i = 0; tf(:) = 0; tl(:) = 0
    DO iT = 1, typeMax
      DO iV = 1, nVar ! sort by type
        ci => me%varData(iV)%p%info
        IF (ci%data_type .EQ. iT) THEN
          i = i + 1
          IF (tf(iT) .EQ. 0) tf(iT) = i
          varReordered(i)%p => me%varData(iV)%p
          nl(i) = - MERGE(ci%used_dimensions(2), 1, ci%ndims .GT. 2)
        END IF
      END DO
      IF (tf(iT) .GT. 0) tl(iT) = i
    END DO
    IF (me%shortcut) THEN ! sort vars of a type by decending nlev -> fewer collect_write cycles
      CALL MOVE_ALLOC(varReordered, me%varData)
    ELSE
      perm(:) = [(j, j = 1, nVar)]
      DO iT = 1, typeMax
        IF (tf(iT) .EQ. 0) CYCLE
        CALL quicksort(nl(tf(iT):tl(iT)), perm(tf(iT):tl(iT)))
      END DO
      DO j = 1, nVar
        me%varData(j)%p => varReordered(perm(j))%p
      END DO
    END IF
    IF (i .NE. nVar) CALL finish(routine, "inconsistency!!!")
    DO iV = 1, nVar
      ci => me%varData(iV)%p%info
      IF (ci%ndims > 3) CALL finish(routine, "ndims > 3 is not supported")
      CALL me%coll%defVar(iV, -nl(iV), ci%data_type, me%description%hmap(ci%hgrid), iOffset)
    END DO
    me%cnkLvs = MERGE(restart_chunk_size, -MINVAL(nl(:),1), restart_chunk_size .GT. 0)
    CALL me%coll%init_win(iOffset, me%shortcut)
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_setup)
  END SUBROUTINE multifilePatchData_createCollectors

  SUBROUTINE multifilePatchData_fileStuff(me, fname, ifile, bWritten, date_dayas)
    CLASS(t_MultifilePatchData), TARGET, INTENT(INOUT) :: me
    CHARACTER(*),                        INTENT(IN)    :: fname
    INTEGER,                             INTENT(IN)    :: ifile
    REAL(dp),                            INTENT(IN)    :: date_dayas
    INTEGER(i8),                         INTENT(INOUT) :: bWritten
    CHARACTER(:), ALLOCATABLE                          :: effectiveFilename
    CLASS(t_restart_patch_description), POINTER        :: desc
    INTEGER :: ncid, iVarId(3), dum, tvid
    CHARACTER(*), PARAMETER :: routine = modname//":fileStuff"

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    desc => me%description
    CALL multifilePayloadPath(fname, desc%id, ifile, effectiveFilename)
    CALL nf(nf90_create(effectiveFilename, NF90_NETCDF4, ncid), routine)
    CALL setup_meta()
    CALL nf(nf90_set_fill(ncid, NF90_NOFILL, dum), routine)
    CALL nf(nf90_enddef(ncid), routine)
    CALL write_glbids(iVarId)
    CALL nf(nf90_put_var(ncid, tvid, date_dayas, [1]), routine)
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
    CALL collect_write()
    CALL nf(nf90_close(ncid), routine)
  CONTAINS

    SUBROUTINE setup_meta()
      INTEGER :: iV
      TYPE(t_rfids) :: rfids
      TYPE(t_var_metadata), POINTER :: ci
      CHARACTER(*), PARAMETER :: routine = modname//":fileStuff:setup_meta"

      CALL rfids%init(ncid, [me%coll%idx(1)%nRecv, me%coll%idx(2)%nRecv, me%coll%idx(3)%nRecv], tvid)
      DO iV = 1, 3
        CALL nf(nf90_def_var(ncid, vNames_glbIdx(iV), NF90_INT, [rfids%gids(iV), rfids%ftid], iVarId(iV)), routine)
      END DO
      DO iV = 1, SIZE(me%varData) !SIZE(vWrNow)
        ci => me%varData(iV)%p%info
        IF (.NOT.has_valid_time_level(ci, me%description%id, &
          &  me%description%nnow, me%description%nnow_rcf)) CYCLE
        CALL rfids%def_ncdfvar(ci, desc%hmap(ci%hgrid))
      END DO
    END SUBROUTINE setup_meta

    SUBROUTINE write_glbids(vids)
      INTEGER, INTENT(IN) :: vids(3)
      INTEGER :: i
      INTEGER, TARGET :: dummy_i(0)
      INTEGER, POINTER :: buf_i(:)

      DO i = 1, 3
        buf_i => dummy_i
        IF (me%coll%idx(i)%nRecv .GT. 0) buf_i => me%coll%glb_idx(i)%p
        CALL nf(nf90_put_var(ncid, vids(i), buf_i, [1,1], [me%coll%idx(i)%nRecv,1]), &
          &     modname//":fileStuff:write_glbids")
        bWritten = bWritten + INT(SIZE(buf_i), i8) * 4_i8
      END DO
    END SUBROUTINE write_glbids

    SUBROUTINE collect_write()
      TYPE(t_var_metadata), POINTER   :: ci
      INTEGER :: vO, i, iV, sub(3), cv, nd, st(3), ct(3)
      TYPE(commonBuf_t) :: rBuf
      INTEGER, ALLOCATABLE :: vS(:), vWrNow(:)
      LOGICAL, DIMENSION(SIZE(me%varData)) :: vWrDone
      INTEGER(i8), ALLOCATABLE :: srcOff(:,:)
      CHARACTER(*), PARAMETER :: routine = modname//":fileStuff:collect_write"

      DO i = 1, 3
        IF (.NOT.ASSOCIATED(me%coll%glb_idx(i)%p)) RETURN
        IF (SIZE(me%coll%glb_idx(i)%p) .LE. 0) RETURN
      END DO
      vWrDone(:) = .false.; sub(:) = 0
      DO WHILE (ANY(.NOT.vWrDone(:)))
        CALL select_vars(me, vWrDone, vWrNow, sub)
        IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
        CALL me%coll%fetch(vWrNow, sub, vS, rBuf, srcOff)
        IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
        IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
        vO = 0
        DO iV = 1, SIZE(vWrNow)
          cv = vWrNow(iV)
          ci => me%varData(cv)%p%info
          nd = ci%ndims
          IF (nd .EQ. 2) THEN
            st(1:2) = [1,1]
            ct(1:2) = [vS(iV),1]
          ELSE IF (sub(3) .EQ. cv) THEN
            st(1:3) = [1,sub(1)+1,1]
            ct(1:3) = [vS(iV)/(sub(2)-sub(1)),sub(2)-sub(1),1]
            IF (sub(2) .EQ. ci%used_dimensions(2)) sub(:) = 0
          ELSE
            st(1:3) = [1,1,1]
            ct(1:3) = [vS(iV)/ci%used_dimensions(2),ci%used_dimensions(2),1]
          END IF
          SELECT CASE(ci%data_type)
          CASE(REAL_T)
            CALL nf(nf90_put_var(ncid, ci%cdiVarId, rBuf%d(1+vO:vO+vS(iV)), st(:nd), ct(:nd)), routine)
          CASE(SINGLE_T)
            CALL nf(nf90_put_var(ncid, ci%cdiVarId, rBuf%s(1+vO:vO+vS(iV)), st(:nd), ct(:nd)), routine)
          CASE(INT_T)
            CALL nf(nf90_put_var(ncid, ci%cdiVarId, rBuf%i(1+vO:vO+vS(iV)), st(:nd), ct(:nd)), routine)
          END SELECT
          vO = vO + vS(iV)
          bWritten = bWritten + INT(vS(iV), i8) * typeByte(typeMap(ci%data_type))
        END DO
        IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
      END DO
    END SUBROUTINE collect_write

  END SUBROUTINE multifilePatchData_fileStuff

  SUBROUTINE multifilePatchData_local_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL me%coll%local_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
  END SUBROUTINE multifilePatchData_local_access

  SUBROUTINE multifilePatchData_remote_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL me%coll%remote_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
  END SUBROUTINE multifilePatchData_remote_access

  SUBROUTINE select_vars(me, vWrDone, vWrNow, sub)
    TYPE(t_MultifilePatchData), INTENT(IN) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":select_vars"
    LOGICAL, INTENT(INOUT) :: vWrDone(SIZE(me%varData))
    INTEGER, ALLOCATABLE, INTENT(OUT) :: vWrNow(:)
    INTEGER, INTENT(INOUT) :: sub(3)
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER :: lFree, iV, lN, tFrst, nW
    LOGICAL :: l_vWrNow(SIZE(me%varData))

    l_vWrNow(:) = .false.
    lFree = me%cnkLvs
    tFrst = -1
    nW = 0
    DO iV = 1, SIZE(me%varData)
      IF (lFree .EQ. 0) EXIT
      IF (vWrDone(iV)) CYCLE
      curInfo => me%varData(iV)%p%info
      IF (.NOT.has_valid_time_level(curInfo, me%description%id, &
        &        me%description%nnow, me%description%nnow_rcf)) THEN
        vWrDone(iV) = .true. ! nothing to do
        CYCLE
      END IF
      IF (tFrst .EQ. -1) THEN
        tFrst = curInfo%data_type
        IF (.NOT.ANY(tFrst .EQ. typeID(:))) CALL finish(routine, &
          & "data type not recognized! Variable" // TRIM(curInfo%name))
      END IF
      IF (tFrst .NE. curInfo%data_type) EXIT
      lN = MERGE(curInfo%used_dimensions(2), 1, curInfo%ndims .EQ. 3)
      IF (lFree .GE. lN .OR. me%shortcut) THEN
        IF (.NOT.me%shortcut) lFree = lFree - lN
        vWrDone(iV) = .true.
        IF (lN .NE. 0) l_vWrNow(iV) = .true.
      ELSE IF (me%cnkLvs .LT. lN .AND. (sub(3) .EQ. iV .OR. sub(3) .EQ. 0)) THEN
        l_vWrNow(iV) = .true.
        sub(1) = sub(2)
        sub(2) = MIN(sub(2) + me%cnkLvs, lN)
        sub(3) = iV
        lFree = lFree - sub(2) + sub(1)
        IF (sub(2) .EQ. lN) vWrDone(iV) = .true.
      END IF
      IF (l_vWrNow(iV)) nW = nW + 1
    END DO
    IF (ALLOCATED(vWrNow)) DEALLOCATE(vWrNow)
    ALLOCATE(vWrNow(nW))
    IF (nW .GT. 0) vWrNow(:) = PACK([(iV, iV = 1, SIZE(me%varData))], l_vWrNow(:))
  END SUBROUTINE select_vars

  SUBROUTINE multifilePatchData_exposeData(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    TYPE(t_var_metadata), POINTER :: ci
    INTEGER :: iV, sub(3), lCnt, cv
    INTEGER(i8) :: offset(3)
    REAL(KIND=dp), POINTER :: r_ptr_3d(:,:,:)
    REAL(KIND=sp), POINTER :: s_ptr_3d(:,:,:)
    INTEGER, POINTER :: i_ptr_3d(:,:,:)
    LOGICAL, DIMENSION(SIZE(me%varData)) :: vWrDone
    INTEGER, ALLOCATABLE :: vWrNow(:)

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    vWrDone(:) = .false.; offset(:) = 0_i8; sub(:) = 0
    DO WHILE (ANY(.NOT.vWrDone(:)))
      CALL select_vars(me, vWrDone, vWrNow, sub)
      DO iV = 1, SIZE(vWrNow)
        cv = vWrNow(iV)
        ci => me%varData(cv)%p%info
        lCnt = MERGE(ci%used_dimensions(2), 1, ci%ndims == 3)
        IF (sub(3) .EQ. cv) THEN
          IF (sub(2) .EQ. lCnt) sub(:) = 0
          IF (sub(3) .NE. 0) CYCLE
        END IF
        SELECT CASE (ci%data_type)
        CASE(REAL_T)
          CALL get_var_3d_ptr(me%varData(cv)%p, r_ptr_3d)
          CALL me%coll%sendField(cv, r_ptr_3d, offset(1))
        CASE(SINGLE_T)
          CALL get_var_3d_ptr(me%varData(cv)%p, s_ptr_3d)
          CALL me%coll%sendField(cv, s_ptr_3d, offset(2))
        CASE(INT_T)
          CALL get_var_3d_ptr(me%varData(cv)%p, i_ptr_3d)
          CALL me%coll%sendField(cv, i_ptr_3d, offset(3))
        END SELECT
      END DO
    END DO
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  END SUBROUTINE multifilePatchData_exposeData

  SUBROUTINE multifilePatchData_destruct(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    CALL me%coll%finalize()
    IF (ALLOCATED(me%varData)) DEALLOCATE(me%varData)
  END SUBROUTINE multifilePatchData_destruct

END MODULE mo_multifile_restart_patch_data
