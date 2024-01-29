! Module for collecting the payload data onto the respective restart
! processes.
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
#include "handle_mpi_error.inc"
#include "omp_definitions.inc"
MODULE mo_multifile_restart_collector

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
#ifndef NOMPI
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr
  HANDLE_MPI_ERROR_USE

  USE mpi, ONLY: addr => MPI_ADDRESS_KIND, MPI_DATATYPE_NULL, MPI_STATUS_IGNORE, &
    & MPI_STATUSES_IGNORE, MPI_TYPECLASS_INTEGER, MPI_COMM_NULL, MPI_INFO_NULL, &
    & MPI_MODE_NOCHECK, MPI_LOCK_EXCLUSIVE, MPI_WIN_NULL, MPI_UNDEFINED, &
    & MPI_Win_lock, MPI_Win_unlock, MPI_Win_free, &
    & MPI_Comm_split, MPI_Comm_size, MPI_Comm_free, &
    & MPI_Type_vector, MPI_Type_create_struct, MPI_Type_match_size, MPI_Type_commit, MPI_Type_free, &
    & MPI_Sizeof
# ifndef NO_MPI_CPTR_ARG
  USE mpi, ONLY: MPI_Alloc_mem
# endif
# ifndef NO_MPI_CHOICE_ARG
  ! MPI_Win_{un,}lock_all don't have choice args but OpenMPI doesn't include an interface in its
  ! use-mpi-tkr variant for compilers that can't ignore argument TKR via directives.
  ! Cray's MPI doesn't export MPI_Waitall and MPI_Waitany.
  USE mpi, ONLY: MPI_Win_create, MPI_Free_mem, MPI_Gather, MPI_Win_lock_all, MPI_Win_unlock_all, &
    & MPI_Waitall, MPI_Waitany
#   ifndef NO_MPI_RGET
  USE mpi, ONLY: MPI_Rget
#   else
  USE mpi, ONLY: MPI_Get
# endif
# endif /* !defined(NO_MPI_CHOICE_ARG) */
  USE mo_impl_constants, ONLY: SINGLE_T, REAL_T, INT_T
  USE mo_mpi, ONLY: p_int, p_barrier
  USE mo_multifile_restart_util, ONLY: mpiDtype, typeByte, typeMap
#else
  USE mo_kind, ONLY: addr => i8
#endif
  USE mo_communication, ONLY: idx_no, blk_no
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
  USE mo_model_domain, ONLY: p_patch
  USE mo_exception, ONLY: finish
  USE mo_kind, ONLY: dp, sp, i8
  USE mo_mpi, ONLY: p_comm_work_restart, p_comm_rank, p_send, p_recv, my_process_is_work
  USE mo_multifile_restart_util, ONLY: iAmRestartWriter, commonBuf_t, &
    & typeMax, typeID, facTtoSP
  USE mo_timer, ONLY: timer_start, timer_stop, timer_restart_collector_setup, &
    & timer_restart_indices_setup, timers_level
  USE mo_fortran_tools, ONLY: t_ptr_1d_int, t_ptr_1d_sp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_MultifileRestartCollector

  TYPE :: t_CollectorIndices
    INTEGER :: nRecv, nSend = 0
    INTEGER, ALLOCATABLE :: nSrcRecv(:), sIdx(:), sBlk(:)
  END TYPE t_CollectorIndices

  INTEGER, PARAMETER :: nOpenReqMax = 16
  REAL(dp), TARGET :: dummy_d(0)
  REAL(sp), TARGET :: dummy_s(0)
  INTEGER,  TARGET :: dummy_i(0)

#ifdef NOMPI
  INTEGER, PARAMETER :: MPI_COMM_NULL = 0, MPI_WIN_NULL = 0
#endif

  TYPE :: t_MultifileRestartCollector
    PRIVATE
    TYPE(commonBuf_t) :: sBuf, rBuf
    INTEGER(addr), ALLOCATABLE :: tOffSv(:)
    INTEGER(addr) :: rBuf_size = 0_addr
    INTEGER :: wComm = MPI_COMM_NULL, win = MPI_WIN_NULL, nVar, destPE, nSrcPE
    LOGICAL :: allocd = .false., wPosted = .false., wStarted = .false., shortcut
    TYPE(t_CollectorIndices), PUBLIC :: idx(3)
    TYPE(t_ptr_1d_int), PUBLIC :: glb_idx(3)
    TYPE(t_ptr_1d_sp) :: wptr
    INTEGER, ALLOCATABLE :: vGrid(:), vType(:), vLevs(:), srcPE(:)
  CONTAINS
    PROCEDURE :: construct     => multifileRestartCollector_construct
    PROCEDURE :: init_win      => multifileRestartCollector_init_win
    PROCEDURE :: defVar        => multifileRestartCollector_defVar
    PROCEDURE :: finalize      => multifileRestartCollector_finalize
    PROCEDURE :: sendField     => multifileRestartCollector_sendField
    PROCEDURE :: fetch         => multifileRestartCollector_fetch
    PROCEDURE :: local_access  => multifileRestartCollector_local_access
    PROCEDURE :: remote_access => multifileRestartCollector_remote_access
  END TYPE t_MultifileRestartCollector

  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_collector"

CONTAINS

  ! On the restart writers, this returns an array with the global
  ! indices of the points collected by this PE. The return value is
  ! not allocated on pure worker procs.
  FUNCTION collectorIndices_init(me, nLoc, decompInfo, destPE, srcPE, lactive) &
    & RESULT(globalIndices)
    INTEGER, POINTER :: globalIndices(:)
    CLASS(t_CollectorIndices),       INTENT(INOUT) :: me
    INTEGER,                         INTENT(IN) :: nLoc, destPE, srcPE(:)
    TYPE(t_grid_domain_decomp_info), INTENT(IN), POINTER :: decompInfo
    LOGICAL,                         INTENT(IN) :: lactive
    INTEGER :: myRank, i, j, nSrcPE
    INTEGER, ALLOCATABLE :: sBuf_i(:)
    IF (timers_level >= 10)  CALL timer_start(timer_restart_indices_setup)
    myRank = p_comm_rank(p_comm_work_restart)
    nSrcPE = SIZE(srcPE)
    IF (my_process_is_work() .AND. lactive) THEN
      DO i = 1, nLoc
        IF(decompInfo%owner_mask(idx_no(i), blk_no(i))) &
          & me%nSend = me%nSend + 1
      END DO
    END IF
    ALLOCATE(me%sIdx(me%nSend), me%sBlk(me%nSend), sBuf_i(me%nSend), &
      & me%nSrcRecv(nSrcPE))
    IF (my_process_is_work() .AND. lactive) THEN
      j = 1
      DO i = 1, nLoc
        IF (decompInfo%owner_mask(idx_no(i), blk_no(i))) THEN
          me%sIdx(j) = idx_no(i)
          me%sBlk(j) = blk_no(i)
          sBuf_i(j) = decompInfo%glb_index(i)
          j = j + 1
        END IF
      END DO
    END IF
    IF (iAmRestartWriter()) THEN
      DO i = 1, nSrcPE
        IF(srcPE(i) == myRank) THEN
          me%nSrcRecv(i) = me%nSend
        ELSE
          CALL p_recv(me%nSrcRecv(i), srcPE(i), 0, comm=p_comm_work_restart)
        END IF
      END DO
    END IF
    IF (destPE /= myRank .AND. lactive) &
      & CALL p_send(me%nSend, destPE, 0, comm=p_comm_work_restart)
    me%nRecv = SUM(me%nSrcRecv(1:nSrcPE))
    ALLOCATE(globalIndices(me%nRecv))
    IF (lactive) THEN
      IF(destPE .NE. myRank .AND. me%nSend .GT. 0) &
          CALL p_send(sBuf_i, destPE, 0, comm=p_comm_work_restart)
      j = 1
      globalIndices(:) = 0._dp
      DO i = 1, nSrcPE
        IF(srcPE(i) .EQ. myRank) THEN
          globalIndices(j:j-1+me%nSrcRecv(i)) = sBuf_i(1:me%nSend)
        ELSE
          IF (me%nSrcRecv(i) > 0) &
            & CALL p_recv(globalIndices(j:j-1+me%nSrcRecv(i)), srcPE(i), 0, &
              &          comm = p_comm_work_restart)
        END IF
        j = j + me%nSrcRecv(i)
      END DO
    END IF
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_indices_setup)
  END FUNCTION collectorIndices_init

  SUBROUTINE MultifileRestartCollector_construct(this, jg, nVar, destPE, srcPE, lactive)
    CLASS(t_MultifileRestartCollector),  INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: jg, nVar, destPE, srcPE(:)
    LOGICAL, INTENT(IN) :: lactive
    TYPE(t_grid_domain_decomp_info), POINTER :: deco
    INTEGER :: i, nElem

    this%nSrcPE = SIZE(srcPE)
    ALLOCATE(this%srcPE(this%nSrcPE))
    this%srcPE(:) = srcPE(:)
    this%destPE = destPE
    IF (timers_level >= 10)  CALL timer_start(timer_restart_collector_setup)
    DO i = 1, 3
      nElem = 0
      NULLIFY(deco)
      IF (my_process_is_work()) THEN
        SELECT CASE(i)
        CASE(1)
          nElem = p_patch(jg)%n_patch_cells
          deco => p_patch(jg)%cells%decomp_info
        CASE(2)
          nElem = p_patch(jg)%n_patch_verts
          deco => p_patch(jg)%verts%decomp_info
        CASE(3)
          nElem = p_patch(jg)%n_patch_edges
          deco => p_patch(jg)%edges%decomp_info
        END SELECT
      END IF
      this%glb_idx(i)%p => collectorIndices_init(this%idx(i), nElem, deco, destPE, srcPE, lactive)
    END DO
    this%sBuf%d => dummy_d
    this%sBuf%s => dummy_s
    this%sBuf%i => dummy_i
    this%rBuf%d => dummy_d
    this%rBuf%s => dummy_s
    this%rBuf%i => dummy_i
    this%nVar = nVar
    ALLOCATE(this%vGrid(nVar), this%vType(nVar), this%vLevs(nVar))
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_collector_setup)
  END SUBROUTINE multifileRestartCollector_construct

  SUBROUTINE multifileRestartCollector_defVar(this, iVar, nLevs, iType, iGrid, iOffset)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: this
    INTEGER,                            INTENT(IN   ) :: iVar, nLevs, iType, iGrid
    INTEGER(KIND=i8),                   INTENT(INOUT) :: iOffset(:)

    IF (timers_level >= 10)  CALL timer_start(timer_restart_collector_setup)
    this%vGrid(iVar) = igrid
    this%vType(iVar) = iType
    this%vLevs(iVar) = nLevs
    ioffset(itype) = ioffset(itype) + this%idx(this%vGrid(iVar))%nSend * nLevs
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_collector_setup)
  END SUBROUTINE multifileRestartCollector_defVar

  SUBROUTINE multifileRestartCollector_finalize(this)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: this
    INTEGER :: i
#ifndef NOMPI
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_finalize"
    INTEGER :: rRank, ierr
#endif

    IF (ALLOCATED(this%srcPE)) DEALLOCATE(this%srcPE)
    DO i = 1, 3
      IF (ALLOCATED(this%idx(i)%nSrcRecv)) &
        & DEALLOCATE(this%idx(i)%nSrcRecv, this%glb_idx(i)%p)
      IF (ALLOCATED(this%idx(i)%sIdx)) DEALLOCATE(this%idx(i)%sIdx, this%idx(i)%sBlk)
    END DO
    IF (ALLOCATED(this%vGrid)) DEALLOCATE(this%vGrid, this%vType, this%vLevs)
    IF (this%allocd) THEN
#ifndef NOMPI
      IF (.NOT.this%shortcut) THEN
        IF (this%wStarted .OR. this%wPosted) &
          & CALL this%local_access()
        IF (my_process_is_work()) THEN
          rrank = p_comm_rank(this%wComm)
          CALL MPI_Win_unlock(rrank, this%win, ierr)
          HANDLE_MPI_ERROR(ierr, 'MPI_Win_unlock')
        END IF
      END IF
      IF (this%win .NE. MPI_WIN_NULL) THEN
        CALL MPI_Win_free(this%win, ierr)
        HANDLE_MPI_ERROR(ierr, 'MPI_Win_free')
      END IF
      CALL MPI_Free_mem(this%wptr%p, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Free_mem')
#else
      DEALLOCATE(this%wptr%p)
#endif
      NULLIFY(this%wptr%p, this%sBuf%d, this%sBuf%s, this%sBuf%i)
      this%allocd = .false.
    END IF
#ifndef NOMPI
    IF (this%wComm .NE. MPI_COMM_NULL) THEN
      CALL MPI_Comm_free(this%wComm, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Comm_free')
      this%wComm = MPI_COMM_NULL
    END IF
    this%win = MPI_WIN_NULL
    IF (this%rBuf_size .GT. 0) THEN
      CALL MPI_Free_mem(this%rBuf%d, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Free_mem')
      NULLIFY(this%rBuf%d, this%rBuf%s, this%rBuf%i)
      this%rBuf_size = 0_addr
    END IF
#endif
    IF (ALLOCATED(this%tOffSv)) DEALLOCATE(this%tOffSv)
    this%wStarted = .false.
    this%wPosted = .false.
  END SUBROUTINE multifileRestartCollector_finalize

  SUBROUTINE multifileRestartCollector_fetch(me, vWrNow, sub, vSize, cBuf, srcOffset)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,           INTENT(IN   ) :: vWrNow(:), sub(3)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: vSize(:)
    TYPE(commonBuf_t), INTENT(  OUT) :: cBuf
    INTEGER(i8), ALLOCATABLE, INTENT(INOUT) :: srcOffset(:,:)
    INTEGER, DIMENSION(SIZE(vWrNow)) :: stride
    INTEGER :: iV, iV_o, nL, nV
#ifndef NOMPI
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_fetch"
    INTEGER(addr) :: t_dsp, fDsp, tBy, bSizeB
    INTEGER :: i, ierr, o_dTy, dTy, n_oReq, iReq, iTy, bSize
#ifndef NO_MPI_RGET
    INTEGER :: getReq(nOpenReqMax)
#endif
    INTEGER(addr), DIMENSION(SIZE(vWrNow)) :: rDsp, rD
    INTEGER, DIMENSION(SIZE(vWrNow)) :: ones, tmp_ddt
    INTEGER :: pct(me%nSrcPE,SIZE(vWrNow)), t_cnt(me%nSrcPE)
    TYPE(c_ptr) :: cptr

    bSize = 0; ones(:) = 1; stride(:) = 0; t_cnt(:) = 0; rDsp(:) = 0_addr
    iTy = me%vType(vWrNow(1))
#endif
    nV = SIZE(vWrNow)
    IF (.NOT.ALLOCATED(srcOffset)) THEN
      ALLOCATE(srcOffset(me%nSrcPE, 3))
      srcOffset(:,:) = 0_i8
    END IF
    IF (ALLOCATED(vSize)) DEALLOCATE(vSize)
    ALLOCATE(vSize(nV))
    vSize(:) = 0
    DO iV_o = 1, nV
      iV = vWrNow(iV_o)
      stride(iV_o) = SUM(me%idx(me%vGrid(iV))%nSrcRecv(:))
      nL = MERGE(me%vLevs(iV), sub(2) - sub(1), sub(3) .NE. iV)
      vSize(iV_o) = stride(iV_o) * nL
#ifndef NOMPI
      pct(:,iV_o) = me%idx(me%vGrid(iV))%nSrcRecv(:)
      rDsp(iV_o) = INT(bSize, addr)
      t_cnt(:) = t_cnt(:) + pct(:,iV_o) * nL
      bSize = bSize + stride(iV_o) * nL
#endif
    END DO
    IF (me%shortcut) THEN
      cBuf%d => me%sBuf%d
      cBuf%s => me%sBuf%s
      cBuf%i => me%sBuf%i
#ifndef NOMPI
    ELSE
      fDsp = facTtoSp(typeMap(iTy))
      dTy = mpiDtype(typeMap(iTy))
      tBy = typeByte(typeMap(iTy))
      bSizeB = INT(bSize, addr) * tBy
      IF (me%rBuf_size .LT. bSizeB) THEN
        IF (me%rBuf_size .GT. 0_addr) THEN
          CALL MPI_Free_mem(me%rBuf%d, ierr)
          HANDLE_MPI_ERROR(ierr, 'MPI_Free_mem')
        END IF
        me%rBuf_size = MAX(bSizeB, MAXVAL(typeByte))
        CALL MPI_Alloc_mem(me%rBuf_size,  MPI_INFO_NULL, cptr, ierr)
        HANDLE_MPI_ERROR(ierr, 'MPI_Alloc_mem')
        CALL C_F_POINTER(cptr, me%rBuf%d, [me%rBuf_size/typeByte(typeMap(REAL_T))])
        CALL C_F_POINTER(cptr, me%rBuf%s, [me%rBuf_size/typeByte(typeMap(SINGLE_T))])
        CALL C_F_POINTER(cptr, me%rBuf%i, [me%rBuf_size/typeByte(typeMap(INT_T))])
      END IF
      n_oReq = 0
      DO i = 1, me%nSrcPE
        o_dTy = MPI_DATATYPE_NULL
        IF (t_cnt(i) > 0) THEN
          DO iV_o = 1, nV
            iV = vWrNow(iV_o)
            nL = MERGE(me%vLevs(iV), sub(2) - sub(1), sub(3) .NE. iV)
            CALL MPI_Type_vector(nL, pct(i,iV_o), stride(iV_o), dTy, tmp_ddt(iV_o), ierr)
            HANDLE_MPI_ERROR(ierr, 'MPI_Type_vector')
            rD(iV_o) = rDsp(iV_o) * tBy
          END DO
          CALL MPI_Type_create_struct(nV, ones, rD, tmp_ddt, o_dTy, ierr)
          HANDLE_MPI_ERROR(ierr, 'MPI_Type_create_struct')
          DO iV_o = 1, nV
            CALL MPI_Type_free(tmp_ddt(iV_o), ierr)
            HANDLE_MPI_ERROR(ierr, 'MPI_Type_free')
          END DO
          CALL MPI_Type_commit(o_dTy, ierr)
          t_dsp = srcOffset(i, typeMap(iTy)) * fDsp + me%tOffSv(4*(i-1)+typeMap(iTy))
          IF (n_oReq .LT. nOpenReqMax) THEN
            n_oReq = n_oReq + 1
            iReq = n_oReq
          ELSE
#ifndef NO_MPI_RGET
            CALL MPI_Waitany(n_oReq, getReq, iReq, MPI_STATUS_IGNORE, ierr)
            HANDLE_MPI_ERROR(ierr, 'MPI_Waitany')
#else
            CALL MPI_Win_unlock(me%srcPE(i-nOpenReqMax), me%win, ierr)
            HANDLE_MPI_ERROR(ierr, 'MPI_Win_unlock')
#endif
          END IF
#ifdef NO_MPI_RGET
          CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, me%srcPE(i), MPI_MODE_NOCHECK, me%win, ierr)
          HANDLE_MPI_ERROR(ierr, 'MPI_Win_lock')
#endif
          SELECT CASE(iTy)
          CASE(REAL_T)
#ifndef NO_MPI_RGET
            CALL MPI_Rget(me%rBuf%d(1:bSize), 1, o_dTy, me%srcPE(i), t_dsp, &
              &           t_cnt(i), dTy, me%win, getReq(iReq), ierr)
#else
            CALL MPI_Get(me%rBuf%d(1:bSize), 1, o_dTy, me%srcPE(i), t_dsp, &
              &          t_cnt(i), dTy, me%win, ierr)
#endif
          CASE(SINGLE_T)
#ifndef NO_MPI_RGET
            CALL MPI_Rget(me%rBuf%s(1:bSize), 1, o_dTy, me%srcPE(i), t_dsp, &
              &           t_cnt(i), dTy, me%win, getReq(iReq), ierr)
#else
            CALL MPI_Get(me%rBuf%s(1:bSize), 1, o_dTy, me%srcPE(i), t_dsp, &
              &          t_cnt(i), dTy, me%win, ierr)
#endif
          CASE(INT_T)
#ifndef NO_MPI_RGET
            CALL MPI_Rget(me%rBuf%i(1:bSize), 1, o_dTy, me%srcPE(i), t_dsp, &
              &           t_cnt(i), dTy, me%win, getReq(iReq), ierr)
#else
            CALL MPI_Get(me%rBuf%i(1:bSize), 1, o_dTy, me%srcPE(i), t_dsp, &
              &          t_cnt(i), dTy, me%win, ierr)
#endif
           END SELECT
#ifndef NO_MPI_RGET
          HANDLE_MPI_ERROR(ierr, 'MPI_Rget')
#else
          HANDLE_MPI_ERROR(ierr, 'MPI_Get')
#endif
          CALL MPI_Type_free(o_dTy, ierr)
          rDsp(:) = rDsp(:) + pct(i,:)
        END IF
      END DO
#ifndef NO_MPI_RGET
      CALL MPI_Waitall(n_oReq, getReq ,MPI_STATUSES_IGNORE, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Waitall')
#else
      DO i = 1, n_oReq
        CALL MPI_Win_unlock(me%srcPE(me%nSrcPE - i + 1), me%win, ierr)
        HANDLE_MPI_ERROR(ierr, 'MPI_Win_unlock')
      END DO
#endif
      srcOffset(:, typeMap(iTy)) = srcOffset(:, typeMap(iTy)) + INT(t_cnt(:), i8)
      cBuf%d => me%rBuf%d
      cBuf%s => me%rBuf%s
      cBuf%i => me%rBuf%i
#endif
    END IF
  END SUBROUTINE multifileRestartCollector_fetch

  SUBROUTINE multifileRestartCollector_sendField(me, iV, input, ioffset)
    CLASS(t_multifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,     INTENT(IN)    :: iV
    CLASS(*),    INTENT(IN)    :: input(:,:,:)
    INTEGER(i8), INTENT(INOUT) :: ioffset
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_sendField"
    INTEGER(i8) :: offset
    INTEGER :: iG, iLev, i, nPnt

    iG = me%vGrid(iV)
    nPnt = me%idx(iG)%nSend
    IF (me%idx(ig)%nSend > 0) THEN
      SELECT TYPE(input)
      TYPE IS (REAL(dp))
        IF (.NOT. ASSOCIATED(me%sBuf%d)) CALL finish(routine, "no send buffer (dp)!")
!ICON_OMP PARALLEL DO COLLAPSE(2) PRIVATE(offset)
        DO iLev = 1, me%vLevs(iV)
          DO i = 1, nPnt
            offset = INT(i, i8) + (INT(iLev, i8) - 1_i8) * INT(nPnt, i8) + ioffset
            me%sBuf%d(offset) = input(me%idx(iG)%sIdx(i), iLev, me%idx(iG)%sBlk(i))
          END DO
        END DO
      TYPE IS (REAL(sp))
        IF (.NOT. ASSOCIATED(me%sBuf%s)) CALL finish(routine, "no send buffer (sp)!")
!ICON_OMP PARALLEL DO COLLAPSE(2) PRIVATE(offset)
        DO iLev = 1, me%vLevs(iV)
          DO i = 1, nPnt
            offset = INT(i, i8) + (INT(iLev, i8) - 1_i8) * INT(nPnt, i8) + ioffset
            me%sBuf%s(offset) = input(me%idx(iG)%sIdx(i), iLev, me%idx(iG)%sBlk(i))
          END DO
        END DO
      TYPE IS (INTEGER)
        IF (.NOT. ASSOCIATED(me%sBuf%i)) CALL finish(routine, "no send buffer (int)!")
!ICON_OMP PARALLEL DO COLLAPSE(2) PRIVATE(offset)
        DO iLev = 1, me%vLevs(iV)
          DO i = 1, nPnt
            offset = INT(i, i8) + (INT(iLev, i8) - 1_i8) * INT(nPnt, i8) + ioffset
            me%sBuf%i(offset) = input(me%idx(iG)%sIdx(i), iLev, me%idx(iG)%sBlk(i))
          END DO
        END DO
      CLASS DEFAULT
        CALL finish(routine, "unrecognized datatype!")
      END SELECT
      ioffset = ioffset + INT(nPnt, i8) * INT(me%vLevs(iV), i8)
    END IF
  END SUBROUTINE multifileRestartCollector_sendField

  SUBROUTINE multifileRestartCollector_init_win(this, isize, shortcut)
    CLASS(t_multifileRestartCollector), INTENT(INOUT) :: this
    INTEGER(i8), INTENT(IN)    :: isize(typeMax)
    LOGICAL, INTENT(IN) :: shortcut
    INTEGER(addr) :: tOffCl(4), wSizes(3)
    CHARACTER(*), PARAMETER :: routine = modname//":CollectorSendBuffer_init_win"
    INTEGER :: i
#ifndef NOMPI
    INTEGER :: ierr, splitKey, wComm_size, j, ii, myRank, &
      &        p_addr, addrBytes
    INTEGER, ALLOCATABLE, DIMENSION(:) :: rank_map
    INTEGER(addr) :: memBytes
    INTEGER(KIND=addr), ALLOCATABLE :: tmpOffSv(:)
    TYPE(c_ptr) :: cMemPtr
#endif

    this%shortcut = shortcut
    NULLIFY(this%wptr%p)
    FORALL(i = 1:3) wSizes(i) = isize(typeID(i))
    tOffCl(:) = 0_addr
    IF (my_process_is_work()) THEN
      tOffCl(2) = tOffCl(1) + wSizes(1) * facTtoSp(1)
      tOffCl(3) = tOffCl(2) + wSizes(2) * facTtoSp(2)
      tOffCl(4) = tOffCl(3) + wSizes(3) * facTtoSp(3)
    END IF
    IF (shortcut) THEN
      ALLOCATE(this%tOffSv(4))
      this%srcPE(1) = 0
      this%tOffSv(:) = tOffCl(:)
#ifndef NOMPI
! this makes "ill" configs work; i.e. using more than 1 writer per worker, but less than 2 -- on average
      CALL MPI_Comm_split(p_comm_work_restart, MPI_UNDEFINED, 0, this%wComm, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Comm_split')
    ELSE
      myRank = p_comm_rank(p_comm_work_restart)
      splitKey = MERGE(1, myRank + 2, iAmRestartWriter())
      CALL MPI_Comm_split(p_comm_work_restart, this%destPE, splitKey, this%wComm, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Comm_split')
      CALL MPI_Comm_size(this%wComm, wComm_size, ierr)
      IF (iAmRestartWriter()) THEN
        ALLOCATE(tmpOffSv(4*wComm_size), rank_map(wComm_size))
      ELSE
        ALLOCATE(tmpOffSv(1), rank_map(1))
      END IF
      CALL MPI_Gather(myRank, 1, p_int, rank_map, 1, p_int, 0, this%wComm, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Gather')
      CALL MPI_Sizeof(1_addr, addrBytes, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Sizeof')
      CALL MPI_Type_match_size(MPI_TYPECLASS_INTEGER, addrBytes, p_addr, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Type_match_size')
      CALL MPI_Gather(tOffCl, 4, p_addr, tmpOffSv, 4, p_addr, 0, this%wComm, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Gather')
      IF (iAmRestartWriter()) THEN
        ALLOCATE(this%tOffSv(4 * this%nSrcPE))
        DO i = 1, this%nSrcPE
          DO j = 1, wComm_size
            IF (rank_map(j) .EQ. this%srcPE(i)) ii = j - 1
          END DO
          this%srcPE(i) = ii
          this%tOffSv((i-1)*4+1:i*4) = tmpOffSv(ii*4+1:(ii+1)*4)
        END DO
      END IF
      DEALLOCATE(rank_map, tmpOffSv)
#else
    ELSE
      CALL finish(routine, "you screwed up!")
#endif
    END IF
    this%destPE = 0
#ifndef NOMPI
    memBytes = tOffCl(4) * typeByte(2)
! as of MPI3.0 standard the following MPI_Alloc_mem interface must be 
! present, if the compiler provides ISO_C_BINDING
    CALL MPI_Alloc_mem(MAX(memBytes, 64_addr), MPI_INFO_NULL, cMemPtr, ierr)
    HANDLE_MPI_ERROR(ierr, 'MPI_Alloc_mem')
    CALL C_F_POINTER(cMemPtr, this%wptr%p, [MAX(1_addr, tOffCl(4))])
    IF (.NOT.shortcut) THEN
      CALL MPI_Win_create(this%wptr%p, memBytes, INT(typeByte(2)), MPI_INFO_NULL, &
        & this%wComm, this%win, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Win_create')
      IF (my_process_is_work()) THEN ! begin with local access
        CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_comm_rank(this%wComm), &
          & MPI_MODE_NOCHECK, this%win, ierr)
        HANDLE_MPI_ERROR(ierr, 'mpi_win_lock')
      END IF
    END IF
#else
    ALLOCATE(this%wptr%p(MAX(1, INT(tOffCl(4)))))
#endif
    tOffCl(1:3) = tOffCl(1:3) + 1_addr
    IF (wSizes(1) .GT. 0_addr) &
      & CALL C_F_POINTER(C_LOC(this%wptr%p(tOffCl(1))), this%sBuf%d, [wSizes(1)])
    IF (wSizes(2) .GT. 0_addr) &
      & CALL C_F_POINTER(C_LOC(this%wptr%p(tOffCl(2))), this%sBuf%s, [wSizes(2)])
    IF (wSizes(3) .GT. 0_addr) &
      & CALL C_F_POINTER(C_LOC(this%wptr%p(tOffCl(3))), this%sBuf%i, [wSizes(3)])
    this%allocd = .true.
  END SUBROUTINE multifileRestartCollector_init_win

  SUBROUTINE multifileRestartCollector_local_access(this)
    CLASS(t_multifileRestartCollector), INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::multifileRestartCollector_local_access"
#ifdef NOMPI
    IF (.NOT.this%allocd) CALL finish(routine, "there is no buffer allocd to fill!")
#else
    INTEGER :: rrank, ierr

    IF (.NOT.this%allocd) CALL finish(routine, "there is no buffer allocd to fill!")
    IF (this%shortcut) RETURN
    IF (this%wStarted) THEN
#ifndef NO_MPI_RGET
      CALL MPI_Win_unlock_all(this%win, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Win_unlock_all')
#endif
      CALL p_barrier(comm=this%wComm)
    END IF
    IF (this%wPosted) THEN
      IF (.NOT.this%wStarted) CALL p_barrier(comm=this%wComm)
      rrank = p_comm_rank(this%wComm)
      CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rrank, &
        & MPI_MODE_NOCHECK, this%win, ierr)
      HANDLE_MPI_ERROR(ierr, 'mpi_win_lock')
    END IF
    this%wPosted = .false.
    this%wStarted = .false.
#endif
  END SUBROUTINE multifileRestartCollector_local_access

  SUBROUTINE multifileRestartCollector_remote_access(this)
    CLASS(t_multifileRestartCollector), INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::multifileRestartCollector_remote_access"
#ifdef NOMPI
    IF (.NOT.this%allocd) CALL finish(routine, "there is no window to expose!")
#else
    INTEGER :: rrank, ierr

    IF (.NOT.this%allocd) CALL finish(routine, "there is no window to expose!")
    IF (this%shortcut) RETURN
    rrank = -1
    IF (my_process_is_work() .AND. .NOT.this%wPosted) THEN
      rrank = p_comm_rank(this%wComm)
      CALL MPI_Win_unlock(rrank, this%win, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Win_unlock')
      CALL p_barrier(this%wComm)
      this%wPosted = .true.
    END IF
    IF (iAmRestartWriter() .AND. .NOT.this%wStarted) THEN
      IF (rrank .EQ. -1) CALL p_barrier(this%wComm)
#ifndef NO_MPI_RGET
      CALL MPI_Win_lock_all(MPI_MODE_NOCHECK, this%win, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Win_lock_all')
#endif
      this%wStarted = .true.
    END IF
#endif
  END SUBROUTINE multifileRestartCollector_remote_access

END MODULE mo_multifile_restart_collector
