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

MODULE mo_multifile_restart_collector

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
  USE mo_kind, ONLY: addr => i8
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

  INTEGER, PARAMETER :: MPI_COMM_NULL = 0, MPI_WIN_NULL = 0

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

    IF (ALLOCATED(this%srcPE)) DEALLOCATE(this%srcPE)
    DO i = 1, 3
      IF (ALLOCATED(this%idx(i)%nSrcRecv)) &
        & DEALLOCATE(this%idx(i)%nSrcRecv, this%glb_idx(i)%p)
      IF (ALLOCATED(this%idx(i)%sIdx)) DEALLOCATE(this%idx(i)%sIdx, this%idx(i)%sBlk)
    END DO
    IF (ALLOCATED(this%vGrid)) DEALLOCATE(this%vGrid, this%vType, this%vLevs)
    IF (this%allocd) THEN
      DEALLOCATE(this%wptr%p)
      NULLIFY(this%wptr%p, this%sBuf%d, this%sBuf%s, this%sBuf%i)
      this%allocd = .false.
    END IF
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
    END DO
    IF (me%shortcut) THEN
      cBuf%d => me%sBuf%d
      cBuf%s => me%sBuf%s
      cBuf%i => me%sBuf%i
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
    ELSE
      CALL finish(routine, "you screwed up!")
    END IF
    this%destPE = 0
    ALLOCATE(this%wptr%p(MAX(1, INT(tOffCl(4)))))
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
    IF (.NOT.this%allocd) CALL finish(routine, "there is no buffer allocd to fill!")
  END SUBROUTINE multifileRestartCollector_local_access

  SUBROUTINE multifileRestartCollector_remote_access(this)
    CLASS(t_multifileRestartCollector), INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::multifileRestartCollector_remote_access"
    IF (.NOT.this%allocd) CALL finish(routine, "there is no window to expose!")
  END SUBROUTINE multifileRestartCollector_remote_access

END MODULE mo_multifile_restart_collector
