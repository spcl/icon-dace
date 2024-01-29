! An implementation of t_scatterPattern that uses MPI_Scatterv to distribute the data.
!
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

MODULE mo_scatter_pattern_scatterv
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_scatter_pattern_base
    USE mo_kind, ONLY: wp, dp, sp, i8
    USE mo_mpi, ONLY: p_real_dp, p_real_sp, p_int, &
    &                 p_gather, p_gatherv, p_scatterv
    USE mo_parallel_config, ONLY: blk_no, idx_no
    USE mo_exception, ONLY: finish

    IMPLICIT NONE

PUBLIC :: t_scatterPatternScatterV

    TYPE, EXTENDS(t_scatterPattern) :: t_scatterPatternScatterV
        !This global description is only created on root rank.
        INTEGER, ALLOCATABLE :: pointCounts(:), displacements(:)    !The recvcounts/sendcounts & displs arrays for
                                                                !MPI_GATHERV/MPI_SCATTERV. Indexes into pointIndices or other
                                                                !arrays that need to be distributed.
        INTEGER, ALLOCATABLE :: pointIndices(:)    !For each point requested by a process, this lists the global index of the point.
        INTEGER :: pointCount !size of pointIndices
    CONTAINS
        PROCEDURE :: construct       => constructScatterPatternScatterV !< override
        PROCEDURE :: distribute_dp   => distributeDataScatterV_dp       !< override
        PROCEDURE :: distribute_spdp => distributeDataScatterV_spdp     !< override
        PROCEDURE :: distribute_sp   => distributeDataScatterV_sp       !< override
        PROCEDURE :: distribute_int  => distributeDataScatterV_int      !< override
        PROCEDURE :: destruct        => destructScatterPatternScatterV  !< override
    END TYPE

PRIVATE

    CHARACTER(*), PARAMETER :: modname = "mo_grid_distribution_scatterv"
    LOGICAL, PARAMETER :: debugModule = .false.

CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------
    !> constructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE constructScatterPatternScatterV(me, jg, loc_arr_len, glb_index, &
         communicator, root_rank)
        CLASS(t_scatterPatternScatterV), TARGET, INTENT(OUT) :: me
        INTEGER, VALUE :: jg, loc_arr_len, communicator
        INTEGER, INTENT(IN) :: glb_index(:)
        INTEGER, OPTIONAL, INTENT(in) :: root_rank

        CHARACTER(*), PARAMETER :: routine = modname//":constructScatterPatternScatterV"
        INTEGER :: comm_size, i, ierr, asize, accum
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%root_rank == me%rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine

        CALL constructScatterPattern(me, jg, loc_arr_len, glb_index, &
             communicator, root_rank)
        comm_size = me%comm_size
        asize = MERGE(comm_size, 1, me%rank == me%root_rank)
        ALLOCATE(me%pointCounts(asize), &
          &      me%displacements(asize), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        CALL p_gather(me%myPointCount, me%pointCounts, me%root_rank, communicator)
        IF (me%rank == me%root_rank) THEN
          accum = 0
          DO i = 1, comm_size
            me%displacements(i) = accum
            accum = accum + me%pointCounts(i)
          END DO
          me%pointCount = accum
          asize = accum
        ELSE
          asize = 1
        END IF
        ALLOCATE(me%pointIndices(asize), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        CALL p_gatherv(glb_index, loc_arr_len, me%pointIndices, &
             me%pointCounts, me%displacements, me%root_rank, communicator)

        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE constructScatterPatternScatterV

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_dp()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_dp(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        REAL(dp), INTENT(INOUT) :: globalArray(:)
        REAL(wp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_dp"
        REAL(dp), ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr, asize
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine

        CALL me%startDistribution()

        asize = MERGE(me%pointCount, 1, me%rank == me%root_rank)
        ALLOCATE(sendArray(asize), &
          &      recvArray(me%myPointCount), stat = ierr)
        IF (ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF (me%rank == me%root_rank) THEN
          DO i = 1, me%pointCount
            sendArray(i) = globalArray(me%pointIndices(i))
          END DO
        END IF

        CALL p_scatterv(sendArray, me%pointCounts, me%displacements, &
          recvArray, me%myPointCount, me%root_rank, me%communicator)

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray, sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 8_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_dp

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_spdp()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_spdp(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        REAL(sp), INTENT(INOUT) :: globalArray(:)
        REAL(wp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_spdp"
        REAL(sp), ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr, asize
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine

        CALL me%startDistribution()

        asize = MERGE(me%pointCount, 1, me%rank == me%root_rank)
        ALLOCATE(sendArray(asize), &
          &      recvArray(me%myPointCount), stat = ierr)
        IF (ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF (me%rank == me%root_rank) THEN
          DO i = 1, me%pointCount
            sendArray(i) = globalArray(me%pointIndices(i))
          END DO
        END IF

        CALL p_scatterv(sendArray, me%pointCounts, me%displacements, &
          recvArray, me%myPointCount, me%root_rank, me%communicator)

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray, sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 4_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_spdp


    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_sp()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_sp(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        REAL(sp), INTENT(INOUT) :: globalArray(:)
        REAL(sp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_sp"
        REAL(sp), ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr, asize
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine

        CALL me%startDistribution()

        asize = MERGE(me%pointCount, 1, me%rank == me%root_rank)
        ALLOCATE(sendArray(asize), &
          &      recvArray(me%myPointCount), stat = ierr)
        IF (ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF (me%rank == me%root_rank) THEN
          DO i = 1, me%pointCount
            sendArray(i) = globalArray(me%pointIndices(i))
          END DO
        END IF

        CALL p_scatterv(sendArray, me%pointCounts, me%displacements, &
          recvArray, me%myPointCount, me%root_rank, me%communicator)

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray, sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 4_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_sp

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_int()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_int(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        INTEGER, INTENT(INOUT) :: globalArray(:)
        INTEGER, INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_sp"
        INTEGER, ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr, asize
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine

        CALL me%startDistribution()

        asize = MERGE(me%pointCount, 1, me%rank == me%root_rank)
        ALLOCATE(sendArray(asize), &
          &      recvArray(me%myPointCount), stat = ierr)
        IF (ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF (me%rank == me%root_rank) THEN
          DO i = 1, me%pointCount
            sendArray(i) = globalArray(me%pointIndices(i))
          END DO
        END IF

        CALL p_scatterv(sendArray, me%pointCounts, me%displacements, &
          recvArray, me%myPointCount, me%root_rank, me%communicator)

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray, sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 4_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_int

    !-------------------------------------------------------------------------------------------------------------------------------
    !> destructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE destructScatterPatternScatterV(me)
        CLASS(t_scatterPatternScatterV), TARGET, INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":destructScatterPatternScatterV"
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine
        DEALLOCATE(me%pointCounts, me%displacements, me%pointIndices)
        CALL destructScatterPattern(me)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE destructScatterPatternScatterV

END MODULE
