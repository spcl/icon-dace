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

! An implementation of t_scatterPattern that uses MPI_Scatter to distribute the data.

MODULE mo_scatter_pattern_scatter
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: wp, dp, sp, i8
    USE mo_scatter_pattern_base
    USE mo_mpi, ONLY: my_process_is_stdio, &
    &                 p_max, p_gather, p_scatter
    USE mo_parallel_config, ONLY: blk_no, idx_no
    USE mo_exception, ONLY: finish

    IMPLICIT NONE

PUBLIC :: t_scatterPatternScatter

    TYPE, EXTENDS(t_scatterPattern) :: t_scatterPatternScatter
        INTEGER :: slapSize    !The count of points sent to each pe, calculated as the maximum of all myPointCount members.
        !> This global description is only created on the root rank.
        !! For each point requested by a process, this lists the
        !! global index of the point.
        INTEGER, ALLOCATABLE :: pointIndices(:,:)
        !> actual length of useful data in pointIndices per rank
        INTEGER, ALLOCATABLE :: point_counts(:)
    CONTAINS
        PROCEDURE :: construct       => constructScatterPatternScatter !< override
        PROCEDURE :: distribute_dp   => distributeDataScatter_dp       !< override
        PROCEDURE :: distribute_spdp => distributeDataScatter_spdp     !< override
        PROCEDURE :: distribute_sp   => distributeDataScatter_sp       !< override
        PROCEDURE :: distribute_int  => distributeDataScatter_int      !< override
        PROCEDURE :: destruct        => destructScatterPatternScatter  !< override
    END TYPE

PRIVATE

    CHARACTER(*), PARAMETER :: modname = "mo_grid_distribution_scatter"
    LOGICAL, PARAMETER :: debugModule = .false.

CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------
    !> constructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE constructScatterPatternScatter(me, jg, loc_arr_len, glb_index, &
         communicator, root_rank)
        CLASS(t_scatterPatternScatter), TARGET, INTENT(OUT) :: me
        INTEGER, VALUE :: jg, loc_arr_len, communicator
        INTEGER, INTENT(IN) :: glb_index(:)
        INTEGER, OPTIONAL, INTENT(in) :: root_rank

        CHARACTER(*), PARAMETER :: routine &
             = modname//":costructScatterPatternScatter"
        INTEGER :: ierr, pt_shape(2)
        INTEGER, ALLOCATABLE :: myIndices(:)
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine
        CALL constructScatterPattern(me, jg, loc_arr_len, glb_index, communicator, root_rank)
        me%slapSize = p_max(me%myPointCount, comm = communicator)
        pt_shape(1) = MERGE(me%slapSize, 1, me%rank == me%root_rank)
        pt_shape(2) = MERGE(me%comm_size, 1, me%rank == me%root_rank)
        ALLOCATE(me%pointIndices(pt_shape(1), pt_shape(2)), &
                 me%point_counts(pt_shape(2)), myIndices(me%slapSize), &
                 stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        myIndices(1:me%myPointCount) = glb_index
        myIndices(me%myPointCount+1:me%slapSize) = -1
        CALL p_gather(me%myPointCount, me%point_counts, me%root_rank, &
             communicator)
        CALL p_gather(myIndices, me%pointIndices, me%root_rank, communicator)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE constructScatterPatternScatter

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_dp
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatter_dp(me, globalArray, localArray, ladd_value)
        CLASS(t_scatterPatternScatter), INTENT(INOUT) :: me
        REAL(dp), INTENT(INOUT) :: globalArray(:)
        REAL(wp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine &
             = modname//":distributeDataScatter_dp"
        REAL(dp), ALLOCATABLE :: sendArray(:,:), recvArray(:)
        INTEGER :: i, j, blk, idx, ierr, send_shape(2)
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine
        CALL me%startDistribution()

        send_shape = SHAPE(me%pointIndices)
        ALLOCATE(sendArray(send_shape(1), send_shape(2)), &
             recvArray(me%slapSize), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF (me%rank == me%root_rank) THEN
!$OMP PARALLEL DO PRIVATE(i,j)
          DO j = 1, me%comm_size
            DO i = 1, me%point_counts(j)
              sendArray(i, j) = globalArray(me%pointIndices(i, j))
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
        CALL p_scatter(sendArray, recvArray, me%root_rank, me%communicator)
        IF(ladd_value) THEN
!$NEC ivdep
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
        CALL me%endDistribution(INT(send_shape(1), i8) * INT(send_shape(2), i8) * 8_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatter_dp

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_spdp
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatter_spdp(me, globalArray, localArray, ladd_value)
        CLASS(t_scatterPatternScatter), INTENT(INOUT) :: me
        REAL(sp), INTENT(INOUT) :: globalArray(:)
        REAL(wp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatter_spdp"
        REAL(sp), ALLOCATABLE :: sendArray(:,:), recvArray(:)
        INTEGER :: i, j, blk, idx, ierr, send_shape(2)
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine

        CALL me%startDistribution()

        send_shape = SHAPE(me%pointIndices)
        ALLOCATE(sendArray(send_shape(1), send_shape(2)), &
             recvArray(me%slapSize), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF (me%rank == me%root_rank) THEN
!$OMP PARALLEL DO PRIVATE(i,j)
          DO j = 1, me%comm_size
            DO i = 1, me%point_counts(j)
              sendArray(i, j) = globalArray(me%pointIndices(i, j))
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
        CALL p_scatter(sendArray, recvArray, me%root_rank, me%communicator)
        IF(ladd_value) THEN
!$NEC ivdep
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
        CALL me%endDistribution(INT(send_shape(1), i8) * INT(send_shape(2), i8) * 4_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatter_spdp

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_sp
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatter_sp(me, globalArray, localArray, ladd_value)
        CLASS(t_scatterPatternScatter), INTENT(INOUT) :: me
        REAL(sp), INTENT(INOUT) :: globalArray(:)
        REAL(sp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine &
             = modname//":distributeDataScatter_sp"
        REAL(sp), ALLOCATABLE :: sendArray(:,:), recvArray(:)
        INTEGER :: i, j, blk, idx, ierr, send_shape(2)
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine

        CALL me%startDistribution()

        send_shape = SHAPE(me%pointIndices)
        ALLOCATE(sendArray(send_shape(1), send_shape(2)), &
             recvArray(me%slapSize), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF(me%rank == me%root_rank) THEN
!$OMP PARALLEL DO PRIVATE(i,j)
            DO j = 1, me%comm_size
              DO i = 1, me%point_counts(j)
                sendArray(i, j) = globalArray(me%pointIndices(i, j))
              END DO
            END DO
!$OMP END PARALLEL DO
        END IF
        CALL p_scatter(sendArray, recvArray, me%root_rank, me%communicator)
        IF(ladd_value) THEN
!$NEC ivdep
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
        CALL me%endDistribution(INT(send_shape(1), i8) * INT(send_shape(2), i8) * 4_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatter_sp

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_int
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatter_int(me, globalArray, localArray, ladd_value)
        CLASS(t_scatterPatternScatter), INTENT(INOUT) :: me
        INTEGER, INTENT(INOUT) :: globalArray(:)
        INTEGER, INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine &
             = modname//":distributeDataScatter_sp"
        INTEGER, ALLOCATABLE :: sendArray(:,:), recvArray(:)
        INTEGER :: i, j, blk, idx, ierr, send_shape(2)
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine
        CALL me%startDistribution()

        send_shape = SHAPE(me%pointIndices)
        ALLOCATE(sendArray(send_shape(1), send_shape(2)), &
             recvArray(me%slapSize), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        IF (me%rank == me%root_rank) THEN
            DO j = 1, me%comm_size
              DO i = 1, me%point_counts(j)
                sendArray(i, j) = globalArray(me%pointIndices(i, j))
              END DO
            END DO
        END IF
        CALL p_scatter(sendArray, recvArray, me%root_rank, me%communicator)
        IF(ladd_value) THEN
!$NEC ivdep
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
        CALL me%endDistribution(INT(send_shape(1), i8) * INT(send_shape(2), i8) * 4_i8)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatter_int

    !-------------------------------------------------------------------------------------------------------------------------------
    !> destructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE destructScatterPatternScatter(me)
        CLASS(t_scatterPatternScatter), TARGET, INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine &
             = modname//":destructScatterPatternScatter"
        LOGICAL :: l_write_debug_info

        l_write_debug_info = debugmodule .AND. me%rank == me%root_rank

        IF (l_write_debug_info) WRITE(0,*) "entering ", routine
        DEALLOCATE(me%pointIndices)
        CALL destructScatterPattern(me)
        IF (l_write_debug_info) WRITE(0,*) "leaving ", routine
    END SUBROUTINE destructScatterPatternScatter

END MODULE
