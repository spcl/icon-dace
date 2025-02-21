! MPI global reshuffle.
!
! This program works on a distributed-memory array, where each PE's
! part is defined by global indices "owner_idx" (local size!).  The
! "reshuffle" operation means that each PE writes a number of "nsend"
! values (where "nsend" may be different for each PE) to the array at
! global indices "glb_idx". This involves, of course, some
! send/receive operations, since the destination indices may be stored
! on a different PE. This implementation of the "reshuffle" operation
! involves no global-size arrays.
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

MODULE mo_reshuffle




  USE mo_exception,          ONLY: finish
  USE mo_util_sort, ONLY: quicksort



  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: reshuffle

  TYPE t_regular_partition
    INTEGER :: n_indices, n_pes, irank
  CONTAINS
    PROCEDURE :: local_size => regular_partition_local_size
    PROCEDURE :: glb2local  => regular_partition_glb2local
    PROCEDURE :: glbidx2pe  => regular_partition_glbidx2pe
  END TYPE t_regular_partition

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_reshuffle'

CONTAINS

  INTEGER FUNCTION regular_partition_local_size(this)
    CLASS(t_regular_partition), INTENT(IN) :: this
    regular_partition_local_size = ((this%n_indices+this%n_pes-1) / this%n_pes)
  END FUNCTION regular_partition_local_size

  INTEGER FUNCTION regular_partition_glb2local(this,iidx)
    CLASS(t_regular_partition), INTENT(IN) :: this
    INTEGER,                    INTENT(IN) :: iidx
    regular_partition_glb2local = iidx - this%local_size()*this%irank 
  END FUNCTION regular_partition_glb2local

  INTEGER FUNCTION regular_partition_glbidx2pe(this,iidx)
    CLASS(t_regular_partition), INTENT(IN) :: this
    INTEGER,                    INTENT(IN) :: iidx

    regular_partition_glbidx2pe = -1
    IF ((iidx > 0) .AND. (iidx <= this%n_indices) .AND. (this%n_indices >= this%n_pes)) THEN
      regular_partition_glbidx2pe = (iidx-1)/this%local_size()
    END IF
  END FUNCTION regular_partition_glbidx2pe

  SUBROUTINE swap_int(a, i,j, permutation)
    INTEGER,  INTENT(INOUT)           :: a(:)           !< array for in-situ sorting
    INTEGER,  INTENT(IN)              :: i,j            !< indices to be exchanged
    INTEGER,  INTENT(INOUT), OPTIONAL :: permutation(:) !< (optional) permutation of indices
    ! local variables
    INTEGER :: t, t_p

    t    = a(i)
    a(i) = a(j)
    a(j) = t
    IF (PRESENT(permutation)) THEN
      t_p            = permutation(i)
      permutation(i) = permutation(j)
      permutation(j) = t_p
    END IF
  END SUBROUTINE swap_int

  SUBROUTINE calc_displs(i_pe, displs)
    INTEGER, INTENT(IN)    :: i_pe(:)    !< i_pe(i) contains the PE where to send/receive index "i" to/from
    INTEGER, INTENT(INOUT) :: displs(0:) !< displs(p) contains the start offset for PE "p"
    INTEGER :: i, counts(0:(SIZE(displs)-1))

    counts(:) = 0
    DO i=1,SIZE(i_pe)
      counts(i_pe(i)) = counts(i_pe(i)) + 1
    END DO
    displs(:) = 0
    DO i=1,(SIZE(displs)-1)
      displs(i) = displs(i-1) + counts(i-1)
    END DO
  END SUBROUTINE calc_displs

  ! -----------------------------------------------------------------------
  ! This subroutine works on a distributed-memory array, where each
  ! PE's part is defined by global indices "owner_idx" (local size!).
  ! The "reshuffle" operation means that each PE writes a number of
  ! values (which may be different for each PE) to the array at global
  ! indices "glb_idx". This involves, of course, some send/receive
  ! operations, since the destination indices may be stored on a
  ! different PE. This implementation of the "reshuffle" operation
  ! involves no global-size arrays.
  !
  ! This routine allows to receive multiple entries for one global
  ! destination index. In this case, these are counted and at most
  ! "ncollisions" distinct values are stored.
  !
  ! MPI_ALLTOALL operations needed:
  ! 1x MPI_ALLTOALL : send/receive counts
  ! 1x MPI_ALLTOALLV: send/receive indices
  ! 2x MPI_ALLTOALLV: send/receive values
  !
  ! IN:  glb_idx(1,...,nsend)    : global indices
  !      values(1,...,nsend)     : values to send
  !      owner_idx(1,...,nlocal) : global indices owned by this PE
  ! OUT: out_values(1,...,ncollisions,1,...,nlocal): buffer for received values
  !      out_count(1,...,ncollisions,1,...,nlocal) : number of received duplicates
  !
  ! Note: Only those entries in out_values are modified which correspond to 
  !       entries in "in_glb_idx" on some PE.
  ! -----------------------------------------------------------------------
  SUBROUTINE reshuffle(description, in_glb_idx, in_values, owner_idx, nglb_indices, communicator, out_values, &
    &                  out_count)
    CHARACTER(LEN=*), INTENT(IN) :: description !< description string (for debugging purposes)
    INTEGER, INTENT(IN), CONTIGUOUS :: in_glb_idx(:)    !< global indices to which "values" correspond
    INTEGER, INTENT(IN), CONTIGUOUS :: in_values(:)     !< values to send
    INTEGER, INTENT(IN), CONTIGUOUS :: owner_idx(:)     !< array indices "owned" by local PE
    INTEGER, INTENT(IN)    :: nglb_indices     !< total size of distributed array
    INTEGER, INTENT(IN)    :: communicator     !< MPI comm.
    INTEGER, INTENT(INOUT), CONTIGUOUS :: out_values(:,:)  !< resulting local part of distributed array
    INTEGER, INTENT(INOUT), CONTIGUOUS :: out_count(:,:)   !< counts, how often an entry was received

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':reshuffle'
    INTEGER                   :: nsend, nlocal, i, j, ncollisions, local_idx, nvals, dst_idx, nerror
    LOGICAL                   :: lfound

    LOGICAL                   :: lcollision

    INTEGER, ALLOCATABLE      :: reg_partition_buf(:,:), reg_partition_modified(:),        &
      &                          reg_partition_count(:,:)

    nlocal = SIZE(owner_idx)
    nsend  = SIZE(in_glb_idx)
    ncollisions = SIZE(out_values,1) ! max no. of concurrent sets

    ! consistency checks:
    IF (SIZE(out_values,2) /= nlocal)  CALL finish(routine, TRIM(description)//" - SIZE(out_values) /= nlocal")
    IF (SIZE(in_values) /= nsend)  CALL finish(routine, TRIM(description)//" - SIZE(in_values) /= nsend")
    IF (MAXVAL(in_glb_idx) > nglb_indices)  THEN
      WRITE (0,*) "MAXVAL(in_glb_idx) = ", MAXVAL(in_glb_idx), "; nglb_indices = ", nglb_indices
      CALL finish(routine, TRIM(description)//" - MAXVAL(in_glb_idx) > nglb_indices")
    END IF
    IF (MINVAL(in_glb_idx) < 1) THEN
      CALL finish(routine, TRIM(description)//" - MINVAL(in_glb_idx) < 1")
    END IF
    IF (MAXVAL(owner_idx) > nglb_indices) THEN
      CALL finish(routine, TRIM(description)//" - MAXVAL(owner_idx) > nglb_indices")
    END IF
    IF (MINVAL(owner_idx) < 1) THEN
      CALL finish(routine, TRIM(description)//" - MINVAL(owner_idx) < 1")
    END IF
    IF ((SIZE(out_values,1) /= SIZE(out_count,1)) .OR. &
      & (SIZE(out_values,2) /= SIZE(out_count,2))) THEN
      CALL finish(routine, TRIM(description)//" - Size of out_count/out_data mismatch")
    END IF


    ! non-MPI mode: local copy
    ALLOCATE(reg_partition_buf(ncollisions, nglb_indices),    &
      &      reg_partition_modified(nglb_indices),            &
      &      reg_partition_count(ncollisions, nglb_indices))
    reg_partition_modified(:) = 0
    reg_partition_buf(:,:)    = 0
    reg_partition_count = 0
    lcollision = .FALSE.
    SEND_LOOP: DO i=1,nsend
      local_idx = in_glb_idx(i)
      ! we try to avoid collisions, when one entry is repeatedly
      ! modified:
      nvals  = reg_partition_modified(local_idx)
      dst_idx = in_values(i)
      LOOP_FOUND : DO j=1,nvals
        IF (reg_partition_buf(j, local_idx) == dst_idx) THEN
          reg_partition_count(j, local_idx) = reg_partition_count(j, local_idx) + 1
          CYCLE SEND_LOOP
        END IF
      END DO LOOP_FOUND
      lcollision = lcollision .OR. nvals >= ncollisions
      nvals = nvals + MERGE(0, 1, lcollision)
      reg_partition_buf(nvals, local_idx)   = dst_idx
      reg_partition_count(nvals, local_idx) = 1
      reg_partition_modified(local_idx)     = nvals
    END DO SEND_LOOP
    IF (lcollision)  CALL finish(routine, TRIM(description)//" - Error! Too many collisions!")
    DO i=1,nlocal
      dst_idx = owner_idx(i)
      IF (reg_partition_modified(dst_idx) > 0) THEN
        out_values(:,i) = reg_partition_buf(:,dst_idx)
        out_count(:,i) = reg_partition_count(:,dst_idx)
      ELSE
        out_count(:,i) = 0
      END IF
    END DO

    ! non-MPI mode: clean up
    DEALLOCATE(reg_partition_modified)

  END SUBROUTINE reshuffle

END MODULE mo_reshuffle
