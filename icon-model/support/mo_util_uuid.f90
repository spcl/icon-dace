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

! This module contains the routines for the calculation of 128bit
! data checksums/fingerprints (which are called "UUID's" in ICON).
!
! The larger part of this functionality is implemented in C routines
! in the "support" subdirectory, while this Fortran module acts
! merely as a wrapper. For *parallel* fingerprint calculation,
! however, the MPI-parallel communication is invoked on the Fortran
! level only.

MODULE mo_util_uuid

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_SIGNED_CHAR, C_NULL_CHAR, &
    &                                    C_DOUBLE, C_INT, C_PTR, C_F_POINTER, C_LOC
  USE mo_util_uuid_types, ONLY: t_uuid
  USE mo_util_sort, ONLY: quicksort
#ifdef __SX__
  USE mo_util_sort,      ONLY: radixsort_int
#endif
#ifndef NOMPI
  USE MPI
#endif

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: uuid_generate
  PUBLIC :: compare_uuid
  PUBLIC :: uuid_parse
  PUBLIC :: uuid_unparse
  PUBLIC :: uuid2char
  PUBLIC :: char2uuid
  PUBLIC :: clear_uuid
  PUBLIC :: OPERATOR(==)
  PUBLIC :: UUID_EQUAL, UUID_EQUAL_LIMITED_ACCURACY, UUID_UNEQUAL

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_uuid'

  TYPE, BIND(C) :: t_context
    INTEGER(C_INT) :: f64(2)
    INTEGER(C_INT) :: f48(2)
    INTEGER(C_INT) :: tl64(2)
    INTEGER(C_INT) :: tl48(2)
    INTEGER(C_INT) :: max_zero
  END TYPE t_context

  ENUM, BIND(c)
    ENUMERATOR :: UUID_EQUAL                  = 0
    ENUMERATOR :: UUID_EQUAL_LIMITED_ACCURACY = 1
    ENUMERATOR :: UUID_UNEQUAL                = 2
  END ENUM

  INTERFACE OPERATOR (==)
    MODULE PROCEDURE uuid_compare
  END INTERFACE OPERATOR (==)

  INTERFACE uuid_generate
    MODULE PROCEDURE uuid_generate_sequential
    MODULE PROCEDURE uuid_generate_parallel
    MODULE PROCEDURE uuid_generate_parallel1D
  END INTERFACE uuid_generate

  INTERFACE
    SUBROUTINE my_uuid_generate(val, nval, uuid) BIND(C,NAME='uuid_generate')
      IMPORT :: C_DOUBLE, C_INT, t_uuid
      REAL(c_double),   INTENT(IN)         :: val(*)
      INTEGER(c_int),   INTENT(IN), VALUE  :: nval
      TYPE(t_uuid),     INTENT(OUT)        :: uuid
    END SUBROUTINE my_uuid_generate
  END INTERFACE

  INTERFACE
    INTEGER(C_INT) FUNCTION my_compare_uuid(uuid_A, uuid_B, min_difference) BIND(C,NAME='compare_UUID')
      IMPORT :: t_uuid, C_INT, C_DOUBLE
      TYPE(t_uuid),     INTENT(IN), VALUE  :: uuid_A, uuid_B
      REAL(c_double),   INTENT(OUT)        :: min_difference
    END FUNCTION my_compare_uuid
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_unparse(uuid_string, uuid) BIND(C,NAME='uuid_unparse')
      IMPORT :: C_CHAR, t_uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: uuid_string
      TYPE(t_uuid),                    INTENT(in)  :: uuid
    END SUBROUTINE my_uuid_unparse
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_parse(uuid, uuid_string) BIND(C,NAME='uuid_parse')
      IMPORT :: C_CHAR, t_uuid
      TYPE(t_uuid),                    INTENT(out) :: uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in)  :: uuid_string
    END SUBROUTINE my_uuid_parse
  END INTERFACE

  INTERFACE
    SUBROUTINE my_encode_uuid(fingerprint, uuid) BIND(C,NAME='encode_uuid')
      IMPORT :: t_uuid, C_PTR
      TYPE(C_PTR), VALUE        :: fingerprint
      TYPE(t_uuid), INTENT(out) :: uuid
    END SUBROUTINE my_encode_uuid
  END INTERFACE

  INTERFACE
    FUNCTION my_uuid_scan_data(val, nval) BIND(C,NAME='uuid_scan_data') RESULT(fingerprint)
      IMPORT :: C_PTR, C_DOUBLE, C_INT
      REAL(C_DOUBLE),     INTENT(IN) :: val(*)
      INTEGER(C_INT), INTENT(IN), VALUE :: nval
      TYPE(C_PTR) :: fingerprint
    END FUNCTION my_uuid_scan_data
  END INTERFACE

  INTERFACE
    FUNCTION my_concat_fingerprints(context0, context1) BIND(C,NAME='concat_fingerprints') RESULT(fingerprint)
      IMPORT :: C_PTR
      TYPE(C_PTR), VALUE :: context0
      TYPE(C_PTR), VALUE :: context1
      TYPE(C_PTR) :: fingerprint
    END FUNCTION my_concat_fingerprints
  END INTERFACE

  INTERFACE
    SUBROUTINE deallocate_c(ptr) BIND(C,NAME='deallocate_fingerprint') 
      IMPORT :: C_PTR
      TYPE(C_PTR), VALUE :: ptr
    END SUBROUTINE deallocate_c
  END INTERFACE

CONTAINS
 
  SUBROUTINE uuid_parse(uuid_string, uuid)
    CHARACTER(len=*), INTENT(in)  :: uuid_string
    TYPE(t_uuid),     INTENT(out) :: uuid
    CALL my_uuid_parse(uuid, TRIM(uuid_string)//C_NULL_CHAR)
  END SUBROUTINE uuid_parse

  SUBROUTINE uuid_unparse(uuid, uuid_string)
    TYPE(t_uuid),     INTENT(in)  :: uuid
    CHARACTER(len=*), INTENT(out) :: uuid_string
    CALL my_uuid_unparse(uuid_string, uuid)
  END SUBROUTINE uuid_unparse

  FUNCTION uuid_compare(uuid1, uuid2)
    LOGICAL :: uuid_compare
    TYPE(t_uuid),     INTENT(in)  :: uuid1
    TYPE(t_uuid),     INTENT(in)  :: uuid2
    uuid_compare = .TRUE.
!CDIR NOVECTOR
    IF (ANY(uuid1%data /= uuid2%data)) THEN
      uuid_compare = .FALSE.
    ENDIF
  END FUNCTION uuid_compare
  
  SUBROUTINE uuid2char(uuid, string)
    TYPE(t_uuid), INTENT(in) :: uuid
    CHARACTER(len=1), INTENT(out) :: string(16)
    string = TRANSFER(uuid%data, string)
  END SUBROUTINE uuid2char

  SUBROUTINE char2uuid(string, uuid)
    CHARACTER(len=1), INTENT(in) :: string(16)
    TYPE(t_uuid), INTENT(out) :: uuid
    uuid%data  = TRANSFER(string, uuid%data)
  END SUBROUTINE char2uuid

  SUBROUTINE uuid_generate_sequential(val, uuid)
    REAL(C_DOUBLE), INTENT(IN)  :: val(:)
    TYPE(t_uuid),   INTENT(out) :: uuid
    CALL my_uuid_generate(val, SIZE(val), uuid)
  END SUBROUTINE uuid_generate_sequential

  INTEGER(C_INT) FUNCTION compare_uuid(uuid_A, uuid_B) 
    TYPE(t_uuid),   INTENT(IN)  :: uuid_A, uuid_B    
    REAL(C_DOUBLE) :: min_difference
    compare_uuid = my_compare_uuid(uuid_A, uuid_B, min_difference) 
  END FUNCTION compare_uuid


  ! Parallel computation of UUID.
  !
  ! In contrast to the sequential version "uuid_generate", here the
  ! data may be spread over different PEs in an unordered way. Without
  ! using a global array, the UUID is then generated by concatenating
  ! independently computed fingerprints.
  !
  SUBROUTINE uuid_generate_parallel(comm, in_val, in_glbidx, glb_nval, uuid)
    INTEGER,        INTENT(IN)  :: comm             !< MPI communicator
    REAL(C_DOUBLE), INTENT(IN)  :: in_val(:,:)      !< input value, dim2: data-parallel
    INTEGER,        INTENT(IN)  :: in_glbidx(:)     !< "in_glbidx(:,i)": global index of "in_val(:,i)"
    INTEGER,        INTENT(IN)  :: glb_nval         !< total no. of global indices
    TYPE(t_uuid),   INTENT(OUT) :: uuid             !< (output:) UUID
    ! local variables
    INTEGER, PARAMETER :: ROOTPE    = 0
    INTEGER, PARAMETER :: dbg_level = 0

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':uuid_generate_parallel'
    TYPE(C_PTR)                          :: fingerprint, fingerprint_merge, new_fingerprint
    TYPE(t_context), TARGET, ALLOCATABLE :: fingerprint_i(:)
    INTEGER                              :: nval, nrow, i, p_error, ipe, npes, target_pe,  &
      &                                     nval_local, chunksize, p_c_double_byte, p_c_double
    REAL(C_DOUBLE)                       :: tmp_c_double
    TYPE(t_context), POINTER             :: ptr
    INTEGER, ALLOCATABLE                 :: sendbuf(:), recvbuf(:,:), permutation(:),      &
      &                                     glbidx_sorted(:), glbidx_local(:), glbidx(:),  &
      &                                     sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
    REAL(C_DOUBLE), ALLOCATABLE          :: val_sorted(:,:), val_local(:,:)
    LOGICAL                              :: lcontiguous

    nrow = SIZE(in_val,1)
    nval = SIZE(in_val,2)

    IF (SIZE(in_glbidx) /= nval) THEN
      WRITE (0,*) routine, ": Size of input fields does not match!" ; RETURN
    END IF

#ifndef NOMPI
    CALL MPI_COMM_SIZE(comm, npes, p_error)
    CALL MPI_COMM_RANK(comm, ipe,  p_error)

    chunksize = INT(FLOOR(REAL(glb_nval)/REAL(npes)))

    ! first, we need to reorder the scattered global indices without
    ! the use of a global-size array. We need to build an
    ! MPI_ALLTOALL, where the global index "glbidx" ends up at PE
    ! glbidx/chunksize:

    ! --- build local copies of "val", "glbidx", sorted by "glbidx":
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'build local copies of "val", "glbidx", sorted by "glbidx":'
    END IF
    ALLOCATE(val_sorted(nrow,nval), permutation(nval), glbidx_sorted(nval), &
      &      glbidx(nval))
    glbidx        = in_glbidx
    permutation   = (/ ( i, i=1,nval ) /)
#ifdef __SX__
    CALL radixsort_int(glbidx, permutation)
#else
    CALL quicksort(glbidx, permutation)
#endif
    val_sorted    = in_val(:,permutation(:))
    glbidx_sorted = glbidx
    DEALLOCATE(permutation)

    ! --- build a "displacement array", i.e. the first index which is
    ! --- sent to a specific PE:
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'build a "displacement array"'
    END IF
    ALLOCATE(sdispls(npes), sendcounts(npes))
    sdispls       = -1
    sendcounts(:) =  0
    DO i=1,SIZE(glbidx_sorted)
      ! every PE gets "chunksize" data entries to work on, except the
      ! last one which may have a slightly higher workload.
      target_pe = MIN( (glbidx_sorted(i)-1)/chunksize + 1, npes)
      sendcounts(target_pe) = sendcounts(target_pe) + 1
      IF (sdispls(target_pe) == -1)  sdispls(target_pe) = (i-1)
    END DO
    WHERE (sdispls(:) < 0)  sdispls(:) = 0

    ! --- each PE tells each other PE how many items it indends to
    ! --- send:
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'each PE tells each other PE how many items'
    END IF
    ALLOCATE(recvcounts(npes), rdispls(npes))
    CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, comm, p_error)

    ! --- build a "displacement array" for the receiver side, i.e. the
    ! --- first index where to receive a block from another PE:
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'build a "displacement array" for the receiver side'
    END IF
    rdispls(:) = 0
    DO i=1,npes-1
      rdispls(i+1) = rdispls(i) + recvcounts(i)
    END DO

    ! --- consistency check: test if our PE is going to receive at
    ! --- most "chunksize" items:
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'consistency check'
    END IF
    nval_local = SUM(recvcounts)
    IF ((nval_local > chunksize) .AND. (ipe < (npes-1))) THEN
      WRITE (0,*) routine, ": PE ", ipe, " - Internal error! #local values=", nval_local, &
        &                  "; chunksize=", chunksize ; RETURN
    END IF
    IF (ANY(sendcounts(:) < 0) .OR. ANY(recvcounts(:) < 0) .OR.  &
      & ANY(sdispls(:) < 0)    .OR. ANY(rdispls(:) < 0)) THEN
      WRITE (0,*) routine, ": PE ", ipe, " - Internal error: counts/displacements!" ; RETURN
    END IF

    ! determine MPI data type for REAL(C_DOUBLE):
    CALL MPI_SIZEOF(tmp_c_double, p_c_double_byte, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, p_c_double_byte, p_c_double, p_error)

    ! --- now perform the ALLTOALLV operation for the data and the
    ! --- global indices:
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'perform the ALLTOALLV operation'
    END IF
    ALLOCATE(glbidx_local(nval_local), val_local(nrow, nval_local))
    CALL MPI_ALLTOALLV(glbidx_sorted, sendcounts, sdispls, MPI_INTEGER, &
      &                glbidx_local,  recvcounts, rdispls, MPI_INTEGER, comm, p_error)
    sendcounts = nrow * sendcounts
    sdispls    = nrow * sdispls
    recvcounts = nrow * recvcounts
    rdispls    = nrow * rdispls
    ! consistency check
    IF (ANY(rdispls+recvcounts > nrow*nval_local)) THEN
      WRITE (0,*) routine, ": PE ", ipe, " - Internal error: receive buffer size!" ; RETURN
    END IF

    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'perform the ALLTOALLV operation II'
    END IF
    CALL MPI_ALLTOALLV(val_sorted, sendcounts, sdispls, p_c_double, &
      &                val_local,  recvcounts, rdispls, p_c_double, comm, p_error)

    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'sort the result'
    END IF
    ALLOCATE(permutation(nval_local))
    permutation  = (/ ( i, i=1,nval_local ) /)
#ifdef __SX__
    CALL radixsort_int(glbidx_local, permutation)
#else
    CALL quicksort(glbidx_local, permutation)
#endif
    val_local    = val_local(:,permutation(:))
    DEALLOCATE(permutation)

    ! --- consistency check: test if this PE has received a contiguous
    ! --- chunk of indices:
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'consistency check II'
    END IF
    lcontiguous = .TRUE.
    DO i=2,nval_local
      IF (glbidx_local(i) /= glbidx_local(i-1)+1) THEN
        lcontiguous = .FALSE.
        EXIT
      END IF
    END DO
    IF (.NOT. lcontiguous) THEN
      WRITE (0,*) routine, ": Internal error, PE has non-contiguous chunk of indices!" ; RETURN
    END IF
    ! check if the every PE except the last one has received
    ! "chunksize" items:
    IF ((ipe < (npes-1)) .AND. (nval_local /= chunksize)) THEN
      WRITE (0,*) routine, ": Internal error, PE ", ipe, " has not the right no. of local values: ", &
        &         "nval_local=",nval_local, "; chunksize=", chunksize
      RETURN
    END IF
    ! check if the last PE has got the index "glb_nval":
    IF (ipe == (npes-1)) THEN
      IF (glbidx_local(nval_local) /= glb_nval) THEN
        WRITE (0,*) routine, ': Internal error, last index not "glb_nval"!' 
        WRITE (0,*) "glb_nval = ", glb_nval, "; npes = ", npes
        WRITE (0,*) "Indices are: ", glbidx_local(:)
        RETURN
      END IF
    END IF

    ! each PE creates "its" part as an independent fingerprint:
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'each PE creates "its" part as an independent fingerprint'
    END IF
    fingerprint = my_uuid_scan_data(val_local, nrow*nval_local)

    IF (npes == 1) THEN

      ! this branch especially holds true in "ptestrun" mode:
      CALL my_encode_uuid(fingerprint, uuid)

    ELSE

      ! gather all fingerprints on PE#0:
      IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
        WRITE (0,*) 'gather all fingerprints on PE#0'
      END IF
      CALL C_F_POINTER(fingerprint, ptr)
      ALLOCATE(recvbuf(9,npes), sendbuf(9))
      sendbuf(:) = (/ ptr%f64(1),  ptr%f64(2),  &
        &             ptr%f48(1),  ptr%f48(2),  &
        &             ptr%tl64(1), ptr%tl64(2), &
        &             ptr%tl48(1), ptr%tl48(2), &
        &             ptr%max_zero /)
      CALL MPI_GATHER(sendbuf, SIZE(sendbuf), MPI_INTEGER, &
        &             recvbuf, SIZE(sendbuf), MPI_INTEGER, ROOTPE, comm, p_error)
      
      IF (ipe == ROOTPE) THEN
        ALLOCATE(fingerprint_i(npes))
        DO i=1,npes
          fingerprint_i(i)%f64(1)   = recvbuf(1,i)
          fingerprint_i(i)%f64(2)   = recvbuf(2,i)
          fingerprint_i(i)%f48(1)   = recvbuf(3,i)
          fingerprint_i(i)%f48(2)   = recvbuf(4,i)
          fingerprint_i(i)%tl64(1)  = recvbuf(5,i)
          fingerprint_i(i)%tl64(2)  = recvbuf(6,i)
          fingerprint_i(i)%tl48(1)  = recvbuf(7,i)
          fingerprint_i(i)%tl48(2)  = recvbuf(8,i)
          fingerprint_i(i)%max_zero = recvbuf(9,i)
        END DO
        DEALLOCATE(sendbuf, recvbuf)
        
        ! concatenate them:
        IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
          WRITE (0,*) 'concatenate fingerprints'
        END IF
        fingerprint_merge = C_LOC(fingerprint_i(1))
        DO i=2,npes
          new_fingerprint   = my_concat_fingerprints(fingerprint_merge, C_LOC(fingerprint_i(i)));
          IF (i>2)  CALL deallocate_c(fingerprint_merge)
          fingerprint_merge = new_fingerprint
        END DO
        
        CALL my_encode_uuid(fingerprint_merge, uuid)
        CALL deallocate_c(fingerprint_merge)
        DEALLOCATE(fingerprint_i)
      END IF

    END IF
    DEALLOCATE(val_sorted, glbidx_sorted, sdispls, sendcounts, &
      &        glbidx_local, val_local, recvcounts, rdispls, glbidx)
    IF ((ipe == 0) .AND. (dbg_level > 0)) THEN
      WRITE (0,*) 'mo_util_uuid done.'
    END IF
#else
    ! non-MPI mode: execute sequential version
    ALLOCATE(permutation(nval))
    permutation(in_glbidx(:)) = (/ (i, i=1,nval) /)
    CALL uuid_generate_sequential(RESHAPE(in_val(:,permutation), (/ SIZE(in_val) /)), uuid)
    DEALLOCATE(permutation)
#endif
  END SUBROUTINE uuid_generate_parallel


  ! Parallel computation of UUID.
  ! Wrapper for 1D array input.
  !
  SUBROUTINE uuid_generate_parallel1D(comm, in_val, in_glbidx, glb_nval, uuid)
    INTEGER,        INTENT(IN)  :: comm
    REAL(C_DOUBLE), INTENT(IN)  :: in_val(:)
    INTEGER,        INTENT(IN)  :: in_glbidx(:)
    INTEGER,        INTENT(IN)  :: glb_nval
    TYPE(t_uuid),   INTENT(OUT) :: uuid
    ! local variables
    REAL(C_DOUBLE), ALLOCATABLE :: tmp_val(:,:)
    ALLOCATE(tmp_val(1,SIZE(in_val)))
    tmp_val(1,:) = in_val(:)
    CALL uuid_generate_parallel(comm, tmp_val, in_glbidx, glb_nval, uuid)
    DEALLOCATE(tmp_val)
  END SUBROUTINE uuid_generate_parallel1D


  SUBROUTINE clear_uuid(uuid)
    TYPE(t_uuid), INTENT(inout) :: uuid
    uuid%data(:) = INT(0, C_SIGNED_CHAR)
  END SUBROUTINE clear_uuid

END MODULE mo_util_uuid
