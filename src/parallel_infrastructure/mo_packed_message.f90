! A class to make sending/receiving of data packets easier.
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

#include "handle_mpi_error.inc"

MODULE mo_packed_message
  USE ISO_C_BINDING, ONLY: c_ptr, C_LOC, C_F_POINTER, C_SIGNED_CHAR
#if defined(__PGI)
  USE ISO_C_BINDING, ONLY: C_SIZE_T, C_SIZEOF
#endif
  USE mo_exception, ONLY: finish, message_text
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_kind, ONLY: sp, dp, i8
#ifndef NOMPI
  HANDLE_MPI_ERROR_USE
  USE mo_mpi, ONLY: p_int, p_get_bcast_role, MPI_SUCCESS
  USE mpi, ONLY: MPI_CHAR, MPI_STATUS_SIZE
# ifndef NO_MPI_CHOICE_ARG
  USE mpi, ONLY: MPI_Send, MPI_Recv, MPI_Bcast
# endif
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_PackedMessage

  INTEGER, PARAMETER, PUBLIC :: kPackOp = 1
  INTEGER, PARAMETER, PUBLIC :: kUnpackOp = 2
#if defined(__PGI)
  INTEGER :: tbytes(7) = -1
#else
! * obtain the storge size of supported data types in units of C_SIGNED_CHAR
! (pgf90 needs some workaround since STORAGE_SIZE appears to be broken)
  INTEGER, PARAMETER :: tbytes(7) = ([ STORAGE_SIZE(1), STORAGE_SIZE(1_i8), &
     & STORAGE_SIZE(1_C_SIGNED_CHAR), STORAGE_SIZE(.false.), &
     & STORAGE_SIZE(1._dp), STORAGE_SIZE(1._sp), STORAGE_SIZE('a') ] + &
     & STORAGE_SIZE(1_C_SIGNED_CHAR) - 1) / STORAGE_SIZE(1_C_SIGNED_CHAR)
#endif
  ! A t_PackedMessage IS used to bundle a number of different values
  ! together into a single message, that can be communicated via a
  ! single CALL.  It IS possible to have ANY number of communication
  ! steps between the packing AND unpacking, including zero (a PE
  ! unpacks its own DATA), AND two (a PE recieves a packed message
  ! AND passes it on, possibly via a different communicator).
  !
  ! If NOMPI IS defined, the communication routines are simply
  ! noops, the packing/unpacking still works as expected.
  !
  ! As an added bonus, this provides packerXXX() routines IN
  ! addition to the packXXX() AND unpackXXX() routines, which allow
  ! folding the packing AND unpacking into the same code. Ie.,
  ! instead of writing a routine containing
  !
  !   message%pack(foo)
  !   message%pack(bar)
  !   message%pack(baz)
  !
  ! AND a second routine containing
  !
  !   message%unpack(foo)
  !   message%unpack(baz) !Error: messed up sequence!
  !   message%unpack(bar) !Error: messed up sequence!
  !
  ! you can WRITE a single routine containing
  !
  !   message%packer(operation, foo)
  !   message%packer(operation, baz)
  !   message%packer(operation, bar)
  !
  ! knowing that it will be simply impossible to mix up the sequence
  ! when setting operation to kUnpackOp to unpack the message.
  !
  TYPE :: t_PackedMessage
    PRIVATE
    INTEGER(C_SIGNED_CHAR), ALLOCATABLE :: messageBuffer(:)
    INTEGER :: messageSize = 0, readPosition = 0
  CONTAINS
    PROCEDURE :: reset => PackedMessage_reset   

    PROCEDURE, PRIVATE :: packBlock => PackedMessage_packBlock

    PROCEDURE, PRIVATE :: packCharacterScalar => PackedMessage_packCharacterScalar
    PROCEDURE, PRIVATE :: packIntScalar => PackedMessage_packIntScalar
    PROCEDURE, PRIVATE :: packLongScalar => PackedMessage_packLongScalar
    PROCEDURE, PRIVATE :: packSingleScalar => PackedMessage_packSingleScalar
    PROCEDURE, PRIVATE :: packDoubleScalar => PackedMessage_packDoubleScalar
    PROCEDURE, PRIVATE :: packLogicalScalar => PackedMessage_packLogicalScalar
    PROCEDURE, PRIVATE :: packIntCcharScalar => PackedMessage_packIntCcharScalar

    PROCEDURE, PRIVATE :: packIntArray => PackedMessage_packIntArray
    PROCEDURE, PRIVATE :: packLongArray => PackedMessage_packLongArray
    PROCEDURE, PRIVATE :: packSingleArray => PackedMessage_packSingleArray
    PROCEDURE, PRIVATE :: packDoubleArray => PackedMessage_packDoubleArray
    PROCEDURE, PRIVATE :: packLogicalArray => PackedMessage_packLogicalArray
    PROCEDURE, PRIVATE :: packIntCcharArray => PackedMessage_packIntCcharArray

    GENERIC :: pack => packCharacterScalar, &
         &             packIntScalar, packLongScalar, &
         &             packSingleScalar, packDoubleScalar, &
         &             packLogicalScalar, packIntCcharScalar, &
         &             packIntArray, packLongArray, &
         &             packSingleArray, packDoubleArray, &
         &             packLogicalArray, packIntCcharArray

    PROCEDURE, PRIVATE :: unpackBlock => PackedMessage_unpackBlock

    PROCEDURE, PRIVATE :: unpackCharacterScalar => PackedMessage_unpackCharacterScalar
    PROCEDURE, PRIVATE :: unpackIntScalar => PackedMessage_unpackIntScalar
    PROCEDURE, PRIVATE :: unpackLongScalar => PackedMessage_unpackLongScalar
    PROCEDURE, PRIVATE :: unpackSingleScalar => PackedMessage_unpackSingleScalar
    PROCEDURE, PRIVATE :: unpackDoubleScalar => PackedMessage_unpackDoubleScalar
    PROCEDURE, PRIVATE :: unpackLogicalScalar => PackedMessage_unpackLogicalScalar
    PROCEDURE, PRIVATE :: unpackIntCcharScalar => PackedMessage_unpackIntCcharScalar

    PROCEDURE, PRIVATE :: unpackIntArray => PackedMessage_unpackIntArray
    PROCEDURE, PRIVATE :: unpackLongArray => PackedMessage_unpackLongArray        
    PROCEDURE, PRIVATE :: unpackSingleArray => PackedMessage_unpackSingleArray
    PROCEDURE, PRIVATE :: unpackDoubleArray => PackedMessage_unpackDoubleArray
    PROCEDURE, PRIVATE :: unpackLogicalArray => PackedMessage_unpackLogicalArray
    PROCEDURE, PRIVATE :: unpackIntCcharArray => PackedMessage_unpackIntCcharArray

    GENERIC :: unpack => unpackCharacterScalar, &
         &               unpackIntScalar, unpackLongScalar, &
         &               unpackSingleScalar, unpackDoubleScalar, &
         &               unpackLogicalScalar, unpackIntCcharScalar, &
         &               unpackIntArray, unpackLongArray, &
         &               unpackSingleArray, unpackDoubleArray, &
         &               unpackLogicalArray, unpackIntCcharArray

    PROCEDURE, PRIVATE :: packerCharacterScalar => PackedMessage_packerCharacterScalar
    PROCEDURE, PRIVATE :: packerIntScalar => PackedMessage_packerIntScalar
    PROCEDURE, PRIVATE :: packerLongScalar => PackedMessage_packerLongScalar
    PROCEDURE, PRIVATE :: packerSingleScalar => PackedMessage_packerSingleScalar
    PROCEDURE, PRIVATE :: packerDoubleScalar => PackedMessage_packerDoubleScalar
    PROCEDURE, PRIVATE :: packerLogicalScalar => PackedMessage_packerLogicalScalar
    PROCEDURE, PRIVATE :: packerIntCcharScalar => PackedMessage_packerIntCcharScalar

    PROCEDURE, PRIVATE :: packerIntArray => PackedMessage_packerIntArray
    PROCEDURE, PRIVATE :: packerLongArray => PackedMessage_packerLongArray        
    PROCEDURE, PRIVATE :: packerSingleArray => PackedMessage_packerSingleArray
    PROCEDURE, PRIVATE :: packerDoubleArray => PackedMessage_packerDoubleArray
    PROCEDURE, PRIVATE :: packerLogicalArray => PackedMessage_packerLogicalArray
    PROCEDURE, PRIVATE :: packerIntCcharArray => PackedMessage_packerIntCcharArray

    GENERIC :: packer => packerCharacterScalar, &
         &               packerIntScalar, packerLongScalar, &
         &               packerSingleScalar, packerDoubleScalar, &
         &               packerLogicalScalar, packerIntCcharScalar, &
         &               packerIntArray, packerLongArray, &
         &               packerSingleArray, packerDoubleArray, &
         &               packerLogicalArray, packerIntCcharArray

    PROCEDURE :: send => PackedMessage_send
    PROCEDURE :: recv => PackedMessage_recv
    PROCEDURE :: bcast => PackedMessage_bcast

    PROCEDURE, PRIVATE :: ensureSpace => PackedMessage_ensureSpace   ! protected, use only in subclasses that implement their own packing

  END TYPE t_PackedMessage

  CHARACTER(len=*), PARAMETER :: modname = "mo_packed_message"

CONTAINS

#if defined(__PGI)
  SUBROUTINE init_tbytes()
    INTEGER, PARAMETER :: crap_i(8) = 1
    INTEGER(i8), PARAMETER :: crap_i8(8) = 1_i8
    INTEGER(C_SIGNED_CHAR), PARAMETER :: crap_Cchar(8) = 1_C_SIGNED_CHAR
    LOGICAL, PARAMETER :: crap_l(8) = .false.
    REAL(dp), PARAMETER :: crap_dp(8) = 1._sp
    REAL(sp), PARAMETER :: crap_sp(8) = 1._dp
    CHARACTER(LEN=1), PARAMETER :: crap_c(8) = 'a'
    INTEGER(C_SIZE_T) :: tmp, schar

    schar = C_SIZEOF(crap_Cchar)
    tmp = C_SIZEOF(crap_i) + schar - 1
    tbytes(1) = INT(tmp) / schar
    tmp = C_SIZEOF(crap_i8) + schar - 1
    tbytes(2) = INT(tmp) / schar
    tmp = C_SIZEOF(crap_Cchar) + schar - 1
    tbytes(3) = INT(tmp) / schar
    tmp = C_SIZEOF(crap_l) + schar - 1
    tbytes(4) = INT(tmp) / schar
    tmp = C_SIZEOF(crap_dp) + schar - 1
    tbytes(5) = INT(tmp) / schar
    tmp = C_SIZEOF(crap_sp) + schar - 1
    tbytes(6) = INT(tmp) / schar
    tmp = C_SIZEOF(crap_c) + schar - 1
    tbytes(7) = INT(tmp) / schar
  END SUBROUTINE init_tbytes
#endif

  SUBROUTINE PackedMessage_reset(me)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me

    me%messageSize = 0
    me%readPosition = 0
  END SUBROUTINE PackedMessage_reset

  SUBROUTINE PackedMessage_ensureSpace(me, requiredSpace)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: requiredSpace
    INTEGER :: ierr, oldSize, newSize
    INTEGER(C_SIGNED_CHAR), ALLOCATABLE :: newBuffer(:)
    CHARACTER(*), PARAMETER :: routine = modname//":ensureSpace"

    oldSize = 0
    IF (ALLOCATED(me%messageBuffer)) oldSize = SIZE(me%messageBuffer)
    IF(oldSize >= me%messageSize + requiredSpace) RETURN
    newSize = MAX(1024, 2*oldSize, me%messageSize + requiredSpace)
    ALLOCATE(newBuffer(newSize), STAT = ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation failed")
    IF (ALLOCATED(me%messageBuffer)) &
      & newBuffer(1:oldSize) = me%messageBuffer(1:oldSize)
    CALL MOVE_ALLOC(newBuffer, me%messageBuffer)
  END SUBROUTINE PackedMessage_ensureSpace

! ********************************************************************
! * new implementation replaces the old implementation via TRANSFER
! * basically this resembles somewhat of a "reinterpret cast" of a 
! * (void*) pointer "cptr" (obtained via C_LOC() in calling routine)
! * to a (char*) pointer "vptr" using C_F_POINTER()
! * then copy the full body of the payload (length=asize) at once 
! * exploiting fortan array syntax
! *
! * avdantages:
! *    - fewer (implicit) copies
! *    - memcpy en bloc instead of copying byte by byte
! *    - avoids fortran CHARACTER-string handling as far as possible
! *    - (gcc >= 7 had problems with the old TRANSFER impemetation)
! *
! * you should not care about the details.
! * do not touch!
! ********************************************************************
  SUBROUTINE PackedMessage_packBlock(me, cptr, asize, tid)
    CLASS(t_PackedMessage), INTENT(INOUT), TARGET :: me
    TYPE(c_ptr), INTENT(INOUT) :: cptr
    INTEGER, INTENT(IN) :: asize, tid
    INTEGER(C_SIGNED_CHAR), POINTER :: vptr(:)
    INTEGER :: bsize

#if defined(__PGI)
    IF (tbytes(1) .EQ. -1) CALL init_tbytes()
#endif
    bsize = asize * tbytes(tid)
    CALL me%ensureSpace(bsize)
    CALL C_F_POINTER(cptr, vptr, [bsize])
    me%messageBuffer(me%messageSize+1:me%messageSize+bsize) = vptr(1:bsize)
    me%messageSize = me%messageSize + bsize
  END SUBROUTINE PackedMessage_packBlock

#define PM_packScalar(__type, __tid) \
    CLASS(t_PackedMessage), INTENT(INOUT) :: me; \
    __type, INTENT(IN), TARGET :: scalar; \
    TYPE(c_ptr) :: cptr; \
    ; \
    cptr = C_LOC(scalar); \
    CALL me%packBlock(cptr, 1, __tid);

  SUBROUTINE PackedMessage_packIntScalar(me, scalar)
  PM_packScalar(INTEGER , 1)
  END SUBROUTINE PackedMessage_packIntScalar

  SUBROUTINE PackedMessage_packLongScalar(me, scalar)
  PM_packScalar(INTEGER(i8) , 2)
  END SUBROUTINE PackedMessage_packLongScalar

  SUBROUTINE PackedMessage_packDoubleScalar(me, scalar)
  PM_packScalar(REAL(dp) , 5)
  END SUBROUTINE PackedMessage_packDoubleScalar

  SUBROUTINE PackedMessage_packSingleScalar(me, scalar)
  PM_packScalar(REAL(sp) , 6)
  END SUBROUTINE PackedMessage_packSingleScalar

  SUBROUTINE PackedMessage_packLogicalScalar(me, scalar)
  PM_packScalar(LOGICAL , 4)
  END SUBROUTINE PackedMessage_packLogicalScalar

  SUBROUTINE PackedMessage_packIntCcharScalar(me, scalar)
  PM_packScalar(INTEGER(C_SIGNED_CHAR) , 3)
  END SUBROUTINE PackedMessage_packIntCcharScalar

  SUBROUTINE PackedMessage_packCharacterScalar(me, scalar)
   CLASS(t_PackedMessage), INTENT(INOUT) :: me
   CHARACTER(*), INTENT(IN), TARGET :: scalar
   TYPE(c_ptr) :: cptr

   cptr = C_LOC(scalar(1:1))
   CALL me%packBlock(cptr, LEN(scalar), 7)
  END SUBROUTINE PackedMessage_packCharacterScalar

#undef PM_packScalar

#define PM_packArray(__type, __tid) \
    CLASS(t_PackedMessage), INTENT(INOUT) :: me; \
    __type, ALLOCATABLE, INTENT(IN), TARGET :: array(:); \
    TYPE(c_ptr) :: cptr; \
    INTEGER, TARGET :: asize; \
    ; \
    asize = 0; \
    IF (ALLOCATED(array)) asize = SIZE(array); \
    cptr = C_LOC(asize); \
    CALL me%packBlock(cptr, 1, 1); \
    IF (asize .GT. 0)  THEN; \
      cptr = C_LOC(array(1)); \
      CALL me%packBlock(cptr, asize, __tid); \
    END IF;

  SUBROUTINE PackedMessage_packIntArray(me, array)
  PM_packArray(INTEGER, 1)
  END SUBROUTINE PackedMessage_packIntArray

  SUBROUTINE PackedMessage_packLongArray(me, array)
  PM_packArray(INTEGER(i8), 2)
  END SUBROUTINE PackedMessage_packLongArray

  SUBROUTINE PackedMessage_packDoubleArray(me, array)
  PM_packArray(REAL(dp), 5)
  END SUBROUTINE PackedMessage_packDoubleArray

  SUBROUTINE PackedMessage_packSingleArray(me, array)
  PM_packArray(REAL(sp), 6)
  END SUBROUTINE PackedMessage_packSingleArray

  SUBROUTINE PackedMessage_packLogicalArray(me, array)
  PM_packArray(LOGICAL, 4)
  END SUBROUTINE PackedMessage_packLogicalArray

  SUBROUTINE PackedMessage_packIntCcharArray(me, array)
  PM_packArray(INTEGER(C_SIGNED_CHAR), 3)
  END SUBROUTINE PackedMessage_packIntCcharArray

#undef PM_packArray

! * see comment above PackedMessage_packBlock
! * you should not care.
! * do not touch!
  SUBROUTINE PackedMessage_unpackBlock(me, cptr, asize, tid)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    TYPE(c_ptr), INTENT(INOUT), TARGET :: cptr
    INTEGER, INTENT(IN) :: asize, tid
    CHARACTER(*), PARAMETER :: routine = modname//"unpackBlock"
    INTEGER(C_SIGNED_CHAR), POINTER :: vptr(:)
    INTEGER :: bsize

#if defined(__PGI)
    IF (tbytes(1) .EQ. -1) CALL init_tbytes()
#endif
    bsize = asize * tbytes(tid)
    IF (bsize + me%readPosition .GT. SIZE(me%messageBuffer)) THEN
      WRITE(message_text, "(4(a,i8))") "out of bounds read at ", &
        & me%readPosition+1, ":", me%readPosition+bsize, &
        & "; bufSize = ", SIZE(me%messageBuffer), " read size = ", bsize
      CALL finish(routine, message_text)
    END IF
    CALL C_F_POINTER(cptr, vptr, [bsize])
    vptr(1:bsize) = me%messageBuffer(me%readPosition+1:me%readPosition+bsize)
    me%readPosition = me%readPosition + bsize
  END SUBROUTINE PackedMessage_unpackBlock

#define PM_unpackScalar(__type, __tid) \
    CLASS(t_PackedMessage), INTENT(INOUT) :: me; \
    __type , INTENT(OUT), TARGET :: scalar; \
    TYPE(c_ptr) :: cptr; \
    ; \
    cptr = C_LOC(scalar); \
    CALL me%unpackBlock(cptr, 1, __tid);

  SUBROUTINE PackedMessage_unpackIntScalar(me, scalar)
  PM_unpackScalar(INTEGER, 1)
  END SUBROUTINE PackedMessage_unpackIntScalar

  SUBROUTINE PackedMessage_unpackLongScalar(me, scalar)
  PM_unpackScalar(INTEGER(i8), 2)
  END SUBROUTINE PackedMessage_unpackLongScalar

  SUBROUTINE PackedMessage_unpackDoubleScalar(me, scalar)
  PM_unpackScalar(REAL(dp), 5)
  END SUBROUTINE PackedMessage_unpackDoubleScalar

  SUBROUTINE PackedMessage_unpackSingleScalar(me, scalar)
  PM_unpackScalar(REAL(sp), 6)
  END SUBROUTINE PackedMessage_unpackSingleScalar

  SUBROUTINE PackedMessage_unpackLogicalScalar(me, scalar)
  PM_unpackScalar(LOGICAL, 4)
  END SUBROUTINE PackedMessage_unpackLogicalScalar

  SUBROUTINE PackedMessage_unpackIntCcharScalar(me, scalar)
  PM_unpackScalar(INTEGER(C_SIGNED_CHAR), 3)
  END SUBROUTINE PackedMessage_unpackIntCcharScalar

  SUBROUTINE PackedMessage_unpackCharacterScalar(me, scalar)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    CHARACTER(*) , INTENT(OUT), TARGET :: scalar
    TYPE(c_ptr) :: cptr

    cptr = C_LOC(scalar(1:1))
    CALL me%unpackBlock(cptr, LEN(scalar), 7)
  END SUBROUTINE PackedMessage_unpackCharacterScalar

#undef PM_unpackScalar

#define PM_unpackArray(__type, __tid) \
    CLASS(t_PackedMessage), INTENT(INOUT) :: me; \
    __type , ALLOCATABLE, INTENT(INOUT) :: array(:); \
    __type , ALLOCATABLE, TARGET :: tmp(:); \
    TYPE(c_ptr) :: cptr; \
    INTEGER, TARGET:: asize; \
    ; \
    cptr = C_LOC(asize); \
    CALL me%unpackBlock(cptr, 1, 1); \
    IF (ALLOCATED(array)) DEALLOCATE(array); \
    IF (asize .GT. 0) THEN; \
      ALLOCATE(tmp(asize)); \
      cptr = C_LOC(tmp(1)); \
      CALL me%unpackBlock(cptr, asize, __tid); \
      CALL MOVE_ALLOC(tmp, array); \
    END IF;

  SUBROUTINE PackedMessage_unpackIntArray(me, array)
  PM_unpackArray(INTEGER, 1)
  END SUBROUTINE PackedMessage_unpackIntArray

  SUBROUTINE PackedMessage_unpackLongArray(me, array)
  PM_unpackArray(INTEGER(i8), 2)
  END SUBROUTINE PackedMessage_unpackLongArray

  SUBROUTINE PackedMessage_unpackDoubleArray(me, array)
  PM_unpackArray(REAL(dp), 5)
  END SUBROUTINE PackedMessage_unpackDoubleArray

  SUBROUTINE PackedMessage_unpackSingleArray(me, array)
  PM_unpackArray(REAL(sp), 6)
  END SUBROUTINE PackedMessage_unpackSingleArray

  SUBROUTINE PackedMessage_unpackLogicalArray(me, array)
  PM_unpackArray(LOGICAL, 4)
  END SUBROUTINE PackedMessage_unpackLogicalArray

  SUBROUTINE PackedMessage_unpackIntCcharArray(me, array)
  PM_unpackArray(INTEGER(C_SIGNED_CHAR), 3)
  END SUBROUTINE PackedMessage_unpackIntCcharArray

#undef PM_unpackArray

#define PM_packerScalar(__type) \
    CLASS(t_PackedMessage), INTENT(INOUT) :: me; \
    INTEGER, INTENT(IN) :: op; \
    __type, INTENT(INOUT) :: scalar; \
    ; \
    IF (op .NE. kUnpackOp) THEN; \
      CALL me%pack(scalar); \
    ELSE; \
      CALL me%unpack(scalar); \
    END IF;

  SUBROUTINE PackedMessage_packerIntScalar(me, op, scalar)
  PM_packerScalar(INTEGER)
  END SUBROUTINE PackedMessage_packerIntScalar

  SUBROUTINE PackedMessage_packerLongScalar(me, op, scalar)
  PM_packerScalar(INTEGER(i8))
  END SUBROUTINE PackedMessage_packerLongScalar

  SUBROUTINE PackedMessage_packerDoubleScalar(me, op, scalar)
  PM_packerScalar(REAL(dp))
  END SUBROUTINE PackedMessage_packerDoubleScalar

  SUBROUTINE PackedMessage_packerSingleScalar(me, op, scalar)
  PM_packerScalar(REAL(sp))
  END SUBROUTINE PackedMessage_packerSingleScalar

  SUBROUTINE PackedMessage_packerLogicalScalar(me, op, scalar)
  PM_packerScalar(LOGICAL)
  END SUBROUTINE PackedMessage_packerLogicalScalar

  SUBROUTINE PackedMessage_packerIntCcharScalar(me, op, scalar)
  PM_packerScalar(INTEGER(C_SIGNED_CHAR))
  END SUBROUTINE PackedMessage_packerIntCcharScalar

  SUBROUTINE PackedMessage_packerCharacterScalar(me, op, scalar)
  PM_packerScalar(CHARACTER(*))
  END SUBROUTINE PackedMessage_packerCharacterScalar

#undef PM_packerScalar

#define PM_packerArray(__type) \
    CLASS(t_PackedMessage), INTENT(INOUT) :: me; \
    INTEGER, INTENT(IN) :: op; \
    __type, ALLOCATABLE, INTENT(INOUT) :: array(:); \
    ; \
    IF (op .NE. kUnpackOp) THEN; \
      CALL me%pack(array); \
    ELSE; \
      CALL me%unpack(array); \
    END IF;

  SUBROUTINE PackedMessage_packerIntArray(me, op, array)
  PM_packerArray(INTEGER)
  END SUBROUTINE PackedMessage_packerIntArray

  SUBROUTINE PackedMessage_packerLongArray(me, op, array)
  PM_packerArray(INTEGER(i8))
  END SUBROUTINE PackedMessage_packerLongArray

  SUBROUTINE PackedMessage_packerDoubleArray(me, op, array)
  PM_packerArray(REAL(dp))
  END SUBROUTINE PackedMessage_packerDoubleArray

  SUBROUTINE PackedMessage_packerSingleArray(me, op, array)
  PM_packerArray(REAL(sp))
  END SUBROUTINE PackedMessage_packerSingleArray

  SUBROUTINE PackedMessage_packerLogicalArray(me, op, array)
  PM_packerArray(LOGICAL)
  END SUBROUTINE PackedMessage_packerLogicalArray

  SUBROUTINE PackedMessage_packerIntCcharArray(me, op, array)
  PM_packerArray(INTEGER(C_SIGNED_CHAR))
  END SUBROUTINE PackedMessage_packerIntCcharArray

#undef PM_packerArray

  ! communication routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! These always send/receive two MPI messages: one with the SIZE of the packed message, AND one with the actual message.
  ! This allows us to ensure that the receiver(s) always have enough memory to store the entire incoming message.

  SUBROUTINE PackedMessage_send(me, dest, tag, comm)
    CLASS(t_PackedMessage), INTENT(IN) :: me
    INTEGER, INTENT(IN) :: dest, tag, comm
#ifndef NOMPI
    INTEGER :: ierr
    CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_send"

    CALL MPI_Send(me%messageSize, 1, p_int, dest, tag, comm, ierr)
    HANDLE_MPI_ERROR(ierr, "MPI_Send")
    CALL MPI_Send(me%messageBuffer, INT(me%messageSize), MPI_CHAR, dest, tag, comm, ierr)
    HANDLE_MPI_ERROR(ierr, "MPI_Send")
#endif
  END SUBROUTINE PackedMessage_send

  SUBROUTINE PackedMessage_recv(me, source, tag, comm)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: source, tag, comm
#ifndef NOMPI
    INTEGER :: ierr, msize, mstat(MPI_STATUS_SIZE)
    CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_recv"

    CALL me%reset()
    CALL MPI_Recv(msize, 1, p_int, source, tag, comm, mstat, ierr)
    HANDLE_MPI_ERROR(ierr, "MPI_Recv")
    CALL me%ensureSpace(msize)
    CALL MPI_Recv(me%messageBuffer, msize, MPI_CHAR, source, tag, comm, mstat, ierr)
    HANDLE_MPI_ERROR(ierr, "MPI_Recv")
    me%messageSize = msize
#endif
  END SUBROUTINE PackedMessage_recv

  SUBROUTINE PackedMessage_bcast(me, root, comm)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: root, comm
#ifndef NOMPI
    INTEGER :: ierr, msize
    LOGICAL :: isSource, isDestination
    CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_bcast"
    INTEGER(C_SIGNED_CHAR) :: dummy_c(0)

    CALL p_get_bcast_role(root, comm, isSource, isDestination)
    msize = me%messageSize
    CALL MPI_Bcast(msize, 1, p_int, root, comm, ierr)
    HANDLE_MPI_ERROR(ierr, "MPI_Bcast")
    IF (isDestination) THEN
      CALL me%reset()
      CALL me%ensureSpace(msize)
    END IF
    IF (isDestination .OR. isSource) THEN
      CALL MPI_Bcast(me%messageBuffer, msize, MPI_CHAR, root, comm, ierr)
    ELSE
      CALL MPI_Bcast(dummy_c, 0, MPI_CHAR, root, comm, ierr)
    END IF
    HANDLE_MPI_ERROR(ierr, "MPI_Bcast")
    IF(isDestination) me%messageSize = msize
#endif
  END SUBROUTINE PackedMessage_bcast

END MODULE mo_packed_message
