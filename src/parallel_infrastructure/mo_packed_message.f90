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


MODULE mo_packed_message
  USE ISO_C_BINDING, ONLY: c_ptr, C_LOC, C_F_POINTER, C_SIGNED_CHAR
  USE mo_exception, ONLY: finish, message_text
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_kind, ONLY: sp, dp, i8


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_PackedMessage

  INTEGER, PARAMETER, PUBLIC :: kPackOp = 1
  INTEGER, PARAMETER, PUBLIC :: kUnpackOp = 2



! * obtain the storge size of supported data types in units of C_SIGNED_CHAR
! (pgf90 needs some workaround since STORAGE_SIZE appears to be broken)
  INTEGER, PARAMETER :: tbytes(7) = ([ STORAGE_SIZE(1), STORAGE_SIZE(1_i8), &
     & STORAGE_SIZE(1_C_SIGNED_CHAR), STORAGE_SIZE(.false.), &
     & STORAGE_SIZE(1._dp), STORAGE_SIZE(1._sp), STORAGE_SIZE('a') ] + &
     & STORAGE_SIZE(1_C_SIGNED_CHAR) - 1) / STORAGE_SIZE(1_C_SIGNED_CHAR)

  ! A t_PackedMessage IS used to bundle a number of different values
  ! together into a single message, that can be communicated via a
  ! single CALL.  It IS possible to have ANY number of communication
  ! steps between the packing AND unpacking, including zero (a PE
  ! unpacks its own DATA), AND two (a PE recieves a packed message
  ! AND passes it on, possibly via a different communicator).
  !
  ! If 1 IS defined, the communication routines are simply
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

    bsize = asize * tbytes(tid)
    CALL me%ensureSpace(bsize)
    CALL C_F_POINTER(cptr, vptr, [bsize])
    me%messageBuffer(me%messageSize+1:me%messageSize+bsize) = vptr(1:bsize)
    me%messageSize = me%messageSize + bsize
  END SUBROUTINE PackedMessage_packBlock


  SUBROUTINE PackedMessage_packIntScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER , INTENT(IN), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%packBlock(cptr, 1,  1);
  END SUBROUTINE PackedMessage_packIntScalar

  SUBROUTINE PackedMessage_packLongScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(i8) , INTENT(IN), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%packBlock(cptr, 1,  2);
  END SUBROUTINE PackedMessage_packLongScalar

  SUBROUTINE PackedMessage_packDoubleScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(dp) , INTENT(IN), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%packBlock(cptr, 1,  5);
  END SUBROUTINE PackedMessage_packDoubleScalar

  SUBROUTINE PackedMessage_packSingleScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(sp) , INTENT(IN), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%packBlock(cptr, 1,  6);
  END SUBROUTINE PackedMessage_packSingleScalar

  SUBROUTINE PackedMessage_packLogicalScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     LOGICAL , INTENT(IN), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%packBlock(cptr, 1,  4);
  END SUBROUTINE PackedMessage_packLogicalScalar

  SUBROUTINE PackedMessage_packIntCcharScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(C_SIGNED_CHAR) , INTENT(IN), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%packBlock(cptr, 1,  3);
  END SUBROUTINE PackedMessage_packIntCcharScalar

  SUBROUTINE PackedMessage_packCharacterScalar(me, scalar)
   CLASS(t_PackedMessage), INTENT(INOUT) :: me
   CHARACTER(*), INTENT(IN), TARGET :: scalar
   TYPE(c_ptr) :: cptr

   cptr = C_LOC(scalar(1:1))
   CALL me%packBlock(cptr, LEN(scalar), 7)
  END SUBROUTINE PackedMessage_packCharacterScalar



  SUBROUTINE PackedMessage_packIntArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, ALLOCATABLE, INTENT(IN), TARGET :: array(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET :: asize;     ;     asize = 0;     IF (ALLOCATED(array)) asize = SIZE(array);     cptr = C_LOC(asize);     CALL me%packBlock(cptr, 1, 1);     IF (asize .GT. 0)  THEN;       cptr = C_LOC(array(1));       CALL me%packBlock(cptr, asize,  1);     END IF;
  END SUBROUTINE PackedMessage_packIntArray

  SUBROUTINE PackedMessage_packLongArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(i8), ALLOCATABLE, INTENT(IN), TARGET :: array(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET :: asize;     ;     asize = 0;     IF (ALLOCATED(array)) asize = SIZE(array);     cptr = C_LOC(asize);     CALL me%packBlock(cptr, 1, 1);     IF (asize .GT. 0)  THEN;       cptr = C_LOC(array(1));       CALL me%packBlock(cptr, asize,  2);     END IF;
  END SUBROUTINE PackedMessage_packLongArray

  SUBROUTINE PackedMessage_packDoubleArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(dp), ALLOCATABLE, INTENT(IN), TARGET :: array(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET :: asize;     ;     asize = 0;     IF (ALLOCATED(array)) asize = SIZE(array);     cptr = C_LOC(asize);     CALL me%packBlock(cptr, 1, 1);     IF (asize .GT. 0)  THEN;       cptr = C_LOC(array(1));       CALL me%packBlock(cptr, asize,  5);     END IF;
  END SUBROUTINE PackedMessage_packDoubleArray

  SUBROUTINE PackedMessage_packSingleArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(sp), ALLOCATABLE, INTENT(IN), TARGET :: array(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET :: asize;     ;     asize = 0;     IF (ALLOCATED(array)) asize = SIZE(array);     cptr = C_LOC(asize);     CALL me%packBlock(cptr, 1, 1);     IF (asize .GT. 0)  THEN;       cptr = C_LOC(array(1));       CALL me%packBlock(cptr, asize,  6);     END IF;
  END SUBROUTINE PackedMessage_packSingleArray

  SUBROUTINE PackedMessage_packLogicalArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     LOGICAL, ALLOCATABLE, INTENT(IN), TARGET :: array(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET :: asize;     ;     asize = 0;     IF (ALLOCATED(array)) asize = SIZE(array);     cptr = C_LOC(asize);     CALL me%packBlock(cptr, 1, 1);     IF (asize .GT. 0)  THEN;       cptr = C_LOC(array(1));       CALL me%packBlock(cptr, asize,  4);     END IF;
  END SUBROUTINE PackedMessage_packLogicalArray

  SUBROUTINE PackedMessage_packIntCcharArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(C_SIGNED_CHAR), ALLOCATABLE, INTENT(IN), TARGET :: array(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET :: asize;     ;     asize = 0;     IF (ALLOCATED(array)) asize = SIZE(array);     cptr = C_LOC(asize);     CALL me%packBlock(cptr, 1, 1);     IF (asize .GT. 0)  THEN;       cptr = C_LOC(array(1));       CALL me%packBlock(cptr, asize,  3);     END IF;
  END SUBROUTINE PackedMessage_packIntCcharArray


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


  SUBROUTINE PackedMessage_unpackIntScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER , INTENT(OUT), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%unpackBlock(cptr, 1,  1);
  END SUBROUTINE PackedMessage_unpackIntScalar

  SUBROUTINE PackedMessage_unpackLongScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(i8) , INTENT(OUT), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%unpackBlock(cptr, 1,  2);
  END SUBROUTINE PackedMessage_unpackLongScalar

  SUBROUTINE PackedMessage_unpackDoubleScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(dp) , INTENT(OUT), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%unpackBlock(cptr, 1,  5);
  END SUBROUTINE PackedMessage_unpackDoubleScalar

  SUBROUTINE PackedMessage_unpackSingleScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(sp) , INTENT(OUT), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%unpackBlock(cptr, 1,  6);
  END SUBROUTINE PackedMessage_unpackSingleScalar

  SUBROUTINE PackedMessage_unpackLogicalScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     LOGICAL , INTENT(OUT), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%unpackBlock(cptr, 1,  4);
  END SUBROUTINE PackedMessage_unpackLogicalScalar

  SUBROUTINE PackedMessage_unpackIntCcharScalar(me, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(C_SIGNED_CHAR) , INTENT(OUT), TARGET :: scalar;     TYPE(c_ptr) :: cptr;     ;     cptr = C_LOC(scalar);     CALL me%unpackBlock(cptr, 1,  3);
  END SUBROUTINE PackedMessage_unpackIntCcharScalar

  SUBROUTINE PackedMessage_unpackCharacterScalar(me, scalar)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    CHARACTER(*) , INTENT(OUT), TARGET :: scalar
    TYPE(c_ptr) :: cptr

    cptr = C_LOC(scalar(1:1))
    CALL me%unpackBlock(cptr, LEN(scalar), 7)
  END SUBROUTINE PackedMessage_unpackCharacterScalar



  SUBROUTINE PackedMessage_unpackIntArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER , ALLOCATABLE, INTENT(INOUT) :: array(:);     INTEGER , ALLOCATABLE, TARGET :: tmp(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET:: asize;     ;     cptr = C_LOC(asize);     CALL me%unpackBlock(cptr, 1, 1);     IF (ALLOCATED(array)) DEALLOCATE(array);     IF (asize .GT. 0) THEN;       ALLOCATE(tmp(asize));       cptr = C_LOC(tmp(1));       CALL me%unpackBlock(cptr, asize,  1);       CALL MOVE_ALLOC(tmp, array);     END IF;
  END SUBROUTINE PackedMessage_unpackIntArray

  SUBROUTINE PackedMessage_unpackLongArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(i8) , ALLOCATABLE, INTENT(INOUT) :: array(:);     INTEGER(i8) , ALLOCATABLE, TARGET :: tmp(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET:: asize;     ;     cptr = C_LOC(asize);     CALL me%unpackBlock(cptr, 1, 1);     IF (ALLOCATED(array)) DEALLOCATE(array);     IF (asize .GT. 0) THEN;       ALLOCATE(tmp(asize));       cptr = C_LOC(tmp(1));       CALL me%unpackBlock(cptr, asize,  2);       CALL MOVE_ALLOC(tmp, array);     END IF;
  END SUBROUTINE PackedMessage_unpackLongArray

  SUBROUTINE PackedMessage_unpackDoubleArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(dp) , ALLOCATABLE, INTENT(INOUT) :: array(:);     REAL(dp) , ALLOCATABLE, TARGET :: tmp(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET:: asize;     ;     cptr = C_LOC(asize);     CALL me%unpackBlock(cptr, 1, 1);     IF (ALLOCATED(array)) DEALLOCATE(array);     IF (asize .GT. 0) THEN;       ALLOCATE(tmp(asize));       cptr = C_LOC(tmp(1));       CALL me%unpackBlock(cptr, asize,  5);       CALL MOVE_ALLOC(tmp, array);     END IF;
  END SUBROUTINE PackedMessage_unpackDoubleArray

  SUBROUTINE PackedMessage_unpackSingleArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     REAL(sp) , ALLOCATABLE, INTENT(INOUT) :: array(:);     REAL(sp) , ALLOCATABLE, TARGET :: tmp(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET:: asize;     ;     cptr = C_LOC(asize);     CALL me%unpackBlock(cptr, 1, 1);     IF (ALLOCATED(array)) DEALLOCATE(array);     IF (asize .GT. 0) THEN;       ALLOCATE(tmp(asize));       cptr = C_LOC(tmp(1));       CALL me%unpackBlock(cptr, asize,  6);       CALL MOVE_ALLOC(tmp, array);     END IF;
  END SUBROUTINE PackedMessage_unpackSingleArray

  SUBROUTINE PackedMessage_unpackLogicalArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     LOGICAL , ALLOCATABLE, INTENT(INOUT) :: array(:);     LOGICAL , ALLOCATABLE, TARGET :: tmp(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET:: asize;     ;     cptr = C_LOC(asize);     CALL me%unpackBlock(cptr, 1, 1);     IF (ALLOCATED(array)) DEALLOCATE(array);     IF (asize .GT. 0) THEN;       ALLOCATE(tmp(asize));       cptr = C_LOC(tmp(1));       CALL me%unpackBlock(cptr, asize,  4);       CALL MOVE_ALLOC(tmp, array);     END IF;
  END SUBROUTINE PackedMessage_unpackLogicalArray

  SUBROUTINE PackedMessage_unpackIntCcharArray(me, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER(C_SIGNED_CHAR) , ALLOCATABLE, INTENT(INOUT) :: array(:);     INTEGER(C_SIGNED_CHAR) , ALLOCATABLE, TARGET :: tmp(:);     TYPE(c_ptr) :: cptr;     INTEGER, TARGET:: asize;     ;     cptr = C_LOC(asize);     CALL me%unpackBlock(cptr, 1, 1);     IF (ALLOCATED(array)) DEALLOCATE(array);     IF (asize .GT. 0) THEN;       ALLOCATE(tmp(asize));       cptr = C_LOC(tmp(1));       CALL me%unpackBlock(cptr, asize,  3);       CALL MOVE_ALLOC(tmp, array);     END IF;
  END SUBROUTINE PackedMessage_unpackIntCcharArray



  SUBROUTINE PackedMessage_packerIntScalar(me, op, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     INTEGER, INTENT(INOUT) :: scalar;     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(scalar);     ELSE;       CALL me%unpack(scalar);     END IF;
  END SUBROUTINE PackedMessage_packerIntScalar

  SUBROUTINE PackedMessage_packerLongScalar(me, op, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     INTEGER(i8), INTENT(INOUT) :: scalar;     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(scalar);     ELSE;       CALL me%unpack(scalar);     END IF;
  END SUBROUTINE PackedMessage_packerLongScalar

  SUBROUTINE PackedMessage_packerDoubleScalar(me, op, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     REAL(dp), INTENT(INOUT) :: scalar;     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(scalar);     ELSE;       CALL me%unpack(scalar);     END IF;
  END SUBROUTINE PackedMessage_packerDoubleScalar

  SUBROUTINE PackedMessage_packerSingleScalar(me, op, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     REAL(sp), INTENT(INOUT) :: scalar;     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(scalar);     ELSE;       CALL me%unpack(scalar);     END IF;
  END SUBROUTINE PackedMessage_packerSingleScalar

  SUBROUTINE PackedMessage_packerLogicalScalar(me, op, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     LOGICAL, INTENT(INOUT) :: scalar;     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(scalar);     ELSE;       CALL me%unpack(scalar);     END IF;
  END SUBROUTINE PackedMessage_packerLogicalScalar

  SUBROUTINE PackedMessage_packerIntCcharScalar(me, op, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     INTEGER(C_SIGNED_CHAR), INTENT(INOUT) :: scalar;     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(scalar);     ELSE;       CALL me%unpack(scalar);     END IF;
  END SUBROUTINE PackedMessage_packerIntCcharScalar

  SUBROUTINE PackedMessage_packerCharacterScalar(me, op, scalar)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     CHARACTER(*), INTENT(INOUT) :: scalar;     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(scalar);     ELSE;       CALL me%unpack(scalar);     END IF;
  END SUBROUTINE PackedMessage_packerCharacterScalar



  SUBROUTINE PackedMessage_packerIntArray(me, op, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     INTEGER, ALLOCATABLE, INTENT(INOUT) :: array(:);     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(array);     ELSE;       CALL me%unpack(array);     END IF;
  END SUBROUTINE PackedMessage_packerIntArray

  SUBROUTINE PackedMessage_packerLongArray(me, op, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     INTEGER(i8), ALLOCATABLE, INTENT(INOUT) :: array(:);     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(array);     ELSE;       CALL me%unpack(array);     END IF;
  END SUBROUTINE PackedMessage_packerLongArray

  SUBROUTINE PackedMessage_packerDoubleArray(me, op, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     REAL(dp), ALLOCATABLE, INTENT(INOUT) :: array(:);     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(array);     ELSE;       CALL me%unpack(array);     END IF;
  END SUBROUTINE PackedMessage_packerDoubleArray

  SUBROUTINE PackedMessage_packerSingleArray(me, op, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     REAL(sp), ALLOCATABLE, INTENT(INOUT) :: array(:);     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(array);     ELSE;       CALL me%unpack(array);     END IF;
  END SUBROUTINE PackedMessage_packerSingleArray

  SUBROUTINE PackedMessage_packerLogicalArray(me, op, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     LOGICAL, ALLOCATABLE, INTENT(INOUT) :: array(:);     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(array);     ELSE;       CALL me%unpack(array);     END IF;
  END SUBROUTINE PackedMessage_packerLogicalArray

  SUBROUTINE PackedMessage_packerIntCcharArray(me, op, array)
  CLASS(t_PackedMessage), INTENT(INOUT) :: me;     INTEGER, INTENT(IN) :: op;     INTEGER(C_SIGNED_CHAR), ALLOCATABLE, INTENT(INOUT) :: array(:);     ;     IF (op .NE. kUnpackOp) THEN;       CALL me%pack(array);     ELSE;       CALL me%unpack(array);     END IF;
  END SUBROUTINE PackedMessage_packerIntCcharArray


  ! communication routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! These always send/receive two MPI messages: one with the SIZE of the packed message, AND one with the actual message.
  ! This allows us to ensure that the receiver(s) always have enough memory to store the entire incoming message.

  SUBROUTINE PackedMessage_send(me, dest, tag, comm)
    CLASS(t_PackedMessage), INTENT(IN) :: me
    INTEGER, INTENT(IN) :: dest, tag, comm
  END SUBROUTINE PackedMessage_send

  SUBROUTINE PackedMessage_recv(me, source, tag, comm)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: source, tag, comm
  END SUBROUTINE PackedMessage_recv

  SUBROUTINE PackedMessage_bcast(me, root, comm)
    CLASS(t_PackedMessage), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: root, comm
  END SUBROUTINE PackedMessage_bcast

END MODULE mo_packed_message
