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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A generic hash table
! A container for any kind of objects that are stored under any kind
! of keys, with O(1) complexity for insertion, lookup, and removal.
! Implementation note:
! The loops look more complex than they are:
! Their bodies are executed only once in most cases, sometimes zero times, and more than once in even less cases.
! This is due to the fact that we grow the hash table to be at least as large as we have entries.
! As such, the lists have at most one entry on average.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_hash_table
  USE ISO_C_BINDING, ONLY: c_int32_t
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  PUBLIC :: t_HashTable, hashTable_make, t_HashIterator
  PRIVATE

  ABSTRACT INTERFACE
    INTEGER(c_int32_t) FUNCTION f_hashFunction(key)
      IMPORT c_int32_t
      CLASS(*), POINTER, INTENT(IN) :: key
    END FUNCTION f_hashFunction

    LOGICAL FUNCTION f_equalKeysFunction(keyA, keyB)
      CLASS(*), POINTER, INTENT(IN) :: keyA, keyB
    END FUNCTION f_equalKeysFunction
  END INTERFACE

  TYPE :: t_HashTable
    PRIVATE
    PROCEDURE(f_hashFunction), NOPASS, POINTER :: getHash
    PROCEDURE(f_equalKeysFunction), NOPASS, POINTER :: equalKeys
    TYPE(t_HashEntryPtr), POINTER :: table(:) ! grows in powers of two
    INTEGER(c_int32_t) :: entryCount = 0 ! the actual entry count
    INTEGER :: hashBits = 5 ! the current count of bits that are used to index the hash table
  CONTAINS
    PROCEDURE :: setEntry => hashTable_setEntry
    PROCEDURE :: removeEntry => hashTable_removeEntry
    PROCEDURE :: getEntry => hashTable_getEntry
    PROCEDURE :: destruct => hashTable_destruct
    PROCEDURE :: getEntryCount => hashTable_getEntryCount

    PROCEDURE, PRIVATE :: findBin => hashTable_findBin
    PROCEDURE, PRIVATE :: growTable => hashTable_growTable
    PROCEDURE, PRIVATE :: removeFromList => hashTable_removeFromList
  END TYPE

  ! provides sequential access to all entries of a hash table
  ! <<< WARNING >>>: iterators are invalidated by setEntry(),
  ! removeEntry(), AND destruct() calls on the corresponding hash
  ! table
  TYPE :: t_HashIterator
    PRIVATE
    TYPE(t_HashTable), POINTER :: table
    INTEGER :: curBin = 0
    TYPE(t_HashEntry), POINTER :: curEntry => NULL()
  CONTAINS
    PROCEDURE :: init => hashIterator_init
    PROCEDURE :: nextEntry => hashIterator_nextEntry ! returns .TRUE. IF the operation was successful
    PROCEDURE :: reset => hashIterator_reset
  END TYPE

  TYPE :: t_HashEntryPtr
    TYPE(t_HashEntry), POINTER :: ptr => NULL()
  END TYPE

  TYPE :: t_HashEntry
    CLASS(*), POINTER :: key, val
    TYPE(t_HashEntryPtr) :: next
    INTEGER(c_int32_t) :: hash
  END TYPE

  CHARACTER(LEN=*), PARAMETER :: modname = "mo_hash_table"

  INTEGER, PARAMETER :: SUCCESS = 0

CONTAINS

  INTEGER FUNCTION hashTable_getEntryCount(me) RESULT(entryCount)
    CLASS(t_HashTable), INTENT(IN) :: me

    entryCount = me%entryCount
  END FUNCTION hashTable_getEntryCount

  FUNCTION hashTable_make(hashFunction, compareFunction) RESULT(resultVar)
    PROCEDURE(f_hashFunction) :: hashFunction
    PROCEDURE(f_equalKeysFunction) :: compareFunction
    TYPE(t_HashTable), POINTER :: resultVar
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":hashTable_make"
    INTEGER :: error

    ALLOCATE (resultVar, STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    resultVar%getHash => hashFunction
    resultVar%equalKeys => compareFunction
    ALLOCATE (resultVar%table(2**resultVar%hashBits), STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END FUNCTION hashTable_make

  FUNCTION hashTable_findBin(me, hash_in) RESULT(res)
    CLASS(t_HashTable), INTENT(IN) :: me
    INTEGER(c_int32_t), INTENT(IN) :: hash_in
    TYPE(t_HashEntryPtr), POINTER :: res
    INTEGER(c_int32_t) :: redHash, i, hash

    hash = MERGE(hash_in, NOT(hash_in), hash_in .LT. 0)
    redHash = 0
    DO i = 1, (31 + me%hashBits - 1)/me%hashBits
      redHash = IEOR(redHash, hash)
      hash = ISHFT(hash, -me%hashBits)
    END DO
    redHash = IAND(redHash, 2**me%hashBits - 1) + 1
    res => me%table(redHash)
  END FUNCTION hashTable_findBin

  SUBROUTINE hashTable_growTable(me)
    CLASS(t_HashTable), INTENT(INOUT) :: me

    CHARACTER(LEN=*), PARAMETER :: routine = modname//":hashTable_growTable"
    TYPE(t_HashEntry), POINTER :: curEntry, nextEntry
    TYPE(t_HashEntryPtr), POINTER :: oldTable(:), bin
    INTEGER :: i, error

    oldTable => me%table

    ! Create a new empty table.
    me%hashBits = me%hashBits + 1
    ALLOCATE (me%table(2**me%hashBits), STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    ! Move over the contents of the old table.
    DO i = 1, SIZE(oldTable, 1)
      curEntry => oldTable(i)%ptr
      DO WHILE (ASSOCIATED(curEntry))
        nextEntry => curEntry%next%ptr

        ! Insert curEntry into the new table.
        bin => me%findBin(curEntry%hash)
        curEntry%next%ptr => bin%ptr
        bin%ptr => curEntry

        curEntry => nextEntry
      END DO
    END DO

    ! Cleanup the old table.
    DEALLOCATE (oldTable)
  END SUBROUTINE hashTable_growTable

  SUBROUTINE hashTable_removeFromList(me, list, key, hash)
    CLASS(t_HashTable), INTENT(INOUT) :: me
    TYPE(t_HashEntryPtr), POINTER, INTENT(INOUT) :: list
    CLASS(*), POINTER, INTENT(IN) :: key
    INTEGER(c_int32_t), INTENT(IN) :: hash
    TYPE(t_HashEntry), POINTER :: curEntry
    TYPE(t_HashEntryPtr), POINTER :: iterator

    iterator => list
    DO WHILE (ASSOCIATED(iterator%ptr))
      curEntry => iterator%ptr
      IF (curEntry%hash == hash) THEN
        IF (me%equalKeys(curEntry%key, key)) THEN
          ! remove from list
          iterator%ptr => curEntry%next%ptr
          me%entryCount = me%entryCount - 1

          ! destroy the entry
          DEALLOCATE (curEntry%key, curEntry%val)
          DEALLOCATE (curEntry)
          CYCLE
        END IF
      END IF
      iterator => curEntry%next ! point to next entry
    END DO
  END SUBROUTINE hashTable_removeFromList

  ! The hash table takes possession of both the key and the val and will DEALLOCATE() them eventually.
  SUBROUTINE hashTable_setEntry(me, key, val)
    CLASS(t_HashTable), INTENT(INOUT) :: me
    CLASS(*), POINTER, INTENT(IN) :: key, val

    CHARACTER(LEN=*), PARAMETER :: routine = modname//":hashTable_setEntry"
    TYPE(t_HashEntry), POINTER :: newEntry
    TYPE(t_HashEntryPtr), POINTER :: bin
    INTEGER :: error

    ! Prepare the new entry.
    ALLOCATE (newEntry, STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    newEntry%key => key
    newEntry%val => val
    newEntry%hash = me%getHash(key)

    ! If there is a prexeisting entry for this key, remove it.
    bin => me%findBin(newEntry%hash)
    CALL me%removeFromList(bin, key, newEntry%hash)

    ! Add the new entry to the list.
    newEntry%next%ptr => bin%ptr
    bin%ptr => newEntry
    me%entryCount = me%entryCount + 1

    ! Check whether we need to grow the table.
    IF (me%entryCount == SIZE(me%table, 1)) CALL me%growTable()
  END SUBROUTINE hashTable_setEntry

  SUBROUTINE hashTable_removeEntry(me, key)
    CLASS(t_HashTable), INTENT(INOUT) :: me
    CLASS(*), POINTER, INTENT(IN) :: key

    INTEGER(c_int32_t) :: hash
    TYPE(t_HashEntryPtr), POINTER :: bin

    hash = me%getHash(key)
    bin => me%findBin(hash)
    CALL me%removeFromList(bin, key, hash)
  END SUBROUTINE hashTable_removeEntry

  FUNCTION hashTable_getEntry(me, key) RESULT(resultVar)
    CLASS(t_HashTable), INTENT(IN) :: me
    CLASS(*), POINTER, INTENT(IN) :: key
    CLASS(*), POINTER :: resultVar

    INTEGER(c_int32_t) :: hash
    TYPE(t_HashEntryPtr), POINTER :: bin
    TYPE(t_HashEntry), POINTER :: curEntry

    resultVar => NULL()
    hash = me%getHash(key)
    bin => me%findBin(hash)
    curEntry => bin%ptr
    DO WHILE (ASSOCIATED(curEntry))
      IF (curEntry%hash == hash) THEN
        IF (me%equalKeys(curEntry%key, key)) THEN
          resultVar => curEntry%val
          RETURN
        END IF
      END IF
      curEntry => curEntry%next%ptr
    END DO
  END FUNCTION hashTable_getEntry

  SUBROUTINE hashTable_destruct(me)
    CLASS(t_HashTable), INTENT(INOUT) :: me

    INTEGER :: i
    TYPE(t_HashEntry), POINTER :: curEntry, nextEntry

    DO i = 1, SIZE(me%table)
      curEntry => me%table(i)%ptr
      DO WHILE (ASSOCIATED(curEntry))
        nextEntry => curEntry%next%ptr
        DEALLOCATE (curEntry%key, curEntry%val)
        DEALLOCATE (curEntry)
        curEntry => nextEntry
      END DO
    END DO
    DEALLOCATE (me%table)
  END SUBROUTINE hashTable_destruct

  SUBROUTINE hashIterator_init(me, table)
    CLASS(t_HashIterator), INTENT(INOUT) :: me
    TYPE(t_HashTable), POINTER, INTENT(IN) :: table

    me%table => table
  END SUBROUTINE hashIterator_init

  SUBROUTINE hashIterator_reset(me)
    CLASS(t_HashIterator), INTENT(INOUT) :: me

    me%curBin = 0
    NULLIFY (me%curEntry)
  END SUBROUTINE hashIterator_reset

  LOGICAL FUNCTION hashIterator_nextEntry(me, key, val) RESULT(resultVar)
    CLASS(t_HashIterator), INTENT(INOUT) :: me
    CLASS(*), POINTER, INTENT(OUT) :: key, val
    LOGICAL :: cont

    NULLIFY (key, val)
    cont = .TRUE.
    ! try to ADVANCE within the current bin
    IF (ASSOCIATED(me%curEntry)) me%curEntry => me%curEntry%next%ptr
    ! search for the next entry
    DO WHILE (cont)
      ! check whether we have found the next entry
      IF (ASSOCIATED(me%curEntry)) THEN
        key => me%curEntry%key
        val => me%curEntry%val
        cont = .FALSE.
      ELSE IF (me%curBin .LT. SIZE(me%table%table, 1)) THEN
        ! check whether we can ADVANCE ANY further ADVANCE to the next bin
        me%curBin = me%curBin + 1
        me%curEntry => me%table%table(me%curBin)%ptr
      ELSE
        cont = .FALSE.
      END IF
    END DO
    resultVar = ASSOCIATED(key)
  END FUNCTION hashIterator_nextEntry

END MODULE mo_hash_table
