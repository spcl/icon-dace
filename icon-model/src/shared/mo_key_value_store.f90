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

! This module provides a simple key value storage based on string comparisons
! using the generic hash tables provided by mo_hash_table.
! As method to calculate hash keys, the DJB algorithm is used.

MODULE mo_key_value_store
  USE mo_exception,             ONLY: finish, message
  USE mo_hash_table,            ONLY: t_HashTable, hashTable_make, t_HashIterator
  USE mo_kind,                  ONLY: wp
  USE mo_mpi,                   ONLY: p_get_bcast_role
  USE mo_packed_message,        ONLY: t_PackedMessage
  USE mo_util_string,           ONLY: tolower
  USE mo_util_texthash,         ONLY: text_hash, text_isEqual, sel_char
#ifdef __PGI
  USE mo_util_texthash,         ONLY: t_char_workaround
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_key_value_store

  TYPE :: t_key_value_store
    PRIVATE
    LOGICAL, PUBLIC :: is_init = .false.
    LOGICAL :: lcase_sensitive = .true.
    TYPE(t_HashTable), POINTER :: table => NULL()
  CONTAINS
    PROCEDURE, PRIVATE :: put_c => key_value_store_put_c
    PROCEDURE, PRIVATE :: put_r => key_value_store_put_r
    PROCEDURE, PRIVATE :: put_i => key_value_store_put_i
    PROCEDURE, PRIVATE :: put_l => key_value_store_put_l
    GENERIC, PUBLIC :: put => put_c, put_r, put_i, put_l
    PROCEDURE, PRIVATE :: put_internal => key_value_store_put_internal
    PROCEDURE, PRIVATE :: get_c => key_value_store_get_c
    PROCEDURE, PRIVATE :: get_r => key_value_store_get_r
    PROCEDURE, PRIVATE :: get_i => key_value_store_get_i
    PROCEDURE, PRIVATE :: get_l => key_value_store_get_l
    GENERIC, PUBLIC :: get => get_c, get_r, get_i, get_l
    PROCEDURE, PRIVATE :: get_internal => key_value_store_get_internal
    PROCEDURE, PRIVATE :: input_pmsg => key_value_store_input_pmsg
    GENERIC, PUBLIC :: input => input_pmsg
    PROCEDURE, PUBLIC :: output => key_value_store_output
    PROCEDURE, PUBLIC :: init => key_value_store_init
    PROCEDURE, PUBLIC :: destruct => key_value_store_destruct
    PROCEDURE, PUBLIC :: bcast => key_value_store_bcast
    PROCEDURE, PUBLIC :: getEntryCount => key_value_store_getEntryCount
    PROCEDURE, PUBLIC :: getIterator => key_value_store_getIterator
  END TYPE t_key_value_store

  CHARACTER(*), PARAMETER :: modname = "mo_key_value_store"
  INTEGER, PARAMETER :: dt_txt = 1, dt_int = 2, dt_log = 3, dt_flt = 4

CONTAINS

  SUBROUTINE key_value_store_init(me, cs)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    LOGICAL, INTENT(IN) :: cs
    
    IF (me%is_init) THEN
      IF (cs .NEQV. me%lcase_sensitive .OR. me%table%getEntryCount() .GT. 0) &
        & CALL key_value_store_destruct(me)
    END IF
    IF (.NOT. me%is_init) THEN
      me%lcase_sensitive = cs
      me%table => hashTable_make(text_hash, text_isEqual)
      me%is_init = .true.
    END IF
  END SUBROUTINE key_value_store_init

  SUBROUTINE key_value_store_bcast(me, root, comm)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: root, comm
    LOGICAL :: isRoot, isDest
    TYPE(t_PackedMessage) :: pmsg

    CALL p_get_bcast_role(root, comm, isRoot, isDest)
    IF (isRoot) CALL key_value_store_output(me, pmsg=pmsg)
    CALL pmsg%bcast(root, comm)
    IF (isDest) CALL key_value_store_input_pmsg(me, pmsg)
  END SUBROUTINE key_value_store_bcast

#ifdef __PGI
#define ALLOC_ASSIGN(__a) \
    __a , POINTER :: val_p; CLASS(*), POINTER :: valObj; \
    ALLOCATE(val_p); val_p = val; valObj => val_p
#define ALLOC_ASSIGN_C \
    TYPE(t_char_workaround), POINTER :: val_p; CLASS(*), POINTER :: valObj; \
    ALLOCATE(val_p); ALLOCATE(CHARACTER(LEN=LEN_TRIM(val)) :: val_p%c); \
    WRITE(val_p%c, "(a)") TRIM(val); valObj => val_p
#else
#define ALLOC_ASSIGN(__a) \
    CLASS(*), POINTER :: valObj; ALLOCATE(valObj , SOURCE=val)
#define ALLOC_ASSIGN_C ALLOC_ASSIGN(xxxx)
#endif

  SUBROUTINE key_value_store_put_internal(me, key, valObj)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    CLASS(*), POINTER, INTENT(IN) :: valObj

    IF (.NOT.me%is_init) CALL finish(modname//"put", "data structure not initialized")
    IF (me%lcase_sensitive) THEN
      CALL put_actual(TRIM(key))
    ELSE
      CALL put_actual(tolower(TRIM(key)))
    END IF
  CONTAINS

    SUBROUTINE put_actual(key)
      CHARACTER(*), INTENT(IN) :: key
      CLASS(*), POINTER :: keyObj
#ifdef __PGI
      TYPE(t_char_workaround), POINTER :: key_p

      ALLOCATE(key_p)
      ALLOCATE(CHARACTER(LEN=LEN(key)) :: key_p%c)
      WRITE(key_p%c, "(a)") key
      keyObj => key_p
#else
      ALLOCATE(keyObj , SOURCE=key)
#endif
      CALL me%table%setEntry(keyObj, valObj)
    END SUBROUTINE put_actual
  END SUBROUTINE key_value_store_put_internal

  SUBROUTINE key_value_store_put_c(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key, val

    ALLOC_ASSIGN_C
    CALL key_value_store_put_internal(me, key, valObj)
  END SUBROUTINE key_value_store_put_c

  SUBROUTINE key_value_store_put_r(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    REAL(wp), INTENT(IN) :: val

    ALLOC_ASSIGN(REAL(wp))
    CALL key_value_store_put_internal(me, key, valObj)
  END SUBROUTINE key_value_store_put_r

  SUBROUTINE key_value_store_put_i(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(IN) :: val

    ALLOC_ASSIGN(INTEGER)
    CALL key_value_store_put_internal(me, key, valObj)
  END SUBROUTINE key_value_store_put_i

  SUBROUTINE key_value_store_put_l(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL, INTENT(IN) :: val

    ALLOC_ASSIGN(LOGICAL)
    CALL key_value_store_put_internal(me, key, valObj)
  END SUBROUTINE key_value_store_put_l

  SUBROUTINE key_value_store_get_internal(me, trimmed_key, opt_err, opt_c, opt_r, opt_i, opt_l)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: trimmed_key
    CHARACTER(:), INTENT(OUT), ALLOCATABLE, OPTIONAL :: opt_c
    REAL(wp), INTENT(OUT), OPTIONAL :: opt_r
    INTEGER, INTENT(OUT), OPTIONAL :: opt_i, opt_err
    LOGICAL, INTENT(OUT), OPTIONAL :: opt_l
    CHARACTER(*), PARAMETER :: routine = modname//"get"
    CLASS(*), POINTER :: valObj

    IF (.NOT.me%is_init) CALL finish(routine, "data structure not initialized")
    IF (me%lcase_sensitive) THEN
      valObj => get_valObj(trimmed_key)
    ELSE
      valObj => get_valObj(tolower(trimmed_key))
    END IF
    IF(ASSOCIATED(valObj)) THEN
      IF (PRESENT(opt_err)) opt_err = 0
      SELECT TYPE(valObj)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        IF (.NOT.PRESENT(opt_c)) &
          & CALL finish(routine, "type mismatch '"//trimmed_key//"', isCHARACTER")
        opt_c = valObj%c
#else
      TYPE IS(CHARACTER(*))
        IF (.NOT.PRESENT(opt_c)) &
          &CALL finish(routine, "type mismatch '"//trimmed_key//"', isCHARACTER")
        opt_c = valObj
#endif
      TYPE IS(REAL(wp))
         IF (.NOT.PRESENT(opt_r)) &
           & CALL finish(routine, "type mismatch '"//trimmed_key//"', is REAL(wp)")
        opt_r = valObj
      TYPE IS(INTEGER)
        IF (.NOT.PRESENT(opt_i)) &
          & CALL finish(routine, "type mismatch '"//trimmed_key//"', is INTEGER")
        opt_i = valObj
      TYPE IS(LOGICAL)
        IF (.NOT.PRESENT(opt_l)) &
          & CALL finish(routine, "type mismatch '"//trimmed_key//"', is LOGICAL")
        opt_l = valObj
      CLASS DEFAULT
        CALL finish(routine, "datatype not supported")
      END SELECT
    ELSE
      IF (PRESENT(opt_err)) THEN
        opt_err = 1
      ELSE
        CALL key_value_store_output(me)
        CALL finish(routine, "key '"//trimmed_key//"' not found")
      END IF
    END IF
  CONTAINS

    FUNCTION get_valObj(key) RESULT(valObj)
      CHARACTER(*), INTENT(IN), TARGET :: key
      CLASS(*), POINTER :: keyObj, valObj

      keyObj => key
      valObj => me%table%getEntry(keyObj)
    END FUNCTION get_valObj
  END SUBROUTINE key_value_store_get_internal

  SUBROUTINE key_value_store_get_c(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: res
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_c = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_c

  SUBROUTINE key_value_store_get_r(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    REAL(wp), INTENT(OUT) :: res
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_r = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_r

  SUBROUTINE key_value_store_get_i(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: res
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_i = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_i

  SUBROUTINE key_value_store_get_l(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL, INTENT(OUT) :: res
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_l = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_l

  SUBROUTINE key_value_store_input_pmsg(me, pmsg)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    TYPE(t_PackedMessage), INTENT(INOUT) :: pmsg
    INTEGER :: natt, i, attType, attLen, namLen, alnLen, alvLen, aInteger
    REAL(wp) :: aDouble
    LOGICAL :: aLogical, lcase_sensitive
    CHARACTER(:), ALLOCATABLE :: attTxt, attName

    alvLen = -1
    alnLen = -1
    CALL pmsg%unpack(lcase_sensitive)
    CALL key_value_store_init(me, lcase_sensitive)
    CALL pmsg%unpack(natt)
    DO i = 1, natt
      CALL pmsg%unpack(namLen)
      IF (namLen .NE. alnLen) THEN
        alnLen = namLen
        IF (ALLOCATED(attName)) DEALLOCATE(attName)
        ALLOCATE(CHARACTER(alnLen) :: attName)
      END IF
      CALL pmsg%unpack(attName)
      CALL pmsg%unpack(attType)
      SELECT CASE(attType)
      CASE(dt_flt)
        CALL pmsg%unpack(aDouble)
        CALL key_value_store_put_r(me, attName, aDouble)
      CASE(dt_int)
        CALL pmsg%unpack(aInteger)
        CALL key_value_store_put_i(me, attName, aInteger)
      CASE(dt_log)
        CALL pmsg%unpack(aLogical)
        CALL key_value_store_put_l(me, attName, aLogical)
      CASE(dt_txt)
        CALL pmsg%unpack(attLen)
        IF (attLen .NE. alvLen) THEN
          alvLen = attLen
          IF (ALLOCATED(attTxt)) DEALLOCATE(attTxt)
          ALLOCATE(CHARACTER(alvLen) :: attTxt)
        END IF
        CALL pmsg%unpack(attTxt)
        CALL key_value_store_put_c(me, attName, attTxt)
      END SELECT
    ENDDO
  END SUBROUTINE key_value_store_input_pmsg

  SUBROUTINE key_value_store_output(me, pmsg, archive, copy, inverse, label)
    CLASS(t_key_value_store), INTENT(IN) :: me
    TYPE(t_PackedMessage), INTENT(INOUT), OPTIONAL :: pmsg
    TYPE(t_key_value_store), INTENT(INOUT), OPTIONAL :: archive, copy, inverse
    CHARACTER(*), INTENT(IN), OPTIONAL :: label
    CHARACTER(*), PARAMETER :: routine = modname//"output"
    TYPE(t_HashIterator) :: iterator
    CLASS(*), POINTER :: curKey, curVal
    CHARACTER(:), POINTER :: ccKey, ccVal
    LOGICAL :: selector(4)

    selector = [PRESENT(pmsg), PRESENT(archive), PRESENT(copy), PRESENT(inverse)]
    IF (COUNT(selector(:)) .GT. 1) CALL finish(routine, "too many arguments")
    IF (selector(1)) THEN
      IF (.NOT.me%is_init) RETURN
      CALL pmsg%pack(me%lcase_sensitive)
      CALL pmsg%pack(me%table%getEntryCount())
    ELSE IF (selector(2)) THEN
      CALL key_value_store_init(archive, me%lcase_sensitive)
    ELSE IF (selector(3)) THEN
      CALL key_value_store_init(copy, me%lcase_sensitive)
    ELSE IF (selector(4)) THEN
      CALL key_value_store_init(inverse, me%lcase_sensitive)
    ELSE IF (.NOT.ANY(selector) .AND. PRESENT(label)) THEN
      CALL message(routine,"START STORAGE DUMP (" // TRIM(label) // ")")
      IF (me%getEntryCount() .EQ. 0) CALL message(routine,"structure is empty")
    END IF
    CALL iterator%init(me%table)
    DO WHILE (iterator%nextEntry(curKey, curVal))
      ccKey => sel_char(curKey, routine, "key: invalid type")
      IF (selector(1)) THEN
        CALL output_pmsg()
      ELSE IF (selector(2)) THEN
        IF (ccKey(1:MIN(4,LEN(ccKey))) /= 'nml_') CYCLE
        ccVal => sel_char(curVal, routine, "nml-archive copy: only character valued pairs!")
        CALL key_value_store_put_c(archive, ccKey, ccVal)
      ELSE IF (selector(3)) THEN
        ccVal => sel_char(curVal, routine, "copy: only character valued pairs!")
        CALL key_value_store_put_c(copy, ccKey, ccVal)
      ELSE IF (selector(4)) THEN
        ccVal => sel_char(curVal, routine, "inversion: only character valued pairs!")
        CALL key_value_store_put_c(inverse, ccVal, ccKey)
      ELSE
        CALL output_stdio()
      END IF
    END DO
    IF (.NOT.ANY(selector) .AND. PRESENT(label)) &
      & CALL message(routine,"END OF STORAGE DUMP (" // label // ")")
  CONTAINS

    SUBROUTINE output_pmsg()
      INTEGER :: clen
  
      clen = LEN_TRIM(ccKey)
      CALL pmsg%pack(clen)
      CALL pmsg%pack(ccKey(1:clen))
      SELECT TYPE(curVal)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        ccVal => curVal%c
#else
      TYPE IS(CHARACTER(*))
        ccVal => curVal
#endif
        CALL pmsg%pack(dt_txt)
        clen = LEN_TRIM(ccVal)
        CALL pmsg%pack(clen)
        CALL pmsg%pack(ccVal(1:clen))
      TYPE IS(REAL(wp))
        CALL pmsg%pack(dt_flt)
        CALL pmsg%pack(curVal)
      TYPE IS(INTEGER)
        CALL pmsg%pack(dt_int)
        CALL pmsg%pack(curVal)
      TYPE IS(LOGICAL)
        CALL pmsg%pack(dt_log)
        CALL pmsg%pack(curVal)
      CLASS DEFAULT
        CALL finish(routine, "val: invalid type")
      END SELECT
    END SUBROUTINE output_pmsg

    SUBROUTINE output_stdio()
      CHARACTER(LEN=576) :: message_text ! message_text from mo_exception is too short

      message_text = ''
      NULLIFY(ccVal)
      SELECT TYPE(curVal)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        ccVal => curVal%c
#else
      TYPE IS(CHARACTER(*))
        ccVal => curVal
#endif
      TYPE IS(REAL(wp))
        WRITE(message_text, "(3a,e12.5,a)") "key = >", ccKey, "< val = >", curVal, "<"
      TYPE IS(INTEGER)
        WRITE(message_text, "(3a,i6,a)") "key = >", ccKey, "< val = >", curVal, "<"
      TYPE IS(LOGICAL)
        WRITE(message_text, "(3a,l1,a)") "key = >", ccKey, "< val = >", curVal, "<"
      CLASS DEFAULT
        CALL finish(routine, "val: invalid type")
      END SELECT
      IF (ASSOCIATED(ccVal)) THEN
        IF (LEN(ccVal) .GT. 448) THEN
          WRITE(message_text, "(5a)") "key = >", ccKey, "< val = >", ccVal(1:400), "< !TRUNCATED!"
        ELSE
          WRITE(message_text, "(5a)") "key = >", ccKey, "< val = >", ccVal, "<"
        END IF
      END IF
      IF (PRESENT(label)) THEN
        CALL message(routine, message_text)
      ELSE
        PRINT "(a)", TRIM(message_text)
      END IF
    END SUBROUTINE output_stdio
  END SUBROUTINE key_value_store_output

  INTEGER FUNCTION key_value_store_getEntryCount(me) RESULT(n)
    CLASS(t_key_value_store), INTENT(IN) :: me

    n = me%table%getEntryCount()
  END FUNCTION key_value_store_getEntryCount

  FUNCTION key_value_store_getIterator(me) RESULT(iterator)
    CLASS(t_key_value_store), INTENT(IN) :: me
    TYPE(t_HashIterator), POINTER :: iterator

    ALLOCATE(iterator)
    CALL iterator%init(me%table)
  END FUNCTION key_value_store_getIterator

  SUBROUTINE key_value_store_destruct(me)
    CLASS(t_key_value_store), INTENT(INOUT) :: me

    IF (ASSOCIATED(me%table)) THEN
      CALL me%table%destruct()
      DEALLOCATE(me%table)
    END IF
    me%lcase_sensitive = .true.
    me%is_init = .false.
  END SUBROUTINE key_value_store_destruct

END MODULE mo_key_value_store
