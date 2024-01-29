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

MODULE mo_util_texthash

  USE ISO_C_BINDING, ONLY: c_int, c_size_t, c_char, c_int32_t
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: text_hash, text_hash_c, text_isEqual, sel_char

#if defined(__PGI) || defined(__FLANG)
  TYPE, PUBLIC :: t_char_workaround
    CHARACTER(:), ALLOCATABLE :: c
  END TYPE t_char_workaround
#endif

  CHARACTER(*), PARAMETER :: modname = "mo_util_texthash"

CONTAINS

  FUNCTION sel_char(key, routine, err_msg) RESULT(ptr)
    CLASS(*), POINTER, INTENT(in) :: key
    CHARACTER(*), INTENT(IN) :: routine, err_msg
    CHARACTER(:), POINTER :: ptr

    SELECT TYPE (key)
#if defined(__PGI) || defined(__FLANG)
    TYPE IS (t_char_workaround)
      ptr => key%c
#endif
    TYPE IS (CHARACTER(*))
      ptr => key
    CLASS DEFAULT
      CALL finish(routine, err_msg)
    END SELECT
  END FUNCTION sel_char

  INTEGER FUNCTION text_hash_c(key) RESULT(hash)
    CHARACTER(*), INTENT(IN), TARGET :: key
    CLASS(*), POINTER :: key_p

    key_p => key
    hash = INT(text_hash(key_p))
  END FUNCTION text_hash_c

  INTEGER(c_int) FUNCTION text_hash(key) RESULT(hash)
    INTERFACE
      FUNCTION util_hashword(text, text_len, inithash) RESULT(hash) BIND(C, NAME='util_hashword')
        IMPORT :: c_char, c_size_t, c_int32_t
        INTEGER(c_int32_t) :: hash
        CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: text
        INTEGER(kind=c_size_t), VALUE, INTENT(IN) :: text_len
        INTEGER(c_int32_t), VALUE, INTENT(IN) :: inithash
      END FUNCTION util_hashword
    END INTERFACE

    CLASS(*), POINTER, INTENT(IN) :: key
    CHARACTER(*), PARAMETER :: routine = modname//":text_hash_cs"
    CHARACTER(:), POINTER :: key_p

    key_p => sel_char(key, routine, "Unknown type for key.")
    hash = INT(util_hashword(key_p, INT(LEN(key_p), c_size_t), 0_c_int32_t), c_int)
  END FUNCTION text_hash

  LOGICAL FUNCTION text_isEqual(keyA, keyB) RESULT(is_equal)
    CLASS(*), POINTER, INTENT(in) :: keyA, keyB
    CHARACTER(*), PARAMETER :: routine = modname//":text_isEqual_cs"
    CHARACTER(:), POINTER :: keyA_p, keyB_p

    keyA_p => sel_char(keyA, routine, "Unknown type for keyA.")
    keyB_p => sel_char(keyB, routine, "Unknown type for keyB.")
    is_equal = keyA_p == keyB_p
  END FUNCTION text_isEqual
END MODULE mo_util_texthash
