!> @file comin_metadata_types.F90
!! @brief Fortran interface to C++ metadata container
!
!  @authors 08/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_metadata_types

  USE iso_c_binding,         ONLY: c_int, c_char, c_double, c_bool, &
                               &   c_ptr, c_loc, c_null_ptr, c_associated
  USE comin_c_utils,         ONLY: convert_c_string, convert_f_string
  USE comin_setup_constants, ONLY: COMIN_ZAXIS_3D, COMIN_HGRID_UNSTRUCTURED_CELL

  PRIVATE

  PUBLIC :: t_comin_var_metadata, t_comin_var_metadata_iterator
  ! The following functions are used to provide a more direct access
  !   for C and Python plugins to the underlying keyval container
  PUBLIC :: comin_keyval_get_char_c
  PUBLIC :: comin_keyval_iterator_begin_c, comin_keyval_iterator_end_c
  PUBLIC :: comin_keyval_iterator_get_key_c
  PUBLIC :: comin_keyval_iterator_compare_c
  PUBLIC :: comin_keyval_iterator_next_c
  PUBLIC :: comin_keyval_iterator_delete_c

#include "comin_global.inc"

  TYPE t_comin_var_metadata_iterator
    TYPE(c_ptr) :: comin_metadata_iterator_current_c = c_null_ptr
    TYPE(c_ptr) :: comin_metadata_iterator_end_c     = c_null_ptr
  CONTAINS
    PROCEDURE :: next    => comin_metadata_iterator_next
    PROCEDURE :: is_end  => comin_metadata_iterator_is_end
    PROCEDURE :: key     => comin_metadata_get_key_from_iterator
    PROCEDURE :: delete  => comin_metadata_iterator_delete
  END TYPE t_comin_var_metadata_iterator

  TYPE t_comin_var_metadata
    TYPE(c_ptr) :: comin_metadata_c = c_null_ptr
  CONTAINS
    PROCEDURE :: create              => comin_metadata_create
    PROCEDURE :: delete              => comin_metadata_delete
    PROCEDURE :: query               => comin_metadata_query
    PROCEDURE :: get_iterator        => comin_metadata_get_iterator
    GENERIC, PUBLIC    :: set        => set_int, set_double, set_char, set_bool
    PROCEDURE, PRIVATE :: set_int    => comin_metadata_set_int
    PROCEDURE, PRIVATE :: set_double => comin_metadata_set_double
    PROCEDURE, PRIVATE :: set_char   => comin_metadata_set_char
    PROCEDURE, PRIVATE :: set_bool   => comin_metadata_set_bool
    GENERIC, PUBLIC    :: get        => get_int, get_double, get_char, get_bool
    PROCEDURE, PRIVATE :: get_int    => comin_metadata_get_int
    PROCEDURE, PRIVATE :: get_double => comin_metadata_get_double
    PROCEDURE, PRIVATE :: get_char   => comin_metadata_get_char
    PROCEDURE, PRIVATE :: get_bool   => comin_metadata_get_bool
  END TYPE t_comin_var_metadata

  INTERFACE
    SUBROUTINE comin_keyval_create_c(comin_keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr) :: comin_keyval_c
    END SUBROUTINE comin_keyval_create_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_delete_c(comin_keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: comin_keyval_c
    END SUBROUTINE comin_keyval_delete_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_query_c(ckey, idx, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_ptr, c_int
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      INTEGER(KIND=c_int), INTENT(out) :: idx
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_query_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_set_int_c(ckey,  val, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_int, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      INTEGER(KIND=c_int), VALUE, INTENT(in) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_int_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_get_int_c(ckey, val, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_int, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      INTEGER(KIND=c_int), INTENT(out) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_int_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_set_double_c(ckey,  val, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_double, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      REAL(KIND=c_double), VALUE, INTENT(in) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_double_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_get_double_c(ckey, val, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_double, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      REAL(KIND=c_double), INTENT(out) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_double_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_set_char_c(ckey,  cval, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      TYPE(c_ptr), VALUE, INTENT(in) :: cval
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_char_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_get_char_c(ckey,  cval, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      TYPE(c_ptr) :: cval
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_char_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_set_bool_c(ckey,  val, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_bool, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_set_bool_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_get_bool_c(ckey, val, keyval_c) BIND(C)
      USE iso_c_binding, ONLY: c_bool, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: ckey
      LOGICAL(KIND=c_bool), INTENT(out) :: val
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
    END SUBROUTINE comin_keyval_get_bool_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_iterator_begin_c(keyval_c, iterator) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
      TYPE(c_ptr) :: iterator
    END SUBROUTINE comin_keyval_iterator_begin_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_iterator_end_c(keyval_c, iterator) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: keyval_c
      TYPE(c_ptr) :: iterator
    END SUBROUTINE comin_keyval_iterator_end_c
  END INTERFACE

  INTERFACE
    FUNCTION comin_keyval_iterator_get_key_c(iterator) RESULT(ckey) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator
      TYPE(c_ptr) :: ckey
    END FUNCTION comin_keyval_iterator_get_key_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_iterator_next_c(iterator) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator
    END SUBROUTINE comin_keyval_iterator_next_c
  END INTERFACE

  INTERFACE
    SUBROUTINE comin_keyval_iterator_delete_c(iterator) BIND(C)
      USE iso_c_binding, ONLY: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator
    END SUBROUTINE comin_keyval_iterator_delete_c
  END INTERFACE

  INTERFACE
    FUNCTION comin_keyval_iterator_compare_c(iterator1, iterator2) BIND(C)
      USE iso_c_binding, ONLY: c_ptr, c_bool
      TYPE(c_ptr), VALUE, INTENT(in) :: iterator1, iterator2
      LOGICAL(KIND=c_bool) :: comin_keyval_iterator_compare_c
    END FUNCTION comin_keyval_iterator_compare_c
  END INTERFACE

CONTAINS

  SUBROUTINE comin_metadata_create(this)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    ! If metadata container is already created, simply ignore this step
    ! (entries can be still be added/modified via set in that case)
    IF ( .NOT. c_associated(this%comin_metadata_c) ) THEN
      CALL comin_keyval_create_c(this%comin_metadata_c)
    ENDIF
  END SUBROUTINE comin_metadata_create

  SUBROUTINE comin_metadata_delete(this)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CALL comin_keyval_delete_c(this%comin_metadata_c)
  END SUBROUTINE comin_metadata_delete

  INTEGER FUNCTION comin_metadata_query(this, key)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    INTEGER(KIND=c_int) :: idx
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)

    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_query_c(C_LOC(ckey), idx, this%comin_metadata_c)

    comin_metadata_query = INT(idx)
  END FUNCTION comin_metadata_query

  SUBROUTINE comin_metadata_set_int(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    INTEGER(KIND=c_int), INTENT(in)  :: val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)

    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_set_int_c(C_LOC(ckey), val, this%comin_metadata_c)

  END SUBROUTINE comin_metadata_set_int

  SUBROUTINE comin_metadata_get_int(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    INTEGER(KIND=c_int), INTENT(out)  :: val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)

    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_get_int_c(C_LOC(ckey), val, this%comin_metadata_c)
  END SUBROUTINE comin_metadata_get_int

  SUBROUTINE comin_metadata_set_double(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    REAL(KIND=c_double), INTENT(in)  :: val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)

    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_set_double_c(C_LOC(ckey), val, this%comin_metadata_c)

  END SUBROUTINE comin_metadata_set_double

  SUBROUTINE comin_metadata_get_double(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    REAL(KIND=c_double), INTENT(out)  :: val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)

    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_get_double_c(C_LOC(ckey), val, this%comin_metadata_c)
  END SUBROUTINE comin_metadata_get_double

  SUBROUTINE comin_metadata_set_char(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key, val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1), cval(LEN(TRIM(val))+1)

    CALL convert_f_string(TRIM(key), ckey)
    CALL convert_f_string(TRIM(val), cval)
    CALL comin_keyval_set_char_c(C_LOC(ckey), C_LOC(cval), this%comin_metadata_c)

  END SUBROUTINE comin_metadata_set_char

  SUBROUTINE comin_metadata_get_char(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)
    TYPE(C_PTR) :: cval

    cval = c_null_ptr
    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_get_char_c(C_LOC(ckey), cval, this%comin_metadata_c)
    val = convert_c_string(cval)
  END SUBROUTINE comin_metadata_get_char

  SUBROUTINE comin_metadata_set_bool(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    LOGICAL, INTENT(in) :: val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)

    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_set_bool_c(C_LOC(ckey), LOGICAL(val,c_bool), this%comin_metadata_c)

  END SUBROUTINE comin_metadata_set_bool

  SUBROUTINE comin_metadata_get_bool(this, key, val)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    CHARACTER(LEN=*), INTENT(in) :: key
    LOGICAL, INTENT(out)  :: val
    CHARACTER(len=1, KIND=c_char), TARGET :: ckey(LEN(TRIM(key))+1)
    LOGICAL(KIND=c_bool) :: cval

    CALL convert_f_string(TRIM(key), ckey)
    CALL comin_keyval_get_bool_c(C_LOC(ckey), cval, this%comin_metadata_c)
    val = LOGICAL(cval)
  END SUBROUTINE comin_metadata_get_bool

  SUBROUTINE comin_metadata_get_iterator(this, iterator)
    CLASS(t_comin_var_metadata), INTENT(inout) :: this
    TYPE(t_comin_var_metadata_iterator), INTENT(inout) :: iterator

    IF (c_associated(iterator%comin_metadata_iterator_current_c)) CALL iterator%delete()
    CALL comin_keyval_iterator_begin_c(this%comin_metadata_c, iterator%comin_metadata_iterator_current_c)
    CALL comin_keyval_iterator_end_c  (this%comin_metadata_c, iterator%comin_metadata_iterator_end_c    )
  END SUBROUTINE comin_metadata_get_iterator

  FUNCTION comin_metadata_get_key_from_iterator(this) RESULT(key)
    CLASS(t_comin_var_metadata_iterator), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: key
    TYPE(C_PTR) :: ckey

    ckey = comin_keyval_iterator_get_key_c(this%comin_metadata_iterator_current_c)
    key = convert_c_string(ckey)
  END FUNCTION comin_metadata_get_key_from_iterator

  SUBROUTINE comin_metadata_iterator_next(this)
    CLASS(t_comin_var_metadata_iterator), INTENT(in) :: this

    CALL comin_keyval_iterator_next_c(this%comin_metadata_iterator_current_c)
  END SUBROUTINE comin_metadata_iterator_next

  SUBROUTINE comin_metadata_iterator_delete(this)
    CLASS(t_comin_var_metadata_iterator), INTENT(inout) :: this

    CALL comin_keyval_iterator_delete_c(this%comin_metadata_iterator_current_c)
    CALL comin_keyval_iterator_delete_c(this%comin_metadata_iterator_end_c    )
    this%comin_metadata_iterator_current_c = c_null_ptr ! Reset c_ptr such that c_associated can be used
    this%comin_metadata_iterator_end_c     = c_null_ptr ! Reset c_ptr such that c_associated can be used
  END SUBROUTINE

  LOGICAL FUNCTION comin_metadata_iterator_is_end(this)
    CLASS(t_comin_var_metadata_iterator), INTENT(inout) :: this

    comin_metadata_iterator_is_end = LOGICAL(comin_keyval_iterator_compare_c( &
      &  this%comin_metadata_iterator_current_c, this%comin_metadata_iterator_end_c) )
  END FUNCTION comin_metadata_iterator_is_end

END MODULE comin_metadata_types
