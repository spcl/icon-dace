!> @file comin_metadata.F90
!! @brief Variable metadata definition.
!
!  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_metadata

  USE iso_c_binding,           ONLY: c_int, c_ptr, c_bool, c_double, c_loc
  USE comin_setup_constants,   ONLY: wp
  USE comin_errhandler_constants, ONLY: COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED,                 &
    &                                   COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR,  &
    &                                   COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR, &
    &                                   COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND,                &
    &                                   COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS,      &
    &                                   COMIN_ERROR_METADATA_KEY_NOT_FOUND,                  &
    &                                   COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE
  USE comin_errhandler,        ONLY: comin_message, comin_error_set, comin_plugin_finish
  USE comin_state,             ONLY: state
  USE comin_c_utils,           ONLY: convert_c_string
  USE comin_variable_types,    ONLY: t_comin_var_descriptor, t_comin_var_item,            &
    &                                t_var_request_list_item, t_comin_var_descriptor_c,   &
    &                                comin_var_descr_match
  USE comin_variable,          ONLY: comin_var_get_from_exposed
  USE comin_descrdata_types,   ONLY: t_comin_descrdata_global
  USE comin_descrdata,         ONLY: comin_descrdata_get_global
  USE comin_metadata_types,    ONLY: t_comin_var_metadata, t_comin_var_metadata_iterator, &
    &                                comin_keyval_get_char_c, &
    &                                comin_keyval_iterator_begin_c, comin_keyval_iterator_end_c, &
    &                                comin_keyval_iterator_get_key_c, comin_keyval_iterator_compare_c, &
    &                                comin_keyval_iterator_next_c, comin_keyval_iterator_delete_c
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_metadata_set_request
  PUBLIC :: comin_metadata_set_host
  PUBLIC :: comin_metadata_get
  PUBLIC :: comin_metadata_get_or
  PUBLIC :: comin_metadata_get_typeid
  PUBLIC :: comin_metadata_get_iterator
  PUBLIC :: comin_metadata_set_integer_c, comin_metadata_set_real_c
  PUBLIC :: comin_metadata_set_character_c, comin_metadata_set_logical_c
  PUBLIC :: comin_metadata_get_integer_c, comin_metadata_get_real_c
  PUBLIC :: comin_metadata_get_character_c, comin_metadata_get_logical_c

#include "comin_global.inc"

  !> Sets metadata for a requested ComIn variable.
  !> **Note:Plugins use the alias `comin_metadata_set`.**
  !! @ingroup plugin_interface
  INTERFACE comin_metadata_set_request
    MODULE PROCEDURE comin_request_set_var_metadata_logical
    MODULE PROCEDURE comin_request_set_var_metadata_integer
    MODULE PROCEDURE comin_request_set_var_metadata_real
    MODULE PROCEDURE comin_request_set_var_metadata_character
  END INTERFACE comin_metadata_set_request

  !> Sets metadata for an exposed variable.
  !> **Note: The host model uses the alias `comin_metadata_set`.**
  !! @ingroup host_interface
  INTERFACE comin_metadata_set_host
    MODULE PROCEDURE comin_metadata_host_set_logical
    MODULE PROCEDURE comin_metadata_host_set_integer
    MODULE PROCEDURE comin_metadata_host_set_real
    MODULE PROCEDURE comin_metadata_host_set_character
  END INTERFACE comin_metadata_set_host

  !> Read-only access to additional information about a given variable.
  !! @ingroup plugin_interface
  INTERFACE comin_metadata_get
    MODULE PROCEDURE comin_metadata_get_logical
    MODULE PROCEDURE comin_metadata_get_integer
    MODULE PROCEDURE comin_metadata_get_real
    MODULE PROCEDURE comin_metadata_get_character
  END INTERFACE comin_metadata_get

  INTERFACE comin_metadata_get_or
    MODULE PROCEDURE comin_metadata_get_or_integer
    MODULE PROCEDURE comin_metadata_get_or_real
    MODULE PROCEDURE comin_metadata_get_or_character
    MODULE PROCEDURE comin_metadata_get_or_logical
  END INTERFACE

CONTAINS

  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_integer(descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    INTEGER,                      INTENT(IN)  :: val !< metadata value
    !
    TYPE(t_comin_var_item), POINTER :: var_item

    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    ELSE
      IF ( ALL(var_item%metadata%query(TRIM(key)) /= (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_INTEGER/) ) ) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
      ENDIF
      CALL var_item%metadata%set(key, val)
    END IF
  END SUBROUTINE comin_metadata_host_set_integer

  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_logical(descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL,                      INTENT(IN)  :: val !< metadata value
    !
    TYPE(t_comin_var_item), POINTER :: var_item

    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    ELSE
      IF ( ALL(var_item%metadata%query(TRIM(key)) /= (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_LOGICAL/) ) ) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
      ENDIF
      CALL var_item%metadata%set(key, val)
    END IF
  END SUBROUTINE comin_metadata_host_set_logical

  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_real(descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    REAL(wp),                     INTENT(IN)  :: val !< metadata value
    !
    TYPE(t_comin_var_item), POINTER :: var_item

    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    ELSE
      IF ( ALL(var_item%metadata%query(TRIM(key)) /= (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_REAL/) ) ) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
      ENDIF
      CALL var_item%metadata%set(key, val)
    END IF
  END SUBROUTINE comin_metadata_host_set_real

  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_character(descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    CHARACTER(LEN=*),             INTENT(IN)  :: val !< metadata value
    !
    TYPE(t_comin_var_item), POINTER :: var_item

    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    ELSE
      IF ( ALL(var_item%metadata%query(TRIM(key)) /= (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_CHARACTER/) ) ) THEN
        CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
      ENDIF
      CALL var_item%metadata%set(key, val)
    END IF
  END SUBROUTINE comin_metadata_host_set_character

  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_integer(var_descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    INTEGER,                      INTENT(OUT) :: val !< metadata value
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF

    IF ( var_item%metadata%query(TRIM(key)) /= COMIN_TYPEID_INTEGER)  THEN
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE)
      RETURN
    ENDIF

    CALL var_item%metadata%get(key, val)
  END SUBROUTINE comin_metadata_get_integer

  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_logical(var_descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL,                      INTENT(OUT) :: val !< metadata value
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF

    IF ( var_item%metadata%query(TRIM(key)) /= COMIN_TYPEID_LOGICAL)  THEN
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE)
      RETURN
    ENDIF

    CALL var_item%metadata%get(key, val)
  END SUBROUTINE comin_metadata_get_logical

  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_real(var_descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    REAL(wp),                     INTENT(OUT) :: val !< metadata value
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF

    IF ( var_item%metadata%query(TRIM(key)) /= COMIN_TYPEID_REAL)  THEN
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE)
      RETURN
    ENDIF

    CALL var_item%metadata%get(key, val)
  END SUBROUTINE comin_metadata_get_real

  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_character(var_descriptor, key, val)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    CHARACTER(LEN=:), ALLOCATABLE,INTENT(OUT) :: val !< metadata value
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF

    IF ( var_item%metadata%query(TRIM(key)) /= COMIN_TYPEID_CHARACTER)  THEN
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE)
      RETURN
    ENDIF

    CALL var_item%metadata%get(key, val)
  END SUBROUTINE comin_metadata_get_character

  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_integer_c(var_descriptor, key, val) &
    & BIND(C, NAME="comin_metadata_get_integer")
    TYPE(t_comin_var_descriptor_c), VALUE,  INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    INTEGER(kind=c_int),            INTENT(OUT) :: val !< metadata value
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    INTEGER :: val_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_metadata_get_integer(var_descriptor_fortran, &
      &  convert_c_string(key), val_fortran)
    val = INT(val_fortran, c_int)
  END SUBROUTINE comin_metadata_get_integer_c

  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_logical_c(var_descriptor, key, val) &
    & BIND(C, NAME="comin_metadata_get_logical")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL(kind=c_bool),           INTENT(OUT) :: val !< metadata value
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    LOGICAL :: val_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_metadata_get_logical(var_descriptor_fortran, &
      &  convert_c_string(key), val_fortran)
    val = LOGICAL(val_fortran, c_bool)
  END SUBROUTINE comin_metadata_get_logical_c

  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_real_c(var_descriptor, key, val) &
    & BIND(C, NAME="comin_metadata_get_real")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    REAL(kind=c_double),            INTENT(OUT) :: val !< metadata value
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    REAL(wp) :: val_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_metadata_get_real(var_descriptor_fortran, &
      &  convert_c_string(key), val_fortran)
    val = REAL(val_fortran, c_double)
  END SUBROUTINE comin_metadata_get_real_c

  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_character_c(var_descriptor, key, val, len) &
    & BIND(C, NAME="comin_metadata_get_character")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    TYPE(c_ptr),                    INTENT(OUT) :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: len !< string length
    ! local
    TYPE(t_comin_var_descriptor)    :: var_descriptor_fortran
    TYPE(t_comin_var_item), POINTER :: var_item

    INTERFACE
      FUNCTION c_strlen(str_ptr) BIND ( C, name = "strlen" ) RESULT(len)
        USE, INTRINSIC :: iso_c_binding
        TYPE(c_ptr), VALUE      :: str_ptr
        INTEGER(kind=c_size_t)  :: len
      END FUNCTION c_strlen
    END INTERFACE

    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF
    ! Create fortran var descriptor and get var_item
    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    var_item => comin_var_get_from_exposed(var_descriptor_fortran)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF

    IF ( var_item%metadata%query(convert_c_string(key)) /= COMIN_TYPEID_CHARACTER)  THEN
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE)
      RETURN
    ENDIF

    ! Not another Fortran detour, get the c_ptr directly
    CALL comin_keyval_get_char_c(key,  val, var_item%metadata%comin_metadata_c)
    len = int(c_strlen(val),c_int)
  END SUBROUTINE comin_metadata_get_character_c

  SUBROUTINE comin_metadata_get_or_integer(metadata, key, val, defaultval)
    TYPE(t_comin_var_metadata), INTENT(inout) :: metadata
    CHARACTER(LEN=*), INTENT(in) :: key
    INTEGER, INTENT(out) :: val
    INTEGER, INTENT(in) :: defaultval

    SELECT CASE ( metadata%query(key) )
    CASE (COMIN_TYPEID_INTEGER)
      CALL metadata%get(key, val)
    CASE (COMIN_TYPEID_UNDEFINED)
      val = defaultval
    CASE DEFAULT
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
    END SELECT
  END SUBROUTINE comin_metadata_get_or_integer

  SUBROUTINE comin_metadata_get_or_real(metadata, key, val, defaultval)
    TYPE(t_comin_var_metadata), INTENT(inout) :: metadata
    CHARACTER(LEN=*), INTENT(in) :: key
    REAL(wp), INTENT(out) :: val
    REAL(wp), INTENT(in) :: defaultval

    SELECT CASE ( metadata%query(key) )
    CASE (COMIN_TYPEID_REAL)
      CALL metadata%get(key, val)
    CASE (COMIN_TYPEID_UNDEFINED)
      val = defaultval
    CASE DEFAULT
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
    END SELECT
  END SUBROUTINE comin_metadata_get_or_real

  SUBROUTINE comin_metadata_get_or_character(metadata, key, val, defaultval)
    TYPE(t_comin_var_metadata), INTENT(inout) :: metadata
    CHARACTER(LEN=*), INTENT(in) :: key
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: val
    CHARACTER(LEN=*), INTENT(in) :: defaultval

    SELECT CASE ( metadata%query(key) )
    CASE (COMIN_TYPEID_CHARACTER)
      CALL metadata%get(key, val)
    CASE (COMIN_TYPEID_UNDEFINED)
      val = defaultval
    CASE DEFAULT
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
    END SELECT
  END SUBROUTINE comin_metadata_get_or_character

  SUBROUTINE comin_metadata_get_or_logical(metadata, key, val, defaultval)
    TYPE(t_comin_var_metadata), INTENT(inout) :: metadata
    CHARACTER(LEN=*), INTENT(in) :: key
    LOGICAL, INTENT(out) :: val
    LOGICAL, INTENT(in) :: defaultval

    SELECT CASE ( metadata%query(key) )
    CASE (COMIN_TYPEID_LOGICAL)
      CALL metadata%get(key, val)
    CASE (COMIN_TYPEID_UNDEFINED)
      val = defaultval
    CASE DEFAULT
      CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
    END SELECT
  END SUBROUTINE comin_metadata_get_or_logical

  SUBROUTINE comin_metadata_set_integer_c(var_descriptor, key, val) &
    &  BIND(C, name="comin_metadata_set_integer")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    INTEGER(kind=c_int), VALUE,     INTENT(IN)  :: val !< metadata value
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_integer(var_descriptor_fortran, &
      &                                         convert_c_string(key), val)
  END SUBROUTINE comin_metadata_set_integer_c

  SUBROUTINE comin_metadata_set_logical_c(var_descriptor, key, val) &
    &  BIND(C, name="comin_metadata_set_logical")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL(C_BOOL), VALUE,         INTENT(IN)  :: val !< metadata value
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_logical(var_descriptor_fortran, &
      &                                         convert_c_string(key), LOGICAL(val))
  END SUBROUTINE comin_metadata_set_logical_c

  SUBROUTINE comin_metadata_set_real_c(var_descriptor, key, val) &
    &  BIND(C, name="comin_metadata_set_real")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    REAL(wp), VALUE,                INTENT(IN)  :: val !< metadata value
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_real(var_descriptor_fortran, &
      &                                      convert_c_string(key), REAL(val, wp))
  END SUBROUTINE comin_metadata_set_real_c

  SUBROUTINE comin_metadata_set_character_c(var_descriptor, key, val) &
    &  BIND(C, name="comin_metadata_set_character")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    TYPE(C_PTR), VALUE,             INTENT(IN)  :: val !< metadata value
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_character(var_descriptor_fortran, &
      &                                           convert_c_string(key), convert_c_string(val))
  END SUBROUTINE comin_metadata_set_character_c

  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_integer(var_descriptor, key, val)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    INTEGER,                       INTENT(IN)  :: val !< metadata value
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global

    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF

    IF (var_descriptor%id == -1) THEN
      comin_global => comin_descrdata_get_global()
      IF (.NOT. ASSOCIATED(comin_global)) CALL comin_plugin_finish("variable ", "global data missing")

      domain_id_start = 1
      domain_id_end   = comin_global%n_dom
    ELSE
      domain_id_start = var_descriptor%id
      domain_id_end   = var_descriptor%id
    END IF

    ! loop over request list
    lfound = .FALSE.
    DO domain_id = domain_id_start, domain_id_end
      var_descriptor_domain    = var_descriptor
      var_descriptor_domain%id = domain_id

      p => state%comin_var_request_list%first()
      DO WHILE (ASSOCIATED(p))
        ASSOCIATE (var_list_request_element => p%item_value)
          IF (comin_var_descr_match(var_list_request_element%descriptor, var_descriptor_domain)) THEN
            lfound = .TRUE.
            IF ( ALL(var_list_request_element%metadata%query(TRIM(key)) /= &
              &      (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_INTEGER/) ) ) THEN
              CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
            ENDIF
            CALL var_list_request_element%metadata%set(key, val)
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND); RETURN
    ENDIF
  END SUBROUTINE comin_request_set_var_metadata_integer

  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_logical(var_descriptor, key, val)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    LOGICAL,                       INTENT(IN)  :: val !< metadata value
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global

    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF

    IF (var_descriptor%id == -1) THEN
      comin_global => comin_descrdata_get_global()
      IF (.NOT. ASSOCIATED(comin_global)) CALL comin_plugin_finish("variable ", "global data missing")

      domain_id_start = 1
      domain_id_end   = comin_global%n_dom
    ELSE
      domain_id_start = var_descriptor%id
      domain_id_end   = var_descriptor%id
    END IF

    ! loop over request list
    lfound = .FALSE.
    DO domain_id = domain_id_start, domain_id_end
      var_descriptor_domain    = var_descriptor
      var_descriptor_domain%id = domain_id

      p => state%comin_var_request_list%first()
      DO WHILE (ASSOCIATED(p))
        ASSOCIATE (var_list_request_element => p%item_value)
          IF (comin_var_descr_match(var_list_request_element%descriptor, var_descriptor_domain)) THEN
            lfound = .TRUE.

            IF ((key == "tracer") .AND. val .AND. (var_descriptor%id /= -1)) THEN
              CALL comin_error_set(COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS); RETURN
            END IF
            IF ( ALL(var_list_request_element%metadata%query(TRIM(key)) /= &
              &      (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_LOGICAL/) ) ) THEN
              CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
            ENDIF
            CALL var_list_request_element%metadata%set(key, val)
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound)  THEN
      CALL comin_error_set(COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND); RETURN
    ENDIF
  END SUBROUTINE comin_request_set_var_metadata_logical

  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_real(var_descriptor, key, val)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    REAL(wp),                      INTENT(IN)  :: val !< metadata value
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global

    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF

    IF (var_descriptor%id == -1) THEN
      comin_global => comin_descrdata_get_global()
      IF (.NOT. ASSOCIATED(comin_global)) CALL comin_plugin_finish("variable ", "global data missing")

      domain_id_start = 1
      domain_id_end   = comin_global%n_dom
    ELSE
      domain_id_start = var_descriptor%id
      domain_id_end   = var_descriptor%id
    END IF

    ! loop over request list
    lfound = .FALSE.
    DO domain_id = domain_id_start, domain_id_end
      var_descriptor_domain    = var_descriptor
      var_descriptor_domain%id = domain_id

      p => state%comin_var_request_list%first()
      DO WHILE (ASSOCIATED(p))
        ASSOCIATE (var_list_request_element => p%item_value)
          IF (comin_var_descr_match(var_list_request_element%descriptor, var_descriptor_domain)) THEN
            lfound = .TRUE.
            IF ( ALL(var_list_request_element%metadata%query(TRIM(key)) /= &
              &      (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_REAL/) ) ) THEN
              CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
            ENDIF
            CALL var_list_request_element%metadata%set(key, val)
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND); RETURN
    ENDIF
  END SUBROUTINE comin_request_set_var_metadata_real

  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_character(var_descriptor, key, val)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    CHARACTER(LEN=*),              INTENT(IN)  :: val !< metadata value
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global

    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF

    IF (var_descriptor%id == -1) THEN
      comin_global => comin_descrdata_get_global()
      IF (.NOT. ASSOCIATED(comin_global)) CALL comin_plugin_finish("variable ", "global data missing")

      domain_id_start = 1
      domain_id_end   = comin_global%n_dom
    ELSE
      domain_id_start = var_descriptor%id
      domain_id_end   = var_descriptor%id
    END IF

    ! loop over request list
    lfound = .FALSE.
    DO domain_id = domain_id_start, domain_id_end
      var_descriptor_domain    = var_descriptor
      var_descriptor_domain%id = domain_id

      p => state%comin_var_request_list%first()
      DO WHILE (ASSOCIATED(p))
        ASSOCIATE (var_list_request_element => p%item_value)
          IF (comin_var_descr_match(var_list_request_element%descriptor, var_descriptor_domain)) THEN
            lfound = .TRUE.
            IF ( ALL(var_list_request_element%metadata%query(TRIM(key)) /= &
              &      (/COMIN_TYPEID_UNDEFINED, COMIN_TYPEID_CHARACTER/) ) ) THEN
              CALL comin_error_set(COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE); RETURN
            ENDIF
            CALL var_list_request_element%metadata%set(key, val)
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND); RETURN
    ENDIF
  END SUBROUTINE comin_request_set_var_metadata_character

  !> Return a ID (integer) describing the the metadata for a given key
  !> string.
  !! @ingroup common
  INTEGER FUNCTION comin_metadata_get_typeid(var_descriptor, key)  RESULT(typeid)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF
    ! Find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF
    typeid = var_item%metadata%query(key)
  END FUNCTION comin_metadata_get_typeid

  INTEGER(KIND=c_int) FUNCTION comin_metadata_get_typeid_c(var_descriptor, key)  &
    & RESULT(typeid) &
    & BIND(C, name="comin_metadata_get_typeid")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE, INTENT(IN)  :: key !< metadata key (name)
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id

    typeid = INT(comin_metadata_get_typeid(var_descriptor_fortran, &
      &                                    convert_c_string(key)), c_int)
  END FUNCTION comin_metadata_get_typeid_c

  !> Return a metadata container iterator
  !! @ingroup common
  SUBROUTINE comin_metadata_get_iterator(var_descriptor, iterator)
    TYPE(t_comin_var_descriptor), INTENT(IN)         :: var_descriptor
    TYPE(t_comin_var_metadata_iterator), INTENT(OUT) :: iterator
    ! local
    TYPE(t_comin_var_item), POINTER                  :: var_item

    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR); RETURN
    END IF
    ! Find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF
    CALL var_item%metadata%get_iterator(iterator)
  END SUBROUTINE comin_metadata_get_iterator

  !> Provide direct access to C iterator begin
  !! @ingroup common
  FUNCTION comin_metadata_get_iterator_begin_c(var_descriptor) &
    & RESULT(iterator) &
    &  BIND(C, name="comin_metadata_get_iterator_begin")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr) :: iterator
    ! local
    TYPE(t_comin_var_descriptor)    :: var_descriptor_fortran
    TYPE(t_comin_var_item), POINTER :: var_item

    ! Create fortran var descriptor and get var_item
    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    var_item => comin_var_get_from_exposed(var_descriptor_fortran)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF

    CALL comin_keyval_iterator_begin_c(var_item%metadata%comin_metadata_c, iterator)
  END FUNCTION comin_metadata_get_iterator_begin_c

  !> Provide direct access to C iterator end
  !! @ingroup common
  FUNCTION comin_metadata_get_iterator_end_c(var_descriptor) &
    &  RESULT(iterator) &
    &  BIND(C, name="comin_metadata_get_iterator_end")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr) :: iterator
    ! local
    TYPE(t_comin_var_descriptor)    :: var_descriptor_fortran
    TYPE(t_comin_var_item), POINTER :: var_item

    ! Create fortran var descriptor and get var_item
    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    var_item => comin_var_get_from_exposed(var_descriptor_fortran)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      CALL comin_error_set(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED); RETURN
    END IF

    CALL comin_keyval_iterator_end_c(var_item%metadata%comin_metadata_c, iterator)
  END FUNCTION comin_metadata_get_iterator_end_c

  FUNCTION comin_metadata_iterator_get_key_c(it) RESULT(key) &
       & BIND(C, NAME="comin_metadata_iterator_get_key")
    TYPE(c_ptr), INTENT(IN), VALUE :: it
    TYPE(c_ptr) :: key
    key = comin_keyval_iterator_get_key_c(it)
  END FUNCTION comin_metadata_iterator_get_key_c

  FUNCTION comin_metadata_iterator_compare_c(it1, it2) RESULT(equal) &
       & BIND(C, NAME="comin_metadata_iterator_compare")
    TYPE(C_PTR), INTENT(IN), VALUE :: it1
    TYPE(C_PTR), INTENT(IN), VALUE :: it2
    LOGICAL(KIND=C_BOOL) :: equal
    equal = comin_keyval_iterator_compare_c(it1, it2)
  END FUNCTION comin_metadata_iterator_compare_c

  SUBROUTINE comin_metadata_iterator_next_c(it) &
       & BIND(C, NAME="comin_metadata_iterator_next")
    TYPE(C_PTR), INTENT(IN), VALUE :: it
    CALL comin_keyval_iterator_next_c(it)
  END SUBROUTINE comin_metadata_iterator_next_c

  SUBROUTINE comin_metadata_iterator_delete_c(it) &
       & BIND(C, NAME="comin_metadata_iterator_delete")
    TYPE(C_PTR), INTENT(IN), VALUE :: it
    CALL comin_keyval_iterator_delete_c(it)
  END SUBROUTINE comin_metadata_iterator_delete_c

END MODULE comin_metadata
