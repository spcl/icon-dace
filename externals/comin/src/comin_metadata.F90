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

  USE iso_c_binding,           ONLY: c_int, c_ptr, C_BOOL, c_double, c_loc
  USE comin_setup_constants,   ONLY: wp
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS,                                       &
    &                                   COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED,                 &
    &                                   COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR,  &
    &                                   COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR, &
    &                                   COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND,                &
    &                                   COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS,      &
    &                                   COMIN_ERROR_METADATA_KEY_NOT_FOUND
  USE comin_errhandler,        ONLY: comin_plugin_finish, comin_message, comin_error
  USE comin_state,             ONLY: state
  USE comin_c_utils,           ONLY: convert_c_string
  USE comin_variable_types,    ONLY: t_comin_var_descriptor, t_comin_var_item,            &
    &                                t_var_request_list_item, t_comin_var_descriptor_c,   &
    &                                t_comin_var_metadata, comin_var_descr_match
  USE comin_variable,          ONLY: comin_var_get_from_exposed
  USE comin_descrdata_types,   ONLY: t_comin_descrdata_global
  USE comin_descrdata,         ONLY: comin_descrdata_get_global
  IMPLICIT NONE

  PUBLIC




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


CONTAINS

  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_integer(descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    INTEGER,                      INTENT(IN)  :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error code
    !
    TYPE(t_comin_var_item), POINTER :: var_item
    CHARACTER(LEN=512)  :: warning_message

    ierr = COMIN_SUCCESS
    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
    ELSE
      CALL comin_metadata_set_integer_core(var_item%metadata, key, val, ierr)
    END IF
    IF (ierr/=COMIN_SUCCESS) THEN
      WRITE(warning_message, '(a,a,a,a)') 'Warning: Metadata ', key, &
        & ' could not be set for field ', descriptor%name
      CALL comin_message(warning_message, 0)
    END IF
  END SUBROUTINE comin_metadata_host_set_integer


  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_logical(descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL,                      INTENT(IN)  :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error code
    !
    TYPE(t_comin_var_item), POINTER :: var_item
    CHARACTER(LEN=512)  :: warning_message

    ierr = COMIN_SUCCESS
    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
    ELSE
      CALL comin_metadata_set_logical_core(var_item%metadata, key, val, ierr)
    END IF
    IF (ierr/=COMIN_SUCCESS) THEN
      WRITE(warning_message, '(a,a,a,a,a)') 'Warning: Metadata ', key, &
        ' could not be set for field ', descriptor%name
      CALL comin_message(warning_message, 0)
    END IF
  END SUBROUTINE comin_metadata_host_set_logical


  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_real(descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    REAL(wp),                     INTENT(IN)  :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error code
    !
    TYPE(t_comin_var_item), POINTER :: var_item
    CHARACTER(LEN=512)  :: warning_message

    ierr = COMIN_SUCCESS
    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
    ELSE
      CALL comin_metadata_set_real_core(var_item%metadata, key, val, ierr)
    END IF
    IF (ierr/=COMIN_SUCCESS) THEN
      WRITE(warning_message, '(a,a,a,a)') 'Warning: Metadata ', key, &
        & ' could not be set for field ', descriptor%name
      CALL comin_message(warning_message, 0)
    END IF
  END SUBROUTINE comin_metadata_host_set_real


  !> Set metadata for item in variable list.
  SUBROUTINE comin_metadata_host_set_character(descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    CHARACTER(LEN=*),             INTENT(IN)  :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error code
    !
    TYPE(t_comin_var_item), POINTER :: var_item
    CHARACTER(LEN=512)  :: warning_message

    ierr = COMIN_SUCCESS
    var_item => comin_var_get_from_exposed(descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
    ELSE
      CALL comin_metadata_set_character_core(var_item%metadata, key, val, ierr)
    END IF
    IF (ierr/=COMIN_SUCCESS) THEN
      WRITE(warning_message, '(a,a,a,a)') 'Warning: Metadata ', key, &
        & ' could not be set for field ', descriptor%name
      CALL comin_message(warning_message, 0)
    END IF
  END SUBROUTINE comin_metadata_host_set_character


  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_integer(var_descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    INTEGER,                      INTENT(OUT) :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ierr = COMIN_SUCCESS
    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR
      RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
      RETURN
    END IF
    CALL comin_metadata_get_integer_core(var_item%metadata, key, val, ierr)
  END SUBROUTINE comin_metadata_get_integer


  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_logical(var_descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL,                      INTENT(OUT) :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ierr = COMIN_SUCCESS
    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR
      RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
      RETURN
    END IF
    CALL comin_metadata_get_logical_core(var_item%metadata, key, val, ierr)
  END SUBROUTINE comin_metadata_get_logical

  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_real(var_descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    REAL(wp),                     INTENT(OUT) :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ierr = COMIN_SUCCESS
    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR
      RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
      RETURN
    END IF
    CALL comin_metadata_get_real_core(var_item%metadata, key, val, ierr)
  END SUBROUTINE comin_metadata_get_real

  !> request the metadata to a variable
  SUBROUTINE comin_metadata_get_character(var_descriptor, key, val, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),             INTENT(IN)  :: key !< metadata key (name)
    CHARACTER(LEN=:), POINTER,    INTENT(OUT) :: val !< metadata value
    INTEGER,                      INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_item), POINTER           :: var_item

    ierr = COMIN_SUCCESS
    ! check if called after primary constructor
    IF (.NOT. state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR
      RETURN
    END IF
    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item)) THEN
      ierr = COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED
      RETURN
    END IF
    CALL comin_metadata_get_character_core(var_item%metadata, key, val, ierr)
  END SUBROUTINE comin_metadata_get_character

  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_integer_c(var_descriptor, key, val, ierr) &
    & BIND(C, NAME="comin_metadata_get_integer")
    TYPE(t_comin_var_descriptor_c), VALUE,  INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    INTEGER(kind=c_int),            INTENT(OUT) :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    INTEGER :: val_fortran, ierr_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    ierr_fortran = COMIN_SUCCESS
    CALL comin_metadata_get_integer(var_descriptor_fortran, &
      &  convert_c_string(key), val_fortran, ierr_fortran)
    val = INT(val_fortran, c_int)
    ierr = INT(ierr_fortran, c_int)
  END SUBROUTINE comin_metadata_get_integer_c


  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_logical_c(var_descriptor, key, val, ierr) &
    & BIND(C, NAME="comin_metadata_get_logical")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL(kind=c_bool),           INTENT(OUT) :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    LOGICAL :: val_fortran
    INTEGER :: ierr_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    ierr_fortran = COMIN_SUCCESS
    CALL comin_metadata_get_logical(var_descriptor_fortran, &
      &  convert_c_string(key), val_fortran, ierr_fortran)
    val = LOGICAL(val_fortran, c_bool)
    ierr = INT(ierr_fortran, c_int)
  END SUBROUTINE comin_metadata_get_logical_c

  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_real_c(var_descriptor, key, val, ierr) &
    & BIND(C, NAME="comin_metadata_get_real")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    REAL(kind=c_double),            INTENT(OUT) :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    REAL(wp) :: val_fortran
    INTEGER  :: ierr_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    ierr_fortran = COMIN_SUCCESS
    CALL comin_metadata_get_real(var_descriptor_fortran, &
      &  convert_c_string(key), val_fortran, ierr_fortran)
    val = REAL(val_fortran, c_double)
    ierr = INT(ierr_fortran, c_int)
  END SUBROUTINE comin_metadata_get_real_c

  !> request the metadata to a variable, C interface
  SUBROUTINE comin_metadata_get_character_c(var_descriptor, key, val, len, ierr) &
    & BIND(C, NAME="comin_metadata_get_character")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    TYPE(c_ptr),                    INTENT(OUT) :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: len !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    ! local
    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    CHARACTER(LEN=:), POINTER :: val_fortran
    INTEGER  :: ierr_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    ierr_fortran = COMIN_SUCCESS
    CALL comin_metadata_get_character(var_descriptor_fortran, &
      &  convert_c_string(key), val_fortran, ierr_fortran)
    val = C_LOC(val_fortran)
    len = LEN_TRIM(val_fortran)
    ierr = INT(ierr_fortran, c_int)
  END SUBROUTINE comin_metadata_get_character_c

  SUBROUTINE comin_metadata_set_integer_c(var_descriptor, key, val, ierr) &
    &  BIND(C, name="comin_metadata_set_integer")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    INTEGER(kind=c_int), VALUE,     INTENT(IN)  :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_integer(var_descriptor_fortran, &
      &                                         convert_c_string(key), val, ierr)
  END SUBROUTINE comin_metadata_set_integer_c


  SUBROUTINE comin_metadata_set_logical_c(var_descriptor, key, val, ierr) &
    &  BIND(C, name="comin_metadata_set_logical")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    LOGICAL(C_BOOL), VALUE,         INTENT(IN)  :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_logical(var_descriptor_fortran, &
      &                                         convert_c_string(key), LOGICAL(val), ierr)
  END SUBROUTINE comin_metadata_set_logical_c

  SUBROUTINE comin_metadata_set_real_c(var_descriptor, key, val, ierr) &
    &  BIND(C, name="comin_metadata_set_real")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    REAL(wp), VALUE,                INTENT(IN)  :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_real(var_descriptor_fortran, &
      &                                      convert_c_string(key), REAL(val, wp), ierr)
  END SUBROUTINE comin_metadata_set_real_c

  SUBROUTINE comin_metadata_set_character_c(var_descriptor, key, val, ierr) &
    &  BIND(C, name="comin_metadata_set_character")
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor !< variable descriptor
    TYPE(c_ptr), VALUE,             INTENT(IN)  :: key !< metadata key (name)
    TYPE(C_PTR), VALUE,             INTENT(IN)  :: val !< metadata value
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_request_set_var_metadata_character(var_descriptor_fortran, &
      &                                           convert_c_string(key), convert_c_string(val), ierr)
  END SUBROUTINE comin_metadata_set_character_c


  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_integer(var_descriptor, key, val, ierr)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    INTEGER,                       INTENT(IN)  :: val !< metadata value
    INTEGER,                       INTENT(OUT) :: ierr !< error status code
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global
    CHARACTER(LEN=512)           :: warning_message

    ierr = COMIN_SUCCESS
    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR
      RETURN
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
            CALL comin_metadata_set_integer_core(var_list_request_element%metadata, key, val, ierr)
            IF (ierr /= COMIN_SUCCESS) THEN
              WRITE(warning_message, '(a,a,a,a)') 'Warning: Metadata ', key, &
                 & ' could not be set for field ', var_descriptor%name
              CALL comin_message(warning_message, 0)
              RETURN
            END IF
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound)  ierr = COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND
  END SUBROUTINE comin_request_set_var_metadata_integer


  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_logical(var_descriptor, key, val, ierr)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    LOGICAL,                       INTENT(IN)  :: val !< metadata value
    INTEGER,                       INTENT(OUT) :: ierr !< error status code
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global
    CHARACTER(LEN=512)           :: warning_message

    ierr = COMIN_SUCCESS
    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR
      RETURN
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
              ierr = COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS
              CALL comin_error(ierr)
              RETURN
            END IF
            CALL comin_metadata_set_logical_core(var_list_request_element%metadata, key, val, ierr)
            IF (ierr /= COMIN_SUCCESS) THEN
              WRITE(warning_message, '(a,a,a,a)') 'Warning: Metadata ', key, &
                  & ' could not be set for field ', var_descriptor%name
              CALL comin_message(warning_message, 0)
              RETURN
            END IF
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound)  ierr = COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND
  END SUBROUTINE comin_request_set_var_metadata_logical

  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_real(var_descriptor, key, val, ierr)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    REAL(wp),                      INTENT(IN)  :: val !< metadata value
    INTEGER,                       INTENT(OUT) :: ierr !< error status code
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global
    CHARACTER(LEN=512)           :: warning_message

    ierr = COMIN_SUCCESS
    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR
      RETURN
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

            IF ((key == "tracer") .AND. (var_descriptor%id /= -1)) THEN
              ierr = COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS
              CALL comin_error(ierr)
              RETURN
            END IF
            CALL comin_metadata_set_real_core(var_list_request_element%metadata, key, val, ierr)
            IF (ierr /= COMIN_SUCCESS) THEN
              WRITE(warning_message, '(a,a,a,a)') 'Warning: Metadata ', key, &
                  & ' could not be set for field ', var_descriptor%name
              CALL comin_message(warning_message, 0)
              RETURN
            END IF
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound)  ierr = COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND
  END SUBROUTINE comin_request_set_var_metadata_real

  !> Sets a specific metadata item (represented by a key-value pair)
  !  for a requested variable. Must be called inside the primary
  !  constructor and the variable must have been previously requested,
  !  otherwise this subroutine aborts with an error status flag.  If
  !  the metadata key does not exist, then this subroutine aborts with
  !  an error status flag.
  !
  SUBROUTINE comin_request_set_var_metadata_character(var_descriptor, key, val, ierr)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor !< variable descriptor
    CHARACTER(LEN=*),              INTENT(IN)  :: key !< metadata key (name)
    CHARACTER(LEN=*),              INTENT(IN)  :: val !< metadata value
    INTEGER,                       INTENT(OUT) :: ierr !< error status code
    ! local
    INTEGER :: domain_id, domain_id_start, domain_id_end
    LOGICAL :: lfound
    TYPE (t_comin_var_descriptor)            :: var_descriptor_domain
    TYPE(t_var_request_list_item),   POINTER :: p
    TYPE(t_comin_descrdata_global),  POINTER :: comin_global
    CHARACTER(LEN=512)           :: warning_message

    ierr = COMIN_SUCCESS
    ! check if called in primary constructor
    IF (state%l_primary_done) THEN
      ierr = COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR
      RETURN
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

            IF ((key == "tracer") .AND. (var_descriptor%id /= -1)) THEN
              ierr = COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS
              CALL comin_error(ierr)
              RETURN
            END IF
            CALL comin_metadata_set_character_core(var_list_request_element%metadata, key, val, ierr)
            IF (ierr /= COMIN_SUCCESS) THEN
              WRITE(warning_message, '(a,a,a,a)') 'Warning: Metadata ', key, &
                 & ' could not be set for field ', var_descriptor%name
              CALL comin_message(warning_message, 0)
              RETURN
            END IF
          END IF
        END ASSOCIATE
        p => p%next()
      END DO
    END DO
    IF (.NOT. lfound)  ierr = COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND
  END SUBROUTINE comin_request_set_var_metadata_character


  !> Return a ID (integer) describing the the metadata for a given key
  !> string.
  !! @ingroup plugin_interface
  INTEGER FUNCTION comin_metadata_get_typeid(key)  RESULT(typeid)
    CHARACTER(LEN=*), INTENT(IN) :: key !< metadata key (name)

    typeid = -1
    IF (key == "restart") THEN
      typeid = 2
    ELSE IF (key == "tracer") THEN
      typeid = 2
    ELSE IF (key == "tracer_turb") THEN
      typeid = 2
    ELSE IF (key == "tracer_conv") THEN
      typeid = 2
    ELSE IF (key == "zaxis_id") THEN
      typeid = 1
    ELSE IF (key == "tracer_vlimit") THEN
      typeid = 1
    ELSE IF (key == "tracer_hlimit") THEN
       typeid = 1
    ELSE IF (key == "tracer_hadv") THEN
       typeid = 1
    ELSE IF (key == "tracer_vadv") THEN
       typeid = 1
    ELSE IF (key == "units") THEN
       typeid = 4
    ELSE IF (key == "standard_name") THEN
       typeid = 4
    ELSE IF (key == "long_name") THEN
       typeid = 4
    ELSE IF (key == "short_name") THEN
       typeid = 4
    END IF
  END FUNCTION comin_metadata_get_typeid


  INTEGER(KIND=c_int) FUNCTION comin_metadata_get_typeid_c(key)  &
    & RESULT(typeid) &
    & BIND(C, name="comin_metadata_get_typeid")
    TYPE(c_ptr), VALUE, INTENT(IN)  :: key !< metadata key (name)
    typeid = INT(comin_metadata_get_typeid(convert_c_string(key)), c_int)
  END FUNCTION comin_metadata_get_typeid_c


  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_set_integer_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),            INTENT(IN)     :: key !< metadata key (name)
    INTEGER,                     INTENT(IN)     :: val !< metadata value
    INTEGER,                     INTENT(OUT)    :: ierr !< error status code

    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "zaxis_id") THEN
      metadata%zaxis_id = val
    ELSE IF (key == "hgrid_id") THEN
      metadata%hgrid_id = val
    ELSE IF (key == "tracer_vlimit") THEN
      metadata%tracer_vlimit = val
    ELSE IF (key == "tracer_hlimit") THEN
      metadata%tracer_hlimit = val
    ELSE IF (key == "tracer_hadv") THEN
      metadata%tracer_hadv = val
    ELSE IF (key == "tracer_vadv") THEN
      metadata%tracer_vadv = val
    ELSE
      ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_set_integer_core


  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_get_integer_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),            INTENT(IN)     :: key !< metadata key (name)
    INTEGER,                     INTENT(OUT)    :: val !< metadata value
    INTEGER,                     INTENT(OUT)    :: ierr !< error status code

    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "zaxis_id") THEN
      val = metadata%zaxis_id
    ELSE IF (key == "hgrid_id") THEN
      val = metadata%hgrid_id
    ELSE IF (key == "tracer_vlimit") THEN
      val = metadata%tracer_vlimit
    ELSE IF (key == "tracer_hlimit") THEN
      val = metadata%tracer_hlimit
    ELSE IF (key == "tracer_hadv") THEN
      val = metadata%tracer_hadv
    ELSE IF (key == "tracer_vadv") THEN
      val = metadata%tracer_vadv
    ELSE
      ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_get_integer_core

  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_set_real_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),            INTENT(IN)     :: key !< metadata key (name)
    REAL(wp),                    INTENT(IN)     :: val !< metadata value
    INTEGER,                     INTENT(OUT)    :: ierr !< error status code

    ! this is currently stupid dummy code to avoid warnings
    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "" .AND. val == 0 .AND. metadata%tracer) THEN
       ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    ELSE
      ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_set_real_core

  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_get_real_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),            INTENT(IN)     :: key !< metadata key (name)
    REAL(wp),                    INTENT(OUT)    :: val !< metadata value
    INTEGER,                     INTENT(OUT)    :: ierr !< error status code

    ! this is currently stupid dummy code to avoid warnings
    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "zaxis_id" .AND. metadata%tracer) THEN
       ierr = -1
       val = 0
    ELSE
      ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_get_real_core

  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_set_logical_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),            INTENT(IN)     :: key !< metadata key (name)
    LOGICAL,                     INTENT(IN)     :: val !< metadata value
    INTEGER,                     INTENT(OUT)    :: ierr !< error status code

    ierr = COMIN_SUCCESS
    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "restart") THEN
      metadata%restart = val
    ELSE IF (key == "tracer") THEN
      metadata%tracer = val
    ELSE IF (key == "tracer_turb") THEN
      metadata%tracer_turb = val
    ELSE IF (key == "tracer_conv") THEN
      metadata%tracer_conv = val
    ELSE
      ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_set_logical_core

  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_get_logical_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),            INTENT(IN)     :: key !< metadata key (name)
    LOGICAL,                     INTENT(OUT)    :: val !< metadata value
    INTEGER,                     INTENT(OUT)    :: ierr !< error status code

    ierr = COMIN_SUCCESS
    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "restart") THEN
      val = metadata%restart
    ELSE IF (key == "tracer") THEN
      val = metadata%tracer
    ELSE IF (key == "tracer_turb") THEN
      val = metadata%tracer_turb
    ELSE IF (key == "tracer_conv") THEN
      val = metadata%tracer_conv
    ELSE
      ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_get_logical_core

  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_set_character_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),            INTENT(IN)     :: key !< metadata key (name)
    CHARACTER(LEN=*),            INTENT(IN)     :: val !< metadata value
    INTEGER,                     INTENT(OUT)    :: ierr !< error status code

    ierr = COMIN_SUCCESS
    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "units") THEN
       metadata%units = val
    ELSE IF (key == "standard_name") THEN
       metadata%standard_name = val
    ELSE IF (key == "long_name") THEN
       metadata%long_name = val
    ELSE IF (key == "short_name") THEN
       metadata%short_name = val
    ELSE
       ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_set_character_core

  ! Auxiliary core function: set data member in internal data structure for
  ! field metadata.
  SUBROUTINE comin_metadata_get_character_core(metadata, key, val, ierr)
    TYPE (t_comin_var_metadata), TARGET, INTENT(INOUT)  :: metadata !< data structure containing metadata
    CHARACTER(LEN=*),                    INTENT(IN)     :: key !< metadata key (name)
    CHARACTER(LEN=:), POINTER,           INTENT(OUT)    :: val !< metadata value
    INTEGER,                             INTENT(OUT)    :: ierr !< error status code

    ierr = COMIN_SUCCESS
    ! later on, the following code passage should be replaced
    ! by a dictionary data structure.
    IF (key == "units") THEN
      val => metadata%units
    ELSE IF (key == "standard_name") THEN
       val => metadata%standard_name
    ELSE IF (key == "long_name") THEN
       val => metadata%long_name
    ELSE IF (key == "short_name") THEN
       val => metadata%short_name
    ELSE
      ierr = COMIN_ERROR_METADATA_KEY_NOT_FOUND
    END IF
  END SUBROUTINE comin_metadata_get_character_core


END MODULE comin_metadata
