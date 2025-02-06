!> @file comin_errhandler.F90
!! @brief Utility functions for error handling.
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_errhandler

  USE mpi
  USE iso_c_binding,              ONLY: c_ptr, c_int, c_char, c_bool
  USE comin_c_utils,              ONLY: convert_c_string, convert_f_string
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS, COMIN_INFO, COMIN_WARNING, &
    &                                   COMIN_ERROR_STATUS, COMIN_ERROR_FATAL, comin_errhandler_get_string
  USE comin_setup_constants,      ONLY: EP_FINISH
  USE comin_state,                ONLY: state, comin_setup_get_verbosity_level, &
    &                                   comin_current_get_ep

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: comin_plugin_finish
  PUBLIC :: comin_message
  PUBLIC :: comin_error_get_message
  PUBLIC :: comin_error_check
  PUBLIC :: comin_error_set, comin_error_get
  PUBLIC :: comin_error_set_errors_return

#include "comin_global.inc"

CONTAINS

  !> Wrapper function for callback to ICON's "finish" routine.
  !! @ingroup plugin_interface
  SUBROUTINE comin_plugin_finish(routine, text)
    CHARACTER(LEN=*), INTENT(IN) :: routine
    CHARACTER(LEN=*), INTENT(IN) :: text
    !
    INTEGER, PARAMETER :: exit_no = 1
    INTEGER :: ierr

    ! skip this routine if the "finish" call was triggered by a plugin
    ! inside the entry point EP_FINISH itself:
    IF (comin_current_get_ep() == EP_FINISH)  RETURN

    IF (ASSOCIATED(state%comin_host_finish)) THEN
      CALL state%comin_host_finish(routine, text)
    ELSE
      WRITE (0,*) routine, "  ", text
      CALL MPI_ABORT(MPI_COMM_WORLD, exit_no, ierr)
      STOP exit_no
    END IF
  END SUBROUTINE comin_plugin_finish

  !> C-wrapper for the `comin_plugin_finish` subroutine.
  !
  SUBROUTINE comin_plugin_finish_c(routine, text) &
    &  BIND(C, name="comin_plugin_finish")
    TYPE(c_ptr), VALUE, INTENT(IN) :: routine
    TYPE(c_ptr), VALUE, INTENT(IN) :: text

    CALL comin_plugin_finish(convert_c_string(routine), convert_c_string(text))
  END SUBROUTINE comin_plugin_finish_c

  !> Prints a message on rank 0 if the global verbosity level larger than lvl
  SUBROUTINE comin_message(message, lvl)
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER, INTENT(IN)   :: lvl

    INTEGER :: iverbosity

    iverbosity = comin_setup_get_verbosity_level()

    IF (lvl < 0) THEN
      CALL comin_plugin_finish("message", "ERROR: Message level must be non-negative.")
    END IF

    IF (state%lstdout .AND. (iverbosity > lvl)) THEN
      WRITE(0, *) TRIM(message)
    END IF
  END SUBROUTINE comin_message

  SUBROUTINE comin_error_get_message(error_code, category, message)
    INTEGER, INTENT(IN)                               :: error_code
    CHARACTER(LEN=11), INTENT(INOUT)                  :: category
    CHARACTER(LEN=MAX_LEN_ERR_MESSAGE), INTENT(INOUT) :: message

    IF (error_code < COMIN_SUCCESS .OR. error_code > COMIN_ERROR_FATAL) THEN
      CALL comin_plugin_finish("error", "ERROR: Unknown error code.")
    END IF

    category = ""
    IF (error_code == COMIN_SUCCESS) THEN
      category = "SUCCESS"
    ELSE IF (error_code < COMIN_WARNING) THEN
      category = "INFO"
    ELSE IF (error_code < COMIN_ERROR_STATUS) THEN
      category = "WARNING"
    ELSE IF (error_code < COMIN_ERROR_FATAL) THEN
      category = "ERROR"
    ELSE
      category = "FATAL ERROR"
    END IF

    message = TRIM(comin_errhandler_get_string(error_code))
  END SUBROUTINE comin_error_get_message

  SUBROUTINE comin_error_get_message_c(error_code, category, message) &
    & BIND(C, name="comin_error_get_message")
    INTEGER(C_INT), VALUE,  INTENT(IN)  :: error_code
    CHARACTER(KIND=C_CHAR), INTENT(OUT) :: category(11)
    CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(MAX_LEN_ERR_MESSAGE)

    CHARACTER(LEN=11)                  :: category_f
    CHARACTER(LEN=MAX_LEN_ERR_MESSAGE) :: message_f
    CALL comin_error_get_message(error_code, category_f, message_f)
    CALL convert_f_string(category_f, category)
    CALL convert_f_string(message_f, message)
  END SUBROUTINE comin_error_get_message_c

  ! check the error code:
  ! does nothing if error_code == COMIN_SUCCESS
  ! prints the corresponding message
  ! calls comin_plugin_finish if the error_code is an error or fatal error
  SUBROUTINE comin_error_check() BIND(C)

    INTEGER :: error_code
    CHARACTER(LEN=11) :: message_prefix
    CHARACTER(LEN=MAX_LEN_ERR_MESSAGE) :: message

    error_code = state%errcode
    IF(error_code == COMIN_SUCCESS) RETURN

    CALL comin_error_get_message(error_code, message_prefix, message)

    IF (error_code < COMIN_ERROR_STATUS) THEN
      CALL comin_message(TRIM(message_prefix) // ": " &
           &// TRIM(message), 0)
    ELSE
      CALL comin_plugin_finish(state%current_plugin%name, &
                               TRIM(message_prefix) // ": " // TRIM(message))
    END IF
  END SUBROUTINE comin_error_check

  SUBROUTINE comin_error_set(errcode)
    INTEGER, INTENT(IN) :: errcode
    state%errcode = errcode
    IF (.NOT. ASSOCIATED(state%current_plugin)) THEN
      CALL comin_error_check()
    ELSE
      IF(.NOT. state%current_plugin%errors_return) THEN
        CALL comin_error_check()
      END IF
    END IF
  END SUBROUTINE comin_error_set

  FUNCTION comin_error_get() &
       & BIND(C)
    INTEGER(C_INT) :: comin_error_get
    comin_error_get = state%errcode
  END FUNCTION comin_error_get

  SUBROUTINE comin_error_set_errors_return(errors_return) BIND(C)
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: errors_return
    state%current_plugin%errors_return = errors_return
  END SUBROUTINE comin_error_set_errors_return

END MODULE comin_errhandler
