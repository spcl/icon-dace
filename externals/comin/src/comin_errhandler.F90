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
  USE iso_c_binding,              ONLY: c_ptr
  USE comin_c_utils,              ONLY: convert_c_string
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS, COMIN_INFO, COMIN_WARNING, &
    &                                   COMIN_ERROR_STATUS, COMIN_ERROR_FATAL,    &
    &                                   COMIN_MSG_STR
  USE comin_setup_constants,      ONLY: EP_FINISH
  USE comin_state,                ONLY: state, comin_setup_get_verbosity_level, &
    &                                   comin_current_get_ep

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: comin_plugin_finish
  PUBLIC :: comin_error
  PUBLIC :: comin_message


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

  !> Prints an error message on rank 0
  SUBROUTINE comin_error(error_code, scope)
    INTEGER, INTENT(IN)   :: error_code !< error code
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: scope !< string of the scope, where the error occured

    CHARACTER(LEN=11) :: message_prefix

    IF (error_code < COMIN_SUCCESS .OR. error_code > COMIN_ERROR_FATAL) THEN
      CALL comin_plugin_finish("error", "ERROR: Unknown error code.")
    END IF

    message_prefix = ""
    IF (error_code == COMIN_SUCCESS) THEN
      message_prefix = "SUCCESS"
    ELSE IF (error_code < COMIN_WARNING) THEN
      message_prefix = "INFO"
    ELSE IF (error_code < COMIN_ERROR_STATUS) THEN
      message_prefix = "WARNING"
    ELSE IF (error_code < COMIN_ERROR_FATAL) THEN
      message_prefix = "ERROR"
    ELSE
      message_prefix = "FATAL ERROR"
    END IF

    IF (PRESENT(scope)) THEN
      CALL comin_message("    " // TRIM(scope) // ": " // TRIM(message_prefix) &
        &// ": " // TRIM(COMIN_MSG_STR(error_code)), 0)
    ELSE
      CALL comin_message(TRIM(message_prefix) // ": " &
        &// TRIM(COMIN_MSG_STR(error_code)), 0)
    END IF

  END SUBROUTINE comin_error

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

END MODULE comin_errhandler
