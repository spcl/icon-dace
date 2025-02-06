!> @file comin_errhandler_constants.F90
!! @brief Constants for error handling.
!
!  @authors 09/2023 :: ICON Community Interface  <icon@dwd.de>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_errhandler_constants

  IMPLICIT NONE

  PUBLIC
#include "comin_global.inc"

  !> define list of error points
  !> COMIN_ERROR_FATAL should always be the last entry
  ENUM, BIND(C)
    ENUMERATOR :: COMIN_SUCCESS = 0,                                        &
      &           COMIN_INFO,                                               &
      &           COMIN_WARNING,                                            &
      &           COMIN_ERROR_STATUS,                                       &
      &           COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR, &
      &           COMIN_ERROR_CALLBACK_COMPLETE,                            &
      &           COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN,                       &
      &           COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC,                    &
      &           COMIN_ERROR_DESCRDATA_FINALIZE,                           &
      &           COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR,      &
      &           COMIN_ERROR_METADATA_KEY_NOT_FOUND,                       &
      &           COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR,       &
      &           COMIN_ERROR_SETUP_FINALIZE,                               &
      &           COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED,              &
      &           COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION,                    &
      &           COMIN_ERROR_PLUGIN_INIT_PRECISION,                        &
      &           COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED,                &
      &           COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED,              &
      &           COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET,                     &
      &           COMIN_ERROR_SETUP_PRECISION_TEST_FAILED,                  &
      &           COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR,         &
      &           COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE,          &
      &           COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE,     &
      &           COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND,                     &
      &           COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED,                      &
      &           COMIN_ERROR_FIELD_NOT_ALLOCATED,                          &
      &           COMIN_ERROR_POINTER_NOT_ASSOCIATED,                       &
      &           COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS,           &
      &           COMIN_ERROR_VAR_SYNC_DEVICE_MEM_NOT_ASSOCIATED,           &
      &           COMIN_ERROR_VAR_GET_OUTSIDE_SECONDARY_CONSTRUCTOR,        &
      &           COMIN_ERROR_VAR_GET_NO_DEVICE,                            &
      &           COMIN_ERROR_VAR_GET_VARIABLE_NOT_FOUND,                   &
      &           COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE,               &
      &           COMIN_ERROR_FATAL
  END ENUM

CONTAINS

  FUNCTION comin_errhandler_get_string(err_code) RESULT(string)
    CHARACTER(LEN=MAX_LEN_ERR_MESSAGE) :: string
    INTEGER, INTENT(IN) :: err_code

    SELECT CASE(err_code)
    CASE(COMIN_SUCCESS)
      string = "Success"
    CASE(COMIN_INFO)
      string = "Info"
    CASE(COMIN_WARNING)
      string = "Warning"
    CASE(COMIN_ERROR_STATUS)
      string = "Error"
    CASE(COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR)
      string = "Callbacks cannot be registered after primary constructor."
    CASE(COMIN_ERROR_CALLBACK_COMPLETE)
      string = "Callback registration could not be completed."
    CASE(COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN)
      string = "Unknown id for ENTRY point."
    CASE(COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC)
      string = "Unable to set descriptive DATA FUNCTION: glb2loc. FUNCTION not associated."
    CASE(COMIN_ERROR_DESCRDATA_FINALIZE)
      string = "Unable to clean up descriptive DATA structure (finalize)."
    CASE(COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR)
      string = "Cannot set metadata outside primary constructor."
    CASE(COMIN_ERROR_METADATA_KEY_NOT_FOUND)
      string = "Metadata key not found."
    CASE(COMIN_ERROR_METADATA_GET_INSIdE_PRIMARYCONSTRUCTOR)
      string = "Unable to get metadata inside primary constructor."
    CASE(COMIN_ERROR_SETUP_FINALIZE)
      string = "Setup finalized failed."
    CASE(COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED)
      string = "ComIn is already initilized."
    CASE(COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION)
      string = "Host and plugin are using incompatible ComIn versions."
    CASE(COMIN_ERROR_PLUGIN_INIT_PRECISION)
      string = "Host and plugin are using incompatible PRECISION (wp)."
    CASE(COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED)
      string = "State is already initialized WITH a different state."
    CASE(COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED)
      string = "The host model error handler PROCEDURE cannot be found."
    CASE(COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET)
      string = "Error handler not set (setup check)."
    CASE(COMIN_ERROR_SETUP_PRECISION_TEST_FAILED)
      string = "PRECISION test failed (setup check)."
    CASE(COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR)
      string = "Variables cannot be requested after primary constructor."
    CASE(COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE)
      string = "Requested variable already exists and is exclusive."
    CASE(COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE)
      string = "Requested exclusive variable already exists."
    CASE(COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND)
      string = "var_descriptor not found."
    CASE(COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED)
      string = "var_item not associated."
    CASE(COMIN_ERROR_FIELD_NOT_ALLOCATED)
      string = "Filed not allocated."
    CASE(COMIN_ERROR_POINTER_NOT_ASSOCIATED)
      string = "POINTER not associated."
    CASE(COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS)
      string = "Traers need to be requested for all domains (id=-1)."
    CASE(COMIN_ERROR_VAR_SYNC_DEVICE_MEM_NOT_ASSOCIATED)
      string = "Sync device memory callback not set (setup check)."
    CASE(COMIN_ERROR_VAR_GET_OUTSIDE_SECONDARY_CONSTRUCTOR)
      string = "Cannot get variable outside of secondary constructor"
    CASE(COMIN_ERROR_VAR_GET_NO_DEVICE)
      string = "No device available"
    CASE(COMIN_ERROR_VAR_GET_VARIABLE_NOT_FOUND)
      string = "Cannot find variable"
    CASE(COMIN_ERROR_FATAL)
      string = "Fatal error"
    END SELECT
  END FUNCTION comin_errhandler_get_string
END MODULE comin_errhandler_constants
