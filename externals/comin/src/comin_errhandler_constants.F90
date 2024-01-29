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

  PRIVATE
  PUBLIC :: COMIN_SUCCESS,                                            &
    &       COMIN_INFO,                                               &
    &       COMIN_WARNING,                                            &
    &       COMIN_ERROR_STATUS,                                       &
    &       COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR, &
    &       COMIN_ERROR_CALLBACK_COMPLETE,                            &
    &       COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN,                       &
    &       COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC,                    &
    &       COMIN_ERROR_DESCRDATA_FINALIZE,                           &
    &       COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR,      &
    &       COMIN_ERROR_METADATA_KEY_NOT_FOUND,                       &
    &       COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR,       &
    &       COMIN_ERROR_SETUP_FINALIZE,                               &
    &       COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED,              &
    &       COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION,                    &
    &       COMIN_ERROR_PLUGIN_INIT_PRECISION,                        &
    &       COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED,                &
    &       COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED,              &
    &       COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET,                     &
    &       COMIN_ERROR_SETUP_PRECISION_TEST_FAILED,                  &
    &       COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR,         &
    &       COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE,          &
    &       COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE,     &
    &       COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND,                     &
    &       COMIN_ERROR_VARIABLE_NOT_TRACER_IN_TURBULENT_TRANSPORT,   &
    &       COMIN_ERROR_VARIABLE_NOT_TRACER_IN_CONVECTIVE_TRANSPORT,  &
    &       COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED,                      &
    &       COMIN_ERROR_FIELD_NOT_ALLOCATED,                          &
    &       COMIN_ERROR_POINTER_NOT_ASSOCIATED,                       &
    &       COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS,           &
    &       COMIN_ERROR_FATAL,                                        &
    &       COMIN_MSG_STR

#include "comin_global.inc"

  !> define list of error points
  !> COMIN_ERROR_FATAL should always be the last entry
  ENUM, BIND(C)
    ENUMERATOR :: COMIN_SUCCESS = 0,                                        &
      ! --- INFO ---
      &           COMIN_INFO,                                               &
      ! --- WARNING ---
      &           COMIN_WARNING,                                            &
      ! --- ERROR ---
      &           COMIN_ERROR_STATUS,                                       &
      &           COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR, &
      &           COMIN_ERROR_CALLBACK_COMPLETE,                            &
      &           COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN,                       &
      &           COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC,                    &
      &           COMIN_ERROR_DESCRDATA_FINALIZE,                           &
      &           COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR,      &
      &           COMIN_ERROR_METADATA_KEY_NOT_FOUND,                       &
      &           COMIN_ERROR_METADATA_GET_INSIdE_PRIMARYCONSTRUCTOR,       &
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
      &           COMIN_ERROR_VARIABLE_NOT_TRACER_IN_TURBULENT_TRANSPORT,   &
      &           COMIN_ERROR_VARIABLE_NOT_TRACER_IN_CONVECTIVE_TRANSPORT,  &
      &           COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED,                      &
      &           COMIN_ERROR_FIELD_NOT_ALLOCATED,                          &
      &           COMIN_ERROR_POINTER_NOT_ASSOCIATED,                       &
      &           COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS,           &
      ! --- FINAL ERROR
      &           COMIN_ERROR_FATAL
  END ENUM

  CHARACTER(LEN=MAX_LEN_ERR_MESSAGE), PARAMETER :: COMIN_MSG_STR(0:COMIN_ERROR_FATAL) = [ &
    &  "Success                                                                                         ", &
    ! --- INFOS ---
    &  "Info                                                                                            ", &
    ! --- WARNINGS ---
    &  "Warning                                                                                         ", &
    ! --- ERRORS ---
    &  "Error                                                                                           ", &
    &  "Callbacks cannot be registered after primary constructor.                                       ", &
    &  "Callback registration could not be completed.                                                   ", &
    &  "Unknown id for entry point.                                                                     ", &
    &  "Unable to set descriptive data function: glb2loc. Function not associated.                      ", &
    &  "Unable to clean up descriptive data structure (finalize).                                       ", &
    &  "Cannot set metadata outside primary constructor.                                                ", &
    &  "Metadata key not found.                                                                         ", &
    &  "Unable to get metadata inside primary constructor.                                              ", &
    &  "Setup finalized failed.                                                                         ", &
    &  "ComIn is already initilized.                                                                    ", &
    &  "Host and plugin are using incompatible ComIn versions.                                          ", &
    &  "Host and plugin are using incompatible precision (wp).                                          ", &
    &  "State is already initialized with a different state.                                            ", &
    &  "The host model error handler procedure cannot be found.                                         ", &
    &  "Error handler not set (setup check).                                                            ", &
    &  "Precision test failed (setup check).                                                            ", &
    &  "Variables cannot be requested after primary constructor.                                        ", &
    &  "Requested variable already exists and is exclusive.                                             ", &
    &  "Requested exclusive variable already exists.                                                    ", &
    &  "var_descriptor not found.                                                                       ", &
    &  "Non-tracer variable cannot be added to turbulent transport.                                     ", &
    &  "Non-tracer variable cannot be added to convective transport.                                    ", &
    &  "var_item not associated.                                                                        ", &
    &  "Filed not allocated.                                                                            ", &
    &  "Pointer not associated.                                                                         ", &
    &  "Traers need to be requested for all domains (id=-1).                                            ", &
    ! --- FINAL ERROR
    &  "Fatal error                                                                                     " ]

END MODULE comin_errhandler_constants
