!> @file comin_state.F90
!! @brief Data shared between host and plugins.
!
!  @authors 10/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_state

  USE iso_c_binding,         ONLY:  c_int, c_ptr, c_double
  USE comin_setup_constants, ONLY:  DOMAIN_UNDEFINED, wp
  USE comin_plugin_types,    ONLY:  t_comin_plugin_info
  USE comin_parallel_types,  ONLY:  t_comin_parallel_info
  USE comin_callback_types,  ONLY:  t_callback_list,                &
    &                               t_comin_callback_context
  USE comin_descrdata_types, ONLY:  t_comin_descrdata_global,       &
    &                               t_comin_descrdata_domain,       &
    &                               t_comin_descrdata_simulation_interval, &
    &                               comin_glb2loc_index_lookup_fct
  USE comin_variable_types,  ONLY:  t_var_descr_list, t_var_list,   &
    &                               t_comin_var_list_context,       &
    &                               t_var_request_list,             &
    &                               comin_var_sync_device_mem_fct
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS
  USE comin_errhandler_types, ONLY: comin_host_errhandler_fct

  PRIVATE
  PUBLIC :: comin_setup_set_verbosity_level
  PUBLIC :: comin_setup_get_verbosity_level
  PUBLIC :: comin_current_get_ep
  PUBLIC :: comin_current_get_domain_id

#include "comin_global.inc"

  TYPE, PUBLIC :: t_comin_state

    ! current error state of comin
    INTEGER :: errcode = COMIN_SUCCESS

    ! verbosity level, related to ICON's `msg_level`.
    ! 0 = silent
    INTEGER :: comin_iverbosity = 0

    LOGICAL :: lstdout = .TRUE.

    ! number of 3rd party plugins associated
    INTEGER :: num_plugins      = 0

    ! plugin infos
    TYPE(t_comin_plugin_info), ALLOCATABLE :: plugin_info(:)

    ! currently active 3rd party plugin
    TYPE(t_comin_plugin_info), POINTER :: current_plugin

    ! (Non-public) data structure
    TYPE(t_comin_parallel_info) :: parallel_info

    !> information about each entry point/callback
    TYPE(t_callback_list)                       :: comin_callback_list
    TYPE(t_comin_callback_context), ALLOCATABLE :: comin_callback_context(:,:)
    INTEGER, ALLOCATABLE                        :: comin_callback_order(:,:)

    !> Create descriptive data structures
    TYPE(t_comin_descrdata_global)                 :: comin_descrdata_global
    TYPE(t_comin_descrdata_domain), ALLOCATABLE    :: comin_descrdata_domain(:)
    TYPE(t_comin_descrdata_simulation_interval)    :: comin_descrdata_simulation_interval
    REAL(wp), ALLOCATABLE                          :: comin_descrdata_timesteplength(:)

    !> List of all variables available (exported) from ICON
    TYPE(t_var_descr_list) :: comin_var_descr_list
    TYPE(t_var_list)       :: comin_var_list

    !> List of variables in context
    TYPE(t_comin_var_list_context), ALLOCATABLE :: comin_var_list_context(:,:)

    !> List of all variables available (exported) from ICON
    TYPE(t_var_request_list) :: comin_var_request_list

    ! translates global cell indices to process-local indices
    PROCEDURE (comin_glb2loc_index_lookup_fct), POINTER, NOPASS :: comin_descrdata_fct_glb2loc_cell

    ! Global variable which contains the callback to ICON's "finish"
    ! routine.
    PROCEDURE(comin_host_errhandler_fct), POINTER, NOPASS :: comin_host_finish => NULL()

    ! current simulation date-time stamp (ISO 8601)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: current_datetime

    INTEGER  :: current_domain_id   = DOMAIN_UNDEFINED
    INTEGER  :: current_ep

    !> Flag tracks if primary constructors were completed
    !> prevents further registration of callbacks and registration of variables
    LOGICAL  :: l_primary_done = .FALSE.

    PROCEDURE (comin_var_sync_device_mem_fct), POINTER, NOPASS :: sync_device_mem

  END TYPE t_comin_state

  TYPE(t_comin_state), PUBLIC, POINTER :: state => NULL()

CONTAINS

  !> Sets verbosity level, related to ICON's `msg_level`.
  !>  0 = silent
  !! @ingroup host_interface
  SUBROUTINE comin_setup_set_verbosity_level(iverbosity)
    INTEGER, INTENT(IN)   :: iverbosity

    state%comin_iverbosity = iverbosity
  END SUBROUTINE comin_setup_set_verbosity_level

  !> Returns verbosity level
  !! @ingroup host_interface
  FUNCTION comin_setup_get_verbosity_level() BIND(C)
    INTEGER(c_int) :: comin_setup_get_verbosity_level
    comin_setup_get_verbosity_level = state%comin_iverbosity
  END FUNCTION comin_setup_get_verbosity_level

  !> Access information on the current entry point being processed by ComIn.
  !! @ingroup plugin_interface
  FUNCTION comin_current_get_ep()  BIND(C)
    INTEGER(c_int) :: comin_current_get_ep
    comin_current_get_ep = INT(state%current_ep, C_INT)
  END FUNCTION comin_current_get_ep

  !> Request information on current domain
  !! @ingroup plugin_interface
  FUNCTION comin_current_get_domain_id()  &
    &  BIND(C)
    INTEGER(c_int) :: comin_current_get_domain_id

    comin_current_get_domain_id = -1
    IF (state%current_domain_id /= DOMAIN_UNDEFINED) THEN
      comin_current_get_domain_id = INT(state%current_domain_id, c_int)
    END IF
  END FUNCTION comin_current_get_domain_id

END MODULE comin_state
