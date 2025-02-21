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

  USE iso_c_binding,         ONLY:  c_int, c_ptr, c_f_pointer, c_double
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
    &                               t_var_request_list
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS,                               &
    &                                   COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED, &
    &                                   COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION,       &
    &                                   COMIN_ERROR_PLUGIN_INIT_PRECISION,           &
    &                                   COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED,   &
    &                                   COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC,       &
    &                                   COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED
  USE comin_errhandler_types, ONLY: comin_host_errhandler_fct
  USE comin_setup_utils,      ONLY: t_comin_setup_version_info, comin_setup_version_compatible

  PRIVATE

  PUBLIC :: comin_setup_init, comin_plugin_init
  PUBLIC :: comin_setup_set_verbosity_level
  PUBLIC :: comin_setup_get_verbosity_level
  PUBLIC :: comin_descrdata_set_fct_glb2loc_cell
  PUBLIC :: comin_setup_errhandler
  PUBLIC :: comin_current_get_ep
  PUBLIC :: comin_current_get_domain_id




  TYPE, PUBLIC :: t_comin_state

    ! verbosity level, related to ICON's `msg_level`.
    ! 0 = silent
    INTEGER :: comin_iverbosity = 0

    LOGICAL :: lstdout = .TRUE.

    ! number of 3rd party plugins associated
    INTEGER :: num_plugins      = 0

    ! currently active 3rd party plugin
    TYPE(t_comin_plugin_info) :: current_plugin

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
    CHARACTER(LEN=32) :: current_datetime

    INTEGER  :: current_domain_id   = DOMAIN_UNDEFINED
    INTEGER  :: current_ep

    !> Flag tracks if primary constructors were completed
    !> prevents further registration of callbacks and registration of variables
    LOGICAL  :: l_primary_done = .FALSE.

  END TYPE t_comin_state

  TYPE(t_comin_state), PUBLIC, POINTER :: state => NULL()

CONTAINS

  !> Initialize the comin state
  !! @ingroup host_interface
  !! This routine needs to be called by the host before any other comin call
  !!
  SUBROUTINE comin_setup_init(lstdout, ierr)
    LOGICAL, INTENT(IN)  :: lstdout     !< do print on stdout or not
    INTEGER, INTENT(OUT) :: ierr        !< error code
    ierr = COMIN_SUCCESS
    IF (ASSOCIATED(state)) THEN
       ! <! cant use comin_message due to circular dependencies
       WRITE(0,*) "ERROR: ComIn is already initalized!"
       ierr = COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED
       RETURN
    END IF
    ALLOCATE(state)
    state%lstdout = lstdout
  END SUBROUTINE comin_setup_init

  !> Initialize the plugin state
  ! This routine is called by the host to set the plugins state
  ! explicitly to the state of the host model. It should by loaded by
  ! the host explicitly from the shared library of the plugin by using
  ! `dlsym`.
  SUBROUTINE comin_plugin_init(state_ptr, host_version, host_wp, ierr) &
       BIND(C)

    ! At the moment assuming only NVHPC and GCC (default way)
    TYPE(C_PTR), VALUE, INTENT(IN)               :: state_ptr
    TYPE(t_comin_setup_version_info), INTENT(IN) :: host_version
    INTEGER(C_INT), INTENT(IN)                   :: host_wp
    INTEGER(C_INT), INTENT(OUT)                  :: ierr

    TYPE(t_comin_state), POINTER :: host_state => NULL()
    ierr = COMIN_SUCCESS

    ! we cant rely on calling methods to obtain the version or wp,
    ! because this might call functions dynamically loaded by the host
    ! or other plugins. Hence we use the version preprocessor and wp
    ! constants driectly.

    IF( .NOT. comin_setup_version_compatible(host_version, &
         & t_comin_setup_version_info(0, 1, 0))) THEN
       WRITE(0,*) "ERROR: host and plugin using incompatible comin versions"
       ierr = COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION
    END IF

    IF( host_wp /= C_DOUBLE ) THEN
       WRITE(0,*) "ERROR: host and plugin using incompatible wp"
       ierr = COMIN_ERROR_PLUGIN_INIT_PRECISION
    END IF

    CALL C_F_POINTER(state_ptr, host_state)

    IF(ASSOCIATED(state) .AND. .NOT. ASSOCIATED(state, host_state)) THEN
       ! <! cant use comin_message due to circular dependencies
       WRITE(0,*) "ERROR: state is already initialized with a different state"
       ierr = COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED
       RETURN
    END IF

    state => host_state
  END SUBROUTINE comin_plugin_init

  !> Sets verbosity level, related to ICON's `msg_level`.
  !>  0 = silent
  !! @ingroup host_interface
  SUBROUTINE comin_setup_set_verbosity_level(iverbosity, ierr)
    INTEGER, INTENT(IN)   :: iverbosity
    INTEGER, INTENT(OUT)  :: ierr

    ierr = COMIN_SUCCESS
    state%comin_iverbosity = iverbosity
  END SUBROUTINE comin_setup_set_verbosity_level


  !> Returns verbosity level
  !! @ingroup host_interface
  FUNCTION comin_setup_get_verbosity_level() BIND(C)
    INTEGER(c_int) :: comin_setup_get_verbosity_level
    comin_setup_get_verbosity_level = state%comin_iverbosity
  END FUNCTION comin_setup_get_verbosity_level


  !> Sets the "global-to-local" index lookup function.
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_set_fct_glb2loc_cell(fct, ierr)
    PROCEDURE(comin_glb2loc_index_lookup_fct) :: fct !< index lookup function
    INTEGER, INTENT(OUT) :: ierr !< error status code
    ierr = COMIN_SUCCESS
    state%comin_descrdata_fct_glb2loc_cell => fct
    IF (.NOT. ASSOCIATED(state%comin_descrdata_fct_glb2loc_cell)) THEN
      ierr = COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC; RETURN
    END IF
  END SUBROUTINE comin_descrdata_set_fct_glb2loc_cell


  !> Sets the global error handler procedure pointer.
  !! @ingroup host_interface
  !!
  !! To be called by the host application (ICON).
  !!
  !! Aborts with an error code if the error handler has already been
  !! set.
  SUBROUTINE comin_setup_errhandler(error_handler, ierr)
    PROCEDURE(comin_host_errhandler_fct) :: error_handler !< error handler
    INTEGER, INTENT(OUT) :: ierr !< return code

    ierr = COMIN_SUCCESS
    state%comin_host_finish => error_handler
    IF (.NOT. ASSOCIATED(state%comin_host_finish)) THEN
      ierr = COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED; RETURN
    END IF
  END SUBROUTINE comin_setup_errhandler


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
