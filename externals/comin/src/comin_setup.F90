!> @file comin_setup.F90
!! @brief Routines to set up ComIn (except for callbacks).
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_setup

  USE comin_callback,        ONLY: comin_callback_complete
  USE comin_variable,        ONLY: comin_var_complete
  USE comin_setup_constants, ONLY: wp
  USE comin_setup_utils,     ONLY: t_comin_setup_version_info, comin_setup_get_version, &
       &                           comin_setup_version_compatible
  USE comin_state,           ONLY: t_comin_state, state
  USE comin_c_utils,         ONLY: convert_c_string, convert_f_string
  USE comin_parallel,        ONLY: comin_parallel_free_mpi_comms
  USE comin_plugin_types,    ONLY: t_comin_plugin_description,        &
    &                              t_comin_plugin_info
  USE comin_errhandler_constants, ONLY: COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET,               &
    &                                   COMIN_ERROR_SETUP_PRECISION_TEST_FAILED,            &
    &                                   COMIN_ERROR_SETUP_FINALIZE,                         &
    &                                   COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED,        &
    &                                   COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC,              &
    &                                   COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION,              &
    &                                   COMIN_ERROR_PLUGIN_INIT_PRECISION,                  &
    &                                   COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED,          &
    &                                   COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED
  USE comin_errhandler,      ONLY: comin_message, comin_plugin_finish, comin_error_set
  USE comin_errhandler_types, ONLY: comin_host_errhandler_fct
  USE comin_descrdata_types, ONLY:  comin_glb2loc_index_lookup_fct

  USE iso_c_binding,         ONLY: c_int, c_ptr, c_funptr, c_f_procpointer, c_null_char, &
    &                              C_ASSOCIATED, c_null_ptr, c_loc, c_char, c_double, c_f_pointer
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_setup_check
  PUBLIC :: comin_plugin_primaryconstructor
  PUBLIC :: comin_setup_init
  PUBLIC :: comin_setup_finalize
  PUBLIC :: comin_current_get_plugin_info
  PUBLIC :: comin_plugin_init
  PUBLIC :: comin_descrdata_set_fct_glb2loc_cell
  PUBLIC :: comin_setup_errhandler

  !> interface for a primary constructor call
  ABSTRACT INTERFACE
    SUBROUTINE comin_plugin_init_fct(state_ptr, host_version, host_wp) BIND(C)
      IMPORT c_ptr, c_int, t_comin_setup_version_info

#ifndef __NVCOMPILER
      TYPE(C_PTR), VALUE, INTENT(IN)               :: state_ptr
#else
      TYPE(C_PTR), INTENT(IN)                      :: state_ptr
#endif

      TYPE(t_comin_setup_version_info), INTENT(IN) :: host_version
      INTEGER(C_INT), INTENT(IN)                   :: host_wp
    END SUBROUTINE comin_plugin_init_fct

    SUBROUTINE comin_primaryconstructor_fct() &
      &  BIND(C)
    END SUBROUTINE comin_primaryconstructor_fct

  END INTERFACE

  INTEGER(c_int), PARAMETER :: rtld_now    =   2 ! (value extracted from the C header file)
  !
  ! interface to linux API
  INTERFACE
    FUNCTION dlopen_ptr(filename,mode) BIND(c,name="dlopen")
      ! void *dlopen(const char *filename, int mode);
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dlopen_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: filename
      INTEGER(c_int), VALUE :: mode
    END FUNCTION dlopen_ptr

    FUNCTION dlsym(handle,name) BIND(c,name="dlsym")
      ! void *dlsym(void *handle, const char *name);
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_funptr) :: dlsym
      TYPE(c_ptr), VALUE :: handle
      CHARACTER(c_char), INTENT(in) :: name(*)
    END FUNCTION dlsym

    FUNCTION dlclose(handle) BIND(c,name="dlclose")
      ! int dlclose(void *handle);
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER(c_int) :: dlclose
      TYPE(c_ptr), VALUE :: handle
    END FUNCTION dlclose

    FUNCTION DLError() RESULT(error) BIND(C,NAME="dlerror")
      ! char *dlerror(void);
      USE ISO_C_BINDING
      TYPE(C_PTR) :: error
    END FUNCTION DLError

  END INTERFACE

  ! list type
  TYPE :: comin_primaryconstructor_ptr
    PROCEDURE(comin_primaryconstructor_fct), POINTER, NOPASS :: fct_ptr
  END TYPE comin_primaryconstructor_ptr

  !> list of primary constructors
  TYPE(c_ptr), ALLOCATABLE :: dl_handles(:)

#include "comin_version.inc"

CONTAINS

  FUNCTION dlopen(filename, mode)
    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER(c_int), VALUE :: mode
    !
    TYPE(c_ptr) :: dlopen
    CHARACTER(len=1, kind=c_char), TARGET :: c_filename(LEN(filename)+1)
    CALL convert_f_string(TRIM(filename), c_filename)
    dlopen = dlopen_ptr(C_LOC(c_filename), mode)
  END FUNCTION dlopen

  !> Execute primary constructors.
 !! @ingroup host_interface
  SUBROUTINE comin_plugin_primaryconstructor(plugin_list)
    !> list of dynamic libs:
    TYPE(t_comin_plugin_description), INTENT(IN), TARGET :: plugin_list(:)
    !
    INTEGER :: i, last_sep_idx
    TYPE(c_funptr) :: setup_fct_c, plugin_init_fct_c
    PROCEDURE(comin_primaryconstructor_fct), BIND(C), POINTER :: setup_fct
    PROCEDURE(comin_plugin_init_fct), BIND(C), POINTER :: plugin_init_fct

    state%num_plugins = SIZE(plugin_list)
    ALLOCATE(dl_handles(state%num_plugins))
    ALLOCATE(state%plugin_info(state%num_plugins))
    DO i=1,state%num_plugins

      IF (TRIM(plugin_list(i)%plugin_library) .EQ. "") THEN
        dl_handles(i) = dlopen_ptr(c_null_ptr, RTLD_NOW)
      ELSE
        dl_handles(i) = dlopen(plugin_list(i)%plugin_library, RTLD_NOW)
      END IF
      IF (.NOT. C_ASSOCIATED(dl_handles(i))) &
        & CALL comin_plugin_finish("comin_plugin_primaryconstructor", &
        &                          "ERROR: Cannot load plugin " // convert_c_string(dlerror()))

      ! We load the symbol comin_plugin_init explicitly from the library
      plugin_init_fct_c = dlsym(dl_handles(i), "comin_plugin_init"//c_null_char)
      IF (.NOT. C_ASSOCIATED(plugin_init_fct_c))  &
        & CALL comin_plugin_finish("comin_plugin_primaryconstructor", &
        &                          "Cannot load 'comin_plugin_init' from plugin: " // convert_c_string(dlerror()))
      CALL C_F_PROCPOINTER( plugin_init_fct_c, plugin_init_fct )
      CALL plugin_init_fct(C_LOC(state), comin_setup_get_version(), wp)

      setup_fct_c = dlsym(dl_handles(i), TRIM(plugin_list(i)%primary_constructor)//c_null_char)
      IF (.NOT. C_ASSOCIATED(setup_fct_c))  &
        & CALL comin_plugin_finish("comin_plugin_primaryconstructor", &
        &                          "Cannot load primary constructor from plugin: " // convert_c_string(dlerror()))

      CALL C_F_PROCPOINTER( setup_fct_c, setup_fct )
      state%current_plugin => state%plugin_info(i)
      state%current_plugin%id = i
      IF (LEN_TRIM(plugin_list(i)%name) > 0) THEN
        state%current_plugin%name = TRIM(plugin_list(i)%name)
      ELSE
        last_sep_idx = SCAN(plugin_list(i)%plugin_library, "/", .TRUE.)
        state%current_plugin%name = TRIM(plugin_list(i)%plugin_library(last_sep_idx+1:)) // &
          & "(" // TRIM(plugin_list(i)%primary_constructor) // ")"
      ENDIF
      state%current_plugin%options = TRIM(plugin_list(i)%options)
      state%current_plugin%comm = TRIM(plugin_list(i)%comm)
      CALL setup_fct()
      NULLIFY(state%current_plugin)

    END DO

    ! after all primary callbacks are done: finalize callback structure and set flag
    CALL comin_callback_complete()
    CALL comin_var_complete()
    state%l_primary_done = .TRUE.
  END SUBROUTINE comin_plugin_primaryconstructor

  !> Performs basic compatibility checks.
  !! @ingroup host_interface
  SUBROUTINE comin_setup_check(plugin_str, wp_check)
    CHARACTER(LEN=*), INTENT(IN)  :: plugin_str !< plugin name
    INTEGER,          INTENT(IN)  :: wp_check !< KIND value for compatibility checks.
    !

    IF (.NOT. ASSOCIATED(state%comin_host_finish)) THEN
      CALL comin_error_set(COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET); RETURN
    END IF

    ! compare floating point precision
    IF (wp /= wp_check) THEN
      CALL comin_error_set(COMIN_ERROR_SETUP_PRECISION_TEST_FAILED); RETURN
    ELSE
      CALL comin_message("     " // plugin_str // ": working precision test successful.", 0)
    END IF
  END SUBROUTINE comin_setup_check

  !> Returns the structure `current_plugin`. It can for example be
  !> used to access the id of the current plugin.
  !! @ingroup plugin_interface
  SUBROUTINE comin_current_get_plugin_info(comin_current_plugin)
    TYPE(t_comin_plugin_info), INTENT(OUT)   :: comin_current_plugin !< plugin info struct

    comin_current_plugin = state%current_plugin

  END SUBROUTINE comin_current_get_plugin_info

  !> request plugin information, C interface
  INTEGER(C_INT) FUNCTION  comin_current_get_plugin_id() &
    & BIND(C, NAME="comin_current_get_plugin_id")
    comin_current_get_plugin_id = INT(state%current_plugin%id, C_INT)
  END FUNCTION comin_current_get_plugin_id

  !> request plugin information, C interface
  SUBROUTINE comin_current_get_plugin_name(val, len) &
    & BIND(C, NAME="comin_current_get_plugin_name")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< plugin name
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%current_plugin%name)
    len = LEN_TRIM(state%current_plugin%name)
  END SUBROUTINE comin_current_get_plugin_name

  !> request plugin information, C interface
  SUBROUTINE comin_current_get_plugin_options(val, len) &
    & BIND(C, NAME="comin_current_get_plugin_options")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< plugin options
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%current_plugin%options)
    len = LEN_TRIM(state%current_plugin%options)
  END SUBROUTINE comin_current_get_plugin_options

  !> request plugin information, C interface
  SUBROUTINE comin_current_get_plugin_comm(val, len) &
    & BIND(C, NAME="comin_current_get_plugin_comm")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< plugin comm. name
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length

    val = C_LOC(state%current_plugin%comm)
    len = LEN_TRIM(state%current_plugin%comm)
  END SUBROUTINE comin_current_get_plugin_comm

  !> Initialize the comin state
  !! @ingroup host_interface
  !! This routine needs to be called by the host before any other comin call
  !!
  SUBROUTINE comin_setup_init(lstdout)
    LOGICAL, INTENT(IN)  :: lstdout     !< do print on stdout or not
    IF (ASSOCIATED(state)) THEN
      ! <! cant use comin_message due to circular dependencies
      CALL comin_error_set(COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED); RETURN
    END IF
    ALLOCATE(state)
    state%lstdout = lstdout
  END SUBROUTINE comin_setup_init

  !> Destructor.
  !! @ingroup host_interface
  SUBROUTINE comin_setup_finalize()
    !
    INTEGER :: i
    INTEGER(c_int) :: ierr_c

    CALL comin_parallel_free_mpi_comms()
    DO i=1,SIZE(dl_handles)
      ierr_c = dlclose(dl_handles(i))
      IF (ierr_c /= 0) THEN
        CALL comin_error_set(COMIN_ERROR_SETUP_FINALIZE)
      END IF
    END DO
    DEALLOCATE(dl_handles)
    DEALLOCATE(state%plugin_info)
  END SUBROUTINE comin_setup_finalize

  !> Initialize the plugin state
  ! This routine is called by the host to set the plugins state
  ! explicitly to the state of the host model. It should by loaded by
  ! the host explicitly from the shared library of the plugin by using
  ! `dlsym`.
  SUBROUTINE comin_plugin_init(state_ptr, host_version, host_wp) &
    BIND(C)

    ! At the moment assuming only NVHPC and GCC (default way)
    TYPE(C_PTR), VALUE, INTENT(IN)               :: state_ptr
    TYPE(t_comin_setup_version_info), INTENT(IN) :: host_version
    INTEGER(C_INT), INTENT(IN)                   :: host_wp

    TYPE(t_comin_state), POINTER :: host_state => NULL()

    ! we cant rely on calling methods to obtain the version or wp,
    ! because this might call functions dynamically loaded by the host
    ! or other plugins. Hence we use the version preprocessor and wp
    ! constants directly.

    IF( .NOT. comin_setup_version_compatible(host_version, &
         & t_comin_setup_version_info(COMIN_VERSION_MAJOR, COMIN_VERSION_MINOR, COMIN_VERSION_PATCH))) THEN
      CALL comin_error_set(COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION); RETURN
    END IF

    IF( host_wp /= C_DOUBLE ) THEN
      CALL comin_error_set(COMIN_ERROR_PLUGIN_INIT_PRECISION); RETURN
    END IF

    CALL C_F_POINTER(state_ptr, host_state)

    IF(ASSOCIATED(state) .AND. .NOT. ASSOCIATED(state, host_state)) THEN
      ! <! cant use comin_message due to circular dependencies
      CALL comin_error_set(COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED); RETURN
    END IF

    state => host_state
  END SUBROUTINE comin_plugin_init

  !> Sets the "global-to-local" index lookup function.
  !! @ingroup host_interface
  SUBROUTINE comin_descrdata_set_fct_glb2loc_cell(fct)
    PROCEDURE(comin_glb2loc_index_lookup_fct) :: fct !< index lookup function
    state%comin_descrdata_fct_glb2loc_cell => fct
    IF (.NOT. ASSOCIATED(state%comin_descrdata_fct_glb2loc_cell)) THEN
      CALL comin_error_set(COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC); RETURN
    END IF
  END SUBROUTINE comin_descrdata_set_fct_glb2loc_cell

  !> Sets the global error handler procedure pointer.
  !! @ingroup host_interface
  !!
  !! To be called by the host application (ICON).
  !!
  !! Aborts with an error code if the error handler has already been
  !! set.
  SUBROUTINE comin_setup_errhandler(error_handler)
    PROCEDURE(comin_host_errhandler_fct) :: error_handler !< error handler

    state%comin_host_finish => error_handler
    IF (.NOT. ASSOCIATED(state%comin_host_finish)) THEN
      CALL comin_error_set(COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED); RETURN
    END IF
  END SUBROUTINE comin_setup_errhandler

END MODULE comin_setup
