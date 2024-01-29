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
  USE comin_setup_utils,     ONLY: t_comin_setup_version_info, comin_setup_get_version
  USE comin_state,           ONLY: state
  USE comin_c_utils,         ONLY: convert_c_string, convert_f_string
  USE comin_parallel,        ONLY: comin_parallel_free_mpi_comms
  USE comin_plugin_types,    ONLY: t_comin_plugin_description,        &
    &                              t_comin_plugin_info, t_comin_plugin_info_c
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS,                                      &
    &                                   COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET,               &
    &                                   COMIN_ERROR_SETUP_PRECISION_TEST_FAILED,            &
    &                                   COMIN_ERROR_SETUP_FINALIZE
  USE comin_errhandler,      ONLY: comin_message, comin_error, comin_plugin_finish

  USE iso_c_binding,         ONLY: c_int, c_ptr, c_funptr, c_f_procpointer, c_null_char, &
    &                              C_ASSOCIATED, c_null_ptr, c_loc, c_char
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_setup_check
  PUBLIC :: comin_plugin_primaryconstructor, comin_setup_finalize

  PUBLIC :: comin_current_get_plugin_info

  !> interface for a primary constructor call
  ABSTRACT INTERFACE
     SUBROUTINE comin_plugin_init_fct(state_ptr, host_version, host_wp, ierr) BIND(C)
      IMPORT c_ptr, c_int, t_comin_setup_version_info
      TYPE(C_PTR), INTENT(IN), VALUE               :: state_ptr
      TYPE(t_comin_setup_version_info), INTENT(IN) :: host_version
      INTEGER(C_INT), INTENT(IN)                   :: host_wp
      INTEGER(C_INT), INTENT(OUT)                  :: ierr
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
 SUBROUTINE comin_plugin_primaryconstructor(plugin_list, ierr)
   !> list of dynamic libs:
   TYPE(t_comin_plugin_description), INTENT(IN), TARGET :: plugin_list(:)
   INTEGER, INTENT(OUT) :: ierr        !< error code
   !
   INTEGER :: i, last_sep_idx
   TYPE(c_funptr) :: setup_fct_c, plugin_init_fct_c
   PROCEDURE(comin_primaryconstructor_fct), BIND(C), POINTER :: setup_fct
   PROCEDURE(comin_plugin_init_fct), BIND(C), POINTER :: plugin_init_fct
   state%current_plugin = t_comin_plugin_info(id=-1)

   ierr  = COMIN_SUCCESS

   state%num_plugins = SIZE(plugin_list)
   ALLOCATE(dl_handles(state%num_plugins))
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
     CALL plugin_init_fct(C_LOC(state), comin_setup_get_version(), wp, ierr)
     IF(ierr /= COMIN_SUCCESS)  &
       & CALL comin_plugin_finish("comin_plugin_primaryconstructor", &
       &                          "Plugin initialization failed")

     setup_fct_c = dlsym(dl_handles(i), TRIM(plugin_list(i)%primary_constructor)//c_null_char)
     IF (.NOT. C_ASSOCIATED(setup_fct_c))  &
       & CALL comin_plugin_finish("comin_plugin_primaryconstructor", &
       &                          "Cannot load primary constructor from plugin: " // convert_c_string(dlerror()))

     CALL C_F_PROCPOINTER( setup_fct_c, setup_fct )
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

   END DO

   ! after all primary callbacks are done: finalize callback structure and set flag
   CALL comin_callback_complete(ierr)
   CALL comin_var_complete(ierr)
   state%l_primary_done = .TRUE.
  END SUBROUTINE comin_plugin_primaryconstructor

  !> Performs basic compatibility checks.
  !! @ingroup host_interface
  SUBROUTINE comin_setup_check(plugin_str, wp_check, ierr)
    CHARACTER(LEN=*), INTENT(IN)  :: plugin_str !< plugin name
    INTEGER,          INTENT(IN)  :: wp_check !< KIND value for compatibility checks.
    INTEGER,          INTENT(OUT) :: ierr !< error status code
    !
    ierr = COMIN_SUCCESS

    IF (.NOT. ASSOCIATED(state%comin_host_finish)) THEN
      ierr = COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET
      CALL comin_error(ierr, plugin_str)
    END IF

    ! compare floating point precision
    IF (wp /= wp_check) THEN
      ierr = COMIN_ERROR_SETUP_PRECISION_TEST_FAILED
      CALL comin_error(ierr, plugin_str)
    ELSE
      CALL comin_message("     " // plugin_str // ": working precision test successful.", 0)
    END IF
  END SUBROUTINE comin_setup_check

  !> Returns the structure `current_plugin`. It can for example be
  !> used to access the id of the current plugin.
  !! @ingroup plugin_interface
  SUBROUTINE comin_current_get_plugin_info(comin_current_plugin, ierr)
     TYPE(t_comin_plugin_info), INTENT(OUT)   :: comin_current_plugin !< plugin info struct
     INTEGER, INTENT(OUT)                     :: ierr                 !< error status flag

     ierr = COMIN_SUCCESS
     comin_current_plugin = state%current_plugin

  END SUBROUTINE comin_current_get_plugin_info


  !> request plugin information, C interface
  INTEGER(C_INT) FUNCTION  comin_current_get_plugin_id() &
    & BIND(C, NAME="comin_current_get_plugin_id")
    comin_current_get_plugin_id = INT(state%current_plugin%id, C_INT)
  END FUNCTION comin_current_get_plugin_id

  !> request plugin information, C interface
  SUBROUTINE comin_current_get_plugin_name(val, len, ierr) &
    & BIND(C, NAME="comin_current_get_plugin_name")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< plugin name
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
    val = C_LOC(state%current_plugin%name)
    len = LEN_TRIM(state%current_plugin%name)
  END SUBROUTINE comin_current_get_plugin_name


  !> request plugin information, C interface
  SUBROUTINE comin_current_get_plugin_options(val, len, ierr) &
    & BIND(C, NAME="comin_current_get_plugin_options")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< plugin options
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
    val = C_LOC(state%current_plugin%options)
    len = LEN_TRIM(state%current_plugin%options)
  END SUBROUTINE comin_current_get_plugin_options


  !> request plugin information, C interface
  SUBROUTINE comin_current_get_plugin_comm(val, len, ierr) &
    & BIND(C, NAME="comin_current_get_plugin_comm")
    TYPE(c_ptr),                    INTENT(OUT) :: val  !< plugin comm. name
    INTEGER(kind=c_int),            INTENT(OUT) :: len  !< string length
    INTEGER(kind=c_int),            INTENT(OUT) :: ierr !< error status code

    ierr = COMIN_SUCCESS ! currently unused
    val = C_LOC(state%current_plugin%comm)
    len = LEN_TRIM(state%current_plugin%comm)
  END SUBROUTINE comin_current_get_plugin_comm


  !> Destructor.
  !! @ingroup host_interface
  SUBROUTINE comin_setup_finalize(ierr)
    INTEGER, INTENT(OUT) :: ierr !< error code
    !
    INTEGER :: i
    INTEGER(c_int) :: ierr_c

    ierr = COMIN_SUCCESS
    CALL comin_parallel_free_mpi_comms()
    DO i=1,SIZE(dl_handles)
      ierr_c = dlclose(dl_handles(i))
      IF (ierr_c /= 0) THEN
        ierr = COMIN_ERROR_SETUP_FINALIZE
        EXIT
      END IF
    END DO
    DEALLOCATE(dl_handles)
  END SUBROUTINE comin_setup_finalize

END MODULE comin_setup
