!> @file comin_callback.F90
!! @brief Routines to handle third party plugin callbacks.
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_callback

  USE ISO_C_BINDING, ONLY: C_INT, C_CHAR
  USE comin_setup_constants,      ONLY: EP_DESTRUCTOR,             &
    &                                   DOMAIN_UNDEFINED, COMIN_FLAG_DEVICE, COMIN_FLAG_WRITE, COMIN_FLAG_READ
  USE comin_state,                ONLY: state
  USE comin_errhandler,           ONLY: comin_message, comin_error_set
  USE comin_errhandler_constants, ONLY: COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR, &
    &                                   COMIN_ERROR_CALLBACK_COMPLETE, COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN
  USE comin_callback_types,       ONLY: comin_callback_routine, t_comin_callback_element,      &
    &                                   t_callback_list_item
  USE comin_variable_types,       ONLY: t_comin_var_context_item, t_var_context_list_item
  USE comin_setup_constants,      ONLY: EP_NAME
  USE comin_c_utils,              ONLY: convert_f_string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_callback_register, comin_callback_context_call, comin_callback_complete
  PUBLIC :: comin_callback_get_ep_name

#include "comin_global.inc"

CONTAINS

  !> Routine to register new callbacks during primary constructor.
  !! @ingroup plugin_interface
  !!
  !! Also stores the currently active 3rd party plugin "current_plugin".
  SUBROUTINE comin_callback_register(entry_point_id, fct_ptr) &
    &  BIND(C)
    INTEGER(kind=C_INT), INTENT(IN), VALUE        :: entry_point_id
    PROCEDURE(comin_callback_routine) :: fct_ptr
    !
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name

    !> callbacks cannot be registered after primary constructor
    IF (state%l_primary_done) THEN
      CALL comin_error_set(COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR)
      RETURN
    ENDIF

    CALL state%comin_callback_list%append( t_callback_list_item(         &
      &  t_comin_callback_element(entry_point_id = entry_point_id, &
      &                           plugin_info    = state%current_plugin, &
      &                           comin_callback = fct_ptr) ))

    CALL comin_callback_get_ep_name(entry_point_id, ep_name)
    CALL comin_message("    registration for '"//ep_name//"' (ep: "//&
         &int2string(entry_point_id,'(i0)')//") associated for 3rd party plugin "//&
         &TRIM(state%current_plugin%name)//" successful.", 12)
  END SUBROUTINE comin_callback_register

  SUBROUTINE comin_callback_complete()
    ! local
    TYPE(t_callback_list_item), POINTER :: p
    INTEGER :: ep_loc, tp_loc
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name
    INTEGER                       :: status

    !> finalize settings made in primary constructor
    CALL comin_message("     Complete primary constructors", 0)

    !> convert callbacks from list to array, since now size known
    ALLOCATE(state%comin_callback_context(1:EP_DESTRUCTOR,1:state%num_plugins), stat=status)
    IF (status /= 0) THEN
      CALL comin_error_set(COMIN_ERROR_CALLBACK_COMPLETE); RETURN
    END IF
    ALLOCATE(state%comin_callback_order(1:EP_DESTRUCTOR,1:state%num_plugins), stat=status)
    IF (status /= 0) THEN
      CALL comin_error_set(COMIN_ERROR_CALLBACK_COMPLETE); RETURN
    END IF
    state%comin_callback_order = 0
    !> go iterate current list of entry points
    p => state%comin_callback_list%first()
    DO WHILE (ASSOCIATED(p))
      ASSOCIATE (var_list_element => p%item_value)

        ep_loc = var_list_element%entry_point_id
        tp_loc = var_list_element%plugin_info%id

        ASSOCIATE (callback_context => state%comin_callback_context(ep_loc, tp_loc))
          IF (.NOT. ASSOCIATED(callback_context%vl)) THEN
            ALLOCATE(callback_context%vl)
          ELSE
            CALL comin_callback_get_ep_name(ep_loc, ep_name)
            CALL comin_message("     WARNING:: Overwrite callback for plugin '"//&
                 &state%current_plugin%name//"' at entry point '"//ep_name//"' (ep: "//&
                 &int2string(ep_loc,'(i0)')//")", 0)
          END IF
          callback_context%vl = t_comin_callback_element(                &
            &              entry_point_id = ep_loc,                      &
            &              plugin_info    = var_list_element%plugin_info, &
            &              comin_callback = var_list_element%comin_callback)
          state%comin_callback_order(ep_loc, tp_loc) = tp_loc
        END ASSOCIATE
      END ASSOCIATE
      p => p%next()
    END DO

    !> delete linked list
    CALL state%comin_callback_list%delete_list()

    !> order/re-order callbacks
    ! Note: re-ordering based on namelist settings will be implemented
    !       in a later version of ComIn
    !       current default: order as at registration

    ! Note:
    ! the adapter library checks for duplicates and exclusiveness of requested
    ! variables and potentially aborts directly in comin_var_request_add
    ! therefore no further check before secondary constructor required

  END SUBROUTINE comin_callback_complete

  !> Routine to find callback routine associated with current entry point
  !! @ingroup host_interface
  RECURSIVE SUBROUTINE comin_callback_context_call(entry_point_id, domain_id, lacc)
    INTEGER, INTENT(IN)            :: entry_point_id
    INTEGER, INTENT(IN)            :: domain_id
    LOGICAL, INTENT(IN)            :: lacc
    TYPE(t_comin_callback_element), POINTER :: cl
    !PROCEDURE(comin_callback_routine), POINTER :: comin_callback_context
    INTEGER   :: thirdpi, loci
    LOGICAL :: lcallbacks_exist
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name

    ! We cant call callbacks before the primary constructors are done
    ! and comin_callback_complete was called. (e.g. EP_FINISH)
    if(.NOT. state%l_primary_done) RETURN

    CALL comin_callback_get_ep_name(entry_point_id, ep_name)
    CALL comin_message("     CONTEXT " // ep_name, 12)

    !> call callback functions given by order in entry_point_order
    lcallbacks_exist = SIZE(state%comin_callback_order,2) > 0
    IF (lcallbacks_exist)  lcallbacks_exist = (SUM(state%comin_callback_order(entry_point_id,:)) /= 0)

    IF (lcallbacks_exist) THEN
      DO thirdpi=1,state%num_plugins
        loci = state%comin_callback_order(entry_point_id, thirdpi)
        NULLIFY(cl)
        IF (loci > 0) cl => state%comin_callback_context(entry_point_id,loci)%vl
        IF (ASSOCIATED(cl)) THEN
          CALL comin_message("     current ep '"//ep_name//"' (ep: "//&
               &int2string(cl%entry_point_id,'(i0)')//") for library: "//cl%plugin_info%name, 0)
          ! undefine domain id for the call of the ComIn destructor
          IF (entry_point_id == EP_DESTRUCTOR) state%current_domain_id = DOMAIN_UNDEFINED
          IF (.NOT. lacc) CALL check_var_no_device(entry_point_id, thirdpi)
          IF (lacc) CALL sync_vars_for_device(entry_point_id, thirdpi, COMIN_FLAG_READ)
          ! set current plugin
          state%current_plugin => cl%plugin_info
          ! set current entry point
          state%current_ep = entry_point_id
          ! set current domain id
          state%current_domain_id = domain_id
          CALL cl%comin_callback
          NULLIFY(state%current_plugin)
          IF (lacc) CALL sync_vars_for_device(entry_point_id, thirdpi, COMIN_FLAG_WRITE)
        ELSE
          CALL comin_message("      entry point '"//ep_name//"' (ep: "//&
               &int2string(entry_point_id,'(i0)')//") not associated", 12)
        END IF
      ENDDO
    ELSE
      CALL comin_message("     no calls associated with entry point '"//&
           &ep_name//"' (ep: "//int2string(entry_point_id,'(i0)')//").", 12)
    END IF
    DEALLOCATE(ep_name)
  END SUBROUTINE comin_callback_context_call

  SUBROUTINE sync_vars_for_device(ep, plugin_id, rw_flag)
    INTEGER, INTENT(IN) :: ep
    INTEGER, INTENT(IN) :: plugin_id
    INTEGER, INTENT(IN) :: rw_flag
    ! locals
    TYPE(t_var_context_list_item), POINTER :: p

    p => state%comin_var_list_context(ep , plugin_id)%var_list%first()
    DO WHILE (ASSOCIATED(p))
      IF (IAND(p%item_value%access_flag, COMIN_FLAG_DEVICE) == 0 .AND. &
 &        IAND(p%item_value%access_flag, rw_flag) /= 0) THEN
        CALL state%sync_device_mem(p%item_value%p, rw_flag == COMIN_FLAG_WRITE)
      END IF
      p => p%next()
    END DO
  END SUBROUTINE sync_vars_for_device

  SUBROUTINE check_var_no_device(ep, plugin_id)
    INTEGER, INTENT(IN) :: ep
    INTEGER, INTENT(IN) :: plugin_id
    ! locals
    TYPE(t_var_context_list_item), POINTER :: p
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name

    CALL comin_callback_get_ep_name(ep, ep_name)

    p => state%comin_var_list_context(ep , plugin_id)%var_list%first()
    DO WHILE (ASSOCIATED(p))
      IF (IAND(p%item_value%access_flag, COMIN_FLAG_DEVICE) /= 0) THEN
        CALL comin_message("WARNING: Device access at a non-ported entrypoint (" // &
             &             ep_name // &
             &             ") requested.", 1)
      END IF
      p => p%next()
    END DO
  END SUBROUTINE check_var_no_device

  ! returns integer n as a string (needed in printing messages)
  FUNCTION int2string(n, opt_fmt)
    CHARACTER(:), ALLOCATABLE :: int2string ! result
    CHARACTER(len=128) :: res
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=128) :: fmt

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(i0)'
    END IF
    WRITE(res,fmt) n
    res = ADJUSTL(res)
    int2string = TRIM(res)
  END FUNCTION int2string

  !> returns entry point name (character string) corresponding to `iep`.
  !! @ingroup plugin_interface
  SUBROUTINE comin_callback_get_ep_name( iep, out_ep_name )
    INTEGER, INTENT(IN)  :: iep   !< entry point ID
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: out_ep_name !< entry point name string

    IF ((iep < 0) .OR. (iep > EP_DESTRUCTOR)) THEN
      CALL comin_error_set(COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN); RETURN
      out_ep_name = "UNKNOWN"
    ELSE
      out_ep_name = TRIM(EP_NAME(iep))
    END IF
  END SUBROUTINE comin_callback_get_ep_name

  !> C-wrapper. Assumes that the `out_ep_name` has been allocated by
  !> the caller.
  !
  SUBROUTINE comin_callback_get_ep_name_c( iep, out_ep_name) &
    &     BIND(C, name="comin_callback_get_ep_name")
    INTEGER(c_int), VALUE, INTENT(IN)  :: iep   !< entry point ID
    CHARACTER(len=1, kind=c_char),DIMENSION(MAX_LEN_EP_NAME+1)  :: out_ep_name
    !
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name

    CALL comin_callback_get_ep_name(iep, ep_name)
    CALL convert_f_string(ep_name, out_ep_name)
  END SUBROUTINE comin_callback_get_ep_name_c

END MODULE comin_callback
