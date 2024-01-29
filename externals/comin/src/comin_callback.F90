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

  USE ISO_C_BINDING, ONLY: C_INT
  USE comin_setup_constants,      ONLY: comin_callback_get_ep_name, EP_DESTRUCTOR,             &
    &                                   DOMAIN_UNDEFINED
  USE comin_state,                ONLY: state
  USE comin_errhandler,           ONLY: comin_message, comin_error
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS,                                            &
    &                                   COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR, &
    &                                   COMIN_ERROR_CALLBACK_COMPLETE
  USE comin_callback_types,       ONLY: comin_callback_routine, t_comin_callback_element,      &
    &                                   t_callback_list_item

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comin_callback_register, comin_callback_context_call, comin_callback_complete

CONTAINS

  !> Routine to register new callbacks during primary constructor.
  !! @ingroup plugin_interface
  !!
  !! Also stores the currently active 3rd party plugin "current_plugin".
  SUBROUTINE comin_callback_register(entry_point_id, fct_ptr, ierr) &
    &  BIND(C)
    INTEGER(kind=C_INT), INTENT(IN), VALUE        :: entry_point_id
    PROCEDURE(comin_callback_routine) :: fct_ptr
    INTEGER(kind=C_INT), INTENT(OUT):: ierr
    !
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name

    ierr = COMIN_SUCCESS

    !> callbacks cannot be registered after primary constructor
    IF (state%l_primary_done) THEN
       ierr = COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR
       CALL comin_error(ierr)
       RETURN
    ENDIF

    CALL state%comin_callback_list%append( t_callback_list_item(         &
      &  t_comin_callback_element(entry_point_id = entry_point_id, &
      &                           thirdparty     = state%current_plugin, &
      &                           comin_callback = fct_ptr) ))

    CALL comin_callback_get_ep_name(entry_point_id, ep_name, ierr)
    IF (ierr /= COMIN_SUCCESS) THEN
       ep_name = "UNKNOWN"
       ierr = COMIN_SUCCESS
    END IF
    CALL comin_message("    registration for '"//ep_name//"' (ep: "//&
         &int2string(entry_point_id,'(i0)')//") associated for 3rd party plugin "//&
         &TRIM(state%current_plugin%name)//" successful.", 12)
  END SUBROUTINE comin_callback_register


  SUBROUTINE comin_callback_complete(ierr)
    INTEGER, INTENT(OUT) :: ierr
    ! local
    TYPE(t_callback_list_item), POINTER :: p
    INTEGER :: ep_loc, tp_loc, ep_name_ierr
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name
    INTEGER                       :: status

    ierr = COMIN_SUCCESS
    ep_name_ierr = COMIN_SUCCESS

    !> finalize settings made in primary constructor
    CALL comin_message("     Complete primary constructors", 0)

    !> convert callbacks from list to array, since now size known
    ALLOCATE(state%comin_callback_context(1:EP_DESTRUCTOR,1:state%num_plugins), stat=status)
    IF (status /= 0) THEN
      ierr = COMIN_ERROR_CALLBACK_COMPLETE; RETURN
    END IF
    ALLOCATE(state%comin_callback_order(1:EP_DESTRUCTOR,1:state%num_plugins), stat=status)
    IF (status /= 0) THEN
      ierr = COMIN_ERROR_CALLBACK_COMPLETE; RETURN
    END IF
    state%comin_callback_order = 0
    !> go iterate current list of entry points
    p => state%comin_callback_list%first()
    DO WHILE (ASSOCIATED(p))
      ASSOCIATE (var_list_element => p%item_value)

        ep_loc = var_list_element%entry_point_id
        tp_loc = var_list_element%thirdparty%id

        ASSOCIATE (callback_context => state%comin_callback_context(ep_loc, tp_loc))
          IF (.NOT. ASSOCIATED(callback_context%vl)) THEN
            ALLOCATE(callback_context%vl)
         ELSE
            CALL comin_callback_get_ep_name(ep_loc, ep_name, ep_name_ierr)
            IF (ep_name_ierr /= COMIN_SUCCESS) ep_name = "UNKNOWN"
            CALL comin_message("     WARNING:: Overwrite callback for plugin '"//&
                 &state%current_plugin%name//"' at entry point '"//ep_name//"' (ep: "//&
                 &int2string(ep_loc,'(i0)')//")", 0)
         END IF
          callback_context%vl = t_comin_callback_element(                &
            &              entry_point_id = ep_loc,                      &
            &              thirdparty     = var_list_element%thirdparty, &
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
  SUBROUTINE comin_callback_context_call(entry_point_id, domain_id)
    INTEGER, INTENT(IN)            :: entry_point_id
    INTEGER, INTENT(IN)            :: domain_id
    TYPE(t_comin_callback_element), POINTER :: cl
    !PROCEDURE(comin_callback_routine), POINTER :: comin_callback_context
    INTEGER   :: thirdpi, loci, ierr
    LOGICAL :: lcallbacks_exist
    CHARACTER(LEN=:), ALLOCATABLE :: ep_name

    ! We cant call callbacks before the primary constructors are done
    ! and comin_callback_complete was called. (e.g. EP_FINISH)
    if(.NOT. state%l_primary_done) RETURN

    CALL comin_callback_get_ep_name(entry_point_id, ep_name, ierr)
    IF (ierr /= COMIN_SUCCESS) THEN
       ep_name = "UNKNOWN"
       ierr = COMIN_SUCCESS
    END IF
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
               &int2string(cl%entry_point_id,'(i0)')//") for library: "//cl%thirdparty%name, 0)
          ! set current plugin
          state%current_plugin = cl%thirdparty
          ! set current entry point
          state%current_ep = entry_point_id
          ! set current domain id
          state%current_domain_id = domain_id
          ! undefine domain id for the call of the ComIn destructor
          IF (entry_point_id == EP_DESTRUCTOR) state%current_domain_id = DOMAIN_UNDEFINED
          CALL cl%comin_callback
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

END MODULE comin_callback
