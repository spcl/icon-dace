!> @file comin_variable
!! @brief Functions to modify and retrieve Variable definition
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_variable

  USE iso_c_binding,           ONLY: c_int, c_ptr, C_LOC, c_null_ptr, C_F_POINTER, C_BOOL
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS,                                          &
    &                                   COMIN_ERROR_POINTER_NOT_ASSOCIATED,                     &
    &                                   COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR,       &
    &                                   COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE,        &
    &                                   COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE,   &
    &                                   COMIN_ERROR_FIELD_NOT_ALLOCATED,                        &
    &                                   COMIN_ERROR_VARIABLE_NOT_TRACER_IN_TURBULENT_TRANSPORT, &
    &                                   COMIN_ERROR_VARIABLE_NOT_TRACER_IN_CONVECTIVE_TRANSPORT
  USE comin_errhandler,        ONLY: comin_plugin_finish, comin_error
  USE comin_setup_constants,   ONLY: wp, EP_SECONDARY_CONSTRUCTOR,                           &
    &                                EP_DESTRUCTOR
  USE comin_state,             ONLY: state
  USE comin_c_utils,           ONLY: convert_f_string, convert_c_string
  USE comin_parallel,          ONLY: comin_parallel_handle_mpi_errcode
  USE comin_variable_types,    ONLY: t_comin_var_descriptor, t_comin_var_item,               &
    &                                t_var_request_list_item, t_comin_var_context_item,      &
    &                                t_comin_var_context_item, t_comin_var_descriptor_c,     &
    &                                t_comin_var_descr_list_item, t_comin_var_context_item,  &
    &                                t_comin_var_ptr, t_var_list_item,                       &
    &                                t_var_context_list_item, t_comin_request_item,          &
    &                                comin_var_descr_match
  USE comin_descrdata_types,   ONLY: t_comin_descrdata_global
  USE comin_descrdata,         ONLY: comin_descrdata_get_global
  IMPLICIT NONE

  PRIVATE

  ! Public procedures, intention: called by host
  PUBLIC :: comin_var_list_finalize, comin_var_list_append
  PUBLIC :: comin_var_update
  PUBLIC :: comin_request_get_list_head
  ! Public procedures, intention: called by host and plugin
  PUBLIC :: comin_var_get_descr_list_head
  ! Public procedures, intention: called by plugin
  PUBLIC :: comin_var_request_add
  PUBLIC :: comin_var_get
  PUBLIC :: comin_var_get_descr_list_next, comin_var_get_descr_list_var_desc
  PUBLIC :: comin_var_to_3d
  ! Public procedures, intention: for internal use
  PUBLIC :: comin_var_complete
  PUBLIC :: comin_var_get_from_exposed




CONTAINS

  ! destructor.
  SUBROUTINE comin_var_item_finalize(this)
    TYPE(t_comin_var_item), INTENT(INOUT) :: this
    DEALLOCATE(this%descriptor%name)
    DEALLOCATE(this%p)
  END SUBROUTINE comin_var_item_finalize


  !> Returns head of variable list.
  !! @ingroup common
  !!
  !! We cannot directly use the `anylist` type-bound procedures from
  !! the user code (C compatibility).
  !!
  FUNCTION comin_var_get_descr_list_head()  RESULT(ptr)
    TYPE(t_comin_var_descr_list_item), POINTER :: ptr
    ptr => state%comin_var_descr_list%first()
  END FUNCTION comin_var_get_descr_list_head

  FUNCTION comin_var_get_descr_list_head_c() RESULT(ptr_c) &
    &  BIND(C, NAME="comin_var_get_descr_list_head")
    TYPE(C_PTR) :: ptr_c
    TYPE(t_comin_var_descr_list_item), POINTER :: ptr
    ptr_c = C_NULL_PTR
    ptr => state%comin_var_descr_list%first()
    IF(.NOT. ASSOCIATED(ptr)) THEN
       CALL comin_plugin_finish("Message of comin_var_get_descr_list_head_c", " ERROR: Pointer not associated.")
    END IF
    ptr_c = C_LOC(ptr)
  END FUNCTION comin_var_get_descr_list_head_c


  FUNCTION comin_var_get_descr_list_next(current) RESULT(ptr_c) BIND(C)
    TYPE(C_PTR), INTENT(IN), VALUE :: current
    TYPE(C_PTR) :: ptr_c
    TYPE(t_comin_var_descr_list_item), POINTER :: ptr, next
    ptr_c = C_NULL_PTR
    CALL C_F_POINTER(current, ptr)
    IF (.NOT. ASSOCIATED(ptr)) THEN
       CALL comin_plugin_finish("Message of comin_var_get_descr_list_next", " ERROR: Pointer not associated.")
    END IF
    next => ptr%next()
    ptr_c = C_LOC(next)
  END FUNCTION comin_var_get_descr_list_next


  SUBROUTINE comin_var_get_descr_list_var_desc(current, var_desc_out, ierr) BIND(C)
    TYPE(C_PTR), INTENT(IN), VALUE :: current
    TYPE(t_comin_var_descriptor_c), INTENT(INOUT) :: var_desc_out
    INTEGER(C_INT), INTENT(OUT) :: ierr
    TYPE(t_comin_var_descr_list_item), POINTER :: ptr => NULL()

    ierr = COMIN_SUCCESS

    CALL C_F_POINTER(current, ptr)
    IF (.NOT. ASSOCIATED(ptr)) THEN
       ierr = COMIN_ERROR_POINTER_NOT_ASSOCIATED
       RETURN
    END IF
    ASSOCIATE (var_item => ptr%item_value)
      var_desc_out%id = var_item%id
      CALL convert_f_string(var_item%name, var_desc_out%name)
    END ASSOCIATE
  END SUBROUTINE comin_var_get_descr_list_var_desc

  !> @return head of variable list.
  !! @ingroup host_interface
  FUNCTION comin_request_get_list_head()  RESULT(ptr)
    TYPE(t_var_request_list_item), POINTER :: ptr
    ptr => state%comin_var_request_list%first()
  END FUNCTION comin_request_get_list_head


  !> Append item to variable list.
  !! @ingroup host_interface
  SUBROUTINE comin_var_list_append(descriptor, p, ierr)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor !< variable descriptor
    TYPE(t_comin_var_ptr),        POINTER     :: p !< new variable  (pointer)
    INTEGER,                      INTENT(OUT) :: ierr !< error status code
    !
    TYPE(t_comin_var_item)                :: var_item
    TYPE(t_comin_var_descriptor), POINTER :: p_descr
    TYPE(t_comin_var_descr_list_item),  POINTER :: last

    ierr = COMIN_SUCCESS ! currently not used

    ! first, add the descriptor to a separate list
    ! (the one that is also exposed to the plugins):
    CALL state%comin_var_descr_list%append(t_comin_var_descr_list_item(descriptor))
    ! get a pointer to the item in the descriptor list
    last => state%comin_var_descr_list%firstptr%prevptr
    p_descr => last%item_value
    ! then add an entry to the other (internal) list of variables
    ! which contains a pointer to the above descriptor
    var_item%descriptor => p_descr
    var_item%p=>p
    CALL state%comin_var_list%append(t_var_list_item(var_item))
  END SUBROUTINE comin_var_list_append


  !> Destruct variable list, deallocate memory.
  !! @ingroup host_interface
  SUBROUTINE comin_var_list_finalize(ierr)
    INTEGER, INTENT(OUT) :: ierr !< error code
    ! local
    TYPE(t_var_list_item), POINTER :: p

    ierr = COMIN_SUCCESS
    p => state%comin_var_list%first()
    DO WHILE (ASSOCIATED(p))
      CALL comin_var_item_finalize(p%item_value)
      p => p%next()
    END DO
    CALL state%comin_var_list%delete_list()
  END SUBROUTINE comin_var_list_finalize


  FUNCTION comin_var_get_c(context_len, context, var_descriptor, flag) &
       & RESULT(var_pointer) &
       & BIND(C, name="comin_var_get")
    INTEGER(c_int),VALUE,           INTENT(IN) :: context_len
    INTEGER(c_int),                 INTENT(IN) :: context(context_len)
    TYPE(t_comin_var_descriptor_c), VALUE, INTENT(IN) :: var_descriptor
    INTEGER(c_int), VALUE,          INTENT(IN) :: flag
    TYPE(c_ptr)                                :: var_pointer
    TYPE(t_comin_var_ptr), POINTER             :: var_pointer_fortran

    TYPE(t_comin_var_descriptor) :: var_descriptor_fortran
    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id

    CALL comin_var_get(context, var_descriptor_fortran, flag, var_pointer_fortran)
    IF(ASSOCIATED(var_pointer_fortran)) THEN
       var_pointer = C_LOC(var_pointer_fortran)
    ELSE
       var_pointer = c_null_ptr
    ENDIF
  END FUNCTION comin_var_get_c


  FUNCTION comin_var_get_ptr(handle) &
       & RESULT(dataptr)                          &
       & BIND(C, NAME="comin_var_get_ptr")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    TYPE(C_PTR)                       :: dataptr
    !
    TYPE(t_comin_var_ptr), POINTER :: p => NULL()
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
       dataptr = C_NULL_PTR
    ELSE
       dataptr  = C_LOC(p%ptr)
    END IF
  END FUNCTION comin_var_get_ptr

  SUBROUTINE comin_var_get_shape(handle, data_shape, ierr) &
       & BIND(C, NAME="comin_var_get_shape")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    INTEGER(C_INT), INTENT(INOUT)       :: data_shape(5)
    INTEGER(C_INT), INTENT(OUT)       :: ierr
    !
    TYPE(t_comin_var_ptr), POINTER :: p => NULL()
    ierr = 0_C_INT
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
       ierr = COMIN_ERROR_POINTER_NOT_ASSOCIATED
       RETURN
    ELSE
       data_shape = SHAPE(p%ptr)
    END IF
  END SUBROUTINE comin_var_get_shape


  SUBROUTINE comin_var_get_pos(handle, pos_jc, pos_jk, pos_jb, pos_jn, ierr) &
       & BIND(C, NAME="comin_var_get_pos")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    INTEGER(C_INT), INTENT(OUT)       :: pos_jc
    INTEGER(C_INT), INTENT(OUT)       :: pos_jk
    INTEGER(C_INT), INTENT(OUT)       :: pos_jb
    INTEGER(C_INT), INTENT(OUT)       :: pos_jn
    INTEGER(C_INT), INTENT(OUT)       :: ierr
    !
    TYPE(t_comin_var_ptr), POINTER :: p => NULL()
    ierr = 0_C_INT
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
       ierr = COMIN_ERROR_POINTER_NOT_ASSOCIATED
       RETURN
    ELSE
       ! Convert to C dimension index
       pos_jc = p%pos_jc - 1
       pos_jk = p%pos_jk - 1
       pos_jb = p%pos_jb - 1
       pos_jn = p%pos_jn - 1
    END IF
  END SUBROUTINE comin_var_get_pos


  SUBROUTINE comin_var_get_ncontained(handle, ncontained, ierr) &
       & BIND(C, NAME="comin_var_get_ncontained")
    TYPE(C_PTR),    INTENT(IN), VALUE :: handle
    INTEGER(C_INT), INTENT(OUT)       :: ncontained
    INTEGER(C_INT), INTENT(OUT)       :: ierr
    !
    TYPE(t_comin_var_ptr), POINTER :: p => NULL()
    ierr = 0_C_INT
    CALL C_F_POINTER(handle, p)
    IF (.NOT. ASSOCIATED(p)) THEN
       ierr = COMIN_ERROR_POINTER_NOT_ASSOCIATED
       RETURN
    ELSE
       ! Convert to C dimension index
       ncontained = p%ncontained - 1
    END IF
  END SUBROUTINE comin_var_get_ncontained


  !> Request a pointer to an ICON variable in context(s).
  !! @ingroup plugin_interface
  SUBROUTINE comin_var_get(context, var_descriptor, flag, var_pointer)
    INTEGER,                      INTENT(IN) :: context(:)
    TYPE(t_comin_var_descriptor), INTENT(IN) :: var_descriptor
    INTEGER,                      INTENT(IN) :: flag
    TYPE(t_comin_var_ptr), POINTER           :: var_pointer
    ! local
    TYPE(t_comin_var_item), POINTER :: var_item
    TYPE(t_comin_var_context_item), POINTER :: var_list_element
    INTEGER :: ic

    var_pointer => NULL()

    IF (.NOT. state%l_primary_done) RETURN
    ! Routine should only be called during secondary constructor
    IF (state%current_ep > EP_SECONDARY_CONSTRUCTOR) RETURN

    ! first find the variable in list of all ICON variables and set the pointer
    var_item => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(var_item))  RETURN
    var_pointer => var_item%p

    ! Register the exposed variable in a second list, which is
    ! restricted to the context.
    IF (.NOT. ALLOCATED(state%comin_var_list_context)) THEN
      ALLOCATE(state%comin_var_list_context(EP_DESTRUCTOR, state%num_plugins))
    END IF
    DO ic = 1, SIZE(context)
      ! ignore EP_SECONDARY_CONSTRUCTOR for var_list
      IF (context(ic) == EP_SECONDARY_CONSTRUCTOR) CYCLE
      var_list_element => comin_var_get_by_context(context(ic), state%current_plugin%id, var_item%descriptor)
      IF (.NOT. ASSOCIATED(var_list_element)) THEN
        ! not in context list: register variable, set access flag
        ASSOCIATE(var_list => state%comin_var_list_context(context(ic) , state%current_plugin%id)%var_list)
          CALL var_list%append( t_var_context_list_item(t_comin_var_context_item( &
            &                                descriptor = var_item%descriptor,    &
            &                                metadata = var_item%metadata,        &
            &                                p = var_item%p,                      &
            &                                access_flag = flag)) )
        END ASSOCIATE
      END IF
    END DO
  END SUBROUTINE comin_var_get


  !> get pointer to a variable exposed by ICON
  FUNCTION comin_var_get_from_exposed(var_descriptor)  RESULT(comin_get_var)
    TYPE(t_comin_var_item), POINTER :: comin_get_var
    TYPE (t_comin_var_descriptor), INTENT(IN) :: var_descriptor
    !
    TYPE(t_var_list_item), POINTER :: p

    comin_get_var => NULL()
    p => state%comin_var_list%first()
    DO WHILE (ASSOCIATED(p))
      IF (comin_var_descr_match(p%item_value%descriptor, var_descriptor)) THEN
        comin_get_var => p%item_value
        EXIT
      END IF
      p => p%next()
    END DO
  END FUNCTION comin_var_get_from_exposed


  !> get pointer to a variable according to context and descriptor
  FUNCTION comin_var_get_by_context(context, plugin_id, var_descriptor)  RESULT(comin_get_var)
    TYPE(t_comin_var_context_item), POINTER   :: comin_get_var
    INTEGER, INTENT(IN)                       :: context, plugin_id
    TYPE (t_comin_var_descriptor), INTENT(IN) :: var_descriptor
    ! local
    TYPE(t_var_context_list_item), POINTER :: p

    comin_get_var => NULL()
    IF (.NOT. ALLOCATED(state%comin_var_list_context)) RETURN
    ASSOCIATE(var_list => state%comin_var_list_context(context, plugin_id)%var_list)
      p => var_list%first()
      DO WHILE (ASSOCIATED(p))
        ! test if already registered for context
        IF (comin_var_descr_match(p%item_value%descriptor, var_descriptor)) THEN
          comin_get_var => p%item_value
          EXIT
        END IF
        p => p%next()
      END DO
    END ASSOCIATE
  END FUNCTION comin_var_get_by_context


  !> subroutine to update a pointer
  !! @ingroup host_interface
  !!
  !! @note This subroutine aborts internally if the variable cannot be found.
  !!
  SUBROUTINE comin_var_update(var_descriptor, data_ptr)
    TYPE (t_comin_var_descriptor), INTENT(IN) :: var_descriptor
    REAL(wp), POINTER        :: data_ptr(:,:,:,:,:)
    ! local
    TYPE(t_comin_var_item), POINTER :: p

    IF (.NOT. ASSOCIATED(data_ptr)) THEN
      CALL comin_plugin_finish("comin_variable::comin_var_update", "Internal error! Unassociated pointer passed in.")
    END IF

    p => comin_var_get_from_exposed(var_descriptor)
    IF (.NOT. ASSOCIATED(p)) THEN
       CALL comin_plugin_finish("comin_variable::comin_var_update", "Internal error! Cannot find variable to update.")
    END IF

    p%p%ptr => data_ptr
  END SUBROUTINE comin_var_update


  SUBROUTINE comin_var_request_add_c(var_descriptor, lmodexclusive, ierr) &
    &  BIND(C, name="comin_var_request_add")
    TYPE (t_comin_var_descriptor_c), VALUE, INTENT(IN)  :: var_descriptor
    LOGICAL(C_BOOL), VALUE,          INTENT(IN)  :: lmodexclusive
    INTEGER(kind=C_INT),             INTENT(OUT) :: ierr
    !
    TYPE (t_comin_var_descriptor) :: var_descriptor_fortran

    var_descriptor_fortran%name = convert_c_string(var_descriptor%name)
    var_descriptor_fortran%id = var_descriptor%id
    CALL comin_var_request_add(var_descriptor_fortran, LOGICAL(lmodexclusive), ierr)
    IF (ierr /= COMIN_SUCCESS)  RETURN
  END SUBROUTINE comin_var_request_add_c


  !> By calling this subroutine inside the primary constructor, 3rd
  !> party plugins may request the creation of additional variables.
  !! @ingroup plugin_interface
  !!
  !!  Note: The lmodexclusive argument provides the information if this
  !!        variable is exclusive to the calling plugin.
  !!
  !!  Note: If a 3rd party plugin requests the creation of a variable
  !!        through this subroutine, it is still not guaranteed that
  !!        this variable is actually created! It might be skipped
  !!        due to inconsistencies, it could be a duplicate
  !!        etc. Therefore, 3rd party plugins still have to evaluate
  !!        the return code of `comin_var_request_add`.
  !!
  SUBROUTINE comin_var_request_add(var_descriptor, lmodexclusive, ierr)
    TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor
    LOGICAL,                       INTENT(IN)  :: lmodexclusive
    INTEGER,                       INTENT(OUT) :: ierr
    ! local
    TYPE(t_comin_descrdata_global),     POINTER :: comin_global
    TYPE (t_comin_var_descriptor)          :: var_descriptor_domain
    INTEGER                                :: domain_id

    comin_global => comin_descrdata_get_global()

    ierr = COMIN_SUCCESS

    IF (state%l_primary_done) THEN
       ierr = COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR
       CALL comin_error(ierr)
       RETURN
    ENDIF

    IF (var_descriptor%id == -1) THEN
      comin_global => comin_descrdata_get_global()
      IF (.NOT. ASSOCIATED(comin_global)) CALL comin_plugin_finish("variable ", "global data missing")

      DO domain_id = 1, comin_global%n_dom
        var_descriptor_domain    = var_descriptor
        var_descriptor_domain%id = domain_id
        CALL comin_var_request_add_element(var_descriptor_domain, lmodexclusive, ierr)
      ENDDO
    ELSE
      CALL comin_var_request_add_element(var_descriptor, lmodexclusive, ierr)
   ENDIF

  CONTAINS

    SUBROUTINE comin_var_request_add_element(var_descriptor, lmodexclusive, ierr)
      TYPE (t_comin_var_descriptor), INTENT(IN)  :: var_descriptor
      LOGICAL,                       INTENT(IN)  :: lmodexclusive
      INTEGER,                       INTENT(OUT) :: ierr
      !
      TYPE(t_var_request_list_item), POINTER :: p
      TYPE(t_comin_request_item):: comin_request_item

      ! check if requested variable already requested or if modexlusive conflicts exist
      ! first find the variable in list of all ICON variables and set the pointer
      p => state%comin_var_request_list%first()
      DO WHILE (ASSOCIATED(p))
        ASSOCIATE (var_list_request_element => p%item_value)
          IF (comin_var_descr_match(var_list_request_element%descriptor, var_descriptor)) THEN
            !> first criterion for abort: variable exists and was requested exclusively
            IF (var_list_request_element%lmodexclusive) THEN
              ierr = COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE; RETURN
              !> second criterion for abort: variable exists and now requested exclusively
            ELSEIF (lmodexclusive) THEN
              ierr = COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE; RETURN
              !> if existing but no conflicts with exclusiveness: expand moduleID information
            ELSE
              IF (.NOT. ALLOCATED(var_list_request_element%moduleID)) THEN
                ! if not allocated something went wrong before (should not happen)
                ierr = COMIN_ERROR_FIELD_NOT_ALLOCATED; RETURN
              ELSE
                var_list_request_element%moduleID = [var_list_request_element%moduleID(:), &
                  &                                  state%current_plugin%id]
              END IF
              ierr = COMIN_SUCCESS; RETURN
            END IF
          END IF
        END ASSOCIATE
        p => p%next()
      END DO

      ! register new variable request
      ASSOCIATE( var_list => state%comin_var_request_list)

        comin_request_item%descriptor =  var_descriptor
        comin_request_item%lmodexclusive = lmodexclusive
        comin_request_item%moduleID   = [state%current_plugin%id]
        CALL var_list%append( t_var_request_list_item(comin_request_item))
      END ASSOCIATE
    END SUBROUTINE comin_var_request_add_element

  END SUBROUTINE comin_var_request_add


  ! Internal subroutine. Consistency checks and similar operations,
  ! done after primary constructors.
  SUBROUTINE comin_var_complete(ierr)
    INTEGER, INTENT(OUT)    :: ierr
    TYPE(t_var_request_list_item),   POINTER :: p

    ierr = COMIN_SUCCESS
    ! consistency checks
    p => state%comin_var_request_list%first()
    DO WHILE (ASSOCIATED(p))
      ASSOCIATE (var_descriptor => p%item_value%descriptor, &
        &        var_metadata   => p%item_value%metadata)
        IF (var_metadata%tracer_turb .AND. .NOT. var_metadata%tracer) THEN
          ierr = COMIN_ERROR_VARIABLE_NOT_TRACER_IN_TURBULENT_TRANSPORT
          CALL comin_error(ierr, var_descriptor%name)
          RETURN
        ENDIF
        IF (var_metadata%tracer_conv .AND. .NOT. var_metadata%tracer) THEN
          ierr = COMIN_ERROR_VARIABLE_NOT_TRACER_IN_CONVECTIVE_TRANSPORT
          CALL comin_error(ierr, var_descriptor%name)
          RETURN
        ENDIF
      END ASSOCIATE
      p => p%next()
    END DO
  END SUBROUTINE comin_var_complete


  !> Convenience operation for accessing 2D/3D fields.
  !! @ingroup plugin_interface
  !!
  !! Assumes that the last dimension is not used!
  FUNCTION comin_var_to_3d(var) RESULT(slice)
    TYPE(t_comin_var_ptr), INTENT(IN)  :: var
    REAL(wp), POINTER :: slice(:,:,:)

    ! this operation is invalid if the field is a container
    IF (var%lcontainer) THEN
       CALL comin_plugin_finish("comin_var_to_3d", " ERROR: Attempt to convert container variable into 3D field.")
    END IF

    SELECT CASE (var%pos_jn)
       CASE(1)
          slice => var%ptr(1, :, :, :, 1)
       CASE(2)
          slice => var%ptr(:, 1, :, :, 1)
       CASE(3)
          slice => var%ptr(:, :, 1, :, 1)
       CASE DEFAULT
          slice => var%ptr(:, :, :, 1, 1)
       END SELECT
  END FUNCTION comin_var_to_3d

END MODULE comin_variable
