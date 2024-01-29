! --------------------------------------------------------------------
!> Example plugin for the ICON Community Interface (ComIn)
!  with basic (not MPI-parallel) callbacks and accessing variables and
!  descriptive data structures.
!
!  Note that in order to demonstrate ComIn's language interoperability,
!  a similary plugin has been implemented in C, see the subdirectory
!  "simple_c".
!
! 
!  @authors 01/2023 :: ICON Community INTERFACE  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
! --------------------------------------------------------------------
MODULE simple_fortran_plugin
  USE mpi
  USE iso_c_binding,           ONLY : C_INT
  USE comin_plugin_interface,  ONLY : comin_callback_register,                                          &
    &                                 comin_var_get, comin_parallel_get_host_mpi_comm,                  &
    &                                 t_comin_var_descriptor, t_comin_var_ptr,                          &
    &                                 comin_var_request_add,                                            &
    &                                 comin_descrdata_get_domain, t_comin_descrdata_domain,             &
    &                                 comin_descrdata_get_global, t_comin_descrdata_global,             &
    &                                 comin_descrdata_get_simulation_interval,                          &
    &                                 t_comin_descrdata_simulation_interval,                            &
    &                                 t_comin_setup_version_info, comin_setup_get_version,              &
    &                                 EP_SECONDARY_CONSTRUCTOR, EP_DESTRUCTOR,                          &
    &                                 EP_ATM_PHYSICS_BEFORE, EP_ATM_WRITE_OUTPUT_BEFORE,                &
    &                                 COMIN_FLAG_READ, COMIN_FLAG_WRITE, COMIN_ZAXIS_2D,                &
    &                                 comin_parallel_get_host_mpi_rank, comin_current_get_domain_id,    &
    &                                 COMIN_DOMAIN_OUTSIDE_LOOP, comin_callback_get_ep_name,            &
    &                                 t_comin_plugin_info, comin_current_get_plugin_info,               &
    &                                 comin_var_get_descr_list_head, t_comin_var_descr_list_item,       &
    &                                 comin_plugin_finish, comin_metadata_set,                          &
    &                                 comin_metadata_get, comin_var_to_3d, comin_descrdata_get_timesteplength

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: pluginname = "simple_fortran_plugin"

  !> working precision (will be compared to ComIn's and ICON's)
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  TYPE(t_comin_setup_version_info) :: version

  TYPE(t_comin_var_ptr),  POINTER :: pres, simple_fortran_var, simple_fortran_tracer
  INTEGER                         :: rank

  !> access descriptive data structures
  TYPE(t_comin_descrdata_domain),     POINTER   :: p_patch
  TYPE(t_comin_descrdata_global),     POINTER   :: p_global
  TYPE(t_comin_descrdata_simulation_interval), POINTER   :: p_simulation_interval
  
  TYPE :: t_comin_var_ptr_list
   TYPE(t_comin_var_ptr), POINTER :: var_ptr
  END TYPE t_comin_var_ptr_list
  
  TYPE(t_comin_var_ptr_list),  ALLOCATABLE :: qv(:)
  
CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main()  BIND(C)
    !
    TYPE(t_comin_plugin_info)     :: this_plugin
    TYPE(t_comin_var_descriptor)  :: simple_fortran_d, simple_fortran_tracer_d
    INTEGER                       :: ierr
    REAL(wp)                      :: dtime_1, dtime_2

    ierr = 0

    rank = comin_parallel_get_host_mpi_rank()
    IF (rank == 0)  WRITE (0,*) "     setup: ", pluginname

    version = comin_setup_get_version()
    IF (version%version_no_major > 1)  THEN
      CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", "incompatible version!")
    END IF

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin, ierr)
    IF (ierr /= 0) THEN
       CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", "Current plugin info could not be obtained!")
    ENDIF
    IF (rank == 0) WRITE (0,'(a,a,a,i4)') "     plugin ", this_plugin%name, " has id: ", this_plugin%id

    !> add requests for additional ICON variables

    ! request host model to add variable simple_fortran_var
    simple_fortran_d = t_comin_var_descriptor(id = 1, name = "simple_fortran_var")
    CALL comin_var_request_add(simple_fortran_d, .FALSE., ierr)
    IF (ierr /= 0) THEN
       CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
            & "Failed for comin_var_request_add to add variable simple_fortran_var!")
    ENDIF
    CALL comin_metadata_set(simple_fortran_d, "tracer", .FALSE., ierr)
    IF (ierr /= 0) THEN
       CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
            & "comin_metadata_set failed for tracer!")
    ENDIF
    CALL comin_metadata_set(simple_fortran_d, "restart", .FALSE., ierr)
    IF (ierr /= 0) THEN
       CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
            & "comin_metadata_set failed for restart!")
    ENDIF

    ! request host model to add tracer simple_fortran_tracer
    simple_fortran_tracer_d = t_comin_var_descriptor( id = -1, name = "simple_fortran_tracer" )
    CALL comin_var_request_add(simple_fortran_tracer_d, .FALSE., ierr)
    IF (ierr /= 0) THEN
       CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
            & "Failed for comin_var_request_add to add tracer simple_fortran_tracer!")
    ENDIF
    CALL comin_metadata_set(simple_fortran_tracer_d, "tracer", .TRUE., ierr)
    IF (ierr /= 0) THEN
       CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
            & "comin_metadata_set failed for tracer simple_fortran_tracer!")
    ENDIF
    CALL comin_metadata_set(simple_fortran_tracer_d, "restart", .FALSE., ierr)
    IF (ierr /= 0) THEN
       CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
            & "comin_metadata_set failed for restart (simple_fortran_tracer)!")
    ENDIF

    ! register callbacks
    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR,   simple_fortran_constructor, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
      & "comin_callback_register failed for simple_fortran_constructor!")
    CALL comin_callback_register(EP_ATM_WRITE_OUTPUT_BEFORE, simple_fortran_diagfct,     ierr)
    IF (ierr /= 0) CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
      & "comin_callback_register failed for simple_fortran_diagfct!")
    CALL comin_callback_register(EP_DESTRUCTOR,              simple_fortran_destructor,  ierr)
    IF (ierr /= 0) CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", &
      & "comin_callback_register failed for simple_fortran_destructor!")

    ! get descriptive data structures
    p_patch    => comin_descrdata_get_domain(1)
    p_global   => comin_descrdata_get_global()
    p_simulation_interval    => comin_descrdata_get_simulation_interval()

    dtime_1 = comin_descrdata_get_timesteplength(1)
    dtime_2 = comin_descrdata_get_timesteplength(2)
    IF (rank == 0) WRITE(0,*) "     timesteplength from comin_descrdata_get_timesteplength", dtime_1, dtime_2
  END SUBROUTINE comin_main


  ! --------------------------------------------------------------------
  ! ComIn secondary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE simple_fortran_constructor()  BIND(C)
    TYPE(t_comin_var_descriptor)         :: var_desc
    TYPE(t_comin_var_descr_list_item), POINTER :: var_descr_list_ptr
    INTEGER                              :: ierr, jg
    LOGICAL                              :: tracer
    REAL(WP), POINTER                    :: tracer_slice(:,:,:)

    IF (rank == 0) THEN
      WRITE (0,*) "     third party callback: secondary constructor."
      WRITE (0,*) "     third party callback: iterate over variable list:"
    END IF
    
    var_descr_list_ptr => comin_var_get_descr_list_head()
    DO WHILE (ASSOCIATED(var_descr_list_ptr))
      ASSOCIATE (descriptor => var_descr_list_ptr%item_value)
        IF (rank == 0)  WRITE (0,*) "Variable found: ", TRIM(descriptor%name), &
          "(", descriptor%id, ")"
      END ASSOCIATE
      var_descr_list_ptr => var_descr_list_ptr%next()
    END DO

    IF (rank == 0)  WRITE (0,*) "     ", pluginname, " - register some variables in some context"
    var_desc%name = 'pres'
    var_desc%id = 1
    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
      &                 var_desc, COMIN_FLAG_READ, pres)
    IF (.NOT. ASSOCIATED(pres)) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")

    IF (ANY([pres%pos_jc, pres%pos_jk, pres%pos_jb] /= [1,2,3])) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Dimension check failed!")

    ALLOCATE(qv(p_global%n_dom))
    DO jg =1, p_global%n_dom
          CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
          &                 t_comin_var_descriptor(name='qv', id=jg), COMIN_FLAG_READ, qv(jg)%var_ptr)
          IF (.NOT. ASSOCIATED(qv(jg)%var_ptr)) &
          &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")
          IF (ANY([qv(jg)%var_ptr%pos_jc, qv(jg)%var_ptr%pos_jk, qv(jg)%var_ptr%pos_jb] /= [1,2,3])) &
          &  CALL comin_plugin_finish("simple_fortran_constructor", "Dimension check failed!")
     END DO

    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
      &                t_comin_var_descriptor(name="simple_fortran_var", id=1), COMIN_FLAG_WRITE, simple_fortran_var)
    IF (.NOT. ASSOCIATED(simple_fortran_var)) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")

    CALL comin_metadata_get(t_comin_var_descriptor(name="simple_fortran_var", id=1), &
      &                     "tracer", tracer, ierr)
    IF (ierr/=0) &
      &  CALL comin_plugin_finish("comin_var_get_metadata_logical", "Internal error!")

    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
      &                t_comin_var_descriptor(name="simple_fortran_tracer", id=1), COMIN_FLAG_WRITE, simple_fortran_tracer)
    IF (.NOT. ASSOCIATED(simple_fortran_tracer)) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")

    ! access the "simple_fortran_tracer" via a 3D array pointer:
    tracer_slice => comin_var_to_3d(simple_fortran_tracer)
    tracer_slice(1,9,1) = 1.0
    WRITE(0,*) SHAPE(tracer_slice)
  END SUBROUTINE simple_fortran_constructor


  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE simple_fortran_diagfct()  BIND(C)
    TYPE(t_comin_plugin_info)       :: this_plugin
    INTEGER                         :: ierr, domain_id, jg, comm, root
    REAL(wp)                        :: local_max, global_max

    IF (rank == 0)  WRITE (0,*) "     third party callback: before output."

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin, ierr)
    IF (ierr /= 0)  CALL comin_plugin_finish("simple_fortran_diagfct (simple_fortran_plugin)", &
         & "Current plugin info could not be obtained!")
    IF (rank == 0)  WRITE (0,'(a,a,a,i4)') "      plugin ", this_plugin%name, " has id: ", this_plugin%id

    domain_id = comin_current_get_domain_id()

    IF (domain_id == COMIN_DOMAIN_OUTSIDE_LOOP) THEN
       IF (rank == 0)  WRITE(0,'(a)')  "      currently not in domain loop"
    ELSE
       IF (rank == 0)  WRITE(0,'(a,i0)')  "      currently on domain ", domain_id
    END IF

    simple_fortran_var%ptr = pres%ptr + 42

    DO jg =1, p_global%n_dom
      comm      = comin_parallel_get_host_mpi_comm()
      local_max = MAXVAL(qv(jg)%var_ptr%ptr)
      root      = 0
      CALL MPI_REDUCE(local_max, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, comm, ierr)
      IF (rank == root)   WRITE(0,*) "domain ", jg, ": global max = ", global_max
    END DO
  END SUBROUTINE simple_fortran_diagfct


  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE simple_fortran_destructor() BIND(C)
    IF (rank == 0)  WRITE (0,*) "     third party callback: destructor."
    NULLIFY(pres, simple_fortran_var, simple_fortran_tracer, p_patch, p_global)
  END SUBROUTINE simple_fortran_destructor

END MODULE simple_fortran_plugin
