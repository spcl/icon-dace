!> @file comin_replay.F90
!! @brief Imitates a host model by reading data from a netcdf file
!! (created with the `comin_run_recorder_plugin`) and plays the data back to comin.
!
!  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

PROGRAM comin_replay
  USE mpi,                     ONLY: MPI_INIT, MPI_COMM_WORLD, MPI_COMM_RANK, MPI_COMM_SIZE,        &
       &                             MPI_FINALIZE, MPI_SUCCESS, MPI_ABORT
  USE comin_host_interface,    ONLY: t_comin_var_ptr, t_comin_setup_version_info,                   &
       &                             t_comin_descrdata_global, t_comin_descrdata_domain,            &
       &                             t_comin_descrdata_simulation_interval,                         &
       &                             t_comin_plugin_description, comin_setup_init,                  &
       &                             comin_setup_errhandler, comin_setup_set_verbosity_level,       &
       &                             comin_setup_check, comin_setup_get_version,                    &
       &                             comin_descrdata_set_global, comin_descrdata_set_domain,        &
       &                             comin_descrdata_set_timesteplength,                            &
       &                             comin_descrdata_set_simulation_interval,                       &
       &                             comin_parallel_mpi_handshake, comin_plugin_primaryconstructor, &
       &                             comin_current_set_datetime, comin_callback_context_call,       &
       &                             t_comin_var_descriptor, comin_request_get_list_head,           &
       &                             comin_var_list_append, comin_metadata_set,                     &
       &                             t_var_request_list_item, COMIN_ZAXIS_2D,                       &
       &                             COMIN_ZAXIS_3D, COMIN_ZAXIS_3D_HALF,mpi_handshake, EP_FINISH,  &
       &                             COMIN_DOMAIN_OUTSIDE_LOOP, EP_ATM_YAC_DEFCOMP_AFTER,           &
       &                             EP_ATM_YAC_SYNCDEF_AFTER, EP_ATM_YAC_ENDDEF_AFTER,             &
       &                             EP_ATM_TIMELOOP_BEFORE, t_comin_var_metadata_iterator,         &
       &                             comin_metadata_get_iterator, COMIN_METADATA_TYPEID_INTEGER,    &
       &                             COMIN_METADATA_TYPEID_REAL, COMIN_METADATA_TYPEID_CHARACTER,   &
       &                             COMIN_METADATA_TYPEID_LOGICAL, comin_metadata_get_or
  USE netcdf,                  ONLY: nf90_open, nf90_get_att, nf90_inq_ncid, nf90_get_var,          &
       &                             nf90_close, NF90_GLOBAL, NF90_NOWRITE
  USE netcdf_utils,            ONLY: nf90, nf90_utils_def_var, nf90_utils_get_shape
  USE comin_descrdata_load,    ONLY: comin_descrdata_load_domain, comin_descrdata_load_global
  USE iso_c_binding,           ONLY: C_DOUBLE, C_NULL_PTR
  USE utils,                   ONLY: int2string
#ifdef ENABLE_YAC
  USE yac,                     ONLY: yac_finit_comm, yac_fdef_calendar, YAC_PROLEPTIC_GREGORIAN, &
       &                             yac_fdef_datetime, yac_fdef_comp, yac_fsync_def, yac_fenddef
#endif

  IMPLICIT NONE

  TYPE :: replay_var
    REAL(C_DOUBLE), POINTER :: data_ptr(:,:,:,:,:)
    TYPE(t_comin_var_ptr) :: comin_var_ptr
  END TYPE replay_var

  INTEGER, PARAMETER :: wp = C_DOUBLE
  INTEGER :: jg, ierr
  INTEGER :: ncid, grp_ncid, varid
  INTEGER :: host_comm, comin_comm, host_rank, host_size
  INTEGER :: file_host_rank, file_host_size
  INTEGER :: i, nplugins
  CHARACTER(len=256) :: prefix = "", filename
  TYPE(t_comin_setup_version_info) :: comin_version
  INTEGER :: file_comin_version(3)

  TYPE(t_comin_descrdata_global)                  :: comin_global
  TYPE(t_comin_descrdata_domain), ALLOCATABLE     :: comin_domain(:)
  TYPE(t_comin_descrdata_simulation_interval)     :: comin_simulation_interval

  TYPE(replay_var), TARGET, ALLOCATABLE :: vars(:)

  TYPE(t_comin_plugin_description) :: plugin_list(16) !< list of dynamic libs (max: 16)
  NAMELIST /comin_nml/ plugin_list

  INTEGER :: shap(2)
  INTEGER, ALLOCATABLE :: current_ep(:)
  INTEGER, ALLOCATABLE :: current_domain_id(:)
  CHARACTER(LEN=32), ALLOCATABLE :: current_datetime(:)
  REAL(C_DOUBLE) :: dt
#ifdef ENABLE_YAC
  INTEGER :: yac_instance_id, yac_comp_id
#endif

  CALL start_mpi()

  IF ( COMMAND_ARGUMENT_COUNT() == 1 ) THEN
    CALL GET_COMMAND_ARGUMENT(1,prefix)
  ENDIF
  filename = TRIM(prefix) // int2string(host_rank) // ".nc"
  IF(host_rank == 0) WRITE (0,*) "reading ", TRIM(filename)
  CALL nf90(nf90_open(TRIM(filename), NF90_NOWRITE, ncid))

  CALL comin_setup_init(host_rank == 0)
  CALL comin_setup_errhandler(finish)
  CALL comin_setup_set_verbosity_level(iverbosity=13)
  CALL comin_setup_check("replay", wp)

  ! check version of file and library
  comin_version = comin_setup_get_version()
  CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "comin_version", file_comin_version))
  IF(host_rank == 0) WRITE (0,*) "using comin version          ", comin_version
  IF(host_rank == 0) WRITE (0,*) "file was written with version ", file_comin_version
  IF(file_comin_version(1) /= comin_version%version_no_major .OR. &
     file_comin_version(2) /= comin_version%version_no_minor) &
    CALL finish("replay", "Incompatible versions")

  ! check host comm
  CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "host_comm_rank", file_host_rank))
  IF(host_rank /= file_host_rank) CALL finish("comin_replay", "Wrong host rank number in file")
  CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "host_comm_size", file_host_size))
  IF(host_size /= file_host_size) CALL finish("comin_replay", "Wrong host comm size in file")

  CALL comin_descrdata_load_global(ncid, comin_global)
#ifdef ENABLE_YAC
  IF (comin_global%yac_instance_id /= -1) THEN
    comin_global%yac_instance_id = yac_instance_id
  ENDIF
#endif
  CALL comin_descrdata_set_global(comin_global)
  ALLOCATE(comin_domain(comin_global%n_dom))
  DO jg = 1,comin_global%n_dom
    CALL nf90(nf90_inq_ncid(ncid, "domain_"//int2string(jg), grp_ncid))
    CALL comin_descrdata_load_domain(grp_ncid, comin_domain(jg))
    CALL comin_descrdata_set_domain(comin_domain)
    CALL nf90(nf90_get_att(grp_ncid, NF90_GLOBAL, "timestep_length", dt))
    CALL comin_descrdata_set_timesteplength(jg, dt)
  END DO

  CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "exp_start", comin_simulation_interval%exp_start))
  CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "exp_stop", comin_simulation_interval%exp_stop))
  CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "run_start", comin_simulation_interval%run_start))
  CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "run_stop", comin_simulation_interval%run_stop))
  CALL comin_descrdata_set_simulation_interval(comin_simulation_interval)
#ifdef ENABLE_YAC
  CALL yac_fdef_datetime(yac_instance_id, &
       & TRIM(comin_simulation_interval%run_start), &
       & TRIM(comin_simulation_interval%run_stop))
#endif

  OPEN(1, file="master.nml")
  READ(1, nml=comin_nml)
  CLOSE(1)
  DO i=1,SIZE(plugin_list)
    ! plugin and primary_constructor are optional. Only if both are not specified the plugin is considered as empty
    IF (len_TRIM(plugin_list(i)%plugin_library) == 0 .AND. &
    &   TRIM(plugin_list(i)%primary_constructor) == "comin_main") EXIT
  END DO

  nplugins = i-1
  CALL comin_parallel_mpi_handshake(comin_comm, plugin_list(1:nplugins)%comm, "atmo") ! currently the minimial example only emulates the atmo

  CALL comin_plugin_primaryconstructor(plugin_list(1:nplugins))

  CALL add_variables(comin_global, comin_domain)

  ! run callback loop
  shap(1:1) = nf90_utils_get_shape(ncid, "current_ep", 1, varid)
  ALLOCATE(current_ep(shap(1)))
  CALL nf90(nf90_get_var(ncid, varid, current_ep))
  shap(1:1) = nf90_utils_get_shape(ncid, "current_domain_id", 1, varid)
  ALLOCATE(current_domain_id(shap(1)))
  CALL nf90(nf90_get_var(ncid, varid, current_domain_id))

  shap = nf90_utils_get_shape(ncid, "current_datetime", 2, varid)
  ALLOCATE(current_datetime(shap(2)))
  CALL nf90(nf90_get_var(ncid, varid, current_datetime))

  DO i=1,SIZE(current_ep)
!    WRITE(0,*) "time is " // TRIM(current_datetime(i)) // &
!         & " calling EP " // int2string(current_ep(i)) // &
    !         & " on domain " // int2string(current_domain_id(i))
    CALL comin_current_set_datetime(trim_null(current_datetime(i)))
    CALL yac_routines(current_ep(i))
    CALL comin_callback_context_call(current_ep(i), current_domain_id(i), .FALSE.)
  END DO

  CALL nf90(nf90_close(ncid))

  CALL MPI_FINALIZE(ierr); CALL handle_mpi_errcode(ierr)

CONTAINS

  SUBROUTINE start_mpi()
    INTEGER :: ierr
#ifdef ENABLE_YAC
    CHARACTER(LEN=256), PARAMETER :: group_names(3) = [&
                                     "replay", &
                                     "comin ", &
                                     "yac   "]
#else
    CHARACTER(LEN=256), PARAMETER :: group_names(2) = [&
                                     "replay", &
                                     "comin "]
#endif
    INTEGER :: group_comms(3)

    CALL MPI_INIT (ierr); CALL handle_mpi_errcode(ierr)
    CALL mpi_handshake(MPI_COMM_WORLD, group_names, group_comms)

#ifdef ENABLE_YAC
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_finit_comm(group_comms(3), yac_instance_id)
#endif

    host_comm = group_comms(1)
    comin_comm = group_comms(2)
    CALL MPI_COMM_RANK(host_comm, host_rank, ierr); CALL handle_mpi_errcode(ierr)

    CALL MPI_COMM_SIZE(host_comm, host_size, ierr); CALL handle_mpi_errcode(ierr)
    IF (host_rank == 0) THEN
      WRITE (0,*) "running with ", host_size, " MPI tasks"
    END IF
  END SUBROUTINE start_mpi

  !> Utility function.
  SUBROUTINE handle_mpi_errcode(errcode)
    INTEGER, INTENT(IN) :: errcode
    IF (errcode .NE. MPI_SUCCESS) THEN
      CALL finish("replay", "Error in MPI program. Terminating.")
    END IF
  END SUBROUTINE handle_mpi_errcode

  SUBROUTINE finish(routine, text)
    CHARACTER(LEN=*), INTENT(IN) :: routine
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text
    INTEGER :: ierr, iexit

    CALL comin_callback_context_call(EP_FINISH, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)

    iexit = -1
    IF (host_rank == 0)  WRITE (0,*) routine, ": ", text, iexit
    CALL MPI_ABORT(MPI_COMM_WORLD, iexit, ierr)
  END SUBROUTINE finish

  SUBROUTINE add_variables(global, domain)
    TYPE(t_comin_descrdata_global), INTENT(IN) :: global
    TYPE(t_comin_descrdata_domain), INTENT(IN) :: domain(:)
    TYPE(t_var_request_list_item), POINTER :: ptr
    TYPE(t_comin_var_descriptor) :: descr
    TYPE(t_comin_var_metadata_iterator) :: metadata_it
    INTEGER :: no_vars = 0, no_tracers = 0, pos_jcjkjb(3), dimshape(5), dom_idx, ncontained
    INTEGER :: tracer_counter(global%n_dom)
    TYPE(t_comin_var_ptr), POINTER :: comin_var_ptr_ptr
    LOGICAL :: tracer, logbuff
    INTEGER :: zaxis_id, intbuff
    REAL(C_DOUBLE) :: realbuff
    CHARACTER(LEN=:),ALLOCATABLE :: current_metadata_key, charbuff

    ptr => comin_request_get_list_head()
    DO WHILE (ASSOCIATED(ptr))
      CALL ptr%item_value%metadata%get("tracer", tracer)
      IF (tracer) THEN
        no_tracers = no_tracers+1
      END IF
      no_vars = no_vars+1
      ptr => ptr%next()
    END DO

    ALLOCATE(vars(no_vars + global%n_dom))  ! plus 1 for the "tracer" container variable (on each domain)

    no_vars = 1
    ! add the "tracer" variables
    DO jg=1,global%n_dom
      ALLOCATE(vars(jg)%data_ptr(global%nproma, domain(jg)%nlev, &
           & domain(jg)%cells%nblks, 1, no_tracers))
      descr%id=jg
      descr%name="tracer"
      vars(no_vars)%comin_var_ptr = t_comin_var_ptr(              &
           &                      descriptor=descr,               &
           &                      ptr=vars(no_vars)%data_ptr,     &
           &                      device_ptr=C_NULL_PTR,          & ! currently not supported
           &                      pos_jc = 1,                     &
           &                      pos_jk = 2,                     &
           &                      pos_jb = 3,                     &
           &                      pos_jn = 5,                     &
           &                      lcontainer = .TRUE.,            &
           &                      ncontained=no_tracers)
      no_vars = no_vars + 1
    END DO

    tracer_counter = 0

    ptr => comin_request_get_list_head()
    DO WHILE (ASSOCIATED(ptr))
      dom_idx = ptr%item_value%descriptor%id
      CALL comin_metadata_get_or(ptr%item_value%metadata, "zaxis_id", zaxis_id, COMIN_ZAXIS_3D)
      IF (zaxis_id == COMIN_ZAXIS_3D) THEN
        dimshape=[global%nproma, domain(dom_idx)%nlev, domain(dom_idx)%cells%nblks,1,1]
        pos_jcjkjb =  [1, 2, 3]
      ELSE IF(zaxis_id == COMIN_ZAXIS_3D_HALF) THEN
        dimshape=[global%nproma, domain(dom_idx)%nlev + 1, domain(dom_idx)%cells%nblks,1,1]
        pos_jcjkjb =  [1, 2, 3]
      ELSE IF(zaxis_id == COMIN_ZAXIS_2D) THEN
        dimshape=[global%nproma,domain(dom_idx)%cells%nblks,1,1,1]
        pos_jcjkjb =  [1, -1, 2]
      ELSE
        CALL finish("comin_replay", "Unknown ZAXIS")
      END IF
      CALL ptr%item_value%metadata%get("tracer", tracer)
      IF (tracer) THEN
        tracer_counter(dom_idx) = tracer_counter(dom_idx) + 1
        ncontained = tracer_counter(dom_idx)
        vars(no_vars)%data_ptr => vars(dom_idx)%data_ptr(:,:,:,:,ncontained:ncontained)
      ELSE
        ALLOCATE(vars(no_vars)%data_ptr(dimshape(1), dimshape(2), &
             & dimshape(3), dimshape(4), dimshape(5)))
        ncontained = 0
      END IF
      vars(no_vars)%comin_var_ptr = t_comin_var_ptr(              &
           &                      descriptor=ptr%item_value%descriptor, &
           &                      ptr=vars(no_vars)%data_ptr,     &
           &                      device_ptr=C_NULL_PTR,          & ! currently not supported
           &                      pos_jc = pos_jcjkjb(1),         &
           &                      pos_jk = pos_jcjkjb(2),         &
           &                      pos_jb = pos_jcjkjb(3),         &
           &                      pos_jn = -1,                    &
           &                      lcontainer = .FALSE., &
           &                      ncontained=ncontained)
      comin_var_ptr_ptr => vars(no_vars)%comin_var_ptr
      CALL comin_var_list_append( &
           &       p          = comin_var_ptr_ptr)

      ! Iterate through request list metadata, forward metadata to var list
      CALL ptr%item_value%metadata%get_iterator(metadata_it)
      DO WHILE(.NOT. metadata_it%is_end())
        current_metadata_key = metadata_it%key()
        SELECT CASE(ptr%item_value%metadata%query(current_metadata_key))
        CASE (COMIN_METADATA_TYPEID_INTEGER)
          CALL ptr%item_value%metadata%get(TRIM(current_metadata_key), intbuff)
          CALL comin_metadata_set(ptr%item_value%descriptor, current_metadata_key, intbuff)
        CASE (COMIN_METADATA_TYPEID_REAL)
          CALL ptr%item_value%metadata%get(current_metadata_key, realbuff)
          CALL comin_metadata_set(ptr%item_value%descriptor, current_metadata_key, realbuff)
        CASE (COMIN_METADATA_TYPEID_CHARACTER)
          CALL ptr%item_value%metadata%get(current_metadata_key, charbuff)
          CALL comin_metadata_set(ptr%item_value%descriptor, current_metadata_key, charbuff)
        CASE (COMIN_METADATA_TYPEID_LOGICAL)
          CALL ptr%item_value%metadata%get(current_metadata_key, logbuff)
          CALL comin_metadata_set(ptr%item_value%descriptor, current_metadata_key, logbuff)
        END SELECT
        CALL metadata_it%next()
      ENDDO
      CALL metadata_it%delete()

      no_vars = no_vars+1
      ptr => ptr%next()
    END DO
  END SUBROUTINE add_variables

  FUNCTION trim_null( string)
    CHARACTER(len=*) :: string
    CHARACTER(len=LEN(string)) :: trim_null

    INTEGER :: pos

    pos = INDEX( string, ACHAR(0) )
    IF ( pos > 0 ) THEN
      trim_null = string(1:pos-1)
    ELSE
      trim_null = string
    ENDIF
  END FUNCTION trim_null

  SUBROUTINE yac_routines(current_ep)
    INTEGER, INTENT(IN) :: current_ep
#ifdef ENABLE_YAC
    IF (current_ep == EP_ATM_YAC_DEFCOMP_AFTER) THEN
      CALL yac_fdef_comp(yac_instance_id, "comin_replay", yac_comp_id)
    ENDIF
    IF (current_ep == EP_ATM_YAC_SYNCDEF_AFTER) THEN
      CALL yac_fsync_def(yac_instance_id)
    ENDIF
    IF (current_ep == EP_ATM_YAC_ENDDEF_AFTER) THEN
      CALL yac_fenddef(yac_instance_id)
    ENDIF
#endif
  END SUBROUTINE yac_routines

END PROGRAM comin_replay
