!> ICON model (program).
!
!  -----------------------------------------------------------------------------
!  This is a mockup (kind of a "scale model") of the ICON source code,
!  with the purpose of demonstrating and testing the interfaces of
!  ICON ComIn (the ICON Community interface).
!  -----------------------------------------------------------------------------
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
PROGRAM minimal_example
  USE mpi,                     ONLY: MPI_COMM_RANK, MPI_COMM_SIZE, MPI_BARRIER, MPI_INIT,     &
    &                                MPI_FINALIZE, MPI_SUCCESS, MPI_COMM_WORLD
  USE OMP_LIB
  USE mo_utilities,            ONLY: wp, finish, mpi_rank, int2string
  USE vars,                    ONLY: t_state, t_var, t_var_item, add_var,                     &
    &                                add_ref, get_var_timelevel, get_var_name,                &
    &                                find_list_element, deallocate_vars
  USE descr_data,         ONLY: setup_descr_data, p_patch, nproma,                  &
    &                                update_exposed_descriptive_data
  USE vars,               ONLY: expose_variables, append_comin_variables,      &
    &                                update_exposed_variables, get_ntracer_comin
  USE descr_data,         ONLY: clear_descr_data, expose_descriptive_data, max_dom,   &
    &                           n_dom, min_rlcell_int, min_rlcell, max_rlcell,        &
    &                           min_rlvert_int, min_rlvert, max_rlvert, &
    &                           min_rledge_int, min_rledge, max_rledge, &
    &                           grf_bdywidth_c,                         &
                                grf_bdywidth_e, dt, exp_start, exp_stop, sim_current,     &
                                run_start, run_stop
  USE comin_host_interface,    ONLY: comin_setup_init,                                        &
    &                                comin_setup_check, comin_callback_context_call,          &
    &                                comin_plugin_primaryconstructor,                         &
    &                                EP_SECONDARY_CONSTRUCTOR, EP_DESTRUCTOR,                 &
    &                                EP_ATM_WRITE_OUTPUT_BEFORE, EP_SECONDARY_CONSTRUCTOR,    &
    &                                EP_ATM_ADVECTION_BEFORE, EP_ATM_PHYSICS_BEFORE,          &
    &                                EP_ATM_YAC_DEFCOMP_BEFORE, EP_ATM_YAC_DEFCOMP_AFTER,     &
    &                                EP_ATM_YAC_SYNCDEF_BEFORE, EP_ATM_YAC_SYNCDEF_AFTER,     &
    &                                EP_ATM_YAC_ENDDEF_BEFORE, EP_ATM_YAC_ENDDEF_AFTER,       &
    &                                EP_ATM_TIMELOOP_START, EP_ATM_TIMELOOP_END,              &
    &                                comin_var_list_finalize, comin_var_request_list_finalize,&
    &                                comin_descrdata_finalize,                                &
    &                                comin_setup_finalize, t_comin_plugin_description,        &
    &                                t_comin_setup_version_info, comin_setup_get_version,     &
    &                                comin_setup_errhandler,                                  &
    &                                comin_parallel_mpi_handshake, mpi_handshake,             &
    &                                COMIN_DOMAIN_OUTSIDE_LOOP, comin_setup_set_verbosity_level

#ifdef ENABLE_YAC
  USE mo_yac_finterface, ONLY: yac_finit_comm, yac_fdef_comp, yac_ffinalize, yac_fdef_datetime, &
    &                          YAC_PROLEPTIC_GREGORIAN, yac_fsync_def, yac_fenddef
#endif

  IMPLICIT NONE

  INTEGER, PARAMETER :: nsteps   =  3      ! number of timesteps
  LOGICAL, PARAMETER :: lrestart = .FALSE. !< flag: simulation restart?
  INTEGER            :: jstep              ! time step counter
  INTEGER            :: jg                 ! counter for domains

  ! ICON variables: data vector with 2 timelevels
  TYPE(t_state) :: state_vector(2)
  INTEGER       :: nnow, nold, tl, mpi_size, nblks_c, nblks_e, nlev, i, ierr, &
    &              comm, nplugins, ntracer, ntracer_comin, comin_comm
  TYPE(t_var_item), POINTER :: vn
  TYPE(t_comin_setup_version_info) :: comin_version
  TYPE(t_comin_plugin_description) :: plugin_list(16) !< list of dynamic libs (max: 16)
  INTEGER :: yac_instance_id = 0
#ifdef ENABLE_YAC
  INTEGER yac_comp_id
#endif

  NAMELIST /comin_nml/ plugin_list

  CALL start_mpi()

  CALL comin_setup_init(mpi_rank == 0)

  ! -------------------------------------------------------------------------
  ! ICON ComIn primary constructor
  ! -------------------------------------------------------------------------
  !
  ! - compatibility checks (library versions)
  ! - requests for additional model variables
  ! - registration of callback functions
  !
  comin_version = comin_setup_get_version()
  IF (mpi_rank == 0)  WRITE (0,*) "     compatibility check."
  IF (mpi_rank == 0)  WRITE (0,'(3(a,i0))') "        linked to ICON Community Interface v", &
       &                        comin_version%version_no_major, ".", comin_version%version_no_minor, &
       &                        ".", comin_version%version_no_patch

  CALL comin_setup_errhandler(finish)
  CALL comin_setup_set_verbosity_level(iverbosity=13)
  CALL comin_setup_check("minimal-example", wp)

  IF (mpi_rank == 0)  WRITE (0,*) "   ICON model setup"

  IF (mpi_rank == 0)  WRITE (0,*) "     ICON // Fill descriptive data structure"

  CALL MPI_COMM_SIZE(comm, mpi_size, ierr); CALL handle_mpi_errcode(ierr)
  CALL setup_descr_data(comm)

  nblks_c = p_patch(1)%nblks_c
  nblks_e = p_patch(1)%nblks_e
  nlev = p_patch(1)%nlev

  ! -------------------------------------------------------------------------
  ! the following code simulates a sequence of ICON "add_var's".
  ! -------------------------------------------------------------------------
  !
  IF (mpi_rank == 0)  WRITE (0,*) "     ICON // Fill variable structure"
  DO tl = 1,2 ! time level

    CALL add_var("u", tl, lcontainer=.FALSE., vc=vn, dimshape=[nproma,nlev,nblks_c,1,1], units="m/s^1")
    state_vector(tl)%u => vn%this%dataarray
    state_vector(tl)%u = REAL(1,wp) *   10._wp
    IF (mpi_rank == 0)  WRITE (0,'(a,a,a,f12.2)') "       model variable '", TRIM(vn%this%name), &
      &        "', filled with value ", vn%this%dataarray(1,1,1,1,1)

    CALL add_var("v", tl, lcontainer=.FALSE., vc=vn, dimshape=[nproma,nlev,nblks_c,1,1], units="m/s^1")
    state_vector(tl)%v => vn%this%dataarray
    state_vector(tl)%v = REAL(2,wp) *   10._wp
    IF (mpi_rank == 0)  WRITE (0,'(a,a,a,f12.2)') "       model variable '", TRIM(vn%this%name), &
      &        "', filled with value ", vn%this%dataarray(1,1,1,1,1)

    CALL add_var("pres", tl, lcontainer=.FALSE., vc=vn, dimshape=[nproma,nlev,nblks_c,1,1], units="Pa")
    state_vector(tl)%pres => vn%this%dataarray
    state_vector(tl)%pres = REAL(3,wp) *   10._wp
    IF (mpi_rank == 0)  WRITE (0,'(a,a,a,f12.2)') "       model variable '", TRIM(vn%this%name), &
      &        "', filled with value ", vn%this%dataarray(1,1,1,1,1)

    CALL add_var("temp", tl, lcontainer=.FALSE., vc=vn, dimshape=[nproma,nlev,nblks_c,1,1], units="K")
    state_vector(tl)%temp => vn%this%dataarray
    state_vector(tl)%temp = REAL(4,wp) *   10._wp
    IF (mpi_rank == 0)  WRITE (0,'(a,a,a,f12.2)') "       model variable '", TRIM(vn%this%name), &
      &        "', filled with value ", vn%this%dataarray(1,1,1,1,1)

    CALL add_var("rho", tl, lcontainer=.FALSE., vc=vn, dimshape=[nproma,nlev,nblks_c,1,1], units="kg/m3")
    state_vector(tl)%rho => vn%this%dataarray
    state_vector(tl)%rho = REAL(5,wp) *   10e-5_wp
    IF (mpi_rank == 0)  WRITE (0,'(a,a,a,f12.2)') "       model variable '", TRIM(vn%this%name), &
         &        "', filled with value ", vn%this%dataarray(1,1,1,1,1)

    CALL add_var("pres_sfc", tl, lcontainer=.FALSE., vc=vn, dimshape=[nproma,nblks_c,1,1,1], units="Pa")
    state_vector(tl)%pres_sfc => vn%this%dataarray
    state_vector(tl)%pres_sfc = REAL(5,wp) *   10._wp
    IF (mpi_rank == 0)  WRITE (0,'(a,a,a,f12.2)') "       model variable '", TRIM(vn%this%name), &
         &        "', filled with value ", vn%this%dataarray(1,1,1,1,1)

    CALL add_var("vn", tl, lcontainer=.FALSE., vc=vn, dimshape=[nproma,nlev,nblks_e,1,1], units='m s-1')
    state_vector(tl)%vn => vn%this%dataarray
    state_vector(tl)%vn = REAL(5,wp) *   10._wp
    IF (mpi_rank == 0)  WRITE (0,'(a,a,a,f12.2)') "       model variable '", TRIM(vn%this%name), &
      &        "', filled with value ", vn%this%dataarray(1,1,1,1,1)

  END DO

  IF (mpi_rank == 0)  WRITE (0,*) "     Associate ICON descriptive data structures and ComIn"
  CALL expose_descriptive_data(n_dom, max_dom, nproma,  &
       &                       min_rlcell_int, min_rlcell, max_rlcell, min_rlvert_int, &
       &                       min_rlvert, max_rlvert, min_rledge_int, min_rledge, max_rledge, &
       &                       grf_bdywidth_c, grf_bdywidth_e, &
       &                       dt, lrestart, exp_start, exp_stop, sim_current, run_start, run_stop, &
       &                       yac_instance_id)

  ! -------------------------------------------------------------------------
  ! ICON ComIn - call to primary constructor
  ! -------------------------------------------------------------------------
  !
  ! - third party modules obtain patch info and request additional variables.
  !
  IF (mpi_rank == 0)  WRITE (0,*) "     call of ComIn primary constructors."

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

  ! loop over the list of requested variables and count the number of
  ! additional tracer variables:
  ntracer_comin = get_ntracer_comin()

  IF (mpi_rank == 0)  WRITE (0,*) "     Add tracer variable (container)."
  ntracer = 6 + ntracer_comin
  CALL add_var("tracer", tl=1, lcontainer=.TRUE., dimshape=[nproma,nlev,nblks_c,ntracer,1], units="")
  CALL add_ref("tracer", "qr", tl=1, is_3d_var=.TRUE., vc=vn, units="kg/kg")
  vn%this%dataarray(:,:,:,vn%this%metadata%ncontained,1) = 1.0
  CALL add_ref("tracer", "qi", tl=1, is_3d_var=.TRUE., vc=vn, units="kg/kg")
  vn%this%dataarray(:,:,:,vn%this%metadata%ncontained,1) = 2.0
  CALL add_ref("tracer", "qv", tl=1, is_3d_var=.TRUE., vc=vn, units="kg/kg")
  vn%this%dataarray(:,:,:,vn%this%metadata%ncontained,1) = 3.0
  CALL add_ref("tracer", "qs", tl=1, is_3d_var=.TRUE., vc=vn, units="kg/kg")
  vn%this%dataarray(:,:,:,vn%this%metadata%ncontained,1) = 4.0
  CALL add_ref("tracer", "qc", tl=1, is_3d_var=.TRUE., vc=vn, units="kg/kg")
  vn%this%dataarray(:,:,:,vn%this%metadata%ncontained,1) = 5.0
  CALL add_ref("tracer", "qg", tl=1, is_3d_var=.TRUE., vc=vn, units="kg/kg")
  vn%this%dataarray(:,:,:,vn%this%metadata%ncontained,1) = 5.0

  CALL append_comin_variables(nproma, nlev, nblks_c)
  CALL expose_variables(nlev)

  ! -------------------------------------------------------------------------
  ! ICON ComIn - call to secondary constructor
  ! -------------------------------------------------------------------------
  !
  ! - third party modules retrieve pointers to data arrays, telling
  !   ICON ComIn about the context where these will be accessed.
  !
  CALL comin_callback_context_call(EP_SECONDARY_CONSTRUCTOR, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)

#ifdef ENABLE_YAC
  CALL yac_fdef_datetime(yac_instance_id, run_start, run_stop)
#endif
  CALL comin_callback_context_call(EP_ATM_YAC_DEFCOMP_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
#ifdef ENABLE_YAC
  CALL yac_fdef_comp(yac_instance_id, "icon", yac_comp_id)
#endif
  CALL comin_callback_context_call(EP_ATM_YAC_DEFCOMP_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
  CALL comin_callback_context_call(EP_ATM_YAC_SYNCDEF_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
#ifdef ENABLE_YAC
  CALL yac_fsync_def(yac_instance_id)
#endif
  CALL comin_callback_context_call(EP_ATM_YAC_SYNCDEF_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
  CALL comin_callback_context_call(EP_ATM_YAC_ENDDEF_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
#ifdef ENABLE_YAC
  CALL yac_fenddef(yac_instance_id)
#endif
  CALL comin_callback_context_call(EP_ATM_YAC_ENDDEF_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)

  ! -------------------------------------------------------------------------
  IF (mpi_rank == 0)  WRITE (0,"(/,a)") "    ICON model time loop"
  ! -------------------------------------------------------------------------

  nnow  = 1
  nold  = 2
  jg    = 1
  DO jstep=1,nsteps
    CALL comin_callback_context_call(EP_ATM_TIMELOOP_START, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
    IF (mpi_rank == 0)  WRITE (0,'(/,a,i0)') "      time step ", jstep
    IF (mpi_rank == 0)  WRITE (0,'(a,/)')    "      ==========="

    !> example of a third-party entry point (callback)
    CALL comin_callback_context_call(EP_ATM_ADVECTION_BEFORE, jg, .FALSE.)

    ! ------------------------------------------------
    ! DUMMY CODE: "ADVECTION"
    !
    IF (mpi_rank == 0)  WRITE (0,'(/,a)') "        <advection>"

    state_vector(nnow)%u    = state_vector(nnow)%u + 2.0_wp
    state_vector(nnow)%v    = state_vector(nnow)%v + 2.0_wp
    state_vector(nnow)%pres = state_vector(nnow)%pres + 2.0_wp
    state_vector(nnow)%temp = state_vector(nnow)%temp + 2.0_wp

    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        u: value    = ", state_vector(nnow)%u(1,1,1,1,1)
    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        v: value    = ", state_vector(nnow)%v(1,1,1,1,1)
    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        pres: value = ", state_vector(nnow)%pres(1,1,1,1,1)
    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        temp: value = ", state_vector(nnow)%temp(1,1,1,1,1)

    IF (mpi_rank == 0)  WRITE (0,*) " "
    ! ------------------------------------------------

    !> example of a third-party entry point (callback)
    CALL comin_callback_context_call(EP_ATM_WRITE_OUTPUT_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)

    ! ------------------------------------------------
    ! DUMMY CODE: "OUTPUT"
    IF (mpi_rank == 0)  WRITE (0,'(/,a,/)') "        <output>"
    !
    ! do nothing
    ! ------------------------------------------------

    !> example of a third-party entry point (callback)
    CALL comin_callback_context_call(EP_ATM_PHYSICS_BEFORE, jg, .FALSE.)

    ! ------------------------------------------------
    ! DUMMY CODE: "PHYS"
    !
    IF (mpi_rank == 0)  WRITE (0,'(/,a)') "        <phys>"

    state_vector(nnow)%u    = state_vector(nnow)%u + 3.0_wp
    state_vector(nnow)%v    = state_vector(nnow)%v + 3.0_wp
    state_vector(nnow)%pres = state_vector(nnow)%pres + 3.0_wp
    state_vector(nnow)%temp = state_vector(nnow)%temp + 3.0_wp

    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        u: value    = ", state_vector(nnow)%u(1,1,1,1,1)
    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        v: value    = ", state_vector(nnow)%v(1,1,1,1,1)
    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        pres: value = ", state_vector(nnow)%pres(1,1,1,1,1)
    IF (mpi_rank == 0)  WRITE (0,'(a,F9.4)') "        temp: value = ", state_vector(nnow)%temp(1,1,1,1,1)

    IF (mpi_rank == 0)  WRITE (0,*) " "
    ! ------------------------------------------------

    ! Switch the time levels in ICON
    nnow = 1 + MOD(jstep,     2)
    nold = 1 + MOD(jstep + 1, 2)
    IF (mpi_rank == 0)  WRITE (0,*) "     nnow, nold = ", nnow, nold

    ! Update all pointers exposed to the ComIn
    sim_current = "1950-01-01T00:"//int2string(jstep,'(i0)')//"0:00"
    CALL update_exposed_descriptive_data(sim_current)
    CALL update_exposed_variables(nnow)

    CALL comin_callback_context_call(EP_ATM_TIMELOOP_END, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
  END DO

  !> example of a third-party entry point (callback)
  !> at the end of the timeloop to initialize clean-up in 3rd party modules
  CALL comin_callback_context_call(EP_DESTRUCTOR, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)

  ! clean-up
  IF (mpi_rank == 0)  WRITE (0,"(/,a)") "    ICON clean-up"

  CALL MPI_BARRIER(comm, ierr); CALL handle_mpi_errcode(ierr)

  CALL comin_var_list_finalize()
  CALL comin_var_request_list_finalize()
  CALL comin_descrdata_finalize()
  CALL comin_setup_finalize()

  CALL clear_descr_data()

  CALL deallocate_vars()
#ifdef ENABLE_YAC
  CALL yac_ffinalize()
#endif
  CALL MPI_FINALIZE(ierr); CALL handle_mpi_errcode(ierr)

CONTAINS

  !> Initialize MPI, this should always be the first call
  SUBROUTINE start_mpi()
    INTEGER :: ierr, p_comm_size
    CHARACTER(LEN=256), PARAMETER :: group_names(3) = ["minimal-example", "yac            ", "comin          "]
    INTEGER :: group_comms(3)

    CALL MPI_INIT (ierr); CALL handle_mpi_errcode(ierr)
    CALL mpi_handshake(MPI_COMM_WORLD, group_names, group_comms)

#ifdef ENABLE_YAC
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_finit_comm(group_comms(2), yac_instance_id)
#endif

    comm = group_comms(1)
    comin_comm = group_comms(3)
    CALL MPI_COMM_RANK(comm, mpi_rank, ierr); CALL handle_mpi_errcode(ierr)

    IF (mpi_rank == 0)  WRITE (0,*) "ICON running"

    CALL MPI_COMM_SIZE(comm, p_comm_size, ierr); CALL handle_mpi_errcode(ierr)
    IF (mpi_rank == 0) THEN
      WRITE (0,*) "   running with ", p_comm_size, &
        &         " MPI tasks and ", OMP_GET_NUM_THREADS(), " thread(s) per task."
    END IF
  END SUBROUTINE start_mpi

  !> Utility function.
  SUBROUTINE handle_mpi_errcode(errcode)
    INTEGER, INTENT(IN) :: errcode
    INTEGER :: ierr
    IF (errcode .NE. MPI_SUCCESS) THEN
      WRITE (0,*) "Error in MPI program. Terminating."
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    END IF
  END SUBROUTINE handle_mpi_errcode

END PROGRAM minimal_example
