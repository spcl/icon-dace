! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Scheduler for internal post-processing.
! ===================================================================
!
! This module manages a "job queue" for internal post-processing
! tasks on the compute PEs.
!
! For example, the interpolation of model variables onto lon-lat
! fields constitutes such a task. Allocating and computing only those
! fields which are required for output (or as intermediate results)
! can save memory and computing time.
!
! Jobs are processed according to user-defined priority levels. In
! principle, all tasks with the same job priority could be processed
! simultaneously (with OpenMP).
!
! List of post-processing tasks implemented so far:
!
! Type of task                        Execution priority
! -------------------------------------------------------------------
!
! z/p/i interpolation setup           HIGH_PRIORITY
! compute vertical velocity           DEFAULT_PRIORITY0
! compute rel. humidity               DEFAULT_PRIORITY0
! compute PV                          DEFAULT_PRIORITY0
! compute mean sea level pressure     DEFAULT_PRIORITY0
! vertical interpolation              DEFAULT_PRIORITY1
! horizontal interpolation edge->cell DEFAULT_PRIORITY2
! horizontal synchronization          DEFAULT_PRIORITY3
! horizontal interpolation            DEFAULT_PRIORITY4
! z/p/i interpolation clean-up        LOW_PRIORITY
!
! ===================================================================
!
! DETAILS: Vertical interpolation of output variables
! ---------------------------------------------------
!
! Vertical interpolation is one of the internal post-processing tasks
! in ICON (see the table above for other post-processing tasks).
! Vertical interpolation of a variable is enabled/disabled simply by
! adding the variable's name to the corresponding namelist parameter
! ("pl_varlist", "hl_varlist", "il_varlist").
!
! Vertical interpolation of a field is only possible if the necessary
! meta-data, e.g. the interpolation method, has been defined in the
! "add_var" call for this variable. A typical example would be
!
!    CALL add_var( p_prog_list, 'rho', p_prog%rho,  &
!      ...
!      vert_interp=create_vert_interp_metadata(                     &
!                    vert_intp_type=vintp_types("P","Z","I"),       &
!                    vert_intp_method=VINTP_METHOD_LIN ) )
!
! In this example, vertical interpolation is enabled for p-, z- and
! i-level interpolation, using a linear interpolation method.
!
! In general, available settings are
!
!  vert_intp_type      : vertical interpolation type, one or more of
!                        mo_var_metadata_types::VINTP_TYPE_LIST.
!                        Default: no vertical interpolation
!
!  vert_intp_method    : vertical interpolation algorithms, listed in mo_impl_constants, and
!                        defined in module mo_nh_vert_interp:
!
!                         VINTP_METHOD_LIN (default) : linear vertical interpolation
!                         VINTP_METHOD_QV            : vertical interpolation of specific humidity,
!                                                      performs cubic interpolation where possible,
!                                                      turning to linear interpolation close to the surface
!                         VINTP_METHOD_PRES          : vertical interpolation of pressure, piecewise
!                                                      analytical integration of the hydrostatic equation
!                         VINTP_METHOD_LIN_NLEVP1    : linear interpolation for half level fields
!                         VINTP_METHOD_VN            : vertical interpolation and extrapolation of horizontal
!                                                      wind, performs cubic interpolation where
!                                                      possible, turning to linear interpolation close to the
!                                                      surface with boundary-layer treatment
!                                                      - Please note than wind fields are treated in a special way,
!                                                        see below. -
!
!  Special treatment of wind fields:
!
!  Wind fields are vertically interpolated as follows: A new
!  z/p/i-variable "vn" is created and a post-processing task for
!  vertical interpolation of the model's "vn" onto this new
!  field. Then, new cell-based variables "u", "v" on the same
!  vertical axis are created and a post-processing task for edge2cell
!  interpolation "vn" -> "u","v". Thus, when requesting "U", "V", we
!  actually get "VN" vertically interpolated.
!  The vertical interpolation method for the wind fields is therefore
!  specified by the "add_var(...)" for the normal velocity component
!  "VN".

!
!  Tuning parameters for the interpolation algorithms
!  (defaults are defined in mo_var_metadata::create_vert_interp_metadata)
!
!     l_hires_intp             : mode for interpolation to (much) finer grid (VINTP_METHOD_VN)
!                                Default: .FALSE.
!     l_restore_fricred        : subtract/restore frictional reduction of wind speed (VINTP_METHOD_VN)
!                                Default: .FALSE.
!     l_loglin                 : setting l_loglin=.TRUE. activates logarithmic interpolation
!                                (only for VINTP_METHOD_LIN)
!                                Default: .FALSE.
!     l_satlimit               : limit input field to water saturation (VINTP_METHOD_QV)
!                                Default: .FALSE.
!     l_restore_pbldev         : restore PBL deviation of QV from extrapolated profile (VINTP_METHOD_QV)
!                                Default: .FALSE.
!
!     l_pd_limit               : Switch for use of positive definite limiter (VINTP_METHOD_LIN)
!                                Default: .FALSE.
!     lower_limit              : Limiter value to avoid negative or unreasonably small values
!                                (VINTP_METHOD_LIN, VINTP_METHOD_QV). For the linear interpolation method,
!                                the "lower_limit" is used only in combination with "l_pd_limit". Default is 0.
!
!     l_extrapol               : Switch for use of downward extrapolation (only for VINTP_METHOD_LIN)
!                                where the gradient between height "zpbl1" and "zpbl2" is used
!                                (see mo_initicon_config). If "l_extrapol==.FALSE.", then for all levels
!                                below the lowermost input level we use the values from the lowermost input
!                                level.
!                                Default: .TRUE.
!                                Note: Logarithmic computation is not used for extrapolation because
!                                      it would be numerically unstable.
!
! ===================================================================
!
! TODO[FP] To increase performance, one should allocate/deallocate
!          the vertical interpolation coefficient tables "vcoeff"
!          only when necessary.
!
! TODO[FP] Interpolation tasks are performed more often than
!          necessary: The activity flag must be adjusted to the output
!          intervals!
! TODO[FP] Do not insert post-processing tasks for patches with 0 cells.
!
! -----------------------------------------------------------------------------------

MODULE mo_pp_scheduler

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_impl_constants,          ONLY: SUCCESS, HINTP_TYPE_NONE, max_var_ml,               &
    &                                   max_var_pl, max_var_hl, max_var_il, TASK_NONE,      &
    &                                   TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I,  &
    &                                   TASK_FINALIZE_IPZ, TASK_INTP_HOR_LONLAT,            &
    &                                   TASK_INTP_VER_PLEV, TASK_INTP_SYNC, TASK_INTP_MSL,  &
    &                                   TASK_COMPUTE_RH, TASK_COMPUTE_PV, TASK_COMPUTE_SMI, &
    &                                   TASK_COMPUTE_SDI2, TASK_COMPUTE_LPI,                &
    &                                   TASK_COMPUTE_HBAS_SC, TASK_COMPUTE_HTOP_SC,         &
    &                                   TASK_COMPUTE_TWATER, TASK_COMPUTE_Q_SEDIM,          &
    &                                   TASK_COMPUTE_DBZ850, TASK_COMPUTE_DBZCMAX,          &
    &                                   TASK_COMPUTE_DBZLMX_LOW,                            &
    &                                   TASK_COMPUTE_CEILING, TASK_COMPUTE_VIS,             &
    &                                   TASK_COMPUTE_WSHEAR_U, TASK_COMPUTE_WSHEAR_V,       &
    &                                   TASK_COMPUTE_LAPSERATE, TASK_COMPUTE_MCONV,         &
    &                                   TASK_COMPUTE_SRH,                                   &
    &                                   TASK_COMPUTE_INVERSION,                             &
    &                                   TASK_INTP_VER_ZLEV,                                 &
    &                                   TASK_INTP_VER_ILEV, TASK_INTP_EDGE2CELL,            &
    &                                   UNDEF_TIMELEVEL, ALL_TIMELEVELS,                    &
    &                                   vname_len, TASK_COMPUTE_OMEGA,                      &
    &                                   TLEV_NNOW, TLEV_NNOW_RCF, HINTP_TYPE_LONLAT_NNB,    &
    &                                   STR_HINTP_TYPE
  USE mo_cdi_constants,           ONLY: GRID_CELL, GRID_UNSTRUCTURED_CELL,                  &
    &                                   GRID_REGULAR_LONLAT, LONLAT_PREFIX
  USE mo_model_domain,            ONLY: p_patch, p_phys_patch
  USE mo_var_list,                ONLY: add_var, find_list_element, t_var_list_ptr
  USE mo_var_list_register,       ONLY: t_vl_register_iter
  USE mo_var, ONLY: t_var, level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_var_metadata_types,      ONLY: t_var_metadata, t_var_metadata_dynamic, t_post_op_meta, &
    &                                   t_hor_interp_meta
  USE mo_var_metadata,            ONLY: create_hor_interp_metadata, vintp_type_id,          &
    &                                   get_timelevel_string, get_var_timelevel, get_var_name
  USE mo_intp_data_strc,          ONLY: p_int_state
  USE mo_intp_lonlat_types,       ONLY: t_lon_lat_intp, lonlat_grids
  USE mo_nonhydro_state,          ONLY: p_nh_state, p_nh_state_lists
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, p_nh_opt_diag
  USE mo_nwp_phy_state,           ONLY: prm_diag
  USE mo_nh_pzlev_config,         ONLY: nh_pzlev_config
  USE mo_name_list_output_config, ONLY: first_output_name_list
  USE mo_name_list_output_types,  ONLY: t_output_name_list, is_grid_info_var, &
    &                                   remap_regular_latlon
  USE mo_parallel_config,         ONLY: nproma
  USE mo_io_config,               ONLY: lnetcdf_flt64_output
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var, grib2_var
  USE mo_util_string,             ONLY: int2string, remove_duplicates,                      &
    &                                   difference, tolower
  USE mo_cdi,                     ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16,    &
    &                                   GRID_UNSTRUCTURED,TSTEP_INSTANT, TSTEP_CONSTANT
  USE mo_zaxis_type,              ONLY: ZA_ALTITUDE, ZA_PRESSURE, ZA_ISENTROPIC, zaxisTypeList
  USE mo_pp_tasks,                ONLY: pp_task_lonlat, pp_task_sync, pp_task_ipzlev_setup, &
    &                                   pp_task_ipzlev, pp_task_compute_field,              &
    &                                   pp_task_intp_msl, pp_task_edge2cell,                & 
    &                                   t_simulation_status, t_job_queue, job_queue,        &
    &                                   HIGH_PRIORITY,                                      &
    &                                   DEFAULT_PRIORITY0, DEFAULT_PRIORITY1,               &
    &                                   DEFAULT_PRIORITY2, DEFAULT_PRIORITY3,               &
    &                                   DEFAULT_PRIORITY4, LOW_PRIORITY, dbg_level,         &
    &                                   t_activity_status
  USE mo_timer,                   ONLY: timers_level, timer_start, timer_stop, timer_opt_diag_atmo
  USE mo_util_texthash,           ONLY: text_hash_c
  USE mo_fortran_tools,           ONLY: assert_acc_device_only

  IMPLICIT NONE
  PRIVATE

  ! functions and subroutines
  PUBLIC :: pp_scheduler_init
  PUBLIC :: pp_scheduler_process
  PUBLIC :: pp_scheduler_finalize
  PUBLIC :: new_simulation_status

  !> module name string
  CHARACTER(*), PARAMETER :: modname = 'mo_pp_scheduler'

CONTAINS

  !--- SCHEDULER ---------------------------------------------------------------------

  !---------------------------------------------------------------
  !> Setup of post-processing job queue.
  !
  ! E.g. setup of optional diagnostic quantities like pz-level
  ! interpolation. Computation of these additional variables is
  ! registered as post-processing tasks.
  !
  ! @note This subroutine must be called after the namelists have been
  !       read and the nh_state has been constructed and _before_
  !       initialization of the name list output.
  !
  SUBROUTINE pp_scheduler_init(l_init_prm_diag)
    LOGICAL, INTENT(IN) :: l_init_prm_diag
    CHARACTER(*), PARAMETER :: routine =  modname//"::pp_scheduler_init"
    INTEGER :: jg, iv
    TYPE(t_var), POINTER :: elem, elem_pres
    TYPE(t_vl_register_iter) :: vl_iter

    if (dbg_level > 5)  CALL message(routine, "Enter")
    !-------------------------------------------------------------
    !--- setup of optional diagnostic fields updated by the 
    !    post-processing scheduler

    !- loop over model level variables
    DO WHILE(vl_iter%next())
      jg = vl_iter%cur%p%patch_id
      DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)%p
        IF (elem%info%l_pp_scheduler_task .EQ. TASK_NONE) CYCLE
        IF (dbg_level > 5)  &
          & CALL message(routine, "Inserting pp task: "//TRIM(elem%info%name))
        SELECT CASE(elem%info%l_pp_scheduler_task)
        CASE (TASK_COMPUTE_RH,      TASK_COMPUTE_OMEGA,   TASK_COMPUTE_PV,      &
          &   TASK_COMPUTE_WSHEAR_U, TASK_COMPUTE_WSHEAR_V, TASK_COMPUTE_SRH,   &
          &   TASK_COMPUTE_LPI,     TASK_COMPUTE_CEILING, TASK_COMPUTE_HBAS_SC, &
          &   TASK_COMPUTE_HTOP_SC, TASK_COMPUTE_TWATER,  TASK_COMPUTE_Q_SEDIM, &
          &   TASK_COMPUTE_DBZ850,  TASK_COMPUTE_DBZCMAX, TASK_COMPUTE_SMI,     &
          &   TASK_COMPUTE_DBZLMX_LOW, TASK_COMPUTE_SDI2, TASK_COMPUTE_VIS,     &
          &   TASK_COMPUTE_INVERSION, TASK_COMPUTE_LAPSERATE, TASK_COMPUTE_MCONV )
          CALL pp_scheduler_register(name=elem%info%name, jg=jg, p_out_var=elem, &
            &    l_init_prm_diag=l_init_prm_diag, job_type=elem%info%l_pp_scheduler_task ) 
        CASE (TASK_INTP_MSL)
          ! mean sea level pressure
          ! find the standard pressure field:
          elem_pres => find_list_element(p_nh_state_lists(jg)%diag_list, 'pres')
          IF (ASSOCIATED (elem_pres)) &
            ! register task for interpolation to z=0:
            & CALL pp_scheduler_register(name=elem%info%name, jg=jg, p_out_var=elem, &
              & job_type=TASK_INTP_MSL, l_init_prm_diag=l_init_prm_diag, opt_p_in_var=elem_pres)
        CASE DEFAULT
          CALL finish(routine, "Unknown pp task type!")
        END SELECT
      ENDDO ! loop over vlist "i"
    ENDDO ! i = 1, SIZE(var_lists)
    !-------------------------------------------------------------
    !--- setup of vertical interpolation onto i/p/z-levels
    CALL pp_scheduler_init_ipz(l_init_prm_diag)
    !-------------------------------------------------------------
    !-- horizontal interpolation regular lon-lat grids
    CALL pp_scheduler_init_lonlat()
    IF (dbg_level > 5)  CALL message(routine, "Done")
  END SUBROUTINE pp_scheduler_init

  !---------------------------------------------------------------
  !> (Internal) Utility routine, add a new "vn" field for a given
  !  axis, based on the meta-data of the standard "vn". - This is done
  !  for all available time levels.
  !
  SUBROUTINE init_vn_horizontal(ll_grid_id, lev_type)
    INTEGER,          INTENT(IN)  :: ll_grid_id      !< lon-lat grid number
    INTEGER,          INTENT(IN)  :: lev_type        !< level type: p/z/i/m
    ! local variables
    CHARACTER(*), PARAMETER :: routine =  modname//"::init_vn_horizontal"
    TYPE(t_var), POINTER :: elem_u, elem_v, elem, new_elem, new_elem2
    INTEGER :: iv, shape3d_ll(3), nblks_lonlat, nlev, jg, tl, vn_hash
    TYPE(t_job_queue),    POINTER :: task
    TYPE(t_var_metadata), POINTER :: info
    REAL(wp),             POINTER :: p_opt_field_r3d(:,:,:)
    CHARACTER(len=vname_len)    :: vname
    CHARACTER(len=4)              :: suffix
    TYPE(t_cf_var)                :: cf
    TYPE(t_grib2_var)             :: grib2
    TYPE(t_post_op_meta)          :: post_op
    TYPE(t_var_list_ptr), POINTER     :: dst_varlist     !< destination variable list
    TYPE (t_lon_lat_intp), POINTER:: ptr_int_lonlat
    CHARACTER(LEN=1)              :: prefix
    TYPE(t_vl_register_iter) :: vl_iter
    TYPE(t_hor_interp_meta)  :: hor_interp_meta

    IF (dbg_level > 5)  CALL message(routine, "Enter")
    
    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_elem => NULL()
    new_elem2 => NULL()
    vn_hash = text_hash_c('vn')
    !- loop over model level variables
    ! Note that there are several "vn" variables with different time
    ! levels, we just add unconditionally all
    DO WHILE(vl_iter%next())
      jg = vl_iter%cur%p%patch_id
      SELECT CASE(lev_type)
      CASE (level_type_ml)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list         
        prefix = "m"
      CASE (level_type_pl)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list_p
        prefix = "p"
      CASE (level_type_hl)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list_z
        prefix = "z"
      CASE (level_type_il)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list_i
        prefix = "i"
      END SELECT
      ! Do not inspect lists which are disabled for output
      IF (.NOT.vl_iter%cur%p%loutput) CYCLE
      ! loop only over model level variables
      IF (vl_iter%cur%p%vlevel_type /= lev_type) CYCLE         
      ! loop only over variables of where domain was requested
      IF (.NOT. lonlat_grids%list(ll_grid_id)%l_dom(jg)) CYCLE

      VAR_LOOP : DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)%p
        ! now, search for "vn" in the variable list:
        info => elem%info
        ! Do not inspect element if it is a container
        IF (info%lcontainer) CYCLE
        ! Do not inspect element if "loutput=.false."
        IF (.NOT. info%loutput) CYCLE
        ! Check for matching name
        IF (vn_hash .NE. vl_iter%cur%p%key_notl(iv)) CYCLE
        ! get time level
        tl = get_var_timelevel(info%name)
        suffix = ''
        IF (tl /= -1)  suffix = get_timelevel_string(tl)
        !- find existing variables "u", "v" (for copying the meta-data):
        elem_u => find_list_element(p_nh_state_lists(jg)%diag_list, "u")
        elem_v => find_list_element(p_nh_state_lists(jg)%diag_list, "v")
        !- predefined array shapes
        nlev = info%used_dimensions(2)
        ptr_int_lonlat => lonlat_grids%list(ll_grid_id)%intp(jg)
        nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
        shape3d_ll = (/ nproma, nlev, nblks_lonlat /)
        !-- create new cell-based variables "u", "v" on lon-lat grid
        !   for the same time level
        IF (dbg_level > 8) &
          CALL message(routine, "horizontal interpolation: create u/v variables on lon-lat grid")
        vname   = LONLAT_PREFIX//TRIM(get_var_name(elem_u%info))//suffix
        cf      = elem_u%info%cf
        grib2   = elem_u%info%grib2
        post_op = elem_u%info%post_op
        hor_interp_meta = elem_u%info%hor_interp
        hor_interp_meta%lonlat_id  = ll_grid_id ! link to the lon-lat grid
        CALL add_var(dst_varlist, TRIM(vname), p_opt_field_r3d, GRID_REGULAR_LONLAT, &
          & info%vgrid, cf, grib2, ldims=shape3d_ll, lrestart=.FALSE.,               &
          & in_group=elem_u%info%in_group, new_element=new_elem, loutput=.TRUE.,     &
          & post_op=post_op, var_class=elem_u%info%var_class, lopenacc=.TRUE.,       &
          & tlev_source=info%tlev_source, hor_interp=hor_interp_meta,                &
          & vert_interp=elem_u%info%vert_interp)
        vname   = LONLAT_PREFIX//TRIM(get_var_name(elem_v%info))//suffix
        cf      = elem_v%info%cf
        grib2   = elem_v%info%grib2
        post_op = elem_v%info%post_op
        hor_interp_meta = elem_v%info%hor_interp
        hor_interp_meta%lonlat_id = ll_grid_id ! link to the lon-lat grid
        CALL add_var(dst_varlist, TRIM(vname), p_opt_field_r3d, GRID_REGULAR_LONLAT, &
          & info%vgrid, cf, grib2, ldims=shape3d_ll, lrestart=.FALSE.,               &
          & in_group=elem_v%info%in_group, new_element=new_elem2, loutput=.TRUE.,    &
          & post_op=post_op, var_class=elem_v%info%var_class, lopenacc=.TRUE.,       &
          & tlev_source=info%tlev_source, hor_interp=hor_interp_meta,                &
          & vert_interp=elem_v%info%vert_interp)
        !-- create and add post-processing task
        task => pp_task_insert(DEFAULT_PRIORITY4)
        WRITE (task%job_name, *) "horizontal interp. ",TRIM(info%name),", ",prefix//"-levels", ", DOM ",jg
        IF (dbg_level > 8) CALL message(routine, task%job_name)
        task%data_input%p_nh_state          => NULL()
        task%data_input%prm_diag            => NULL()
        task%data_input%nh_pzlev_config     => NULL()
        task%data_input%jg                  =  jg           
        task%data_input%p_patch             => p_patch(jg)
        task%data_input%p_nh_opt_diag       => p_nh_opt_diag(jg)
        task%data_input%p_int_state         => p_int_state(jg)
        task%job_type                       =  TASK_INTP_HOR_LONLAT
        task%activity                       =  new_activity_status(l_output_step=.TRUE.)
        task%activity%check_dom_active      =  .TRUE.
        task%activity%i_timelevel           =  tl
        task%data_input%var                 => elem      ! set input variable
        task%data_output%var                => new_elem  ! set output variable "u"
        task%data_output%var_2              => new_elem2 ! set output variable "v"
      END DO VAR_LOOP
    END DO
    if (dbg_level > 5)  CALL message(routine, "Done")
  END SUBROUTINE init_vn_horizontal

  !---------------------------------------------------------------
  !> Setup of lon-lat interpolation tasks.
  !
  SUBROUTINE pp_scheduler_init_lonlat

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  modname//"::pp_scheduler_init_lonlat"
    INTEGER :: iv, jg, ndom, ierrstat, ivar, i, j, nvars_ll, nblks_lonlat, &
      & ilev_type, max_var, ilev, n_uv_hrz_intp, var_shape(5)
    LOGICAL                               :: l_horintp
    TYPE (t_output_name_list), POINTER    :: p_onl
    TYPE(t_job_queue),         POINTER    :: task
    TYPE(t_var_list_ptr),      POINTER    :: p_opt_diag_list
    REAL(wp), POINTER                     :: p_opt_field_r3d(:,:,:)
    INTEGER,  POINTER                     :: p_opt_field_i3d(:,:,:)
    TYPE(t_var), POINTER :: elem, new_elem
    CHARACTER(LEN=vname_len),  POINTER    :: varlist(:)
    INTEGER, ALLOCATABLE :: ll_vargrid(:), ll_varlevs(:)
    CHARACTER(LEN=vname_len), ALLOCATABLE :: ll_varlist(:)
    CHARACTER(LEN=vname_len)              :: vname
    TYPE(t_var_metadata),POINTER          :: info
    TYPE(t_var_metadata_dynamic),POINTER  :: info_dyn
    TYPE (t_lon_lat_intp),     POINTER    :: ptr_int_lonlat
    INTEGER :: uv_hrz_intp_grid(4*lonlat_grids%ngrids), &
         uv_hrz_intp_levs(4*lonlat_grids%ngrids)
    CHARACTER(LEN=1)                      :: prefix
    TYPE(t_vl_register_iter) :: vl_iter
    TYPE(t_hor_interp_meta)  :: hor_interp_meta

    if (dbg_level > 5)  CALL message(routine, "Enter")

    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_elem => NULL()

    !-------------------------------------------------------------
    !-- horizontal interpolation regular lon-lat grids

    ndom = SIZE(p_nh_opt_diag)
    ALLOCATE(ll_varlist(ndom*max_var_ml), ll_vargrid(ndom*max_var_ml), &
      &      ll_varlevs(ndom*max_var_ml), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! flag. will be set if horizontal interpolation tasks have been
    ! created
    l_horintp = .FALSE.

    ! loop over the output definitions, collect pairs of variable
    ! names/lon-lat interpolation requests

    p_onl => first_output_name_list
    nvars_ll            =  0
    n_uv_hrz_intp       =  0
    uv_hrz_intp_grid(:) = -1
    uv_hrz_intp_levs(:) = -1

    NML_LOOP : DO WHILE (ASSOCIATED(p_onl))

      ! Selection criterion:
      ! - lon-lat interpolation is requested
      IF (p_onl%remap == remap_regular_latlon) THEN

        ! do the same for the three level types: model levels (ml),
        ! height levels (hl), pressure levels (pl) and isentropic levels (il):
        DO ilev=1,4

          SELECT CASE(ilev)
          CASE (1)
            varlist => p_onl%ml_varlist
            ilev_type  =  level_type_ml
            max_var    =  max_var_ml
          CASE (2)
            varlist => p_onl%pl_varlist
            ilev_type  =  level_type_pl
            max_var    =  max_var_pl
          CASE (3)
            varlist => p_onl%hl_varlist
            ilev_type  =  level_type_hl
            max_var    =  max_var_hl
          CASE (4)
            varlist => p_onl%il_varlist
            ilev_type  =  level_type_il
            max_var    =  max_var_il
          END SELECT
          IF (varlist(1) == ' ') CYCLE


          ivar_loop: DO ivar=1,max_var
            IF (varlist(ivar) == ' ')             CYCLE
            IF (is_grid_info_var(varlist(ivar)))  CYCLE
            ! check, if have not yet registered this variable:
            DUPLICATE_LOOP : DO i=1,nvars_ll
              IF ((ll_varlist(i) == varlist(ivar))   .AND. &
                & (ll_vargrid(i) == p_onl%lonlat_id) .AND. &
                & (ll_varlevs(i) == ilev_type)) THEN
                CYCLE ivar_loop
              END IF
            END DO DUPLICATE_LOOP

            nvars_ll=nvars_ll+1
            ll_varlist(nvars_ll) = varlist(ivar)
            ll_vargrid(nvars_ll) = p_onl%lonlat_id
            ll_varlevs(nvars_ll) = ilev_type
            ! return a special flag, if var name matches "u" or "v":
            IF (ll_varlist(nvars_ll) == "u" &
              & .OR. ll_varlist(nvars_ll) == "v") THEN
              ! check if this lon-lat grid has not yet been
              ! registered (because we may specify "u" AND "v"):
              DO j=1,n_uv_hrz_intp
                IF ((uv_hrz_intp_grid(n_uv_hrz_intp) == p_onl%lonlat_id) .AND. &
                  & (uv_hrz_intp_levs(n_uv_hrz_intp) == ilev_type)) THEN
                  CYCLE ivar_loop
                END IF
              END DO
              n_uv_hrz_intp = n_uv_hrz_intp + 1
              uv_hrz_intp_grid(n_uv_hrz_intp) = p_onl%lonlat_id
              uv_hrz_intp_levs(n_uv_hrz_intp) = ilev_type
            END IF
          END DO ivar_loop

        END DO
      END IF
      p_onl => p_onl%next
    END DO NML_LOOP

    !-- take care of the special case of "u", "v", where we take "vn"
    !   as the source for horizontal interpolation:
    IF (n_uv_hrz_intp > 0) THEN
      DO i=1,n_uv_hrz_intp
        ! insert post-processing tasks for "vn":
        CALL init_vn_horizontal(uv_hrz_intp_grid(i), uv_hrz_intp_levs(i))
      END DO
    END IF
    !-- loop over requested variables, add variables ("add_var") and
    !-- register interpolation tasks for all requested domains:
    DO ivar=1,nvars_ll
      IF (dbg_level > 8) &
        CALL message(routine, "horizontal interpolation: "&
        &     //"Looking for input var '"//TRIM(ll_varlist(ivar))//"'")
      vname = ll_varlist(ivar)
      !- loop over model level variables
      ! Note that there may be several variables with different time levels,
      ! we just add unconditionally all
      LIST_LOOP : DO WHILE(vl_iter%next())
        ! Do not inspect lists which are disabled for output
        IF (.NOT. vl_iter%cur%p%loutput) CYCLE
        ! Do not inspect lists if vertical level type does not
        ! match (p/z/i/model levels):
        IF (vl_iter%cur%p%vlevel_type/=ll_varlevs(ivar)) CYCLE LIST_LOOP
        ! loop only over variables on requested domains:
        jg = vl_iter%cur%p%patch_id
        IF (.NOT. lonlat_grids%list(ll_vargrid(ivar))%l_dom(jg)) CYCLE LIST_LOOP
        VAR_LOOP : DO iv = 1, vl_iter%cur%p%nvars
          elem => vl_iter%cur%p%vl(iv)%p
          info     => elem%info
          info_dyn => elem%info_dyn
          ! Do not inspect element if it is a container
          IF (info%lcontainer) CYCLE VAR_LOOP
          ! Do not inspect element if "loutput=.false."
          IF (.NOT. info%loutput) CYCLE VAR_LOOP
          ! Do not inspect element if it does not support horizontal
          ! interpolation
          IF (info%hor_interp%hor_intp_type==HINTP_TYPE_NONE) CYCLE VAR_LOOP
          ! Check for matching name
          vname = ll_varlist(ivar)
          ! "u", "v" are processed separately, see above.
          IF ((TRIM(vname) == "u") .OR. (TRIM(vname) == "v"))  CYCLE VAR_LOOP
          ! Check for matching name (take care of suffix of
          ! time-dependent variables):
          IF (TRIM(vname) /= tolower(get_var_name(info))) CYCLE VAR_LOOP
          IF (info%hgrid /= GRID_UNSTRUCTURED_CELL)  CYCLE VAR_LOOP
          ! Found it, add it to the variable list of optional
          ! diagnostics       
          SELECT CASE(ll_varlevs(ivar))
          CASE (level_type_ml)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list         
            prefix = "m"
          CASE (level_type_pl)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_p
            prefix = "p"
          CASE (level_type_hl)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_z
            prefix = "z"
          CASE (level_type_il)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_i
            prefix = "i"
          END SELECT
          ! set local values for "nblks" and "npromz"
          ptr_int_lonlat => lonlat_grids%list(ll_vargrid(ivar))%intp(jg)
          nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
          var_shape      =  info%used_dimensions(:)
          IF (zaxisTypeList%is_2d(info%vgrid) .AND. (info%ndims /= 2)) THEN
            CALL finish(routine, "Inconsistent dimension info: "//TRIM(info%name)//"!")
          END IF
          IF (zaxisTypeList%is_2d(info%vgrid)) THEN
            var_shape(2:3)   =  (/ 1, nblks_lonlat /)
          ELSE
            var_shape(3)     =  nblks_lonlat
          END IF
          ! cross-check: some output fields set "undefined" values as
          ! a fixed value and communicate this to the output module to
          ! create a bit mask (GRIB). However, for lon-lat
          ! interpolated output products this "missval" is preserved
          ! exactly only when nearest-neighbor interpolation is
          ! applied.
          IF (info%lmiss .AND. (info%hor_interp%hor_intp_type /= HINTP_TYPE_LONLAT_NNB)) THEN
            CALL finish(routine, "User tried to interpolate field with missing value!")
          END IF

          hor_interp_meta=create_hor_interp_metadata(                &
            &               hor_intp_type = HINTP_TYPE_NONE,         &
            &               lonlat_id = ll_vargrid(ivar) )

          SELECT CASE (info%hgrid)
          CASE (GRID_UNSTRUCTURED_CELL)
            !--- REAL fields
            IF (ASSOCIATED(elem%r_ptr)) THEN
              CALL add_var( p_opt_diag_list, LONLAT_PREFIX//info%name, p_opt_field_r3d,          &
                &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
                &           ldims=var_shape, lrestart=.FALSE.,                    &
                &           tracer_info=info_dyn%tracer,                          &
                &           loutput=.TRUE., new_element=new_elem,              &
                &           isteptype=info%isteptype,                             &
                &           hor_interp=hor_interp_meta,                           &
                &           vert_interp=info%vert_interp,                         &
                &           post_op=info%post_op,                                 &
                &           lmiss=info%lmiss,                                     &
                &           missval=info%missval%rval, var_class=info%var_class,  &
                &           tlev_source=info%tlev_source, lopenacc=.TRUE.)
            END IF
            !--- INTEGER fields
            IF (ASSOCIATED(elem%i_ptr)) THEN
              CALL add_var( p_opt_diag_list, LONLAT_PREFIX//info%name, p_opt_field_i3d,          &
                &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
                &           ldims=var_shape, lrestart=.FALSE.,                    &
                &           loutput=.TRUE., new_element=new_elem,              &
                &           isteptype=info%isteptype,                             &
                &           hor_interp=hor_interp_meta,                           &
                &           vert_interp=info%vert_interp,                         &
                &           post_op=info%post_op,                                 &
                &           lmiss=info%lmiss,                                     &
                &           missval=info%missval%ival, var_class=info%var_class,  &
                &           tlev_source=info%tlev_source, lopenacc=.TRUE.)
            END IF
            ! SINGLE PRECISION FLOAT fields
            IF (ASSOCIATED(elem%s_ptr)) THEN
              CALL add_var( p_opt_diag_list, LONLAT_PREFIX//info%name, p_opt_field_r3d,          &
                &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
                &           ldims=var_shape, lrestart=.FALSE.,                    &
                &           tracer_info=info_dyn%tracer,                          &
                &           loutput=.TRUE., new_element=new_elem,              &
                &           isteptype=info%isteptype,                             &
                &           hor_interp=hor_interp_meta,                           &
                &           vert_interp=info%vert_interp,                         &
                &           post_op=info%post_op,                                 &
                &           lmiss=info%lmiss,                                     &
                &           missval=REAL(info%missval%sval,wp),                   &
                &           var_class=info%var_class,                             &
                &           tlev_source=info%tlev_source, lopenacc=.TRUE.)
            END IF
            ! LOGICAL fields
            IF (ASSOCIATED(elem%l_ptr)) THEN
              CALL finish(routine, "Regular-grid output of LOGICAL field "//TRIM(info%name)//" unsupported!")
            END IF
          CASE DEFAULT
            CALL finish(routine, "Unsupported grid type!")
          END SELECT
          ! If actions have been defined for the source field, we define those actions 
          ! for the target field (lat-lon-field) as well.
          ! Strictly speaking this is only necessary for the RESET action, but for the 
          ! the time being, we copy all (currently only the RESET action exists).
          ! By this, we assure that for statistically processed fields the  start and end 
          ! interval in the GRIB message is correct.
          ! The drawback is, that the reset action is performed for lon-lat fields as well, 
          ! even though it is not necessary. 
          IF (info%action_list%n_actions > 0 ) new_elem%info%action_list = info%action_list
          !-- create and add post-processing task
          task => pp_task_insert(DEFAULT_PRIORITY4)
          WRITE (task%job_name, *) "horizontal interp. ",TRIM(info%name),", DOM ",jg, &
            &                      " on ", prefix, "-levels, intp TYPE: ", &
            &                      TRIM(STR_HINTP_TYPE(elem%info%hor_interp%hor_intp_type))
          IF (dbg_level > 8) CALL message(routine, task%job_name)
          task%data_input%p_nh_state      => NULL()
          task%data_input%prm_diag        => NULL()
          task%data_input%nh_pzlev_config => NULL()
          task%data_input%jg              =  jg           
          task%data_input%p_patch         => p_patch(jg)
          task%data_input%p_nh_opt_diag   => p_nh_opt_diag(jg)
          task%data_input%p_int_state     => p_int_state(jg)
          task%job_type                   =  TASK_INTP_HOR_LONLAT
          task%activity                   =  new_activity_status(l_output_step=.TRUE.)
          task%activity%check_dom_active  =  .TRUE.
          task%activity%i_timelevel       =  get_var_timelevel(elem%info%name)
          task%data_input%var             => elem       ! set input variable
          task%data_output%var            => new_elem   ! set output variable
          ! Flag. Denotes that at least one interpolation task has
          ! been created.
          l_horintp = .TRUE.
        ENDDO VAR_LOOP ! loop over vlist "i"
      ENDDO LIST_LOOP ! i = 1, SIZE(var_lists)
    END DO ! ivar
    DEALLOCATE(ll_varlist, ll_vargrid, ll_varlevs, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    ! If at least one interpolation task has been created, we add a
    ! setup task which synchronizes the halo regions:
    IF (l_horintp) THEN
      IF (dbg_level >= 10) &
        CALL message(routine, "Creating synchronization task for horizontal interpolation.")
      task => pp_task_insert(DEFAULT_PRIORITY3)
      WRITE (task%job_name, *) "horizontal interp. SYNC"
      IF (dbg_level > 8) CALL message(routine, task%job_name)
      task%job_type = TASK_INTP_SYNC
      task%activity = new_activity_status(l_output_step=.TRUE.)
      task%activity%check_dom_active = .FALSE. ! i.e. no domain-wise (in-)activity 
      task%activity%i_timelevel      = ALL_TIMELEVELS
    END IF
    IF (dbg_level > 5)  CALL message(routine, "Done")
  END SUBROUTINE pp_scheduler_init_lonlat

  !---------------------------------------------------------------
  !> (Internal) Utility routine, collecting variable names from output
  !  namelist.
  SUBROUTINE collect_output_variables(jg, vintp, max_var, nvars, &
    &                                 l_intp, var_names, l_uv_vertical_intp)
    INTEGER, INTENT(IN)      :: jg                              !< current domain
    INTEGER, INTENT(IN) :: vintp                  !< level_type_{pl|hl|il}
    INTEGER, INTENT(IN)      :: max_var                         !< maximum no. of variables
    INTEGER, INTENT(OUT)     :: nvars                           !< actual no. of variables
    LOGICAL, INTENT(OUT)     :: l_intp                          !< Flag. .FALSE. if there is no variable
    CHARACTER(LEN=vname_len), INTENT(OUT) :: var_names(:)     !< list of variable names (strings)
    LOGICAL, INTENT(OUT)     :: l_uv_vertical_intp              !< Flag. .TRUE., if "u" or "v" contained
    ! local variables
    CHARACTER(*), PARAMETER :: routine =  modname//"::collect_output_variables"
    TYPE (t_output_name_list), POINTER :: p_onl
    LOGICAL :: l_jg_active
    INTEGER :: ivar
    CHARACTER(LEN=vname_len), POINTER :: nml_varlist(:)         !< varlist (hl/ml/pl/il) in output_nml namelist

    l_uv_vertical_intp = .FALSE.
    p_onl => first_output_name_list
    nvars = 0
    l_intp = .FALSE.
    IF (.NOT.ANY(vintp .EQ. [level_type_hl, level_type_pl, level_type_il])) &
      & CALL finish(routine, "Internal error!")
    NML_LOOP : DO WHILE (ASSOCIATED(p_onl))
      SELECT CASE (vintp)
      CASE (level_type_hl)
        nml_varlist => p_onl%hl_varlist
      CASE (level_type_pl)
        nml_varlist => p_onl%pl_varlist
      CASE (level_type_il)
        nml_varlist => p_onl%il_varlist
      END SELECT

      IF (dbg_level >= 21)  WRITE (0,*) nml_varlist 

      l_jg_active = (jg == p_phys_patch(p_onl%dom)%logical_id)
      
      ! Selection criteria: 
      ! - domain is requested
      ! - "Z"/"P"/"I"-level interpolation is requested
      IF (l_jg_active .AND. (nml_varlist(1) /= ' ')) THEN
        l_intp = .TRUE.
        DO ivar=1,max_var
          IF (nml_varlist(ivar) == ' ')             CYCLE
          IF (is_grid_info_var(nml_varlist(ivar)))  CYCLE
          nvars=nvars+1
          var_names(nvars) = nml_varlist(ivar)
          ! return a special flag, if var name matches "u" or "v":
          l_uv_vertical_intp = l_uv_vertical_intp &
            & .OR. var_names(nvars) == "u" .OR. var_names(nvars) == "v"
        END DO
      END IF
      p_onl => p_onl%next
    END DO NML_LOOP
    DO ivar = nvars+1,SIZE(var_names)
      var_names(ivar) = " "
    END DO
  END SUBROUTINE collect_output_variables


  !---------------------------------------------------------------
  !> (Internal) Utility routine, add a new variable field for a given
  !  axis, based on the meta-data of an existing variable field (which
  !  is defined on model/half levels).
  !
  SUBROUTINE copy_variable(vname, src_vl, dst_axis, shape3d, ptr, dst_vl)
    CHARACTER(*), INTENT(IN) :: vname        !< name of variable
    TYPE(t_var_list_ptr), INTENT(IN) :: src_vl !< source variable list
    INTEGER, INTENT(IN) :: dst_axis, shape3d(3)
    REAL(wp), POINTER, INTENT(INOUT) :: ptr(:,:,:)  !< reference to field
    TYPE(t_var_list_ptr), INTENT(INOUT) :: dst_vl !< destination variable list
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::copy_variable"
    TYPE(t_var), POINTER :: e
     
    ! find existing variable
    e => find_list_element(src_vl, TRIM(vname))
    IF (.NOT.ASSOCIATED(e)) CALL finish(routine, "Variable not found!")
    ! add new variable, copy the meta-data from the existing variable
    CALL add_var(dst_vl, TRIM(vname), ptr, e%info%hgrid, dst_axis, e%info%cf,     &
      & e%info%grib2, ldims=shape3d, tracer_info=e%info_dyn%tracer,               &
      & post_op=e%info%post_op, loutput=.TRUE., var_class=e%info%var_class,       &
      & tlev_source=e%info%tlev_source, hor_interp=e%info%hor_interp,             &
      & lrestart=.FALSE., vert_interp=e%info%vert_interp, lopenacc=.TRUE.)
  END SUBROUTINE copy_variable

  !---------------------------------------------------------------
  !> (Internal) Utility routine, add a new "vn" field for a given
  !  axis, based on the meta-data of the standard "vn". - This is done
  !  for all available time levels.
  !
  SUBROUTINE init_vn_vertical(jg, job_type, prefix, l_init_prm_diag, nlev, dst_axis, dst_vl)
    INTEGER,          INTENT(IN)  :: jg              !< domain number
    INTEGER,          INTENT(IN)  :: job_type        !< vertical interpolation type
    CHARACTER(LEN=*), INTENT(IN)  :: prefix          !< job name prefix
    LOGICAL,          INTENT(IN)  :: l_init_prm_diag !< Flag. If .TRUE., then prm_diag data structure is available
    INTEGER,          INTENT(IN)  :: nlev            !< number of vertical levels
    INTEGER,          INTENT(IN)  :: dst_axis        !< destination axis
    TYPE(t_var_list_ptr), POINTER :: dst_vl     !< destination variable list
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::init_vn_vertical"
    TYPE(t_var), POINTER :: elem_u, elem_v, elem, vn_elem, new_elem, new_elem2
    INTEGER :: iv, shape3d_c(3), shape3d_e(3), nblks_c, nblks_e, tl, vn_hash
    TYPE(t_job_queue),    POINTER :: task
    TYPE(t_var_metadata), POINTER :: info
    REAL(wp),             POINTER :: p_opt_field_r3d(:,:,:)
    CHARACTER(len=vname_len)    :: name
    CHARACTER(len=4)              :: suffix
    TYPE(t_cf_var)                :: cf
    TYPE(t_grib2_var)             :: grib2
    TYPE(t_post_op_meta)          :: post_op
    TYPE(t_vl_register_iter) :: vl_iter

    if (dbg_level > 5)  CALL message(routine, "Enter")
    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_elem  => NULL()
    new_elem2 => NULL()
    !- find existing variables "u", "v" (for copying the meta-data):
    elem_u => find_list_element(p_nh_state_lists(jg)%diag_list, "u")
    elem_v => find_list_element(p_nh_state_lists(jg)%diag_list, "v")
    !- predefined array shapes
    nblks_c   = p_patch(jg)%nblks_c
    nblks_e   = p_patch(jg)%nblks_e
    shape3d_c = (/ nproma, nlev, nblks_c /)
    shape3d_e = (/ nproma, nlev, nblks_e /)
    vn_hash = text_hash_c('vn')
    !- loop over model level variables
    ! Note that there may be several variables with different time levels,
    ! we just add unconditionally all
    DO WHILE(vl_iter%next())
      ! Do not inspect lists which are disabled for output
      IF (.NOT.vl_iter%cur%p%loutput) CYCLE
      ! loop only over model level variables
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE         
      ! loop only over variables of current domain
      IF (vl_iter%cur%p%patch_id /= jg) CYCLE
      DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)%p
        info => elem%info
        ! Do not inspect element if it is a container
        IF (info%lcontainer) CYCLE
        ! Do not inspect element if "loutput=.false."
        IF (.NOT. info%loutput) CYCLE
        ! Check for matching name (take care of suffix of
        ! time-dependent variables):
        IF (vn_hash .NE.  vl_iter%cur%p%key_notl(iv)) CYCLE
        ! get time level
        tl = get_var_timelevel(info%name)
        suffix = ''
        IF (tl /= -1) suffix = get_timelevel_string(tl)
        !-- create a new z/p/i-variable "vn":
        CALL add_var(dst_vl, TRIM(info%name), p_opt_field_r3d, elem%info%hgrid, dst_axis,         &
          & info%cf, info%grib2, ldims=shape3d_e, new_element=vn_elem, post_op=info%post_op,      &
          & lrestart=.FALSE., var_class=info%var_class, tlev_source=elem%info%tlev_source,        &
          & hor_interp=info%hor_interp, vert_interp=info%vert_interp, lopenacc=.TRUE.)
        !-- create a post-processing task for vertical interpolation of "vn"
        task => pp_task_insert(DEFAULT_PRIORITY1)
        task%job_name        =  &
          &  TRIM(prefix)//" interp. "//TRIM(info%name)//", DOM "//TRIM(int2string(jg))
        IF (dbg_level > 8) CALL message(routine, task%job_name)
        task%job_type                    =  job_type
        task%activity                    =  new_activity_status(l_output_step=.TRUE.)
        task%activity%check_dom_active   =  .TRUE.
        task%activity%i_timelevel        =  tl
        task%data_input%jg               =  jg 
        task%data_input%p_patch          => p_patch(jg)          
        task%data_input%p_int_state      => p_int_state(jg)
        task%data_input%p_nh_state       => p_nh_state(jg)
        task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
        task%data_input%p_nh_opt_diag    => p_nh_opt_diag(jg)
        IF (l_init_prm_diag) THEN
          task%data_input%prm_diag       => prm_diag(jg)
        ELSE
          task%data_input%prm_diag       => NULL() 
        END IF
        task%data_input%var  => elem      ! set input variable
        task%data_output%var => vn_elem   ! set output variable
        !-- create new cell-based variables "u", "v" on the same vertical axis
        name    = TRIM(get_var_name(elem_u%info))//suffix
        cf      = elem_u%info%cf
        grib2   = elem_u%info%grib2
        post_op = elem_u%info%post_op
        CALL add_var(dst_vl, TRIM(name), p_opt_field_r3d, GRID_UNSTRUCTURED_CELL,         &
          & dst_axis, cf, grib2, ldims=shape3d_c, lrestart=.FALSE., lopenacc=.TRUE.,      &
          & in_group=elem_u%info%in_group, new_element=new_elem, post_op=post_op,         &
          & var_class=elem_u%info%var_class, tlev_source=elem_u%info%tlev_source,         &
          & hor_interp=elem_u%info%hor_interp, vert_interp=elem_u%info%vert_interp )
        name    = TRIM(get_var_name(elem_v%info))//suffix
        cf      = elem_v%info%cf
        grib2   = elem_v%info%grib2
        post_op = elem_v%info%post_op
        CALL add_var(dst_vl, TRIM(name), p_opt_field_r3d, GRID_UNSTRUCTURED_CELL,         &
          & dst_axis, cf, grib2, ldims=shape3d_c, lrestart=.FALSE., lopenacc=.TRUE.,      &
          & in_group=elem_v%info%in_group, new_element=new_elem2, post_op=post_op,        &
          & var_class=elem_v%info%var_class, tlev_source=elem_v%info%tlev_source,         &
          & hor_interp=elem_v%info%hor_interp, vert_interp=elem_v%info%vert_interp)
        !-- create a post-processing task for edge2cell interpolation "vn" -> "u","v"
        task => pp_task_insert(DEFAULT_PRIORITY2)
        WRITE (task%job_name, *) "edge2cell interp. ",TRIM(info%name),", DOM ",jg
        IF (dbg_level > 8) CALL message(routine, task%job_name)
        task%data_input%p_nh_state      => NULL()
        task%data_input%prm_diag        => NULL()
        task%data_input%nh_pzlev_config => NULL()
        task%data_input%jg              =  jg           
        task%data_input%p_patch         => p_patch(jg)
        task%data_input%p_nh_opt_diag   => p_nh_opt_diag(jg)
        task%data_input%p_int_state     => p_int_state(jg)
        task%job_type                   =  TASK_INTP_EDGE2CELL
        task%activity                   =  new_activity_status(l_output_step=.TRUE.)
        task%activity%check_dom_active  =  .TRUE.
        task%activity%i_timelevel       =  tl
        task%data_input%var             => vn_elem    ! set input variable
        task%data_output%var            => new_elem   ! set output variable
        task%data_output%var_2          => new_elem2  ! set Y-component
      END DO
    END DO
    if (dbg_level > 5)  CALL message(routine, "Done")
  END SUBROUTINE init_vn_vertical

  !---------------------------------------------------------------
  !> Setup of i/p/z-level interpolation tasks.
  !  - collects lists of variables for i/p/z interpolation
  !  - adds variable fields where interpolation results will be stored
  !  - creates "post-processing tasks", i.e. entries in a list 
  !    which is regularly traversed during the model run.
  !
  ! See SUBROUTINE pp_scheduler_init for further details.
  !
  SUBROUTINE pp_scheduler_init_ipz(l_init_prm_diag)
    LOGICAL, INTENT(IN) :: l_init_prm_diag !< Flag. If .TRUE., then prm_diag data structure is available

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  modname//"::pp_scheduler_init_ipz"
    INTEGER,      PARAMETER :: init_tasks(3) = &
      &  (/ TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I /)
    CHARACTER,    PARAMETER :: init_names(3) = &
      &  (/ 'z', 'p', 'i' /)
    INTEGER                            :: &
      &  jg, ndom, ibits, ierrstat, ivar, i, isteptype, iv, &
      &  iaxis, vgrid, nlev, nvars_pl, nvars_hl, nvars_il, nvars,   &
      &  job_type, shape3d(3), datatype_flt, intp_id
    LOGICAL                            :: &
      &  l_intp_p, l_intp_z, l_intp_i, found, &
      &  l_uv_vertical_intp_z, l_uv_vertical_intp_p, l_uv_vertical_intp_i, &
      &  l_uv_vertical_intp
    TYPE(t_job_queue),         POINTER :: task
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_var_list_ptr),          POINTER :: p_opt_diag_list_p, p_opt_diag_list_z, &
      &                                   p_opt_diag_list_i, p_opt_diag_list
    REAL(wp),                  POINTER :: p_opt_field_r3d(:,:,:)
    TYPE(t_var), POINTER :: elem, new_elem
    ! variable lists (for all domains + output name lists):
    CHARACTER(LEN=vname_len), TARGET, ALLOCATABLE  :: &
         &                                pl_varlist(:), hl_varlist(:), il_varlist(:)
    CHARACTER(LEN=vname_len),  POINTER   :: varlist(:)
    CHARACTER(LEN=10)                    :: prefix
    TYPE(t_var_metadata),POINTER         :: info
    TYPE(t_var_metadata_dynamic),POINTER :: info_dyn
    TYPE(t_cf_var)                       :: cf_desc
    TYPE(t_grib2_var)                    :: grib2_desc
    TYPE(t_vl_register_iter) :: vl_iter

    ! define NetCDF output precision
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_elem => NULL()
    if (dbg_level > 5)  CALL message(routine, "Enter")
    ndom = SIZE(p_nh_opt_diag)
    ! loop over the domains and output definitions
    ! - check for which domain p- or z-level interpolation has been
    !   requested
    ! - register setup routine for pz-level interpolation (must be
    !   executed ahead of interpolation tasks)
    ! - for each variable: register post-processing task

    ALLOCATE(pl_varlist(ndom*max_var_pl), hl_varlist(ndom*max_var_hl), &
         &   il_varlist(ndom*max_var_il), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! list indices for the vertical interpolation types:
    shape3d(1) = nproma
    DOM_LOOP : DO jg=1,ndom
      IF (dbg_level > 8)  CALL message(routine, "DOM "//int2string(jg))
      !-- check if any output name list requests p- or z- or i-level 
      !-- interpolation for this domain, collect the list of variables
      ! loop in search of pressure-level interpolation      
      CALL collect_output_variables(jg, level_type_pl, max_var_pl, &                      ! in
        &                           nvars_pl, l_intp_p, pl_varlist, l_uv_vertical_intp_p) ! out
      ! loop in search of height-level interpolation
      CALL collect_output_variables(jg, level_type_hl, max_var_hl, &                      ! in
        &                           nvars_hl, l_intp_z, hl_varlist, l_uv_vertical_intp_z) ! out
      ! loop in search of i-level interpolation
      CALL collect_output_variables(jg, level_type_il, max_var_il, &                      ! in
        &                           nvars_il, l_intp_i, il_varlist, l_uv_vertical_intp_i) ! out
      ! now, we have total variables lists "hl_varlist(1:nvars_hl)"
      ! and "pl_varlist(1:nvars_pl)" and "il_varlist(1:nvars_il)"
      ! some debugging output...
      IF (dbg_level > 8)  THEN
        DO i=1,nvars_pl
          WRITE (message_text,*) "p var list: ", pl_varlist(i)
          CALL message(routine, message_text)
        END DO
        DO i=1,nvars_hl
          WRITE (message_text,*) "h var list: ", hl_varlist(i)
          CALL message(routine, message_text)
        END DO
        DO i=1,nvars_il
          WRITE (message_text,*) "i var list: ", il_varlist(i)
          CALL message(routine, message_text)
        END DO
      END IF
      ! skip domain if no p/z-interpolation requested:
      IF (.NOT. (l_intp_z .OR. l_intp_p .OR. l_intp_i)) CYCLE DOM_LOOP
      ! remove duplicates from variable lists
      IF (l_intp_z) CALL remove_duplicates(hl_varlist, nvars_hl)
      IF (l_intp_p) CALL remove_duplicates(pl_varlist, nvars_pl)
      IF (l_intp_i) CALL remove_duplicates(il_varlist, nvars_il)
      !-- First, add some diagnostic variables which are essential for
      !-- p/z-level interpolation, e.g. temp_z, pres_z:
      p_diag_pz         => p_nh_opt_diag(jg)%diag_pz
      p_opt_diag_list_z => p_nh_opt_diag(jg)%opt_diag_list_z
      p_opt_diag_list_p => p_nh_opt_diag(jg)%opt_diag_list_p
      p_opt_diag_list_i => p_nh_opt_diag(jg)%opt_diag_list_i
      ibits     = DATATYPE_PACK16   ! "entropy" of horizontal slice
      ! predefined array shapes
      shape3d(3) = p_patch(jg)%nblks_c
      ! add new variable fields for the z/p/i-axis, based on the
      ! meta-data of an existing variable field (which is defined on
      ! model/half levels):
      IF (l_intp_z) THEN
        shape3d(2) = nh_pzlev_config(jg)%zlevels%nvalues
        CALL copy_variable("temp", p_nh_state_lists(jg)%diag_list, ZA_ALTITUDE, shape3d, &
          &                p_diag_pz%z_temp, p_opt_diag_list_z)
        CALL copy_variable("pres", p_nh_state_lists(jg)%diag_list, ZA_ALTITUDE, shape3d, &
          &                p_diag_pz%z_pres, p_opt_diag_list_z)
      END IF
      IF (l_intp_p) THEN
        shape3d(2) = nh_pzlev_config(jg)%plevels%nvalues
        cf_desc    = t_cf_var('gh', 'm', 'geopotential height', datatype_flt)
        grib2_desc = grib2_var(0, 3, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'gh', p_diag_pz%p_gh,                  &
          & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc, grib2_desc,           &
          & ldims=shape3d, lrestart=.FALSE., lopenacc=.TRUE.)
        !$ACC ENTER DATA ATTACH(p_diag_pz%p_gh)
        CALL copy_variable("temp",   p_nh_state_lists(jg)%diag_list,    ZA_PRESSURE, shape3d, &
          &                p_diag_pz%p_temp, p_opt_diag_list_p)
      END IF
      IF (l_intp_i) THEN
        shape3d(2) = nh_pzlev_config(jg)%ilevels%nvalues
        cf_desc    = t_cf_var('gh', 'm', 'geopotential height', datatype_flt)
        grib2_desc = grib2_var(0, 3, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_opt_diag_list_i, 'gh', p_diag_pz%i_gh,                  &
          & GRID_UNSTRUCTURED_CELL, ZA_ISENTROPIC, cf_desc, grib2_desc,         &
          & ldims=shape3d, lrestart=.FALSE., lopenacc=.TRUE.)
        !$ACC ENTER DATA ATTACH(p_diag_pz%i_gh)
        CALL copy_variable("temp",   p_nh_state_lists(jg)%diag_list,    ZA_ISENTROPIC, shape3d, &
          &                p_diag_pz%i_temp, p_opt_diag_list_i)
      END IF

      !-- register interpolation setup as post-processing task(s)
      DO i=1,SIZE(init_tasks)
        IF ((init_tasks(i) == TASK_INIT_VER_Z) .AND. .NOT. l_intp_z) CYCLE
        IF ((init_tasks(i) == TASK_INIT_VER_P) .AND. .NOT. l_intp_p) CYCLE
        IF ((init_tasks(i) == TASK_INIT_VER_I) .AND. .NOT. l_intp_i) CYCLE

        task => pp_task_insert(HIGH_PRIORITY)
        task%data_input%p_int_state      => p_int_state(jg)
        task%data_input%jg               =  jg           
        task%data_input%p_patch          => p_patch(jg)
        task%data_input%p_nh_state       => p_nh_state(jg)
        task%data_input%p_nh_opt_diag    => p_nh_opt_diag(jg)
        IF (l_init_prm_diag) THEN
          task%data_input%prm_diag       => prm_diag(jg)
        ELSE
          task%data_input%prm_diag       => NULL() 
        END IF
        task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
        task%activity     = new_activity_status(l_output_step=.TRUE.)
        task%activity%check_dom_active   = .TRUE.
        task%activity%i_timelevel        = ALL_TIMELEVELS
        task%job_type                    = init_tasks(i)
        task%job_name                    = "Init: "//init_names(i)//"-level interpolation, DOM "//TRIM(int2string(jg))
        IF (dbg_level > 8) CALL message(routine, task%job_name)
      END DO

      !-- register clean-up routine as a post-processing task
      task => pp_task_insert(LOW_PRIORITY)
      task%activity     = new_activity_status(l_output_step=.TRUE.)
      task%activity%check_dom_active   =  .TRUE.
      task%activity%i_timelevel        =  ALL_TIMELEVELS
      task%data_input%p_nh_opt_diag => p_nh_opt_diag(jg)
      task%job_name     = "Clean-up: ipz-level interpolation, level "//TRIM(int2string(jg))
      IF (dbg_level > 8) CALL message(routine, task%job_name)
      task%job_type     = TASK_FINALIZE_IPZ
      task%data_input%p_int_state      => NULL()
      task%data_input%jg               =  jg           
      task%data_input%p_patch          => p_patch(jg)
      task%data_input%p_nh_state       => p_nh_state(jg)
      task%data_input%prm_diag         => NULL() 
      task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)

      ! remove already defined variables from list of requested output
      ! fields:
      IF (l_intp_z) CALL difference(hl_varlist, nvars_hl, &
        &                           (/ "temp  ", "pres  ", "u     ", "v     " /), 4)
      IF (l_intp_p) CALL difference(pl_varlist, nvars_pl, &
        &                           (/ "gh    ", "temp  ", "u     ", "v     " /), 4)
      IF (l_intp_i) CALL difference(il_varlist, nvars_il, &
        &                           (/ "gh    ", "temp  ", "u     ", "v     " /), 4)

      !-- loop over requested p-, z-, and i-level variables, add variables
      !-- ("add_var") and register interpolation tasks:
      DO iaxis = 1, 3
        SELECT CASE(iaxis)
        CASE(1)
          prefix  =  "z-level"
          varlist => hl_varlist
          nvars   =  nvars_hl
          nlev    =  nh_pzlev_config(jg)%zlevels%nvalues
          vgrid   =  ZA_ALTITUDE
          p_opt_diag_list => p_opt_diag_list_z
          job_type = TASK_INTP_VER_ZLEV
          l_uv_vertical_intp = l_uv_vertical_intp_z
!CDIR NOIEXPAND
          intp_id = vintp_type_id("Z")
        CASE(2)
          prefix  =  "p-level"
          varlist => pl_varlist
          nvars   =  nvars_pl
          nlev    =  nh_pzlev_config(jg)%plevels%nvalues
          vgrid   =  ZA_PRESSURE
          p_opt_diag_list => p_opt_diag_list_p
          job_type = TASK_INTP_VER_PLEV
          l_uv_vertical_intp = l_uv_vertical_intp_p
!CDIR NOIEXPAND
          intp_id = vintp_type_id("P")
        CASE(3)
          prefix  =  "i-level"
          varlist => il_varlist
          nvars   =  nvars_il
          nlev    =  nh_pzlev_config(jg)%ilevels%nvalues
          vgrid   =  ZA_ISENTROPIC
          p_opt_diag_list => p_opt_diag_list_i
          job_type = TASK_INTP_VER_ILEV
          l_uv_vertical_intp = l_uv_vertical_intp_i
!CDIR NOIEXPAND
          intp_id = vintp_type_id("I")
        END SELECT
        shape3d(2) = nlev
        !-- if "u", "v" appear in the variable list...
        IF (l_uv_vertical_intp) &
          ! for each time level, create a new z/p/i-variable "vn",
          ! create new cell-based variables "u", "v" on the same
          ! vertical axis, create a post-processing task for vertical
          ! interpolation of "vn", create a post-processing task for
          ! edge2cell interpolation "vn" -> "u","v":
          & CALL init_vn_vertical(jg, job_type, prefix, l_init_prm_diag, &
            &                   nlev, vgrid, p_opt_diag_list)

        DO ivar=1,nvars
          IF (dbg_level > 8) &
            CALL message(routine, TRIM(prefix)//": Looking for input var '"//TRIM(varlist(ivar))//"'")
          found = .FALSE.
        
          !- loop over model level variables
          ! Note that there may be several variables with different time levels,
          ! we just add unconditionally all 
          ! TODO(HB): this is n(2) complexity... mabe try to overcome...
          DO WHILE(vl_iter%next())
            ! Do not inspect lists which are disabled for output
            IF (.NOT. vl_iter%cur%p%loutput) CYCLE
            ! loop only over model level variables
            IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE         
            ! loop only over variables of current domain
            IF (vl_iter%cur%p%patch_id /= jg) CYCLE

            DO iv = 1, vl_iter%cur%p%nvars
              elem     => vl_iter%cur%p%vl(iv)%p
              info     => elem%info
              info_dyn => elem%info_dyn
              ! Do not inspect element if it is a container
              IF (info%lcontainer) CYCLE
              ! Do not inspect element if "loutput=.false."
              IF (.NOT. info%loutput) CYCLE
              ! Inspect element only if vertical interpolation matches
              IF (.NOT. info%vert_interp%vert_intp_type(intp_id)) CYCLE
              ! Check for matching name (take care of suffix of
              ! time-dependent variables):
              IF (TRIM(varlist(ivar)) /= tolower(get_var_name(info))) CYCLE
              ! Found it, add it to the variable list of optional
              ! diagnostics

              IF (ANY(info%used_dimensions(1:3:2) .NE. shape3d(1:3:2))) &
                & CALL finish(routine, "Unexpected field size!")
              ! fields interpolated to pressure levels are time
              ! dependent, rather than constant in time:
              isteptype = MERGE(TSTEP_INSTANT, info%isteptype, &
                & info%isteptype .EQ. TSTEP_CONSTANT .AND. (l_intp_p .OR. l_intp_i))

              CALL add_var( p_opt_diag_list, info%name, p_opt_field_r3d,    &
                &           info%hgrid, vgrid, info%cf, info%grib2,         &
                &           ldims=shape3d, lrestart=.FALSE.,                &
                &           tracer_info=info_dyn%tracer,                    &
                &           loutput=.TRUE., new_element=new_elem,        &
                &           isteptype=isteptype,                            &
                &           post_op=info%post_op, var_class=info%var_class, &
                &           tlev_source=info%tlev_source,                   &
                &           hor_interp=info%hor_interp,                     &
                &           vert_interp=info%vert_interp, lopenacc=.TRUE.)

              !-- add post-processing task for interpolation

              task => pp_task_insert(DEFAULT_PRIORITY1)
              task%job_name        =  &
                &  TRIM(prefix)//" interp. "//TRIM(info%name)  &
                &  //", DOM "//TRIM(int2string(jg))            &
                & //", vintp type: "//TRIM(int2string(elem%info%vert_interp%vert_intp_method))
              IF (dbg_level > 8) CALL message(routine, task%job_name)

              task%job_type                    =  job_type
              task%activity                    =  new_activity_status(l_output_step=.TRUE.)
              task%activity%check_dom_active   =  .TRUE.
              task%activity%i_timelevel        =  get_var_timelevel(info%name)
              task%data_input%jg               =  jg 
              task%data_input%p_patch          => p_patch(jg)          
              task%data_input%p_int_state      => p_int_state(jg)
              task%data_input%p_nh_state       => p_nh_state(jg)
              task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
              task%data_input%p_nh_opt_diag    => p_nh_opt_diag(jg)
              IF (l_init_prm_diag) THEN
                task%data_input%prm_diag       => prm_diag(jg)
              ELSE
                task%data_input%prm_diag       => NULL()
              END IF
              task%data_input%var  => elem       ! set input variable
              task%data_output%var => new_elem   ! set output variable
              found = .TRUE.
            ENDDO ! loop over vlist "i"
          ENDDO ! i = 1, SIZE(var_lists)
          ! Check that at least one element with this name has been found
          IF(.NOT. found) &
            CALL finish(routine,TRIM(prefix)//" interpolation: No feasible variable found: "&
            &                   //TRIM(varlist(ivar)))
        END DO ! ivar
      END DO ! height/pressure/isentropic axis
    END DO DOM_LOOP  ! jg
    DEALLOCATE(pl_varlist, hl_varlist, il_varlist,STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    IF (dbg_level > 5)  CALL message(routine, "Done")
  END SUBROUTINE pp_scheduler_init_ipz

  !---------------------------------------------------------------
  !> Register a new post-processing task for computing additional 
  !  diagnostic fields (not for interpolation).
  !
  !  The new task will be added to the dynamic list of post-processing
  !  jobs and called on a regular basis.
  !
  SUBROUTINE pp_scheduler_register(name, jg, p_out_var,             &
    &                              l_init_prm_diag, job_type,       &
    &                              opt_priority, opt_l_output_step, &
    &                              opt_l_accumulation_step,         &
    &                              opt_p_in_var)

    CHARACTER(LEN=*)                   , INTENT(IN) :: name
    INTEGER                            , INTENT(IN) :: jg
    TYPE (t_var), POINTER                  :: p_out_var
    LOGICAL                            , INTENT(IN) :: l_init_prm_diag
    INTEGER                            , INTENT(IN) :: job_type
    INTEGER, OPTIONAL                  , INTENT(IN) :: opt_priority
    LOGICAL, OPTIONAL                  , INTENT(IN) :: opt_l_output_step
    LOGICAL, OPTIONAL                  , INTENT(IN) :: opt_l_accumulation_step
    TYPE (t_var), POINTER, OPTIONAL        :: opt_p_in_var
    ! local variables
    LOGICAL                    :: l_output_step
    LOGICAL                    :: l_accumulation_step
    INTEGER                    :: priority
    TYPE(t_job_queue), POINTER :: task
    
    ! set default values
    l_output_step = .TRUE.
    l_accumulation_step = .TRUE.
    priority      = DEFAULT_PRIORITY0
    ! assign optional parameters:
    IF (PRESENT(opt_priority)) priority = opt_priority
    IF (PRESENT(opt_l_output_step)) l_output_step = opt_l_output_step
    IF (PRESENT(opt_l_accumulation_step)) l_accumulation_step = opt_l_accumulation_step
    ! create a post-processing task and fill its input/output data:
    task => pp_task_insert(priority)
    WRITE (task%job_name, *) TRIM(name),", DOM ",jg
    IF (PRESENT(opt_p_in_var)) THEN
      task%data_input%var            => opt_p_in_var
    END IF
    task%data_input%jg               =  jg           
    task%data_input%p_nh_state       => p_nh_state(jg)
    task%data_input%p_patch          => p_patch(jg)
    task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
    IF (l_init_prm_diag) THEN
      task%data_input%prm_diag       => prm_diag(jg)
    END IF
    task%data_output%var             => p_out_var
    task%job_type                    =  job_type
    task%activity                    =  new_activity_status( &
      &   l_output_step=l_output_step, &
      &   l_accumulation_step=l_accumulation_step &
      )
    task%activity%check_dom_active   =  .TRUE.
    task%activity%i_timelevel        =  ALL_TIMELEVELS
  END SUBROUTINE pp_scheduler_register


  !---------------------------------------------------------------
  !> Loop over job queue, call active tasks.
  !
  SUBROUTINE pp_scheduler_process(simulation_status, lacc)
    TYPE(t_simulation_status), INTENT(IN) :: simulation_status
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_scheduler_process"
    TYPE(t_job_queue), POINTER :: ptr_task
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    CALL assert_acc_device_only(routine, lacc)
    IF (dbg_level >= 10) THEN
      CALL message(routine,"Processing task list...")
    END IF
    ptr_task => job_queue
    ! loop over job queue
    LOOP_JOB : DO
      IF (.NOT. ASSOCIATED(ptr_task)) EXIT
      IF (.NOT. pp_task_is_active(ptr_task, simulation_status)) THEN
        IF (dbg_level > 20) THEN
          WRITE(message_text,*) "Skipping task '", TRIM(ptr_task%job_name), "'"
          CALL message(routine, TRIM(message_text))
        END IF
        ptr_task => ptr_task%next
        CYCLE LOOP_JOB
      END IF

      IF (dbg_level > 5) THEN
        WRITE(message_text,*) "Staging task '", TRIM(ptr_task%job_name), "'"
        CALL message(routine, TRIM(message_text))
      END IF

      SELECT CASE ( ptr_task%job_type )

        ! initialize vertical interpolation:
      CASE ( TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I, &
           & TASK_FINALIZE_IPZ)
        CALL pp_task_ipzlev_setup(ptr_task, lacc=.TRUE.)

        ! perform horizontal interpolation:
      CASE ( TASK_INTP_HOR_LONLAT )
        CALL pp_task_lonlat(ptr_task)

        ! perform vertical interpolation:
      CASE ( TASK_INTP_VER_PLEV, TASK_INTP_VER_ZLEV, TASK_INTP_VER_ILEV )
        CALL pp_task_ipzlev(ptr_task, lacc=.TRUE.)

        ! synchronize halo regions:
      CASE ( TASK_INTP_SYNC )
        CALL pp_task_sync(simulation_status)

        ! compute mean sea level pressure:
      CASE ( TASK_INTP_MSL )
        CALL pp_task_intp_msl(ptr_task, lacc=.TRUE.)

        ! compute relative humidty, vertical velocity, potential vorticity, ...
      CASE ( TASK_COMPUTE_RH, TASK_COMPUTE_OMEGA, TASK_COMPUTE_PV, TASK_COMPUTE_SDI2,              &
        &    TASK_COMPUTE_LPI, TASK_COMPUTE_CEILING, TASK_COMPUTE_HBAS_SC, TASK_COMPUTE_HTOP_SC,   &
        &    TASK_COMPUTE_TWATER, TASK_COMPUTE_Q_SEDIM, TASK_COMPUTE_DBZ850, TASK_COMPUTE_DBZCMAX, &
        &    TASK_COMPUTE_DBZLMX_LOW, TASK_COMPUTE_SMI, TASK_COMPUTE_WSHEAR_U,                     &
        &    TASK_COMPUTE_WSHEAR_V, TASK_COMPUTE_SRH, TASK_COMPUTE_LAPSERATE, TASK_COMPUTE_VIS,    &
        &    TASK_COMPUTE_INVERSION, TASK_COMPUTE_MCONV)
        IF (timers_level >= 5) CALL timer_start(timer_opt_diag_atmo)
        CALL pp_task_compute_field(ptr_task)
        IF (timers_level >= 5) CALL timer_stop(timer_opt_diag_atmo)

        ! vector reconstruction on cell centers:
      CASE ( TASK_INTP_EDGE2CELL )
        CALL pp_task_edge2cell(ptr_task)

      CASE DEFAULT
        CALL finish(routine, "Unknown post-processing job.")
      END SELECT

      ptr_task => ptr_task%next
    END DO LOOP_JOB

  END SUBROUTINE pp_scheduler_process


  !---------------------------------------------------------------
  !> Destruction of post-processing job queue
  !
  SUBROUTINE pp_scheduler_finalize()
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_scheduler_finalize"
    INTEGER                    :: ierrstat
    TYPE(t_job_queue), POINTER :: tmp
    
    CALL message(routine, "")
    ! destroy linked list
    DO
      IF (.NOT. ASSOCIATED(job_queue)) EXIT
      ! remove list head
      tmp => job_queue%next
      DEALLOCATE(job_queue, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      job_queue => tmp
    END DO
    
  END SUBROUTINE pp_scheduler_finalize


  !--- UTILITY ROUTINES --------------------------------------------------------------

  !---------------------------------------------------------------
  !> @return .TRUE. if given post-processing task is in active state.
  ! 
  ! Tasks may be inactive, e.g. outside the output intervals.
  FUNCTION pp_task_is_active(ptr_task, sim_status)
    LOGICAL :: pp_task_is_active
    TYPE(t_job_queue), POINTER :: ptr_task
    TYPE(t_simulation_status),  INTENT(IN) :: sim_status
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_is_active"
    INTEGER :: jg, tlev_source, timelevel

    ! compare simulation status to post-processing tasks activity
    ! flags, then check if any of the activity conditions is
    ! fulfilled:
    pp_task_is_active = .FALSE.
    IF (ANY(ptr_task%activity%status_flags(:)  .AND.  &
      &     sim_status%status_flags(:))) pp_task_is_active = .TRUE.

    ! check, if current task applies only to domains which are
    ! "active":
    IF (ptr_task%activity%check_dom_active) THEN
      jg          = ptr_task%data_input%jg

      IF (.NOT. sim_status%ldom_active(jg)) THEN
        pp_task_is_active = .FALSE.
      END IF

      IF  (ptr_task%activity%i_timelevel /= ALL_TIMELEVELS) THEN
         tlev_source = ptr_task%data_input%var%info%tlev_source

         SELECT CASE (tlev_source)
         CASE(TLEV_NNOW);     timelevel = sim_status%i_timelevel_dyn(jg)
         CASE(TLEV_NNOW_RCF); timelevel = sim_status%i_timelevel_phy(jg)
         CASE DEFAULT
            CALL finish(routine, 'Unsupported tlev_source')
         END SELECT
         
         ! check, if current task matches the variable time level (TL1,
         ! TL2, ...) of the simulation status:
         IF  (ptr_task%activity%i_timelevel /= timelevel) THEN
            pp_task_is_active = .FALSE.
         END IF
      END IF
    END IF

  END FUNCTION pp_task_is_active  


  !---------------------------------------------------------------
  !> Insert task into job queue.
  FUNCTION pp_task_insert(job_priority) RESULT(element)
    INTEGER, INTENT(IN) :: job_priority
    TYPE(t_job_queue), POINTER :: element
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_insert"
    TYPE(t_job_queue), POINTER :: tmp, nb_left
    INTEGER                    :: ierrstat

    IF (dbg_level > 5)  CALL message(routine, "Inserting pp task")
    
    ! find the correct position in list:
    tmp     => job_queue
    nb_left => NULL()
    DO
      IF (.NOT. ASSOCIATED(tmp))                EXIT
      IF (tmp%job_priority > job_priority) EXIT
      nb_left => tmp
      tmp     => tmp%next
    END DO
    ! insert element into linked list
    IF (ASSOCIATED(nb_left)) THEN
      ALLOCATE(nb_left%next, stat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      element => nb_left%next
    ELSE
      ALLOCATE(job_queue, stat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      element => job_queue
    END IF
    ! fill new list element with data
    element%job_priority            =  job_priority
    element%next                    => tmp
    element%activity%i_timelevel    =  UNDEF_TIMELEVEL
    element%data_input%jg           =  -1
  END FUNCTION pp_task_insert


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for "t_simulation_status" variables
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION new_simulation_status(l_output_step, l_first_step, l_last_step, l_accumulation_step,        &
    &                            l_dom_active, i_timelevel_dyn, i_timelevel_phy)  &
    RESULT(sim_status)

    TYPE(t_simulation_status) :: sim_status
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_output_step, l_first_step, l_last_step, l_accumulation_step
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_dom_active(:)
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  i_timelevel_dyn(:), i_timelevel_phy(:)
    ! local variables
    INTEGER :: ndom

    ! set default values
    sim_status%status_flags(:) = (/ .FALSE., .FALSE., .FALSE. , .FALSE./)

    ! supersede with user definitions
    IF (PRESENT(l_output_step)) sim_status%status_flags(1) = l_output_step
    IF (PRESENT(l_first_step)) sim_status%status_flags(2) = l_first_step
    IF (PRESENT(l_last_step)) sim_status%status_flags(3) = l_last_step
    IF (PRESENT(l_accumulation_step)) sim_status%status_flags(4) = l_accumulation_step

    ! as a default, all domains are "inactive", i.e. the activity
    ! flags are not considered:
    sim_status%ldom_active(:) = .FALSE.
    IF  (PRESENT(l_dom_active)) THEN
      ndom = SIZE(l_dom_active)
      sim_status%ldom_active(1:ndom) = l_dom_active(1:ndom)
    END IF

    ! as a default, no special timelevel is set for (in-)activity:
    sim_status%i_timelevel_dyn(:) = ALL_TIMELEVELS
    sim_status%i_timelevel_phy(:) = ALL_TIMELEVELS
    IF (PRESENT(i_timelevel_dyn)) sim_status%i_timelevel_dyn = i_timelevel_dyn
    IF (PRESENT(i_timelevel_phy)) sim_status%i_timelevel_phy = i_timelevel_phy

  END FUNCTION new_simulation_status


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for "t_simulation_status" variables
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION new_activity_status(l_output_step, l_first_step, l_last_step, l_accumulation_step, &
    &                          check_dom_active, i_timelevel)  &
    RESULT(activity_status)

    TYPE(t_activity_status) :: activity_status
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_output_step, l_first_step, l_last_step, l_accumulation_step
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  check_dom_active
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  i_timelevel

    ! set default values
    activity_status%status_flags(:) = (/ .FALSE., .FALSE., .FALSE., .FALSE. /)

    ! supersede with user definitions
    IF (PRESENT(l_output_step)) activity_status%status_flags(1) = l_output_step
    IF (PRESENT(l_first_step)) activity_status%status_flags(2) = l_first_step
    IF (PRESENT(l_last_step)) activity_status%status_flags(3) = l_last_step
    IF (PRESENT(l_accumulation_step)) activity_status%status_flags(4) = l_accumulation_step

    ! as a default, all domains are "inactive", i.e. the activity
    ! flags are not considered:
    activity_status%check_dom_active = .FALSE.
    IF (PRESENT(check_dom_active)) activity_status%check_dom_active = check_dom_active

    ! as a default, no special timelevel is set for (in-)activity:
    activity_status%i_timelevel = ALL_TIMELEVELS
    IF (PRESENT(i_timelevel)) activity_status%i_timelevel = i_timelevel

  END FUNCTION new_activity_status

END MODULE mo_pp_scheduler
