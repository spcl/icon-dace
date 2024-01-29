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

MODULE mo_run_nml

  USE mo_run_config, ONLY: config_ltestcase       => ltestcase,       &
                         & config_ldynamics       => ldynamics,       &
                         & config_iforcing        => iforcing,        &
                         & config_ltransport      => ltransport,      &
                         & config_ntracer         => ntracer,         &
                         & config_lart            => lart,            &
                         & config_ldass_lhn       => ldass_lhn,       &
                         & config_luse_radarfwo   => luse_radarfwo,   &
                         & config_radarnmlfile    => radarnmlfile,    &
                         & config_lvert_nest      => lvert_nest,      &
                         & config_nlev            => nlev,            &
                         & config_num_lev         => num_lev,         &
                         & config_nshift          => nshift,          &
                         & config_nsteps          => nsteps,          &
                         & config_dtime           => dtime,           &
                         & config_ltimer          => ltimer,          &
                         & config_timers_level    => timers_level,    &
                         & config_activate_sync_timers => activate_sync_timers, &
                         & config_msg_level       => msg_level,       &
                         & config_logmaxrss       => logmaxrss,       &
                         & config_logmaxrss_all   => logmaxrss_all,   &
                         & config_output          => output,          &
                         & config_output_mode     => output_mode,     &
                         & config_test_mode       => test_mode,       &
                         & t_output_mode, max_output_modes,           &
                         & config_debug_check_level => debug_check_level, &
                         & config_restart_filename  => restart_filename, &
                         & config_profiling_output => profiling_output, &
                         & config_check_uuid_gracefully => check_uuid_gracefully, &
                         & cfg_modelTimeStep => modelTimeStep
#ifdef HAVE_RADARFWO
  USE radar_data_namelist, ONLY: radar_config_radarnmlfile => radarnmlfile
#endif
  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish, set_msg_timestamp
  USE mo_impl_constants, ONLY: max_dom, max_ntracer, inoforcing, IHELDSUAREZ,     &
                               INWP,iaes,ILDF_ECHAM,IMPIOM,INOFORCING,ILDF_DRY, &
                               MAX_CHAR_LENGTH
  USE mo_io_units,       ONLY: nnml, nnml_output
  USE mo_namelist,       ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,            ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_util_string,    ONLY: one_of
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings

  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,   &
       &                            open_and_restore_namelist, close_tmpfile
  USE mtime,                  ONLY: max_timedelta_str_len
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_run_namelist

  CHARACTER(len=*), PARAMETER :: modname = 'mo_run_nml'

CONTAINS
  !>
  !!
  SUBROUTINE read_run_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    CHARACTER(len=*), PARAMETER :: routine = modname//'::read_run_namelist'

    !------------------------------------------------------------------------
    ! Namelist variables
    !------------------------------------------------------------------------

    LOGICAL :: ltestcase       ! if .TRUE. then
    ! - compute analytical initial state,
    !   depending on the specified test case,
    ! - compute analytical boundary conditions,
    ! - if applicable, compute analytical forcing

    LOGICAL :: ldynamics       ! if .TRUE., switch on adiabatic dynamics
    INTEGER :: iforcing        ! adiabatic forcing

    LOGICAL :: ltransport      ! if .TRUE., switch on large-scale tracer transport
    INTEGER :: ntracer         ! number of advected tracers
    !> switch for ICON-ART (Treatment of Aerosols and Trace Gases)
    LOGICAL :: lart

    !> switch for assimilation of radar data using latent heat nudging
    LOGICAL :: ldass_lhn

    LOGICAL :: luse_radarfwo(max_dom)  !< switch for radar forward operator EMVORADO
    CHARACTER(LEN=255) :: radarnmlfile !< name of the file containing the radar namelist
    LOGICAL :: lvert_nest         ! if .TRUE., switch on vertical nesting
    INTEGER :: num_lev(max_dom)   ! number of full levels for each domain
    INTEGER :: nshift (max_dom)   ! half level of parent domain which coincides 
    ! with the upper boundary of the current domain jg

    INTEGER  :: nsteps            ! number of time steps
    REAL(wp) :: dtime             ! [s] length of a time step

    LOGICAL :: ltimer        ! if .TRUE., wallclock timers are switched on
    INTEGER :: timers_level  ! what level of timers to run
    LOGICAL :: activate_sync_timers

    LOGICAL :: logmaxrss     ! log maxrss for three mpi ranks
    LOGICAL :: logmaxrss_all ! log maxrss for all mpi ranks
    
    INTEGER :: msg_level     ! how much printout is generated during runtime
    LOGICAL :: msg_timestamp ! If .TRUE.: Precede output messages by time stamp.


    INTEGER :: test_mode  ! if =0 then run the standard version,
    ! otherwise special setup for (performance) tests, see Namelist_overview
    INTEGER :: debug_check_level

    CHARACTER(len=max_timedelta_str_len) :: modelTimeStep

    !> output mode (logicals)
    !  one or multiple of "none", "nml", "totint"
    CHARACTER(len=32) :: output(max_output_modes)

    INTEGER :: profiling_output  !< switch defining the kind of timer output

    LOGICAL :: check_uuid_gracefully !< Flag. If .TRUE. then we give only warnings for non-matching UUIDs

    !> file name for restart/checkpoint files (containg keyword
    !> substition patterns)
    CHARACTER(len=MAX_CHAR_LENGTH) :: restart_filename

    NAMELIST /run_nml/ ltestcase, ldynamics, iforcing, ltransport,     &
      &                ntracer, lart, ldass_lhn, luse_radarfwo,        &
      &                radarnmlfile, ltimer,                           &
      &                lvert_nest, num_lev, nshift, nsteps, dtime,     &
      &                timers_level, activate_sync_timers, logmaxrss,  &
      &                logmaxrss_all, msg_level, test_mode, output,    &
      &                msg_timestamp, debug_check_level,               &
      &                restart_filename, profiling_output,             &
      &                check_uuid_gracefully, modelTimeStep


    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------

    ltestcase       = .TRUE.
    ldynamics       = .TRUE.
    iforcing        = inoforcing

    ltransport      = .FALSE.
    ntracer         = 0
    lart            = .FALSE.
    ldass_lhn       = .FALSE.

    luse_radarfwo(:) = .FALSE.

    radarnmlfile = ''

    lvert_nest = .FALSE. ! no vertical nesting
    num_lev(:) = 31      ! number of full levels for each domain
    nshift(:)  = 0       ! please do not change the default.
                         ! otherwise the initialization of 
                         ! p_patch(jg)%nshift in "import patches" 
                         ! will not work properly.

    nsteps = -999

    ! Note: The default needs to be empty, since there exist
    ! concurrent namelist parameters to specify these values:
    modelTimeStep = ''
    dtime         = 0._wp
    
    ltimer               = .TRUE.
    timers_level         = 1
    activate_sync_timers = .FALSE.
    logmaxrss            = .FALSE.
    logmaxrss_all        = .FALSE.
    msg_level            = 10
    msg_timestamp        = .FALSE.
    test_mode            = 0
    debug_check_level    = 0

    output(:) = " "
    output(1) = "default"

    restart_filename = "<gridfile>_restart_<mtype>_<rsttime>.<extension>"
    profiling_output = config_profiling_output
    check_uuid_gracefully = .FALSE.

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('run_nml')
      READ(funit,NML=run_nml)
      CALL close_tmpfile(funit)
    END IF

    !----------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !----------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('run_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, run_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, run_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, run_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! Sanity check
    !----------------------------------------------------
    SELECT CASE (iforcing)                                                     
    CASE(INOFORCING,IHELDSUAREZ,INWP,iaes,ILDF_DRY,ILDF_ECHAM,IMPIOM)
    CASE DEFAULT
      CALL finish( TRIM(routine),'wrong value for iforcing')
    END SELECT

    IF ((ntracer<0).OR.(ntracer>max_ntracer)) CALL finish( TRIM(routine), &
      & 'wrong number of tracers. Valid range: 0 <= ntracer <= max_ntracer')

    IF (ANY(num_lev < 0)) CALL finish(TRIM(routine),'"num_lev" must be positive')
    IF (ANY(nshift  < 0)) CALL finish(TRIM(routine),'"nshift" must be positive')

    IF (nsteps < 0 .AND. nsteps /= -999) CALL finish(TRIM(routine),'"nsteps" must not be negative')

    IF (.NOT. ltimer) timers_level = 0

    IF (lart) THEN
#ifndef __ICON_ART
      CALL finish( TRIM(routine),'model set to run with ART but compiled with --disable-art')
#endif
    END IF

    !----------------------------------------------------
    ! Fill part of the configuration state
    !----------------------------------------------------

    config_ltestcase       = ltestcase 
    config_ldynamics       = ldynamics 
    config_iforcing        = iforcing 

    config_ltransport      = ltransport 
    config_ntracer         = ntracer 
    config_lart            = lart
    config_ldass_lhn       = ldass_lhn

    config_luse_radarfwo(:) = luse_radarfwo(:)
    config_radarnmlfile    = radarnmlfile
#ifdef HAVE_RADARFWO
    IF(LEN_TRIM(config_radarnmlfile)>0) radar_config_radarnmlfile = config_radarnmlfile
#endif

    config_lvert_nest      = lvert_nest
    config_nlev            = num_lev(1)
    config_num_lev(:)      = num_lev(:)
    config_nshift(:)       = nshift(:)

    config_nsteps          = nsteps  
    config_dtime           = dtime 

    config_ltimer          = ltimer
    config_timers_level    = timers_level
    config_activate_sync_timers = activate_sync_timers

    config_logmaxrss       = logmaxrss
    config_logmaxrss_all   = logmaxrss_all
    config_msg_level       = msg_level
    CALL set_msg_timestamp(msg_timestamp)
    config_test_mode    = test_mode
    config_debug_check_level = debug_check_level

    config_restart_filename = restart_filename
    config_profiling_output = profiling_output

    config_check_uuid_gracefully = check_uuid_gracefully

    cfg_modelTimeStep       = modelTimeStep

    IF (TRIM(output(1)) /= "default") THEN
      config_output(:) = output(:)
    ELSE
      config_output(:) = " "
      config_output(1) = "nml"
      config_output(2) = "totint"
    END IF
    CALL parse_output_mode(config_output, config_output_mode)
    
    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=run_nml)
      CALL store_and_close_namelist(funit, 'run_nml')
    ENDIF
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=run_nml)

  END SUBROUTINE read_run_namelist
  !-------------


  !> Sets LOGICAL values in "output_mode" according to user input
  !
  SUBROUTINE parse_output_mode(output, om)
    CHARACTER(len=32)   , INTENT(IN)  :: output(max_output_modes)
    TYPE (t_output_mode), INTENT(OUT) :: om
    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = &
      &  modname//'::parse_output_mode'
    CHARACTER(len=8), PARAMETER :: valid_names(4) &
         = (/ "none    ", "nml     ", "totint  ", "maxwinds" /)
    INTEGER, PARAMETER :: midx_none = 1, midx_nml = 2, midx_totint = 3, &
         midx_maxwinds = 4
    INTEGER :: i, match_idx(max_output_modes)

    ! check if user input is in list of valid names:
    DO i=1,max_output_modes
      IF (LEN_TRIM(output(i)) /= 0) THEN
        match_idx(i) = one_of(output(i), valid_names)
        IF  (match_idx(i) == -1) THEN
          CALL finish(routine, "Syntax error: unknown output mode.")
        END IF
      ELSE
        match_idx(i) = -1
      END IF
    END DO

    ! for each logical of type t_output_mode, check if the
    ! corresponding keyword is in the list of strings
    om%l_nml      = ANY(match_idx == midx_nml)
    om%l_totint   = ANY(match_idx == midx_totint)
    om%l_maxwinds = ANY(match_idx == midx_maxwinds)
    om%l_none     = ANY(match_idx == midx_none) &
      .OR. .NOT. (om%l_nml .OR. om%l_totint .OR. om%l_maxwinds)

    ! consistency check:
    ! error: "none" has been chosen in combination with others:
    IF (om%l_none .AND. (om%l_nml .OR. om%l_totint .OR. om%l_maxwinds)) &
      CALL finish(routine, "Syntax error when setting output to 'none'.")

  END SUBROUTINE parse_output_mode

END MODULE mo_run_nml
