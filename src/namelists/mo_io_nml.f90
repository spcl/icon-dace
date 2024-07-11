! Contains the setup of the variables for io.
!
!
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

MODULE mo_io_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_char_length, max_ntracer, max_dom, &
    &                              PRES_MSL_METHOD_GME, RH_METHOD_WMO, max_echotop, max_wshear, max_srh
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio, p_n_work
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_nml_and_att,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_io_config,          ONLY: config_lkeep_in_sync           => lkeep_in_sync          , &
                                 & config_dt_diag                 => dt_diag                , &
                                 & config_gust_interval           => gust_interval          , &
                                 & config_ff10m_interval          => ff10m_interval         , &
                                 & config_celltracks_interval     => celltracks_interval    , &
                                 & config_dt_lpi                  => dt_lpi                 , &
                                 & config_dt_hailcast             => dt_hailcast            , &
                                 & config_wdur_min_hailcast       => wdur_min_hailcast      , &
                                 & config_dt_celltracks           => dt_celltracks          , &
                                 & config_dt_radar_dbz            => dt_radar_dbz           , &
                                 & config_echotop_meta            => echotop_meta           , &
                                 & t_echotop_meta                                           , &
                                 & config_precip_interval         => precip_interval        , &
                                 & config_totprec_d_interval      => totprec_d_interval     , &
                                 & config_runoff_interval         => runoff_interval        , &
                                 & config_itype_dursun            => itype_dursun           , &
                                 & config_itype_convindices       => itype_convindices      , &
                                 & config_itype_hzerocl           => itype_hzerocl          , &
                                 & config_sunshine_interval       => sunshine_interval      , &
                                 & config_melt_interval           => melt_interval          , &
                                 & config_maxt_interval           => maxt_interval          , &
                                 & config_dt_checkpoint           => dt_checkpoint          , &
                                 & config_inextra_2d              => inextra_2d             , &
                                 & config_inextra_3d              => inextra_3d             , &
                                 & config_lflux_avg               => lflux_avg              , &
                                 & config_itype_pres_msl          => itype_pres_msl         , &
                                 & config_output_nml_dict         => output_nml_dict        , &
                                 & config_linvert_dict            => linvert_dict           , &
                                 & config_netcdf_dict             => netcdf_dict            , &
                                 & config_lnetcdf_flt64_output    => lnetcdf_flt64_output   , &
                                 & config_itype_rh                => itype_rh               , &
                                 & config_restart_file_type       => restart_file_type      , &
                                 & config_write_initial_state     => write_initial_state    , &
                                 & config_write_last_restart      => write_last_restart     , &
                                 & config_timeSteps_per_outputStep  => timeSteps_per_outputStep, &
                                 & config_lmask_boundary            => lmask_boundary          , &
                                 & config_restart_write_mode        => restart_write_mode   , &
                                 & config_nrestart_streams          => nrestart_streams     , &
                                 & config_checkpoint_on_demand      => checkpoint_on_demand , &
                                 & config_wshear_uv_heights       => wshear_uv_heights      , &
                                 & config_srh_heights             => srh_heights

  USE mo_exception,        ONLY: finish
  USE mo_util_string,      ONLY: tolower
  USE mtime,               ONLY: max_timedelta_str_len


  IMPLICIT NONE
  PUBLIC :: read_io_namelist

  ! module name
  CHARACTER(*), PARAMETER :: modname = "mo_io_nml"

  
CONTAINS
  !! Read Namelist for I/O.
  !!
  !! This subroutine
  !! - reads the Namelist for I/O
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_io_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename

    CHARACTER(*), PARAMETER :: routine = modname//":read_io_namelist"
    INTEGER                        :: istat, funit
    INTEGER                        :: iunit, jg

    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------

    LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
    REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
    REAL(wp):: gust_interval(max_dom)     ! time interval over which maximum wind gusts are taken
    REAL(wp):: ff10m_interval(max_dom)    ! time interval over which ff10m is averaged
    REAL(wp):: celltracks_interval(max_dom)  ! time interval over which extrema of cell track vars are taken
                                             !  (LPI_MAX, UH_MAX, VORW_CTMAX, W_CTMAX, DBZ_CTMAX)
    TYPE(t_echotop_meta) :: echotop_meta(max_dom) ! meta data for echotops (ECHOTOP, ECHOTOPinM)
    CHARACTER(len=max_timedelta_str_len) :: precip_interval(max_dom)   ! time interval over which precipitation variables are accumulated
    CHARACTER(len=max_timedelta_str_len) :: totprec_d_interval(max_dom)! time interval over which the special tot_prec_d is accumulated
    CHARACTER(len=max_timedelta_str_len) :: runoff_interval(max_dom)   ! time interval over which runoff variables are accumulated
    CHARACTER(len=max_timedelta_str_len) :: sunshine_interval(max_dom) ! time interval over which sunshine duration is accumulated
    CHARACTER(len=max_timedelta_str_len) :: melt_interval(max_dom)     ! time interval over which snow melt is accumulated
    CHARACTER(len=max_timedelta_str_len) :: maxt_interval(max_dom)     ! time interval for tmax_2m and tmin_2m 
    REAL(wp):: dt_lpi                     ! calling frequency [seconds] of lpi diagnosis for hourly maximum calculation
    REAL(wp):: dt_hailcast                ! calling frequency [seconds] of hail diagnosis for hourly maximum calculation
    REAL(wp):: wdur_min_hailcast          ! minimal updraft persistence [seconds] for hailcast to be activated
    REAL(wp):: dt_celltracks              ! calling frequency [seconds] of celltrack diagnosis for hourly maximum calculation
                                          ! this pertains to the following variables: tcond_max/tcond10_max, uh_max, vorw_ctmax, w_ctmax
    REAL(wp):: dt_radar_dbz               ! calling frequency [seconds] of radar reflectivity diagnosis for hourly maximum calculation

    REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

    INTEGER :: inextra_2d                 ! number of extra output fields for debugging
    INTEGER :: inextra_3d                 ! number of extra output fields for debugging
    INTEGER :: itype_dursun               ! if 0 the sunshine duration is counted if >120W/m^2
                                          ! if 1 the sunshine duration is counted only
                                          !    if direct radiation > 200 W/m^2 and relative sunshine duration in % is computed
    INTEGER :: itype_convindices          ! if 1 CAPE_MU/CIN_MU are approximated via the CAPE/CIN of the parcel with maximum equivalent temperature
                                          ! if 2 the full computation is done
    INTEGER :: itype_hzerocl              ! Specifies height of freezing level if T < 0 Celsius in the whole atmospheric column
                                          ! 1: set hzerocl to orography height (default)
                                          ! 2: set hzerocl to -999.0_wp (undef)
                                          ! 3: set hzerocl to extrapolated value below ground (assuming -6.5 K/km)
    LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated
                                          !  from the beginning of the run
                                          ! if .TRUE. the output fluxex are average values
                                          !  from the beginning of the run, except of
                                          !  TOT_PREC that would be accumulated
    INTEGER :: itype_pres_msl             ! Specifies method for computation of mean sea level pressure
                                          ! 1: GME-type extrapolation
                                          ! 2: stepwise analytical integration
                                          ! 3: IFS method
                                          ! 4: IFS method with consistency correction
    INTEGER :: itype_rh                   ! Specifies method for computation of relative humidity
                                          ! 1: WMO: water only (e_s=e_s_water)
                                          ! 2: IFS: mixed phases (e_s=a*e_s_water + b*e_s_ice)


    CHARACTER(LEN=filename_max) :: &
      &        output_nml_dict,    &      !< maps variable names onto the internal ICON names.
      &        netcdf_dict                !< maps internal variable names onto names in output file (NetCDF only).

    LOGICAL :: linvert_dict               !< inverts columns in output_nml_dict (allows using the same dictionary file as for input)

    LOGICAL :: lnetcdf_flt64_output       !< if .TRUE. floating point valued NetCDF output
                                          !  is written in 64-bit instead of 32-bit accuracy


    INTEGER :: restart_file_type

    LOGICAL :: write_initial_state

    LOGICAL :: write_last_restart

    INTEGER :: timeSteps_per_outputStep

    LOGICAL :: lmask_boundary(max_dom) ! flag: true, if interpolation zone should be masked *in output*

    CHARACTER(LEN = 256) :: restart_write_mode

    ! When using the restart write mode "dedicated proc mode", it is
    ! possible to split the restart output into several files, as if
    ! "nrestart_streams" * "num_io_procs" restart processes were
    ! involved. This speeds up the read-in process, since all the
    ! files may then be read in parallel.
    INTEGER :: nrestart_streams

    ! Allows checkpointing (followed by stopping) during runtime triggered by a file named 'stop_icon' in the workdir
    LOGICAL :: checkpoint_on_demand

    REAL(wp) :: wshear_uv_heights(1:max_wshear)

    REAL(wp) :: srh_heights(1:max_srh)

    NAMELIST/io_nml/ lkeep_in_sync, dt_diag, dt_checkpoint,               &
      &              inextra_2d, inextra_3d,                              &
      &              lflux_avg, itype_pres_msl, itype_rh,                 &
      &              output_nml_dict, netcdf_dict, linvert_dict,          &
      &              lnetcdf_flt64_output,                                &
      &              restart_file_type, write_initial_state,              &
      &              write_last_restart, timeSteps_per_outputStep,        &
      &              lmask_boundary, gust_interval, restart_write_mode,   &
      &              nrestart_streams, celltracks_interval, echotop_meta, &
      &              precip_interval, totprec_d_interval, runoff_interval,&
      &              maxt_interval, checkpoint_on_demand,                 &
      &              nrestart_streams, dt_lpi, dt_celltracks,             &
      &              dt_hailcast, wdur_min_hailcast,                      &
      &              dt_radar_dbz, sunshine_interval, itype_dursun,       &
      &              itype_convindices, itype_hzerocl, melt_interval,     &
      &              wshear_uv_heights, srh_heights, ff10m_interval

    !-----------------------
    ! 1. default settings
    !-----------------------
    lkeep_in_sync           = .FALSE.

    dt_diag                 = 86400._wp    !  1 day

    ! Note: The default needs to be empty, since there exists a
    ! concurrent namelist parameter to specify this value:
    dt_checkpoint           = 0._wp  ! unspecified

    gust_interval(:)        = 3600._wp     ! 1 hour
    ff10m_interval(:)       = 600._wp      ! 10 min
    celltracks_interval(:)  = 3600._wp     ! 1 hour
    DO jg=1, max_dom 
      ! echotop_meta(jg)%nechotop will be re-computed later in mo_nml_crosscheck.f90
      echotop_meta(jg)%nechotop                  = 0
      echotop_meta(jg)%time_interval             = 3600._wp     ! 1 hour
      echotop_meta(jg)%dbzthresh(1:max_echotop)  = -999.99_wp   ! missing value
    END DO
    precip_interval(:)      = "P01Y"       ! 1 year
    totprec_d_interval(:)   = "PT01H"      ! 1 hour
    runoff_interval(:)      = "P01Y"       ! 1 year
    sunshine_interval(:)    = "P01Y"       ! 1 year
    melt_interval(:)        = "P01Y"       ! 1 year
    maxt_interval(:)        = "PT06H"      ! 6 hours
    dt_lpi                  = 180._wp      ! 3 minutes
    dt_celltracks           = 120._wp      ! 2 minutes
    dt_radar_dbz            = 120._wp      ! 2 minutes
    dt_hailcast             = 180._wp      ! 3 minutes
    wdur_min_hailcast       = 900._wp      ! 15 minutes
    inextra_2d              = 0     ! no extra output 2D fields
    inextra_3d              = 0     ! no extra output 3D fields
    itype_dursun            = 0
    itype_convindices       = 1
    itype_hzerocl           = 1
    lflux_avg               = .TRUE.
    itype_pres_msl          = PRES_MSL_METHOD_GME
    itype_rh                = RH_METHOD_WMO       ! WMO: water only
    output_nml_dict         = ' '
    netcdf_dict             = ' '
    linvert_dict            = .FALSE.
    lnetcdf_flt64_output    = .FALSE.

    restart_file_type       = config_restart_file_type
    write_initial_state     = config_write_initial_state
    write_last_restart      = config_write_last_restart
    timeSteps_per_outputStep        = config_timeSteps_per_outputStep

    lmask_boundary(:)       = .FALSE.

    restart_write_mode = ""
    nrestart_streams   = 1
    checkpoint_on_demand = .FALSE.

    wshear_uv_heights(:) = -999.99_wp  ! missing value
    wshear_uv_heights(1:3) = (/ 1000.0_wp, 3000.0_wp, 6000.0_wp /)
    
    srh_heights(:) = -999.99_wp  ! missing value
    srh_heights(1:2) = (/ 1000.0_wp, 3000.0_wp /)
    

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('io_nml')
      READ(funit,NML=io_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('io_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, io_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, io_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, io_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill echotop_meta with defaults, if nothing
    !    has been specified in the namelist
    !----------------------------------------------------

    DO jg=1, max_dom
      ! echotop_meta(jg)%nechotop will be computed later in mo_nml_crosscheck.f90
      IF ( ALL(echotop_meta(jg)%dbzthresh(1:max_echotop) < -900.0_wp) ) THEN
        echotop_meta(jg)%dbzthresh(1)              = 18.0_wp
        echotop_meta(jg)%dbzthresh(2:max_echotop)  = -999.99_wp ! missing value
      END IF
    END DO

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_lkeep_in_sync           = lkeep_in_sync
    config_dt_diag                 = dt_diag
    config_gust_interval(:)        = gust_interval(:)
    config_ff10m_interval(:)       = ff10m_interval(:)
    config_celltracks_interval(:)  = celltracks_interval(:)
    config_echotop_meta(:)         = echotop_meta(:)
    config_precip_interval(:)      = precip_interval(:)
    config_totprec_d_interval(:)   = totprec_d_interval(:)
    config_runoff_interval(:)      = runoff_interval(:)
    config_itype_dursun            = itype_dursun
    config_itype_convindices       = itype_convindices
    config_itype_hzerocl           = itype_hzerocl
    config_sunshine_interval(:)    = sunshine_interval(:)
    config_melt_interval(:)        = melt_interval(:)
    config_maxt_interval(:)        = maxt_interval(:)
    config_dt_checkpoint           = dt_checkpoint
    config_dt_lpi                  = dt_lpi
    config_dt_hailcast             = dt_hailcast
    config_wdur_min_hailcast       = wdur_min_hailcast
    config_dt_celltracks           = dt_celltracks
    config_dt_radar_dbz            = dt_radar_dbz
    config_inextra_2d              = inextra_2d
    config_inextra_3d              = inextra_3d
    config_lflux_avg               = lflux_avg
    config_itype_pres_msl          = itype_pres_msl
    config_itype_rh                = itype_rh
    config_output_nml_dict         = output_nml_dict
    config_netcdf_dict             = netcdf_dict
    config_linvert_dict            = linvert_dict
    config_lnetcdf_flt64_output    = lnetcdf_flt64_output
    config_restart_file_type       = restart_file_type
    config_write_initial_state     = write_initial_state
    config_timeSteps_per_outputStep= timeSteps_per_outputStep
    config_write_last_restart      = write_last_restart
    config_lmask_boundary          = lmask_boundary
    config_restart_write_mode      = tolower(restart_write_mode)
    config_nrestart_streams        = nrestart_streams
    config_checkpoint_on_demand    = checkpoint_on_demand
    config_wshear_uv_heights       = wshear_uv_heights
    config_srh_heights             = srh_heights
    
    ! --- consistency check:

    ! Each work can send its data only to one restart PE. Therefore it
    ! is not possible to have more restart files than source PEs.
    IF ((nrestart_streams < 0) .OR. (nrestart_streams > p_n_work)) THEN
      CALL finish(routine, "Invalid choice of parameter value: nrestart_streams!")
    END IF
    IF((itype_convindices < 0) .OR. (itype_convindices > 2)) THEN
          CALL finish(routine, 'itype_convindices can only have the value 0 (no MU computation), or 1 or 2 when requesting computation.')
    ENDIF

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=io_nml)
      CALL store_and_close_namelist(funit, 'io_nml')
    ENDIF
    !-----------------------------------------------------
    ! 6. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=io_nml)
    END IF


  END SUBROUTINE read_io_namelist

END MODULE mo_io_nml
