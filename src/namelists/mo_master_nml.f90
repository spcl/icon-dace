! Master namelist.
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

MODULE mo_master_nml

  USE mo_exception,      ONLY: finish, message, warning
  USE mo_io_units,       ONLY: filename_max, nnml
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED
  USE mo_util_string,    ONLY: t_keyword_list, associate_keyword, with_keywords, tolower
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings
  USE mo_mpi,            ONLY: my_process_is_stdio
  USE mtime,             ONLY: max_calendar_str_len, setCalendar,                            &
       &                       proleptic_gregorian, year_of_365_days,  year_of_360_days,     &
       &                       max_datetime_str_len, max_timedelta_str_len,                  &
       &                       datetime, newDatetime, deallocateDatetime,                    &
       &                       timedelta, newTimedelta, deallocateTimedelta,                 &
       &                       datetimeToString, OPERATOR(+), register_print_mtime_procedure 
  USE mo_master_config,  ONLY: master_component_models, addModel, noOfModels, maxNoOfModels, &
       &                       setInstitution, setRestart, setReadRestartNamelists,          &
       &                       setRestartWriteLast, setModelBaseDir,                         &
       &                       cfg_experimentReferenceDate => experimentReferenceDate,       &
       &                       cfg_experimentStartDate     => experimentStartDate,           &
       &                       cfg_experimentStopDate      => experimentStopDate,            &
       &                       cfg_calendar                => calendar_str,                  &
       &                       cfg_checkpointTimeIntval => checkpointTimeIntval,             &
       &                       cfg_restartTimeIntval => restartTimeIntval

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: read_master_namelist
    
CONTAINS

  !>
  !! Initialization of variables that contain general information
  !! about the coupled model run. The configuration is read from
  !! namelist 'master_nml'.
  !!
  INTEGER FUNCTION read_master_namelist(namelist_filename)
    
    CHARACTER(len=*), INTENT(in) :: namelist_filename

    ! Local variables
    
    ! Namelist variables

    CHARACTER(len=256) :: institute = ''

    !> Flag: True, if model run is initialized from restart state.
    LOGICAL :: lRestart             = .FALSE.
    !> Flag: True, if model run should create restart at experiment end.
    !  This is independent from the settings of the restart interval.
    LOGICAL :: lrestart_write_last  = .TRUE.

    ! if is restart, the flag indecates if we should read the restart namelists
    LOGICAL :: read_restart_namelists  = .TRUE.

    CHARACTER(len=filename_max) :: model_base_dir          = ''
    CHARACTER(len=132)          :: model_name              = ''
    CHARACTER(len=filename_max) :: model_namelist_filename = ''
    
    INTEGER :: model_type 
    INTEGER :: model_min_rank
    INTEGER :: model_max_rank
    INTEGER :: model_inc_rank
    INTEGER :: model_rank_group_size
    
    CHARACTER(len=max_calendar_str_len) :: calendar                 = ''
    CHARACTER(len=max_datetime_str_len) :: experimentReferenceDate  = ''   
    CHARACTER(len=max_datetime_str_len) :: experimentStartDate      = ''
    CHARACTER(len=max_datetime_str_len) :: experimentStopDate       = ''

    CHARACTER(len=max_timedelta_str_len) :: forecastLeadTime        = ''
    
    CHARACTER(len=max_timedelta_str_len) :: checkpointTimeIntval    = ''
    CHARACTER(len=max_timedelta_str_len) :: restartTimeIntval       = ''

    TYPE(datetime), POINTER :: experiment_start_date, experiment_stop_date
    TYPE(timedelta), POINTER :: forecast_lead_time
    
    NAMELIST /master_nml/              &
         &    institute,               &
         &    lRestart,                &
         &    lrestart_write_last,     &
         &    model_base_dir,          &
         &    read_restart_namelists  
    
    NAMELIST /master_time_control_nml/ &
         &    calendar,                &
         &    experimentReferenceDate, &   
         &    experimentStartDate,     &
         &    experimentStopDate,      &
         &    forecastLeadTime,        &
         &    checkpointTimeIntval,    &
         &    restartTimeIntval        
    
    NAMELIST /master_model_nml/        &
         &    model_name,              &
         &    model_namelist_filename, &
         &    model_type,              &
         &    model_min_rank,          &
         &    model_max_rank,          &
         &    model_inc_rank,          &
         &    model_rank_group_size

    INTEGER :: icalendar
    INTEGER :: istat
    INTEGER :: iunit
    LOGICAL :: lrewind
    
    CHARACTER(len=*), PARAMETER :: routine = 'mo_master_nml:read_master_namelist'
    
    TYPE (t_keyword_list), POINTER :: keywords         => NULL()
    
    ! Read  master_nml (done so far by all MPI processes)
    
    istat = 0
    CALL open_nml(namelist_filename, lwarn=.TRUE., istat=istat)
    IF (istat /= 0) THEN
      read_master_namelist = -1
      RETURN
    ENDIF

    ! --------------------------------------------------------------------------------
    ! MASTER_NML
    ! --------------------------------------------------------------------------------
    
    CALL position_nml('master_nml', STATUS=istat)

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, master_nml)  ! write defaults to temporary text file
    END IF

    IF (istat == POSITIONED) THEN
      READ (nnml, master_nml)
    ENDIF

    IF (my_process_is_stdio()) THEN
      iunit = temp_settings()
      WRITE(iunit, master_nml)  ! write settings to temporary text file
    END IF

    SELECT CASE (institute)
    CASE ('DWD')
      CALL setInstitution('Deutscher Wetterdienst')
    CASE ('MPIM')  
      CALL setInstitution('Max Planck Institute for Meteorology')
    CASE ('KIT')  
      CALL setInstitution('Karlsruhe Institute of Technology')
    CASE ('CSCS')
      CALL setInstitution('Swiss National Supercomputing Centre')      
    CASE DEFAULT
      CALL setInstitution('Max Planck Institute for Meteorology/Deutscher Wetterdienst')
    END SELECT
    
    ! save namelist variables in configuration

    CALL setRestart(lRestart)
    CALL setRestartWriteLast(lrestart_write_last)
    CALL setModelBaseDir(model_base_dir)
    CALL setReadRestartNamelists(read_restart_namelists)
   

    ! --------------------------------------------------------------------------------
    ! MASTER_TIME_CONTROL_NML
    ! --------------------------------------------------------------------------------

    ! Note: The default needs to be empty, since there exist
    ! concurrent namelist parameters to specify these values:
    calendar                = ""
    experimentReferenceDate = ""
    experimentStartDate     = ""
    experimentStopDate      = ""
    forecastLeadTime        = ""
    checkpointTimeIntval    = ""
    restartTimeIntval       = ""

    CALL position_nml('master_time_control_nml', STATUS=istat)

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, master_time_control_nml)  ! write defaults to temporary text file
    END IF

    IF (istat == POSITIONED) THEN
      READ (nnml, master_time_control_nml)
    ENDIF
    
    IF (my_process_is_stdio()) THEN
      iunit = temp_settings()
      WRITE(iunit, master_time_control_nml)  ! write settings to temporary text file
    END IF

    cfg_experimentReferenceDate = experimentReferenceDate
    cfg_experimentStartDate     = experimentStartDate
    cfg_experimentStopDate      = experimentStopDate
    cfg_checkpointTimeIntval    = checkpointTimeIntval
    cfg_restartTimeIntval       = restartTimeIntval
    cfg_calendar                = calendar

    ! set calendar (singleton, so do not change later!)
    
    SELECT CASE (toLower(calendar))
    CASE ('proleptic gregorian')
      icalendar  = proleptic_gregorian
    CASE ('365 day year')  
      icalendar = year_of_365_days
    CASE ('360 day year')  
      icalendar = year_of_360_days
    CASE default
      icalendar  = proleptic_gregorian
      CALL message('','No calendar selected! Use default proleptic Gregorian.')
    END SELECT

    CALL setCalendar(icalendar)
    CALL register_print_mtime_procedure(warning)
    
    IF (experimentStartDate /= "") THEN
      IF (experimentStopDate == "") THEN
        IF (forecastLeadTime /= "") THEN
          forecast_lead_time => newTimedelta(TRIM(forecastLeadTime))
          experiment_start_date => newDatetime(TRIM(experimentStartDate))
          experiment_stop_date => newDatetime(TRIM(experimentStartDate))
          experiment_stop_date = experiment_start_date + forecast_lead_time
          CALL datetimeToString(experiment_stop_date, experimentStopDate)
          cfg_experimentStopDate = experimentStopDate
          CALL deallocateDatetime(experiment_stop_date)
          CALL deallocateDatetime(experiment_start_date)
          CALL deallocateTimedelta(forecast_lead_time)      
        ELSE
          CALL finish('','Need forecastLeadTime AND experimentStartDate set in master_time_control_nml namelist.')
        ENDIF
      ENDIF
    ENDIF


    ! --------------------------------------------------------------------------------
    ! MASTER_MODEL_NML
    ! --------------------------------------------------------------------------------
    
    ! for positioning to the first entry of the namelist 
    lrewind = .TRUE.
    DO
      CALL position_nml('master_model_nml', lrewind=lrewind, status=istat)
      IF ( istat /= POSITIONED ) EXIT

      IF (noOfModels() >= maxNoOfModels) THEN
        CALL finish(routine, 'no of models >= max no of models')
      ENDIF
      
      ! change to be able to allow fetching the next namelist entry
      lrewind = .FALSE.
      
      ! default values

      model_name              = ''
      model_namelist_filename = ''
      model_type              = -1
      model_min_rank          =  0
      model_max_rank          = -1 
      model_inc_rank          = 1
      model_rank_group_size   = 1
      
      IF (my_process_is_stdio()) THEN
        iunit = temp_defaults()
        WRITE(iunit, master_model_nml)  ! write defaults to temporary text file
      END IF
      
      READ (nnml, master_model_nml)     ! overwrite default settings

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, master_model_nml)  ! write settings to temporary text file
      END IF

      CALL addModel()
      master_component_models(noOfModels())%model_name = model_name

      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)

      master_component_models(noOfModels())%model_namelist_filename = &
        &  TRIM(with_keywords(keywords, model_namelist_filename))

      master_component_models(noOfModels())%model_type = model_type

      master_component_models(noOfModels())%model_min_rank = model_min_rank
      master_component_models(noOfModels())%model_max_rank = model_max_rank
      master_component_models(noOfModels())%model_inc_rank = model_inc_rank
      master_component_models(noOfModels())%model_rank_group_size = model_rank_group_size
     
    ENDDO
      
    CLOSE (nnml, IOSTAT=istat)

    read_master_namelist = 0

  END FUNCTION read_master_namelist

END MODULE mo_master_nml
