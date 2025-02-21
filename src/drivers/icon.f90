! This is the master program of the ICON model.
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
PROGRAM icon

  USE mo_exception,           ONLY: message_text, message, finish, enable_logging
  USE mo_io_units,            ONLY: filename_max
  USE mo_mpi,                 ONLY: start_mpi , stop_mpi, my_process_is_global_root,    &
    &                               my_process_is_stdio
  USE mo_master_init,         ONLY: init_master_control
  USE mo_master_control,      ONLY: get_my_namelist_filename, get_my_process_type,      &
    &                               atmo_process, ocean_process, ps_radiation_process,  &
    &                               hamocc_process, jsbach_process, icon_output_process,&
    &                               wave_process
  USE mo_master_control,      ONLY: testbed_process
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: OPERATOR(>)
  USE mo_util_vcs,            ONLY: show_version

  USE mo_ocean_model,         ONLY: ocean_model
  USE mo_hamocc_model,        ONLY: hamocc_model

  USE mo_wave_model,          ONLY: wave_model

  USE mo_icon_testbed,        ONLY: icon_testbed

  USE mo_atmo_model,          ONLY: atmo_model

  USE mo_jsbach_model,        ONLY: jsbach_model

  USE mo_icon_output_model, ONLY: icon_output_driver



  IMPLICIT NONE

  INTEGER                     :: master_control_status, my_process_component
  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=filename_max) :: master_namelist_filename="icon_master.namelist"



  ! handling of comand-line arguments:
  TYPE t_cmdarg_option
    CHARACTER(len=1024) :: arg   !< (case-sensitive) option
    CHARACTER(len=1024) :: help  !< help string
  END TYPE t_cmdarg_option

  ENUM, BIND(C)
    ENUMERATOR :: ARG_UNKNOWN = 0, &
      &           ARG_HELP,        &
      &           ARG_VERSION
  END ENUM

  TYPE(t_cmdarg_option), PARAMETER :: cmdarg_options(2) =            &
    & [                                                              &
    &  t_cmdarg_option("--help",    "print this help message."),     &
    &  t_cmdarg_option("--version", "print version info and exit.")  &
    & ]

  INTEGER :: i,j
  LOGICAL :: lmatch, lmatch_ij, lcmdarg(0:SIZE(cmdarg_options))
  CHARACTER(len=1024) :: arg


!--------------------------------------------------------------------


  !-------------------------------------------------------------------
  ! Initialize MPI, this should always be the first call
  CALL start_mpi('ICON')

  !-------------------------------------------------------------------
  !set up signal trapping on IBM: export USE_SIGNAL_HANDLING=yes


  ! print info on the current version:
  CALL show_version()

  ! When executing ICON, it is now possible to invoke command-line
  ! arguments (logical switches).
  i = 1
  lcmdarg(:) = .FALSE.
  DO
    CALL get_command_argument(i, arg)
    IF (LEN_TRIM(arg) == 0) EXIT

    lmatch = .FALSE.
    DO j=1,SIZE(cmdarg_options)
      lmatch_ij  = (TRIM(arg) == TRIM(cmdarg_options(j)%arg))
      lcmdarg(j) = lcmdarg(j) .OR. lmatch_ij
      lmatch     = lmatch     .OR. lmatch_ij
    END DO
    IF (.NOT. lmatch) THEN
      lcmdarg(ARG_UNKNOWN) = .TRUE.
      CALL message("", "command-line argument '"//TRIM(arg)//"' unknown!")
    END IF
    i = i+1
  END DO

  ! print a list of available options and exit
  IF (ANY([lcmdarg(ARG_UNKNOWN), lcmdarg(ARG_HELP)])) THEN
    CALL message("", "")
    CALL message("", "list of available command-line options:")
    DO j=1,SIZE(cmdarg_options)
      WRITE(message_text,'(a,A20,a)') "    ", [CHARACTER(50) :: TRIM(cmdarg_options(j)%arg)], TRIM(cmdarg_options(j)%help)
      CALL message('', message_text)
    END DO
    CALL message("", "")
  END IF

  ! for '--version' ICON prints the same information as usually at the
  ! beginning of stdout and aborts then.
  IF (ANY([lcmdarg(ARG_UNKNOWN), lcmdarg(ARG_HELP), lcmdarg(ARG_VERSION)])) THEN
    CALL stop_mpi ! Shut down MPI
    STOP 'icon'
  END IF


  !-------------------------------------------------------------------
  ! Initialize the master control

  master_control_status = init_master_control(TRIM(master_namelist_filename))



  my_namelist_filename = get_my_namelist_filename()
  my_process_component = get_my_process_type()

  CALL enable_logging(my_process_is_stdio())

  SELECT CASE (my_process_component)

  CASE (atmo_process)
    CALL atmo_model  (my_namelist_filename, TRIM(master_namelist_filename))

  CASE (ocean_process)
    CALL ocean_model (my_namelist_filename, TRIM(master_namelist_filename))

  CASE (hamocc_process)
    CALL hamocc_model (my_namelist_filename, TRIM(master_namelist_filename))

  CASE (wave_process)
    CALL wave_model (my_namelist_filename, TRIM(master_namelist_filename))

  CASE (jsbach_process)
    CALL jsbach_model (my_namelist_filename, TRIM(master_namelist_filename))

  CASE (testbed_process)
    CALL icon_testbed(my_namelist_filename, TRIM(master_namelist_filename))

  CASE (icon_output_process)
    CALL icon_output_driver(my_namelist_filename, TRIM(master_namelist_filename))

  CASE default
    CALL finish("icon","my_process_component is unknown")

  END SELECT

  IF (ASSOCIATED(time_config%tc_exp_stopdate) .AND. ASSOCIATED(time_config%tc_stopdate)) THEN
    ! write the control.status file
    IF (my_process_is_global_root()) THEN
      OPEN (500, FILE="finish.status")
      IF ((time_config%tc_exp_stopdate > time_config%tc_stopdate) .AND. time_config%tc_write_restart) THEN
        WRITE(500,*) "RESTART"
      ELSE
        WRITE(500,*) "OK"
      ENDIF
      CLOSE(500)
    END IF
  END IF

  ! Shut down MPI
  CALL stop_mpi


END PROGRAM icon
