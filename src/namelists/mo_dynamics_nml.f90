! Namelist variables shared by the hydrostatic and nonhydrostatic
! dynamical cores
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

MODULE mo_dynamics_nml

  USE mo_dynamics_config,     ONLY: config_iequations     => iequations,     &
                                  & config_divavg_cntrwgt => divavg_cntrwgt, &
                                  & config_lcoriolis      => lcoriolis,      &
                                  & config_lmoist_thdyn   => lmoist_thdyn,   &
                                  & config_ldeepatmo      => ldeepatmo

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_physical_constants,  ONLY: grav
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio 

  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,   &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_dynamics_namelist

  !---------------------------------------------------------------
  ! Namelist variables 
  !---------------------------------------------------------------
  ! time stepping scheme 

  INTEGER  :: iequations

  REAL(wp) :: divavg_cntrwgt ! weight of central cell for divergence averaging

  LOGICAL  :: lcoriolis      ! if .TRUE.,  the Coriolis force is switched on

  LOGICAL  :: lmoist_thdyn   ! if .TRUE., moisture terms included in first law

  LOGICAL  :: ldeepatmo      ! if .TRUE., deep-atmosphere modification is applied 
                             ! to the governing equations, on which the dynamical core is based

  NAMELIST/dynamics_nml/ iequations, divavg_cntrwgt, &
                         lcoriolis, lmoist_thdyn, ldeepatmo

CONTAINS
  !>
  !!
  SUBROUTINE read_dynamics_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_nml:read_dynamics_namelist'

    !------------------------------------------------------------
    ! Set up the default values
    !------------------------------------------------------------
    iequations     = -999
    divavg_cntrwgt = 0.5_wp
    lcoriolis      = .TRUE.
    lmoist_thdyn   = .TRUE.
    ldeepatmo      = .FALSE.

    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by 
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('dynamics_nml')
      READ(funit,NML=dynamics_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('dynamics_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, dynamics_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, dynamics_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, dynamics_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    ! Temporary sanity check, until iequations gets removed completely
    IF (iequations /= -999) THEN
      WRITE(message_text,'(a)') 'WARNING: The Namelist variable iequations is obsolete and will be removed soon.'
      CALL message(routine, message_text)
    ENDIF

    !-----------------------------------------------------
    ! 4. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=dynamics_nml)
      CALL store_and_close_namelist(funit, 'dynamics_nml')
    ENDIF

    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=dynamics_nml)

    !-----------------------------------------------------
    ! 5. Fill configuration state
    !-----------------------------------------------------

    config_divavg_cntrwgt = divavg_cntrwgt
    config_lcoriolis      = lcoriolis
    config_lmoist_thdyn   = lmoist_thdyn
    config_ldeepatmo      = ldeepatmo

  END SUBROUTINE read_dynamics_namelist
  !-------------

END MODULE mo_dynamics_nml
