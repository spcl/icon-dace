! Contains the variables to set up the coupling.
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

MODULE mo_coupling_nml

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------

  USE mo_impl_constants,  ONLY: max_char_length
  USE mo_io_units,        ONLY: nnml
  USE mo_namelist,        ONLY: open_nml, close_nml, position_nml, POSITIONED
  USE mo_exception,       ONLY: message, finish
  USE mo_coupling_config, ONLY: config_coupled_to_ocean, config_coupled_to_waves,    &
    &                           config_coupled_to_hydrodisc, config_coupled_to_atmo, &
    &                           config_coupled_to_output, config_coupled_to_aero,    &
    &                           config_coupled_to_o3
  USE mo_coupling_utils,  ONLY: cpl_config_file_exists
  USE mo_master_control,  ONLY: get_my_process_type, get_my_process_name,           &
    &                           atmo_process, ocean_process, ps_radiation_process,  &
    &                           hamocc_process, jsbach_process, icon_output_process,&
    &                           wave_process, testbed_process

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_coupling_namelist

CONTAINS

  !!  Initialization of variables that contain general information.
  !!
  !!               Initialization of variables that contain general information
  !!               about the coupled model run. The configuration is read from
  !!               namelist 'icon_cpl'.
  !!

  SUBROUTINE read_coupling_namelist (namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    !
    ! Local variables
    !
    LOGICAL :: coupled_to_ocean, can_couple_to_ocean, &
               coupled_to_waves, can_couple_to_waves, &
               coupled_to_atmo, can_couple_to_atmo, &
               coupled_to_hydrodisc, can_couple_to_hydrodisc, &
               coupled_to_output, can_couple_to_output, &
               coupled_to_aero, can_couple_to_aero, &
               coupled_to_o3, can_couple_to_o3

    LOGICAL :: coupled_mode
    INTEGER :: istat

    INTEGER :: my_process_component

    CHARACTER(len=max_char_length), PARAMETER :: &
         &   routine = 'mo_coupling_nml:read_coupling_namelist'

    NAMELIST /coupling_mode_nml/ coupled_to_ocean, coupled_to_waves, coupled_to_atmo, &
         coupled_to_hydrodisc, coupled_to_output, coupled_to_aero, coupled_to_o3

    !--------------------------------------------------------------------
    ! 1. Set default values
    !--------------------------------------------------------------------

    coupled_to_ocean     = .FALSE.
    coupled_to_waves     = .FALSE.
    coupled_to_atmo      = .FALSE.
    coupled_to_hydrodisc = .FALSE.
    coupled_to_output    = .FALSE.
    coupled_to_aero      = .FALSE.
    coupled_to_o3        = .FALSE.

    can_couple_to_ocean     = .FALSE.
    can_couple_to_waves     = .FALSE.
    can_couple_to_atmo      = .FALSE.
    can_couple_to_hydrodisc = .FALSE.
    can_couple_to_output    = .FALSE.
    can_couple_to_aero      = .FALSE.
    can_couple_to_o3        = .FALSE.

    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (done so far by all MPI processes)
    !--------------------------------------------------------------------

#ifdef YAC_coupling

    CALL open_nml (TRIM(namelist_filename))

    CALL position_nml('coupling_mode_nml',STATUS=istat)
    IF (istat==POSITIONED) THEN
      READ (nnml, coupling_mode_nml)
    ENDIF

    CALL close_nml

#endif

    config_coupled_to_ocean     = coupled_to_ocean
    config_coupled_to_waves     = coupled_to_waves
    config_coupled_to_atmo      = coupled_to_atmo
    config_coupled_to_hydrodisc = coupled_to_hydrodisc
    config_coupled_to_output    = coupled_to_output
    config_coupled_to_aero      = coupled_to_aero
    config_coupled_to_o3        = coupled_to_o3

    coupled_mode = ANY((/coupled_to_ocean,     &
                         coupled_to_waves,     &
                         coupled_to_atmo,      &
                         coupled_to_hydrodisc, &
                         coupled_to_output,    &
                         coupled_to_aero,      &
                         coupled_to_o3/))

    !----------------------------------------------------
    ! 3. Sanity checks
    !----------------------------------------------------

#ifndef YAC_coupling

    if (coupled_mode) &
      CALL finish( &
        routine, "(coupled_mode == .TRUE.) " // &
        "but not compiled coupling support")
#endif

    my_process_component = get_my_process_type()

    ! set supported coupling types

    SELECT CASE (my_process_component)
      CASE (atmo_process)
        can_couple_to_ocean = .TRUE.
        can_couple_to_waves = .TRUE.
        can_couple_to_hydrodisc = .TRUE.
        can_couple_to_output = .TRUE.
        can_couple_to_aero = .TRUE.
        can_couple_to_o3 = .TRUE.
      CASE (ocean_process)
        can_couple_to_atmo = .TRUE.
        can_couple_to_hydrodisc = .TRUE.
        can_couple_to_output = .TRUE.
      CASE (hamocc_process)
        can_couple_to_atmo = .TRUE.
      CASE (wave_process)
        can_couple_to_atmo = .TRUE.
      CASE (jsbach_process)
        can_couple_to_ocean = .TRUE.
      CASE (testbed_process)
      CASE (icon_output_process)
      CASE default
        CALL finish(routine, "my_process_component is unsupported")
    END SELECT

    ! checks for unsupported values in coupling_nml
    IF (coupled_to_ocean .AND. .NOT. can_couple_to_ocean) THEN
      CALL finish( &
        routine, 'Component ' // TRIM(get_my_process_name()) // &
        ' does not support coupling to ocean')
    ENDIF

    IF (coupled_to_waves .AND. .NOT. can_couple_to_waves) THEN
      CALL finish( &
        routine, 'Component ' // TRIM(get_my_process_name()) // &
        ' does not support coupling to waves')
    ENDIF

    IF (coupled_to_atmo .AND. .NOT. can_couple_to_atmo) THEN
      CALL finish( &
        routine, 'Component ' // TRIM(get_my_process_name()) // &
        ' does not support coupling to atmo')
    ENDIF

    IF (coupled_to_hydrodisc .AND. .NOT. can_couple_to_hydrodisc) THEN
      CALL finish( &
        routine, 'Component ' // TRIM(get_my_process_name()) // &
        ' does not support coupling to hydrodisc')
    ENDIF

    IF (coupled_to_output .AND. .NOT. can_couple_to_output) THEN
      CALL finish( &
        routine, 'Component ' // TRIM(get_my_process_name()) // &
        ' does not support coupling to output')
    ENDIF

    IF (coupled_to_aero .AND. .NOT. can_couple_to_aero) THEN
      CALL finish( &
        routine, 'Component ' // TRIM(get_my_process_name()) // &
        ' does not support coupling to aero')
    ENDIF

    IF (coupled_to_o3 .AND. .NOT. can_couple_to_o3) THEN
      CALL finish( &
        routine, 'Component ' // TRIM(get_my_process_name()) // &
        ' does not support coupling to o3')
    ENDIF

    IF (coupled_mode .AND. .NOT. cpl_config_file_exists()) THEN
      CALL message( &
        routine, &
        'run is configured to be coupled, but coupler configuration files are not available')
    END IF

  END SUBROUTINE read_coupling_namelist

END MODULE mo_coupling_nml
