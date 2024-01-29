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

MODULE mo_comin_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,   &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_master_config,       ONLY: isRestart
  USE mo_comin_config,        ONLY: comin_config

#ifndef __NO_ICON_COMIN__
  USE comin_host_interface,   ONLY: t_comin_plugin_description
#endif


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_comin_namelist

#ifndef __NO_ICON_COMIN__
  ! Namelist variables
  TYPE(t_comin_plugin_description) :: plugin_list(16) !< list of dynamic libs (max: 16)

  NAMELIST /comin_nml/ plugin_list
#endif


CONTAINS
  !> Read namelist settings for ICON ComIn
  !
  !  Sets defaults for various `comin_nml` user settings, then opens
  !  and reads the namelist file `filename`, overwriting some
  !  defaults.
  !
  SUBROUTINE read_comin_namelist( filename )
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !
#ifndef __NO_ICON_COMIN__
    ! CHARACTER(LEN=*),PARAMETER :: routine='mo_comin_nml:read_comin_namelist'
    INTEGER :: istat, funit, iunit, i

    ! --- set up the default values
    !
    DO i=1,SIZE(plugin_list)
      ! name of primary constructor.
      plugin_list(i)%primary_constructor = "comin_main"

      ! set default of plugin_list(1)%plugin_name to "", such that ICON
      ! expects static linking if not set otherwise.
      plugin_list(i)%plugin_library = ""

      ! name of MPI communicator. left as an empty string if the
      ! application does not require a communicator for this plugin.
      plugin_list(i)%comm = ""

      !  character string containing user-defined options
      !  (e.g. a python script filename) for a specific plugin.
      plugin_list(i)%options = ""
    END DO

    ! --- if this is a resumed integration, overwrite the defaults
    !     above by values in the restart file.
    !
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('comin_nml')
      READ (funit,NML=comin_nml)
      CALL close_tmpfile(funit)
    END IF

    ! --- read user's (new) specifications.
    !
    CALL open_nml(TRIM(filename))
    CALL position_nml ('comin_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE (iunit, comin_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, comin_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE (iunit, comin_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    ! --- store the namelist for restart.
    !
    IF (my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE (funit,NML=comin_nml)
      CALL store_and_close_namelist(funit, 'comin_nml')
    ENDIF

    ! --- write the contents of the namelist to an ASCII file.
    IF (my_process_is_stdio()) WRITE (nnml_output,nml=comin_nml)

    ! --- fill configuration state.
    !
    DO i=1,SIZE(plugin_list)
       IF (TRIM(plugin_list(i)%primary_constructor) == "comin_main" .AND. &
            LEN_TRIM(plugin_list(i)%plugin_library) == 0) EXIT
    END DO
    comin_config%nplugins  = i-1
    comin_config%plugin_list   = plugin_list

    ! --- store the namelist for restart
    !
    IF (my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE (funit,NML=comin_nml)
      CALL store_and_close_namelist(funit, 'comin_nml')
    ENDIF

    ! --- write the contents of the namelist to an ASCII file
    !
    IF (my_process_is_stdio()) WRITE(nnml_output,nml=comin_nml)
#endif
  END SUBROUTINE read_comin_namelist

END MODULE mo_comin_nml
