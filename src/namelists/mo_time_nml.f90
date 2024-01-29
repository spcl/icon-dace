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

MODULE mo_time_nml

  USE mo_kind,                  ONLY: wp
  USE mo_time_config,           ONLY: cfg_dt_restart          => dt_restart,          &
    &                                 cfg_icalendar           => icalendar,           &
    &                                 cfg_is_relative_time    => is_relative_time,    &
    &                                 cfg_ini_datetime_string => ini_datetime_string, &
    &                                 cfg_end_datetime_string => end_datetime_string, &
    &                                 restart_ini_datetime_string,                    &
    &                                 restart_end_datetime_string,                    &
    &                                 restart_calendar
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_master_control,        ONLY: isRestart
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                   ONLY: my_process_is_stdio 
  USE mo_restart_nml_and_att,   ONLY: open_and_restore_namelist, close_tmpfile, &
                                    & open_tmpfile, store_and_close_namelist
  USE mo_nml_annotate,          ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_time_namelist

  !--------------------
  ! Namelist variables
  !--------------------

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_time_nml"

CONTAINS

  !-------------------------------------------------------------------------
  !!  Initialization of variables that contain general information.
  !!
  !!  Initialization of variables that contain general information
  !!  about the model run. The configuration is read from
  !!  namelist 'time_nml'.
  !!
  SUBROUTINE read_time_namelist( filename )

   CHARACTER(LEN=*), INTENT(IN)        :: filename
   ! local variables
   CHARACTER(len=*), PARAMETER         :: routine = modname//'::read_time_namelist'

   INTEGER            ::  calendar         ! calendar type
   REAL(wp)           ::  dt_restart       ! [s] length of a restart cycle 

   CHARACTER(len=32)  ::  ini_datetime_string, end_datetime_string

   !> LOGICAL is_relative_time: .TRUE., if time loop shall start with
   !> step 0 regardless whether we are in a standard run or in a
   !> restarted run (which means re-initialized run):
   LOGICAL            ::  is_relative_time

   NAMELIST /time_nml/ calendar,            &
     &                 ini_datetime_string, &
     &                 end_datetime_string, &
     &                 dt_restart,          &
     &                 is_relative_time

   INTEGER                             :: istat, funit
   INTEGER                             :: iunit

   !------------------------------------------------------------------------
   ! Default values
   !------------------------------------------------------------------------

   ! Note: The default needs to be empty, since there exist
   ! concurrent namelist parameters to specify these values:
   calendar  = -1   ! unspecified

   ! Initialize start and end date by empty strings, which means
   ! "undefined". The same information may be set through other
   ! namelist parameters from "master_time_control_nml".
   ini_datetime_string = ""     ! initial date and time
   end_datetime_string = ""     ! end date and time
   dt_restart          = 0._wp  ! length of restart cycle (unspecified)

   is_relative_time = .FALSE.

   IF (isRestart()) THEN

     ! 2.1 Overwrite the defaults above by values in the restart file
     funit = open_and_restore_namelist('time_nml')
     READ(funit,NML=time_nml)
     CALL close_tmpfile(funit) 

     ! store the namelist settings originating from the restart file:
     restart_calendar            = calendar
     restart_ini_datetime_string = ini_datetime_string
     restart_end_datetime_string = end_datetime_string
     
   END IF
    
   !------------------------------------------------------------------------
   !  Read user's (new) specifications. (Done so far by all MPI processes)
   !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('time_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, time_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, time_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, time_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    cfg_icalendar = calendar

    cfg_ini_datetime_string = ini_datetime_string
    cfg_end_datetime_string = end_datetime_string

    cfg_dt_restart          = dt_restart
    cfg_is_relative_time    = is_relative_time

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=time_nml)
      CALL store_and_close_namelist(funit, 'time_nml')
    ENDIF
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=time_nml)

 END SUBROUTINE read_time_namelist

END MODULE mo_time_nml
