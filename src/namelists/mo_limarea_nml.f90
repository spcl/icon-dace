! Contains the setup of variables related to the boundary relaxation of limited
! area models from an external dataset
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

MODULE mo_limarea_nml

  USE mo_kind,                ONLY: wp, i8
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist     , &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mtime,                  ONLY: MAX_TIMEDELTA_STR_LEN, newTimedelta, getPTStringFromMS

  IMPLICIT NONE
  PRIVATE

  PUBLIC read_limarea_namelist

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_limarea_nml'

CONTAINS
  !>
  !!
  SUBROUTINE read_limarea_namelist( filename )
    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::configure_latbc"

    INTEGER                              :: istat, funit
    INTEGER                              :: iunit, errno
    REAL(wp)                             :: dtime_latbc_in_ms
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: dtime_latbc_str

    !------------------------------------------------------------------------
    ! Namelist variables
    !------------------------------------------------------------------------
    !> type of limited area boundary nudging
    INTEGER                         :: itype_latbc
    !> dt between two consequtive external latbc files
    REAL(wp)                        :: dtime_latbc
    !> number of vertical levels in boundary data
    INTEGER                         :: nlev_latbc
    !> prefix of latbc files
    CHARACTER(LEN=filename_max)     :: latbc_filename
    !> directory containing external latbc files
    !! FIXME: since latbc_path refers to a path, using max_char_length
    !! seems an odd choice
    CHARACTER(LEN=MAX_CHAR_LENGTH)  :: latbc_path
    !> grid file defining the lateral boundary
    CHARACTER(LEN=FILENAME_MAX)     :: latbc_boundary_grid
    !> if set to TRUE, qi and qc are read from latbc data
    LOGICAL                         :: latbc_contains_qcqi
    !> take initial lateral boundary conditions from first guess
    LOGICAL                         :: init_latbc_from_fg
    !> use hydrostatic pressure for lateral boundary nudging
    LOGICAL                         :: nudge_hydro_pres
    !> factor for pressure bias correction of latbc data
    REAL(wp)                        :: fac_latbc_presbiascor

    ! dictionary which maps internal variable names onto
    ! GRIB2 shortnames or NetCDF var names used for lateral boundary nudging.
    CHARACTER(LEN=filename_max) :: latbc_varnames_map_file

    !> if LatBC data is unavailable: number of retries
    INTEGER                         :: nretries

    !> if LatBC data is unavailable: idle wait seconds between retries
    INTEGER                         :: retry_wait_sec


    NAMELIST /limarea_nml/ itype_latbc, dtime_latbc, nlev_latbc,                         &
      &                     latbc_filename, latbc_path, latbc_boundary_grid,             &
      &                     latbc_varnames_map_file, init_latbc_from_fg,                 &
      &                     nudge_hydro_pres, latbc_contains_qcqi,                       &
      &                     nretries, retry_wait_sec, fac_latbc_presbiascor

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------
    itype_latbc         = 0

    dtime_latbc         = -1._wp
    nlev_latbc          = -1

    latbc_filename      = "prepiconR<nroot>B<jlev>_<y><m><d><h>.nc"
    latbc_path          = "./"
    latbc_boundary_grid = ""  ! empty string means: whole domain is read for lateral boundary
    latbc_varnames_map_file = " "
    latbc_contains_qcqi = .TRUE.
    init_latbc_from_fg  = .FALSE.
    nudge_hydro_pres    = .TRUE.
    fac_latbc_presbiascor = 0._wp

    nretries            = 0
    retry_wait_sec      = 10

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('limarea_nml')
      READ(funit,NML=limarea_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('limarea_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, limarea_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, limarea_nml)
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, limarea_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! check for deprecated namelist parameters
    !----------------------------------------------------

    IF ((nlev_latbc /= -1) .AND. my_process_is_stdio()) THEN
      WRITE (0,*) " "
      WRITE (0,*) "---------------------------------------------"
      WRITE (0,*) "DEPRECATED NAMELIST PARAMETER: nlev_latbc !!!"
      WRITE (0,*) "---------------------------------------------"
      WRITE (0,*) " "
    END IF


    !----------------------------------------------------
    ! sanity checks
    !----------------------------------------------------

    IF (dtime_latbc > 86400._wp) THEN
      CALL finish(routine, "Namelist setting of limarea_nml/dtime_latbc too large for mtime conversion!")
    END IF

    IF ( (itype_latbc > 0) .AND. (dtime_latbc <= 0._wp) ) THEN
      CALL finish(routine, "Illegal setting of limarea_nml/dtime_latbc. Must be > 0!")
    END IF


    !----------------------------------------------------
    ! Fill the configuration state
    !----------------------------------------------------

    latbc_config%itype_latbc         = itype_latbc
    latbc_config%dtime_latbc         = dtime_latbc
    latbc_config%latbc_filename      = latbc_filename
    latbc_config%latbc_path          = TRIM(latbc_path)//'/'
    latbc_config%latbc_boundary_grid = latbc_boundary_grid
    latbc_config%latbc_contains_qcqi = latbc_contains_qcqi
    latbc_config%lsparse_latbc       = (LEN_TRIM(latbc_boundary_grid) > 0)
    latbc_config%latbc_varnames_map_file = latbc_varnames_map_file
    latbc_config%init_latbc_from_fg  = init_latbc_from_fg
    latbc_config%nudge_hydro_pres    = nudge_hydro_pres
    latbc_config%fac_latbc_presbiascor = fac_latbc_presbiascor
    latbc_config%nretries            = nretries
    latbc_config%retry_wait_sec      = retry_wait_sec

    ! convert dtime_latbc into mtime object
    dtime_latbc_in_ms = 1000._wp * dtime_latbc
    CALL getPTStringFromMS(NINT(dtime_latbc_in_ms,i8), dtime_latbc_str)
    latbc_config%dtime_latbc_mtime => newTimedelta(dtime_latbc_str, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in initialization of dtime_latbc time delta.")


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      funit = open_tmpfile()
      WRITE(funit,NML=limarea_nml)
      CALL store_and_close_namelist(funit, 'limarea_nml')
    ENDIF

    !-----------------------------------------------------
    !write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=limarea_nml)

  END SUBROUTINE read_limarea_namelist

END MODULE mo_limarea_nml
