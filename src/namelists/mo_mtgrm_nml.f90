! Contains the setup of the variables for meteogram output.
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

MODULE mo_meteogram_nml

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_exception,          ONLY: finish
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_nml_and_att,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_meteogram_config,   ONLY: t_station_list, t_meteogram_output_config, &
    &                              meteogram_output_config, MAX_NVARS, &
    &                              MAX_NAME_LENGTH, MAX_NUM_STATIONS, FTYPE_NETCDF
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PUBLIC :: read_meteogram_namelist

CONTAINS

  !! Read Namelist for meteogram output.
  !!
  !! This subroutine
  !! - reads the Namelist for meteogram output
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_meteogram_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    ! local variables
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_mtgrm_nml::read_meteogram_namelist'

    INTEGER                        :: istat, funit, idom, istation, &
      &                               nstations, idx
    INTEGER                        :: iunit

    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------

    TYPE t_list_of_stations
      REAL(wp)              :: lat, lon
      CHARACTER (LEN=48)    :: zname
    END TYPE t_list_of_stations

    !> Flag. True, if meteogram output is enabled
    LOGICAL                         :: lmeteogram_enabled(max_dom)
    !> string with file name prefix for output file
    CHARACTER (len=MAX_NAME_LENGTH) :: zprefix(max_dom)
    !> file type (NetCDF, ...)
    INTEGER                         :: ftype(max_dom)
    !> Flag. Separate files for each PE
    LOGICAL                         :: ldistributed(max_dom)
    !> initial time step for meteogram output
    INTEGER                         :: n0_mtgrm(max_dom)
    !> output interval (in time steps)
    INTEGER                         :: ninc_mtgrm(max_dom)
    !> activate tile specific output
    LOGICAL                         :: loutput_tiles

    ! same for all patches:
    !> list of meteogram stations
    TYPE(t_list_of_stations) :: stationlist_tot(MAX_NUM_STATIONS)

    ! Positive-list of variables (optional). Only variables contained in
    ! this list are included in this meteogram. If the default list is
    ! not changed by user input, then all available variables are added
    ! to the meteogram
    CHARACTER(len=MAX_NAME_LENGTH)    :: var_list(MAX_NVARS)

    !> maximum number of time stamps stored before flush
    INTEGER :: max_time_stamps(max_dom)
    !> flush silently when time stamp buffer is exhausted or warn user?
    LOGICAL :: silent_flush(max_dom)
    !> instead of clobbering existing file, append instead
    LOGICAL :: append_if_exists(max_dom)

    !> Namelist for meteogram output
    NAMELIST/meteogram_output_nml/ lmeteogram_enabled, zprefix, ldistributed,            &
         &                         loutput_tiles, n0_mtgrm, ninc_mtgrm, stationlist_tot, &
         &                         var_list, max_time_stamps, silent_flush,   &
         &                         append_if_exists

    CHARACTER(len=*), PARAMETER :: &
         nml_name = 'meteogram_output_nml'

    !-----------------------
    ! 1. default settings
    !-----------------------

    lmeteogram_enabled(:)    =          .FALSE.
    zprefix(:)               =     "METEOGRAM_"
    ftype(:)                 =     FTYPE_NETCDF
    ldistributed(:)          =           .TRUE.
    loutput_tiles            =          .FALSE.
    n0_mtgrm(:)              =               0
    ninc_mtgrm(:)            =               1
    var_list(:)              =             " "
    stationlist_tot(:)%lon   = 0._wp
    stationlist_tot(:)%lat   = 0._wp
    stationlist_tot(:)%zname = ""
    ! To determine the length of the provided list, we set the default
    ! value of the longitude to "-1."
    stationlist_tot(:)%lon   = -1._wp
    stationlist_tot(1)%zname = "Hamburg"
    stationlist_tot(1)%lon   = 53.633_wp
    stationlist_tot(1)%lat   =  9.983_wp

    ! "max_time_stamps" : Maximum no. of time stamps which can be
    !                     sampled in the meteogram buffer.
    !
    ! Note that this parameter directly the affects the memory
    ! consumption on the PE which samples the meteograms! The
    ! following formula calculates the buffer size and gives typical
    ! values:
    !
    !    max_nlevs        = 91
    !    nvars            = 27
    !    nsfcvars         = 61
    !    
    !    p_real_dp_byte   = 8
    !    MAX_HEADER_SIZE  = 60
    !    MAX_DATE_LEN     = 16
    !    MAX_NUM_STATIONS = 60
    !    
    !    max_sfcvar_size  = max_time_stamps*p_real_dp_byte
    !    max_var_size     = (max_time_stamps+1)*p_real_dp_byte*max_nlevs
    !    max_buf_size     = MAX_HEADER_SIZE*p_real_dp_byte + max_time_stamps*(MAX_DATE_LEN+4)
    !                       + nvars*max_var_size  + nsfcvars*max_sfcvar_size              
    !    total            = max_buf_size * MAX_NUM_STATIONS
    !
    ! This means, for example, that sampling 10000 time stamps at once
    ! in the buffer would allocate 12 GiB of memory! Be careful!
    !
    max_time_stamps          =    200

    silent_flush             = .FALSE.
    append_if_exists         = .FALSE.
    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist(nml_name)
      READ(funit,NML=meteogram_output_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(filename)
    CALL position_nml (nml_name, status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, meteogram_output_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, meteogram_output_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, meteogram_output_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    ! determine length of station list:
    DO nstations = 0, MAX_NUM_STATIONS - 1
      IF (stationlist_tot(nstations + 1)%lon == -1._wp) EXIT
    END DO

    ! fill in values for each model domain:
    DO idom=1,max_dom
      meteogram_output_config(idom)%lenabled     = lmeteogram_enabled(idom)
      meteogram_output_config(idom)%zprefix      = zprefix(idom)
      meteogram_output_config(idom)%ftype        = ftype(idom)
#ifdef NOMPI
      meteogram_output_config(idom)%ldistributed = .TRUE.
#else
      meteogram_output_config(idom)%ldistributed = ldistributed(idom)
#endif
      meteogram_output_config(idom)%loutput_tiles= loutput_tiles
      meteogram_output_config(idom)%n0_mtgrm     = n0_mtgrm(idom)
      meteogram_output_config(idom)%ninc_mtgrm   = ninc_mtgrm(idom)
      meteogram_output_config(idom)%nstations    = nstations
      meteogram_output_config(idom)%var_list     = var_list

      meteogram_output_config(idom)%max_time_stamps = max_time_stamps(idom)
      meteogram_output_config(idom)%silent_flush = silent_flush(idom)
      meteogram_output_config(idom)%append_if_exists = append_if_exists(idom)

      ALLOCATE(meteogram_output_config(idom)%station_list(nstations))

      meteogram_output_config(idom)%station_list(:)%zname        = &
        &  stationlist_tot(1:nstations)%zname
      meteogram_output_config(idom)%station_list(:)%location%lon = &
        &  stationlist_tot(1:nstations)%lon
      meteogram_output_config(idom)%station_list(:)%location%lat = &
        &  stationlist_tot(1:nstations)%lat

    END DO

    ! consistency check
    idx = 0
    DO idom=1,max_dom
      IF (meteogram_output_config(idom)%lenabled) THEN
        idx = idom
        EXIT
      END IF
    END DO
    DO idom=(idx+1),max_dom
      IF (.NOT. meteogram_output_config(idom)%lenabled) CYCLE
      IF (meteogram_output_config(idom)%ldistributed .NEQV. meteogram_output_config(idx)%ldistributed) THEN
        CALL finish( TRIM(routine), "Inconsistent namelist setting for domains (ldistributed)!")
      END IF
    END DO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=meteogram_output_nml)
      CALL store_and_close_namelist(funit, nml_name)
    ENDIF

    !-----------------------------------------------------
    ! 6. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=meteogram_output_nml)
    END IF

  END SUBROUTINE read_meteogram_namelist

END MODULE mo_meteogram_nml
