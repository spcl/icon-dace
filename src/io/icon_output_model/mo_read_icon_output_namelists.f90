! Read namelists, make sanity checks specific to each namelist and make
! a cross check once all namelists of a component are available.
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

MODULE mo_icon_output_read_namelists

  USE mo_master_control,    ONLY: get_my_process_name
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish, warning
  USE mo_parallel_config,   ONLY: check_parallel_configuration, p_test_run, l_fast_sum, &
      &                           use_dp_mpi2io
  USE mo_run_config,        ONLY: nsteps, dtime, nlev, configure_run
  USE mo_time_config,       ONLY: time_config, dt_restart
  USE mo_io_config,         ONLY: dt_checkpoint, write_initial_state, lnetcdf_flt64_output
  USE mo_grid_config,       ONLY: grid_rescale_factor, use_duplicated_connectivity, init_grid_configuration
  USE mo_master_config,     ONLY: isRestart
  USE mo_time_management,   ONLY: compute_timestep_settings,                        &
    &                             compute_restart_settings,                         &
    &                             compute_date_settings

  USE mo_mpi                 ,ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist            ,ONLY: open_nml_output, close_nml_output, position_nml, positioned, open_nml, close_nml
  USE mo_nml_annotate        ,ONLY: log_nml_settings, temp_defaults, temp_settings

  USE mo_time_nml            ,ONLY: read_time_namelist

  USE mo_parallel_nml        ,ONLY: read_parallel_namelist
  USE mo_run_nml             ,ONLY: read_run_namelist
  USE mo_io_nml              ,ONLY: read_io_namelist
  USE mo_gribout_nml         ,ONLY: read_gribout_namelist
  USE mo_dbg_nml             ,ONLY: read_dbg_namelist

  USE mo_grid_nml            ,ONLY: read_grid_namelist
  USE mo_dynamics_nml        ,ONLY: read_dynamics_namelist
  ! USE mo_extpar_nml          ,ONLY: read_extpar_namelist

  USE mo_name_list_output_init,ONLY: read_name_list_output_namelists
! #ifndef __NO_ICON_ATMO__
  USE mo_coupling_nml        ,ONLY: read_coupling_namelist
! #endif

  USE mo_icon_output_variables, ONLY: zlevels, dz_full_level
  USE mo_ocean_nml,             ONLY: n_zlev, dzlev_m

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_icon_output_namelists
  
  NAMELIST/vertical_levels_nml/ &
    & zlevels, dz_full_level
  

CONTAINS

  !---------------------------------------------------------------------
  !>
  !! Read namelists for ocean models
  !!
!<Optimize:inUse>
  SUBROUTINE read_icon_output_namelists(icon_output_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: icon_output_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
    
    INTEGER :: status
    INTEGER :: iunit
    CHARACTER(len=*), PARAMETER :: method_name = "read_icon_output_namelists"

    !-----------------------------------------------------------------
    ! Create a new file in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    !-----------------------------------------------------------------

    IF(my_process_is_stdio()) CALL open_nml_output(TRIM(icon_output_namelist_filename)//'_output')

    !-----------------------------------------------------------------
    ! Read namelists that are shared by all components of the model.
    ! This means that the same namelists with the same values are
    ! read by all components of a coupled system.
    !-----------------------------------------------------------------

    CALL read_time_namelist           (TRIM(shr_namelist_filename))

    !-----------------------------------------------------------------
    ! Read namelists that are specific to the oce model.
    ! In case of a coupled simulation, the atmosphere model may also
    ! read some of these namelists, but probably from a different
    ! ASCII file containing different values.
    !-----------------------------------------------------------------

    ! General
    !
    CALL read_parallel_namelist       (TRIM(icon_output_namelist_filename))
    CALL read_run_namelist            (TRIM(icon_output_namelist_filename))
    CALL read_io_namelist             (TRIM(icon_output_namelist_filename))
    CALL read_name_list_output_namelists (TRIM(icon_output_namelist_filename))
    CALL read_dbg_namelist            (TRIM(icon_output_namelist_filename))

    ! Grid
    !
    CALL read_grid_namelist           (TRIM(icon_output_namelist_filename))
!    CALL read_interpol_namelist       (TRIM(icon_output_namelist_filename))

    ! Dynamics, transport and physics
    ! (still using the old setup)
    !
    CALL read_dynamics_namelist       (TRIM(icon_output_namelist_filename))

    ! Boundary conditions
    !
    ! CALL read_extpar_namelist         (TRIM(icon_output_namelist_filename))

    !
    ! GRIB2 output
    CALL read_gribout_namelist        (TRIM(icon_output_namelist_filename))

    ! Coupling
    !
    CALL read_coupling_namelist       (TRIM(icon_output_namelist_filename))

    
    CALL open_nml(TRIM(icon_output_namelist_filename))
    !==================================================================
    ! NOTE: DO NOT USE STATUS FLAG in READ(nnml) WITHOUT CHECKING IT  !
    ! This will result undetected unread namelists                    !
    !==================================================================
    CALL position_nml ('vertical_levels_nml', status=status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, vertical_levels_nml) ! write defaults to temporary text file
    END IF
    SELECT CASE (status)
    CASE (positioned)
      READ (nnml, vertical_levels_nml)                         ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, vertical_levels_nml) ! write settings to temporary text file
      END IF
    END SELECT

    ! set the patch-related nlev variable to the ocean setup n_ zlev
    nlev = zlevels
    n_zlev = zlevels
    dzlev_m(1:zlevels) = dz_full_level(1:zlevels)
    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------
    IF (my_process_is_stdio()) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (my_process_is_stdio()) CALL log_nml_settings(TRIM(icon_output_namelist_filename)//".log")
    
    CALL  close_nml
    
    !-----------------------------------------------------------------
    ! Do some checks
    
    !-----------------------------------------------------------------
    CALL check_parallel_configuration()

    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings(TRIM(get_my_process_name()), dt_restart, nsteps)

    CALL init_grid_configuration

    !--------------------------------------------------------------------
    ! checking the meanings of the io settings
    !--------------------------------------------------------------------
    IF (lnetcdf_flt64_output) THEN
       CALL message(TRIM(method_name),'NetCDF output of floating point variables will be in 64-bit accuracy')
       IF (.NOT. use_dp_mpi2io) THEN
          use_dp_mpi2io = .TRUE.
          CALL message(TRIM(method_name),'--> use_dp_mpi2io is changed to .TRUE. to allow 64-bit accuracy in the NetCDF output.')
       END IF
    ELSE
       CALL message(TRIM(method_name),'NetCDF output of floating point variables will be in 32-bit accuracy')
    END IF

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run
  
    
  END SUBROUTINE read_icon_output_namelists
  !-------------------------------------------------------------------------


END MODULE mo_icon_output_read_namelists
