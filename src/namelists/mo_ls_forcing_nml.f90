! Namelist file for large-scale forcing terms
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

MODULE mo_ls_forcing_nml

  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_kind,                ONLY: wp
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_run_config,          ONLY: ltestcase
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_grid_config,         ONLY: is_plane_torus
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_ls_forcing_namelist, is_ls_forcing, is_subsidence_moment, is_subsidence_heat, &
            is_advection, is_advection_uv,is_advection_tq,is_geowind, is_rad_forcing, is_theta,&
            is_nudging, is_nudging_uv, is_nudging_tq, nudge_start_height, nudge_full_height,   &
            dt_relax, is_sim_rad

  LOGICAL  :: is_ls_forcing         !true if any forcing is on
  LOGICAL  :: is_subsidence_moment  !true if subsidence is on for u and v
  LOGICAL  :: is_subsidence_heat    !true if subsidence is on for thermodyn. variables
  LOGICAL  :: is_advection          !true if horizontal advective forcing is on for any variable
  LOGICAL  :: is_advection_uv       !true if horizontal advective forcing is on for u and v
  LOGICAL  :: is_advection_tq       !true if horizontal advective forcing is on for temperature and moisture
  LOGICAL  :: is_geowind            !true if geostophic wind is set 
  LOGICAL  :: is_rad_forcing        !true if radiative forcing is on
  LOGICAL  :: is_theta              !true is forcings are in terms of theta
  LOGICAL  :: is_nudging            !true if nudging applied
  LOGICAL  :: is_nudging_uv         !true if nudging applied to u and v
  LOGICAL  :: is_nudging_tq         !true if nudging applied to temperature and moisture
  LOGICAL  :: is_sim_rad            !true if simplified radiation scheme should be used
  REAL(wp) :: nudge_start_height    !height where nudging starts                [m]
  REAL(wp) :: nudge_full_height     !height where nudging reaches full strength [m]
  REAL(wp) :: dt_relax              !time scale for nudging                     [s]
 
  NAMELIST/ls_forcing_nml/ is_subsidence_moment, is_subsidence_heat, is_advection,is_advection_uv,is_advection_tq, &
                           is_geowind, is_rad_forcing, is_theta, is_nudging,is_nudging_uv, is_nudging_tq, &
                           nudge_start_height, nudge_full_height, dt_relax, is_sim_rad

CONTAINS
  !-------------------------------------------------------------------------
  !! Read Namelist for LS forcing
  !!
  !! This subroutine 
  !! - reads the Namelist 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state 
  !!
  SUBROUTINE read_ls_forcing_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_ls_forcing_nml: read_ls_forcing_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    is_ls_forcing        = .FALSE.
    is_subsidence_moment = .FALSE.
    is_subsidence_heat   = .FALSE.
    is_advection         = .FALSE.
    is_advection_uv      = .TRUE.
    is_advection_tq      = .TRUE.
    is_geowind           = .FALSE.
    is_rad_forcing       = .FALSE.
    is_theta             = .FALSE.
    is_nudging           = .FALSE.
    is_nudging_uv        = .TRUE.
    is_nudging_tq        = .TRUE.
    is_sim_rad           = .FALSE.
    nudge_start_height   = 1000.0_wp
    nudge_full_height    = 2000.0_wp
    dt_relax             = 3600.0_wp

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    ! ltestcase is added here because it was causing trouble for AMIP runs
    ! restarting with a file generated long ago
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN 
      funit = open_and_restore_namelist('ls_forcing_nml')
      READ(funit,NML=ls_forcing_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ls_forcing_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ls_forcing_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ls_forcing_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ls_forcing_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !4. checks
    !If any atmospheric forcing is ON, turn on is_ls_forcing
    ! see configure_model/mo_nml_crosscheck.f90 for surface forcing in SCM setup 
    IF(is_subsidence_moment .OR. is_subsidence_heat .OR. is_advection .OR. &
      & is_geowind .OR. is_rad_forcing .OR. is_nudging .OR. is_sim_rad) THEN
        is_ls_forcing = .TRUE.
    END IF






    IF(is_ls_forcing .AND. .NOT.ltestcase) &
        CALL message(routine,'ltestcase is turned ON because is_ls_forcing is ON!')

    !Check for testcases with large-scale forcing
    ! fixme: move to nml crosscheck
    IF (is_rad_forcing .AND. atm_phy_nwp_config(1)%inwp_radiation>0) &
      CALL finish(routine, 'both inwp_rad and rad_forcing are turned on!')

    IF (is_geowind .AND. .NOT.is_plane_torus) &
      CALL finish(routine, 'is_geowind is only applicable if is_plane_torus is turned on!')

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=ls_forcing_nml)                    
      CALL store_and_close_namelist(funit,'ls_forcing_nml') 
    ENDIF

    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=ls_forcing_nml)

  END SUBROUTINE read_ls_forcing_namelist

END MODULE mo_ls_forcing_nml
