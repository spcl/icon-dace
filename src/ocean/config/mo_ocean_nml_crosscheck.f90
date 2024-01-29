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

MODULE mo_ocean_nml_crosscheck

  USE mo_master_control,    ONLY: get_my_process_name
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish, warning
  USE mo_grid_config,       ONLY: init_grid_configuration
  USE mo_parallel_config,   ONLY: check_parallel_configuration, p_test_run, l_fast_sum, &
      &                           use_dp_mpi2io, proc0_shift
  USE mo_run_config,        ONLY: nsteps, dtime, nlev
  USE mo_time_config,       ONLY: time_config, dt_restart
  USE mo_io_config,         ONLY: dt_checkpoint, write_initial_state, lnetcdf_flt64_output
  USE mo_grid_config,       ONLY: grid_rescale_factor, use_duplicated_connectivity
  USE mo_ocean_nml
  USE mo_master_config,     ONLY: isRestart
  USE mo_time_management,   ONLY: compute_timestep_settings,                        &
    &                             compute_restart_settings,                         &
    &                             compute_date_settings

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ocean_crosscheck

CONTAINS
  

!<Optimize:inUse>
  SUBROUTINE check_thicknesses

    ! ensure, that all used thicknesses are non-zero in case of non-shallow-water.
    !For shallow-water this test makes no sense since dzlev==0 
    IF(iswm_oce==0)THEN
      IF (MINVAL(dzlev_m(1:n_zlev)) <= 0.0_wp) THEN
        CALL finish("check_thicknesses","Found zero or negative thicknesses")
      END IF
    ENDIF  
  END SUBROUTINE check_thicknesses

!<Optimize:inUse>
  SUBROUTINE ocean_crosscheck()

    CHARACTER(len=*), PARAMETER :: method_name =  'mo_ocean_nml_crosscheck:ocean_crosscheck'

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

    ! set the patch-related nlev variable to the ocean setup n_ zlev
    nlev = n_zlev

    IF (p_test_run .AND. l_fast_sum ) THEN                      
       CALL warning(method_name, "p_test_run sets l_fast_sum=.false.")
       l_fast_sum = .false.                                     
    ENDIF                                                       

    SELECT CASE (select_solver)
      CASE (select_gmres, select_gmres_r, select_gmres_mp_r, select_cg, select_cg_mp, select_cgj, select_bcgs, &
        & select_legacy_gmres, select_mres)
      CASE default
        CALL finish(method_name, "Unknown solver type")
    END SELECT

    IF (no_tracer < 1) THEN
      CALL warning("ocean_crosscheck", "no_tracer < 1, use_constant_mixing")
      PPscheme_type = PPscheme_Constant_type
    ENDIF

    CALL check_thicknesses

    IF  (RichardsonDiffusion_threshold < convection_InstabilityThreshold) &
      CALL finish (method_name, "RichardsonDiffusion_threshold < convection_InstabilityThreshold")

     
    IF (l_rigid_lid .AND. iswm_oce /= 1) THEN
      CALL finish(method_name, "l_rigid_lid .AND. iswm_oce /= 1")
    ENDIF

    IF (use_duplicated_connectivity) THEN
      use_duplicated_connectivity = .FALSE.
      CALL message(method_name, "Set use_duplicated_connectivity to FALSE")
    ENDIF
    
    IF (isRestart() .AND. write_initial_state) THEN
      CALL warning(method_name, "write_initial_state is disabled for restarts")
      write_initial_state = .false.
    ENDIF

    IF ((VelocityDiffusion_order == 21 .or. VelocityDiffusion_order == 213) .and. .not. laplacian_form == 1) &
      CALL finish(method_name,"harmonic+biharmonic velocity diffusion requires curl-curl form")

    IF (vert_mix_type /=1) &
       PPscheme_type = -1

    IF (proc0_shift > 0 .AND. select_solver /= select_cg) THEN
      CALL finish(method_name, "proc0_shift only works with the CG solver (select_cg=4)")
    END IF

  END SUBROUTINE ocean_crosscheck


END MODULE mo_ocean_nml_crosscheck
