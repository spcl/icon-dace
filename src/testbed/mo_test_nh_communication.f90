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

MODULE mo_test_nh_communication

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, proc_split, push_glob_comm, pop_glob_comm
  USE mo_timer,               ONLY: ltimer, print_timer

  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_icon_testbed_config, ONLY: testbed_iterations

  USE mo_model_domain,        ONLY: p_patch  
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  
  USE mo_icon_comm_lib
  USE mo_atmo_nonhydrostatic, ONLY: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic
  USE mo_async_latbc_types,   ONLY: t_latbc_data

  !-------------------------------------------------------------------------
  ! for the nh run
  USE mo_nonhydro_types,      ONLY: t_nh_state
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: n_dom, n_dom_start
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  
  !-------------------------------------------------------------------------

!-------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE

PUBLIC :: test_nh_communication

CONTAINS
  

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_nh_communication(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    ! 2D variables    

    
    INTEGER :: patch_no, i

    TYPE(t_latbc_data) :: latbc !< data structure for async latbc prefetching

    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_nh_communication"

    !---------------------------------------------------------------------
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    !---------------------------------------------------------------------

    ltimer = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_nonhydrostatic(latbc)
        
    CALL work_mpi_barrier()
    
    !---------------------------------------------------------------------
    DO i=1,testbed_iterations
!       CALL integrate_nh_test_comm(p_nh_state, p_patch, p_int_state, datetime, p_grf_state, &
!         &               1, jstep, dtime, sim_time, 1,                            &
!         &               l_compute_diagnostic_quants                              )
    ENDDO
    !---------------------------------------------------------------------
    

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic(latbc)
    CALL destruct_atmo_model()
!     CALL destruct_icon_communication()
    CALL message(TRIM(method_name),'clean-up finished')

    !---------------------------------------------------------------------
    ! print the timers
!    IF (my_process_is_stdio()) THEN
      CALL message("===================", "=======================")
      WRITE(message_text,*) "Communication Iterations=", testbed_iterations
      CALL message(method_name, TRIM(message_text))
      CALL print_timer()
!    ENDIF
    !---------------------------------------------------------------------
     

  END SUBROUTINE test_nh_communication
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
#ifdef __NO_NESTING__
   SUBROUTINE integrate_nh_test_comm (p_nh_state, p_patch, p_int_state,  &
  &        p_grf_state, jg, nstep_global, dt_loc, dtadv_loc, sim_time,   &
  &        num_steps, l_compute_diagnostic_quants)
#else
  RECURSIVE SUBROUTINE integrate_nh_test_comm (p_nh_state, p_patch, p_int_state,  &
  &        p_grf_state, jg, nstep_global, dt_loc, dtadv_loc, sim_time, &
  &        num_steps, l_compute_diagnostic_quants)
#endif

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_stepping:integrate_nh'

    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch(n_dom_start:n_dom)    !< patch
    TYPE(t_int_state),TARGET,INTENT(in)  :: p_int_state(n_dom_start:n_dom)!< interpolation state
    TYPE(t_nh_state), TARGET, INTENT(inout) :: p_nh_state(n_dom) !< nonhydrostatic state
    TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state(n_dom_start:n_dom)!< gridref state

    INTEGER , INTENT(IN)    :: jg           !< current grid level
    INTEGER , INTENT(IN)    :: nstep_global !< counter of global time step
    INTEGER , INTENT(IN)    :: num_steps    !< number of time steps to be executed
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level
    REAL(wp), INTENT(IN)    :: dtadv_loc    !< advective time step applicable to 
                                            !< local grid level
    REAL(wp), INTENT(INOUT) :: sim_time(n_dom) !< elapsed simulation time on each
                                               !< grid level
    LOGICAL, INTENT(IN) :: l_compute_diagnostic_quants    !< computation of diagnostic quantities

   END SUBROUTINE integrate_nh_test_comm 
  !-------------------------------------------------------------------------

END MODULE mo_test_nh_communication

