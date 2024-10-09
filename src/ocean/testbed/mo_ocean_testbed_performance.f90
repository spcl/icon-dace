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

! @brief Testbed for ocean

MODULE mo_testbed_ocean_performance

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer, activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method
  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no

  ! USE mo_icon_testbed_config, ONLY: testbed_iterations

  USE mo_math_types,                ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, sea
  USE mo_math_constants,            ONLY: pi
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,                                                  &
    &                                     threshold_min_T, threshold_max_T, threshold_min_S, threshold_max_S, &
    &                                     type_3dimRelax_Temp, para_3dimRelax_Temp,                           &
    &                                     type_3dimRelax_Salt, para_3dimRelax_Salt,                           &
    &                                     iswm_oce, l_edge_based,                                             &
    &                                     FLUX_CALCULATION_HORZ, FLUX_CALCULATION_VERT,                       &
    &                                     forcing_enable_freshwater
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, test_mode
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state
  USE mo_model_domain,              ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_sea_ice_types,             ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff

  USE mo_grid_config,         ONLY: n_dom

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: ocean_test_performance

INTEGER :: testbed_iterations = 100

CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE ocean_test_performance( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, external_data, physics_parameters,             &
    & oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    CHARACTER(*), PARAMETER :: method_name = "ocean_test_performance"

    SELECT CASE (test_mode)  !  1000 - 1100
      CASE (1000)
        CALL test_prepare_tracer_transport( patch_3d, ocean_state, operators_coefficients)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT



  END SUBROUTINE ocean_test_performance

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_prepare_tracer_transport( patch_3d, ocean_state, operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    INTEGER :: timer_prep_trace_trans_0, timer_prep_trace_trans_1, timer_prep_trace_trans_2, &
      & timer_prep_trace_trans_3
    INTEGER :: iter
    CHARACTER(*), PARAMETER :: method_name = "mo_testbed_ocean_performance:ocean_test_performance"

    !---------------------------------------------------------------------
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    !---------------------------------------------------------------------
        
    !---------------------------------------------------------------------
    ! DO the tests
    timer_prep_trace_trans_0  = new_timer("prep_trace_trans_0")
    timer_prep_trace_trans_1  = new_timer("prep_trace_trans_1")
    timer_prep_trace_trans_2  = new_timer("prep_trace_trans_2")
    timer_prep_trace_trans_3  = new_timer("prep_trace_trans_3")

!     !---------------------------------------------------------------------
!     ! measure the original implementation
!     CALL work_mpi_barrier()
!     DO iter=1, testbed_iterations
!       CALL timer_start(timer_prep_trace_trans_0)
!       CALL prepare_tracer_transport_0( patch_3D, ocean_state(1), operators_coefficients)
!       CALL timer_stop(timer_prep_trace_trans_0)
!     ENDDO
!     !---------------------------------------------------------------------
!     ! measure 1st implementation
!     CALL work_mpi_barrier()
!     DO iter=1, testbed_iterations
!       CALL timer_start(timer_prep_trace_trans_1)
!       CALL prepare_tracer_transport_1( patch_3D, ocean_state(1), operators_coefficients)
!       CALL timer_stop(timer_prep_trace_trans_1)
!     ENDDO
!     !---------------------------------------------------------------------
!     ! measure 2nd implementation
!     CALL work_mpi_barrier()
!     DO iter=1, testbed_iterations
!       CALL timer_start(timer_prep_trace_trans_2)
!       CALL prepare_tracer_transport_2( patch_3D, ocean_state(1), operators_coefficients)
!       CALL timer_stop(timer_prep_trace_trans_2)
!     ENDDO
!     !---------------------------------------------------------------------
!     ! measure 3d implementation
!     CALL work_mpi_barrier()
!     DO iter=1, testbed_iterations
!       CALL timer_start(timer_prep_trace_trans_3)
!       CALL prepare_tracer_transport_3( patch_3D, ocean_state(1), operators_coefficients)
!       CALL timer_stop(timer_prep_trace_trans_3)
!     ENDDO
!     !---------------------------------------------------------------------
! 

    !---------------------------------------------------------------------
    ! print the timers
    CALL work_mpi_barrier()
    CALL message("===================", "=======================")
    CALL message(method_name, TRIM(message_text))
    CALL print_timer()
    !---------------------------------------------------------------------
    

  END SUBROUTINE test_prepare_tracer_transport
  !-------------------------------------------------------------------------


END MODULE mo_testbed_ocean_performance

