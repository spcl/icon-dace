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

MODULE mo_emvorado_interface

  USE mo_model_domain,    ONLY: p_patch
  USE mo_timer,           ONLY: timer_start, timer_stop   , &
       &                        timer_radar_tot           , &
       &                        timer_radar_asynio_barrier, &
       &                        timer_radar_asynio        , &
       &                        timer_radar_out           , &
       &                        timer_radar_comm          , &
       &                        timer_radar_barrier       , &
       &                        timer_radar_composites    , &
                                timer_radar_bubbles                          
  USE mo_real_timer,      ONLY: timer_report_short_gen
  USE mo_kind,            ONLY: sp, dp, wp
  USE mo_run_config,      ONLY: ltimer, nsteps, msg_level
  USE mo_exception,       ONLY: message
  USE mtime,              ONLY: datetime
  USE mo_util_mtime,      ONLY: getElapsedSimTimeInSeconds

  USE mo_emvorado_config, ONLY: config_emvorado
  
!==============================================================================

  IMPLICIT NONE

!==============================================================================

  PUBLIC :: radar_mpi_barrier, emvorado_radarfwo
  
!==============================================================================

CONTAINS

!==============================================================================


  ! This is to be called on the workers ( my_process_is_work() )
  SUBROUTINE emvorado_radarfwo (mtime_current, ntime_dyn, ntime_qx, n_dom_model, radar_flag_doms_model, jstep, endstep)

    TYPE(datetime), INTENT(in) :: mtime_current
    INTEGER, INTENT(in)        :: n_dom_model                            ! Number of model domains
    INTEGER, INTENT(in)        :: ntime_dyn (1:n_dom_model)              ! time level for u, v, w, t, p, rho
    INTEGER, INTENT(in)        :: ntime_qx  (1:n_dom_model)              ! time level for qv, qc, qr, qi, ...
    LOGICAL, INTENT(in)        :: radar_flag_doms_model (1:n_dom_model)  ! Switch for running EMVORADO for each domain
    INTEGER, INTENT(in)        :: jstep, endstep
    
    REAL(wp)  :: sim_time     !< elapsed simulation time
    INTEGER   :: sendtag_radar = 0
    INTEGER   :: jg, jg_list(n_dom_model)
    CHARACTER(len=12) :: csimtime

    CALL timer_start (timer_radar_tot)
      
    sim_time = getElapsedSimTimeInSeconds(mtime_current)
    csimtime(:) = ' '
    WRITE (csimtime, '(f0.1)') sim_time

    CALL timer_stop  (timer_radar_tot)

    CALL config_emvorado (n_dom_model, radar_flag_doms_model (1:n_dom_model))
      
    CALL timer_start (timer_radar_tot)

      
    CALL timer_stop (timer_radar_tot)

  END SUBROUTINE emvorado_radarfwo

  SUBROUTINE radar_mpi_barrier

    INTEGER :: p_error
    

  END SUBROUTINE radar_mpi_barrier

END MODULE mo_emvorado_interface
