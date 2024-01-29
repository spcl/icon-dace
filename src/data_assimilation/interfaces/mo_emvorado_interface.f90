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

#ifndef NOMPI
# ifdef HAVE_RADARFWO
  USE mpi,                ONLY: MPI_Barrier
# endif
# ifndef NO_MPI_CHOICE_ARG
  USE mpi,                ONLY: MPI_Isend, MPI_Recv
# endif
  USE mo_mpi,             ONLY: stop_mpi, abort_mpi, mpi_success, get_my_global_mpi_id, &
       &                        MPI_STATUS_SIZE, p_real, p_int, p_bcast, &
       &                        my_process_is_mpi_radarioroot, process_mpi_all_workroot_id, &
       &                        my_process_is_mpi_workroot, process_mpi_all_radarioroot_id, &
       &                        my_process_is_radario, &
       &                        process_mpi_radario_size, get_my_mpi_all_communicator, get_my_mpi_all_id
  USE mo_io_units,        ONLY: nerr
#endif
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
#ifdef HAVE_RADARFWO
  USE radar_data,         ONLY: icomm_radar, icomm_radario, &
       &                        num_radario, my_radario_id, time_mod_sec
  USE radar_organize,     ONLY: any_radar_activity, organize_radar
#endif
  
!==============================================================================

  IMPLICIT NONE

!==============================================================================

  PUBLIC :: radar_mpi_barrier, emvorado_radarfwo
#ifndef NOMPI
  PUBLIC :: detach_emvorado_io, exchg_with_detached_emvorado_io
#endif
  
!==============================================================================

CONTAINS

!==============================================================================

#ifndef NOMPI

  ! This is to be called on the combined radar communicator PEs ( my_process_is_radar() )
  SUBROUTINE exchg_with_detached_emvorado_io (grid_Center, grid_Subcenter)

    INTEGER, INTENT(inout) :: grid_Center(:), grid_Subcenter(:)

#ifdef HAVE_RADARFWO
    CALL p_bcast (grid_Center   , 0, icomm_radar)
    CALL p_bcast (grid_Subcenter, 0, icomm_radar)
#endif
    
  END SUBROUTINE exchg_with_detached_emvorado_io

  ! This is to be called on the radar I/O PEs only ( my_process_is_radario() )
  SUBROUTINE detach_emvorado_io (n_dom_model, radar_flag_doms_model)

    INTEGER, INTENT(in) :: n_dom_model                            ! Number of model domains
    LOGICAL, INTENT(in) :: radar_flag_doms_model (1:n_dom_model)  ! Switch for running EMVORADO for each domain

    
    INTEGER :: mpierr,                                &
               mpistatus(MPI_STATUS_SIZE),            &
               sendtag_radar = 0   ! gets SAVE attribute automatically
    REAL(kind=wp) :: sim_time

    INTEGER :: jg
    INTEGER, ALLOCATABLE   :: jg_list(:)

    ! Detach the separate radar IO procs by calling organize_radar('compute') on pure radar IO PEs:
    !  organize_radar('compute') does its job for one domain and one timestep, so we need
    !  an endless timeloop with an MPI-commun. to trigger computations and final EXIT.

    ! This involves communication among the whole icomm_radar and synchronizes
    !  itself with the workers in a separate call to config_emvorado in mo_atmo_model.f90:
    CALL config_emvorado (n_dom_model, radar_flag_doms_model(1:n_dom_model))
          
    ALLOCATE (jg_list(n_dom_model))
 
#ifdef HAVE_RADARFWO

    IF (msg_level >= 11 .AND. my_process_is_mpi_radarioroot()) THEN
      CALL message('detach_emvorado_io(): ', 'detaching radar-IO procs from workers', &
                    all_print=.TRUE.)
    END IF

    radario_async: DO

      CALL timer_start (timer_radar_asynio_barrier)

      IF (my_process_is_mpi_radarioroot()) THEN
        CALL mpi_recv (sim_time, 1, p_real, &
             process_mpi_all_workroot_id, sendtag_radar, get_my_mpi_all_communicator(), mpistatus, mpierr )
        CALL mpi_recv (jg_list, n_dom_model, p_int, &
             process_mpi_all_workroot_id, sendtag_radar+1, get_my_mpi_all_communicator(), mpistatus, mpierr )
        sendtag_radar = sendtag_radar + 2
      END IF
      CALL p_bcast (sim_time, 0, icomm_radario)
      CALL p_bcast (jg_list , 0, icomm_radario)

      CALL timer_stop  (timer_radar_asynio_barrier)
      CALL timer_start (timer_radar_asynio)

      IF (sim_time < -9999.0_wp) THEN
        EXIT radario_async
      ELSE
        DO jg = 1, n_dom_model
          IF (jg_list(jg) > -1) THEN
            IF (msg_level >= 13 .AND. my_process_is_mpi_radarioroot()) THEN
              CALL message('detach_emvorado_io(): ', 'calling organize_radar(''compute'')', &
                            all_print=.TRUE.)
            END IF
            time_mod_sec = sim_time   ! internal timer of EMVORADO
            CALL organize_radar(action='compute', itimelevel_dyn=1, itimelevel_qx=1, idom_in=jg)
            ! itimelevel_dyn=1, itimelevel_qx=1 are dummies here, but jg is needed because of the
            ! radario_dom(jg) communicators
          END IF
        END DO
        CALL timer_stop (timer_radar_asynio)
      END IF

    END DO radario_async
#endif

    DEALLOCATE (jg_list)

    CALL timer_stop (timer_radar_asynio)

#ifndef __SCT__
    IF (ltimer) THEN

      CALL radar_mpi_barrier ()

#ifdef HAVE_RADARFWO
      CALL timer_report_short_gen('EMVORADO async IO'           , &
                                   (/ timer_radar_asynio_barrier, &
                                   timer_radar_asynio        , &
                                   timer_radar_out           , &
                                   timer_radar_comm          , &
                                   timer_radar_barrier       , &
                                   timer_radar_composites    , &
                                   timer_radar_bubbles    /) , &
                                   icomm_radario, num_radario, my_radario_id, 0)
#endif
    ENDIF
#endif

#ifdef HAVE_RADARFWO
    ! Shut down async. radar I/O:
    WRITE (*,'(a,f0.1,a,i5)') 'Finished with asynchronous radar IO at time ',time_mod_sec, &
         ' sec on proc ', get_my_mpi_all_id()
    CALL stop_mpi()
    ! Info for Luis: This STOP is not an error stop but it ends the program
    ! regularly on an asynchronous PE. It is mandatory and cannot be replaced by
    ! finish()!
    STOP
#endif

  END SUBROUTINE detach_emvorado_io
#endif

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
#ifndef NOMPI
    INTEGER   :: mpierr, mpirequest
#endif

    CALL timer_start (timer_radar_tot)
      
    sim_time = getElapsedSimTimeInSeconds(mtime_current)
    csimtime(:) = ' '
    WRITE (csimtime, '(f0.1)') sim_time
#ifdef HAVE_RADARFWO
    time_mod_sec = sim_time   ! internal timer of EMVORADO
#endif

    CALL timer_stop  (timer_radar_tot)

    CALL config_emvorado (n_dom_model, radar_flag_doms_model (1:n_dom_model))
      
    CALL timer_start (timer_radar_tot)

#ifdef HAVE_RADARFWO
    DO jg = 1, n_dom_model
      IF ( radar_flag_doms_model (jg) .AND. p_patch(jg)%ldom_active .AND. any_radar_activity( idom_in=jg ) ) THEN
        jg_list(jg) = jg
      ELSE
        jg_list(jg) = -1
      END IF
    END DO

#ifndef NOMPI
    IF ( ANY(jg_list > -1) .AND. my_process_is_mpi_workroot() .AND. process_mpi_radario_size > 0) THEN
      CALL mpi_isend (sim_time, 1, p_real, &
                      process_mpi_all_radarioroot_id, sendtag_radar, get_my_mpi_all_communicator(), mpirequest, mpierr )
      CALL mpi_isend (jg_list, n_dom_model, p_int, &
                      process_mpi_all_radarioroot_id, sendtag_radar+1, get_my_mpi_all_communicator(), mpirequest, mpierr )
      sendtag_radar = sendtag_radar + 2
    END IF
#endif
      
    DO jg = 1, n_dom_model
      IF ( jg_list(jg) > -1 ) THEN
        IF (msg_level >= 11 .AND. my_process_is_mpi_workroot()) THEN
          CALL message('emvorado_radarfwo(): ', 'calling organize_radar(''compute'') for sim_time = '//TRIM(csimtime)//' s', &
                        all_print=.TRUE.)
        END IF
        CALL organize_radar('compute', ntime_dyn(jg), ntime_qx(jg), jg)    
      END IF
    END DO

#ifndef NOMPI
    ! send an END signal to the radario procs after last organize_radar:
    IF ( jstep == endstep ) THEN
      IF (my_process_is_mpi_workroot() .AND. process_mpi_radario_size > 0) THEN
        CALL mpi_isend (-9999.999_wp, 1, p_real, &
                        process_mpi_all_radarioroot_id, sendtag_radar, get_my_mpi_all_communicator(), mpirequest, mpierr )
        jg_list(:) = -999
        CALL mpi_isend (jg_list, n_dom_model, p_int, &
                        process_mpi_all_radarioroot_id, sendtag_radar+1, get_my_mpi_all_communicator(), mpirequest, mpierr )
        sendtag_radar = sendtag_radar + 2
      END IF
    END IF
#endif
#endif
      
    CALL timer_stop (timer_radar_tot)

  END SUBROUTINE emvorado_radarfwo

  SUBROUTINE radar_mpi_barrier

    INTEGER :: p_error
    
#ifndef NOMPI
#ifdef HAVE_RADARFWO
    CALL MPI_BARRIER (icomm_radar, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' global_mpi_barrier on ', get_my_global_mpi_id(), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE radar_mpi_barrier

END MODULE mo_emvorado_interface
