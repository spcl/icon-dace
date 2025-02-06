!
! Provides interface to ART-routines dealing with coagulation
!
! This module provides an interface to the ART-routines.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
!
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

MODULE mo_art_coagulation_interface

  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: min_rlcell_int
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_parallel_config,               ONLY: nproma
  USE mo_exception,                     ONLY: finish
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_run_config,                    ONLY: lart,iqr,iqnr,iqs
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config,         ONLY: kstart_tracer
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                              timer_art, timer_art_coagInt

! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_pollen, t_fields_volc, &
                                          &   t_mode_fields
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_clipping,                  ONLY: art_clip_lt
! Coagulation Routines
  USE mo_art_coagulation,               ONLY: art_calc_coag_coefficients
  
  USE omp_lib

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_coagulation_interface'

  PUBLIC  :: art_coagulation_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_coagulation_interface(pt_prog,pt_diag, dtime, p_patch, &
              &                  prm_diag, p_metrics, tracer)
!>
!! Interface for ART-routines dealing with coagulation
!!
!! This interface calls the ART-routines, if ICON has been 
!! built including the ART-package. Otherwise, this is simply a dummy 
!! routine.
!!
  TYPE(t_nh_prog), TARGET, INTENT(inout) :: & 
    &  pt_prog                           !< Prognostic variables
  TYPE(t_nh_diag), TARGET, INTENT(inout) :: &
    &  pt_diag                           !< Diagnostic variables
  REAL(wp), INTENT(IN)              :: &
    &  dtime                             !< Time step (fast physics)
  TYPE(t_patch), TARGET, INTENT(IN) :: &
    &  p_patch                           !< Patch on which computation is performed
  TYPE(t_nwp_phy_diag), INTENT(IN)  :: &
    &  prm_diag                          !< Diagnostic variables (Physics)
  TYPE(t_nh_metrics), INTENT(IN)    :: &
    &  p_metrics                         !< Metrics (dz, ...)
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:,:,:)                   !< Tracer mixing ratios [kg/kg]
  ! Local Variables
  INTEGER                 :: & 
    &  jg, jb,               & !< patch id and loop counters
    &  i_startblk, i_endblk, & !< Start and end of block loop
    &  istart, iend,         & !< Start and end of nproma loop
    &  i_rlstart, i_rlend,   & !< Relaxation start and end
    &  i_nchdom,             & !< Number of child domains
    &  kstart,               & !< top level for all aerosol processes (transport, dynamics, ...)
    &  nlev,                 & !< Number of levels (equals index of lowest full level)
    &  num_radioact,         & !< counter for number of radionuclides
    &  coag_iteration          !< Iteration counter for entire coagulation process 
                               !  (prepare, calc_coeff, finalize)
  REAL(wp),POINTER        :: &
    &  rho(:,:,:)              !< Pointer to air density [kg m-3]
  TYPE(t_mode), POINTER   :: this_mode
  CLASS(t_mode_fields), POINTER :: fields
  !-----------------------------------------------------------------------
    
  rho => pt_prog%rho
  
  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  jg         = p_patch%id
  nlev       = p_patch%nlev
  i_rlstart  = grf_bdywidth_c+1
  i_rlend    = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

  IF(lart .AND. TRIM(art_config(jg)%cart_coag_xml)/='') THEN 
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_coagInt)

    IF (art_config(jg)%lart_aerosol) THEN
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        CALL art_air_properties(pt_diag%pres(:,:,jb),pt_diag%temp(:,:,jb), &
          &                     istart,iend,1,nlev,p_art_data(jg)%air_prop%art_free_path(:,:,jb), &
          &                     p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb))
      ENDDO
      
      num_radioact = 0
      kstart = MINVAL(kstart_tracer(jg,p_art_data(jg)%aero%itr_start:p_art_data(jg)%aero%itr_end))

      !this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
      
!$omp parallel do default(shared) private(jb, istart, iend, coag_iteration, this_mode, fields)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        DO coag_iteration=1,3
          this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
          DO WHILE(ASSOCIATED(this_mode))
            ! Select type of mode
            fields=>this_mode%fields
            SELECT TYPE (fields)
              CLASS IS (t_fields_2mom)
                ! Is coagulation enabled for this mode?
                IF(fields%do_coag <= 0) THEN
                  this_mode => this_mode%next_mode
                  CYCLE
                END IF
              
                SELECT CASE (coag_iteration)
                  CASE (1)
                    ! Before coagulation, the modal parameters have to be calculated
                    CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),        &
                      &                     istart, iend, kstart, nlev, jb, tracer(:,:,jb,:))
                  CASE (2)
                    ! Calculate coagulation coefficients
                    !  - results will be stored in fields%coag_util%p_coag%coagcoeff0 and 
                    !    fields%coag_util%p_coag%coagcoeff3
                    CALL art_calc_coag_coefficients(fields, tracer(:,:,:,:),rho(:,:,jb),         &
                      &                             pt_diag%temp(:,:,jb),                         &
                      &                             p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb), &
                      &                             istart, iend,                                 &
                      &                             kstart, nlev, jb)
                  CASE (3)
                    ! Updating fields
                    CALL fields%update_number(tracer(:,:,:,:), jb, dtime,                         &
                      &                       istart, iend, kstart, nlev, rho(:,:,jb))
                    CALL fields%update_mass(tracer(:,:,:,:), jb, dtime,                           &
                      &                     istart, iend, kstart, nlev, rho(:,:,jb))
                  CASE DEFAULT
                    CALL finish(TRIM(routine)//':art_coagulation_interface', &
                      &         'ART: Too many coag_iterations')
                END SELECT ! coag_iteration 
              CLASS IS (t_fields_pollen)
                ! Do nothing
              CLASS DEFAULT
                CALL finish(TRIM(routine)//':art_coagulation_interface', &
                    &      'ART: Unknown mode field type')
            END SELECT

            this_mode => this_mode%next_mode
          END DO !associated(this_mode)
          NULLIFY(this_mode)
        END DO !coag_iteration
      END DO !< jb  
!$omp end parallel do
    
      ! ----------------------------------
      ! --- Clip the tracers
      ! ----------------------------------
!$omp parallel do default(shared) private(jb, istart, iend)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        CALL art_clip_lt(tracer(istart:iend,1:nlev,jb,:),0.0_wp)
      ENDDO
!$omp end parallel do
      
    ENDIF !lart_aerosol

    IF (timers_level > 3) CALL timer_stop(timer_art_coagInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF !lart

END SUBROUTINE art_coagulation_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_coagulation_interface

