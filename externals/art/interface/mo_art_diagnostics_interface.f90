!
! This module provides and interface to routines calculating diagnostic
! properties of aerosols and trace gases
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_art_diagnostics_interface

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_run_config,                    ONLY: lart
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_diagInt
  USE mo_run_config,                    ONLY: iforcing
  USE mo_impl_constants,                ONLY: inwp
  USE mo_fortran_tools,                 ONLY: set_acc_host_or_device

! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_aero_optical_props,        ONLY: art_calc_aod, art_calc_bsc, art_calc_aodvar
  USE mo_art_diag_state,                ONLY: art_create_diagnostics
  USE mo_art_diagnostics,               ONLY: art_volc_diagnostics, art_radio_diagnostics, &
                                          &   art_dust_diagnostics, art_seas_diagnostics,  &
                                          &   art_soot_diagnostics, art_radio_diagnostics_dt_phy
  USE mo_art_clipping,                  ONLY: art_clip_lt
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_diagnostics_interface_init, art_diagnostics_interface
    
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_diagnostics_interface_init(jg, this_list, p_prog_list)

  INTEGER,INTENT(in)             :: &
    &   jg                            !< patch id
  TYPE(t_var_list_ptr),INTENT(INOUT) :: &
    &   this_list                     !< current list: prognostic or phys. tend.  
  TYPE(t_var_list_ptr),INTENT(IN)            :: &
    &   p_prog_list                   !< current prognostic state list (needed in diag-case when this_list=p_diag_list)
 
  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_diagInt)

    CALL art_create_diagnostics(jg, this_list, p_prog_list)

    IF (timers_level > 3) CALL timer_stop(timer_art_diagInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF ! lart

END SUBROUTINE art_diagnostics_interface_init
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_diagnostics_interface(rho, pres, p_trac, dz, hml, jg, &
  &                                  dt_phy_jg, p_sim_time, lacc)
  
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:,:),          & !< Air density
    &  pres(:,:,:),         & !< Air pressure
    &  dz(:,:,:),           & !< Layer thickness
    &  hml(:,:,:)             !< Height of main layer
  REAL(wp),INTENT(inout) :: &
    &  p_trac(:,:,:,:)        !< Tracer mixing ratios [kg/kg]
  INTEGER, INTENT(in)    :: &
    &  jg                     !< Patch id
  REAL(wp), INTENT(IN), OPTIONAL   :: &
    &  dt_phy_jg(:)           !< time interval for all physics packages on domain jg
  REAL(wp), INTENT(IN), OPTIONAL   :: &
    &  p_sim_time
  LOGICAL, OPTIONAL, INTENT(IN) :: &
    &  lacc

  REAL(wp), POINTER ::       &
        & cmd(:,:)             !< Pointer to current median diameter
  REAL(wp) ::                &
        & ini_cmd              !< Initial diameter
  INTEGER                ::  &
    &  jb,                   & !< Loop index: jb -> blocks
    &  istart, iend,         & !< Start and end index of nproma loop
    &  var_med_dia             !< control variable for varying median diameter (1=varying median 
                               !< dia, 0=constant median dia)
  LOGICAL                ::  &
    &  l_output_step,        &
    &  l_init_aod
  TYPE(t_mode), POINTER ::   &
    &  this_mode               !< Pointer to current mode
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  LOGICAL :: &                 !< OpenACC flag
    &  lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  IF ( PRESENT(dt_phy_jg) .AND. PRESENT(p_sim_time) ) THEN
    l_output_step = .FALSE.
  ELSE
    l_output_step = .TRUE.
  END IF

  IF (lart) THEN 
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_diagInt)

    art_atmo => p_art_data(jg)%atmo
    ! --- Get the loop indizes

    var_med_dia = 1

    IF (l_output_step) THEN

      !$ACC DATA PRESENT(art_atmo, p_trac) IF(lzacc)
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, istart, iend)

        ! ----------------------------------
        ! --- Clipping as convection might produce negative values
        ! ----------------------------------

        IF (iforcing == inwp) THEN
          CALL art_clip_lt(p_trac(istart:iend,1:art_atmo%nlev,jb,:),0.0_wp,lacc=lzacc)
        ENDIF

        ! --------------------------------------
        ! --- Calculate aerosol optical depths
        ! --------------------------------------
        IF (art_config(jg)%iart_dust > 0    .OR. &
          & art_config(jg)%iart_seasalt > 0 .OR. &
          & art_config(jg)%iart_fire > 0 .OR.    &
          & art_config(jg)%iart_volcano > 0) THEN
          IF (var_med_dia == 1) THEN
            l_init_aod = .TRUE.
            this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
            DO WHILE(ASSOCIATED(this_mode))
              ! Select type of mode
              SELECT TYPE (fields => this_mode%fields)
              CLASS is (t_fields_2mom)
                cmd => fields%diameter(:,:,jb)
                ini_cmd = fields%info%diameter_ini_nmb
                CALL art_calc_aodvar(rho(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb),           &
                  &                  istart, iend, art_atmo%nlev, jb, jg, p_art_data(jg), &
                  &                  cmd, ini_cmd, TRIM(fields%name),                     &
                  &                  l_init_aod)
                IF (fields%name(:4) == 'dust') l_init_aod = .FALSE.
              CLASS DEFAULT
                  ! No ARI for monodisperse particles
              END SELECT
              this_mode => this_mode%next_mode
            ENDDO !associated(this_mode)
          ELSE
            CALL art_calc_aod(rho(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb),            &
              &               istart, iend, art_atmo%nlev, jb, jg, p_art_data(jg))
          ENDIF
          CALL art_calc_bsc(rho(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb),              &
            &               istart, iend, art_atmo%nlev, jb, jg, p_art_data(jg))
        ENDIF

        ! -------------------------------------
        ! --- Calculate mineral dust products
        ! -------------------------------------
        IF (art_config(jg)%iart_dust > 0) THEN
          CALL art_dust_diagnostics( rho(:,:,jb), pres(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb),  &
            &                        istart, iend, art_atmo%nlev, jb, p_art_data(jg) )
        END IF

        ! -------------------------------------
        ! --- Calculate sea salt products
        ! -------------------------------------
        IF (art_config(jg)%iart_seasalt > 0) THEN
          CALL art_seas_diagnostics( rho(:,:,jb), p_trac(:,:,jb,:),                     &
            &                        istart, iend, art_atmo%nlev, jb, p_art_data(jg) )
        END IF

        ! -------------------------------------
        ! --- Calculate volcanic ash products
        ! -------------------------------------
        IF (art_config(jg)%iart_volcano > 0) THEN
          CALL art_volc_diagnostics( rho(:,:,jb), pres(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb), hml(:,:,jb),  &
            &                        istart, iend, art_atmo%nlev, jb, p_art_data(jg),                       &
            &                        art_config(jg)%iart_volcano )
        END IF

        ! -------------------------------------
        ! --- Calculate sea salt products
        ! -------------------------------------
        IF (art_config(jg)%iart_fire > 0) THEN
          CALL art_soot_diagnostics( rho(:,:,jb), p_trac(:,:,jb,:),                     &
            &                        istart, iend, art_atmo%nlev, jb, p_art_data(jg) )
        END IF

        ! ---------------------------------------------------------
        ! --- Calculate average and maximum radionuclide products 
        ! ---------------------------------------------------------
        IF (art_config(jg)%iart_radioact == 1) THEN
          CALL art_radio_diagnostics( jg, pres(:,:,jb), p_trac(:,:,jb,:),                &
            &                         istart, iend, art_atmo%nlev, jb, p_art_data(jg) )
        END IF

      ENDDO

      !$ACC END DATA

    ELSE

!$OMP PARALLEL
      IF ( p_sim_time > 1.e-6_wp) THEN

!$OMP DO PRIVATE(jb,istart,iend) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, istart, iend)

          ! ---------------------------------------------------------
          ! --- Calculate average and maximum radionuclide products 
          ! ---------------------------------------------------------
          IF (art_config(jg)%iart_radioact == 1) THEN
            CALL art_radio_diagnostics_dt_phy( jg, rho(:,:,jb), p_trac(:,:,jb,:),                &
              &                                istart, iend, art_atmo%nlev, jb, p_art_data(jg),  &
              &                                dt_phy_jg(:), p_sim_time )
          END IF

        ENDDO
!$OMP END DO NOWAIT

      END IF
!$OMP END PARALLEL

    END IF

    IF (timers_level > 3) CALL timer_stop(timer_art_diagInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  END IF !lart

END SUBROUTINE art_diagnostics_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_diagnostics_interface
