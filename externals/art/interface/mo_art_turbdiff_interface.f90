!
! Provides interface to ART-routines dealing with turbulent diffusion
! using the COSMO turbulence scheme of Matthias Raschendorfer.
!
! This module provides an interface to the ART-routines.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
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

MODULE mo_art_turbdiff_interface


  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_nonhydro_types,                ONLY: t_nh_metrics, t_nh_diag, t_nh_prog
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE turb_data,                        ONLY: modvar
  USE mo_run_config,                    ONLY: lart
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_turbdiffInt
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_diag_types,                ONLY: t_art_diag
  USE mo_art_surface_value,             ONLY: art_surface_value
  USE mo_art_config,                    ONLY: art_config
  USE mo_parallel_config,               ONLY: nproma
  USE mo_fortran_tools,                 ONLY: assert_acc_device_only
#ifdef _OPENACC
  USE mo_exception,                     ONLY: warning
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_turbdiff_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_turbdiff_interface( defcase,  & !>in
    &          p_patch,                      & !>in
    &          p_prog_rcf,                   & !>in
    &          prm_nwp_tend,                 & !>in
    &          ncloud_offset,                & !>in
    &          ptr,                          & !>out
    &          idx_nturb_tracer,             & !>out, optional
    &          p_rho,                        & !>in, optional
    &          p_metrics, p_diag, prm_diag,  & !>in, optional
    &          jb,                           & !>in, optional
    &          opt_sv, opt_fc,               & !>in, optional
    &          i_st, i_en, dt,               & !>in, optional
    &          lacc )

!! Interface for ART-routines dealing with turbdiff
!!
!! This interface calls the ART-routines, if ICON has been
!! built including the ART-package. Otherwise, this is simply a dummy
!! routine.
!!
  CHARACTER(len=*), INTENT(in)              :: &
    &  defcase                           !< definition of case
  TYPE(t_patch), TARGET, INTENT(in)         :: &
    &  p_patch                           !< patch on which computation is performed
  TYPE(t_nh_prog),  INTENT(in)              :: &
    &  p_prog_rcf                        !< the prog vars
  TYPE(t_nwp_phy_tend), INTENT(in)          :: &
    &  prm_nwp_tend                      !< atm phys tendencies
  INTEGER, INTENT(IN)                       :: &
    &  ncloud_offset                     !< index offset due to additional cloud variables 
                                         !< to be diffused.
  TYPE(modvar), INTENT(inout)               :: &
    &  ptr(:)                            !< passive tracer pointer type structure for diffusion
  REAL(wp), INTENT(in)                      :: &
    &  dt

  INTEGER, INTENT(inout), OPTIONAL          :: &
    &  idx_nturb_tracer(:)               !< indices of the turbulent tracers in the prognostic list
  REAL(wp), INTENT(in), OPTIONAL            :: &
    &  p_rho(:,:,:)                      !< air density
  TYPE(t_nh_metrics), INTENT(in), OPTIONAL  :: &
    &  p_metrics                         !< metrical fields
  TYPE(t_nh_diag), INTENT(in), OPTIONAL     :: &
    &  p_diag                            !< the diag vars
  TYPE(t_nwp_phy_diag), INTENT(in), OPTIONAL:: &
    &  prm_diag                          !< atm phys vars
  INTEGER, INTENT(in), OPTIONAL             :: &
    &  jb, i_st, i_en
  REAL(wp), INTENT(out), TARGET, OPTIONAL   :: &
    &  opt_sv(:,:,:)                     !< surface value according to transfer coeff.
  LOGICAL, INTENT(in), OPTIONAL             :: &
    &  opt_fc
  LOGICAL, INTENT(in), OPTIONAL             :: &
    &  lacc


!Local variables
  REAL(wp), POINTER         :: &
    &  sv(:,:,:),              & !< surface value of tracer
    &  vdep(:,:,:)               !< deposition velocity of tracer
  TYPE(t_art_diag), POINTER :: &
    &  art_diag                  !< Pointer to ART diagnostic fields
  INTEGER                   :: &
    &  jg, idx_trac, jk, jc,   & !< loop indices
    &  nlev,                   & !<
    &  idx_tot                   !< counter for total number of fields 
                                 !< (add. cloud vars + tracer vars) to be diffused
  REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS :: &
    & dvar_av, dvar_at
  
  CALL assert_acc_device_only("art_turbdiff_interface", lacc)
  
  jg  = p_patch%id
  IF ( lart ) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_turbdiffInt)

    SELECT CASE(TRIM(defcase))
    
    CASE('setup_ptr')

      nlev = p_patch%nlev

      sv => p_art_data(jg)%turb_fields%sv
      vdep => p_art_data(jg)%turb_fields%vdep

      art_diag => p_art_data(jg)%diag

      IF ( .NOT. PRESENT(opt_sv) ) THEN
        CALL art_surface_value( p_patch, p_prog_rcf, p_metrics, p_diag, prm_diag, &
          &                     art_diag, p_rho(:,:,:), dt, jg, jb, vdep, sv, lacc=lacc )
      END IF

      DO idx_trac = 1, art_config(jg)%nturb_tracer

        idx_tot = idx_trac + ncloud_offset

        ! set up pointer to tracer type structure for diffusion
        ptr(idx_tot)%av => p_prog_rcf%turb_tracer(jb,idx_trac)%ptr
        ptr(idx_tot)%at => prm_nwp_tend%turb_tracer_tend(jb,idx_trac)%ptr
        dvar_at => ptr(idx_tot)%at

        !$ACC DATA PRESENT(dvar_at)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = 1, nproma
            dvar_at(jc,jk) = 0._wp
          END DO
        END DO
        !$ACC END PARALLEL

        !$ACC END DATA

        IF ( PRESENT(opt_sv) ) THEN
          ptr(idx_tot)%sv => opt_sv(:,jb,idx_trac)
          IF ( PRESENT(opt_fc) ) THEN
            ptr(idx_tot)%fc = opt_fc
          ELSE
            ptr(idx_tot)%fc = .FALSE.
          END IF
        ELSE
          ptr(idx_tot)%sv => sv(:,jb,idx_trac)
          ptr(idx_tot)%fc = .FALSE.
        END IF

        IF (art_config(jg)%lart_emiss_turbdiff) THEN
          IF (p_art_data(jg)%emiss%exists(p_prog_rcf%turb_tracer(jb,idx_trac)%idx_tracer)) THEN

            ptr(idx_tot)%fc = .TRUE.
          END IF
        END IF

        ! save the index of the current turbulent tracer in the prognostic list
        IF (PRESENT(idx_nturb_tracer)) THEN
           idx_nturb_tracer(idx_trac) = p_prog_rcf%turb_tracer(jb,idx_trac)%idx_tracer
        END IF

      END DO

    CASE('update_ptr')

      nlev = p_patch%nlev !< Number of vertical full levels
      ! update tracers due to diffusion

      !$ACC DATA COPYIN(i_en, i_st)

      DO idx_trac = 1, art_config(jg)%nturb_tracer

        idx_tot = idx_trac + ncloud_offset

        dvar_av => ptr(idx_tot)%av
        dvar_at => ptr(idx_tot)%at

        !$ACC DATA PRESENT(dvar_at, dvar_av)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
!NEC$ IVDEP
          DO jc = i_st, i_en
            dvar_av(jc,jk) = MAX( 0._wp, dvar_av(jc,jk)     &
              &             + dt * dvar_at(jc,jk) )
          END DO
        END DO
        !$ACC END PARALLEL

        !$ACC END DATA

      END DO

      !$ACC WAIT
      !$ACC END DATA

    END SELECT

    IF (timers_level > 3) CALL timer_stop(timer_art_turbdiffInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  END IF

END SUBROUTINE art_turbdiff_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_turbdiff_interface

