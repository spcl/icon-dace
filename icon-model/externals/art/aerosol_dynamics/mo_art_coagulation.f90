!
! mo_art_coagulation
! This module provides subroutines for modal coagulation.
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

MODULE mo_art_coagulation
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_physical_constants,            ONLY: ak
  USE mo_exception,                     ONLY: finish
! ART
  USE mo_art_modes,                     ONLY: t_fields_2mom
  USE mo_art_modes_linked_list,         ONLY: t_mode

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_coagulation'

  PUBLIC :: art_calc_coag_coefficients

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_coag_coefficients(this_fields, tracer, rho, &
  &                                   temp, dyn_vis,            &
  &                                   istart, iend,             &
  &                                   kstart, kend,             &
  &                                   jb)
!<
! SUBROUTINE art_calc_coag_coefficients
! This subroutine calculates the coagulation coefficents which are to be used
! when >this_mode< coagulates with its partner-modes (Note: no multiplication with 0th moments yet)
! Part of Module: mo_art_coagulation
! Based on: COSMO-ART code
! Author: Sven Werchner, KIT
! Initial Release: 2017-MM-DD
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - 
!>
  TYPE(t_fields_2mom), INTENT(INOUT)::   &
    & this_fields
  REAL(wp), INTENT(INOUT)            :: &
    & tracer(:,:,:,:),                  & !< mass/number - concentration [kg kg-1 | # kg-1]
    & rho(:,:)                            !< Air density [kg m-3]
  REAL(wp),INTENT(IN)            ::     &
    & temp(:,:),                        & !< Temperature [(jc,jk)]
    & dyn_vis(:,:)                        !< Dynamic viscosity [(jc,jk)] 
                                          !  (p_art_data(jg)%air_prop%art_dyn_visc)
  INTEGER, INTENT(IN)            ::     &
    & istart, iend,                     & !< Start and end index of nproma-loop (jc)
    & kstart, kend,                     & !< Start and end index of vertical-loop (jk)
    & jb                                  !< index of block-loop
  
  REAL(wp)                       ::     &
    & k_nc,                             & !< Near-continuum prefactor
    & k_fm,                             & !< Free-molecular prefactor
    & beta0nc, beta0fm,                 & !< Coagulation coefficients 0th moment
                                          !<   near continuum and free molecular
    & beta3nc, beta3fm,                 & !< Coagulation coefficients 3rd moment
                                          !<   near continuum and free molecular
!    & kn_t, kn_p,                       & !< Knudsen-Numbers (for this_fields and partner_fields)
!    & dia_t, dia_p,                     & !< Diameter (for this_fields and partner_fields)
    & dia_r, dia_ri,                    & !< Ratio (and inverse ratio) of diameter 
                                          !    [dia_t/dia_p (dia_p/dia_t)]
    & exp_aero_t  , exp_aero_p  ,       & !< = exp(1/8*ln^2(sigma)) 
                                          !    (for this_fields and partner_fields)
    & exp_aero_t4 , exp_aero_p4 ,       & !< = exp_aero_{t|p}^4
    & exp_aero_t9 , exp_aero_p9 ,       & !< = exp_aero_{t|p}^9
    & exp_aero_t16, exp_aero_p16,       & !< = exp_aero_{t|p}^16
    & exp_aero_t25, exp_aero_t36,       & !< = exp_aero_t^{25|36}
    & exp_aero_t49, exp_aero_t64,       & !< = exp_aero_t^{49|64}
    & exp_aero_t100,                    & !< = exp_aero_t^100
    & b0,                               & !< = bm0i or bm0 depending on bi- or uni-modal coag.
    & self                                !< = 1 or 0.5 depending on bi- or uni-modal coag.
  REAL(wp), ALLOCATABLE          ::     &   
    & kn_t(:,:), kn_p(:,:),             & !< Knudsen-Numbers (for this_fields and partner_fields)
    & dia_t(:,:), dia_p(:,:),           & !< Diameter (for this_fields and partner_fields)
    & dens_t(:,:), dens_p(:,:),         & !< Density (of this_fields and partner_fields)
    & thrdmom_t(:,:), thrdmom_p(:,:)      !< Third moment (for this_fields and partner_fields)
  INTEGER                        ::     &
    & jc, jk, jn,                       & !< Loop-Indices
    & this_itr0, partner_itr0             !< xxx%itr0
  
  
  ! Constants
  REAL(wp), PARAMETER            ::     &
    & bm0  = 0.8e0_wp,                  &
    & bm0i = 0.9e0_wp,                  &
    & bm3i = 0.9e0_wp,                  &
    & a    = 1.246e0_wp
   
  ALLOCATE( kn_t(istart:iend,kstart:kend) )
  ALLOCATE( kn_p(istart:iend,kstart:kend) )
  ALLOCATE( dia_t(istart:iend,kstart:kend) )
  ALLOCATE( dia_p(istart:iend,kstart:kend) )
  ALLOCATE( dens_t(istart:iend,kstart:kend) )
  ALLOCATE( dens_p(istart:iend,kstart:kend) )
  ALLOCATE( thrdmom_t(istart:iend,kstart:kend) )
  ALLOCATE( thrdmom_p(istart:iend,kstart:kend) )

  ! constants for this_fields
  exp_aero_t   = this_fields%info%exp_aero
  exp_aero_t4  = exp_aero_t**4
  exp_aero_t9  = exp_aero_t**9
  exp_aero_t16 = exp_aero_t**16
  exp_aero_t25 = exp_aero_t**25
  exp_aero_t36 = exp_aero_t**36
  exp_aero_t49 = exp_aero_t**49
  exp_aero_t64 = exp_aero_t**64
  exp_aero_t100= exp_aero_t**100
  
  this_itr0 = this_fields%itr0

  ! Knudsen numbers
  kn_t(:,:) = this_fields%knudsen_nr(istart:iend,kstart:kend,jb)
  ! Diameters
  dia_t(:,:)  = this_fields%diameter(istart:iend,kstart:kend,jb)
  ! Density
  dens_t(:,:)  = this_fields%density(istart:iend,kstart:kend,jb)
  ! Third moment
  thrdmom_t(:,:) = this_fields%third_moment(istart:iend,kstart:kend,jb)

  DO jn = 1,this_fields%coag_util%n_modes
    SELECT TYPE(partner_fields => this_fields%coag_util%p_coag(jn)%p_coagulateWith)
      CLASS IS(t_fields_2mom)
        !partner_fields = this_fields%coag_util%p_coag(jn)%p_coagulateWith
        partner_itr0 = partner_fields%itr0
        ! constants for partner_fields
        exp_aero_p   = partner_fields%info%exp_aero
        exp_aero_p4  = exp_aero_p**4
        exp_aero_p9  = exp_aero_p**9
        exp_aero_p16 = exp_aero_p**16

        ! Knudsen numbers
        kn_p(:,:) = partner_fields%knudsen_nr(istart:iend,kstart:kend,jb)
        ! Diameters
        dia_p(:,:) = partner_fields%diameter(istart:iend,kstart:kend,jb)
        ! Density
        dens_p(:,:) = partner_fields%density(istart:iend,kstart:kend,jb)
        ! Third moment
        thrdmom_p(:,:) = partner_fields%third_moment(istart:iend,kstart:kend,jb)

        b0 = bm0i
        self = 1.0_wp
        
        !Self-Coagulation: change numbers correction factor AND no change in mass AND take half
        IF (TRIM(ADJUSTL(this_fields%name)) == TRIM(ADJUSTL(partner_fields%name))) THEN
          b0 = bm0
          self = 0.5_wp
          this_fields%coag_util%p_coag(jn)%coagcoeff3(:,:,jb) = 0.0_wp
        END IF

        DO jk = kstart,kend
!NEC$ ivdep
          DO jc = istart,iend
            IF ( thrdmom_t(jc,jk) <= 0.0_wp                .OR.                    &
               & thrdmom_p(jc,jk) <= 0.0_wp                .OR.                    &
               & (tracer(jc,jk,jb,this_itr0   ) <  1.0_wp  .AND.                   &
               & thrdmom_t(jc,jk) <  10.e-25_wp)           .OR.                    &
               & (tracer(jc,jk,jb,partner_itr0) <  1.0_wp  .AND.                   &
               & thrdmom_p(jc,jk) <  10.e-25_wp)         ) THEN
              this_fields%coag_util%p_coag(jn)%coagcoeff0(jc,jk,jb) = 0.0_wp
              this_fields%coag_util%p_coag(jn)%coagcoeff3(jc,jk,jb) = 0.0_wp
              CYCLE
            END IF
            
            k_nc = 2.0_wp/3.0_wp * ak * temp(jc,jk) / dyn_vis(jc,jk)
            IF (dyn_vis(jc,jk) == 0.0_wp) THEN
              CALL finish('art_calc_coag_coefficients:','dyn_vis is zero!')
            END IF
            k_fm = SQRT(6.0_wp * ak * temp(jc,jk) / (dens_t(jc,jk) + dens_p(jc,jk)))
            IF (dens_t(jc,jk) + dens_p(jc,jk) == 0.0_wp) THEN
              CALL finish('art_calc_coag_coefficients:','sum of density is zero!')
            END IF
            dia_r  = dia_t(jc,jk) / dia_p(jc,jk)
            dia_ri = dia_p(jc,jk) / dia_t(jc,jk)
            IF (dia_p(jc,jk) == 0) THEN
              CALL finish('art_calc_coag_coefficients:','dia_p is zero!')
            END IF
            IF (dia_t(jc,jk) == 0) THEN
              CALL finish('art_calc_coag_coefficients:','dia_t is zero!')
            END IF
            ! 0th moment
            beta0nc = self * k_nc * ( 2.0_wp + a * kn_t(jc,jk)                                    &
              &     * ( exp_aero_t4 + dia_ri * exp_aero_t16 * exp_aero_p4 )                       &
              &     + a * kn_p(jc,jk) * ( exp_aero_p4 + dia_r  * exp_aero_p16 * exp_aero_t4 )     &
              &     + (dia_r + dia_ri) * exp_aero_p4 * exp_aero_t4 )
            beta0fm = self * k_fm * b0 * SQRT(dia_t(jc,jk)) * (exp_aero_t + SQRT(dia_ri) * exp_aero_p&
              &                               + 2.0_wp * dia_ri * exp_aero_t * exp_aero_p4        &
              &                               + dia_ri**2.0_wp * exp_aero_t9 * exp_aero_p16       &
              &                               + dia_r**(3.0/2.0) * exp_aero_t16 * exp_aero_p9     &
              &                               + 2.0_wp * SQRT(dia_r) * exp_aero_t4 * exp_aero_p )
            this_fields%coag_util%p_coag(jn)%coagcoeff0(jc,jk,jb) = beta0nc * beta0fm             &
              &                                                   / ( beta0nc + beta0fm )         &
              &                                         * tracer(jc,jk,jb,this_itr0) * rho(jc,jk) &
              &                                         * tracer(jc,jk,jb,partner_itr0) * rho(jc,jk)
            IF(beta0nc + beta0fm == 0.0_wp) THEN
              CALL finish('art_calc_coag_coefficients:','beta0nc + beta0fm is zero!')
            END IF
            ! 3rd moment if not self coagulation
            IF (self == 1.0_wp) THEN
              beta3nc = self * k_nc * dia_t(jc,jk)**3.0_wp * (2.0_wp*exp_aero_t36 + a * kn_t(jc,jk)&
                &     * ( exp_aero_t16 + dia_ri * exp_aero_t4 * exp_aero_p4 ) + a * kn_p(jc,jk)   &
                &     * ( exp_aero_t36 * exp_aero_p4 + dia_r * exp_aero_t64 * exp_aero_p16 )      &
                &     + dia_ri * exp_aero_t16 * exp_aero_p4 + dia_r * exp_aero_t64 * exp_aero_p4 )
              beta3fm = self * k_fm * bm3i * dia_t(jc,jk)**(7.0/2.0)                              &
                &     * ( exp_aero_t49 + SQRT(dia_ri) * exp_aero_t36 * exp_aero_p                 &
                &     + 2.0_wp * dia_ri * exp_aero_t25 * exp_aero_p4                              &
                &     + dia_ri**2.0_wp * exp_aero_t9 * exp_aero_p16                               &
                &     + dia_r**(3.0/2.0) * exp_aero_t100 * exp_aero_p9                            &
                &     + 2.0_wp * SQRT(dia_r) * exp_aero_t64 * exp_aero_p )
              this_fields%coag_util%p_coag(jn)%coagcoeff3(jc,jk,jb) = beta3nc * beta3fm           &
                &                                                   / ( beta3nc + beta3fm )       &
                &                                        * tracer(jc,jk,jb,this_itr0) * rho(jc,jk)&
                &                                        * tracer(jc,jk,jb,partner_itr0) * rho(jc,jk)
              IF(beta3nc + beta3fm == 0.0_wp) THEN
                CALL finish('art_calc_coag_coefficients:','beta3nc + beta3fm is zero!')
              END IF
            END IF
        
          END DO !jc
        END DO !jk
      CLASS DEFAULT
        WRITE(0,*) "Wrong type"
    END SELECT
  END DO !jn

  ! ----------------------------------
  ! --- Clean up
  ! ----------------------------------

  DEALLOCATE(kn_t)
  DEALLOCATE(kn_p)
  DEALLOCATE(dia_t)
  DEALLOCATE(dia_p)
  DEALLOCATE(dens_t)
  DEALLOCATE(dens_p)
  DEALLOCATE(thrdmom_t)
  DEALLOCATE(thrdmom_p)

END SUBROUTINE art_calc_coag_coefficients
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_coagulation
