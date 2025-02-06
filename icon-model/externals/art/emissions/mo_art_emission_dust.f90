!
! mo_art_emission_dust
! This module provides subroutines for the emission of mineral dust aerosol.
! Based on Vogel et al. (2006) - A model of dust transport applied to the Dead Sea area
! Extended for global application as described in Rieger (2016) - Der Einfluss von
! natuerlichem Aerosol auf Wolken ueber Mitteleuropa, Dissertation an der Fakultaet fuer
! Physik des Karlsruher Instituts fuer Technologie (KIT)
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

MODULE mo_art_emission_dust
  ! ICON Routines
  USE mo_kind,                          ONLY: wp
  USE mo_physical_constants,            ONLY: grav
  USE mo_math_constants,                ONLY: pi
  ! ART Routines
  USE mo_art_aerosol_utilities,         ONLY: calc_ustar
  USE mo_art_external_types,            ONLY: t_art_soil_table,t_art_soil_properties
    
    
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: art_prepare_emission_dust
  PUBLIC :: art_emission_dust

  INTEGER, PARAMETER ::   integ_steps = 1000  !< Number of numerical integration steps
  
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_emission_dust(jb, istart, iend, u, v, rho, tcm, dzsoil, &
  &                                  soilwater, soilice, soil_prop, hsnow,     &
  &                                  ddust1, ddust2, ddust3, lart_diag_out,    &
  &                                  ustart_min_diag, ustar_diag)
!<
! SUBROUTINE art_prepare_emission_dust
! This subroutine calculates the emission fluxes of mineral dust aerosol
! and updates the mineral dust aerosol concentration. Equation numbers are related 
! to Vogel et al., 2006.
! Based on: Vogel et al. (2006) - A model of dust transport applied to the Dead Sea Area
!           COSMO-ART code by C. HOOSE and H. VOGEL
! Part of Module: mo_art_emission_dust
! Author: Daniel Rieger, KIT
! Initial Release: 2013-05-16
! Modifications:
! 2014-11-11: Daniel Rieger, KIT
! - Unified ICON/COSMO-ART emission routine.
!>
  REAL(wp), INTENT(IN)  :: &
    &  u(:), v(:),         & !< Wind velocities in lowest model layer [m s-1]
    &  rho(:),             & !< Air density [kg m-3]
    &  tcm(:),             & !< Transfer coefficient for momentum
    &  soilwater(:),       & !< Soil water content
    &  soilice(:),         & !< Soil ice content
    &  dzsoil,             & !< Thickness of uppermost soil layer
    &  hsnow(:),           & !< snow height [m]
    &  ddust1,             & !<  1.5E-6_wp median diameter dust mode a
    &  ddust2,             & !<  6.7E-6_wp median diameter dust mode b
    &  ddust3                !< 14.2E-6_wp median diameter dust mode c
  LOGICAL, INTENT(IN)   :: &
    &  lart_diag_out         !< diagnostic output?
  INTEGER, INTENT(IN)   :: &
    &  jb,                 & !< Block index
    &  istart, iend          !< Start and end index of nproma loop
  TYPE(t_art_soil_properties),INTENT(INOUT) :: &
    &  soil_prop             !< Storage container for all fields related to soil
  REAL(wp),INTENT(INOUT)  :: & !< collect for output:
    &  ustart_min_diag(:), & !< minimum value of ustar_t (Eq. 3.4) at dp_min
    &  ustar_diag(:)         !< friction velocity
  ! Local Variables
  TYPE(t_art_soil_table),POINTER  :: &
    &  soil_tab                   !< Pointer to table for soil properties of one specific soil type
  REAL(wp)                        :: &
    &    dstart, dend,               & !< Lower and upper bound of integration (diameter)
    &    dlnd,                       & !< Logarithmical numerical integration step of ln(diameter)
    &    dint(integ_steps+1),        & !< diameter at integration
    &    log_dint(integ_steps+1),    & !< its logarithm
    &    abdint(integ_steps+1),      & !< dint with factors a and b
    &    sqrt_abdint(integ_steps+1), & !< sqrt(dint)
    &    abdint_p15(integ_steps+1),  & !< dint**1.5
    &    fac33a,fac33b,              & !< Factors needed for Eq. 3.3
    &    an33,gamma33,               & !< Factors needed for Eq. 3.3
    &    waterdens,                  & !< Density of water (=1000)
    &    cwhite,                     & !< Constant C in White et. al, 1979 Paper
    &    cstrich,                    & !< Modified cwhite: see Eq. below
    &    wbgrav,                     & !< Gravimetric water content
    &    f_eta, f_eta_wb,            & !< Factor for u* for non-ideal conditions (Eq. 3.4)
    &    q1dint, q2dint, q3dint,     & !< Integral of Eq. 3.8
    &    de3,                        & !< Current diameter at integration step integ
    &    de2,                        & !< 
    &    logsigma_min,logdgm_min,    & !< ln(sigma) and ln(dgm) for minimal dispersed
    &    logsigma_ful,logdgm_ful,    & !< ln(sigma) and ln(dgm) for fully dispersed
    &    mass_m, mass_f,             & !< mass of minimal/fully dispersed
    &    stot,                       & !< Cross sectional area of all dust particles on the ground
    &    ustar,                      & !< friction velocity
    &    dp_min,                     & !< diameter of the minimum value of ustar_t (Eq.3.4)
    &    ustart_min,                 & !< minimum value of ustar_t (Eq. 3.4) at dp_min
    &    ustart2,                    & !< ustart * ustart (square of Eq. 3.3)
    &    gamma_disperse,             & !< factor which characterizes the actual dispersion of the 
                                       !    soil =f(ustar)
    &    intfactor1,                 & !< factor of the integration of Eq. 3.7
    &    intfactor2,                 & !<
    &    ekin,                       & !< Kinetic energy
    &    ekin_e1(integ_steps+1),     & !< Difference kinetic energy and binding energy 1
    &    ekin_e2(integ_steps+1),     & !< Difference kinetic energy and binding energy 2
    &    ekin_e3(integ_steps+1),     & !< Difference kinetic energy and binding energy 3
    &    e3fac(integ_steps+1),       & !< Factor (=0 if ekin_e3<0, =1 if ekin_e3>0)
    &    ekinfac,                    & !< Prefactor to calculate kinetic energy
    &    p1,p2,p3,                   & !< 
    &    p1fac,p2fac,p3fac             !<

  REAL(wp),ALLOCATABLE              :: &
    &    emiss_a(:),                   & !< Vertical emission flux of mode a
    &    emiss_b(:),                   & !< Vertical emission flux of mode b
    &    emiss_c(:)                      !< Vertical emission flux of mode c
                                          
  REAL(wp),PARAMETER                :: & 
    &    e1       = 3.61E-7_wp,        & !< binding energy 1
    &    e2       = 3.52E-7_wp,        & !< binding energy 2
    &    e3       = 3.46E-7_wp,        & !< binding energy 3
    &    beta     = 163._wp,           & !< proportionality factor between hor. and vert. dust flux
    &    bulkdens = 1.5E3_wp,          & !< Bulk soil density
    &    rho_p    = 2.65E3_wp            !< Bulk density of particles
                                          
  INTEGER                           :: & 
    &    integ,                        & !< Counter of integration steps
    &    i,js,jn                         !< Loop indizes
    
  ! ----------------------------------
  ! --- Initializations
  ! ----------------------------------
  
  dend        = 2.E-3_wp
  an33        = 0.0123_wp
  gamma33     = 3.E-4_wp
  waterdens   = 1000._wp
  cwhite      = 0.375_wp
  
  dp_min = SQRT(gamma33/(rho_p*grav))
  
  ALLOCATE(emiss_a(soil_prop%nsoil_types))
  ALLOCATE(emiss_b(soil_prop%nsoil_types))
  ALLOCATE(emiss_c(soil_prop%nsoil_types))
  
  ! ----------------------------------
  ! --- Loop over all horizontal gridpoints
  ! ----------------------------------
  DO i = istart, iend
    emiss_a(:) = 0.0_wp
    emiss_b(:) = 0.0_wp
    emiss_c(:) = 0.0_wp
    
    IF (soil_prop%dust_mask(i,jb) .AND. &
      & hsnow  (i) == 0.0_wp      .AND. &
      & soilice(i) == 0.0_wp) THEN
      
      ! ----------------------------------
      ! --- Calculations
      ! ----------------------------------
      ! Friction velocity
      CALL calc_ustar(ustar,tcm(i),u(i),v(i))

      ! size-dependent threshold value of the friction velocity: factors needed
      ! for Eq. 3.3
      fac33a = an33*rho_p*grav/rho(i)
      fac33b = an33*gamma33/rho(i)

      wbgrav = 100._wp*(soilwater(i)/dzsoil)*bulkdens/waterdens
      ! Calculate impact of soil humidity (Eq. 3.6)
      f_eta_wb = SQRT(1.0_wp + 1.21_wp          &
        &      * (MAX(wbgrav - soil_prop%wstrich(i,jb)*1.0_wp,0.0_wp))**0.68_wp )

      ! Factor for u* for non-ideal conditions (Eq. 3.4)
      f_eta = MAX(1.0_wp,f_eta_wb*soil_prop%f_z0(i,jb))

      ! calculate minimum value of ustart (Eq. 3.3 and 3.4)
      ustart_min = f_eta*SQRT((fac33a*dp_min) &
        &        + (fac33b/dp_min))
      
      ! collect ustar and ustart_min for diagnostic output
      IF (lart_diag_out) THEN
        ustart_min_diag(i) = ustart_min
        ustar_diag     (i) = ustar
      ENDIF

      IF (ustar > ustart_min) THEN
        ! factor which characterizes the actual dispersion of the soil
        gamma_disperse  = EXP( -0.5_wp*((ustar-ustart_min)**3) )

        DO js = 1, soil_prop%nsoil_types

          soil_tab => soil_prop%soil_type(js)

          ! We neglect emissions of soil types with a very low fraction
          IF (soil_tab%fr_soil(i,jb) > 1.e-5_wp) THEN
            ! Calculate cross sectional area of the current soil type with
            ! current gamma_disperse
            stot = gamma_disperse*soil_tab%stot_min + (1._wp-gamma_disperse)*soil_tab%stot_ful

            ! prefactor from Eq. 3.7 and enumerator of Eq. 3.8 with 3.9 without
            ! m_s and 1/d_p
            IF (stot > 0.0_wp) THEN
              cstrich = pi*cwhite*rho(i)*ustar**3 &
                &     / (4.0_wp*stot*grav)        &
                &     *  6.0_wp/(pi*rho_p)
            ELSE
              cstrich = 0.0_wp
            ENDIF

            de3     = (e3/(pi/12.0_wp*rho_p*17.0_wp**2*ustar**2))**(1.0_wp/3.0_wp)
            ustart2 = fac33a*de3 + fac33b/de3
            de2     = (e2/(pi/12.0_wp*rho_p*17.0_wp**2*ustar**2))**(1.0_wp/3.0_wp)
            dstart  = de2
            dlnd    = (LOG(dend)-LOG(dstart))/integ_steps

            ekinfac = 289._wp* (pi/12._wp)*rho_p*(ustar**2)

            ! Precomputation of quantities that do not depend on the soil mode
            DO integ = 1,integ_steps+1

              dint(integ)     = dstart*EXP(REAL(integ-1,wp)*dlnd)
              log_dint(integ) = LOG(dint(integ))

              abdint(integ) = fac33a*dint(integ)+fac33b/dint(integ)
              sqrt_abdint(integ) = SQRT(abdint(integ))
              abdint_p15(integ)  = sqrt_abdint(integ)**3

              ekin           = ekinfac*(dint(integ)**3)

              ekin_e1(integ) = MAX(0._wp,ekin - e1)
              ekin_e2(integ) = MAX(0._wp,ekin - e2)
              ekin_e3(integ) = MAX(1.e-100_wp,ekin - e3)
              e3fac(integ)   = MERGE(0._wp,1._wp,ekin-e3 <= 0._wp)
            ENDDO

            DO jn = 1, soil_tab%nsoil_modes
              ! first Peak (p3) outside of num. integration (Chap. 2.3.4 Dipl.Thesis C. Hoose)
              IF (soil_tab%fr_mode(jn,1) > 0.0_wp .AND. soil_tab%fr_mode(jn,2) > 0.0_wp) THEN
                logsigma_min = LOG(soil_tab%std_dev(jn,1))
                logdgm_min   = LOG(soil_tab%diam_med(jn,1))
                logsigma_ful = LOG(soil_tab%std_dev(jn,2))
                logdgm_ful   = LOG(soil_tab%diam_med(jn,2))

                ! Total mass density for minimal and fully dispersed
                mass_m = soil_tab%fr_mode(jn,1) &
                  &    * EXP(-(LOG(de3) - logdgm_min)**2/(2.0_wp*(logsigma_min)**2))/logsigma_min
                mass_f = soil_tab%fr_mode(jn,2) &
                  &    * EXP(-(LOG(de3) - logdgm_ful)**2/(2.0_wp*(logsigma_ful)**2))/logsigma_ful

                ! Eq. 3.7 the two brackets; the missing 1/d_p from Eq. 3.8/3.9
                ! and lognormal distribution without constant prefactors
                q3dint = MAX(0.0_wp,                            &
                  &         (1.0_wp + SQRT(ustart2)*f_eta/ustar &
                  &         - (ustart2)*f_eta**2/(ustar**2)     &
                  &         - (ustart2)**1.5_wp*f_eta**3        &
                  &         / ustar**3)/de3                     &
                  &         * (gamma_disperse*mass_m + (1.0_wp - gamma_disperse)*mass_f))

                ! dust flux of mode Mode 3
                intfactor1 = cstrich*bulkdens    &
                  &        / (SQRT(2._wp*pi))    &
                  &        * pi / 6._wp*rho_p*beta

                emiss_c(js)= emiss_c(js)                 &
                  &        + intfactor1                  &
                  &        *(ddust3**3)/e3               &
                  &        *(LOG(de2)-LOG(de3))*q3dint
              

                p1fac   = intfactor1*ddust1**3/e1
                p2fac   = intfactor1*ddust2**3/e2
                p3fac   = intfactor1*ddust3**3/e3
              
                DO integ = 1,integ_steps+1
                    
                  ! Total mass density for minimal and fully dispersed
                  mass_m = soil_tab%fr_mode(jn,1) * EXP(-(log_dint(integ)-logdgm_min)**2 & 
                    &    / (2._wp*(logsigma_min)**2)) / logsigma_min
                  mass_f = soil_tab%fr_mode(jn,2) * EXP(-(log_dint(integ)-logdgm_ful)**2 &
                    &    / (2._wp*(logsigma_ful)**2)) / logsigma_ful
                    
                  intfactor2 = (1._wp                             &
                    &        + sqrt_abdint(integ)*f_eta/(ustar)   &
                    &        - abdint(integ)*f_eta**2/(ustar**2)  &
                    &        - abdint_p15(integ)*f_eta**3         &
                    &        / (ustar**3))/dint(integ)            &
                    &        * (gamma_disperse*mass_m + (1._wp-gamma_disperse)*mass_f) &
                    &        * dlnd

                  p1 = ekin_e1(integ)/ekin_e3(integ)
                  q1dint = p1 * p1fac * intfactor2

                  p2 = ekin_e2(integ)/ekin_e3(integ) * (1._wp-p1)
                  q2dint = p2 * p2fac * intfactor2

                  p3 = e3fac(integ) -p1 -p2
                  q3dint = p3 * p3fac * intfactor2

                  emiss_a(js) = emiss_a(js) + MAX(q1dint,0._wp)
                  emiss_b(js) = emiss_b(js) + MAX(q2dint,0._wp)
                  emiss_c(js) = emiss_c(js) + MAX(q3dint,0._wp)
                ENDDO !integ_steps
              ENDIF !fr_mode
            ENDDO !nsoil_modes
          ENDIF !fr_soil
        ENDDO !nsoil_types
      ENDIF ! ustar

      soil_prop%emiss_rate_a(i,jb) = 0.0_wp
      soil_prop%emiss_rate_b(i,jb) = 0.0_wp
      soil_prop%emiss_rate_c(i,jb) = 0.0_wp
      DO js = 1, soil_prop%nsoil_types
        soil_tab => soil_prop%soil_type(js)
        soil_prop%emiss_rate_a(i,jb) = soil_prop%emiss_rate_a(i,jb) &
          &                          + emiss_a(js)*soil_tab%fr_soil(i,jb)
        soil_prop%emiss_rate_b(i,jb) = soil_prop%emiss_rate_b(i,jb) &
          &                          + emiss_b(js)*soil_tab%fr_soil(i,jb)
        soil_prop%emiss_rate_c(i,jb) = soil_prop%emiss_rate_c(i,jb) &
          &                          + emiss_c(js)*soil_tab%fr_soil(i,jb)
      ENDDO
    ELSE
      soil_prop%emiss_rate_a(i,jb) = 0.0_wp
      soil_prop%emiss_rate_b(i,jb) = 0.0_wp
      soil_prop%emiss_rate_c(i,jb) = 0.0_wp
    ENDIF !dust_mask,hsnow, soilice
  ENDDO !i

  DEALLOCATE(emiss_a,emiss_b,emiss_c)
  
END SUBROUTINE art_prepare_emission_dust
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_dust(dz, lu_frac_shrub_eg, lu_frac_shrub, lu_frac_grass, lu_frac_bare, &
  &                          lu_frac_sparse, jb, istart, iend, soil_prop, emiss_rate)
!<
! SUBROUTINE art_emission_dust
! This subroutine updates the mineral dust concentration with the according emission fluxes 
! Based on: Vogel et al. (2006) - A model of dust transport applied to the Dead Sea Area
!           COSMO-ART code by C. HOOSE and H. VOGEL
! Part of Module: mo_art_emission_dust
! Author: Daniel Rieger, KIT
! Initial Release: 2013-05-16
! Modifications:
! 2014-11-11: Daniel Rieger, KIT
! - Unified ICON/COSMO-ART emission routine.
!>
  REAL(wp),INTENT(IN)   :: &
    &  dz(:),              & !< layer thickness of bottommost layer [m]
    &  lu_frac_shrub_eg(:),& !< Shrub cover Grassland/Forest // Evergreen
    &  lu_frac_shrub(:),   & !< Closed to open Shrubland (deciduous)
    &  lu_frac_grass(:),   & !< Grassland//herbaceous
    &  lu_frac_bare(:),    & !< Landuse class fraction bare soil
    &  lu_frac_sparse(:)     !< Landuse class fraction sparse cover
  INTEGER,INTENT(IN)    :: &
    &  jb,                 & !< Block index
    &  istart, iend          !< Start and end index of nproma loop
  TYPE(t_art_soil_properties),INTENT(IN) :: &
    &  soil_prop             !< Storage container for all fields related to soil
  REAL(wp),INTENT(INOUT)    :: &
    &  emiss_rate(istart:iend,3) !< Output emission rate [ug m-3 s-1]
  ! Local variables
  INTEGER               :: &
    &  i                     !< Loop index
  
  ! ----------------------------------
  ! --- General calculations (no need to loop around)
  ! ----------------------------------
  ! Calculation of vertical flux into the lowest model layer
  ! Bestimmung der einzelnen Fluesse durch Gewichtung nach 
  ! Anteil des Bodentypes an der Gitterbox (pro)
  ! Multiplikation mit Anteil erodierbarer Flaeche (psoil)
  ! Umrechnung kg/m3 -> ug/m3
  DO i = istart, iend
!!  CASE ('dusta')
    IF (soil_prop%dust_mask(i,jb)) THEN
      emiss_rate(i,1) = soil_prop%emiss_rate_a(i,jb)/dz(i)*1.0e9_wp       &
        &               * ( lu_frac_bare(i) + 0.4_wp * lu_frac_sparse(i)  &
        &                   + 0.125_wp * (lu_frac_shrub_eg(i) + lu_frac_shrub(i) + lu_frac_grass(i)) )
      emiss_rate(i,1) = MAX(emiss_rate(i,1),0.0_wp)
    ENDIF

!!  CASE ('dustb')
    IF (soil_prop%dust_mask(i,jb)) THEN
      emiss_rate(i,2) = soil_prop%emiss_rate_b(i,jb)/dz(i)*1.0e9_wp       &
        &               * ( lu_frac_bare(i) + 0.4_wp * lu_frac_sparse(i)  &
        &                   + 0.125_wp * (lu_frac_shrub_eg(i) + lu_frac_shrub(i) + lu_frac_grass(i)) )
      emiss_rate(i,2) = MAX(emiss_rate(i,2),0.0_wp)
    ENDIF

!!  CASE ('dustc')
    IF (soil_prop%dust_mask(i,jb)) THEN
      emiss_rate(i,3) = soil_prop%emiss_rate_c(i,jb)/dz(i)*1.0e9_wp       &
        &               * ( lu_frac_bare(i) + 0.4_wp * lu_frac_sparse(i)  &
        &                   + 0.125_wp * (lu_frac_shrub_eg(i) + lu_frac_shrub(i) + lu_frac_grass(i)) )
      emiss_rate(i,3) = MAX(emiss_rate(i,3),0.0_wp)
    ENDIF
  ENDDO !i

END SUBROUTINE art_emission_dust
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_dust
