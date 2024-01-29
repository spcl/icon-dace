!
! mo_art_aerosol_utilities
! This module provides utility subroutines for modal aerosol.
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

MODULE mo_art_aerosol_utilities
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_physical_constants,            ONLY: rgrav
  USE mo_math_constants,                ONLY: pi, sqrt2
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_exception,                     ONLY: finish,message,message_text
  USE mo_var_list,                      ONLY: get_tracer_info_dyn_by_idx, &
    &                                         t_var_list_ptr
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic
  USE turb_data,                        ONLY: akt
  USE mo_fortran_tools,                 ONLY: set_acc_host_or_device
!! ART
!  USE mo_art_data,                      ONLY: akt, rakt
  
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_aerosol_utilities'

  PUBLIC :: art_air_properties
  PUBLIC :: calc_ustar
  PUBLIC :: calc_aerodynamic_resistance
  PUBLIC :: art_modal_parameters
  PUBLIC :: art_calc_number_from_mass
  PUBLIC :: art_modeshift

  INTERFACE art_calc_number_from_mass
    MODULE PROCEDURE art_calc_number_from_mass_scal
    MODULE PROCEDURE art_calc_number_from_mass_1d
    MODULE PROCEDURE art_calc_number_from_mass_2d
  END INTERFACE art_calc_number_from_mass

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_air_properties(pres,temp,istart,iend,kstart,kend,free_path,dyn_visc,lacc)
!<
! SUBROUTINE art_air_properties
! This subroutine calculates the mean free path and the dynamic
! viscosity of the air
! Based on:
! - viscosity: W. Sutherland (1893), Philosophical Magazine Series 5
! Part of Module: mo_art_aerosol_utilities
! Author: Daniel Rieger, KIT
! Initial Release: 2013-09-30
! Modifications:
! 2016-11-03: Daniel Rieger, KIT
! - Modularization and update of loop structure to general ART physics coupling concept
!>
  REAL(wp), INTENT(in)           :: &
    &  pres(:,:),                   & !< pressure (Pa)
    &  temp(:,:)                      !< temperature (K)
  INTEGER, INTENT(in)            :: &
    &  istart,iend,                 & !< Start and end index of nproma loop
    &  kstart,kend                    !< Start and end index of vertical loop
  REAL(wp),INTENT(inout)         :: &
    &  free_path(:,:),              & !< mean free path
    &  dyn_visc(:,:)                  !< dynamic viscosity
  LOGICAL, OPTIONAL, INTENT(in)  :: &
    &  lacc
  !Local variables
  REAL(wp)                       :: &
    &  p_ref,                       & !< reference pressure, 1013.25 hPa
    &  t_ref                          !< reference temperature
  INTEGER                        :: &
    &  jc, jk                         !< loop indizes
  LOGICAL                        :: & !< OpenACC flag
    &  lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  ! ----------------------------------
  ! --- set initial/reference values
  ! ----------------------------------

  p_ref     =  101325.0_wp
  t_ref     =  288.15_wp

  ! ----------------------------------
  ! --- Get loop indizes to loop over the whole domain
  ! ----------------------------------

  !$ACC DATA PRESENT(dyn_visc, free_path, pres, temp) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jk=kstart,kend
!NEC$ ivdep
    DO jc=istart,iend
! ---------------------------------
! --- Calculate mean free path [ m ]
! ----------------------------------
      free_path(jc,jk) = 6.6328e-8_wp * p_ref * temp(jc,jk) &
         &             / (t_ref * pres(jc,jk))
      ! *** 6.6328d-8 is the sea level values given in Table I.2.8
      ! *** on page 10 of U.S. Standard Atmosphere 1962
! ----------------------------------
! --- Calculate dynamic viscosity [ kg m**-1 s**-1 ]
! ----------------------------------
      dyn_visc(jc,jk) = 1.458e-6_wp * temp(jc,jk) * SQRT(temp(jc,jk))  &
         &            / (temp(jc,jk) + 110.4_wp)
      ! *** U.S. Standard Atmosphere 1962 page 14 expression
      ! *** for dynamic viscosity is:
      ! *** dynamic viscosity =  beta * T * sqrt(T) / ( T + S)
      ! *** where beta = 1.458d-6 [ kg m^-1 sec^-1 K**-0.5 ], s = 110.4 [ K ].
    ENDDO !jc
  ENDDO !jk
  !$ACC END PARALLEL

  !$ACC END DATA

END SUBROUTINE art_air_properties
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE calc_ustar(ustar,tcm,u,v)
!<
! SUBROUTINE calc_ustar
! This subroutine calculates the friction velocity (ustar)
! Based on: COSMO-ART code by Heike Vogel.
! Part of Module: mo_art_aerosol_utilities
! Author: Daniel Rieger, KIT
! Initial Release: 2014-04-17
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - ...
!>
  REAL(wp), INTENT(out) :: &
    &  ustar                 !< friction velocity
  REAL(wp), INTENT(in)  :: &
    &  tcm,                & !< transfer coefficient (momentum)
    &  u,                  & !< wind velocity
    &  v                     !< wind velocity

  ustar    = 0.0_wp

  ! ----------------------------------
  ! --- Calculation of friction velocity
  ! ----------------------------------
  IF(tcm /=0.0_wp) THEN
    ustar = SQRT(u*u + v*v)*SQRT(tcm)
  ELSE
    ustar = SQRT(u*u + v*v)*SQRT(1.0e-11_wp)
  ENDIF

END SUBROUTINE calc_ustar
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE calc_aerodynamic_resistance(ustar, theta_v,gz0, qv, qc, qr,t,t_s,zref,tch,tcm,ra)
!<
! SUBROUTINE calc_aerodynamic_resistance
! This subroutine calculates the aerodynamic resistance
! Based on: COSMO-ART code by Heike Vogel.
! Part of Module: mo_art_aerosol_utilities
! Author: Daniel Rieger, KIT
! Initial Release: 2014-04-17
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - ...
!>
  REAL (KIND=wp), INTENT(in)  :: &
    &  ustar,                    & !< friction velocity
    &  theta_v,                  & !< virtual potential temperature
    &  gz0,                      & !< surface roughness length * gravity constant
    &  qv,                       & !< water vapor mass mixing ratio
    &  qc,                       & !< cloud water mass mixing ratio
    &  qr,                       & !< rain water mass mixing ratio
    &  t,                        & !< temperature in lowest model layer
    &  t_s,                      & !< surface temperature
    &  zref,                     & !< height of layer
    &  tch,                      & !< transfer coefficient (heat)
    &  tcm                         !< transfer coefficient (momentum)

  REAL (KIND=wp), INTENT(out) :: &
    &  ra                          !< aerodynamic resistance
! local variables
  REAL (KIND=wp)  :: &
    &  theta,        &             !< potential temperature
    &  tstar,        &             !<
    &  sr1,          &             !<
    &  sr2,          &             !<
    &  lstar_inv,    &             !<
    &  z0                          !< surface roughness length
  REAL(wp) ::        &
    &  rakt                        !< 1/akt (akt=Karman constant)

  rakt = 1._wp/akt



  ! ----------------------------------
  ! --- Initializations
  ! ----------------------------------

  lstar_inv  = 0.0_wp
  z0 = gz0 * rgrav

  ! ----------------------------------
  ! --- Calculate potential temperature out of virtual potential temperature
  ! ----------------------------------

  theta = theta_v / (1+0.61*qv-(qc+qr))

  ! ----------------------------------
  ! ---
  ! ----------------------------------

  IF(tcm /= 0.0_wp) THEN
    tstar      = tch*(t-t_s)/sqrt(tcm)
  ELSE
    tstar      = tch*(t-t_s)/sqrt(1.0e-11_wp)
  ENDIF

  lstar_inv = akt * rgrav * tstar/(theta * ustar**2)

  z0 = MAX(z0,0.0001_wp)
  z0 = MIN(z0,0.1_wp)

  ! ----------------------------------
  ! --- Aerodynamic resistance for convective conditions
  ! ----------------------------------

  IF ( lstar_inv < 0.0e0_wp ) THEN
    sr1 = SQRT(1._wp-9._wp*zref*lstar_inv)
    sr2 = SQRT(1._wp-9._wp*z0*lstar_inv)
     !*** Check for ln(0), Formula from Taylor's theorem
     !*** Threshold: -9.*Z0*LS(I,J) = 0.1
     !*** For Zr = 10.*Z0
    IF((-9._wp*z0*lstar_inv) <= 0.1_wp) THEN
      ra=0.74_wp*rakt/ustar * LOG(8.2e0_wp)
    ELSE
      ra= 0.74_wp*rakt/ustar* &
       &  (LOG((sr1-1._wp)/(sr1+1._wp))-  &
       &   LOG((sr2-1._wp)/(sr2+1._wp)))
    ENDIF
  ELSE
  ! ----------------------------------
  ! --- Aerodynamic resistance for neutral and stable stratification
  ! ----------------------------------
    ra = 1.*rakt/ustar*(0.74_wp*LOG(zref/z0)+ &
     &   4.7_wp*lstar_inv*(zref-z0))
  ENDIF

END SUBROUTINE calc_aerodynamic_resistance
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_modal_parameters(rho, third_moment, total_mass, density, diameter,        &
  &                             knudsen_nr, exp_aero, dg0e, free_path,                   &
  &                             init_nmb_conc, init_mass_conc,                           &
  &                             itr0, itr3, ntr, istart, iend, kstart, nlev, tracer)
!<
! SUBROUTINE art_modal_parameters
! Calculates modal parameters and derived variables,
! log-squared of std deviation, mode mean size, Knudsen number)
! based on current values of moments for the modes.
! Part of Module: mo_art_aerosol_utilities
! based on code of F. BINKOWSKI
! modified by N. RIEMER, H. VOGEL, B. VOGEL, C. HOOSE
! rewritten for ICON-ART by D. RIEGER
! Author: Daniel Rieger, KIT
! Initial Release: 2013-09-02
! Modifications:
! 2014-12-08: Daniel Rieger, KIT
! - Modifications for unified ICON-/COSMO-ART physical parameterizations
! 2018-07-30: Lukas Muser, KIT
! - Modification of infamous "to avoid problems" statement
!>
  REAL(wp),INTENT(inout)  :: &
    &  third_moment(:,:),    & !< 3rd moment of aerosol distribution
    &  total_mass(:,:),      & !< Total mass of current aerosol mode (kg kg-1)
    &  density(:,:),         & !< Density of aerosol
    &  diameter(:,:),        & !< Current diameter of aerosol
    &  knudsen_nr(:,:)         !< Knudsen number
  REAL(wp),INTENT(in)     :: &
    &  rho(:),               & !< Density of species
    &  init_nmb_conc,        & !< Initial number concentration of current mode
    &  init_mass_conc(:),    & !< Initial mass concentration of species
    &  exp_aero,             & !< aerosol exponent
    &  dg0e,                 & !< Count median diameter at emission -> at initialization? !!!
    &  free_path(:,:)          !< Mean free path
  INTEGER, INTENT(in)     :: &
    &  itr0, itr3(:),        & !< Indices of number (0) and mass mixing ratios (3) 
                               !    in tracer container
    &  ntr,                  & !< Number of species contained in mode
    &  istart, iend,         & !< Start and end indices of nproma loop
    &  kstart, nlev            !< Number of vertical (full) levels
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:)           !< Tracer field (jc,jk,ntracer)
! Local Variables
  REAL(wp)                :: &
    &  f6dpim9,              & !< 1.0e-9 * 6/PI
    &  one_third               !< 1/3
  REAL(wp)                :: &
    &  spfac(ntr-1)            !< Species dependent factor for conversion
  INTEGER                 :: &
    &  ierror,               & !< Error return value
    &  jc, jk, jsp             !< loop indizes
  INTEGER                 :: &
    &  mask(istart:iend,kstart:nlev)

  f6dpim9   = 1.0E-9_wp * 6.0_wp / pi
  one_third = 1.0_wp/3.0_wp

  DO jsp = 1, ntr-1
    spfac(jsp)  = f6dpim9 / rho(jsp)
  ENDDO

  ! Set initial values moved here, to be nproma check safe
  total_mass(:,:)   = 0._wp
  density(:,:)      = 0._wp
  diameter(:,:)     = dg0e
  ! diameter(:,:)     = 1.E-09_wp
  knudsen_nr(:,:)   = 0._wp
  third_moment(:,:) = 0._wp

  DO jk = kstart, nlev
!NEC$ ivdep
    DO jc = istart, iend
      knudsen_nr(jc,jk)   = 2.0_wp * free_path(jc,jk) / dg0e
    END DO
  END DO

  ! Diagnose aerosol 3rd moment and total mass of mode
  DO jsp = 1, ntr-1 ! loop only over mass tracers
    DO jk = kstart, nlev
!NEC$ ivdep
      DO jc = istart, iend
!MARKER        tracer(jc,jk,itr3(jsp)) = MAX(tracer(jc,jk,itr3(jsp)),1.e-24_wp)

!      IF (tracer(jc,jk,itr0)<1._wp) tracer(jc,jk,itr0) = 1._wp

      ! Diagnose aerosol 3rd moment and total mass of mode
        IF (tracer(jc,jk,itr3(jsp))<1.e-24_wp) tracer(jc,jk,itr3(jsp)) = 1.e-24_wp
        third_moment(jc,jk) = third_moment(jc,jk) + spfac(jsp) * tracer(jc,jk,itr3(jsp))
        total_mass(jc,jk) = total_mass(jc,jk) + tracer(jc,jk,itr3(jsp))
      ENDDO
    END DO
  END DO

  DO jk = kstart, nlev
!NEC$ ivdep
    DO jc = istart, iend
      IF (tracer(jc,jk,itr0) > init_nmb_conc) THEN
        IF (third_moment(jc,jk) > 0._wp) THEN
          mask(jc,jk) = 1
        ELSE
          mask(jc,jk) = 2
        END IF
      ELSE
        IF (third_moment(jc,jk) >= 10.e-25_wp) THEN
          mask(jc,jk) = 3
        ELSE
          mask(jc,jk) = 4
        END IF
      END IF
    END DO
  END DO
  

  DO jk = kstart, nlev
!NEC$ ivdep
    DO jc =  istart, iend
      SELECT CASE(mask(jc,jk))
        CASE(1)
          ! only do this if the number concentration is above initial number concentration
          ! 3rd moment is > 0

          ! Diagnose density and diameter of mode
          density(jc,jk)  = f6dpim9 * total_mass(jc,jk) / third_moment(jc,jk)
          diameter(jc,jk) = (third_moment(jc,jk)          &
            &             / (tracer(jc,jk,itr0) * (exp_aero**36) ) )**one_third

          CALL art_check_diameters(diameter(jc,jk), dg0e)   ! !!! should be discussed

          ! Calculate Knudsen Number
          knudsen_nr(jc,jk)   = 2.0_wp * free_path(jc,jk) / diameter(jc,jk)

        CASE(2)
          ! only do this if the number concentration is above initial number concentration
          ! 3rd moment is <= 0
          third_moment(jc,jk)  = init_nmb_conc * dg0e**3 * exp_aero**36 ! 0.0_wp
          tracer(jc,jk,itr0)   = init_nmb_conc
        CASE DEFAULT
          ! Nothing
      END SELECT
    END DO
  END DO


  DO jsp = 1, ntr-1
    DO jk = kstart, nlev
!NEC$ ivdep
      DO jc =  istart, iend
        SELECT CASE(mask(jc,jk))
          CASE(2)
            ! only do this if the number concentration is above initial number concentration
            ! 3rd moment is <= 0
            tracer(jc,jk,itr3(jsp)) = init_mass_conc(jsp)
            total_mass(jc,jk) = total_mass(jc,jk) + tracer(jc,jk,itr3(jsp))
          CASE DEFAULT
            ! Nothing
        END SELECT
      END DO
    END DO
  END DO

  DO jk = kstart, nlev
!NEC$ ivdep
    DO jc =  istart, iend
      SELECT CASE(mask(jc,jk))
        CASE(2)
          ! only do this if the number concentration is above initial number concentration
          ! 3rd moment is <= 0
          density(jc,jk)  = f6dpim9 * total_mass(jc,jk) / third_moment(jc,jk)
        CASE(3) 
          ! This is the case if the number conc. is lower than init_nmb_conc
          ! is 3rd moment over 10e-25? -> then precision is high enough to calculate
          ! number density !!! should be discussed
          tracer(jc,jk,itr0) = third_moment(jc,jk) / ((dg0e**3)*(exp_aero**36))
          density(jc,jk)  = f6dpim9 * total_mass(jc,jk) / third_moment(jc,jk)
        CASE(4)
          ! This is the case if the number conc. is lower than init_nmb_conc
          ! 3rd moment below 10e-25 -> we set concentrations back to initial values in order to
          ! avoid problems
          third_moment(jc,jk)  = init_nmb_conc * dg0e**3 * exp_aero**36 ! 0.0_wp
          tracer(jc,jk,itr0)   = init_nmb_conc
        CASE DEFAULT
          !Nothing
      END SELECT
    END DO
  END DO

  DO jsp = 1, ntr-1
    DO jk = kstart, nlev
!NEC$ ivdep
      DO jc =  istart, iend
        SELECT CASE(mask(jc,jk))
          CASE(4)
            ! This is the case if the number conc. is lower than init_nmb_conc
            ! 3rd moment below 10e-25 -> we set concentrations back to initial values in order to
            ! avoid problems
            tracer(jc,jk,itr3(jsp)) = init_mass_conc(jsp)
            total_mass(jc,jk) = total_mass(jc,jk) + tracer(jc,jk,itr3(jsp))
          CASE DEFAULT
            ! Nothing
        END SELECT
      END DO
    END DO
  END DO

  DO jk = kstart, nlev
!NEC$ ivdep
    DO jc =  istart, iend
      SELECT CASE(mask(jc,jk))
        CASE(4)
          ! This is the case if the number conc. is lower than init_nmb_conc
          ! 3rd moment below 10e-25 -> we set concentrations back to initial values in order to
          ! avoid problems
          density(jc,jk)  = f6dpim9 * total_mass(jc,jk) / third_moment(jc,jk)
        CASE DEFAULT
          ! Nothing
      END SELECT
    ENDDO !jc
  ENDDO !jk

END SUBROUTINE art_modal_parameters
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_check_diameters(p_diam,p_idia)
!<
! SUBROUTINE art_check_diameters
! This subroutine checks the diameters and - if thresholds
! are exceeded - sets them to given values
! Part of Module: mo_art_aerosol_utilities
! Based on: COSMO-ART code
! Author: Daniel Rieger, KIT
! Initial Release: 2013-09-02
! Modifications:
! 2017-02-28: Daniel Rieger, KIT
! - moved to mo_art_aerosol_utilities
!>
  REAL(wp),INTENT(inout)      :: &
    &  p_diam
  REAL(wp),INTENT(in)         :: &
    &  p_idia

!  IF (    TRIM(mode_name) == 'AC'  &       ! Accumulation Mode
!   & .OR. TRIM(mode_name) == 'ACm' &       ! Mixed Accumulation Mode
!   & .OR. TRIM(mode_name) == 'NU'  &       ! Nucleation Mode
!   & .OR. TRIM(mode_name) == 'NUm' &       ! Mixed Nucleation Mode
!   & .OR. TRIM(mode_name) == 'SOOT'&       ! Soot Mode
!   & .OR. TRIM(mode_name) == 'COAR') THEN  ! Coarse Mode
!    p_diam   = MIN(p_diam,   100.e-06_wp)
!  ELSE !Dust, Seasalt
    p_diam   = MIN(p_diam,   18.e-06_wp) !This check might need modification for SEASALT
    p_diam   = MAX(p_diam,    1.e-09_wp)
!  ENDIF

!  IF (    TRIM(mode_name) == 'ACm' &       ! Mixed Accumulation Mode
!   & .OR. TRIM(mode_name) == 'NUm') THEN   ! Mixed Nucleation Mode
!    IF (p_diam<p_idia) p_diam =  p_idia
!  ENDIF

END SUBROUTINE art_check_diameters
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_number_from_mass_scal(dg3, sigmag, this_list, itr3, &
  &                                       mass, number)
!<
! SUBROUTINE art_calc_number_from_mass_scal
! This subroutine calculates a number concentration, number rate or
! specific number from a mass concentration, mass rate or mass mixing
! ratio
! Part of Module: mo_art_aerosol_utilities
! Based on: COSMO-ART code
! Author: Daniel Rieger, KIT
!         Lukas Muser, KIT
! Initial Release: 2020-06-08
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! -
!>
  TYPE(t_var_list_ptr), INTENT(in) :: &
    &  this_list               !< Prognostic var_list
  REAL(wp), INTENT(in)    :: &
    &  dg3, sigmag,          & !< 3rd moment median diameter (m) and standard deviation
    &  mass                    !< Mass to be converted from
  INTEGER, INTENT(in)     :: &
    &  itr3                    !< Index of mass tracer to be converted from in container
  REAL(wp), INTENT(out)   :: &
    &  number                  !< Number to be converted to
! Local variables
  TYPE(t_var_metadata_dynamic)      :: &
    &  info_dyn                !< tracer metadata
  REAL(wp)                :: &
    &  factnum,              & !< Factor to convert 3rd moment to 0th moment
    &  fac,                  & !< Factor to convert 3rd moment to mass
    &  rho                     !< Density of particle
  INTEGER                 :: &
    &  ierror                  !< Error return value


  CALL get_tracer_info_dyn_by_idx(this_list, itr3, info_dyn)
  CALL info_dyn%tracer%opt_meta%get('rho', rho, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_calc_number_from_mass_scal',          &
      &                                'rho not available for tracer '                            &
      &                                //TRIM(info_dyn%tracer%name)//'.')

  factnum = EXP( 4.5_wp * (LOG( sigmag )**2)) / (dg3 ** 3)
  fac     = 1.0E-9_wp * 6.0_wp / pi / rho

  number  = factnum * fac * mass

END SUBROUTINE art_calc_number_from_mass_scal
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_number_from_mass_1d(dg3, sigmag, this_list, itr3, &
  &                                     istart, iend, mass, number)
!<
! SUBROUTINE art_calc_number_from_mass_1d
! This subroutine calculates a number concentration, number rate or
! specific number from a mass concentration, mass rate or mass mixing
! ratio
! Part of Module: mo_art_aerosol_utilities
! Based on: COSMO-ART code
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-20
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! -
!>
  TYPE(t_var_list_ptr), INTENT(in) :: &
    &  this_list               !< Prognostic var_list
  REAL(wp), INTENT(in)    :: &
    &  dg3, sigmag,          & !< 3rd moment median diameter (m) and standard deviation
    &  mass(:)                 !< Mass to be converted from
  INTEGER, INTENT(in)     :: &
    &  itr3,                 & !< Index of mass tracer to be converted from in container
    &  istart, iend            !< Start and end index of loop
  REAL(wp), INTENT(out)   :: &
    &  number(:)               !< Number to be converted to
! Local variables
  TYPE(t_var_metadata_dynamic)      :: &
    &  info_dyn                !< tracer metadata
  REAL(wp)                :: &
    &  factnum,              & !< Factor to convert 3rd moment to 0th moment
    &  fac,                  & !< Factor to convert 3rd moment to mass
    &  rho                     !< Density of particle
  INTEGER                 :: &
    &  ierror,               & !< Error return value
    &  jc                      !< loop index


  CALL get_tracer_info_dyn_by_idx(this_list, itr3, info_dyn)
  CALL info_dyn%tracer%opt_meta%get('rho', rho, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_calc_number_from_mass_1d', &
      &                         'rho not available for tracer '//TRIM(info_dyn%tracer%name)//'.')

  factnum = EXP( 4.5_wp * (LOG( sigmag )**2)) / (dg3 ** 3)
  fac     = 1.0E-9_wp * 6.0_wp / pi / rho

  DO jc = istart, iend
    number(jc) = factnum * fac * mass(jc)
  ENDDO

END SUBROUTINE art_calc_number_from_mass_1d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_number_from_mass_2d(dg3, sigmag, this_list, itr3, &
  &                                     istart, iend, nlev, mass, number)
!<
! SUBROUTINE art_calc_number_from_mass_2d
! This subroutine calculates a number concentration, number rate or
! specific number from a mass concentration, mass rate or mass mixing
! ratio
! Part of Module: mo_art_aerosol_utilities
! Based on: COSMO-ART code
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-20
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! -
!>
  TYPE(t_var_list_ptr), INTENT(in) :: &
    &  this_list               !< Prognostic var_list
  REAL(wp), INTENT(in)    :: &
    &  dg3, sigmag,          & !< 3rd moment median diameter (m) and standard deviation
    &  mass(:,:)               !< Mass to be converted from
  INTEGER, INTENT(in)     :: &
    &  itr3,                 & !< Index of mass tracer to be converted from in container
    &  istart, iend, nlev      !< Start and end index of loop, number of vertical levels
  REAL(wp), INTENT(out)   :: &
    &  number(:,:)             !< Number to be converted to
! Local variables
  TYPE(t_var_metadata_dynamic) :: &
    &  info_dyn                !< tracer metadata
  REAL(wp)                :: &
    &  factnum,              & !< Factor to convert 3rd moment to 0th moment
    &  fac,                  & !< Factor to convert 3rd moment to mass
    &  rho                     !< Density of particle
  INTEGER                 :: &
    &  ierror,               & !< Error return value
    &  jc, jk                  !< loop index
  
  CALL get_tracer_info_dyn_by_idx(this_list, itr3, info_dyn)
  CALL info_dyn%tracer%opt_meta%get('rho', rho, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_calc_number_from_mass_1d', &
      &                         'rho not available for tracer '//TRIM(info_dyn%tracer%name)//'.')

  factnum = EXP( 4.5_wp * (LOG( sigmag )**2)) / (dg3 ** 3)
  fac     = 1.0E-9_wp * 6.0_wp / pi / rho

  DO jc = istart, iend
    DO jk = 1, nlev
      number(jc,jk) = factnum * fac * mass(jc,jk)
    ENDDO
  ENDDO

END SUBROUTINE art_calc_number_from_mass_2d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_modeshift(sigma_s, sigma_l, diameter_s, diameter_l, &
  &                      diameter_thresh, itr0, itr3, njsp,       &
  &                      kstart, nlev, istart, iend, tracer)
!<
! SUBROUTINE art_modeshift
! Generic routine for performing modeshift
! Part of Module: mo_art_aerosol_utilities
! based on code of F. BINKOWSKI
! modified by N. RIEMER, H. VOGEL, B. VOGEL, C. HOOSE
! rewritten for ICON-ART by D. RIEGER
! Author: Heike Vogel, KIT
! Initial Release: ****-**-**
! Modifications:
! 2017-27-04: Simon, Gruber, KIT
! - Modifications for unified ICON-/COSMO-ART physical parameterizations
!>
  REAL(wp),INTENT(in)     :: &
    &  sigma_s,              & !< Standard deviation of smaller mode
    &  sigma_l,              & !< Standard deviation of larger mode
    &  diameter_s(:,:),      & !< Diameter of smaller mode
    &  diameter_l(:,:),      & !< Diameter of larger mode
    &  diameter_thresh         !< Threshold diameter for shifting
  INTEGER, INTENT(in)     :: &
    &  itr0(:),              & !< Index pair number (dim1=2: shift from/to)
    &  itr3(:,:),            & !< Index pairs mass (dim1=2: shift from/to, dim2=njsp3)
    &  njsp,                 & !< Number of species for mass mixing ratios
    &  istart, iend,         & !< Start and end indices of nproma loop
    &  kstart, nlev            !< Number of vertical (full) levels
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:)           !< Tracer field (jc,jk,ntracer)
! Local Variables
  REAL(wp)                :: &
    &  phnum,                & !< fraction of number dist. with diameter < d_thresh
    &  fnum,                 & !< fraction of number dist. with diameter > d_thresh
    &  phm3,                 & !< fraction of mass dist. with diameter < d_thresh
    &  fm3,                  & !< fraction of mass dist. with diameter > d_thresh
    &  num,                  & !< number in intersection
    &  xm3,                  & !< dummy
    &  m3,                   & !< dummy
    &  aaa,                  & !< dummy
    &  log_sig_s,            & !< LOG(sigma_s)
    &  log_sig_l               !< LOG(sigma_l)
  INTEGER                 :: &
    &  jc, jk, jsp             !< loop indizes

!!TO DO:
!! ! aitken -> accum
!! IF (growth_rate_of_3rd_moment_aitken > growth_rate_of_3rd_moment_accum &
!!
!! ! mixed aitken -> mixed accum
!! IF (growth_rate_of_3rd_moment_aitken > growth_rate_of_3rd_moment_accum &

  log_sig_s = LOG(sigma_s)
  log_sig_l = LOG(sigma_l)
  xm3       = 3.0_wp*log_sig_s/sqrt2 ! Factor used in error function call below

  DO jk = kstart, nlev
!NEC$ ivdep
    DO jc = istart, iend
      ! Check whether shifting necessary
      IF (diameter_s(jc,jk) > diameter_thresh &
      & .AND. tracer(jc,jk,itr0(1)) > tracer(jc,jk,itr0(2))) THEN
        ! aaa is the value of ln(dd/DGNUC)/(SQRT2*XXLSGN), where
        ! dd is the diameter at which the smaller and the larger mode number
        ! distributions intersect (overlap).
        aaa = getaf(tracer(jc,jk,itr0(1)), tracer(jc,jk,itr0(2)), &
          &         diameter_s(jc,jk), diameter_l(jc,jk),         &
          &         log_sig_s, log_sig_l, sqrt2)
        ! Do not let num become negative because this means that no more than one half
        ! of the total mode number may be transferred per call.
        num = MAX(aaa,xm3)
        m3  = num - xm3 ! set up for 3rd moment and mass transfers

        IF(m3 > 0.0_wp) THEN
          ! fnum and fm3 are the fractions of the number and 3rd moment
          ! distributions with  diameters greater than dd respectively.
          ! phnum and phm3 are the fractions of the number and 3rd moment
          ! distributions with diameters less than dd.
          phnum = 0.5_wp*(1.0_wp + SQRT(1.0_wp - EXP(-4.0_wp*num*num/pi)))
          fnum  = 0.5_wp*(1.0_wp - SQRT(1.0_wp - EXP(-4.0_wp*num*num/pi)))
          phm3  = 0.5_wp*(1.0_wp + SQRT(1.0_wp - EXP(-4.0_wp*m3*m3/pi)))
          fm3   = 0.5_wp*(1.0_wp - SQRT(1.0_wp - EXP(-4.0_wp*m3*m3/pi)))
          ! Update number tracers
          tracer(jc,jk,itr0(2)) = tracer(jc,jk,itr0(2)) + tracer(jc,jk,itr0(1))*fnum
          tracer(jc,jk,itr0(1)) = tracer(jc,jk,itr0(1))*phnum
          ! Update mass tracers
          DO jsp = 1, njsp
            tracer(jc,jk,itr3(2,jsp)) = tracer(jc,jk,itr3(2,jsp)) + tracer(jc,jk,itr3(1,jsp))*fm3
            tracer(jc,jk,itr3(1,jsp)) = tracer(jc,jk,itr3(1,jsp))*phm3
          ENDDO ! jsp
        ENDIF ! m3 > 0.0_wp
      ENDIF ! dia_s > thresh .AND. number conc. small > number conc. large
    ENDDO ! jc
  ENDDO ! jk

END SUBROUTINE art_modeshift
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION getaf(ni,nj,dgni,dgnj,xlsgi,xlsgj,sqrt2)
!<
! FUNCTION getaf
! Function for set up new processor for renaming
! of particles from i to j modes
! Part of Module: mo_art_aerosol_utilities
! based on code of F. BINKOWSKI
! modified by N. RIEMER, H. VOGEL, B. VOGEL, C. HOOSE
! rewritten for ICON-ART by D. RIEGER
! Author: Heike Vogel, KIT
! Initial Release: ****-**-**
! Modifications:
! 2017-27-04: Simon, Gruber, KIT
! - Modifications for unified ICON-/COSMO-ART physical parameterizations
!>
  REAL(wp), INTENT(in)    :: & !<
    &  ni, nj,               & !<
    &  dgni, dgnj,           & !<
    &  xlsgi, xlsgj,         & !<
    &  sqrt2                   !<
! Local Variables
  REAL(wp)                :: & !<
    &  aa, bb, cc,           & !<
    &  disc, qq,             & !<
    &  alfa, l, yji            !<

    alfa = xlsgi/xlsgj
    yji  = LOG(dgnj/dgni)/(sqrt2*xlsgi)
    aa   = 1.0_wp - alfa*alfa
    l    = LOG(alfa*nj/ni)
    bb   = 2.0_wp*yji*alfa*alfa
    cc   = l - yji*yji*alfa*alfa
    disc = bb*bb - 4.0_wp*aa*cc
    IF( disc < 0.0_wp) THEN
      getaf = - 5.0_wp    ! error in intersection
      RETURN
    ENDIF
    qq    = -0.5_wp*(bb + SIGN(1.0e0_wp,bb)*SQRT(disc))
    getaf = cc/qq

END FUNCTION getaf
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_aerosol_utilities
