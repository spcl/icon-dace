!
! Description:
!  This module contains variables that are used in the grid scale 
!  parameterizations (Microphysics). 
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

MODULE gscp_data

!==============================================================================

USE, INTRINSIC :: iso_fortran_env, ONLY: wp => real64, i4 => int32

USE mo_math_constants    , ONLY: pi

USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 r_d   => rd    , & !! gas constant for dry air
                                 lh_s  => als   , & !! latent heat of sublimation
                                 t0    => tmelt !! melting temperature of ice/snow

USE mo_exception,          ONLY: finish, message, message_text

!==============================================================================

IMPLICIT NONE

PUBLIC

!==============================================================================

! Hardcoded switches to select autoconversion and n0s calculation
! ---------------------------------------------------------------

INTEGER,  PARAMETER ::  &
  iautocon       = 1,   &
  isnow_n0temp   = 2


! Epsilons and thresholds
! -----------------------

REAL (KIND=wp), PARAMETER ::  &
  zqmin = 1.0E-15_wp, & ! threshold for computations
  zeps  = 1.0E-15_wp    ! small number

     
REAL (KIND=wp), PARAMETER :: zxiconv = 1.0E-09_wp ! mean crystal mass for convectively generated ice (gscp3 only)

! Variables which are (mostly) initialized in gscp_set_coefficients
! -----------------------------------------------------------------

  REAL (KIND=wp)     ::           &
    ccsrim,    & !
    ccsagg,    & !
    ccsdep,    & !
    ccsvel,    & !
    ccsvxp,    & !
    ccslam,    & !
    ccslxp,    & !
    ccsaxp,    & !
    ccsdxp,    & !
    ccshi1,    & !
    ccdvtp,    & !
    ccidep,    & !
    ccswxp,    & !
    zconst,    & !
    zcev,      & !
    zbev,      & !
    zcevxp,    & !
    zbevxp,    & !
    zvzxp,     & !
    zvz0r,     & !
    vtxexp,    & !
    kc_c1,     & !
    kc_c2,     & !
    zn0r,      & ! N0_rain
    zar,       & !
    zceff_min, & ! Minimum value for sticking efficiency
    v0snow,    & ! factor in the terminal velocity for snow
    zvz0i,     & ! Terminal fall velocity of ice  (original value of Heymsfield+Donner 1990: 3.29)
    icesedi_exp,&! exponent for density correction for coud ice sedimentation
    cloud_num = 200.00e+06_wp      ! cloud droplet number concentration

! More variables
! --------------

  REAL (KIND=wp)     ::             &
    rain_n0_factor =  1.0_wp,       & ! COSMO_EU default
    mu_rain        =  0.0_wp,       & ! COSMO_EU default
    mu_snow        =  0.0_wp,       & ! COSMO_EU default
    ageo_snow                         ! global zams, will be zams_ci or zams_gr

! Even more variables (currently only used in Carmen's scheme for KC05 fall speed)
! --------------

  REAL (KIND=wp)     ::             &
    kc_alpha       =  0.5870086_wp, & !..alf, CGS is 0.00739
    kc_beta        =  2.45_wp,      & !..exponent  in mass-size relation
    kc_gamma       =  0.120285_wp,  & !..gam, CGS is 0.24
    kc_sigma       =  1.85_wp,      & !..exponent  in area-size relation
    do_i           =  5.83_wp,      & ! coefficients for drag correction
    co_i           =  0.6_wp          ! coefficients for turbulence correction


! Parameters for autoconversion of cloud water and cloud ice 
! ----------------------------------------------------------

REAL    (KIND=wp   ), PARAMETER ::  &
  zccau  = 4.0E-4_wp,    & ! autoconversion coefficient (cloud water to rain)
  zciau  = 1.0E-3_wp,    & ! autoconversion coefficient (cloud ice   to snow)

  zkcau  = 9.44e+09_wp,  & ! kernel coeff for SB2001 autoconversion
  zkcac  = 5.25e+00_wp,  & ! kernel coeff for SB2001 accretion
  zcnue  = 2.00e+00_wp,  & ! gamma exponent for cloud distribution
  zxcmin = 5.24e-13_wp,  & ! minimum mass of cloud droplets (10 micron diameter)
  zxstar = 2.60e-10_wp,  & ! separating mass between cloud and rain
  zkphi1 = 6.00e+02_wp,  & ! constant in phi-function for autoconversion
  zkphi2 = 0.68e+00_wp,  & ! exponent in phi-function for autoconversion
  zkphi3 = 5.00e-05_wp,  & ! exponent in phi-function for accretion

  zhw    = 2.270603_wp,  & ! Howell factor
  zecs   = 0.9_wp,       & ! Collection efficiency for snow collecting cloud water

  zadi   = 0.217_wp,     & ! Formfactor in the size-mass relation of ice particles
  zbdi   = 0.302_wp,     & ! Exponent in the size-mass relation of ice particles
  zams_ci= 0.069_wp,     & ! Formfactor in the mass-size relation of snow particles for cloud ice scheme
  zams_gr= 0.069_wp,     & ! Formfactor in the mass-size relation of snow particles for graupel scheme
  zbms   = 2.000_wp,     & ! Exponent in the mass-size relation of snow particles
                           ! (do not change this, exponent of 2 is hardcoded for isnow_n0temp=2)
  zv1s   = 0.50_wp,      & ! Exponent in the terminal velocity for snow

  zami   = 130.0_wp,     & ! Formfactor in the mass-size relation of cloud ice
  zn0s0  = 8.0E5_wp,     & ! 
  zn0s1  = 13.5_wp * 5.65E5_wp, & ! parameter in N0S(T)
  zn0s2  = -0.107_wp,    & ! parameter in N0S(T), Field et al
  zcac   = 1.72_wp,      & ! (15/32)*(PI**0.5)*(ECR/RHOW)*V0R*AR**(1/8)
  zcicri = 1.72_wp,      & ! (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
  zcrcri = 1.24E-3_wp,   & ! (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
  zcsmel = 1.48E-4_wp,   & ! 4*LHEAT*N0S*AS**(-2/3)/(RHO*lh_f)
  zbsmel = 20.32_wp,     & !        0.26*sqrt(    RHO*v0s/eta)*Gamma(21/8)*AS**(-5/24)
  zasmel = 2.43E3_wp,    & ! DIFF*lh_v*RHO/LHEAT

  zcrfrz = 1.68_wp,      & ! coefficient for raindrop freezing
  zcrfrz1= 9.95e-5_wp,   & !FR: 1. coefficient for immersion raindrop freezing: alpha_if
  zcrfrz2= 0.66_wp,      & !FR: 2. coefficient for immersion raindrop freezing: a_if

  zrho0  = 1.225e+0_wp,  & ! reference air density
  zrhow  = 1.000e+3_wp,  & ! density of liquid water

  zdv    = 2.22e-5_wp,   & ! molecular diffusion coefficient for water vapour
  zlheat = 2.40E-2_wp,   & ! thermal conductivity of dry air
  zeta   = 1.75e-5_wp      ! kinematic viscosity of air 


! Additional parameters
! ---------------------

REAL    (KIND=wp   ), PARAMETER ::  &
  zthet  = 248.15_wp,       & ! temperature for het. nuc. of cloud ice
  zthn   = 236.15_wp,       & ! temperature for hom. freezing of cloud water
  ztrfrz = 271.15_wp,       & ! threshold temperature for heterogeneous freezing of raindrops
  ztmix  = 250.15_wp,       & ! threshold temperature for mixed-phase cloud freezing of cloud drops (Forbes 2012)
  znimax_Thom = 250.E+3_wp, & ! FR: maximal number of ice crystals 
  zmi0   = 1.0E-12_wp,      & ! initial crystal mass for cloud ice nucleation
  zmimax = 1.0E-9_wp,       & ! maximum mass of cloud ice crystals   
  zmsmin = 3.0E-9_wp,       & ! initial mass of snow crystals        
  zbvi   = 0.16_wp,         & ! v = zvz0i*rhoqi^zbvi
!
  v_sedi_rain_min    = 0.7_wp, & ! in m/s; minimum terminal fall velocity of rain    particles (applied only near the ground)
  v_sedi_snow_min    = 0.1_wp, & ! in m/s; minimum terminal fall velocity of snow    particles (applied only near the ground)
  v_sedi_graupel_min = 0.4_wp    ! in m/s; minimum terminal fall velocity of graupel particles (applied only near the ground)


! Constant exponents in the transfer rate equations
! -------------------------------------------------

REAL    (KIND=wp   ), PARAMETER ::  &
  x1o12  =  1.0_wp/12.0_wp, & !
  x3o16  =  3.0_wp/16.0_wp, & !
  x7o8   =  7.0_wp/ 8.0_wp, & !
  x2o3   =  2.0_wp/ 3.0_wp, & !
  x5o24  =  5.0_wp/24.0_wp, & !
  x1o8   =  1.0_wp/ 8.0_wp, & !
  x13o8  = 13.0_wp/ 8.0_wp, & !
  x13o12 = 13.0_wp/12.0_wp, & !
  x27o16 = 27.0_wp/16.0_wp, & !
  x1o3   =  1.0_wp/ 3.0_wp, & !
  x1o2   =  1.0_wp/ 2.0_wp, & !
  x3o4   =  0.75_wp,        & !
  x7o4   =  7.0_wp/ 4.0_wp    !

  REAL (KIND=wp),     PARAMETER ::  &
    mma(10) = (/5.065339_wp, -0.062659_wp, -3.032362_wp, 0.029469_wp, -0.000285_wp, &
                0.312550_wp,  0.000204_wp,  0.003199_wp, 0.000000_wp, -0.015952_wp /), &
    mmb(10) = (/0.476221_wp, -0.015896_wp,  0.165977_wp, 0.007468_wp, -0.000141_wp, &
                0.060366_wp,  0.000079_wp,  0.000594_wp, 0.000000_wp, -0.003577_wp /)


! Parameters relevant to support supercooled liquid water (SLW), sticking efficiency, ...
! ---------------------------------------------------------------------------------------

REAL    (KIND=wp   ), PARAMETER ::  &
  dist_cldtop_ref  = 500.0_wp,  & ! Reference length for distance from cloud top (Forbes 2012)
  reduce_dep_ref   = 0.1_wp,    & ! lower bound on snow/ice deposition reduction
  zceff_fac        = 3.5E-3_wp, & ! Scaling factor [1/K] for temperature-dependent cloud ice sticking efficiency
  tmin_iceautoconv = 188.15_wp    ! Temperature at which cloud ice autoconversion starts



!=======================================================================

CONTAINS

!==============================================================================
!==============================================================================
!>  Module procedure "hydci_pp_init" in "gscp" to initialize some
!!  coefficients which are used in "hydci_pp"
!------------------------------------------------------------------------------

SUBROUTINE gscp_set_coefficients (igscp, idbg, tune_zceff_min, tune_v0snow, tune_zvz0i, &
  &                               tune_mu_rain, tune_rain_n0_factor, tune_icesedi_exp)

!------------------------------------------------------------------------------
!> Description:
!!   Calculates some coefficients for the microphysics schemes. 
!!   Usually called only once at model startup.
!------------------------------------------------------------------------------

  INTEGER  ,INTENT(IN)           ::  igscp

  INTEGER  ,INTENT(IN) ,OPTIONAL ::  idbg              !! debug level
  REAL(wp) ,INTENT(IN) ,OPTIONAL ::  tune_zceff_min
  REAL(wp) ,INTENT(IN) ,OPTIONAL ::  tune_v0snow
  REAL(wp) ,INTENT(IN) ,OPTIONAL ::  tune_zvz0i
  REAL(wp) ,INTENT(IN) ,OPTIONAL ::  tune_icesedi_exp
  REAL(wp) ,INTENT(IN) ,OPTIONAL ::  tune_mu_rain
  REAL(wp) ,INTENT(IN) ,OPTIONAL ::  tune_rain_n0_factor
  
! Local variable
  REAL(wp) :: zams  ! local value of zams
  
!------------------------------------------------------------------------------
!>  Initial setting of local and global variables
!------------------------------------------------------------------------------


  IF (igscp == 2) THEN
    zams = zams_gr         ! default for graupel scheme
  ELSE
    zams = zams_ci         ! default for cloud ice scheme
  END IF

  ageo_snow = zams           ! zams is local, but ageo_snow will survive

  IF (PRESENT(tune_zceff_min)) THEN
    zceff_min = tune_zceff_min
  ELSE
    zceff_min = 0.075_wp     ! COSMO default
  ENDIF

  IF (PRESENT(tune_v0snow)) THEN
    IF (tune_v0snow <= 0._wp) THEN
      IF (igscp == 2) THEN
        v0snow = 20.0_wp     ! default for graupel scheme
      ELSE
        v0snow = 25.0_wp     ! default for cloud ice scheme
      END IF
    ELSE
      v0snow = tune_v0snow   ! use ICON namelist value
    END IF
  ELSE
    v0snow = 20.0_wp         ! COSMO default
  ENDIF

  IF (PRESENT(tune_zvz0i)) THEN
    zvz0i = tune_zvz0i
  ELSE
    zvz0i = 1.25_wp          ! COSMO default
  ENDIF

  IF (PRESENT(tune_icesedi_exp)) THEN
    icesedi_exp = tune_icesedi_exp
  ELSE
    icesedi_exp = 0.33_wp
  ENDIF

  IF (PRESENT(tune_mu_rain)) THEN
    mu_rain = tune_mu_rain
  ELSE
    mu_rain = 0.0_wp         ! COSMO-EU default
  ENDIF

  IF (PRESENT(tune_rain_n0_factor)) THEN
    rain_n0_factor = tune_rain_n0_factor
  ELSE
    rain_n0_factor = 1.0_wp         ! COSMO-EU default
  ENDIF

  zconst = zkcau / (20.0_wp*zxstar) * (zcnue+2.0_wp)*(zcnue+4.0_wp)/(zcnue+1.0_wp)**2
  ccsrim = 0.25_wp*pi*zecs*v0snow*GAMMA(zv1s+3.0_wp)
  ccsagg = 0.25_wp*pi*v0snow*GAMMA(zv1s+3.0_wp)
  ccsdep = 0.26_wp*GAMMA((zv1s+5.0_wp)/2.0_wp)*SQRT(1.0_wp/zeta)
  ccsvxp = -(zv1s/(zbms+1.0_wp)+1.0_wp)
  ccsvel = zams*v0snow*GAMMA(zbms+zv1s+1.0_wp)      &
          *(zams*GAMMA(zbms+1.0_wp))**ccsvxp
  ccsvxp = ccsvxp + 1.0_wp
  ccslam = zams*GAMMA(zbms+1.0_wp)
  ccslxp = 1.0_wp / (zbms+1.0_wp)
  ccswxp = zv1s*ccslxp
  ccsaxp = -(zv1s+3.0_wp)
  ccsdxp = -(zv1s+1.0_wp)/2.0_wp
  ccshi1 = lh_s*lh_s/(zlheat*r_v)
  ccdvtp = 2.22E-5_wp * t0**(-1.94_wp) * 101325.0_wp
  ccidep = 4.0_wp * zami**(-x1o3)
  zn0r   = 8.0E6_wp * EXP(3.2_wp*mu_rain) * (0.01_wp)**(-mu_rain)  ! empirical relation adapted from Ulbrich (1983)
  zn0r   = zn0r * rain_n0_factor                                   ! apply tuning factor to zn0r variable
  zar    = pi*zrhow/6.0_wp * zn0r * GAMMA(mu_rain+4.0_wp)      ! pre-factor in lambda of rain
  zcevxp = (mu_rain+2.0_wp)/(mu_rain+4.0_wp)
  zcev   = 2.0_wp*pi*zdv/zhw*zn0r*zar**(-zcevxp) * GAMMA(mu_rain+2.0_wp)
  zbevxp = (2.0_wp*mu_rain+5.5_wp)/(2.0_wp*mu_rain+8.0_wp)-zcevxp
  zbev   =  0.26_wp * SQRT(    zrho0*130.0_wp/zeta)*zar**(-zbevxp) &
           * GAMMA((2.0_wp*mu_rain+5.5_wp)/2.0_wp) / GAMMA(mu_rain+2.0_wp)

  zvzxp  = 0.5_wp/(mu_rain+4.0_wp)
  zvz0r  = 130.0_wp*GAMMA(mu_rain+4.5_wp)/GAMMA(mu_rain+4.0_wp)*zar**(-zvzxp)

  IF (PRESENT(idbg)) THEN
    IF (idbg > 10) THEN
      CALL message('gscp_set_coefficients',': Initialized coefficients for microphysics')
      WRITE (message_text,'(A,E10.3)') '      zams   = ',zams   ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsvel = ',ccsvel ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsrim = ',ccsrim ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsagg = ',ccsagg ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsdep = ',ccsdep ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccslxp = ',ccslxp ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccidep = ',ccidep ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      mu_r   = ',mu_rain; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zn0r   = ',zn0r   ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zbev   = ',zbev   ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zbevxp = ',zbevxp ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zcev   = ',zcev   ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zcevxp = ',zcevxp ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zvzxp  = ',zvzxp  ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zvz0r  = ',zvz0r  ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '   zceff_min = ',zceff_min ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      v0snow = ',v0snow ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '       zvz0i = ',zvz0i  ; CALL message('',message_text)
    ENDIF
  ENDIF

  CALL message('gscp_set_coefficients','microphysical values initialized')

END SUBROUTINE gscp_set_coefficients

!==============================================================================

END MODULE gscp_data
