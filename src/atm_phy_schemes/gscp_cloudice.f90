!
! Description of *gscp_cloudice*:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!   processes related to the formation of grid scale precipitation. This
!   includes the sedimentation of rain and snow. The precipitation fluxes at
!   the surface are also calculated here.
!
! Method:
!   Prognostic one-moment bulk microphysical parameterization.
!   The sedimentation of rain and snow is computed implicitly.
!
! Reference   This is an adaption of subroutine cloudice in file src_gscp.f90
!  of the COSMO-Model. Equation numbers refer to
!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!
!------------------------------------------------------------------------------
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

MODULE gscp_cloudice

!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------

USE mo_kind,               ONLY: wp         , &
                                 i4
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 lh_v  => alv   , & !! latent heat of vapourization
                                 lh_s  => als   , & !! latent heat of sublimation
!                                lh_f  => alf   , & !! latent heat of fusion
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 b3    => tmelt , & !! melting temperature of ice/snow
                                 t0    => tmelt     !! melting temperature of ice/snow

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b4w   => c4les     !!               -- " --
USE mo_satad,              ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
                                 sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
                                 latent_heat_vaporization, &
                                 latent_heat_sublimation
USE mo_exception,          ONLY: message, message_text
USE mo_run_config,         ONLY: ldass_lhn

!------------------------------------------------------------------------------

USE gscp_data, ONLY: &          ! all variables are used here

    ccsrim,    ccsagg,    ccsdep,    ccsvel,    ccsvxp,    ccslam,       &
    ccslxp,    ccsaxp,    ccsdxp,    ccshi1,    ccdvtp,    ccidep,       &
    ccswxp,    zconst,    zcev,      zbev,      zcevxp,    zbevxp,       &
    zvzxp,     zvz0r,                                                    &
    v0snow,                                                              &
    x13o8,     x1o2,      x27o16,    x7o4,      x7o8,                    &
    zbvi,      zcac,      zccau,     zciau,     zcicri,                  &
    zcrcri,    zcrfrz,    zcrfrz1,   zcrfrz2,   zeps,      zkcac,        &
    zkphi1,    zkphi2,    zkphi3,    zmi0,      zmimax,    zmsmin,       &
    zn0s0,     zn0s1,     zn0s2,     znimax_thom,          zqmin,        &
    zrho0,     zthet,     zthn,      ztmix,     ztrfrz,                  &
    zvz0i,     x2o3,      x5o24,     zams => zams_ci, zasmel,            &
    zbsmel,    zcsmel,    icesedi_exp,                                   &
    iautocon,  isnow_n0temp, dist_cldtop_ref,   reduce_dep_ref,          &
    tmin_iceautoconv,     zceff_fac, zceff_min,                          &
    mma, mmb, v_sedi_rain_min, v_sedi_snow_min

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: cloudice

LOGICAL, PARAMETER :: &
  lorig_icon   = .TRUE.  , &  ! switch for original ICON setup (only for cloudice)
  lred_depgrowth = .TRUE., &  ! switch for reduced depositional growth near tops of stratus clouds
                              ! (combined with increased 'ztmix' parameter in order not to degrade T2M in Siberian winter)

  lsedi_ice    = .FALSE. , &  ! switch for sedimentation of cloud ice (Heymsfield & Donner 1990 *1/3)
  lstickeff    = .FALSE. , &  ! switch for sticking coeff. (work from Guenther Zaengl)
  lsuper_coolw = .TRUE.       ! switch for supercooled liquid water (work from Felix Rieper)

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------


!==============================================================================

CONTAINS

#ifdef _OPENACC
! GPU code can't flush to zero double precision denormals
! So to avoid CPU-GPU differences we'll do it manually
FUNCTION make_normalized(v)
  !$ACC ROUTINE SEQ
  REAL(wp) :: v, make_normalized

  IF (ABS(v) <= 2.225073858507201e-308_wp) THEN
    make_normalized = 0.0_wp
  ELSE
    make_normalized = v
  END IF
END FUNCTION
#else
FUNCTION make_normalized(v)
  REAL(wp) :: v, make_normalized
    make_normalized = v
END FUNCTION
#endif

!==============================================================================
!> Module procedure "cloudice" in "gscp_cloudice" for computing effects of 
!!  grid scale precipitation including cloud water, cloud ice, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE cloudice (                &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,qnc,        & !! prognostic variables
  !xxx: this should become a module variable, e.g. in a new module mo_gscp_data.f90
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
  prr_gsp,prs_gsp,pri_gsp,           & !! surface precipitation rates
  qrsflux,                           & !  total precipitation flux
  l_cv,                              &
  ithermo_water,                     & !  water thermodynamics
  ldiag_ttend,     ldiag_qtend     , &
  ddt_tend_t     , ddt_tend_qv     , &
  ddt_tend_qc    , ddt_tend_qi     , & !> ddt_tend_xx are tendencies
  ddt_tend_qr    , ddt_tend_qs     , & !!    necessary for dynamics
  ddt_diag_au    , ddt_diag_ac     , & !!
  ddt_diag_ev    , ddt_diag_nuc    , & !! ddt_diag_xxx are optional
  ddt_diag_idep  , ddt_diag_sdep   , & !!   diagnostic tendencies of all
  ddt_diag_agg   , ddt_diag_rim    , & !!   individual microphysical processes
  ddt_diag_rcri  , ddt_diag_icri   , & !!
  ddt_diag_dau   , ddt_diag_iau    , & !!
  ddt_diag_imelt , ddt_diag_smelt  , & !!
  ddt_diag_cfrz  , ddt_diag_rfrz   , & !!
  ddt_diag_shed                      ) !!

!------------------------------------------------------------------------------
! Description:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, cloud ice, water vapor, rain and snow due to cloud
!   microphysical processes related to the formation of grid scale
!   precipitation.
!   The variables are updated in this subroutine. Rain and snow are
!   prognostic variables. The precipitation fluxes at the surface are stored
!   on the corresponding global fields.
!
! Method:
!   The sedimentation of rain and snow is computed implicitly.
!
! Vectorization:
!   Most computations in this routine are grouped in IF-clauses. But the IFs
!   inside DO-loops often hinder or even disables vectorization. 
!   For the big IF-chunks, the condition is now checked at the beginning of 
!   the subroutine and the corresponding indices are stored in extra index
!   arrays. The following loops then are running only over these indices,
!   avoiding at least some IF-clauses inside the loops.
!
!------------------------------------------------------------------------------
!! Declarations:
!!
!------------------------------------------------------------------------------
!! Modules used: These are declared in the module declaration section
!! -------------

!! Subroutine arguments:
!! --------------------

  INTEGER, INTENT(IN) ::  &
    nvec          ,    & !> number of horizontal points
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    ivstart   ,    & !> optional start index for horizontal direction
    ivend     ,    & !! optional end index   for horizontal direction
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=wp), INTENT(IN) :: &
    zdt             ,    & !> time step for integration of microphysics     (  s  )
    qi0,qc0                !> cloud ice/water threshold for autoconversion

  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

  INTEGER, INTENT(IN), OPTIONAL :: &
    ithermo_water          !! water thermodynamics

  LOGICAL, INTENT(IN), OPTIONAL :: &
    ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
    ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs                     !! specific snow content                         (kg/kg)

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
       qrsflux        ! total precipitation flux (nudg)

  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
    pri_gsp,             & !! precipitation rate of cloud ice, grid-scale   (kg/(m2*s))
    qnc                    !! cloud number concentration

  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs         !! tendency qs                                      ( 1/s )

  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
    ddt_diag_au     , & !> optional output autoconversion rate cloud to rain           ( 1/s )
    ddt_diag_ac     , & !! optional output accretion rate cloud to rain                ( 1/s )
    ddt_diag_ev     , & !! optional output evaporation of rain                         ( 1/s )
    ddt_diag_nuc    , & !! optional output mass nucleation of cloud ice                ( 1/s )
    ddt_diag_idep   , & !! optional output depositional growth of cloud ice            ( 1/s )
    ddt_diag_sdep   , & !! optional output depositional growth of snow                 ( 1/s )
    ddt_diag_agg    , & !! optional output aggregation snow collecting cloud ice       ( 1/s )
    ddt_diag_rim    , & !! optional output riming of snow by cloud water               ( 1/s )
    ddt_diag_rcri   , & !! optional output cloud ice + rain -> snow (rcri is sink qr)  ( 1/s )
    ddt_diag_icri   , & !! optional output cloud ice + rain -> snow (icri is sink qi)  ( 1/s )
    ddt_diag_dau    , & !! optional output depositional cloud ice autoconversion       ( 1/s )
    ddt_diag_iau    , & !! optional output aggregational cloud ice autoconversion      ( 1/s )
    ddt_diag_imelt  , & !! optional output melting of cloud ice                        ( 1/s )
    ddt_diag_smelt  , & !! optional output melting of snow                             ( 1/s )
    ddt_diag_cfrz   , & !! optional output freezing of cloud water                     ( 1/s )
    ddt_diag_rfrz   , & !! optional output rainwater freezing                          ( 1/s )
    ddt_diag_shed       !! optional output shedding                                    ( 1/s )

  !! Local parameters: None, parameters are in module header, gscp_data or data_constants
  !! ----------------
  
  !> Local scalars:
  !! -------------
  
  INTEGER (KIND=i4   ) ::  &
    iv, k             !> loop indices

  REAL    (KIND=wp   ) :: nnr

  REAL    (KIND=wp   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL    (KIND=wp   ) ::  &
    fpvsw,             & ! name of statement function
    fxna ,             & ! statement function for ice crystal number
    fxna_cooper ,      & ! statement function for ice crystal number, Cooper(1986) 
    ztx  ,             & ! dummy argument for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    znimix,            & ! number of ice crystals at ztmix -> threshold temp for mixed-phase clouds
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zdtr ,             & ! reciprocal of timestep for integration
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    zsssum, zssmax,    & ! terms for limiting snow depletion
    znin,              & ! number of cloud ice crystals at nucleation
    fnuc,              & !FR: coefficient needed for Forbes (2012) SLW layer parameterization 
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,& ! terms for limiting total cloud ice depletion
    zqvsw,             & ! sat. specitic humidity at ice and water saturation
    zqvsidiff,         & ! qv-zqvsi
    ztfrzdiff,         & ! ztrfrz-t  
    zztau, zxfac, zx1, zx2, ztt, &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc, zeff, zlog_10

  REAL    (KIND=wp   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst      ! layer tendency of snow

  REAL    (KIND=wp   ) ::  &
    tg,ppg,rhog,qvg,qcg,qig,qrg,qsg,     & ! copies of the respective state variable for each cell
    temp_c,                              & ! temperature in deg. Cesius
    zlnqrk,zlnqsk,zlnqik,zlnlogmi,       &
    ccswxp_ln1o2,zvzxp_ln1o2,zbvi_ln1o2, &
    alf,bet,m2s,m3s,hlp,maxevap

  LOGICAL :: &
    llqr,llqs,llqc,llqi  !   switch for existence of qr, qs, qc, qi

  LOGICAL :: lldiag_ttend, lldiag_qtend

  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in                     !! specific snow content                         (kg/kg)


!! Local (automatic) arrays:
!! -------------------------

  REAL    (KIND=wp   ) ::   &
    zqvsi             ,     & !> sat. specitic humidity at ice and water saturation
    zqvsw_up    (nvec),     & ! sat. specitic humidity at ice and water saturation in previous layer
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zvzi        (nvec),     & ! terminal fall velocity of ice
    zpkr        (nvec),     & ! precipitation flux of rain
    zpks        (nvec),     & ! precipitation flux of snow
    zpki        (nvec),     & ! precipitation flux of ice
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
    zprvi       (nvec)


  REAL    (KIND=wp   ) ::   &
    zcsdep            ,     & !
    zcidep            ,     & !
    zvz0s             ,     & !
    zcrim             ,     & !
    zcagg             ,     & !
    zbsdep            ,     & !
    zcslam            ,     & !
    zn0s              ,     & !
    zimr              ,     & !
    zims              ,     & !
    zimi              ,     & !
    zzar              ,     & !
    zzas              ,     & !
    zzai              ,     & !
    zqrk              ,     & ! rain water per volume of air (kg/m3)
    zqsk              ,     & ! snow water per volume of air (kg/m3)
    zqik              ,     & ! ice water per volume of air (kg/m3)
    zdtdh             ,     & !
    z1orhog           ,     & ! 1/rhog
    zrho1o2           ,     & ! (rho0/rhog)**1/2
    zrhofac_qi        ,     & ! (rho0/rhog)**icesedi_exp
    zeln7o8qrk        ,     & !
    zeln7o4qrk        ,     & ! FR new
    zeln27o16qrk      ,     & !
    zeln13o8qrk       ,     & !
    zeln5o24qsk       ,     & !
    zeln2o3qsk          


  REAL    (KIND=wp   ) ::  &
    scau   , & ! transfer rate due to autoconversion of cloud water
    scac   , & ! transfer rate due to accretion of cloud water
    snuc   , & ! transfer rate due nucleation of cloud ice
    scfrz  , & ! transfer rate due homogeneous freezing of cloud water
    simelt , & ! transfer rate due melting of cloud ice
    sidep  , & ! transfer rate due depositional growth of cloud ice
    ssdep  , & ! transfer rate due depositional growth of snow
    sdau   , & ! transfer rate due depositional cloud ice autoconversion
    srim   , & ! transfer rate due riming of snow
    sshed  , & ! transfer rate due shedding
    sicri  , & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  , & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   , & ! transfer rate due aggregation of snow and cloud ice
    siau   , & ! transfer rate due autoconversion of cloud ice
    ssmelt , & ! transfer rate due melting of snow
    sev    , & ! transfer rate due evaporation of rain
    srfrz  , & ! transfer rate due to rainwater freezing
    reduce_dep,&!FR: coefficient: reduce deposition at cloud top (Forbes 2012)
    dist_cldtop(nvec) !FR: distance from cloud top layer

#ifdef __LOOP_EXCHANGE
   REAL (KIND = wp )  ::  zlhv(ke), zlhs(ke)
#else
   REAL (KIND = wp )  ::  zlhv(nvec), zlhs(nvec) ! Latent heat if vaporization and sublimation
#endif

  LOGICAL :: lvariable_lh   ! Use constant latent heat (default .true.)

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine cloudice
!------------------------------------------------------------------------------

!> Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )

! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_wp * EXP(0.2_wp * (t0 - ztx))
  fxna_cooper(ztx) = 5.0E+0_wp * EXP(0.304_wp * (t0 - ztx))   ! FR: Cooper (1986) used by Greg Thompson(2008)

! Define reciprocal of heat capacity of dry air (at constant pressure vs at constant volume)

  IF (PRESENT(l_cv)) THEN
    IF (l_cv) THEN
      z_heat_cap_r = cvdr
    ELSE
      z_heat_cap_r = cpdr
    ENDIF
  ELSE
    z_heat_cap_r = cpdr
  ENDIF

  IF (PRESENT(ithermo_water)) THEN
     lvariable_lh = (ithermo_water .NE. 0)
  ELSE  ! Default themodynamic is constant latent heat
     lvariable_lh = .false.
  END IF

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local and global variables
!------------------------------------------------------------------------------
  ! Input data
  !$ACC DATA &
  !$ACC   PRESENT(dz, t, p, rho, qv, qc, qi, qr, qs, qnc) &
  !$ACC   PRESENT(prr_gsp, prs_gsp, qrsflux, pri_gsp) &
  ! automatic arrays
  !$ACC   CREATE(zvzr, zvzs, zvzi) &
  !$ACC   CREATE(zpkr, zpks, zpki) &
  !$ACC   CREATE(zprvr, zprvs, zprvi, zqvsw_up) &
  !$ACC   CREATE(dist_cldtop, zlhv, zlhs)

! Some constant coefficients
  IF( lsuper_coolw) THEN
    znimax = znimax_Thom         !znimax_Thom = 250.E+3_wp,
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSEIF(lorig_icon) THEN
    znimax = 150.E+3_wp     ! from previous ICON code 
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSE
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    znimix = fxna(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  END IF

  zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0
  zlog_10 = LOG(10._wp) ! logarithm of 10

  ! Precomputations for optimization
  ccswxp_ln1o2   = EXP (ccswxp * LOG (0.5_wp))
  zvzxp_ln1o2    = EXP (zvzxp * LOG (0.5_wp))
  zbvi_ln1o2     = EXP (zbvi * LOG (0.5_wp))

! Optional arguments

  IF (PRESENT(ivstart)) THEN
    iv_start = ivstart
  ELSE
    iv_start = 1
  END IF
  IF (PRESENT(ivend)) THEN
    iv_end = ivend
  ELSE
    iv_end = nvec
  END IF
  IF (PRESENT(kstart)) THEN
    k_start = kstart
  ELSE
    k_start = 1
  END IF
  IF (PRESENT(idbg)) THEN
    izdebug = idbg
  ELSE
    izdebug = 0
  END IF
  IF (PRESENT(ldiag_ttend)) THEN
    lldiag_ttend = ldiag_ttend
  ELSE
    lldiag_ttend = .FALSE.
  ENDIF
  IF (PRESENT(ldiag_qtend)) THEN
    lldiag_qtend = ldiag_qtend
  ELSE
    lldiag_qtend = .FALSE.
  ENDIF

  !$ACC DATA CREATE(t_in) IF(lldiag_ttend)
  !$ACC DATA CREATE(qv_in, qc_in, qi_in, qr_in, qs_in) IF(lldiag_qtend)

  ! save input arrays for final tendency calculation
  IF (lldiag_ttend) THEN
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    t_in  = t
    !$ACC END KERNELS
  ENDIF
  IF (lldiag_qtend) THEN
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    qv_in = qv
    qc_in = qc
    qi_in = qi
    qr_in = qr
    qs_in = qs
    !$ACC END KERNELS
  END IF

! timestep for calculations
  zdtr  = 1.0_wp / zdt

! output for various debug levels
  IF (izdebug > 15) THEN
    CALL message('gscp_cloudice: ','Start of cloudice')
    WRITE (message_text,'(A,E10.3)') '      zams   = ',zams   ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
  END IF
  IF (izdebug > 20) THEN
    WRITE (message_text,*) '   nvec = ',nvec       ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke           ; CALL message('',message_text)
    WRITE (message_text,*) '   ivstart = ',ivstart ; CALL message('',message_text)
    WRITE (message_text,*) '   ivend   = ',ivend   ; CALL message('',message_text)
  END IF
  IF (izdebug > 50) THEN
#if defined( _OPENACC )
    CALL message('gscp_cloudice','GPU-info : update host before cloudice')
#endif
    !$ACC UPDATE HOST(dz, t, p, rho, qv, qc, qi, qr, qs) ASYNC(1)
    !$ACC WAIT(1)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN dz  = ',MAXVAL(dz),MINVAL(dz)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN T   = ',MAXVAL(t),MINVAL(t)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN p   = ',MAXVAL(p),MINVAL(p)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN rho = ',MAXVAL(rho),MINVAL(rho)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qv  = ',MAXVAL(qv),MINVAL(qv)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qc  = ',MAXVAL(qc),MINVAL(qc)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qr  = ',MAXVAL(qr),MINVAL(qr)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qi  = ',MAXVAL(qi),MINVAL(qi)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qs  = ',MAXVAL(qs),MINVAL(qs)
    CALL message('',message_text)
  ENDIF

  ! Delete precipitation fluxes from previous timestep
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR
  DO iv = iv_start, iv_end
    prr_gsp (iv) = 0.0_wp
    prs_gsp (iv) = 0.0_wp
    zpkr(iv)     = 0.0_wp
    zpks(iv)     = 0.0_wp
    zpki(iv)     = 0.0_wp
    zprvr(iv)    = 0.0_wp
    zprvs(iv)    = 0.0_wp
    zprvi(iv)    = 0.0_wp
    zvzr(iv)     = 0.0_wp
    zvzs(iv)     = 0.0_wp
    zvzi(iv)     = 0.0_wp
    dist_cldtop(iv) = 0.0_wp
    zqvsw_up(iv) = 0.0_wp
    pri_gsp (iv) = 0.0_wp
#ifndef __LOOP_EXCHANGE
    zlhv(iv)     = lh_v
    zlhs(iv)     = lh_s    
#endif
  END DO
  !$ACC END PARALLEL

  ! Initialize latent heats to constant values
#ifdef __LOOP_EXCHANGE
  zlhv(:) = lh_v
  zlhs(:) = lh_s
#endif

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! *********************************************************************

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP SEQ
#ifdef __LOOP_EXCHANGE
  DO iv = iv_start, iv_end  !loop over horizontal domain

    ! Calculate latent heats if necessary
    IF ( lvariable_lh ) THEN
      DO  k = k_start, ke  ! loop over levels
        tg      = make_normalized(t(iv,k))
        zlhv(k) = latent_heat_vaporization(tg)
        zlhs(k) = latent_heat_sublimation(tg)
      END DO
    END IF

    DO  k = k_start, ke  ! loop over levels
#else
  DO  k = k_start, ke  ! loop over levels

    ! Calculate latent heats if necessary
    IF ( lvariable_lh ) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tg)
      DO  iv = iv_start, iv_end  !loop over horizontal domain
        tg      = make_normalized(t(iv,k))
        zlhv(iv) = latent_heat_vaporization(tg)
        zlhs(iv) = latent_heat_sublimation(tg)
      END DO
    END IF

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(alf, bet, fnuc, llqc, llqi, llqr) &
    !$ACC   PRIVATE(llqs, m2s, m3s, maxevap, nnr, ppg, qcg) &
    !$ACC   PRIVATE(rhog, qig, qrg, qsg, qvg, reduce_dep) &
    !$ACC   PRIVATE(sagg, scac, scau, scfrz, sdau, sev, siau) &
    !$ACC   PRIVATE(sicri, sidep, simelt, snuc, srcri, srfrz) &
    !$ACC   PRIVATE(srim, ssdep, sshed, ssmelt, temp_c) &
    !$ACC   PRIVATE(tg, z1orhog, zbsdep, zcagg, zcidep, zcorr) &
    !$ACC   PRIVATE(zcrim, zcsdep, zcslam, zdtdh, zdvtp, zeff) &
    !$ACC   PRIVATE(zeln13o8qrk, zeln27o16qrk, zeln5o24qsk) &
    !$ACC   PRIVATE(zeln2o3qsk, zeln7o4qrk, zeln7o8qrk) &
    !$ACC   PRIVATE(zhi, zimi, zimr, zims, hlp) &
    !$ACC   PRIVATE(zlnlogmi, zlnqik, zlnqrk, zlnqsk) &
    !$ACC   PRIVATE(zmi, zn0s, znid, znin, zphi, zqct) &
    !$ACC   PRIVATE(zqik, zqit, zqrk, zqrt, zqsk, zqst) &
    !$ACC   PRIVATE(zqvsi, zqvsidiff, zqvsw, zqvsw0) &
    !$ACC   PRIVATE(zqvt, zrho1o2, zrhofac_qi, zscmax, zscsum) &
    !$ACC   PRIVATE(zsimax, zsisum, zsrmax, zsrsum) &
    !$ACC   PRIVATE(zssmax, zsssum, zsvidep, zsvisub, zsvmax) &
    !$ACC   PRIVATE(ztau, ztc, ztfrzdiff, ztt, zvz0s, zx1, zx2) &
    !$ACC   PRIVATE(zxfac, zzai, zzar, zzas, zztau)
    DO iv = iv_start, iv_end  !loop over horizontal domain
#endif

      ! add part of latent heating calculated in subroutine graupel to model latent
      ! heating field: subtract temperature from model latent heating field
      IF (ldass_lhn) THEN
        qrsflux(iv,k) = 0.0_wp
      ENDIF

      !----------------------------------------------------------------------------
      ! Section 2: Check for existence of rain and snow
      !            Initialize microphysics and sedimentation scheme
      !----------------------------------------------------------------------------

      zcrim  = 0.0_wp
      zcagg  = 0.0_wp
      zbsdep = 0.0_wp
      zvz0s  = 0.0_wp
      zn0s   = zn0s0
      reduce_dep = 1.0_wp  !FR: Reduction coeff. for dep. growth of rain and ice  

      !----------------------------------------------------------------------------
      ! 2.1: Preparations for computations and to check the different conditions
      !----------------------------------------------------------------------------

      qrg  = make_normalized(qr(iv,k))
      qsg  = make_normalized(qs(iv,k))
      qvg  = make_normalized(qv(iv,k))
      qcg  = make_normalized(qc(iv,k))
      qig  = make_normalized(qi(iv,k))
      tg   = t(iv,k)
      ppg  = p(iv,k)
      rhog = rho(iv,k)

      !..for density correction of fall speeds
      z1orhog = 1.0_wp/rhog
      hlp     = LOG(zrho0*z1orhog)
      zrho1o2 = EXP(hlp*x1o2)
      zrhofac_qi = EXP(hlp*icesedi_exp)

      zqrk = qrg * rhog
      zqsk = qsg * rhog
      zqik = qig * rhog

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqi = zqik > zqmin

      zdtdh = 0.5_wp * zdt / dz(iv,k)

      zzar = zqrk/zdtdh + zprvr(iv) + zpkr(iv)
      zzas = zqsk/zdtdh + zprvs(iv) + zpks(iv)
      zzai = zqik/zdtdh + zprvi(iv) + zpki(iv)

      zpkr(iv) = 0.0_wp
      zpks(iv) = 0.0_wp
      zpki(iv) = 0.0_wp

      !-------------------------------------------------------------------------
      ! qs_prepare:
      !-------------------------------------------------------------------------
      IF (llqs) THEN
        IF (isnow_n0temp == 1) THEN
          ! Calculate n0s using the temperature-dependent
          ! formula of Field et al. (2005)
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
          zn0s = zn0s1*EXP(zn0s2*ztc)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSEIF (isnow_n0temp == 2) THEN
          ! Calculate n0s using the temperature-dependent moment
          ! relations of Field et al. (2005) who assume bms=2.0
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)

          nnr  = 3._wp
          hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
              + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
              + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
          alf = EXP(hlp*zlog_10) ! 10.0_wp**hlp
          bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
              + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
              + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3

          ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
          m2s = qsg * rhog / zams   ! UB rho added as bugfix
          m3s = alf*EXP(bet*LOG(m2s))

          hlp  = zn0s1*EXP(zn0s2*ztc)
          zn0s = 13.50_wp * m2s * (m2s / m3s)**3
          zn0s = MAX(zn0s,0.5_wp*hlp)
          zn0s = MIN(zn0s,1.0E2_wp*hlp)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSE
          ! Old constant n0s
          zn0s = 8.0E5_wp
        ENDIF
        zcrim  = ccsrim*zn0s
        zcagg  = ccsagg*zn0s
        zbsdep = ccsdep*SQRT(v0snow)
        zvz0s  = ccsvel*EXP(ccsvxp * LOG(zn0s))
        zlnqsk = zvz0s * EXP (ccswxp * LOG (zqsk)) * zrho1o2
        ! Prevent terminal fall speed of snow from being zero at the surface level
        IF ( k == ke ) zlnqsk = MAX( zlnqsk, v_sedi_snow_min )
        zpks(iv) = zqsk * zlnqsk
        IF (zvzs(iv) == 0.0_wp) THEN
          zvzs(iv) = zlnqsk * ccswxp_ln1o2
        ENDIF
      ENDIF ! qs_prepare
    
      ! sedimentation fluxes

      !-------------------------------------------------------------------------
      ! qr_sedi:
      !-------------------------------------------------------------------------

      IF (llqr) THEN
        zlnqrk = zvz0r * EXP (zvzxp * LOG (zqrk)) * zrho1o2
        ! Prevent terminal fall speed of rain from being zero at the surface level
        IF ( k == ke ) zlnqrk = MAX( zlnqrk, v_sedi_rain_min )
        zpkr(iv) = zqrk * zlnqrk
        IF (zvzr(iv) == 0.0_wp) THEN
          zvzr(iv) = zlnqrk * zvzxp_ln1o2
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      ! qi_sedi:
      !-------------------------------------------------------------------------

      IF (llqi) THEN
        zlnqik = zvz0i * EXP (zbvi * LOG (zqik)) * zrhofac_qi
        zpki(iv) = zqik * zlnqik
        IF (zvzi(iv) == 0.0_wp) THEN
          zvzi(iv) = zlnqik * zbvi_ln1o2
        ENDIF
      ENDIF  ! qi_sedi

      ! Prevent terminal fall speeds of precip hydrometeors from being zero at the surface level
      IF ( k == ke ) THEN
        zvzr(iv) = MAX( zvzr(iv), v_sedi_rain_min )
        zvzs(iv) = MAX( zvzs(iv), v_sedi_snow_min )
      ENDIF

      !--------------------------------------------------------------------------
      ! 2.3: Second part of preparations
      !--------------------------------------------------------------------------

      zeln7o8qrk    = 0.0_wp
      zeln7o4qrk    = 0.0_wp !FR
      zeln27o16qrk  = 0.0_wp
      zeln13o8qrk   = 0.0_wp
      zeln5o24qsk   = 0.0_wp
      zeln2o3qsk    = 0.0_wp
      zsrmax        = 0.0_wp
      zssmax        = 0.0_wp

      !FR old
      !   zcsdep    = 3.2E-2_wp
      zcsdep        = 3.367E-2_wp
      zcidep        = 1.3E-5_wp
      zcslam        = 1e10_wp

      scau          = 0.0_wp
      scac          = 0.0_wp
      snuc          = 0.0_wp
      scfrz         = 0.0_wp
      simelt        = 0.0_wp
      sidep         = 0.0_wp
      ssdep         = 0.0_wp
      sdau          = 0.0_wp
      srim          = 0.0_wp
      sshed         = 0.0_wp
      sicri         = 0.0_wp
      srcri         = 0.0_wp
      sagg          = 0.0_wp
      siau          = 0.0_wp
      ssmelt        = 0.0_wp
      sev           = 0.0_wp
      srfrz         = 0.0_wp

      hlp = 1._wp/(rhog * r_v * tg)
      zqvsi = sat_pres_ice(tg) * hlp
      zqvsw = sat_pres_water(tg) * hlp

      zpkr(iv)   = MIN( zpkr(iv) , zzar )
      zpks(iv)   = MIN( zpks(iv) , zzas )
      zpki(iv)   = MIN( zpki(iv) , zzai )

      zzar   = zdtdh * (zzar-zpkr(iv))
      zzas   = zdtdh * (zzas-zpks(iv))
      zzai   = zdtdh * (zzai-zpki(iv))

      zimr   = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh)
      zims   = 1.0_wp / (1.0_wp + zvzs(iv) * zdtdh)
      zimi   = 1.0_wp / (1.0_wp + zvzi(iv) * zdtdh)

      zqrk   = zzar*zimr
      zqsk   = zzas*zims
      zqik   = zzai*zimi

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqi =  qig > zqmin 
      llqc =  qcg > zqmin

      !!----------------------------------------------------------------------------
      !! 2.4: IF (llqr): ic1
      !!----------------------------------------------------------------------------

      IF (llqr) THEN
        zlnqrk   = LOG (zqrk)
        IF ( qig+qcg > zqmin ) THEN
          zeln7o8qrk   = EXP (x7o8   * zlnqrk)
        ENDIF
        IF ( tg < ztrfrz ) THEN
          zeln7o4qrk   = EXP (x7o4   * zlnqrk) !FR new
          zeln27o16qrk = EXP (x27o16 * zlnqrk)
        ENDIF
        IF (llqi) THEN
          zeln13o8qrk  = EXP (x13o8  * zlnqrk)
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.5: IF (llqs): ic2
      !!----------------------------------------------------------------------------

      IF (llqs) THEN

        zlnqsk   = LOG (zqsk)
        zeln5o24qsk  = EXP (x5o24  * zlnqsk)
        zeln2o3qsk   = EXP (x2o3   * zlnqsk)

      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.6: -- This section is only used in the graupel scheme
      !!----------------------------------------------------------------------------

      !!----------------------------------------------------------------------------
      !! 2.7:  slope of snow PSD and coefficients for depositional growth (llqi,llqs; ic3)
      !!----------------------------------------------------------------------------    
llqi =  zqik > zqmin 
      IF (llqi .OR. llqs) THEN
        zdvtp  = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg
        zhi    = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
        hlp    = zdvtp / (1.0_wp + zhi)
        zcidep = ccidep * hlp

        IF (llqs) THEN
          zcslam = EXP(ccslxp * LOG(ccslam * zn0s / zqsk ))
          zcslam = MIN(zcslam,1.0E15_wp)
          zcsdep = 4.0_wp * zn0s * hlp
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.8: Deposition nucleation for low temperatures below a threshold (llqv)
      !!----------------------------------------------------------------------------    

      IF (( tg < zthet .AND. qvg >  8.E-6_wp &
                       .AND. qig <= 0.0_wp )) THEN
        IF( qvg > zqvsi ) THEN
          IF( lsuper_coolw .OR. lorig_icon) THEN
            znin  = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin  = MIN( fxna(tg), znimax )
          END IF
          snuc = zmi0 * z1orhog * znin * zdtr
        ENDIF
      ENDIF

      !!--------------------------------------------------------------------------
      !! Section 3: Search for cloudy grid points with cloud water and
      !!            calculation of the conversion rates involving qc (ic6)
      !!--------------------------------------------------------------------------
!X!
      IF (llqc) THEN

        zscmax = qcg*zdtr
        IF( tg > zthn ) THEN
          IF (iautocon == 0) THEN
            ! Kessler (1969) autoconversion rate
            scau = zccau * MAX( qcg - qc0, 0.0_wp )
            scac = zcac  * qcg * zeln7o8qrk
          ELSEIF (iautocon == 1) THEN
            ! Seifert and Beheng (2001) autoconversion rate
            ! with constant cloud droplet number concentration qnc
            IF (qcg > 1.0E-6_wp) THEN
              ztau = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
              ztau = MAX(ztau,1.E-30_wp)
              hlp  = EXP(zkphi2*LOG(ztau))
              zphi = zkphi1 * hlp * (1.0_wp - hlp)**3
              scau = zconst * qcg*qcg*qcg*qcg/(qnc(iv)*qnc(iv)) &
                   * (1.0_wp + zphi/(1.0_wp - ztau)**2)
              zphi = (ztau/(ztau+zkphi3))**4
              scac = zkcac * qcg * qrg * zphi
            ELSE
              scau = 0.0_wp
              scac = 0.0_wp
            ENDIF
          ENDIF
          IF (llqs) THEN
            srim = zcrim*  EXP(ccsaxp * LOG(zcslam)) * qcg
          ENDIF

          IF( tg >= t0 ) THEN
            sshed = srim
            srim  = 0.0_wp
          ENDIF
          ! Check for maximum depletion of cloud water and adjust the
          ! transfer rates accordingly
          zscsum = scau + scac + srim + sshed
          zcorr  = zscmax / MAX( zscmax, zscsum )
          scau   = zcorr*scau
          scac   = zcorr*scac
          srim   = zcorr*srim
          sshed  = zcorr*sshed

        ELSE !tg >= zthn: ! hom. freezing of cloud and rain water
          scfrz = zscmax
        ENDIF
        ! Calculation of heterogeneous nucleation of cloud ice.
        ! This is done in this section, because we require water saturation
        ! for this process (i.e. the existence of cloud water) to exist.
        ! Heterogeneous nucleation is assumed to occur only when no
        ! cloud ice is present and the temperature is below a nucleation
        ! threshold.
        IF( tg <= 267.15_wp .AND. qig <= 0.0_wp ) THEN   
          IF (lsuper_coolw .OR. lorig_icon) THEN
            znin  = MIN( fxna_cooper(tg), znimax )
            snuc = zmi0 * z1orhog * znin * zdtr
          ELSE
            znin = MIN( fxna(tg), znimax )
            snuc = zmi0 * z1orhog * znin * zdtr
          END IF
        ENDIF

        ! Calculation of in-cloud rainwater freezing
        IF ( tg < ztrfrz .AND. qrg > 0.1_wp*qcg ) THEN
          IF (lsuper_coolw) THEN
            srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
          ELSE
            ztfrzdiff = ztrfrz-tg
            srfrz = zcrfrz*ztfrzdiff*SQRT(ztfrzdiff)* zeln27o16qrk
          ENDIF
        ENDIF

        ! Calculation of reduction of depositional growth at cloud top (Forbes 2012)
        IF( k>k_start .AND. k<ke .AND. lred_depgrowth ) THEN
          znin = MIN(fxna_cooper(tg), znimax )
          fnuc = MIN(znin/znimix, 1.0_wp)

          !! distance from cloud top
          IF( qv(iv,k-1) + qc(iv,k-1) < zqvsw_up(iv) .AND. qi(iv,k-1) + qs(iv,k-1) < zqmin ) THEN ! upper cloud layer
            dist_cldtop(iv) = 0.0_wp    ! reset distance to upper cloud layer
          ELSE
            dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
          END IF

          ! with asymptotic behaviour dz -> 0 (xxx)
          !        reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
          !                             dist_cldtop(iv)/dist_cldtop_ref + &
          !                             (1.0_wp-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_wp)

          ! without asymptotic behaviour dz -> 0
          reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
                        dist_cldtop(iv)/dist_cldtop_ref), 1.0_wp)

        END IF ! Reduction of dep. growth of snow/ice 

      ENDIF

      !------------------------------------------------------------------------
      ! Section 4: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi and qs
      !------------------------------------------------------------------------

      IF ( llqi .OR. llqs ) THEN
        llqs =  zqsk > zqmin !zqsk > zqmin
        llqi =   qig > zqmin

        IF (tg<=t0) THEN           ! cold case 

          zqvsidiff = qvg-zqvsi
          zsvmax    = zqvsidiff * zdtr

          IF( lsuper_coolw .OR. lorig_icon) THEN
            znin   = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin   = MIN( fxna(tg), znimax )
          END IF
          ! Change in sticking efficiency needed in case of cloud ice sedimentation
          ! (based on Guenther Zaengls work)
          IF (lstickeff .OR. lorig_icon) THEN
            zeff     = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
            zeff     = MAX(zeff, zceff_min, zceff_fac*(tg-tmin_iceautoconv)) 
          ELSE !original sticking efficiency of cloud ice
            zeff     = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
            zeff     = MAX(zeff,0.2_wp)
          END IF
          zmi       = MIN( rhog*qig/znin, zmimax )
          zmi       = MAX( zmi0, zmi )
          sagg      = zcagg * EXP(ccsaxp*LOG(zcslam)) * qig * zeff
          siau      = zciau * MAX( qig - qi0, 0.0_wp ) * zeff
          znid      = rhog * qig/zmi
          IF (llqi) THEN
            zlnlogmi  = LOG (zmi)
            sidep     = zcidep * znid * EXP(0.33_wp * zlnlogmi) * zqvsidiff
          ELSE
            sidep = 0.0_wp
          ENDIF
          zsvidep   = 0.0_wp
          zsvisub   = 0.0_wp
          ! for sedimenting quantities the maximum 
          ! allowed depletion is determined by the predictor value. 
          IF (lsedi_ice .OR. lorig_icon) THEN
            zsimax  = zzai*z1orhog*zdtr
          ELSE
            zsimax  = qig*zdtr
          ENDIF
          IF( sidep > 0.0_wp ) THEN
            IF (lred_depgrowth ) THEN
              sidep = sidep * reduce_dep  !FR new: depositional growth reduction
            END IF
            zsvidep = MIN( sidep, zsvmax )
          ELSEIF ( sidep < 0.0_wp ) THEN
            IF (k < ke) THEN
              zsvisub  =   MAX (   sidep,  zsvmax)
              zsvisub  = - MAX ( zsvisub, -zsimax)
            ELSE
              zsvisub = - MAX(sidep, zsvmax )
              IF (zsvisub > zsimax) THEN
                zzai = zsvisub/(z1orhog*zdtr)
                zpki(iv) = MIN( zpki(iv) , zzai )
                zqik = zzai*zimi
                zsimax   = zsvisub
              ENDIF
            ENDIF
          ENDIF

          IF (llqi) THEN
            zlnlogmi   = LOG  (zmsmin/zmi)
            zztau      = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
            sdau       = zsvidep/MAX(zztau,zeps)
          ELSE
            sdau    =  0.0_wp
          ENDIF

          sicri      = zcicri * qig * zeln7o8qrk
          ! Allow growth of snow only if the existing amount of snow is sufficiently large 
          ! for a meaningful distiction between snow and cloud ice
          IF (qsg > 1.e-7_wp) srcri = zcrcri * (qig/zmi) * zeln13o8qrk

          zxfac = 1.0_wp + zbsdep * EXP(ccsdxp*LOG(zcslam))
          ssdep = zcsdep * zxfac * zqvsidiff / (zcslam+zeps)**2
          ! Check for maximal depletion of vapor by sdep
          IF (ssdep > 0.0_wp) THEN
            !FR new: depositional growth reduction
            IF (lred_depgrowth ) THEN
              ssdep = ssdep*reduce_dep
            END IF
            ! GZ: This limitation is crucial for numerical stability in the tropics!
            ssdep = MIN(ssdep, zsvmax-zsvidep)
          END IF


          ! Suppress depositional growth of snow if the existing amount is too small for a
          ! a meaningful distiction between cloud ice and snow
          IF (qsg <= 1.e-7_wp) ssdep = MIN(ssdep, 0.0_wp)

          ! Check for maximal depletion of cloud ice
          zsisum = siau + sdau + sagg + sicri + zsvisub
          zcorr  = 0.0_wp
          IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
          sidep  = zsvidep - zcorr*zsvisub
          sdau   = zcorr*sdau
          siau   = zcorr*siau
          sagg   = zcorr*sagg
          sicri  = zcorr*sicri

        ELSE ! tg > 0 - warm case

          !------------------------------------------------------------------------
          ! Section 5: Search for warm grid points with cloud ice and/or snow and
          !            calculation of the melting rates of qi and ps
          !------------------------------------------------------------------------

          ! cloud ice melts instantaneously
          IF (lsedi_ice .OR. lorig_icon) THEN
            simelt = zzai*z1orhog*zdtr
          ELSE
            simelt = qig*zdtr
          ENDIF

          zqvsw0     = zpvsw0/(rhog * r_v * t0)
          zx1        = (tg - t0) + zasmel*(qvg - zqvsw0)
          zx2        = 1.0_wp + zbsmel * zeln5o24qsk
          ssmelt     = zcsmel * zx1 * zx2 * zeln2o3qsk
          ssmelt     = MAX( ssmelt, 0.0_wp )

        ENDIF ! tg

      ENDIF

      !--------------------------------------------------------------------------
      ! Section 6: Search for grid points with rain in subsaturated areas
      !            and calculation of the evaporation rate of rain
      !--------------------------------------------------------------------------

      IF( (llqr) .AND. (qvg+qcg <= zqvsw)) THEN

        zlnqrk   = LOG (zqrk)
        zx1      = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
        ! Limit evaporation rate in order to avoid overshoots towards supersaturation
        ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
        temp_c = tg - t0
        maxevap     = (0.61_wp-0.0163_wp*temp_c+1.111e-4_wp*temp_c**2)*(zqvsw-qvg)/zdt
        sev    = MIN(zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk), maxevap)

        ! Calculation of below-cloud rainwater freezing
        IF ( tg < ztrfrz ) THEN
          IF (lsuper_coolw) THEN
            !FR new: reduced rain freezing rate
            srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
          ELSE
            srfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk
          ENDIF
        ENDIF


      ENDIF

      !--------------------------------------------------------------------------
      ! Section 7: Calculate the total tendencies of the prognostic variables.
      !            Update the prognostic variables in the interior domain.
      !--------------------------------------------------------------------------
      zsrmax = zzar*z1orhog*zdtr

      zssmax   = zzas * z1orhog * zdtr  

      zsrsum = sev + srfrz + srcri
      zcorr  = 1.0_wp
      IF(zsrsum > 0._wp) THEN
        zcorr  = zsrmax / MAX( zsrmax, zsrsum )
      ENDIF
      sev   = zcorr*sev
      srfrz = zcorr*srfrz
      srcri = zcorr*srcri
      ssmelt = MIN(ssmelt, zssmax)

      IF (ssdep < 0.0_wp ) THEN
        ssdep = MAX(ssdep, - zssmax)
      ENDIF      

      zqvt =   sev    - sidep  - ssdep  - snuc 
      zqct =   simelt - scau   - scfrz  - scac   - sshed  - srim 
      zqit =   snuc   + scfrz  - simelt - sicri  + sidep  - sdau   - sagg   - siau
      zqrt =   scau   + sshed  + scac   + ssmelt - sev    - srcri  - srfrz
      zqst =   siau   + sdau   + sagg   - ssmelt + sicri  + srcri  + srim   + ssdep + srfrz

#ifdef __LOOP_EXCHANGE      
      ztt = z_heat_cap_r*( zlhv(k)*(zqct+zqrt) + zlhs(k)*(zqit+zqst) )
#else
      ztt = z_heat_cap_r*( zlhv(iv)*(zqct+zqrt) + zlhs(iv)*(zqit+zqst) )
#endif

      ! Update variables and add qi to qrs for water loading
      IF (lsedi_ice .OR. lorig_icon) THEN
        qig = MAX ( 0.0_wp, (zzai*z1orhog + zqit*zdt)*zimi)
      ELSE
        qig = MAX ( 0.0_wp, qig + zqit*zdt)
      END IF
      qrg = MAX ( 0.0_wp, (zzar*z1orhog + zqrt*zdt)*zimr)
      qsg = MAX ( 0.0_wp, (zzas*z1orhog + zqst*zdt)*zims)

      !----------------------------------------------------------------------
      ! Section 8: Store satuaration specitic humidity at ice and water
      !            saturation of this layer 
      !----------------------------------------------------------------------
      zqvsw_up(iv) = zqvsw ! to be available in the layer below (layer k-1)
      
      !----------------------------------------------------------------------
      ! Section 10: Complete time step
      !----------------------------------------------------------------------

      IF ( k /= ke) THEN
        ! Store precipitation fluxes and sedimentation velocities 
        ! for the next level
        zprvr(iv) = qrg*rhog*zvzr(iv)
        zprvs(iv) = qsg*rhog*zvzs(iv)
        zprvi(iv) = qig*rhog*zvzi(iv)

        IF (zprvr(iv) <= zqmin) zprvr(iv)=0.0_wp
        IF (zprvs(iv) <= zqmin) zprvs(iv)=0.0_wp
        IF (zprvi(iv) <= zqmin) zprvi(iv)=0.0_wp

        ! for the latent heat nudging
        IF (ldass_lhn) THEN
          IF (lsedi_ice .OR. lorig_icon) THEN
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvi(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpki(iv))
          ELSE 
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv))
          END IF
        ENDIF

        IF (qrg+qr(iv,k+1) <= zqmin) THEN
          zvzr(iv)= 0.0_wp
        ELSE
          zvzr(iv)= zvz0r * EXP(zvzxp*LOG((qrg+qr(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qsg+qs(iv,k+1) <= zqmin) THEN
          zvzs(iv)= 0.0_wp
        ELSE
          zvzs(iv)= zvz0s * EXP(ccswxp*LOG((qsg+qs(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qig+qi(iv,k+1) <= zqmin ) THEN
          zvzi(iv)= 0.0_wp
        ELSE
          zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_wp*rhog)) * zrhofac_qi
        ENDIF
          
      ELSE
        ! Precipitation fluxes at the ground
        prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
        IF (lsedi_ice .OR. lorig_icon) THEN
          prs_gsp(iv) = 0.5_wp * (rhog*qsg*zvzs(iv) + zpks(iv))
          pri_gsp(iv) = 0.5_wp * (rhog*qig*zvzi(iv) + zpki(iv))
        ELSE
          prs_gsp(iv) = 0.5_wp * (qsg*rhog*zvzs(iv) + zpks(iv))
        END IF

        ! for the latent heat nudging
        IF (ldass_lhn) THEN
          qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)
        ENDIF


      ENDIF

      ! Update of prognostic variables or tendencies
      qr (iv,k) = qrg
      qs (iv,k) = qsg
      qi (iv,k) = qig
      t  (iv,k) = t (iv,k) + ztt*zdt 
      qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvt*zdt )
      qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqct*zdt )

#if !defined (_OPENACC) && !defined (__SX__)
      IF (izdebug > 15) THEN
        ! Check for negative values
        IF (qr(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qv'
          CALL message('',message_text)
        ENDIF
      ENDIF
#endif

    END DO  !loop over iv

  END DO ! loop over levels
  !$ACC END PARALLEL

!------------------------------------------------------------------------------
! final tendency calculation for ICON
!
! Note: as soon as we have a new satad subroutine in ICON, this tendency
! calculation will be done in the k-loop and the original 3D variables wont
! be used to store the new values. Then we wont need the _in variables anymore.
!------------------------------------------------------------------------------

! calculated pseudo-tendencies

  IF ( lldiag_ttend ) THEN
    !$ACC DATA &
    !$ACC   PRESENT(ddt_tend_t, t, t_in)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_t (iv,k) = (t (iv,k) - t_in (iv,k))*zdtr
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA
  ENDIF

  IF ( lldiag_qtend ) THEN
    !$ACC DATA &
    !$ACC   PRESENT(ddt_tend_qv, ddt_tend_qc, ddt_tend_qr, ddt_tend_qs) &
    !$ACC   PRESENT(ddt_tend_qi, qv_in, qc_in, qr_in, qs_in, qi_in)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_qv(iv,k) = MAX(-qv_in(iv,k)*zdtr,(qv(iv,k) - qv_in(iv,k))*zdtr)
        ddt_tend_qc(iv,k) = MAX(-qc_in(iv,k)*zdtr,(qc(iv,k) - qc_in(iv,k))*zdtr)
        ddt_tend_qi(iv,k) = MAX(-qi_in(iv,k)*zdtr,(qi(iv,k) - qi_in(iv,k))*zdtr)
        ddt_tend_qr(iv,k) = MAX(-qr_in(iv,k)*zdtr,(qr(iv,k) - qr_in(iv,k))*zdtr)
        ddt_tend_qs(iv,k) = MAX(-qs_in(iv,k)*zdtr,(qs(iv,k) - qs_in(iv,k))*zdtr)
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !$ACC END DATA
    
  ENDIF

  IF (izdebug > 15) THEN
#ifdef _OPENACC
   CALL message('gscp_cloudice', 'GPU-info : update host after cloudice')
#endif
   !$ACC UPDATE HOST(t, qv, qc, qi, qr, qs) ASYNC(1)
   !$ACC WAIT(1)
   CALL message('gscp_cloudice', 'UPDATED VARIABLES')
   WRITE(message_text,'(A,2E20.9)') 'cloudice  T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
   CALL message('', TRIM(message_text))
  ENDIF

  !$ACC END DATA ! IF(lldiag_qtend)
  !$ACC END DATA ! IF(lldiag_ttend)
  !$ACC END DATA ! general

!------------------------------------------------------------------------------
! End of subroutine cloudice
!------------------------------------------------------------------------------

END SUBROUTINE cloudice

!==============================================================================

END MODULE gscp_cloudice
