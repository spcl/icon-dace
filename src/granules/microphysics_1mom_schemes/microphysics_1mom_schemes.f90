! Public interfaces to the 1moment microhpysics-schemes
! This includes routines to
!   - initialize the microphysics
!   - run the microphysics
!   - getters and setters for variables private to microphysics
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
MODULE microphysics_1mom_schemes

  USE, INTRINSIC :: iso_fortran_env, ONLY: wp => real64
  USE gscp_graupel, ONLY: graupel
  USE gscp_cloudice, ONLY: cloudice
  USE gscp_kessler, ONLY: kessler
  USE gscp_ice, ONLY: cloudice2mom
  USE gscp_data, ONLY: & 
      zvz0i, znimax_Thom, zthn, isnow_n0temp, mu_rain, zami, &
      zams_ci, zams_gr, zbms, zcnue, &
      mma, mmb, gscp_set_coefficients, zmi0, zmimax, zn0r, ageo_snow,zmsmin, &
      zn0s1, zn0s2, zxiconv, cloud_num
  

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: microphysics_1mom_init
  PUBLIC :: graupel_run, cloudice_run, kessler_run, cloudice2mom_run
  PUBLIC :: set_terminal_fall_velocity_ice, get_terminal_fall_velocity_ice
  PUBLIC :: get_params_for_dbz_calculation
  PUBLIC :: get_params_for_reff_coefficients
  PUBLIC :: get_params_for_ncn_calculation
  PUBLIC :: get_mean_crystal_mass
  PUBLIC :: get_cloud_number
  PUBLIC :: get_snow_temperature
  

CONTAINS
  ! Returns the cloud number from module gscp_data
  SUBROUTINE get_cloud_number(cloud_number)
    REAL(wp), INTENT(OUT) :: cloud_number
    cloud_number = cloud_num
  END SUBROUTINE get_cloud_number

  ! Returns the mean crystal mass from module gscp_data
  SUBROUTINE get_mean_crystal_mass(mean_crystal_mass)
    REAL(wp), INTENT(OUT) :: mean_crystal_mass
    mean_crystal_mass = zxiconv
  END SUBROUTINE get_mean_crystal_mass

  ! Returns the parameters required for the reff coefficients from module gscp_data
  SUBROUTINE get_params_for_reff_coefficients(zami_arg, zmi0_arg, zmimax_arg, zn0r_arg, mu_rain_arg, &
    ageo_snow_arg, zbms_arg, zmsmin_arg)
    REAL(wp), INTENT(OUT) :: zami_arg, zmi0_arg, zmimax_arg, zn0r_arg, mu_rain_arg, &
                            ageo_snow_arg, zbms_arg, zmsmin_arg
    zami_arg = zami
    zmi0_arg = zmi0
    zmimax_arg = zmimax
    zn0r_arg = zn0r
    mu_rain_arg = mu_rain
    ageo_snow_arg = ageo_snow
    zbms_arg = zbms
    zmsmin_arg = zmsmin
  END SUBROUTINE get_params_for_reff_coefficients

  ! Returns the parameters required for the ncn calculation from module gscp_data
  SUBROUTINE get_params_for_ncn_calculation(isnow_n0temp_arg, ageo_snow_arg,&
                                          zn0s1_arg, zn0s2_arg, &
                                          znimax_Thom_arg, mma_arg, mmb_arg)
    REAL(wp), INTENT(OUT) :: ageo_snow_arg, zn0s1_arg, zn0s2_arg, &
                             znimax_Thom_arg, mma_arg(10),  &
                             mmb_arg(10)
    INTEGER, INTENT(OUT)  :: isnow_n0temp_arg
    isnow_n0temp_arg = isnow_n0temp
    znimax_Thom_arg = znimax_Thom
    zn0s1_arg = zn0s1
    zn0s2_arg = zn0s2
    mma_arg(:) = mma(:)
    mmb_arg(:) = mmb(:)
    ageo_snow_arg = ageo_snow
  END SUBROUTINE get_params_for_ncn_calculation

  ! Returns the parameters required for the dbz calculation from module gscp_data
  SUBROUTINE get_params_for_dbz_calculation(isnow_n0temp_arg, zami_arg, mu_rain_arg, zams_ci_arg, zams_gr_arg, zbms_arg, &
                                          znimax_Thom_arg, zthn_arg, mma_arg, mmb_arg, zcnue_arg)
    REAL(wp), INTENT(OUT) :: zami_arg, mu_rain_arg, zams_ci_arg, zams_gr_arg,  &
                             zbms_arg, znimax_Thom_arg, zthn_arg, mma_arg(10),  &
                             mmb_arg(10), zcnue_arg
    INTEGER, INTENT(OUT)  :: isnow_n0temp_arg
    isnow_n0temp_arg = isnow_n0temp
    zami_arg = zami
    mu_rain_arg = mu_rain
    zams_ci_arg = zams_ci
    zams_gr_arg = zams_gr
    zbms_arg = zbms
    znimax_Thom_arg = znimax_Thom
    zthn_arg = zthn
    mma_arg(:) = mma(:)
    mmb_arg(:) = mmb(:)
    zcnue_arg = zcnue
  END SUBROUTINE get_params_for_dbz_calculation

  ! Set the terminal fall velocity for ice from gscp_data
  SUBROUTINE set_terminal_fall_velocity_ice(new_terminal_fall_velocity)
    REAL(kind=wp), INTENT(IN) :: new_terminal_fall_velocity
    zvz0i = new_terminal_fall_velocity
  END SUBROUTINE

  ! Returns the terminal fall velocity for ice from gscp_data
  SUBROUTINE get_terminal_fall_velocity_ice(terminal_fall_velocity)
    REAL(kind=wp), INTENT(OUT) :: terminal_fall_velocity
    terminal_fall_velocity = zvz0i
  END SUBROUTINE

  ! Returns the snow temperature from gscp_data
  SUBROUTINE get_snow_temperature(snow_temperature)
    INTEGER, INTENT(OUT) :: snow_temperature
    snow_temperature = isnow_n0temp
  END SUBROUTINE

  ! Initializes the microphysics with the given coefficients stored in gscp_data
  SUBROUTINE microphysics_1mom_init( &
    igscp, &
    tune_zceff_min,  &
    tune_v0snow,  &
    tune_zvz0i, &
    tune_mu_rain,  &
    tune_icesedi_exp, &
    tune_rain_n0_factor)

    INTEGER  ,INTENT(IN) ::  igscp
    REAL(wp) ,INTENT(IN) ::  tune_zceff_min
    REAL(wp) ,INTENT(IN) ::  tune_v0snow
    REAL(wp) ,INTENT(IN) ::  tune_zvz0i
    REAL(wp) ,INTENT(IN) ::  tune_mu_rain
    REAL(wp) ,INTENT(IN) ::  tune_icesedi_exp
    REAL(wp) ,INTENT(IN) ::  tune_rain_n0_factor

    CALL gscp_set_coefficients( igscp    = igscp, &
      &                        tune_zceff_min   = tune_zceff_min, &
      &                        tune_v0snow      = tune_v0snow, &
      &                        tune_zvz0i       = tune_zvz0i, &
      &                        tune_icesedi_exp = tune_icesedi_exp, &
      &                        tune_mu_rain        = tune_mu_rain,&
      &                        tune_rain_n0_factor = tune_rain_n0_factor)

  END SUBROUTINE microphysics_1mom_init


  ! Calls the graupel microphysics scheme
  SUBROUTINE graupel_run(             &
    nvec,ke,                           & !> array dimensions
    ivstart,ivend, kstart,             & !! optional start/end indicies
    idbg,                              & !! optional debug level
    zdt, dz,                           & !! numerics parameters
    t,p,rho,qv,qc,qi,qr,qs,qg,qnc,     & !! prognostic variables
    qi0,qc0, zninc,                    & !! cloud ice/water threshold for autoconversion
    prr_gsp,prs_gsp,pri_gsp,prg_gsp,   & !! surface precipitation rates
    qrsflux,                           & !  total precipitation flux
    l_cv,                              &
    ithermo_water,                     & !  water thermodynamics
    ldass_lhn,                         &
    ldiag_ttend,     ldiag_qtend     , &
    ddt_tend_t     , ddt_tend_qv     , &
    ddt_tend_qc    , ddt_tend_qi     , & !> ddt_tend_xx are tendencies
    ddt_tend_qr    , ddt_tend_qs)!!    necessary for dynamics

    INTEGER, INTENT(IN) :: nvec          ,    & !> number of horizontal points
      ke                     !! number of grid points in vertical direction

    INTEGER, INTENT(IN) ::  ivstart   ,    & !> optional start index for horizontal direction
      ivend     ,    & !! optional end index   for horizontal direction
      kstart    ,    & !! optional start index for the vertical index
      idbg             !! optional debug level

    REAL(KIND=wp), INTENT(IN) :: zdt             ,    & !> time step for integration of microphysics     (  s  )
      qi0,qc0 !> cloud ice/water threshold for autoconversion

    REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: dz              ,    & !> layer thickness of full levels                (  m  )
      rho             ,    & !! density of moist air                          (kg/m3)
      p                      !! pressure                                      ( Pa  )

    LOGICAL, INTENT(IN):: l_cv, &                   !! if true, cv is used instead of cp
      ldass_lhn

    INTEGER, INTENT(IN):: ithermo_water          !! water thermodynamics

    LOGICAL, INTENT(IN):: ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
      ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

    REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::  t               ,    & !> temperature                                   (  K  )
      qv              ,    & !! specific water vapor content                  (kg/kg)
      qc              ,    & !! specific cloud water content                  (kg/kg)
      qi              ,    & !! specific cloud ice   content                  (kg/kg)
      qr              ,    & !! specific rain content                         (kg/kg)
      qs              ,    & !! specific snow content                         (kg/kg)
      qg              ,    & !! specific graupel content                      (kg/kg)
      zninc                  !! number of cloud ice crystals at nucleation

    REAL(KIND=wp), INTENT(INOUT) :: qrsflux(:,:)       ! total precipitation flux (nudg)

    REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::  prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
      prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
      prg_gsp,             & !! precipitation rate of graupel, grid-scale     (kg/(m2*s))
      qnc                    !! cloud number concentration

    REAL(KIND=wp), DIMENSION(:), INTENT(INOUT)::   pri_gsp                !! precipitation rate of ice, grid-scale        (kg/(m2*s))

    REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT)::   ddt_tend_t      , & !> tendency T                                       ( 1/s )
      ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
      ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
      ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
      ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
      ddt_tend_qs         !! tendency qs                                      ( 1/s )

    CALL graupel (                                     &
      & nvec   =nvec                            ,    & !> in:  actual array size
      & ke     =ke                              ,    & !< in:  actual array size
      & ivstart=ivstart                        ,    & !< in:  start index of calculation
      & ivend  =ivend                          ,    & !< in:  end index of calculation
      & kstart =kstart                  ,    & !< in:  vertical start index
      & zdt    =zdt                     ,    & !< in:  timestep
      & qi0    =qi0        ,    & 
      & qc0    =qc0        ,    & 
      & dz     =dz     ,    & !< in:  vertical layer thickness
      & t      =t           ,    & !< in:  temp,tracer,...
      & p      =p           ,    & !< in:  full level pres
      & rho    =rho          ,    & !< in:  density
      & qv     =qv    ,    & !< in:  spec. humidity
      & qc     =qc    ,    & !< in:  cloud water
      & qi     =qi    ,    & !< in:  cloud ice
      & qr     =qr    ,    & !< in:  rain water
      & qs     =qs    ,    & !< in:  snow
      & qg     =qg    ,    & !< in:  graupel
      & qnc    = qnc                            ,    & !< cloud number concentration
      & zninc  = zninc                          ,    & !< number of cloud ice crystals at nucleation
      & prr_gsp=prr_gsp     ,    & !< out: precipitation rate of rain
      & prs_gsp=prs_gsp     ,    & !< out: precipitation rate of snow
      & pri_gsp=pri_gsp      ,    & !< out: precipitation rate of cloud ice
      & prg_gsp=prg_gsp  ,    & !< out: precipitation rate of graupel
      & qrsflux= qrsflux       ,    & !< out: precipitation flux
      & ldiag_ttend = ldiag_ttend                 ,    & !< in:  if temp. tendency shall be diagnosed
      & ldiag_qtend = ldiag_qtend                 ,    & !< in:  if moisture tendencies shall be diagnosed
      & ddt_tend_t  = ddt_tend_t                  ,    & !< out: tendency temperature
      & ddt_tend_qv = ddt_tend_qv                 ,    & !< out: tendency QV
      & ddt_tend_qc = ddt_tend_qc                 ,    & !< out: tendency QC
      & ddt_tend_qi = ddt_tend_qi                 ,    & !< out: tendency QI
      & ddt_tend_qr = ddt_tend_qr                 ,    & !< out: tendency QR
      & ddt_tend_qs = ddt_tend_qs                 ,    & !< out: tendency QS
      & idbg=idbg                          ,    &
      & l_cv=l_cv                               ,    &
      & ldass_lhn = ldass_lhn                     ,    &
      & ithermo_water=ithermo_water)!< in: latent heat choice

  END SUBROUTINE graupel_run

  ! Calls the cloudice microphysics scheme
  SUBROUTINE cloudice_run(             &
    nvec,ke,                           & !> array dimensions
    ivstart,ivend, kstart,             & !! optional start/end indicies
    idbg,                              & !! optional debug level
    zdt, dz,                           & !! numerics parameters
    t,p,rho,qv,qc,qi,qr,qs,qnc,zninc,  & !! prognostic variables
    qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
    prr_gsp,prs_gsp,pri_gsp,   & !! surface precipitation rates
    qrsflux,                           & !  total precipitation flux
    l_cv,                              &
    ithermo_water,                     & !  water thermodynamics
    ldass_lhn,                         &
    ldiag_ttend,     ldiag_qtend     , &
    ddt_tend_t     , ddt_tend_qv     , &
    ddt_tend_qc    , ddt_tend_qi     , & !> ddt_tend_xx are tendencies
    ddt_tend_qr    , ddt_tend_qs)!!    necessary for dynamics

    INTEGER, INTENT(IN) :: nvec          ,    & !> number of horizontal points
      ke                     !! number of grid points in vertical direction

    INTEGER, INTENT(IN) ::  ivstart   ,    & !> optional start index for horizontal direction
      ivend     ,    & !! optional end index   for horizontal direction
      kstart    ,    & !! optional start index for the vertical index
      idbg             !! optional debug level

    REAL(KIND=wp), INTENT(IN) :: zdt             ,    & !> time step for integration of microphysics     (  s  )
      qi0,qc0 !> cloud ice/water threshold for autoconversion

    REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: dz              ,    & !> layer thickness of full levels                (  m  )
      rho             ,    & !! density of moist air                          (kg/m3)
      p                      !! pressure                                      ( Pa  )

    LOGICAL, INTENT(IN):: l_cv, &                   !! if true, cv is used instead of cp
      ldass_lhn

    INTEGER, INTENT(IN):: ithermo_water          !! water thermodynamics

    LOGICAL, INTENT(IN):: ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
      ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

    REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::  t               ,    & !> temperature                                   (  K  )
      qv              ,    & !! specific water vapor content                  (kg/kg)
      qc              ,    & !! specific cloud water content                  (kg/kg)
      qi              ,    & !! specific cloud ice   content                  (kg/kg)
      qr              ,    & !! specific rain content                         (kg/kg)
      qs              ,    & !! specific snow content                         (kg/kg)
      zninc                  !! number of cloud ice crystals at nucleation

    REAL(KIND=wp), INTENT(INOUT) :: qrsflux(:,:)       ! total precipitation flux (nudg)

    REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::  prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
      prs_gsp,             &  !! precipitation rate of snow, grid-scale        (kg/(m2*s))
      pri_gsp,             &
      qnc                     !! cloud number concentration


    REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT)::   ddt_tend_t      , & !> tendency T                                       ( 1/s )
      ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
      ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
      ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
      ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
      ddt_tend_qs         !! tendency qs                                      ( 1/s )

    CALL cloudice (                                     &
      & nvec   =nvec                            ,    & !> in:  actual array size
      & ke     =ke                              ,    & !< in:  actual array size
      & ivstart=ivstart                        ,    & !< in:  start index of calculation
      & ivend  =ivend                          ,    & !< in:  end index of calculation
      & kstart =kstart                  ,    & !< in:  vertical start index
      & zdt    =zdt                     ,    & !< in:  timestep
      & qi0    =qi0        ,    & 
      & qc0    =qc0        ,    & 
      & dz     =dz     ,    & !< in:  vertical layer thickness
      & t      =t           ,    & !< in:  temp,tracer,...
      & p      =p           ,    & !< in:  full level pres
      & rho    =rho          ,    & !< in:  density
      & qv     =qv    ,    & !< in:  spec. humidity
      & qc     =qc    ,    & !< in:  cloud water
      & qi     =qi    ,    & !< in:  cloud ice
      & qr     =qr    ,    & !< in:  rain water
      & qs     =qs    ,    & !< in:  snow
      & qnc    = qnc                            ,    & !< cloud number concentration
      & zninc  = zninc                          ,    & !< number of cloud ice crystals at nucleation
      & prr_gsp=prr_gsp     ,    & !< out: precipitation rate of rain
      & prs_gsp=prs_gsp     ,    & !< out: precipitation rate of snow
      & pri_gsp=pri_gsp      ,    & !< out: precipitation rate of cloud ice
      & qrsflux= qrsflux       ,    & !< out: precipitation flux
      & ldiag_ttend = ldiag_ttend                 ,    & !< in:  if temp. tendency shall be diagnosed
      & ldiag_qtend = ldiag_qtend                 ,    & !< in:  if moisture tendencies shall be diagnosed
      & ddt_tend_t  = ddt_tend_t                  ,    & !< out: tendency temperature
      & ddt_tend_qv = ddt_tend_qv                 ,    & !< out: tendency QV
      & ddt_tend_qc = ddt_tend_qc                 ,    & !< out: tendency QC
      & ddt_tend_qi = ddt_tend_qi                 ,    & !< out: tendency QI
      & ddt_tend_qr = ddt_tend_qr                 ,    & !< out: tendency QR
      & ddt_tend_qs = ddt_tend_qs                 ,    & !< out: tendency QS
      & idbg=idbg                          ,    &
      & l_cv=l_cv                               ,    &
      & ldass_lhn = ldass_lhn                     ,    &
      & ithermo_water=ithermo_water)!< in: latent heat choice

  END SUBROUTINE cloudice_run

  ! Calls the kessler microphysics scheme
  SUBROUTINE kessler_run(             &
    nvec,ke,                           & !> array dimensions
    ivstart,ivend, kstart,             & !! optional start/end indicies
    idbg,                              & !! optional debug level
    zdt, dz,                           & !! numerics parameters
    t,p,rho,qv,qc,qr,     & !! prognostic variables
    qc0,                           & !! cloud ice/water threshold for autoconversion
    prr_gsp,   & !! surface precipitation rates
    qrsflux,                           & !  total precipitation flux
    l_cv,                              &
    ldass_lhn,                         &
    ldiag_ttend,     ldiag_qtend     , &
    ddt_tend_t     , ddt_tend_qv     , &
    ddt_tend_qc, & !> ddt_tend_xx are tendencies
    ddt_tend_qr)!!    necessary for dynamics

    INTEGER, INTENT(IN) :: nvec          ,    & !> number of horizontal points
      ke                     !! number of grid points in vertical direction

    INTEGER, INTENT(IN) ::  ivstart   ,    & !> optional start index for horizontal direction
      ivend     ,    & !! optional end index   for horizontal direction
      kstart    ,    & !! optional start index for the vertical index
      idbg             !! optional debug level

    REAL(KIND=wp), INTENT(IN) :: zdt             ,    & !> time step for integration of microphysics     (  s  )
      qc0 !> cloud ice/water threshold for autoconversion

    REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: dz              ,    & !> layer thickness of full levels                (  m  )
      rho             ,    & !! density of moist air                          (kg/m3)
      p                      !! pressure                                      ( Pa  )

    LOGICAL, INTENT(IN):: l_cv, &                   !! if true, cv is used instead of cp
      ldass_lhn

    LOGICAL, INTENT(IN):: ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
      ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

    REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::  t               ,    & !> temperature                                   (  K  )
      qv              ,    & !! specific water vapor content                  (kg/kg)
      qr, &
      qc                     !! specific cloud water content                  (kg/kg)

    REAL(KIND=wp), INTENT(INOUT) :: qrsflux(:,:)       ! total precipitation flux (nudg)

    REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::  prr_gsp !> precipitation rate of rain, grid-scale        (kg/(m2*s))


    REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT)::   ddt_tend_t      , & !> tendency T                                       ( 1/s )
      ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
      ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
      ddt_tend_qr         !! tendency qr                                      ( 1/s )

    CALL kessler (                                     &
      & nvec   =nvec                            ,    & !> in:  actual array size
      & ke     =ke                              ,    & !< in:  actual array size
      & ivstart=ivstart                        ,    & !< in:  start index of calculation
      & ivend  =ivend                          ,    & !< in:  end index of calculation
      & kstart =kstart                  ,    & !< in:  vertical start index
      & zdt    =zdt                     ,    & !< in:  timestep
      & qc0    =qc0        ,    & 
      & dz     =dz     ,    & !< in:  vertical layer thickness
      & t      =t           ,    & !< in:  temp,tracer,...
      & p      =p           ,    & !< in:  full level pres
      & rho    =rho          ,    & !< in:  density
      & qv     =qv    ,    & !< in:  spec. humidity
      & qc     =qc    ,    & !< in:  cloud water
      & qr     =qr    ,    & 
      & prr_gsp=prr_gsp     ,    & !< out: precipitation rate of rain
      & qrsflux= qrsflux       ,    & !< out: precipitation flux
      & ldiag_ttend = ldiag_ttend                 ,    & !< in:  if temp. tendency shall be diagnosed
      & ldiag_qtend = ldiag_qtend                 ,    & !< in:  if moisture tendencies shall be diagnosed
      & ddt_tend_t  = ddt_tend_t                  ,    & !< out: tendency temperature
      & ddt_tend_qv = ddt_tend_qv                 ,    & !< out: tendency QV
      & ddt_tend_qc = ddt_tend_qc                 ,    & !< out: tendency QC
      & ddt_tend_qr = ddt_tend_qr                 ,    & !< out: tendency QR
      & idbg=idbg                          ,    &
      & l_cv=l_cv                               ,    &
      & ldass_lhn = ldass_lhn)

  END SUBROUTINE kessler_run

  ! Calls the cloudice2mom microphysics scheme
  SUBROUTINE cloudice2mom_run(             &
    nvec,ke,                           & !> array dimensions
    ivstart,ivend, kstart,             & !! optional start/end indicies
    idbg,                              & !! optional debug level
    zdt, dz,                           & !! numerics parameters
    w,t,p,rho,qv,qc,qi,qr,qs,qnc,     & !! prognostic variables
    qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
    qni, ninact, &
    prr_gsp,prs_gsp,pri_gsp,   & !! surface precipitation rates
    qrsflux,                           & !  total precipitation flux
    l_cv,                              &
    ithermo_water,                     & !  water thermodynamics
    ldass_lhn,                         &
    ldiag_ttend,     ldiag_qtend     , &
    ddt_tend_t     , ddt_tend_qv     , &
    ddt_tend_qc    , ddt_tend_qi     , & !> ddt_tend_xx are tendencies
    ddt_tend_qr    , ddt_tend_qs)!!    necessary for dynamics

    INTEGER, INTENT(IN) :: nvec          ,    & !> number of horizontal points
      ke                     !! number of grid points in vertical direction

    INTEGER, INTENT(IN) ::  ivstart   ,    & !> optional start index for horizontal direction
      ivend     ,    & !! optional end index   for horizontal direction
      kstart    ,    & !! optional start index for the vertical index
      idbg             !! optional debug level

    REAL(KIND=wp), INTENT(IN) :: zdt             ,    & !> time step for integration of microphysics     (  s  )
      qi0,qc0 !> cloud ice/water threshold for autoconversion

    REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: dz              ,    & !> layer thickness of full levels                (  m  )
      rho             ,    & !! density of moist air                          (kg/m3)
      p                      !! pressure                                      ( Pa  )

    LOGICAL, INTENT(IN):: l_cv, &                   !! if true, cv is used instead of cp
      ldass_lhn

    INTEGER, INTENT(IN):: ithermo_water          !! water thermodynamics

    LOGICAL, INTENT(IN):: ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
      ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

    REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::  t               ,    & !> temperature                                   (  K  )
      qv              ,    & !! specific water vapor content                  (kg/kg)
      qc              ,    & !! specific cloud water content                  (kg/kg)
      qi              ,    & !! specific cloud ice   content                  (kg/kg)
      qr              ,    & !! specific rain content                         (kg/kg)
      qs              ,    & !! specific snow content                         (kg/kg)
      qni             ,    & !< cloud ice number
      w               ,    &
      ninact                 !< activated ice nuclei

    REAL(KIND=wp), INTENT(INOUT) :: qrsflux(:,:)       ! total precipitation flux (nudg)

    REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::  prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
      prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
      pri_gsp, &
      qnc                    !! cloud number concentration


    REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT)::   ddt_tend_t      , & !> tendency T                                       ( 1/s )
      ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
      ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
      ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
      ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
      ddt_tend_qs         !! tendency qs                                      ( 1/s )

    CALL cloudice2mom(                                     &
      & nvec   =nvec                            ,    & !> in:  actual array size
      & ke     =ke                              ,    & !< in:  actual array size
      & ivstart=ivstart                        ,    & !< in:  start index of calculation
      & ivend  =ivend                          ,    & !< in:  end index of calculation
      & kstart =kstart                  ,    & !< in:  vertical start index
      & zdt    =zdt                     ,    & !< in:  timestep
      & qi0    =qi0        ,    & 
      & qc0    =qc0        ,    & 
      & dz     =dz     ,    & !< in:  vertical layer thickness
      & t      =t           ,    & !< in:  temp,tracer,...
      & p      =p           ,    & !< in:  full level pres
      & rho    =rho          ,    & !< in:  density
      & w = w, &
      & qv     =qv    ,    & !< in:  spec. humidity
      & qc     =qc    ,    & !< in:  cloud water
      & qi     =qi    ,    & !< in:  cloud ice
      & qr     =qr    ,    & !< in:  rain water
      & qs     =qs    ,    & !< in:  snow
      & qni    = qni  ,    & !< inout:  cloud ice number
      & ninact = ninact,   & !< inout:  activated ice nuclei
      & qnc    = qnc                            ,    & !< cloud number concentration
      & prr_gsp=prr_gsp     ,    & !< out: precipitation rate of rain
      & prs_gsp=prs_gsp     ,    & !< out: precipitation rate of snow
      & pri_gsp=pri_gsp      ,    & !< out: precipitation rate of cloud ice
      & qrsflux= qrsflux       ,    & !< out: precipitation flux
      & ldiag_ttend = ldiag_ttend                 ,    & !< in:  if temp. tendency shall be diagnosed
      & ldiag_qtend = ldiag_qtend                 ,    & !< in:  if moisture tendencies shall be diagnosed
      & ddt_tend_t  = ddt_tend_t                  ,    & !< out: tendency temperature
      & ddt_tend_qv = ddt_tend_qv                 ,    & !< out: tendency QV
      & ddt_tend_qc = ddt_tend_qc                 ,    & !< out: tendency QC
      & ddt_tend_qi = ddt_tend_qi                 ,    & !< out: tendency QI
      & ddt_tend_qr = ddt_tend_qr                 ,    & !< out: tendency QR
      & ddt_tend_qs = ddt_tend_qs                 ,    & !< out: tendency QS
      & idbg=idbg                          ,    &
      & l_cv=l_cv                               ,    &
      & ldass_lhn = ldass_lhn                     ,    &
      & ithermo_water=ithermo_water)!< in: latent heat choice

  END SUBROUTINE cloudice2mom_run

END MODULE microphysics_1mom_schemes

