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

!
! Data module for variables of the turbulence parameterization!
!
! Description of *turb_data*:
!  This module contains parameters that are used in the turbulence
!  parameterizations. With some of these parameters a tuning of the schemes
!  is possible.

MODULE turb_data

!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

USE mo_kind,                ONLY: wp           ! KIND-type parameter for real variables

USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

USE mo_turbdiff_config,     ONLY: turbdiff_config

!==============================================================================

IMPLICIT NONE

PUBLIC

!==============================================================================
! Fixed configuration parameters:
! ----------------------------------------------------------------------------
INTEGER, PARAMETER :: &
   ntmax=3,          & !maxmal number of time-levels for TKE
!
!  Indices for the two discriminated variable-types:
!
   mom=1,       & !momentum variables
   sca=mom+1,   & !scalar   variables
   ntyp=sca             !related number of variable-types ('mom' and 'sca')
!
!  Indices associated to paricular model variables:
!
INTEGER, PARAMETER :: &
   u_m=1,       & !zonal velocity-component at the mass center
   v_m=u_m+1,   & !meridional ,,      ,,    ,, ,,   ,,    ,,
   nvel=v_m,          & !number of velocity-components active for turbulece ('u_m', 'v_m')

   tem=nvel+1,  & !                       temperature
   tem_l=tem,   & !liquid-water           temperature
   tet=tem,     & !             potential temperature
   tet_l=tet,   & !liquid-water potential temperature 

   vap=tem+1,   & !water vapor (specific humidity 'qv')
   h2o_g=vap,   & !total water ('qv+qc')
   nred=h2o_g,        & !number of progn. variables being active for turbulence after reduction by local sat.adj.
   ninv=nred-nvel,    & !number of scalar variables being conserved during 'vap'<->'liq'-transistions ('tet_l', 'h2o_g')

   liq=vap+1,   & !liquid water (mass fraction 'qc')
   nmvar=liq,         & !number of progn. variables beding active for turbulence, which is equal to the
                        !number of variables being included into the single-column turbulence statistics
   nscal=liq-nvel,    & !number of scalar variables being active for turbulence ('tem', 'vap', 'liq')

   w_m=nmvar,   & !vertical velocity-component (only used for optional signle-column turbulence statistics)

!  Notice that: 'u_m, v_m'        are within [1, nvel];
!  but:         'tet_l, h2o_g'    are within [nvel+1, nred]
!  and:         'tem (tet), vap, liq' within [nvel+1, nmvar],
!  while:       'w_m'            is equal to 'nred+1=nmvar=liq'.

   naux=5,            & !(positive) number of auxilary variables

   ndim=MAX(nmvar,naux) !(positive) limit of last dimension used for 'zaux' and 'zvari'

!  Note:
!  The index "0" for the last dimension of 'zaux' and 'zvari' is also used, and it refers to 
!   'zaux' : saturation fraction (cl_cv) 
!   'zvari': potential available energy per volume (pressure) on half levels or
!            vertical gradient of effectively availabel kinetic energy (per mass)
!             due to near surface thermal inhomogeneity, being called here Circulation Kinetic Energy (CKE) 

LOGICAL, PARAMETER :: &
   ldynimp=.FALSE., &   !dynamical calculation of implicit weights for semi-implicit vertical diffusion
   lprecnd=.FALSE., &   !preconditioning of tridiagonal matrix      ,,     ,,          ,,        ,,

   lporous=.FALSE., &   !Vertically Resolved Roughness Layer (VRRL) representing a porous atmospheric medium
   !Note: The VRRL-treatment is not yet complete!

   ltst2ml =.FALSE., &   !test required, whether  2m-level is above the lowest main-level
   ltst10ml=.FALSE.      !test required, whether 10m-level is above the lowest main-level

   !Attention: 
   !So far, the  2m-level is assumed to be always below the lowest main-level,
   !    and the 10m-level is assumed to be always below the lowest half-level!!
   
!==============================================================================
! Parameters that may be used for tuning and special configurations:
! ----------------------------------------------------------------------------

! Attention:
! The given initializations are default settings of the boundary layer
! parameters. Some of these initial parameter values may be changed afterwards
! by model input NAMELISTs, respectively in SUB 'get_turbdiff_param'!

! 1. Numerical parameters:
!-----------------------------

REAL (KIND=wp)     ::        &
  impl_s       =  1.20_wp,   & ! implicit weight near the surface (maximal value)
  impl_t       =  0.75_wp,   & ! implicit weight near top of the atmosphere (maximal value)

  ! Minimal diffusion coefficients in [m^2/s] for vertical
  tkhmin       =  0.75_wp,   & ! scalar (heat) transport
  tkmmin       =  0.75_wp,   & ! momentum transport
  tkhmin_strat =  0.75_wp,   & ! scalar (heat) transport, enhanced value for stratosphere
  tkmmin_strat =  4.00_wp,   & ! momentum transport,      enhanced value for stratosphere

  ditsmot      =  0.00_wp,   & ! smoothing factor for direct time-step iteration

  tndsmot      =  0.00_wp,   & ! vertical smoothing factor for diffusion tendencies
  frcsmot      =  0.00_wp,   & ! vertical smoothing factor for TKE forcing (in ICON only in the tropics)
  tkesmot      =  0.15_wp,   & ! time smoothing factor for TKE and diffusion coefficients
  stbsmot      =  0.00_wp,   & ! time smoothing factor for stability function
  frcsecu      =  1.00_wp,   & ! security factor for TKE-forcing       (<=1)
  tkesecu      =  1.00_wp,   & ! security factor in  TKE equation      (out of [0; 1])
  stbsecu      =  0.00_wp,   & ! security factor in stability function (out of [0; 1])
  prfsecu      =  0.50_wp,   & ! relat. secur. fact. for prof. funct.  (out of ]0; 1[)

  epsi         =  1.0E-6_wp    ! relative limit of accuracy for comparison of numbers

INTEGER            ::        &
  it_end       =  1            ! number of initialization iterations (>=0)

! 2. Parameters describing physical properties of the lower boundary 
!    of the atmosphere:
!------------------------------------------

REAL (KIND=wp)     ::        &
  rlam_mom     =  0.0_wp,    & ! scaling factor of the laminar boundary layer for momentum
  rlam_heat    = 10.0_wp,    & ! scaling factor of the laminar boundary layer for heat

  rat_lam      =  1.0_wp,    & ! vapour/heat ratio of laminar scaling factors (over land)
  rat_sea      =  0.8_wp,    & ! sea/land ratio of laminar scaling factors for heat (and vapor)
  rat_glac     =  3.0_wp,    & ! glacier/land ratio of laminar scaling factors for heat (and vapor)

  rat_can      =  1.0_wp,    & ! ratio of canopy height over sai*z0m

  ! scaling factor for additional shear-forcing by Non-Turbulent subgrid Circulations (NTCs)
  ! or via Lower Limits of Diffusion-Coefficients (LLDCs) in the surface layer:
  rsur_sher    =  0.0_wp,    & ! (so far deactivated)

  z0m_dia      =  0.2_wp,    & ! roughness length of a typical synoptic station [m]

  alpha0       =  0.0123_wp, & ! Charnock-parameter
  alpha0_max   =  0.0335_wp, & ! upper limit of velocity-dependent Charnock-parameter
  alpha0_pert  =  0.0_wp,    & ! additive ensemble perturbation of Charnock-parameter

  alpha1       =  0.7500_wp    ! parameter scaling the molecular roughness of water waves

  !$ACC DECLARE COPYIN(alpha0, alpha0_max, alpha0_pert)

! 3. Parameters that should be external parameter fields being not yet 
!    available:
!------------------------------------------

REAL (KIND=wp)     ::        &
  c_lnd        = 2.0_wp,     & ! surface area density of the roughness elements over land
  c_sea        = 1.5_wp,     & ! surface area density of the waves over sea
  c_soil       = 1.0_wp,     & ! surface area density of the (evaporative) soil surface
  c_stm        = 0.0_wp,     & ! (so far deactivated)
! c_stm        = 2.5_wp,     & ! surface area density of stems and branches at the plant-covered part of the surface
  !Note:
  !"c_stm=0" matches with the previous not consistant scaling, considering "c_stm=2.5" only as an additional surface 
  ! part that can hold interception water.
  e_surf       = 1.0_wp        ! exponent to get the effective surface area


! 4. Parameters that should be dynamical fields being not yet available:
!------------------------------------------

REAL (KIND=wp)     ::        &
  z0_ice       =  0.001_wp     ! roughness length of sea ice


! 5. Parameters for modelling turbulent diffusion:
!------------------------------------------

REAL (KIND=wp)     ::        &
  tur_len      = 500.0_wp,   & ! asymptotic maximal turbulent distance [m]
  pat_len      = 100.0_wp,   & ! effective global length scale of subscale surface patterns over land [m]
                               ! (should be dependent on location)
  len_min      =  1.0E-6_wp, & ! minimal turbulent length scale [m]

  vel_min      =  0.01_wp,   & ! minimal velocity scale [m/s]

  akt          =  0.4_wp,    & ! von Karman-constant

  ! Length scale factors for pressure destruction of turbulent
  a_heat       =  0.74_wp,   & ! scalar (heat) transport
  a_mom        =  0.92_wp,   & ! momentum transport

  ! Length scale factors for dissipation of turbulent
  d_heat       =  10.1_wp,   & ! scalar (temperature) variance
  d_mom        =  16.6_wp,   & ! momentum variance

  ! Length scale factor for turbulent transport (vertical diffusion) of TKE
  c_diff       =  0.20_wp,   & ! (including turb. pressure-transport)

  ! Length scale factor for separate horizontal shear circulations 
  a_hshr       =  1.00_wp,   & ! contributing to shear-production of TKE

  ! Length scale factor for the stability correction
  a_stab       =  0.00_wp,   & ! applied to integral turbulent length-scale

  ! Dimensionless parameters used in the sub grid scale condensation scheme
  ! (statistical cloud scheme):

  clc_diag     =  0.5_wp,    & ! cloud cover at saturation
  q_crit       =  1.6_wp,    & ! critical value for normalized super-saturation

  c_scld       =  1.0_wp       ! shape-factor (0<=c_scld) applied to pure 'cl_cv' at the moist correct. 
                               !  by turbulent phase-transit., providing an eff. 'cl_cv', which scales 
                               !  the implicit liquid-water flux under turbulent sat.-adj.:
                               !  <1: small eff. 'cl_cv' even at large pure 'cl_cv'
                               !  =1:       eff. 'cl_cv' just equals   pure 'cl_cv'
                               !  >1: large eff. 'cl_cv' even at small pure 'cl_cv'

!==============================================================================
! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

LOGICAL :: &
  loldtur       =.FALSE., & ! use settings to simulate old ijk turbulence version
                            ! if .TRUE.: new ICON-like settings are used
  ltkecon       =.FALSE., & ! consider convective buoyancy production in TKE-equation
  ltkeshs       =.TRUE. , & ! consider separ. horiz. shear production in TKE-equation
  loutshs       =.TRUE. , & ! consider separ. horiz. shear production of TKE for output
  ltkesso       =.TRUE. , & ! consider mechanical SSO-wake production in TKE-equation
  loutsso       =.TRUE. , & ! consider mechanical SSO-wake production of TKE for output
  ltkenst       =.TRUE. , & ! consider produc. by near-surf. thermals in TKE-equation
  loutnst       =.FALSE., & ! consider produc. by near-surf. thermals of TKE for output
!    Output of additional TKE production terms:
  loutbms       =.FALSE., & ! consider TKE-production by turbulent buoyancy, total mechanical shear
                            !  or grid-scale mechanical shear for additional output

  lnonloc       =.FALSE., & ! nonlocal calculation of vertical gradients used for turbul. diff.
  lprfcor       =.FALSE., & ! using the profile values of the lowest main level instead of
                            ! the mean value of the lowest layer for surface flux calulations

  ltmpcor       =.FALSE., & ! consideration minor turbulent sources in the enthalpy budget 
  lcpfluc       =.FALSE.    ! consideration of fluctuations of the heat capacity of air

LOGICAL :: &
  lexpcor       =.FALSE., & ! explicit warm-cloud correct. of implicitly calculated turbul. diff.
  lsflcnd       =.TRUE. , & ! lower flux condition for vertical diffusion calculation
  lcirflx       =.FALSE., & ! consideration of non-turbulent fluxes related to near-surface circulations
  lfreeslip     =.FALSE.    ! free-slip lower boundary condition (use for idealized runs only!)

! Notice that the following switches are provided by the parameter-list of 
! SUB 'turbdiff' or 'turbtran':

! lstfnct                   :calculation of stability function required
! lnsfdia                   :calculation of (synoptical) near-surface variables required
!! lmomdif (obsolete)       :calculation of complete gradient diffusion of horizontal momenum
! lum_dif                   :calculation of vertical gradient diffusion of horizontal u-momenum
! lum_dif                   :calculation of vertical gradient diffusion of horizontal v-momenum
! lscadif                   :calculation of complete gradient diffusion of scalar properties
! lturatm                   :running turbulence model between atmosph. layers (updating diffusion coefficients)
!!ltursrf (obsolete)        :running turbulence model at the surface layer (updating transfer coefficients
! lsfluse                   :use explicit heat flux densities at the suface
! ltkeinp                   :TKE present as input (at level k=ke1 for current time level 'ntur')
! lgz0inp                   :gz0 present as input

!==============================================================================
! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

INTEGER :: &
  imode_tran    =0, & ! mode of TKE-equation in transfer scheme             (compare 'imode_turb')
  imode_turb    =1, & ! mode of TKE-equation in turbulence scheme
                      !  0: diagnostic equation
                      !  1: prognostic equation (default)
                      !  2: prognostic equation (implicitly positive definit)
  icldm_tran    =2, & ! mode of cloud representation in transfer parametr.  (compare 'icldm_turb')
  icldm_turb    =2, & ! mode of cloud representation in turbulence parametr.
                      ! -1: ignoring cloud water completely (pure dry scheme)
                      !  0: no clouds considered (all cloud water is evaporated)
                      !  1: only grid scale condensation possible
                      !  2: also sub grid (turbulent) condensation considered
  itype_wcld    =2, & ! type of water cloud diagnosis within the turbulence scheme:
                      ! 1: employing a scheme based on relative humitidy
                      ! 2: employing a statistical saturation adjustment
  itype_sher    =0    ! type of mean shear-production for TKE
                      ! 0: only vertical shear of horizontal wind
                      ! 1: previous plus horizontal shear correction
                      ! 2: previous plus shear from vertical velocity

! To reproduce the old ijk turbulence settings as good as possible, all these switches 
! have to be set to 1. If loldtur=.TRUE., the re-setting is done in organize_physics.

! These are the settings for the ICON-like setup of the physics
INTEGER :: &
! imode_stbcorr =2, & ! mode of correcting the stability function (related to 'stbsecu')
  imode_stbcorr =1, & ! mode of correcting the stability function (related to 'stbsecu')
                      ! 1: always for strict.-non-stb. strat. using a restr. 'gama' in terms of prev. forc.
                      ! 2: only to avoid non-physic. solution or if current 'gama' is too large
  ilow_def_cond =2, & ! type of the default condition at the lower boundary
                      ! 1: zero surface gradient 
                      ! 2: zero surface value
  imode_calcirc =2, & ! mode of treating the raw "circulation term" (related to 'pat_len', imode_pat_len')
                      ! 1: explicit calculation of the flux convergence
                      ! 2: quasi implicit treatment by calculation of effective TKE-gradients

  !Note: 
  !The theoretical background of the "circulation term" has meanwhile been fundamentally revised.
  !Accordingly, it is going to be substituted by two complementary approaches:
  !i) a thermal SSO parameterization and ii) a new "circulation term" due to thermal surface patterns.
  !Hence in the following, the still active raw parameterization is referred to as raw "circulation term".

  imode_pat_len =2, & ! mode of determining a length scale of surface patterns used for the "circulation-term" 
                      !  and additional roughness by tile-variation of land-use:
                      ! 1: employing the constant value 'pat_len' only
                      !    - raw "circulation term" considered as to be due to thermal surface-patterns.
                      ! 2: using the standard deviat. of SGS orography as a lower limit 
                      !    - raw "circulation term" considered as to be due to thermal SSO effect,
  imode_frcsmot =2, & ! if "frcsmot>0", apply smoothing of TKE source terms 
                      ! 1: globally or 
                      ! 2: in the tropics only (if 'trop_mask' is present) 
  imode_shshear =2, & ! mode of calculat. the separated horizontal shear mode (related to 'ltkeshs', 'a_hshr')
                      ! 0: with a constant lenght scale and based on 3D-shear and incompressibility
                      ! 1: with a constant lenght scale and considering the trace constraint for the 2D-strain tensor
                      ! 2: with a Ri-number depend. length sclale correct. and the trace constraint for the 2D-strain tensor
  imode_tkesso  =1, & ! mode of calculat. the SSO source term for TKE production
                      ! 1: original implementation
                      ! 2: with a Ri-dependent reduction factor for Ri>1
                      ! 3: as "2", but additional reduction for mesh sizes < 2 km
  imode_tkvmini =2, & ! mode of calculating the minimal turbulent diff. coeffecients
                      ! 1: with a constant value
                      ! 2: with a stability dependent correction
  imode_rat_sea =1, & ! mode of scaling the laminar resistance for heat over sea (related to 'rat_sea')
                      ! 1: constant ratio 'rat_sea' to land 
                      ! 2: with a correction for a strongly overheated SST
  imode_vel_min =2, & ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                      ! 1: with a constant value
                      ! 2: with a stability dependent correction of 'vel_min'
  imode_snowsmot=1, & ! mode to treating the aerodynamic surface-smoothing by snow
                      ! 0: no smoothing active at all
                      ! 1: no impact on SAI, but full smoothing of z0 (G. Zaengl's approach)
                      ! 2: "1", but with full smoothing of SAI: full smoothing of z0 and SAI
                      ! 3: dynamical smoothing of z0 and SAI dependent on snow- and roughness height
  imode_charpar =2    ! mode of estimating the Charnock-Parameter
                      ! 1: use a constant value 
                      ! 2: use a wind-dependent value with a constant upper bound
                      ! 3: as "2", but with reduction at wind speeds above 25 m/s for more realistic TC wind speeds

  !$ACC DECLARE COPYIN(imode_charpar)

INTEGER :: &

  itype_2m_diag =1, & ! type of 2m-diagnostics for temperature and -dewpoint
                      ! 1: Considering a fictive surface roughness of a SYNOP lawn
                      ! 2: Considering the mean surface roughness of a grid box
                      !    and using an exponential roughness layer profile
  imode_2m_diag =2, & ! mode of 2m-diagnostics of temperature and dew-point (related to 'itype_2m_diag')
                      ! (-)1: direct interpolation of temperature and specific humidity
                      ! (-)2: interpol. of conserved quantities and subsequent statistical saturation adjustm.,
                      !       allowing particularly for the diagnostic of cloud water at the 2m-level (fog)
                      !  > 0: extra pressure-calculat. at 2m-level
                      !  < 0: surface pressure applied at 2m-level
  imode_nsf_wind=1, & ! mode of local wind-definition at near-surface levels (applied for 10m wind
                      !  diagnostics as well as for calculation of sea-surface roughness)
                      ! 1: ordinary wind speed (magnitude of grid-scale averaged wind-vector)
                      ! 2: including relative wind-speed amplification by NTCs  (at "rsur_shear>0")
                      !                                        or even by LLDCs (at "imode_suradap=3")
  imode_qvsatur =2, & ! mode of calculating the saturat. humidity
                      ! 1: old version using total pressure
                      ! 2: new version using partial pressure of dry air
  imode_stadlim =2, & ! mode of limitting statist. saturation adjustment (SUB 'turb_cloud')
                      ! 1: only absolut upper limit of stand. dev. of local super-satur. (sdsd)
                      ! 2: relative limit of sdsd and upper limit of cloud-water
  imode_trancnf =2, & ! mode of configuring the transfer-scheme (SUB 'turbtran')
                      ! 1: old version: start. with lamin. diffus.; with a lamin. correct. for profile-funct.;
                      !    interpol. T_s rather then Tet_l onto zero-level; calcul. only approx. Tet_l-grads.;
                      !    using an upper bound for TKE-forcing; without transmit. skin-layer depth to turbul.
                      ! 2: 1-st ConSAT: start. with estim. Ustar, without a laminar correct. for prof.-funct.;
                      !    interpol. Tet_l onto zero-level; calcul. Tet_l-gradients directly; 
                      !    without an upper bound for TKE-forcing; with transmit. skin-layer depth to turbul.
                      ! 3: 2-nd ConSAT: as "2", but with a hyperbolic interpol. of profile function
                      !    for stable stratification
                      ! 4: 3-rd ConSAT: as "3", but without using an upper interpolation node
  imode_lamdiff =1, & ! mode of considering laminar diffusion at surface layer
                      ! 1: not applied for the profile functions in case of "imode_trancnf.GE.2"
                      ! 2: always applied
  imode_tkemini =1, & ! mode of adapting q=2TKE**2 and the TMod. to Lower Limits for Diff. Coeffs. (LLDCs)
                      ! 1: LLDC treated as corrections of stability length without any further adaptation
                      ! 2: TKE adapted to that part of LLDC representing so far missing shear forcing, while the
                      !     assumed part of LLDC representing missing drag-forces has no feedback to the TMod.
                      !Notice:
                      !Above the surface layer, the TKE-adaption of "imode_tkemini=2" is always applied
                      ! in case of "rsur_sher>0"!
  imode_suradap =0, & ! mode of adapting surface-layer profile-functions to Lower Limits of Diff. Coeffs. (LLDCs)
                      ! 0: no adaptations at all
                      ! 1: removing the artific. drag contrib. by the LLDC for momentum at level "k=ke"
                      ! 2: "1" and also removing  shear contrib. by LLDCs at level "k=ke" 
                      ! 3: "1" and employing shear contrib. by LLDCs at surface-level "k=ke1"
                      !Notice:
                      !Any shear contrib. by NTCs or LLDCs is only considered at surf.-lev., if "rsur_sher.GT.0".
  imode_tkediff =2, & ! mode of implicit TKE-Diffusion (related to 'c_diff')
                      ! 1: in terms of q=SQRT(2*TKE)) 
                      ! 2: in terms of TKE=0.5*q**2
  imode_adshear =2    ! mode of considering addit. shear by scale interaction (realt. to 'ltkesso', 'ltkeshs',
                      ! 'ltkecon', 'ltkenst')
                      ! 1: not  considered for stability functions
                      ! 2: also considered for stability functions

! Notice that the following selectors are provided by the parameter-list of 
! SUB 'turb_diffusion' or 'turb_transfer':

! iini                :type of initialization (0: no, 1: separate before the time loop
!                                                   , 2: within the first time step)
! itnd                :type of tendency cons. (0: no, 1: in implicit vertical diffusion equation
!                                                     2: by adding to current profile before vertical diffusion
!                                                     3: by using corrected virtual vertical profiles

!==============================================================================
! Declarations of utility variables:

! Turbulence parameters which are computed during model run
!-------------------------------------------------------------------------------

REAL (KIND=wp)     ::        &
  ! do we need it as TARGET?
  ! these variables are set in SR turb_param
  c_tke,tet_g,rim, &
  c_m,c_h, b_m,b_h,  sm_0, sh_0, &
  d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
  a_3,a_5,a_6,                 &

  ! these parameters are physical constants which are either taken as
  ! they are or set to 0.0 in the turbulence for special applications
  tur_rcpv,         & ! cp_v/cp_d - 1
  tur_rcpl            ! cp_l/cp_d - 1 (where cp_l=cv_l)

! Definition of used data types
!-------------------------------

TYPE modvar !model variable
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
             av(:,:) => NULL(), & !atmospheric values
             sv(:)   => NULL(), & !surface     values (concentration or flux density)
             at(:,:) => NULL()    !atmospheric time tendencies
     LOGICAL                                 ::         &
             fc                   !surface values are flux densities
     INTEGER                                 ::         &
             kstart  = 1          !start level for vertical diffusion
END TYPE modvar

TYPE turvar !turbulence variables
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
             tkv(:,:) => NULL(), & !turbulent coefficient for vert. diff.
             tsv(:)   => NULL()    !turbulent velocity at the surface
END TYPE turvar

TYPE varprf !variable profile
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
             bl(:,:), & !variable at boundary model levels
             ml(:,:)    !variable at main     model levels
END TYPE varprf


! 7. Switches from the COSMO-Model: must only be defined for ICON
! ---------------------------------------------------------------

REAL(KIND=wp), POINTER :: &
  impl_weight(:)  ! implicit weights for tridiagonal solver

! Switches controlling turbulent diffusion:
! ------------------------------------------


!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE get_turbdiff_param (jg)

   ! This subroutines overwrites the above given initial values of those parameters, switches or selectors,
   !  which belong to the ICON NAMELIST 'turbdiff_nml', with the given NAMELIST settings.
   ! It is called in the SUBs 'nwp_turbtrans' and 'nwp_turbdiff' before calling the parameterization schemes
   !  'turbtran' and 'turbdiff' respectively.

   !Attention: 
   !Without this setting, perturbations applied to any of these parameters are not effective 
   ! in 'turbdiff' or 'turbtran'!

   INTEGER, INTENT(IN) :: jg !patch index

   impl_weight   => turbdiff_config(jg)%impl_weight

   imode_tran   = turbdiff_config(jg)%imode_tran
   icldm_tran   = turbdiff_config(jg)%icldm_tran
   imode_turb   = turbdiff_config(jg)%imode_turb
   icldm_turb   = turbdiff_config(jg)%icldm_turb
   itype_wcld   = turbdiff_config(jg)%itype_wcld
   itype_sher   = turbdiff_config(jg)%itype_sher
   imode_shshear= turbdiff_config(jg)%imode_shshear
   imode_frcsmot= turbdiff_config(jg)%imode_frcsmot
   imode_tkesso = turbdiff_config(jg)%imode_tkesso

   ltkesso      = turbdiff_config(jg)%ltkesso
   ltkecon      = turbdiff_config(jg)%ltkecon
   ltkeshs      = turbdiff_config(jg)%ltkeshs
   lexpcor      = turbdiff_config(jg)%lexpcor
   ltmpcor      = turbdiff_config(jg)%ltmpcor
   lprfcor      = turbdiff_config(jg)%lprfcor
   lnonloc      = turbdiff_config(jg)%lnonloc
   lfreeslip    = turbdiff_config(jg)%lfreeslip
   lcpfluc      = turbdiff_config(jg)%lcpfluc
   lsflcnd      = turbdiff_config(jg)%lsflcnd

   tur_len      = turbdiff_config(jg)%tur_len
   pat_len      = turbdiff_config(jg)%pat_len
   a_stab       = turbdiff_config(jg)%a_stab
   a_hshr       = turbdiff_config(jg)%a_hshr
   !Attention(MR): Without this setting, perturbations applied to 'a_hshr' are not effective in 'turbdiff'!

   impl_s         = turbdiff_config(jg)%impl_s
   impl_t         = turbdiff_config(jg)%impl_t
   c_diff         = turbdiff_config(jg)%c_diff
   tkhmin         = turbdiff_config(jg)%tkhmin
   tkmmin         = turbdiff_config(jg)%tkmmin
   tkhmin_strat   = turbdiff_config(jg)%tkhmin_strat
   tkmmin_strat   = turbdiff_config(jg)%tkmmin_strat

   tkesmot        = turbdiff_config(jg)%tkesmot
   frcsmot        = turbdiff_config(jg)%frcsmot

   imode_snowsmot = turbdiff_config(jg)%imode_snowsmot
   imode_charpar  = turbdiff_config(jg)%imode_charpar
   alpha0         = turbdiff_config(jg)%alpha0
   alpha0_max     = turbdiff_config(jg)%alpha0_max
   alpha0_pert    = turbdiff_config(jg)%alpha0_pert
   alpha1         = turbdiff_config(jg)%alpha1

   !$ACC UPDATE DEVICE(imode_charpar, alpha0, alpha0_max, alpha0_pert) ASYNC(1)

   rlam_heat      = turbdiff_config(jg)%rlam_heat
   rlam_mom       = turbdiff_config(jg)%rlam_mom
   rat_lam        = turbdiff_config(jg)%rat_lam
   rat_sea        = turbdiff_config(jg)%rat_sea
   rat_glac       = turbdiff_config(jg)%rat_glac

   q_crit         = turbdiff_config(jg)%q_crit

END SUBROUTINE get_turbdiff_param
!==============================================================================

END MODULE turb_data
