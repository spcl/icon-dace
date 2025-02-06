!------------------------------------------------------------------------------
!    'cmn_fjx_mod.f90'  for fast-JX code v 7.4+ (prather 8/15)
!         note that module and file begin with 'cmn_"
!         small changes: LREF=51 instead of hardwired, JVMAP replaces JMAP
!         MASFAC param added, also cloud params - see below
!------------------------------------------------------------------------------
!
! NB - ALL of these common variables are set paramters,
!    They are NOT to be used as variables for a local solution
!    Thus this entire set is 'in' only after it is initialized
!-----------------------------------------------------------------------
!
! !INTERFACE:
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

      MODULE FJX_CMN_MOD
        USE mo_kind,                 ONLY: wp
        USE mo_physical_constants,   ONLY: avo, amd, grav, earth_radius
        USE mo_math_constants,       ONLY: pi, pi2, deg2rad

      implicit none
      public

!-----------------------------------------------------------------------
      ! Turn on/off RRTM-G absorption cross sections
      logical, parameter::  LRRTMG= .true.
!-----------------------------------------------------------------------

      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
     ! integer::  JXL_, JXL1_
      integer, parameter ::  JXL_=150, JXL1_=JXL_+1
      ! JXL2_: 2*JXL_ + 2 = mx no.levels in basic FJX grid (mid-level)
      !integer ::  JXL2_
      integer, parameter ::  JXL2_=2*JXL_+2
      ! W_   = dim = no. of Fast-J Wavelength bins:  currenly only 18, 
      ! TROP-ONLY is done by zeroing FL fluxes
      integer, parameter ::  W_=18
      ! S_   = dim = number of wavelength bins INLCUDING the Solar-J extensions (RRTMG value = 27)
      integer, parameter ::  S_=27
      ! X_   = dim = max no. of X-section data sets (input data)
      integer, parameter ::  X_=73
      ! A_   = dim = max no. of Aerosol Mie sets (input data) not including clouds and SSA
      integer, parameter ::  A_=40
      ! SSA_   = dim = no. of strat sulfate aerosol types (input data)
      integer, parameter ::  SSA_=18
      ! C_   = dim = no. of cld-data sets (input data) ///currently just 4xLiq and 2xIce
      integer, parameter ::  C_=6
      ! N_  = no. of levels in Mie scattering arrays
      !     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
      integer, parameter ::  N_=601
      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      integer, parameter ::  M_=4
      ! M2_ = 2*M_ = 8, replaces MFIT
      integer, parameter ::  M2_=2*M_

!-----------------------------------------------------------------------
      ! 4 Gauss pts = 8-stream
      real(wp), DIMENSION(M_), parameter  ::  &
                         EMU = [.06943184420297d0, .33000947820757d0, &
                                .66999052179243d0, .93056815579703d0]
      real(wp), DIMENSION(M_), parameter  :: &
                         WT  = [.17392742256873d0, .32607257743127d0, &
                                .32607257743127d0, .17392742256873d0]
!-----------------------------------------------------------------------

      ! MASFAC: Conversion factor for pressure to column density
      real(wp), parameter   ::  &
                         MASFAC = 100.d0*avo/(amd*grav*10.d0)
      ! ZZHT: scale height (cm) used above top of CTM ZHL(LPAR+1)
      real(wp), parameter   :: ZZHT = 5.d5
      ! RAD: Radius of Earth (cm)
      real(wp), parameter   :: RAD = earth_radius*100._wp
      ! ATAU: heating rate (factor increase from one layer to the next)
      real(wp), parameter   :: ATAU = 1.120d0
      ! ATAU0: minimum heating rate
      real(wp), parameter   :: ATAU0 = 0.010d0
      ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
      !ainteger :: JTAUMX 
      integer, parameter  :: JTAUMX = (N_ - 4*JXL_)/2
      character(LEN=25), dimension(8), parameter :: TITCLD =  &
      ['clear sky - no clouds    ', &
       'avg cloud cover          ', &
       'avg cloud cover^3/2      ', &
       'ICAs - avg direct beam   ', &
       'ICAs - random N ICAs     ', &
       'QCAs - midpt of bins     ', &
       'QCAs - avg clouds in bins', &
       'ICAs - use all ICAs***   ']

!---- Variables in file 'FJX_spec.dat' (RD_XXX)
      ! WL: Centres of wavelength bins - 'effective wavelength'  (nm)
      real(wp)  WL(S_)
      ! WBIN: Boundaries of wavelength bins                  (microns)
      real(wp)  WBIN(S_+1)
      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      ! FW: Solar flux in W/m2
      ! FP: PAR quantum action spectrum
      real(wp)  FL(S_),FW(S_),FP(S_)
      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      real(wp)  QRAYL(S_)
      ! SJSUB:  intended for breakdown of the super-bins (1:27) into smaller sub-bins.
      real(wp)  SJSUB(S_,15)
      real(wp)  QO2(W_,3)   ! QO2: O2 cross-sections
      real(wp)  QO3(W_,3)   ! QO3: O3 cross-sections
      real(wp)  Q1D(W_,3)   ! Q1D: O3 => O(1D) quantum yield

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      real(wp)  QQQ(W_,3,X_)
      ! TQQ: Temperature for supplied cross sections
      real(wp)  TQQ(3,X_)
      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      integer LQQ(X_)

      ! TITLEJX: Title (short & long) for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*6  TITLEJX(X_)
      CHARACTER*16 TITLEJL(X_)
      ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1  SQQ(X_)

!---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)
      ! TITLAA: Aerosol Mie Titles
      character*12  TITLAA(A_)
      ! QAA: Aerosol scattering phase functions
      real(wp)  QAA(5,A_)
      ! WAA: 5 Wavelengths for the supplied phase functions
      real(wp)  WAA(5,A_)
      ! PAA: Phase function: first 8 terms of expansion
      real(wp)  PAA(8,5,A_)
      ! RAA: Effective radius associated with aerosol type
      real(wp)  RAA(A_)
      ! SAA: Single scattering albedo
      real(wp)  SAA(5,A_)
      ! DAA: density (g/cm^3)
      real(wp)  DAA(A_)
      ! NAA: Number of categories for scattering phase functions
      integer NAA

!---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)
      ! NCC: Number of categories for cloud scattering phase functions
      integer NCC
      ! TITLCC: Cloud type titles
      character*12  TITLCC(C_)
      ! RCC: Effective radius associated with cloud type
      real(wp)  RCC(C_)
      ! GCC: Effective geometric cross section
      real(wp)  GCC(C_)
      ! DCC: density (g/cm^3)
      real(wp)  DCC(C_)
      ! QCC: Cloud Q-ext
      real(wp)  QCC(S_,C_)
      ! WCC: Wavelengths for supplied phase functions
      real(wp)  WCC(S_,C_)
      ! SCC: Single scattering albedo
      real(wp)  SCC(S_,C_)
      ! PCC: Phase function: first 8 terms of expansion
      real(wp)  PCC(8,S_,C_)

!---- Variables in file 'FJX_scat-ssa.dat' (RD_SSA)
      ! NSS: Number of categories for Strat Sulf Aerosol scattering phase functions
      integer NSS
      ! TITLSS: Cloud type titles
      character*12  TITLSS(SSA_)
      ! RSS: Effective radius associated with cloud type
      real(wp)  RSS(SSA_)
      ! GSS: Effective geometric cross section
      real(wp)  GSS(SSA_)
      ! DSS: density (g/cm^3)
      real(wp)  DSS(SSA_)
      ! TSS: temperature (K)
      real(wp)  TSS(SSA_)
      ! WSS: weight percent sulfuric acid (%)
      real(wp)  WSS(SSA_)
      ! QSS: Q-ext      ----begin wavelength dependent quantities
      real(wp)  QSS(S_,SSA_)
      ! SSS: Single scattering albedo
      real(wp)  SSS(S_,SSA_)
      ! PSS: Phase function: first 8 terms of expansion
      real(wp)  PSS(8,S_,SSA_)

!---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)
      ! WMM: U Michigan aerosol wavelengths
      real(wp)  WMM(6)
      ! UMAER: U Michigan aerosol data sets
      real(wp)  UMAER(3,6,21,33)

!---- Variables in file 'atmos_std.dat' (RD_PROF)
      integer, parameter ::  LREF=51   ! layer dim. in reference profiles
      integer, parameter ::  JREF=18   ! latitude dim. in reference profiles

!----- T and O3, H2O, CH4, reference profiles added underscore _ because TREF used in RRTMG_SW
      real(wp), DIMENSION(LREF,JREF,12) :: T_REF, O_REF, H2O_REF, CH4_REF
      integer NJX,NW1,NW2,NS1,NS2

!-----------------NEW for FJX72 parameters for cloud grid now here------
      integer, parameter :: &
            LPAR= 149, LWEPAR=34  &   !this can be set by CTM code
            !LPAR= 57, LWEPAR=34  &   !this can be set by CTM code
           ,L_=LPAR, L1_=L_+1 &   ! L_ = number of CTM layers
           ,L2_=2*L_+2 &        ! no. levels in the Fast-JX grid that
                       ! includes both layer edges and layer mid-points
           ,JVL_=LPAR &  ! vertical(levels) dim for J-values sent to CTM
           ,JVN_=101 &  ! max no. of J-values
           ,AN_=25     ! # FJX aerosols in layer (needs NDX for each)

!-----------------------------------------------------------------------
      ! variables used to map fast-JX J's onto CTM J's
!-----------------------------------------------------------------------
      real(wp)  JFACTA(JVN_)  ! multiplication factor for fast-JX calculated J
      integer JIND(JVN_)    ! index arrays that map Jvalue(j) onto rates
      integer NRATJ         ! number of Photolysis reactions in CTM chemistry, NRATJ <= JVN_
      integer NWBIN         ! used to shut off SFlux (and Fast-J calc) for strat wavelengths
                            ! only values that invoke action are 8 & 12 (mimic W_=8 or 12)
      integer NSBIN         ! likewise used to zero FL, if NSBIN = 18 then FL(19:27)=0.0

      character(LEN=6) JVMAP(JVN_) !label of J-value used to match w/FJX J's
      character(LEN=50) JLABEL(JVN_) ! label of J-value used in the chem model

! Cloud Cover parameters
!-----------------------------------------------------------------------
! NB CBIN_ was set at 20, but with NRG=6 groups, 10 gives a more reasonable number of ICAs
      integer, parameter :: CBIN_ = 10     ! # of quantized cloud fraction bins
      integer, parameter :: ICA_ = 20000   ! Max # of indep colm atmospheres
      integer, parameter :: NQD_ = 4       ! # of cloud fraction bins (4)

      real(wp),  parameter ::  CPI    = pi
      real(wp),  parameter ::  C2PI   = pi2
      real(wp),  parameter ::  CPI180 = deg2rad
      real(wp),  parameter ::  G0     = grav
      real(wp),  parameter ::  G100 = 100.d0/G0
!-------data to set up the random number sequence for use in cloud-JX
      integer, parameter :: NRAN_ = 10007  ! dimension for random number
      real(wp)   RAN4(NRAN_)      ! Random number set

! Solar J parameters - for Cloud-J v7.4 just do one sub-bin per RRTMg superbin
!-----------------------------------------------------------------------
      integer, parameter:: W_r = S_-W_ ! Cloud-J fix, W_r = 82 for RRTMg in rrsw_fasj_cmn.f90
      real(wp), dimension(L1_,W_r)  ::  TAUG_RRTMG    ! Cloud-J fix, is defined in rrsw_fasj_cmn.f90

      integer,  parameter, dimension(S_) ::   NGC =        &
        [  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,          &
           1,  1,  1,  1,  1,  1,  1,  1,  1,  1,          &
           1,  1,  1,  1,  1,  1,  1  ]
! actual RRMTMg values
!      integer,  parameter, dimension(S_) ::   NGC =       &
!       [   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,         &
!           1,  1,  1,  1,  1,  1,  1,  5, 10,  2,         &
!          10, 10,  8,  8, 12,  6, 12]

      END MODULE FJX_CMN_MOD
