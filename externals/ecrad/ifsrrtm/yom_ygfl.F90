! # 1 "ifsrrtm/yom_ygfl.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yom_ygfl.f90"
module yom_ygfl

use parkind1 , only : jpim, jprb

implicit none
save

!-------------------------------------------------------------------------
! contains the descriptors of gfl arrays
!-------------------------------------------------------------------------

! jpgfl : max number of gfl fields
! jpnamed_gfl : number of currently pre-defined components of gfl
! jpghg : number of greenhouse gas fields
! jpgrg : number of reactive gas fields
! jpchem : number of chemical species
! jpaero : number of active aerosol fields
! jpaerout: number of output aerosol fields
! jpuvp : number of output from uv processor
! jptrac : number of tracers for diagnostics
! jpera40 : number of era40 diagnostic fields
! jpch4s  : number of added fields related to methane
! jpnogw  : number of diagnostic fields for noro gwd scheme
! jpsldia : number of sl dynamics diagnostic fields
! jpchem_assim : maximum number of assimilated of chemical species
!-------------------------------------------------------------------------

integer(kind=jpim), parameter :: jpgfl=2163
integer(kind=jpim), parameter :: jpnamed_gfl=27
integer(kind=jpim), parameter :: jpghg=3
integer(kind=jpim), parameter :: jptrac=10
integer(kind=jpim), parameter :: jpchem=66
integer(kind=jpim), parameter :: jpgrg=5       
integer(kind=jpim), parameter :: jpchem_assim=5
integer(kind=jpim), parameter :: jpaero=16
integer(kind=jpim), parameter :: jpforc=1800
integer(kind=jpim), parameter :: jpera40=14
integer(kind=jpim), parameter :: jpsldia=7
integer(kind=jpim), parameter :: jpezdiag=50
integer(kind=jpim), parameter :: jpch4s=2
integer(kind=jpim), parameter :: jpnogw=2
integer(kind=jpim), parameter :: jpaerout=17
integer(kind=jpim), parameter :: jpuvp=2
integer(kind=jpim), parameter :: jpphys=8   
integer(kind=jpim), parameter :: grib_code_gfl_phys=81  ! ajgdb hopefully harmless

type type_gfl_comp ! individual field descriptor

sequence ! daand: necessary to avoid memory corruption with gfortran 4.3.3

character(len=16)  :: cname     = ''        ! arpege field name 
integer(kind=jpim) :: igrbcode  = -999      ! grib code
logical            :: ladv      = .false.   ! field advected or not
logical            :: ladv5     = .false.   ! field advected without wind increments
logical            :: ltdiablin = .false.   ! diabatic tendency is interpolated by lin. int.
logical            :: lhorturb  = .false.   ! horizontal part affected by 3d turbulence
integer(kind=jpim) :: nreqin    = 0         ! 1 if field requiered in input, 0 if not, -1 if initialised
                                            ! with a reference value refvali
logical            :: lreqout   = .false.   ! t if field requiered in output
logical            :: lgpingp   = .true.    ! gp field input as gp
logical            :: lgp       = .false.   ! field exists and of grid-point type
logical            :: lsp       = .false.   ! field exists and of spectral type
logical            :: lcders    = .false.   ! derivatives required (spectral only)
logical            :: lactive   = .false.   ! field in use
logical            :: lthermact = .false.   ! field thermodynamically active
real(kind=jprb)    :: r         = 0.0_jprb
real(kind=jprb)    :: rcp       = 0.0_jprb
logical            :: lt9       = .false.   ! field in t-dt gfl
logical            :: lt1       = .false.   ! field in t+dt gfl
logical            :: lt5       = .false.   ! field in trajectory gfl
logical            :: lphy      = .false.   ! field in physics gfl
logical            :: lpt       = .false.   ! field in pc phy. tend. gfl (gflpt)
logical            :: ltrajio   = .false.   ! field written to and from trajectory structure
logical            :: ldiag     = .false.   ! field is "diagnostic" at t; e.g. cloud fraction 
logical            :: lpc       = .false.   ! field in predictor/corrector time stepping (gflpc)
real(kind=jprb)    :: refvali   = 0.0_jprb  ! reference value for init, used in case nreqin==-1
! lam specific attributes (arome/aladin)
logical            :: ladjust0  = .false.   ! true if field is thermodynamically adjusted at t
                                            ! (immediatly after inverse spectral transforms)
logical            :: ladjust1  = .false.   ! true if field is thermodynamically adjusted at t+dt
                                            ! (after sl interpolations and nl residuals)
integer(kind=jpim) :: ncoupling = 0         ! 1 if field is coupled by davies relaxation, 0 if not,
                                            ! -1 if coupled with reference value for coupling refvalc
real(kind=jprb)    :: refvalc   = 0.0_jprb  ! reference value for coupling, used in case ncoupling==-1
logical            :: lbiper    = .false.   ! true if field must be biperiodised inside the transforms
! end lam specific attributes (arome/aladin)
character(len=12)  :: cslint    = ''        ! s.l interpolaion "type"
integer(kind=jpim) :: mp        = -99999999 ! basic field "pointer"
integer(kind=jpim) :: mpl       = -99999999 ! zonal derivative "pointer"
integer(kind=jpim) :: mpm       = -99999999 ! meridional derivative "pointer"
integer(kind=jpim) :: mp9       = -99999999 ! basic field "pointer" t-dt
integer(kind=jpim) :: mp9_ph    = -99999999 ! basic field "pointer" for physics
integer(kind=jpim) :: mp1       = -99999999 ! basic field "pointer" t+dt
integer(kind=jpim) :: mp5       = -99999999 ! basic field "pointer" trajectory
integer(kind=jpim) :: mp5l      = -99999999 ! zonal derivative "pointer" trajectory
integer(kind=jpim) :: mp5m      = -99999999 ! meridional derivative "pointer" trajectory
integer(kind=jpim) :: mpslp     = -99999999 ! basic field "pointer" physics
integer(kind=jpim) :: mpsp      = -99999999 ! basic field "pointer" spectral space
integer(kind=jpim) :: mp_spl    = -99999999 ! basic field "pointer" spline interpolation
integer(kind=jpim) :: mp_sl1    = -99999999 ! basic field "pointer" in slbuf1
integer(kind=jpim) :: mp_slx    = -99999999 ! basic field "pointer" in slbuf1 for cpg_pt
integer(kind=jpim) :: mppt      = -99999999 ! physics tendency "pointer"
integer(kind=jpim) :: mppc      = -99999999 ! predictor/corrector auxiliary array "pointer"

! daand: intflex attributes
logical            :: lwater                ! true for water species
logical            :: lprecip               ! true for precipitating water species
real(kind=jprb)    :: rlzer                 ! latent heat change at 0k

! gems nl ext
integer(kind=jpim) :: ncouplo4              ! coupled to ctm by oasis4 intefrace
logical            :: lassim                ! use as control variable (either monitored or assimilated)
integer(kind=jpim) :: igribdv               ! grib code of deposition velocity 
integer(kind=jpim) :: igribtc               ! grib code of total column
integer(kind=jpim) :: igribsfc              ! grib code of surface flux 
logical            :: ldiff                 ! diffusion  on
logical            :: lconv                 ! convection on
real(kind=jprb)    :: rmolmass              ! molar mass 
real(kind=jprb)    :: refold                ! efolding decay time 
real(kind=jprb)    :: henrya                ! henry constant a 
real(kind=jprb)    :: henryb                ! henry constant b 
logical            :: lnegfix               ! cut off negative values in sugridug an
logical            :: lmassfix              ! correct mass error of sl advection in gpmodel (if ltrcmfix)
type(type_gfl_comp),pointer :: previous     ! pointer to previously def. field

end type type_gfl_comp

type type_gfl_naml ! individual field descriptor for namelist input

sequence ! daand: necessary to avoid memory corruption with gfortran 4.3.3

character(len=16)  :: cname     ! arpege field name 
integer(kind=jpim) :: igrbcode  ! grib code
integer(kind=jpim) :: nreqin    ! 1 if field required in input, 0 if not, -1 if initialised
                                ! with a reference value refvali
real(kind=jprb) :: refvali      ! reference value for initialisation, used in case nreqin==-1
logical :: lreqout              ! t if field requiered in output
logical :: lgpingp              ! gp field input as gp
logical :: lgp                  ! field exists and of grid-point type
logical :: lsp                  ! field exists and of spectral type
logical :: lcders               ! derivatives required (spectral only)
logical :: lt9                  ! field in t-dt gfl
logical :: lt1                  ! field in t+dt gfl
logical :: lt5                  ! field in trajectory gfl
logical :: lphy                 ! field with physics tendencies gfl
logical :: lpt                  ! field in pc physics tendency gflpt
logical :: ltrajio              ! field written to and from trajectory structure
logical :: ldiag                ! field is "diagnostic" at t; e.g. cloud fraction 
logical :: lpc                  ! field in predictor/corrector time stepping gflpc
logical :: ladv                 ! field advected or not
logical :: ladv5                ! field advected without wind increments
logical :: lintlin              ! linear interpolation for field
logical :: ltdiablin            ! diabatic tendency is interpolated by linear int.
logical :: lhorturb             ! horizontal part affected by 3d turbulence
logical :: lqm                  ! quasi-monotonous interpolation for field
logical :: lqmh                 ! quasi-monotonous interpolation in horizontal for field
logical :: lqm3d                ! quasi-monotone interpolation applied directly in 3 dimensions
logical :: lslhd                ! semi-lagrangian horizontal diffusion used for field
logical :: lcomad               ! comad weights used for sl interpolation of field
logical :: lhv                  ! hermite vertical interpolation used for field (only ozone sofar)
logical :: lvsplip              ! vertical spline interpolation used for field (only ozone sofar)
integer(kind=jpim) :: ncoupling ! 1 if field is coupled by davies relaxation, 0 if not,
                                ! -1 if coupled with reference value for coupling refvalc
real(kind=jprb) :: refvalc      ! reference value for coupling, used in case 
                                ! ncoupling==-1
! gems nl ext
integer(kind=jpim)  :: ncouplo4 ! coupled to ctm by oasis4 intefrace =1 input,=2 in&output,=-1 none
logical             :: lassim   ! use as control variable (either monitored or assimilated)
integer(kind=jpim)  :: igribdv  ! grib code of deposition velocity 
integer(kind=jpim)  :: igribtc  ! grib code of total column
integer(kind=jpim)  :: igribsfc ! grib code of surface flux 
logical             :: ldiff    ! diffusion  on
logical             :: lconv    ! convection on
logical             :: lnegfix  ! cut off negative values in sugridug and callpar
logical             :: lmassfix ! correct mass error of sl advection in gpmodel (if ltrcmfix)
real(kind=jprb)     :: rmolmass ! molar mass 
real(kind=jprb)     :: refold   ! efolding  decay time 
real(kind=jprb)     :: henrya   ! henry constant a 
real(kind=jprb)     :: henryb   ! henry constant b 

end type type_gfl_naml

!-------------------------------------------------------------------------
! derived types for describing the gfl structure.
!-------------------------------------------------------------------------
! modifications:
! 03/07/09 c. fischer - add arome/aladin attributes
! 03/10/01 c. moussy  - add arome/aladin attributes coupling
! 03/10/31 m. tudor   - add physics tendencies for predictor-corrector
! 05/10/10 j. haseler - switch for i/o to trajectory structure
! 2004-nov f. vana    - update of cslint attribute
! 20-feb-2005 vivoda  - 3tl eul pc scheme (gflpc)
! 07/06/27 e. holm    - tl/ad advection without wind increments ladv5
! 12/04/08 j. flemming - gfl attribute extention for gems 
! 22-feb-11 f. vana   - ltdiablin and lhorturb
! spring 2011 ecmwf   - lintlin
! nov. 2013           - lcomad
! 2013-11, d. degrauwe - intflex attributes

type type_gfld

sequence ! daand: necessary to avoid memory corruption with gfortran 4.3.3

! overall descriptor,dimensioning etc.
integer(kind=jpim) :: numflds     = 0  ! number of gfl fields
integer(kind=jpim) :: nders       = 0  ! number of horizontal derivatives fields
integer(kind=jpim) :: numspflds   = 0  ! number of spectrally represented gfl fields
integer(kind=jpim) :: numgpflds   = 0  ! number of grid-point gfl fields
integer(kind=jpim) :: numflds9    = 0  ! number of gfl fields in (t-dt) part
integer(kind=jpim) :: numflds1    = 0  ! number of gfl fields in (t+dt) array
integer(kind=jpim) :: numspflds1  = 0  ! number of spectrally represented gfl fields (t+dt)
integer(kind=jpim) :: numflds5    = 0  ! number of gfl fields (trajectory)
integer(kind=jpim) :: numfldsphy  = 0  ! number of gfl fields (phys.)
integer(kind=jpim) :: numflds_spl = 0  ! number of gfl fields (s.l. spline interpolation)
integer(kind=jpim) :: numflds_sl1 = 0  ! number of gfl fields in s.l. buffer 1
integer(kind=jpim) :: numfldspc   = 0  ! number of gfl fields (predictor/corrector)
integer(kind=jpim) :: ndim        = 0  ! dimension of main array holding gfl fields(gfl)
integer(kind=jpim) :: numfldspt   = 0  ! number of gfl fields (phy. tend.)
integer(kind=jpim) :: ndim0       = 0  ! dimension of t0 part of gfl
integer(kind=jpim) :: ndim9       = 0  ! dimension of t-dt part of gfl
integer(kind=jpim) :: ndim1       = 0  ! dimension of t+dt array (gflt1)
integer(kind=jpim) :: ndim5       = 0  ! dimension of traj. gfl array (gfl5)
integer(kind=jpim) :: ndimslp     = 0  ! diminsion of s.l. phys. gfl array (gflslp)
integer(kind=jpim) :: ndim_spl    = 0  ! dim. of arrays holding gfl fields (s.l.spline int.)
integer(kind=jpim) :: ndimpt      = 0  ! dimension of phy. tend. gfl array (gflpt)
integer(kind=jpim) :: ndimpc      = 0  ! dimension of iterative scheme auxiliary array (gflpc)

integer(kind=jpim) :: ngfl_ext
integer(kind=jpim) :: ngfl_forc
integer(kind=jpim) :: ngfl_ezdiag
integer(kind=jpim) :: nghg
integer(kind=jpim) :: ntrac
integer(kind=jpim) :: ngrg
integer(kind=jpim) :: ngrg_cplo4
integer(kind=jpim) :: ngrg_assim
integer(kind=jpim) :: naero
integer(kind=jpim) :: nactaero
integer(kind=jpim) :: nddhaero
integer(kind=jpim) :: nera40
integer(kind=jpim) :: nnogw
integer(kind=jpim) :: naerout
integer(kind=jpim) :: nuvp
integer(kind=jpim) :: nsldia
integer(kind=jpim) :: nsldiagp
integer(kind=jpim) :: ngfl_phys
logical :: lco2sfc
logical :: lch4sfc
logical :: laerosfc
logical :: lfire
logical :: laerodiu
logical :: ltrcmfix       ! activates tracer mass fixer
logical :: ltrcmfix_ps    ! adjust pressure to conserve dry mass in mass fixer calculations
logical :: laerout
logical :: luvpout
logical :: lchem

integer(kind=jpim) :: ngems   ! the total number of "gems" fields.
integer(kind=jpim) :: nchem
integer(kind=jpim) :: nchem_assim
integer(kind=jpim) :: nchem_flx 
integer(kind=jpim) :: nchem_dv
integer(kind=jpim) :: nchem_tc
integer(kind=jpim) :: nchem_scv

!     ------------------------------------------------------------------
!      mass fixers
!     ------------------------------------------------------------------
integer(kind=jpim) :: nnegafix     ! num of fields to apply -ve fixer
integer(kind=jpim) :: noptnegfix   ! 1: simple negative fixer (reset to 0)
                                   ! 2: reset to local minimum

logical :: lqm3dcons      ! bermejo & staniforth quasi-monotone limiter with improved
                          ! conservation option. when true, applied to all gfl s.t. lqm3d=true
logical :: ladvnegfix              ! activates negative fixer for advection
logical :: ltrcmfbc                ! activate bermejo & conde if true
logical :: ltrcmfpr                ! activate priestley algorithm if true
logical :: ltrcmfmg                ! activate mac gregor's algorithm if true
logical :: lextradf                ! extra diagnostics 


integer(kind=jpim) :: nfldsfix     ! number of fields to be fixed
integer(kind=jpim) :: noptmfix     ! bermejo & conde fixer option for calculating its weight
integer(kind=jpim) :: noptvfe      ! use vertical fe in calculation of column mass total
integer(kind=jpim) :: npmfix       ! parameter used in weight calculation
integer(kind=jpim) :: nmfdiaglev   ! determines global diagnostic output level for fixer:
                                   ! 0 - nothing, 1 - norms printed, 2 - norms + monotonicity
integer(kind=jpim) :: nmfixflds(jpnamed_gfl+jpghg+jpgrg+jpchem+jpaero+jptrac) 
                                   ! index of fields to be corrected by mass fixers
integer(kind=jpim) :: nnegflds(jpnamed_gfl+jpghg+jpgrg+jpchem+jpaero+jptrac)  
                                   ! index of fields to be corrected by sl -ve fixer
real(kind=jprb)    :: zmfixeps     ! threshold for mass fixing scheme

type(type_gfl_comp) :: ycomp(jpgfl)    ! general descriptor of all components

type(type_gfl_comp),pointer  :: yq          => null() ! specific humidity
type(type_gfl_comp),pointer  :: yi          => null() ! ice water
type(type_gfl_comp),pointer  :: yl          => null() ! liquid water
type(type_gfl_comp),pointer  :: ylconv      => null() ! liquid water (conv. part)
type(type_gfl_comp),pointer  :: yiconv      => null() ! ice    water (conv. part)
type(type_gfl_comp),pointer  :: yrconv      => null() ! rain         (conv. part)
type(type_gfl_comp),pointer  :: ysconv      => null() ! snow         (conv. part)
type(type_gfl_comp),pointer  :: yirad       => null() ! radiative cloud ice water
type(type_gfl_comp),pointer  :: ylrad       => null() ! radiative cloud liquid water
type(type_gfl_comp),pointer  :: ys          => null() ! snow
type(type_gfl_comp),pointer  :: yr          => null() ! rain
type(type_gfl_comp),pointer  :: yg          => null() ! graupel
type(type_gfl_comp),pointer  :: yh          => null() ! hail
type(type_gfl_comp),pointer  :: ytke        => null() ! turbulent kinetic energy
type(type_gfl_comp),pointer  :: ytte        => null() ! turbulent total energy
type(type_gfl_comp),pointer  :: yefb1       => null() ! first variable efb scheme
type(type_gfl_comp),pointer  :: yefb2       => null() ! second variable efb scheme
type(type_gfl_comp),pointer  :: yefb3       => null() ! third variable efb scheme
type(type_gfl_comp),pointer  :: ya          => null() ! cloud fraction
type(type_gfl_comp),pointer  :: yo3         => null() ! ozone
type(type_gfl_comp),pointer  :: ysrc        => null() ! second-order flux for arome s'rc'/2sigma_s2 multiplied by lambda_3
type(type_gfl_comp),pointer  :: ymxl        => null() ! prognostic mixing length
type(type_gfl_comp),pointer  :: yscc2       => null() ! saturation deficit^2 for tompkins
type(type_gfl_comp),pointer  :: ygcca       => null() ! skewness for tompkins
type(type_gfl_comp),pointer  :: ycpf        => null() ! convective precipitation flux
type(type_gfl_comp),pointer  :: yspf        => null() ! stratiform precipitation flux
type(type_gfl_comp),pointer  :: ycvgq       => null() ! moisture convergence for french physics
type(type_gfl_comp),pointer  :: yqva        => null() ! total humidity variation
type(type_gfl_comp),pointer  :: yghg(:)     => null() ! greenhouse gases
type(type_gfl_comp),pointer  :: ygrg(:)     => null() ! reactive gases
type(type_gfl_comp),pointer  :: ychem(:)    => null() ! chemistry
type(type_gfl_comp),pointer  :: ygrgtend(:) => null() ! reactive gases tendecies
type(type_gfl_comp),pointer  :: yaero(:)    => null() ! aerosols
type(type_gfl_comp),pointer  :: ytrac(:)    => null() ! tracers for diagnostics
type(type_gfl_comp),pointer  :: ylrch4      => null() ! ch4 loss rate (instantaneous field)
type(type_gfl_comp),pointer  :: ych4s       => null() ! ch4 atmospheric sink (accumulated field)
type(type_gfl_comp),pointer  :: yforc(:)    => null() ! large scale forcing
type(type_gfl_comp),pointer  :: yezdiag(:)  => null() ! easy diagnostics
type(type_gfl_comp),pointer  :: yera40(:)   => null() ! era40 diagnostic fields
type(type_gfl_comp),pointer  :: ynogw(:)    => null() ! noro gwd scheme
type(type_gfl_comp),pointer  :: ysldia(:)   => null() ! sl dynamics diagnostics
type(type_gfl_comp),pointer  :: yaerout(:)  => null() ! aerosol outputs
type(type_gfl_comp),pointer  :: yuvp(:)     => null() ! uv-processor output
type(type_gfl_comp),pointer  :: yphys(:)    => null() ! phys output


type(type_gfl_comp),pointer  :: ysdsat      => null() ! standard deviation of the
                                                      ! saturation depression (sigma_s) 
type(type_gfl_comp),pointer  :: ycvv        => null() ! convective vertical velocity
type(type_gfl_comp),pointer  :: yrkth       => null() ! rasch-kristjansson h tendency
type(type_gfl_comp),pointer  :: yrktqv      => null() ! rasch-kristjansson qv tendency
type(type_gfl_comp),pointer  :: yrktqc      => null() ! rasch-kristjansson qc tendency

! prognostic convection variables: add 6 named components
type(type_gfl_comp),pointer  :: yuom        => null() ! updraught vert velocity
type(type_gfl_comp),pointer  :: yual        => null() ! updraught mesh fraction
type(type_gfl_comp),pointer  :: ydom        => null() ! downdraught vert velocity
type(type_gfl_comp),pointer  :: ydal        => null() ! downdraught mesh fraction
type(type_gfl_comp),pointer  :: yuen        => null() ! updraught entrainment
type(type_gfl_comp),pointer  :: yunebh      => null() ! pseudo-historic convective

! extra fields

type(type_gfl_comp),pointer  :: yext(:)     => null() ! extra fields

type(type_gfl_naml)  :: yq_nl                 ! specific humidity
type(type_gfl_naml)  :: yi_nl                 ! ice water
type(type_gfl_naml)  :: yl_nl                 ! liquid water
type(type_gfl_naml)  :: ylconv_nl             ! liquid water (conv. part)
type(type_gfl_naml)  :: yiconv_nl             ! ice    water (conv. part)
type(type_gfl_naml)  :: yrconv_nl             ! rain         (conv. part)
type(type_gfl_naml)  :: ysconv_nl             ! snow         (conv. part)
type(type_gfl_naml)  :: yirad_nl              ! radiative cloud ice water
type(type_gfl_naml)  :: ylrad_nl              ! radiative cloud liquid water
type(type_gfl_naml)  :: ys_nl                 ! snow
type(type_gfl_naml)  :: yr_nl                 ! rain
type(type_gfl_naml)  :: yg_nl                 ! graupels
type(type_gfl_naml)  :: yh_nl                 ! hail
type(type_gfl_naml)  :: ytke_nl               ! turbulent kinetic energy
type(type_gfl_naml)  :: ytte_nl               ! turbulent total energy
type(type_gfl_naml)  :: yefb1_nl              ! first variable efb scheme
type(type_gfl_naml)  :: yefb2_nl              ! second variable efb scheme
type(type_gfl_naml)  :: yefb3_nl              ! third variable efb scheme
type(type_gfl_naml)  :: ya_nl                 ! cloud fraction
type(type_gfl_naml)  :: yo3_nl                ! ozone
type(type_gfl_naml)  :: ysrc_nl               ! second-order flux for arome
                                              ! s'rc'/2sigma_s2
                                              ! multiplied by lambda_3
type(type_gfl_naml)  :: ymxl_nl               ! prognostic mixing length
type(type_gfl_naml)  :: yscc2_nl              ! saturation deficit^2 for tompkins
type(type_gfl_naml)  :: ygcca_nl              ! skewness for tompkins
type(type_gfl_naml)  :: ycpf_nl               ! convective precipitation flux
type(type_gfl_naml)  :: yspf_nl               ! stratiform precipitation flux
type(type_gfl_naml)  :: ycvgq_nl              ! moisture convergence for french physics
type(type_gfl_naml)  :: yqva_nl               ! total humidity variation

type(type_gfl_naml)  :: yghg_nl(jpghg)        ! greenhouse gases
type(type_gfl_naml)  :: ygrg_nl(jpgrg)        ! reactive gases
type(type_gfl_naml)  :: ychem_nl(jpchem)      ! chemical species
type(type_gfl_naml)  :: ygrgtend_nl(jpgrg)    ! reactive gases tendecies
type(type_gfl_naml)  :: yaero_nl(jpaero)      ! aerosol fields
type(type_gfl_naml)  :: ytrac_nl(jptrac)      ! tracers for diagnostics
type(type_gfl_naml)  :: yera40_nl(jpera40)    ! era40 diagnostic fields
type(type_gfl_naml)  :: ynogw_nl(jpnogw)      ! noro gwd scheme
type(type_gfl_naml)  :: ysldia_nl(jpsldia)    ! sl dynamics diagnostics
type(type_gfl_naml)  :: ylrch4_nl             ! ch4 loss rate
type(type_gfl_naml)  :: ych4s_nl              ! ch4 atmospheric sink
type(type_gfl_naml)  :: yaerout_nl(jpaerout)  ! aerosol outputs
type(type_gfl_naml)  :: yuvp_nl(jpuvp)        ! uv-processor outputs
type(type_gfl_naml)  :: yrkth_nl              ! rasch-kristjansson h tendency
type(type_gfl_naml)  :: yrktqv_nl             ! rasch-kristjansson qv tendency
type(type_gfl_naml)  :: yrktqc_nl             ! rasch-kristjansson qc tendency
type(type_gfl_naml)  :: yphys_nl(jpphys)      ! phys outputs 

! extra fields
type(type_gfl_naml)  :: ysdsat_nl             ! standard deviation of the
                                              ! saturation depression (sigma_s) 
type(type_gfl_naml)  :: ycvv_nl               ! convective vertical velocity
type(type_gfl_naml)  :: yforc_nl(jpforc)      ! forcing precursor
type(type_gfl_naml)  :: yezdiag_nl(jpezdiag)  ! easy diagnostics
type(type_gfl_naml)  :: yext_nl(jpgfl-jpnamed_gfl-jpghg-jpgrg-jpforc-jpezdiag-jpaero-jptrac-jpera40-&
 &                              jpnogw-jpsldia-jpch4s-jpaerout-jpuvp-jpchem-jpphys) ! extra fields

! prognostic convection variables: 6 more namelist components
type(type_gfl_naml)  :: yuom_nl               ! updraught vert velocity
type(type_gfl_naml)  :: yual_nl               ! updraught mesh fraction
type(type_gfl_naml)  :: ydom_nl               ! downdraught vert velocity
type(type_gfl_naml)  :: ydal_nl               ! downdraught mesh fraction
type(type_gfl_naml)  :: yuen_nl               ! updraught entrainment
type(type_gfl_naml)  :: yunebh_nl             ! pseudi hist conv cloud fraction

end type type_gfld

! gfl general descriptor
type(type_gfld), pointer :: ygfl => null()

end module yom_ygfl
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define __sizeof_int__ 4
! #define __sizeof_pointer__ 8
! #define __gfortran__ 1
! #define __gfc_real_16__ 1
! #define __stdc_hosted__ 0
! #define __no_math_errno__ 1
! #define __sizeof_float__ 4
! #define __pic__ 2
! #define _language_fortran 1
! #define __sizeof_long__ 8
! #define __gfc_int_8__ 1
! #define __dynamic__ 1
! #define __sizeof_short__ 2
! #define __gnuc__ 13
! #define __sizeof_long_double__ 16
! #define __biggest_alignment__ 16
! #define __atomic_relaxed 0
! #define _lp64 1
! #define __ecrad_little_endian 1
! #define __gfc_int_1__ 1
! #define __order_big_endian__ 4321
! #define __byte_order__ __order_little_endian__
! #define __sizeof_size_t__ 8
! #define __pic__ 2
! #define __sizeof_double__ 8
! #define __atomic_consume 1
! #define __gnuc_minor__ 3
! #define __gfc_int_16__ 1
! #define __lp64__ 1
! #define __atomic_seq_cst 5
! #define __sizeof_long_long__ 8
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

