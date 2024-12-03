! # 1 "ifsrrtm/yoerad.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoerad.f90"
module yoerad

use parkind1  ,only : jpim     ,jprb

implicit none

save

!     ------------------------------------------------------------------
!*    ** *yoerad* - control options for radiation configuration
!     ------------------------------------------------------------------

integer(kind=jpim) :: naer
integer(kind=jpim) :: nmode
integer(kind=jpim) :: nozocl
integer(kind=jpim) :: nradfr
integer(kind=jpim) :: nradpfr
integer(kind=jpim) :: nradpla
integer(kind=jpim) :: nradint
integer(kind=jpim) :: nradres
integer(kind=jpim) :: nradnfr
integer(kind=jpim) :: nradsfr
integer(kind=jpim) :: nrade1h, nrade3h
integer(kind=jpim) :: nradelg
integer(kind=jpim) :: novlp
integer(kind=jpim) :: nrproma
integer(kind=jpim) :: nsw
integer(kind=jpim) :: nswnl
integer(kind=jpim) :: nswtl
integer(kind=jpim) :: ntsw
integer(kind=jpim) :: nuv
integer(kind=jpim) :: ncsradf
integer(kind=jpim) :: niceopt
integer(kind=jpim) :: nliqopt
integer(kind=jpim) :: nradip
integer(kind=jpim) :: nradlp
integer(kind=jpim) :: ninhom
integer(kind=jpim) :: nlayinh
integer(kind=jpim) :: nlngr1h
integer(kind=jpim) :: npertaer
integer(kind=jpim) :: npertoz
integer(kind=jpim) :: nscen
integer(kind=jpim) :: nhincsol
integer(kind=jpim) :: nmcica
integer(kind=jpim) :: nghgrad
integer(kind=jpim) :: ndecolat
integer(kind=jpim) :: nminice
integer(kind=jpim) :: nvolcvert
integer(kind=jpim) :: nredglw
integer(kind=jpim) :: nredgsw
integer(kind=jpim) :: nspmapl(16), nspmaps(14)

logical :: lerad1h
logical :: lepo3ra
logical :: lonewsw
logical :: lecsrad
logical :: lrrtm
logical :: lsrtm
logical :: ldiffc
logical :: lhvolca
logical :: lnewaer
logical :: lnotroaer
logical :: lrayl
logical :: loptrproma
logical :: leco2var
logical :: lhghg
logical :: lemodal
logical :: leso4his
logical :: letracgms
logical :: laerclim, laervisi
logical :: lvolcspec
logical :: lvolcdamp
logical :: ldiagforcing
logical :: lapproxlwupdate
logical :: lapproxswupdate
logical :: lcentredtimesza

character (len = 256) ::  crtabledir
character (len = 32) ::   crtablefil
logical :: lccnl
logical :: lccno
logical :: lperpet

real(kind=jprb) :: raovlp , rbovlp
real(kind=jprb) :: rccnlnd, rccnsea
logical :: ledbug
real(kind=jprb) :: rpertoz, rre2de
real(kind=jprb) :: rlwinhf, rswinhf
real(kind=jprb) :: rminice
real(kind=jprb) :: rvolcspec(3)
real(kind=jprb) :: rns, rsigair

!        * e.c.m.w.f. physics package *

!     j.-j. morcrette       e.c.m.w.f.      89/07/14
! modifications
!    r j hogan 20 may  2014: added lapproxlwupdate
!    r j hogan 19 june 2014: added lapproxswupdate
!    r j hogan 19 nov  2014: added lcentredtimesza

!  name     type     purpose
!  ----  :  ----   : ---------------------------------------------------
! lerad1h: logical : .t. to allow more frequent radiation calculations
!                  : during first n hours of forecast
! nlngr1h: integer : number forecast hours during which more frequent
!                    radiation calculations are required
! lepo3ra: logical : .t. if prognostic ozone (ec) is passed to radiation
! naer   : integer : configuration index for aerosols
! nmode  : integer : configuration for radiation code: flux vs. radiance
! nozocl : integer : choice of ozone climatology (0 old, 1 new)
! nradfr : integer : frequency of full radiation computations
!                    if(nradfr.gt.0): rad every 'nradfr' time-steps
!                    if(nradfr.lt.0): rad every '-nradfr' hours
! nradpfr: integer : print frequency for rad.statistics (in rad.t.steps)
! nradpla: integer : print rad.statistics every 'nradpla' rows
! nradint: integer : radiation interpolation method
!                  : 1 = spectral transform interpolation
!                  : 2 =  4 point horizontal interpolation
!                  : 3 = 12 point horizontal interpolation
! nradres: integer : radiation grid spectral resolution
! nradnfr: integer : normal   frequency of radiation steps
! nradsfr: integer : start-up frequency of radiation steps
! nrade1h: integer : start-up frequency of radiation steps for eps
! nrade3h: integer : subsequent frequency of radiation steps for eps
! nradelg: integer : length in hours during which the frequency of radiation is increased for eps
! novlp  : integer : cloud overlap configuration
! nrproma: integer : vector length for radiation calculations
! nsw    : integer : number of shortwave spectral intervals
! nswnl  : integer : number of shortwave spectral intervals in nl model
! nswtl  : integer : number of shortwave spectral intervals in tl model
! ntsw   : integer : maximum possible number of sw spectral intervals 
! nuv    : integer : number of uv spectral intervals for the uv processor   
! loptrproma:logical: .t. nrproma will be optimised
!                   : .f. nrproma will not be optimised (forced
!                   :         by negative nrproma in namelist)

! nradip : integer : index for diagnosis of ice cloud effective radius
!          0=ebcu/smsh  1=ebcu/ebcu  2=fuli/fuli  3=fu/fu&al
! nradlp : integer : index for diagnosis of liq. cloud effective radius
!          0=yf/smsh    1=asl/hsa    2=asl/lili
! niceopt: integer : index for ice cloud optical properties
!          0=40u        1=40-130     2=30-60      3=sun'01
! nliqopt: integer : index for liquid water cloud optical properties
!          0=f(p)       1=10/13      2=martin_et_al

! lonewsw: logical : .t. if new sw code is active
! lecsrad: logical : .t. if clear-sky radiation is archived as pextr2
! ncsradf: integer : 1 if accumulated, 2 if instantaneous
! lrrtm  : logical : .t. if rrtm140mr is used for lw radiation transfer

! lhvolca: logical : .t. if giss history of volcanic aerosols is on
! lnewaer: logical : .t. if aerosol monthly distributions are used
! lnotroaer:logical: .t. if no tropospheric aerosols
! crtabledir: char : if nradint > 0 specifies directory path for radiation
!                  : grid rtable namelist
! crtablefil: char : if nradint > 0 specifies file name of radiation 
!                  : grid rtable namelist
! lrayl  : logical : .t. new rayleigh for sw-6 version

! raovlp : real    : coefficients for alpha1 factor in hogan & 
! rbovlp : real    : illingworth's parametrization

! lccnl  : logical : .t. if ccn concentration over land is diagnosed
! lccno  : logical : .t. if ccn concentration over ocean is diagnosed
! rccnlnd: real    : number concentration (cm-3) of ccns over land
! rccnsea: real    : number concentration (cm-3) of ccns over sea

! ldiffc : logical : .t. if savijarvi's diffusivity correction is on

! ninhom : integer : 0 if no inhomogeneity scaling effect 
!                    1 if simple 0.7 scaling
!                    2 if barker, 3 if cairns et al.
! rlwinhf: real    : inhomog. scaling factor for cloud lw optical thickness
! rswinhf: real    : inhomog. scaling factor for cloud sw optical thickness

! npertaer : interger : percentage of perturbation for aerosol   
! npertozone : integer : percentage of perturbation for ozone 
! nhincsol:integer :
!        = 0 no variability of solar constant is accounted for 
!        = 1 if year-to-year variability of solar constant is accounted for 
!        = 2 if month-to-month variability of solar constant is accounted for 
!        = 3 if year-to-year variability of solar constant is accounted for according to cmip5 recommendations
! leco2var: logical: .t. if era-40/amip2 variability of ghg is on
! lhghg  : logical : .t. if variability of greenhouse gases (including co2) is on
! n.b.: lhghg supercedes leco2var and allows using better specification of trace gases
! nscen  : integer : 21st century scenario for ghg (1=a1b, 2=a2, 3=b1)
! rre2de : real    : conversion factor betwwen effective radius and particle size
! rminice: real    : minimum size for ice particles (um)
!                    for ice
! nminice: integer : 1-6 minimum ice particle size depends on latitude, 0=independent of latitude
! ndecolat:integer : decorrelation length for cf and cw 
!                     0: specified independent of latitude, 1: shonk-hogan, 2: improved
! nmcica : integer :  0: no mcica
!                     1: mcica w maximum-random in cloud generator
!                     2: mcica w generalized overlap in cloud generator
! leso4his: logical:.t.: use historical/projected so4 data per decade and month
! nghgrad: integer : configuration of 3d ghg climatologies accounted for in radiation
!                     0: global values
!                     1: co2       2: ch4    3: n2o    4: no2    5:cfc11   6:cfc12
!                    12: co2+ch4  13: co2+ch4+n2o     
!                    16: co2+ch4+n2o+cfc11+cfc12
! letracgms: logical : f=cariolle climatol. t=gems-derived clim for co2, ch4, o3
! laerclim : logical : .t. for output of the climatological aerosol optical depth at 550 nm
! laervisi : logical : .t. for output of the visibility (from diagnsotic or prognostic aerosols)
! nvolcvert: integer : vertical distribution of volcanic aerosol
!                       0: original profile, diagnosed from t
!                       1: original profile, but upper boundary at 10hpa
!                       2: lower boundary diagnosed from ozone, upper boundary at 10hpa
! lvolcspec: logical : t for specified volcanic aerosol
! lvolcdamp: logical : t for damping of specified volcanic aerosol from initial value
! rvolcspec: real    : specified volcanic aerosol (total optical depth) in nh/tropics/sh
! rns                : derived from avogadro
! rsigair: invariant terms in expression of rayleigh scattering cross-section
! nredgsw  : logical : 0 full resolution for rrtm_sw (224)
!                      1 ecmwf high resolution model configuration (_sw: 112)
!                      2 ecmwf eps configuration (_sw: 56)
! nredglw  : logical : 0 full resolution for rrtm_lw (256)
!                      1 ecmwf high resolution model configuration (_lw: 140)
!                      2 ecmwf eps configuration (_lw: 70)
! ldiagforcing : logical : t write input ozone, ghg and aerosol forcing to 3d fields 
!                            to be used for diagnostics only; do not use in production runs
! lapproxlwupdate : logical : update the longwave upwelling flux every
!                             timestep/gridpoint using the stored rate
!                             of change of the fluxes with respect to
!                             the surface upwelling longwave flux
! lapproxswupdate : logical : update the shortwave upwelling flux
!                             every gridpoint to account for the local
!                             value of surface albedo, and every
!                             timestep using manners et al. (2009)
!                             correction for solar zenith angle change
! lcentredtimesza : logical : compute solar zenith angle in radiation
!                             scheme half way between calls to
!                             radiation scheme (rather than previous
!                             behaviour, which is half way between
!                             calls plus half a model timestep)
! ------------------------------------------------------------------
end module yoerad
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

