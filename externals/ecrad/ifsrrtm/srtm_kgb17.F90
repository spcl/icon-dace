! # 1 "ifsrrtm/srtm_kgb17.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb17.f90"
subroutine srtm_kgb17

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o, co2)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  , only : jprb
use ecradhook   , only : lhook, dr_hook, jphook
use yomlun    , only : nulrad
use yommp0    , only : nproc, myproc
use mpl_module, only : mpl_broadcast
use yomtag    , only : mtagrad
use yoesrta17 , only : ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr, &
  &   ka_d, kb_d

!     ------------------------------------------------------------------

implicit none

! kurucz
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 30 "ifsrrtm/srtm_kgb17.f90" 2

if (lhook) call dr_hook('srtm_kgb17',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb17:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb17:')
endif

sfluxref(:,1) = (/ &
 & 3.15613_jprb  ,  3.03449_jprb  ,  2.92069_jprb  ,  2.63874_jprb   , &
 & 2.34581_jprb  ,  2.06999_jprb  ,  1.70906_jprb  ,  1.29085_jprb   , &
 & 0.874851_jprb ,  0.0955392_jprb,  0.0787813_jprb,  0.0621951_jprb , &
 & 0.0459076_jprb,  0.0294129_jprb,  0.0110387_jprb,  0.00159668_jprb /)  
  
sfluxref(:,2) = (/ &
 & 2.83147_jprb  ,  2.95919_jprb  ,  2.96674_jprb  ,  2.77677_jprb   , &
 & 2.46826_jprb  ,  2.11481_jprb  ,  1.73243_jprb  ,  1.30279_jprb   , &
 & 0.882714_jprb ,  0.0962350_jprb,  0.0802122_jprb,  0.0636194_jprb , &
 & 0.0472620_jprb,  0.0299051_jprb,  0.0110785_jprb,  0.00159668_jprb /)  
  
sfluxref(:,3) = (/ &
 & 2.82300_jprb  ,  2.94845_jprb  ,  2.95887_jprb  ,  2.77593_jprb   , &
 & 2.47096_jprb  ,  2.12596_jprb  ,  1.73847_jprb  ,  1.30796_jprb   , &
 & 0.884395_jprb ,  0.0966936_jprb,  0.0801996_jprb,  0.0640199_jprb , &
 & 0.0472803_jprb,  0.0300515_jprb,  0.0112366_jprb,  0.00160814_jprb /)  
  
sfluxref(:,4) = (/ &
 & 2.81715_jprb  ,  2.93789_jprb  ,  2.95091_jprb  ,  2.77046_jprb   , &
 & 2.47716_jprb  ,  2.13591_jprb  ,  1.74365_jprb  ,  1.31277_jprb   , &
 & 0.887443_jprb ,  0.0967016_jprb,  0.0803391_jprb,  0.0642442_jprb , &
 & 0.0472909_jprb,  0.0300720_jprb,  0.0114817_jprb,  0.00161875_jprb /)  
  
sfluxref(:,5) = (/ &
 & 2.82335_jprb  ,  2.93168_jprb  ,  2.91455_jprb  ,  2.75213_jprb   , &
 & 2.49168_jprb  ,  2.14408_jprb  ,  1.75726_jprb  ,  1.32401_jprb   , &
 & 0.893644_jprb ,  0.0969523_jprb,  0.0805197_jprb,  0.0639936_jprb , &
 & 0.0475099_jprb,  0.0305667_jprb,  0.0115372_jprb,  0.00161875_jprb /)  
     
!     rayleigh extinction coefficient at v = 3625 cm-1.
rayl = 6.86e-10_jprb

strrat = 0.364641_jprb

layreffr = 30

!     ------------------------------------------------------------------

!     the array ka contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  the first 
!     index in the array, js, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  for instance, 
!     js=1 refers to dry air, js = 2 corresponds to the paramter value 1/8, 
!     js = 3 corresponds to the parameter value 2/8, etc.  the second index
!     in the array, jt, which runs from 1 to 5, corresponds to different
!     temperatures.  more specifically, jt = 3 means that the data are for
!     the reference temperature tref for this  pressure level, jt = 2 refers
!     to tref-15, jt = 1 is for tref-30, jt = 4 is for tref+15, and jt = 5
!     is for tref+30.  the third index, jp, runs from 1 to 13 and refers
!     to the jpth reference pressure level (see taumol.f for these levels
!     in mb).  the fourth index, ig, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
!     the array kb contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. the first 
!     index in the array, jt, which runs from 1 to 5, corresponds to 
!     different temperatures.  more specifically, jt = 3 means that the 
!     data are for the reference temperature tref for this pressure 
!     level, jt = 2 refers to the temperature tref-15, jt = 1 is for
!     tref-30, jt = 4 is for tref+15, and jt = 5 is for tref+30.  
!     the second index, jp, runs from 13 to 59 and refers to the jpth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  the third index, ig, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.
!     -----------------------------------------------------------------


forref(:, 1) = (/ 0.553258e-03_jprb, 0.555486e-03_jprb, 0.601339e-03_jprb, 0.708280e-03_jprb /)
forref(:, 2) = (/ 0.158558e-02_jprb, 0.162957e-02_jprb, 0.204991e-02_jprb, 0.475881e-02_jprb /)
forref(:, 3) = (/ 0.772542e-02_jprb, 0.784562e-02_jprb, 0.111979e-01_jprb, 0.229016e-01_jprb /)
forref(:, 4) = (/ 0.255097e-01_jprb, 0.256272e-01_jprb, 0.270691e-01_jprb, 0.259505e-01_jprb /)
forref(:, 5) = (/ 0.323263e-01_jprb, 0.324495e-01_jprb, 0.305535e-01_jprb, 0.263993e-01_jprb /)
forref(:, 6) = (/ 0.346920e-01_jprb, 0.348255e-01_jprb, 0.323586e-01_jprb, 0.276357e-01_jprb /)
forref(:, 7) = (/ 0.366509e-01_jprb, 0.366412e-01_jprb, 0.344434e-01_jprb, 0.319223e-01_jprb /)
forref(:, 8) = (/ 0.378451e-01_jprb, 0.375341e-01_jprb, 0.374369e-01_jprb, 0.320334e-01_jprb /)
forref(:, 9) = (/ 0.407348e-01_jprb, 0.396203e-01_jprb, 0.393988e-01_jprb, 0.318343e-01_jprb /)
forref(:,10) = (/ 0.433035e-01_jprb, 0.426488e-01_jprb, 0.408085e-01_jprb, 0.332749e-01_jprb /)
forref(:,11) = (/ 0.428254e-01_jprb, 0.441151e-01_jprb, 0.408887e-01_jprb, 0.327077e-01_jprb /)
forref(:,12) = (/ 0.443226e-01_jprb, 0.446690e-01_jprb, 0.404676e-01_jprb, 0.350492e-01_jprb /)
forref(:,13) = (/ 0.466103e-01_jprb, 0.460809e-01_jprb, 0.401286e-01_jprb, 0.370427e-01_jprb /)
forref(:,14) = (/ 0.483928e-01_jprb, 0.477284e-01_jprb, 0.380684e-01_jprb, 0.387940e-01_jprb /)
forref(:,15) = (/ 0.506987e-01_jprb, 0.490016e-01_jprb, 0.467069e-01_jprb, 0.368998e-01_jprb /)
forref(:,16) = (/ 0.510836e-01_jprb, 0.522771e-01_jprb, 0.500130e-01_jprb, 0.483406e-01_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

selfref(:, 1) = (/ &
 & 0.160537e-01_jprb, 0.149038e-01_jprb, 0.138363e-01_jprb, 0.128452e-01_jprb, 0.119251e-01_jprb, &
 & 0.110709e-01_jprb, 0.102779e-01_jprb, 0.954175e-02_jprb, 0.885829e-02_jprb, 0.822379e-02_jprb /)  
selfref(:, 2) = (/ &
 & 0.365753e-01_jprb, 0.342267e-01_jprb, 0.320288e-01_jprb, 0.299720e-01_jprb, 0.280474e-01_jprb, &
 & 0.262463e-01_jprb, 0.245609e-01_jprb, 0.229837e-01_jprb, 0.215078e-01_jprb, 0.201267e-01_jprb /)  
selfref(:, 3) = (/ &
 & 0.127419e+00_jprb, 0.118553e+00_jprb, 0.110304e+00_jprb, 0.102629e+00_jprb, 0.954883e-01_jprb, &
 & 0.888442e-01_jprb, 0.826624e-01_jprb, 0.769107e-01_jprb, 0.715593e-01_jprb, 0.665802e-01_jprb /)  
selfref(:, 4) = (/ &
 & 0.378687e+00_jprb, 0.348961e+00_jprb, 0.321568e+00_jprb, 0.296325e+00_jprb, 0.273064e+00_jprb, &
 & 0.251629e+00_jprb, 0.231876e+00_jprb, 0.213674e+00_jprb, 0.196901e+00_jprb, 0.181444e+00_jprb /)  
selfref(:, 5) = (/ &
 & 0.472822e+00_jprb, 0.435018e+00_jprb, 0.400236e+00_jprb, 0.368236e+00_jprb, 0.338794e+00_jprb, &
 & 0.311706e+00_jprb, 0.286783e+00_jprb, 0.263854e+00_jprb, 0.242757e+00_jprb, 0.223348e+00_jprb /)  
selfref(:, 6) = (/ &
 & 0.505620e+00_jprb, 0.465050e+00_jprb, 0.427736e+00_jprb, 0.393416e+00_jprb, 0.361849e+00_jprb, &
 & 0.332815e+00_jprb, 0.306111e+00_jprb, 0.281550e+00_jprb, 0.258959e+00_jprb, 0.238181e+00_jprb /)  
selfref(:, 7) = (/ &
 & 0.530488e+00_jprb, 0.487993e+00_jprb, 0.448902e+00_jprb, 0.412943e+00_jprb, 0.379864e+00_jprb, &
 & 0.349434e+00_jprb, 0.321443e+00_jprb, 0.295694e+00_jprb, 0.272007e+00_jprb, 0.250218e+00_jprb /)  
selfref(:, 8) = (/ &
 & 0.540222e+00_jprb, 0.497746e+00_jprb, 0.458610e+00_jprb, 0.422551e+00_jprb, 0.389327e+00_jprb, &
 & 0.358716e+00_jprb, 0.330511e+00_jprb, 0.304524e+00_jprb, 0.280580e+00_jprb, 0.258519e+00_jprb /)  
selfref(:, 9) = (/ &
 & 0.565727e+00_jprb, 0.522899e+00_jprb, 0.483313e+00_jprb, 0.446724e+00_jprb, 0.412905e+00_jprb, &
 & 0.381646e+00_jprb, 0.352753e+00_jprb, 0.326048e+00_jprb, 0.301365e+00_jprb, 0.278550e+00_jprb /)  
selfref(:,10) = (/ &
 & 0.610122e+00_jprb, 0.562337e+00_jprb, 0.518295e+00_jprb, 0.477702e+00_jprb, 0.440289e+00_jprb, &
 & 0.405806e+00_jprb, 0.374023e+00_jprb, 0.344730e+00_jprb, 0.317730e+00_jprb, 0.292846e+00_jprb /)  
selfref(:,11) = (/ &
 & 0.645176e+00_jprb, 0.588957e+00_jprb, 0.537636e+00_jprb, 0.490788e+00_jprb, 0.448022e+00_jprb, &
 & 0.408982e+00_jprb, 0.373344e+00_jprb, 0.340812e+00_jprb, 0.311114e+00_jprb, 0.284004e+00_jprb /)  
selfref(:,12) = (/ &
 & 0.651737e+00_jprb, 0.596547e+00_jprb, 0.546031e+00_jprb, 0.499792e+00_jprb, 0.457469e+00_jprb, &
 & 0.418730e+00_jprb, 0.383272e+00_jprb, 0.350816e+00_jprb, 0.321108e+00_jprb, 0.293916e+00_jprb /)  
selfref(:,13) = (/ &
 & 0.661086e+00_jprb, 0.607954e+00_jprb, 0.559093e+00_jprb, 0.514159e+00_jprb, 0.472836e+00_jprb, &
 & 0.434834e+00_jprb, 0.399886e+00_jprb, 0.367747e+00_jprb, 0.338191e+00_jprb, 0.311011e+00_jprb /)  
selfref(:,14) = (/ &
 & 0.692554e+00_jprb, 0.635574e+00_jprb, 0.583282e+00_jprb, 0.535293e+00_jprb, 0.491251e+00_jprb, &
 & 0.450834e+00_jprb, 0.413741e+00_jprb, 0.379701e+00_jprb, 0.348461e+00_jprb, 0.319791e+00_jprb /)  
selfref(:,15) = (/ &
 & 0.714646e+00_jprb, 0.657179e+00_jprb, 0.604334e+00_jprb, 0.555737e+00_jprb, 0.511049e+00_jprb, &
 & 0.469954e+00_jprb, 0.432164e+00_jprb, 0.397412e+00_jprb, 0.365455e+00_jprb, 0.336068e+00_jprb /)  
selfref(:,16) = (/ &
 & 0.782126e+00_jprb, 0.710682e+00_jprb, 0.645764e+00_jprb, 0.586776e+00_jprb, 0.533177e+00_jprb, &
 & 0.484473e+00_jprb, 0.440219e+00_jprb, 0.400007e+00_jprb, 0.363468e+00_jprb, 0.330266e+00_jprb /)  

if (lhook) call dr_hook('srtm_kgb17',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb17:error reading file radsrtm")

end subroutine srtm_kgb17
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

