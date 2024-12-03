! # 1 "ifsrrtm/srtm_kgb19.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb19.f90"
subroutine srtm_kgb19

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 16:  4650-5150 cm-1 (low - h2o,co2; high - co2)
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
use yoesrta19 , only : ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr, &
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
! # 30 "ifsrrtm/srtm_kgb19.f90" 2

if (lhook) call dr_hook('srtm_kgb19',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb19:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb19:')
endif

sfluxref(:,1) = (/ &
 & 3.25791_jprb    , 3.29697_jprb    , 3.16031_jprb    , 2.96115_jprb    , &
 & 2.69238_jprb    , 2.33819_jprb    , 1.92760_jprb    , 1.44918_jprb    , &
 & 0.979764_jprb   , 0.107336_jprb   , 8.94523e-02_jprb, 6.98325e-02_jprb, &
 & 5.12051e-02_jprb, 3.23645e-02_jprb, 1.23401e-02_jprb, 1.71339e-03_jprb /)  
sfluxref(:,2) = (/ &
 & 3.22769_jprb    , 3.28817_jprb    , 3.16687_jprb    , 2.97662_jprb    , &
 & 2.69495_jprb    , 2.34392_jprb    , 1.92900_jprb    , 1.45391_jprb    , &
 & 0.982522_jprb   , 0.107638_jprb   , 8.92458e-02_jprb, 6.99885e-02_jprb, &
 & 5.09679e-02_jprb, 3.23789e-02_jprb, 1.22673e-02_jprb, 1.56040e-03_jprb /)  
sfluxref(:,3) = (/ &
 & 3.22294_jprb    , 3.27780_jprb    , 3.17424_jprb    , 2.97143_jprb    , &
 & 2.69785_jprb    , 2.34993_jprb    , 1.93155_jprb    , 1.45196_jprb    , &
 & 0.985329_jprb   , 0.108027_jprb   , 8.93552e-02_jprb, 6.99937e-02_jprb, &
 & 5.11678e-02_jprb, 3.24846e-02_jprb, 1.20636e-02_jprb, 1.56040e-03_jprb /)  
sfluxref(:,4) = (/ &
 & 3.22445_jprb    , 3.26113_jprb    , 3.18438_jprb    , 2.96921_jprb    , &
 & 2.69579_jprb    , 2.35586_jprb    , 1.93454_jprb    , 1.44949_jprb    , &
 & 0.987347_jprb   , 0.108611_jprb   , 8.91643e-02_jprb, 7.02236e-02_jprb, &
 & 5.12980e-02_jprb, 3.25282e-02_jprb, 1.21189e-02_jprb, 1.56040e-03_jprb /)  
sfluxref(:,5) = (/ &
 & 3.22497_jprb    , 3.25109_jprb    , 3.18741_jprb    , 2.96970_jprb    , &
 & 2.69460_jprb    , 2.36020_jprb    , 1.93301_jprb    , 1.45224_jprb    , &
 & 0.988564_jprb   , 0.108255_jprb   , 8.93830e-02_jprb, 7.03655e-02_jprb, &
 & 5.13017e-02_jprb, 3.29414e-02_jprb, 1.21189e-02_jprb, 1.56040e-03_jprb /)  
sfluxref(:,6) = (/ &
 & 3.22632_jprb    , 3.24174_jprb    , 3.18524_jprb    , 2.97402_jprb    , &
 & 2.69807_jprb    , 2.35742_jprb    , 1.93377_jprb    , 1.45621_jprb    , &
 & 0.988132_jprb   , 0.108344_jprb   , 8.93188e-02_jprb, 7.04907e-02_jprb, &
 & 5.17938e-02_jprb, 3.31465e-02_jprb, 1.21155e-02_jprb, 1.56040e-03_jprb /)  
sfluxref(:,7) = (/ &
 & 3.22793_jprb    , 3.23589_jprb    , 3.17720_jprb    , 2.97869_jprb    , &
 & 2.70293_jprb    , 2.35436_jprb    , 1.93557_jprb    , 1.45868_jprb    , &
 & 0.988654_jprb   , 0.108198_jprb   , 8.93375e-02_jprb, 7.09790e-02_jprb, &
 & 5.24733e-02_jprb, 3.31298e-02_jprb, 1.21126e-02_jprb, 1.56040e-03_jprb /)  
sfluxref(:,8) = (/ &
 & 3.22966_jprb    , 3.24087_jprb    , 3.15676_jprb    , 2.98171_jprb    , &
 & 2.70894_jprb    , 2.34975_jprb    , 1.93855_jprb    , 1.46354_jprb    , &
 & 0.988544_jprb   , 0.108574_jprb   , 9.02522e-02_jprb, 7.12908e-02_jprb, &
 & 5.24844e-02_jprb, 3.31084e-02_jprb, 1.21060e-02_jprb, 1.56040e-03_jprb /)  
sfluxref(:,9) = (/ &
 & 3.27240_jprb    , 3.24666_jprb    , 3.13886_jprb    , 2.95238_jprb    , &
 & 2.70190_jprb    , 2.34460_jprb    , 1.93948_jprb    , 1.47111_jprb    , &
 & 0.990821_jprb   , 0.108730_jprb   , 9.01625e-02_jprb, 7.13261e-02_jprb, &
 & 5.24813e-02_jprb, 3.31083e-02_jprb, 1.21126e-02_jprb, 1.56040e-03_jprb /)  

!     rayleigh extinction coefficient at v = 4900 cm-1.
rayl = 2.29e-09_jprb

strrat = 5.49281_jprb

layreffr = 3

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

  
forref(:, 1) = (/ 0.106275e-05_jprb, 0.104185e-05_jprb, 0.420154e-05_jprb /)
forref(:, 2) = (/ 0.154343e-05_jprb, 0.653193e-05_jprb, 0.174596e-04_jprb /)
forref(:, 3) = (/ 0.348917e-05_jprb, 0.108420e-04_jprb, 0.540849e-04_jprb /)
forref(:, 4) = (/ 0.145822e-04_jprb, 0.156027e-04_jprb, 0.881263e-04_jprb /)
forref(:, 5) = (/ 0.220204e-04_jprb, 0.819892e-04_jprb, 0.817937e-04_jprb /)
forref(:, 6) = (/ 0.447840e-04_jprb, 0.121116e-03_jprb, 0.932635e-04_jprb /)
forref(:, 7) = (/ 0.166516e-03_jprb, 0.147640e-03_jprb, 0.754029e-04_jprb /)
forref(:, 8) = (/ 0.234756e-03_jprb, 0.145934e-03_jprb, 0.771734e-04_jprb /)
forref(:, 9) = (/ 0.289207e-03_jprb, 0.146768e-03_jprb, 0.677806e-04_jprb /)
forref(:,10) = (/ 0.334959e-03_jprb, 0.125513e-03_jprb, 0.636648e-04_jprb /)
forref(:,11) = (/ 0.333755e-03_jprb, 0.136575e-03_jprb, 0.593651e-04_jprb /)
forref(:,12) = (/ 0.340042e-03_jprb, 0.116259e-03_jprb, 0.595192e-04_jprb /)
forref(:,13) = (/ 0.422470e-03_jprb, 0.148691e-03_jprb, 0.630266e-04_jprb /)
forref(:,14) = (/ 0.440655e-03_jprb, 0.461917e-04_jprb, 0.108222e-04_jprb /)
forref(:,15) = (/ 0.486207e-03_jprb, 0.428458e-03_jprb, 0.108086e-04_jprb /)
forref(:,16) = (/ 0.657463e-03_jprb, 0.657446e-03_jprb, 0.126190e-04_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).
     
selfref(:, 1) = (/ &
 & 0.331728e-03_jprb, 0.287480e-03_jprb, 0.249135e-03_jprb, 0.215904e-03_jprb, 0.187106e-03_jprb, &
 & 0.162149e-03_jprb, 0.140520e-03_jprb, 0.121777e-03_jprb, 0.105534e-03_jprb, 0.914573e-04_jprb /)  
selfref(:, 2) = (/ &
 & 0.882628e-03_jprb, 0.698914e-03_jprb, 0.553439e-03_jprb, 0.438244e-03_jprb, 0.347026e-03_jprb, &
 & 0.274795e-03_jprb, 0.217598e-03_jprb, 0.172306e-03_jprb, 0.136442e-03_jprb, 0.108042e-03_jprb /)  
selfref(:, 3) = (/ &
 & 0.115461e-02_jprb, 0.937203e-03_jprb, 0.760730e-03_jprb, 0.617486e-03_jprb, 0.501215e-03_jprb, &
 & 0.406837e-03_jprb, 0.330231e-03_jprb, 0.268049e-03_jprb, 0.217576e-03_jprb, 0.176607e-03_jprb /)  
selfref(:, 4) = (/ &
 & 0.103450e-02_jprb, 0.960268e-03_jprb, 0.891360e-03_jprb, 0.827397e-03_jprb, 0.768024e-03_jprb, &
 & 0.712911e-03_jprb, 0.661754e-03_jprb, 0.614267e-03_jprb, 0.570188e-03_jprb, 0.529272e-03_jprb /)  
selfref(:, 5) = (/ &
 & 0.289040e-02_jprb, 0.240129e-02_jprb, 0.199495e-02_jprb, 0.165737e-02_jprb, 0.137692e-02_jprb, &
 & 0.114392e-02_jprb, 0.950351e-03_jprb, 0.789535e-03_jprb, 0.655933e-03_jprb, 0.544938e-03_jprb /)  
selfref(:, 6) = (/ &
 & 0.361772e-02_jprb, 0.306611e-02_jprb, 0.259861e-02_jprb, 0.220239e-02_jprb, 0.186659e-02_jprb, &
 & 0.158198e-02_jprb, 0.134077e-02_jprb, 0.113634e-02_jprb, 0.963078e-03_jprb, 0.816234e-03_jprb /)  
selfref(:, 7) = (/ &
 & 0.329878e-02_jprb, 0.318245e-02_jprb, 0.307021e-02_jprb, 0.296194e-02_jprb, 0.285749e-02_jprb, &
 & 0.275671e-02_jprb, 0.265950e-02_jprb, 0.256571e-02_jprb, 0.247522e-02_jprb, 0.238793e-02_jprb /)  
selfref(:, 8) = (/ &
 & 0.293562e-02_jprb, 0.300077e-02_jprb, 0.306737e-02_jprb, 0.313544e-02_jprb, 0.320503e-02_jprb, &
 & 0.327615e-02_jprb, 0.334886e-02_jprb, 0.342318e-02_jprb, 0.349915e-02_jprb, 0.357680e-02_jprb /)  
selfref(:, 9) = (/ &
 & 0.281453e-02_jprb, 0.295894e-02_jprb, 0.311076e-02_jprb, 0.327038e-02_jprb, 0.343818e-02_jprb, &
 & 0.361459e-02_jprb, 0.380006e-02_jprb, 0.399504e-02_jprb, 0.420002e-02_jprb, 0.441553e-02_jprb /)  
selfref(:,10) = (/ &
 & 0.239488e-02_jprb, 0.262487e-02_jprb, 0.287696e-02_jprb, 0.315325e-02_jprb, 0.345607e-02_jprb, &
 & 0.378798e-02_jprb, 0.415176e-02_jprb, 0.455048e-02_jprb, 0.498749e-02_jprb, 0.546647e-02_jprb /)  
selfref(:,11) = (/ &
 & 0.271001e-02_jprb, 0.292235e-02_jprb, 0.315134e-02_jprb, 0.339826e-02_jprb, 0.366453e-02_jprb, &
 & 0.395167e-02_jprb, 0.426131e-02_jprb, 0.459521e-02_jprb, 0.495527e-02_jprb, 0.534354e-02_jprb /)  
selfref(:,12) = (/ &
 & 0.206702e-02_jprb, 0.232254e-02_jprb, 0.260966e-02_jprb, 0.293226e-02_jprb, 0.329475e-02_jprb, &
 & 0.370204e-02_jprb, 0.415969e-02_jprb, 0.467391e-02_jprb, 0.525169e-02_jprb, 0.590090e-02_jprb /)  
selfref(:,13) = (/ &
 & 0.227023e-02_jprb, 0.257331e-02_jprb, 0.291685e-02_jprb, 0.330626e-02_jprb, 0.374766e-02_jprb, &
 & 0.424799e-02_jprb, 0.481511e-02_jprb, 0.545794e-02_jprb, 0.618660e-02_jprb, 0.701253e-02_jprb /)  
selfref(:,14) = (/ &
 & 0.851078e-03_jprb, 0.111512e-02_jprb, 0.146109e-02_jprb, 0.191439e-02_jprb, 0.250832e-02_jprb, &
 & 0.328653e-02_jprb, 0.430617e-02_jprb, 0.564215e-02_jprb, 0.739261e-02_jprb, 0.968616e-02_jprb /)  
selfref(:,15) = (/ &
 & 0.742711e-02_jprb, 0.721347e-02_jprb, 0.700598e-02_jprb, 0.680446e-02_jprb, 0.660873e-02_jprb, &
 & 0.641863e-02_jprb, 0.623400e-02_jprb, 0.605468e-02_jprb, 0.588052e-02_jprb, 0.571137e-02_jprb /)  
selfref(:,16) = (/ &
 & 0.107170e-01_jprb, 0.101913e-01_jprb, 0.969138e-02_jprb, 0.921599e-02_jprb, 0.876392e-02_jprb, &
 & 0.833402e-02_jprb, 0.792521e-02_jprb, 0.753646e-02_jprb, 0.716677e-02_jprb, 0.681522e-02_jprb /)  

if (lhook) call dr_hook('srtm_kgb19',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb19:error reading file radsrtm")

end subroutine srtm_kgb19
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

