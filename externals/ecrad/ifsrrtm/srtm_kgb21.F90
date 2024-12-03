! # 1 "ifsrrtm/srtm_kgb21.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb21.f90"
subroutine srtm_kgb21

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
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
use yoesrta21 , only : ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr  ,&
  &  ka_d, kb_d

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
! # 30 "ifsrrtm/srtm_kgb21.f90" 2

if (lhook) call dr_hook('srtm_kgb21',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb21:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb21:')
endif

sfluxref(:, 1) = (/ &
 & 16.1643_jprb , 15.5806_jprb, 14.7254_jprb    , 13.5541_jprb    , &
 & 11.9519_jprb ,10.44410_jprb, 8.37884_jprb    , 6.26384_jprb    , &
 & 4.28435_jprb ,0.465228_jprb, 0.385095_jprb   ,0.304226_jprb    , &
 & 0.222479_jprb,0.143286_jprb, 5.58046e-02_jprb, 7.84856e-03_jprb /)  
sfluxref(:, 2) = (/ &
 & 15.6451_jprb , 15.3170_jprb, 14.6987_jprb    , 13.7350_jprb    , &
 & 12.2267_jprb ,10.51646_jprb, 8.47150_jprb    , 6.38873_jprb    , &
 & 4.33536_jprb ,0.470610_jprb,0.389426_jprb    ,0.306461_jprb    , &
 & 0.223537_jprb,0.143273_jprb, 5.58179e-02_jprb, 7.84856e-03_jprb /)  
sfluxref(:, 3) = (/ &
 & 15.6092_jprb , 15.3293_jprb, 14.6881_jprb    , 13.6693_jprb    , &
 & 12.2342_jprb ,10.52010_jprb, 8.49442_jprb    , 6.42138_jprb    , &
 & 4.35865_jprb ,0.473349_jprb,0.391349_jprb    ,0.308861_jprb    , &
 & 0.224666_jprb,0.144799_jprb, 5.58176e-02_jprb, 7.84881e-03_jprb /)  
sfluxref(:, 4) = (/ &
 & 15.5786_jprb , 15.3422_jprb, 14.6894_jprb    , 13.6040_jprb    , &
 & 12.2567_jprb ,10.49400_jprb, 8.53521_jprb    , 6.44427_jprb    , &
 & 4.37208_jprb ,0.475709_jprb,0.392956_jprb    ,0.309737_jprb    , &
 & 0.226274_jprb,0.146483_jprb, 5.59325e-02_jprb, 7.84881e-03_jprb /)  
sfluxref(:, 5) = (/ &
 & 15.5380_jprb , 15.3826_jprb, 14.6575_jprb    , 13.5722_jprb    , &
 & 12.2646_jprb ,10.47672_jprb, 8.57158_jprb    , 6.46343_jprb    , &
 & 4.38259_jprb ,0.477647_jprb,0.393982_jprb    ,0.310686_jprb    , &
 & 0.227620_jprb,0.148376_jprb, 5.60398e-02_jprb, 7.83925e-03_jprb /)  
sfluxref(:, 6) = (/ &
 & 15.5124_jprb , 15.3986_jprb, 14.6240_jprb    , 13.5535_jprb    , &
 & 12.2468_jprb ,10.48891_jprb, 8.60434_jprb    , 6.47985_jprb    , &
 & 4.39448_jprb ,0.478267_jprb,0.395618_jprb    ,0.311043_jprb    , &
 & 0.230927_jprb,0.148774_jprb, 5.61189e-02_jprb, 7.83925e-03_jprb /)  
sfluxref(:, 7) = (/ &
 & 15.4910_jprb , 15.4028_jprb, 14.5772_jprb    , 13.5507_jprb    , &
 & 12.2122_jprb ,10.52735_jprb, 8.62650_jprb    , 6.49644_jprb    , &
 & 4.41173_jprb ,0.478627_jprb,0.396433_jprb    ,0.314199_jprb    ,  &
 & 0.233125_jprb,0.149052_jprb, 5.62309e-02_jprb, 7.83925e-03_jprb /)  
sfluxref(:, 8) = (/ &
 & 15.4562_jprb , 15.3928_jprb, 14.5510_jprb    , 13.5122_jprb    , &
 & 12.1890_jprb , 10.5826_jprb, 8.65842_jprb    , 6.51558_jprb    , &
 & 4.42747_jprb ,0.480669_jprb,0.400143_jprb    ,0.318144_jprb    , &
 & 0.233937_jprb,0.149119_jprb, 5.62309e-02_jprb, 7.83925e-03_jprb /)  
sfluxref(:, 9) = (/ &
 & 15.0069_jprb , 15.1479_jprb, 14.7802_jprb    , 13.6085_jprb    , &
 & 12.2793_jprb , 10.6929_jprb, 8.72723_jprb    , 6.57114_jprb    , &
 & 4.46330_jprb ,0.486724_jprb,0.401446_jprb    ,0.318879_jprb    , &
 & 0.233959_jprb,0.149119_jprb, 5.62309e-02_jprb, 7.83925e-03_jprb /)  

!     rayleigh extinction coefficient at v = 6925 cm-1.
rayl = 9.41e-09_jprb

strrat = 0.0045321_jprb

layreffr = 8

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


forref(:, 1) = (/ 0.110008e-06_jprb, 0.630912e-06_jprb, 0.363159e-05_jprb, 0.616892e-05_jprb /)
forref(:, 2) = (/ 0.429709e-05_jprb, 0.789174e-05_jprb, 0.217416e-04_jprb, 0.639393e-04_jprb /)
forref(:, 3) = (/ 0.436283e-04_jprb, 0.526247e-04_jprb, 0.116341e-03_jprb, 0.205616e-03_jprb /)
forref(:, 4) = (/ 0.215627e-03_jprb, 0.234522e-03_jprb, 0.280497e-03_jprb, 0.838668e-03_jprb /)
forref(:, 5) = (/ 0.529283e-03_jprb, 0.620848e-03_jprb, 0.935561e-03_jprb, 0.171252e-02_jprb /)
forref(:, 6) = (/ 0.212267e-02_jprb, 0.218564e-02_jprb, 0.222227e-02_jprb, 0.199650e-02_jprb /)
forref(:, 7) = (/ 0.291120e-02_jprb, 0.281168e-02_jprb, 0.259543e-02_jprb, 0.210159e-02_jprb /)
forref(:, 8) = (/ 0.316249e-02_jprb, 0.310695e-02_jprb, 0.279501e-02_jprb, 0.208076e-02_jprb /)
forref(:, 9) = (/ 0.354993e-02_jprb, 0.336989e-02_jprb, 0.298930e-02_jprb, 0.180424e-02_jprb /)
forref(:,10) = (/ 0.397729e-02_jprb, 0.367409e-02_jprb, 0.328982e-02_jprb, 0.177807e-02_jprb /)
forref(:,11) = (/ 0.408831e-02_jprb, 0.398792e-02_jprb, 0.352727e-02_jprb, 0.192470e-02_jprb /)
forref(:,12) = (/ 0.433926e-02_jprb, 0.420667e-02_jprb, 0.383894e-02_jprb, 0.220836e-02_jprb /)
forref(:,13) = (/ 0.436397e-02_jprb, 0.433769e-02_jprb, 0.425752e-02_jprb, 0.237343e-02_jprb /)
forref(:,14) = (/ 0.440525e-02_jprb, 0.449018e-02_jprb, 0.451881e-02_jprb, 0.269169e-02_jprb /)
forref(:,15) = (/ 0.491350e-02_jprb, 0.481760e-02_jprb, 0.475799e-02_jprb, 0.362666e-02_jprb /)
forref(:,16) = (/ 0.561641e-02_jprb, 0.524553e-02_jprb, 0.512473e-02_jprb, 0.493802e-02_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

selfref(:, 1) = (/ &
 & 0.115887e-03_jprb, 0.926537e-04_jprb, 0.740783e-04_jprb, 0.592270e-04_jprb, 0.473530e-04_jprb, &
 & 0.378596e-04_jprb, 0.302694e-04_jprb, 0.242010e-04_jprb, 0.193491e-04_jprb, 0.154700e-04_jprb /)  
selfref(:, 2) = (/ &
 & 0.459557e-03_jprb, 0.381962e-03_jprb, 0.317469e-03_jprb, 0.263866e-03_jprb, 0.219313e-03_jprb, &
 & 0.182283e-03_jprb, 0.151505e-03_jprb, 0.125924e-03_jprb, 0.104662e-03_jprb, 0.869904e-04_jprb /)  
selfref(:, 3) = (/ &
 & 0.166821e-02_jprb, 0.151103e-02_jprb, 0.136866e-02_jprb, 0.123970e-02_jprb, 0.112290e-02_jprb, &
 & 0.101710e-02_jprb, 0.921266e-03_jprb, 0.834463e-03_jprb, 0.755839e-03_jprb, 0.684623e-03_jprb /)  
selfref(:, 4) = (/ &
 & 0.460175e-02_jprb, 0.421372e-02_jprb, 0.385842e-02_jprb, 0.353307e-02_jprb, 0.323516e-02_jprb, &
 & 0.296236e-02_jprb, 0.271257e-02_jprb, 0.248385e-02_jprb, 0.227440e-02_jprb, 0.208262e-02_jprb /)  
selfref(:, 5) = (/ &
 & 0.101589e-01_jprb, 0.924742e-02_jprb, 0.841772e-02_jprb, 0.766247e-02_jprb, 0.697497e-02_jprb, &
 & 0.634917e-02_jprb, 0.577951e-02_jprb, 0.526096e-02_jprb, 0.478893e-02_jprb, 0.435926e-02_jprb /)  
selfref(:, 6) = (/ &
 & 0.328043e-01_jprb, 0.300853e-01_jprb, 0.275917e-01_jprb, 0.253048e-01_jprb, 0.232075e-01_jprb, &
 & 0.212839e-01_jprb, 0.195198e-01_jprb, 0.179020e-01_jprb, 0.164182e-01_jprb, 0.150574e-01_jprb /)  
selfref(:, 7) = (/ &
 & 0.405936e-01_jprb, 0.376032e-01_jprb, 0.348331e-01_jprb, 0.322671e-01_jprb, 0.298901e-01_jprb, &
 & 0.276883e-01_jprb, 0.256486e-01_jprb, 0.237591e-01_jprb, 0.220089e-01_jprb, 0.203876e-01_jprb /)  
selfref(:, 8) = (/ &
 & 0.448362e-01_jprb, 0.413811e-01_jprb, 0.381923e-01_jprb, 0.352492e-01_jprb, 0.325329e-01_jprb, &
 & 0.300259e-01_jprb, 0.277121e-01_jprb, 0.255766e-01_jprb, 0.236056e-01_jprb, 0.217866e-01_jprb /)  
selfref(:, 9) = (/ &
 & 0.479741e-01_jprb, 0.445389e-01_jprb, 0.413497e-01_jprb, 0.383889e-01_jprb, 0.356400e-01_jprb, &
 & 0.330880e-01_jprb, 0.307188e-01_jprb, 0.285191e-01_jprb, 0.264770e-01_jprb, 0.245812e-01_jprb /)  
selfref(:,10) = (/ &
 & 0.519308e-01_jprb, 0.484130e-01_jprb, 0.451335e-01_jprb, 0.420761e-01_jprb, 0.392259e-01_jprb, &
 & 0.365687e-01_jprb, 0.340916e-01_jprb, 0.317822e-01_jprb, 0.296293e-01_jprb, 0.276222e-01_jprb /)  
selfref(:,11) = (/ &
 & 0.572039e-01_jprb, 0.527780e-01_jprb, 0.486945e-01_jprb, 0.449270e-01_jprb, 0.414510e-01_jprb, &
 & 0.382439e-01_jprb, 0.352849e-01_jprb, 0.325549e-01_jprb, 0.300361e-01_jprb, 0.277122e-01_jprb /)  
selfref(:,12) = (/ &
 & 0.601046e-01_jprb, 0.554411e-01_jprb, 0.511395e-01_jprb, 0.471716e-01_jprb, 0.435116e-01_jprb, &
 & 0.401356e-01_jprb, 0.370215e-01_jprb, 0.341490e-01_jprb, 0.314994e-01_jprb, 0.290554e-01_jprb /)  
selfref(:,13) = (/ &
 & 0.616595e-01_jprb, 0.567145e-01_jprb, 0.521662e-01_jprb, 0.479826e-01_jprb, 0.441346e-01_jprb, &
 & 0.405951e-01_jprb, 0.373395e-01_jprb, 0.343450e-01_jprb, 0.315906e-01_jprb, 0.290571e-01_jprb /)  
selfref(:,14) = (/ &
 & 0.647916e-01_jprb, 0.592493e-01_jprb, 0.541811e-01_jprb, 0.495465e-01_jprb, 0.453083e-01_jprb, &
 & 0.414326e-01_jprb, 0.378885e-01_jprb, 0.346475e-01_jprb, 0.316837e-01_jprb, 0.289735e-01_jprb /)  
selfref(:,15) = (/ &
 & 0.694231e-01_jprb, 0.637703e-01_jprb, 0.585777e-01_jprb, 0.538079e-01_jprb, 0.494265e-01_jprb, &
 & 0.454019e-01_jprb, 0.417050e-01_jprb, 0.383091e-01_jprb, 0.351897e-01_jprb, 0.323244e-01_jprb /)  
selfref(:,16) = (/ &
 & 0.761764e-01_jprb, 0.701815e-01_jprb, 0.646584e-01_jprb, 0.595700e-01_jprb, 0.548820e-01_jprb, &
 & 0.505629e-01_jprb, 0.465838e-01_jprb, 0.429178e-01_jprb, 0.395403e-01_jprb, 0.364286e-01_jprb /)  

if (lhook) call dr_hook('srtm_kgb21',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb21:error reading file radsrtm")

end subroutine srtm_kgb21
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

