! # 1 "ifsrrtm/srtm_kgb18.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb18.f90"
subroutine srtm_kgb18

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
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
use yoesrta18 , only : ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr, &
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
! # 30 "ifsrrtm/srtm_kgb18.f90" 2

if (lhook) call dr_hook('srtm_kgb18',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb18:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb18:')
endif

sfluxref(:,1) = (/ &
 & 3.65840_jprb    , 3.54375_jprb    , 3.34481_jprb    , 3.10534_jprb    , &
 & 2.79879_jprb    , 2.42841_jprb    , 1.98748_jprb    , 1.49377_jprb    , &
 & 1.00196_jprb    , 0.108342_jprb   , 8.95099e-02_jprb, 7.05199e-02_jprb, &
 & 5.16432e-02_jprb, 3.27635e-02_jprb, 1.25133e-02_jprb, 1.73001e-03_jprb /)    
sfluxref(:,2) = (/ &
 & 3.86372_jprb    , 3.48521_jprb    , 3.30790_jprb    , 3.08103_jprb    , &
 & 2.77552_jprb    , 2.40722_jprb    , 1.97307_jprb    , 1.48023_jprb    , &
 & 0.993055_jprb   , 0.107691_jprb   , 8.84430e-02_jprb, 6.99354e-02_jprb, &
 & 5.07881e-02_jprb, 3.24121e-02_jprb, 1.19442e-02_jprb, 1.57612e-03_jprb /)  
sfluxref(:,3) = (/ &
 & 3.90370_jprb    , 3.50657_jprb    , 3.30629_jprb    , 3.06046_jprb    , &
 & 2.76982_jprb    , 2.39907_jprb    , 1.96358_jprb    , 1.47458_jprb    , &
 & 0.988475_jprb   , 0.106698_jprb   , 8.75242e-02_jprb, 6.85898e-02_jprb, &
 & 5.04798e-02_jprb, 3.13718e-02_jprb, 1.09533e-02_jprb, 1.57612e-03_jprb /)  
sfluxref(:,4) = (/ &
 & 3.93165_jprb    , 3.52058_jprb    , 3.31346_jprb    , 3.04944_jprb    , &
 & 2.76074_jprb    , 2.39433_jprb    , 1.95556_jprb    , 1.46712_jprb    , &
 & 0.984056_jprb   , 0.105885_jprb   , 8.73062e-02_jprb, 6.84054e-02_jprb, &
 & 4.87443e-02_jprb, 2.99295e-02_jprb, 1.09533e-02_jprb, 1.57612e-03_jprb /)  
sfluxref(:,5) = (/ &
 & 3.94082_jprb    , 3.55221_jprb    , 3.31863_jprb    , 3.04730_jprb    , &
 & 2.74918_jprb    , 2.38328_jprb    , 1.95212_jprb    , 1.45889_jprb    , &
 & 0.978888_jprb   , 0.105102_jprb   , 8.65732e-02_jprb, 6.74563e-02_jprb, &
 & 4.76592e-02_jprb, 2.91017e-02_jprb, 1.09533e-02_jprb, 1.57612e-03_jprb /)  
sfluxref(:,6) = (/ &
 & 3.94198_jprb    , 3.58743_jprb    , 3.32106_jprb    , 3.05866_jprb    , &
 & 2.74115_jprb    , 2.36939_jprb    , 1.94305_jprb    , 1.45180_jprb    , &
 & 0.971784_jprb   , 1.04045e-01_jprb, 8.53731e-02_jprb, 6.60654e-02_jprb, &
 & 4.63228e-02_jprb, 2.91016e-02_jprb, 1.09552e-02_jprb, 1.57612e-03_jprb /)  
sfluxref(:,7) = (/ &
 & 3.93596_jprb    , 3.63366_jprb    , 3.33144_jprb    , 3.06252_jprb    , &
 & 2.74054_jprb    , 2.35492_jprb    , 1.92769_jprb    , 1.44300_jprb    , &
 & 0.961809_jprb   , 1.02867e-01_jprb, 8.34164e-02_jprb, 6.41005e-02_jprb, &
 & 4.61826e-02_jprb, 2.91006e-02_jprb, 1.09553e-02_jprb, 1.57612e-03_jprb /)  
sfluxref(:,8) = (/ &
 & 3.92520_jprb    , 3.69078_jprb    , 3.35656_jprb    , 3.07055_jprb    , &
 & 2.73862_jprb    , 2.34430_jprb    , 1.90187_jprb    , 1.42242_jprb    , &
 & 0.946676_jprb   , 9.96302e-02_jprb, 8.14421e-02_jprb, 6.38622e-02_jprb, &
 & 4.61794e-02_jprb, 2.91017e-02_jprb, 1.09553e-02_jprb, 1.57612e-03_jprb /)  
sfluxref(:,9) = (/ &
 & 3.80721_jprb    , 3.74437_jprb    , 3.50205_jprb    , 3.18009_jprb    , &
 & 2.75757_jprb    , 2.29188_jprb    , 1.84382_jprb    , 1.35694_jprb    , &
 & 0.914040_jprb   , 9.86811e-02_jprb, 8.14321e-02_jprb, 6.38541e-02_jprb, &
 & 4.61795e-02_jprb, 2.90960e-02_jprb, 1.09613e-02_jprb, 1.57612e-03_jprb /)  

!     rayleigh extinction coefficient at v = 4325 cm-1.
rayl = 1.39e-09_jprb

strrat = 38.9589_jprb

layreffr = 6

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

  
forref(:, 1) = (/ 0.860560e-06_jprb, 0.130439e-05_jprb, 0.382378e-05_jprb /)
forref(:, 2) = (/ 0.817926e-06_jprb, 0.158599e-05_jprb, 0.658771e-04_jprb /)
forref(:, 3) = (/ 0.129369e-05_jprb, 0.824406e-05_jprb, 0.952778e-04_jprb /)
forref(:, 4) = (/ 0.438918e-05_jprb, 0.375356e-04_jprb, 0.119111e-03_jprb /)
forref(:, 5) = (/ 0.306057e-04_jprb, 0.622798e-04_jprb, 0.100740e-03_jprb /)
forref(:, 6) = (/ 0.891934e-04_jprb, 0.856393e-04_jprb, 0.635583e-04_jprb /)
forref(:, 7) = (/ 0.171959e-03_jprb, 0.173431e-03_jprb, 0.611721e-04_jprb /)
forref(:, 8) = (/ 0.357795e-03_jprb, 0.247261e-03_jprb, 0.488864e-04_jprb /)
forref(:, 9) = (/ 0.326623e-03_jprb, 0.289471e-03_jprb, 0.548834e-04_jprb /)
forref(:,10) = (/ 0.345103e-03_jprb, 0.320898e-03_jprb, 0.633214e-04_jprb /)
forref(:,11) = (/ 0.392567e-03_jprb, 0.325153e-03_jprb, 0.744479e-04_jprb /)
forref(:,12) = (/ 0.349277e-03_jprb, 0.345610e-03_jprb, 0.916479e-04_jprb /)
forref(:,13) = (/ 0.425161e-03_jprb, 0.348452e-03_jprb, 0.125788e-03_jprb /)
forref(:,14) = (/ 0.407594e-03_jprb, 0.435836e-03_jprb, 0.287583e-03_jprb /)
forref(:,15) = (/ 0.521605e-03_jprb, 0.486596e-03_jprb, 0.483511e-03_jprb /)
forref(:,16) = (/ 0.773790e-03_jprb, 0.737247e-03_jprb, 0.665939e-03_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).
     
selfref(:, 1) = (/ &
 & 0.750370e-03_jprb, 0.644938e-03_jprb, 0.554321e-03_jprb, 0.476436e-03_jprb, 0.409494e-03_jprb, &
 & 0.351957e-03_jprb, 0.302505e-03_jprb, 0.260002e-03_jprb, 0.223470e-03_jprb, 0.192071e-03_jprb /)  
selfref(:, 2) = (/ &
 & 0.136135e-02_jprb, 0.113187e-02_jprb, 0.941076e-03_jprb, 0.782440e-03_jprb, 0.650546e-03_jprb, &
 & 0.540885e-03_jprb, 0.449709e-03_jprb, 0.373902e-03_jprb, 0.310874e-03_jprb, 0.258471e-03_jprb /)  
selfref(:, 3) = (/ &
 & 0.333950e-02_jprb, 0.256391e-02_jprb, 0.196845e-02_jprb, 0.151129e-02_jprb, 0.116030e-02_jprb, &
 & 0.890824e-03_jprb, 0.683934e-03_jprb, 0.525093e-03_jprb, 0.403143e-03_jprb, 0.309515e-03_jprb /)  
selfref(:, 4) = (/ &
 & 0.793392e-02_jprb, 0.589865e-02_jprb, 0.438548e-02_jprb, 0.326048e-02_jprb, 0.242408e-02_jprb, &
 & 0.180223e-02_jprb, 0.133991e-02_jprb, 0.996186e-03_jprb, 0.740636e-03_jprb, 0.550642e-03_jprb /)  
selfref(:, 5) = (/ &
 & 0.828169e-02_jprb, 0.703139e-02_jprb, 0.596984e-02_jprb, 0.506856e-02_jprb, 0.430335e-02_jprb, &
 & 0.365366e-02_jprb, 0.310206e-02_jprb, 0.263374e-02_jprb, 0.223612e-02_jprb, 0.189852e-02_jprb /)  
selfref(:, 6) = (/ &
 & 0.834190e-02_jprb, 0.780225e-02_jprb, 0.729750e-02_jprb, 0.682541e-02_jprb, 0.638386e-02_jprb, &
 & 0.597087e-02_jprb, 0.558460e-02_jprb, 0.522332e-02_jprb, 0.488541e-02_jprb, 0.456936e-02_jprb /)  
selfref(:, 7) = (/ &
 & 0.119082e-01_jprb, 0.112566e-01_jprb, 0.106406e-01_jprb, 0.100583e-01_jprb, 0.950785e-02_jprb, &
 & 0.898755e-02_jprb, 0.849571e-02_jprb, 0.803080e-02_jprb, 0.759132e-02_jprb, 0.717590e-02_jprb /)  
selfref(:, 8) = (/ &
 & 0.144004e-01_jprb, 0.141762e-01_jprb, 0.139554e-01_jprb, 0.137381e-01_jprb, 0.135241e-01_jprb, &
 & 0.133135e-01_jprb, 0.131062e-01_jprb, 0.129021e-01_jprb, 0.127011e-01_jprb, 0.125033e-01_jprb /)  
selfref(:, 9) = (/ &
 & 0.186171e-01_jprb, 0.175281e-01_jprb, 0.165027e-01_jprb, 0.155373e-01_jprb, 0.146284e-01_jprb, &
 & 0.137726e-01_jprb, 0.129670e-01_jprb, 0.122084e-01_jprb, 0.114942e-01_jprb, 0.108218e-01_jprb /)  
selfref(:,10) = (/ &
 & 0.209396e-01_jprb, 0.195077e-01_jprb, 0.181737e-01_jprb, 0.169309e-01_jprb, 0.157731e-01_jprb, &
 & 0.146945e-01_jprb, 0.136897e-01_jprb, 0.127535e-01_jprb, 0.118814e-01_jprb, 0.110689e-01_jprb /)  
selfref(:,11) = (/ &
 & 0.203661e-01_jprb, 0.193311e-01_jprb, 0.183487e-01_jprb, 0.174163e-01_jprb, 0.165312e-01_jprb, &
 & 0.156911e-01_jprb, 0.148937e-01_jprb, 0.141368e-01_jprb, 0.134184e-01_jprb, 0.127365e-01_jprb /)  
selfref(:,12) = (/ &
 & 0.226784e-01_jprb, 0.210210e-01_jprb, 0.194848e-01_jprb, 0.180608e-01_jprb, 0.167409e-01_jprb, &
 & 0.155174e-01_jprb, 0.143834e-01_jprb, 0.133322e-01_jprb, 0.123579e-01_jprb, 0.114547e-01_jprb /)  
selfref(:,13) = (/ &
 & 0.221773e-01_jprb, 0.210306e-01_jprb, 0.199433e-01_jprb, 0.189122e-01_jprb, 0.179344e-01_jprb, &
 & 0.170071e-01_jprb, 0.161278e-01_jprb, 0.152939e-01_jprb, 0.145032e-01_jprb, 0.137533e-01_jprb /)  
selfref(:,14) = (/ &
 & 0.275920e-01_jprb, 0.252595e-01_jprb, 0.231241e-01_jprb, 0.211693e-01_jprb, 0.193797e-01_jprb, &
 & 0.177415e-01_jprb, 0.162417e-01_jprb, 0.148687e-01_jprb, 0.136117e-01_jprb, 0.124610e-01_jprb /)  
selfref(:,15) = (/ &
 & 0.288687e-01_jprb, 0.269968e-01_jprb, 0.252462e-01_jprb, 0.236092e-01_jprb, 0.220783e-01_jprb, &
 & 0.206466e-01_jprb, 0.193078e-01_jprb, 0.180559e-01_jprb, 0.168851e-01_jprb, 0.157902e-01_jprb /)  
selfref(:,16) = (/ &
 & 0.371842e-01_jprb, 0.347595e-01_jprb, 0.324929e-01_jprb, 0.303741e-01_jprb, 0.283934e-01_jprb, &
 & 0.265419e-01_jprb, 0.248112e-01_jprb, 0.231933e-01_jprb, 0.216809e-01_jprb, 0.202671e-01_jprb /)  

if (lhook) call dr_hook('srtm_kgb18',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb18:error reading file radsrtm")

end subroutine srtm_kgb18
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

