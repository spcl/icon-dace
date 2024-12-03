! # 1 "ifsrrtm/srtm_kgb20.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb20.f90"
subroutine srtm_kgb20

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
!     reformatted for f90 by jjmorcrette, ecmwf
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
use yoesrta20 , only : ka, kb, selfref, forref, sfluxref, rayl, absch4, layreffr  ,&
  & ka_d, kb_d

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
! # 29 "ifsrrtm/srtm_kgb20.f90" 2

if (lhook) call dr_hook('srtm_kgb20',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb20:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb20:')
endif

sfluxref = (/ &
 & 9.34081_jprb , 8.93720_jprb    , 8.19346_jprb    , 7.39196_jprb    , &
 & 6.12127_jprb , 5.23956_jprb    , 4.24941_jprb    , 3.20013_jprb    , &
 & 2.16047_jprb , 0.234509_jprb   , 0.194593_jprb   , 0.151512_jprb   , &
 & 0.110315_jprb, 7.09959e-02_jprb, 2.70573e-02_jprb, 3.36042e-03_jprb /)  
  
absch4 = (/   &
 & 1.01381e-03_jprb,6.33692e-03_jprb,1.94185e-02_jprb,4.83210e-02_jprb, &
 & 2.36574e-03_jprb,6.61973e-04_jprb,5.64552e-04_jprb,2.83183e-04_jprb, &
 & 7.43623e-05_jprb,8.90159e-07_jprb,6.98728e-07_jprb,6.51832e-08_jprb, &
 & 2.96619e-08_jprb,         0._jprb,         0._jprb,         0._jprb /)  

!     rayleigh extinction coefficient at v = 5670 cm-1.
rayl = 4.12e-09_jprb

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


forref(:, 1) = (/ 0.214504e-06_jprb, 0.460418e-06_jprb, 0.357608e-05_jprb, 0.192037e-05_jprb /)
forref(:, 2) = (/ 0.142576e-05_jprb, 0.364463e-05_jprb, 0.117033e-04_jprb, 0.112085e-04_jprb /)
forref(:, 3) = (/ 0.101536e-04_jprb, 0.124096e-04_jprb, 0.509190e-04_jprb, 0.565282e-04_jprb /)
forref(:, 4) = (/ 0.143394e-03_jprb, 0.154700e-03_jprb, 0.466498e-03_jprb, 0.918829e-03_jprb /)
forref(:, 5) = (/ 0.251631e-02_jprb, 0.241729e-02_jprb, 0.240057e-02_jprb, 0.350408e-02_jprb /)
forref(:, 6) = (/ 0.410309e-02_jprb, 0.416851e-02_jprb, 0.390925e-02_jprb, 0.383694e-02_jprb /)
forref(:, 7) = (/ 0.445387e-02_jprb, 0.448657e-02_jprb, 0.432310e-02_jprb, 0.370739e-02_jprb /)
forref(:, 8) = (/ 0.458150e-02_jprb, 0.460014e-02_jprb, 0.450245e-02_jprb, 0.336718e-02_jprb /)
forref(:, 9) = (/ 0.465423e-02_jprb, 0.465595e-02_jprb, 0.467006e-02_jprb, 0.368061e-02_jprb /)
forref(:,10) = (/ 0.493955e-02_jprb, 0.490181e-02_jprb, 0.481941e-02_jprb, 0.367577e-02_jprb /)
forref(:,11) = (/ 0.511876e-02_jprb, 0.490981e-02_jprb, 0.493303e-02_jprb, 0.357423e-02_jprb /)
forref(:,12) = (/ 0.509845e-02_jprb, 0.511556e-02_jprb, 0.504031e-02_jprb, 0.355915e-02_jprb /)
forref(:,13) = (/ 0.523822e-02_jprb, 0.530473e-02_jprb, 0.523811e-02_jprb, 0.414259e-02_jprb /)
forref(:,14) = (/ 0.551133e-02_jprb, 0.535831e-02_jprb, 0.546702e-02_jprb, 0.473875e-02_jprb /)
forref(:,15) = (/ 0.609781e-02_jprb, 0.589859e-02_jprb, 0.561187e-02_jprb, 0.528981e-02_jprb /)
forref(:,16) = (/ 0.644958e-02_jprb, 0.631718e-02_jprb, 0.625201e-02_jprb, 0.600448e-02_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

selfref(:, 1) = (/ &
 & 0.217058e-03_jprb, 0.176391e-03_jprb, 0.143342e-03_jprb, 0.116486e-03_jprb, 0.946614e-04_jprb, &
 & 0.769257e-04_jprb, 0.625131e-04_jprb, 0.508007e-04_jprb, 0.412828e-04_jprb, 0.335481e-04_jprb /)  
selfref(:, 2) = (/ &
 & 0.598055e-03_jprb, 0.484805e-03_jprb, 0.393000e-03_jprb, 0.318580e-03_jprb, 0.258252e-03_jprb, &
 & 0.209348e-03_jprb, 0.169705e-03_jprb, 0.137569e-03_jprb, 0.111518e-03_jprb, 0.904008e-04_jprb /)  
selfref(:, 3) = (/ &
 & 0.102691e-02_jprb, 0.930281e-03_jprb, 0.842740e-03_jprb, 0.763437e-03_jprb, 0.691596e-03_jprb, &
 & 0.626516e-03_jprb, 0.567560e-03_jprb, 0.514152e-03_jprb, 0.465769e-03_jprb, 0.421940e-03_jprb /)  
selfref(:, 4) = (/ &
 & 0.388569e-02_jprb, 0.365098e-02_jprb, 0.343045e-02_jprb, 0.322324e-02_jprb, 0.302854e-02_jprb, &
 & 0.284561e-02_jprb, 0.267372e-02_jprb, 0.251222e-02_jprb, 0.236047e-02_jprb, 0.221789e-02_jprb /)  
selfref(:, 5) = (/ &
 & 0.349845e-01_jprb, 0.326678e-01_jprb, 0.305045e-01_jprb, 0.284845e-01_jprb, 0.265982e-01_jprb, &
 & 0.248369e-01_jprb, 0.231921e-01_jprb, 0.216563e-01_jprb, 0.202222e-01_jprb, 0.188831e-01_jprb /)  
selfref(:, 6) = (/ &
 & 0.613705e-01_jprb, 0.562676e-01_jprb, 0.515890e-01_jprb, 0.472994e-01_jprb, 0.433665e-01_jprb, &
 & 0.397606e-01_jprb, 0.364545e-01_jprb, 0.334233e-01_jprb, 0.306442e-01_jprb, 0.280961e-01_jprb /)  
selfref(:, 7) = (/ &
 & 0.656981e-01_jprb, 0.602660e-01_jprb, 0.552830e-01_jprb, 0.507120e-01_jprb, 0.465190e-01_jprb, &
 & 0.426726e-01_jprb, 0.391443e-01_jprb, 0.359077e-01_jprb, 0.329387e-01_jprb, 0.302153e-01_jprb /)  
selfref(:, 8) = (/ &
 & 0.671782e-01_jprb, 0.616461e-01_jprb, 0.565695e-01_jprb, 0.519110e-01_jprb, 0.476361e-01_jprb, &
 & 0.437132e-01_jprb, 0.401134e-01_jprb, 0.368100e-01_jprb, 0.337787e-01_jprb, 0.309970e-01_jprb /)  
selfref(:, 9) = (/ &
 & 0.675902e-01_jprb, 0.620888e-01_jprb, 0.570351e-01_jprb, 0.523928e-01_jprb, 0.481284e-01_jprb, &
 & 0.442110e-01_jprb, 0.406125e-01_jprb, 0.373069e-01_jprb, 0.342703e-01_jprb, 0.314809e-01_jprb /)  
selfref(:,10) = (/ &
 & 0.708308e-01_jprb, 0.651419e-01_jprb, 0.599099e-01_jprb, 0.550981e-01_jprb, 0.506728e-01_jprb, &
 & 0.466030e-01_jprb, 0.428600e-01_jprb, 0.394176e-01_jprb, 0.362517e-01_jprb, 0.333401e-01_jprb /)  
selfref(:,11) = (/ &
 & 0.698445e-01_jprb, 0.646584e-01_jprb, 0.598573e-01_jprb, 0.554128e-01_jprb, 0.512982e-01_jprb, &
 & 0.474892e-01_jprb, 0.439630e-01_jprb, 0.406986e-01_jprb, 0.376766e-01_jprb, 0.348791e-01_jprb /)  
selfref(:,12) = (/ &
 & 0.743921e-01_jprb, 0.682057e-01_jprb, 0.625337e-01_jprb, 0.573334e-01_jprb, 0.525655e-01_jprb, &
 & 0.481942e-01_jprb, 0.441863e-01_jprb, 0.405118e-01_jprb, 0.371428e-01_jprb, 0.340540e-01_jprb /)  
selfref(:,13) = (/ &
 & 0.775758e-01_jprb, 0.709818e-01_jprb, 0.649484e-01_jprb, 0.594277e-01_jprb, 0.543764e-01_jprb, &
 & 0.497544e-01_jprb, 0.455253e-01_jprb, 0.416556e-01_jprb, 0.381149e-01_jprb, 0.348751e-01_jprb /)  
selfref(:,14) = (/ &
 & 0.776545e-01_jprb, 0.714761e-01_jprb, 0.657894e-01_jprb, 0.605550e-01_jprb, 0.557372e-01_jprb, &
 & 0.513026e-01_jprb, 0.472209e-01_jprb, 0.434639e-01_jprb, 0.400058e-01_jprb, 0.368229e-01_jprb /)  
selfref(:,15) = (/ &
 & 0.855675e-01_jprb, 0.787337e-01_jprb, 0.724456e-01_jprb, 0.666598e-01_jprb, 0.613360e-01_jprb, &
 & 0.564374e-01_jprb, 0.519301e-01_jprb, 0.477827e-01_jprb, 0.439666e-01_jprb, 0.404552e-01_jprb /)  
selfref(:,16) = (/ &
 & 0.934781e-01_jprb, 0.855190e-01_jprb, 0.782376e-01_jprb, 0.715761e-01_jprb, 0.654819e-01_jprb, &
 & 0.599065e-01_jprb, 0.548058e-01_jprb, 0.501394e-01_jprb, 0.458704e-01_jprb, 0.419648e-01_jprb /)  
     
!     -----------------------------------------------------------------

if (lhook) call dr_hook('srtm_kgb20',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb20:error reading file radsrtm")

end subroutine srtm_kgb20
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

