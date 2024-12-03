! # 1 "ifsrrtm/rrtm_kgb12.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb12.f90"
subroutine rrtm_kgb12

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use yommp0    , only : nproc, myproc
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto12, only : kao, kao_d, selfrefo, forrefo, fracrefao

!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle



! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 28 "ifsrrtm/rrtm_kgb12.f90" 2

if (lhook) call dr_hook('rrtm_kgb12',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d
  kao = real(kao_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb12:')
endif

! planck fraction mapping level : p = 174.1640 mbar, t= 215.78 k
      fracrefao(:, 1) = (/ &
     &  1.3984e-01_jprb,1.6809e-01_jprb,1.8072e-01_jprb,1.5400e-01_jprb,1.2613e-01_jprb,9.6959e-02_jprb, &
     &  5.9713e-02_jprb,3.8631e-02_jprb,2.6937e-02_jprb,3.1711e-03_jprb,2.3458e-03_jprb,1.4653e-03_jprb, &
     &  1.0567e-03_jprb,6.6504e-04_jprb,2.4957e-04_jprb,3.5172e-05_jprb/)
      fracrefao(:, 2) = (/ &
     &  1.2745e-01_jprb,1.6107e-01_jprb,1.6568e-01_jprb,1.5436e-01_jprb,1.3183e-01_jprb,1.0166e-01_jprb, &
     &  6.4506e-02_jprb,4.7756e-02_jprb,3.4472e-02_jprb,3.7189e-03_jprb,2.9349e-03_jprb,2.1469e-03_jprb, &
     &  1.3746e-03_jprb,7.1691e-04_jprb,2.8057e-04_jprb,5.6242e-05_jprb/)
      fracrefao(:, 3) = (/ &
     &  1.2181e-01_jprb,1.5404e-01_jprb,1.6540e-01_jprb,1.5255e-01_jprb,1.3736e-01_jprb,9.8856e-02_jprb, &
     &  6.8927e-02_jprb,5.1385e-02_jprb,3.7046e-02_jprb,4.0302e-03_jprb,3.0949e-03_jprb,2.3772e-03_jprb, &
     &  1.6538e-03_jprb,8.9641e-04_jprb,4.6991e-04_jprb,1.1251e-04_jprb/)
      fracrefao(:, 4) = (/ &
     &  1.1794e-01_jprb,1.4864e-01_jprb,1.6316e-01_jprb,1.5341e-01_jprb,1.3986e-01_jprb,9.6656e-02_jprb, &
     &  7.2478e-02_jprb,5.5061e-02_jprb,3.8886e-02_jprb,4.3398e-03_jprb,3.3576e-03_jprb,2.4891e-03_jprb, &
     &  1.7674e-03_jprb,1.0764e-03_jprb,7.7689e-04_jprb,1.1251e-04_jprb/)
      fracrefao(:, 5) = (/ &
     &  1.1635e-01_jprb,1.4342e-01_jprb,1.5924e-01_jprb,1.5670e-01_jprb,1.3740e-01_jprb,9.7087e-02_jprb, &
     &  7.6250e-02_jprb,5.7802e-02_jprb,4.0808e-02_jprb,4.4113e-03_jprb,3.6035e-03_jprb,2.6269e-03_jprb, &
     &  1.7586e-03_jprb,1.6498e-03_jprb,7.7689e-04_jprb,1.1251e-04_jprb/)
      fracrefao(:, 6) = (/ &
     &  1.1497e-01_jprb,1.3751e-01_jprb,1.5587e-01_jprb,1.5904e-01_jprb,1.3140e-01_jprb,1.0159e-01_jprb, &
     &  7.9729e-02_jprb,6.1475e-02_jprb,4.2382e-02_jprb,4.5291e-03_jprb,3.8161e-03_jprb,2.7683e-03_jprb, &
     &  1.9899e-03_jprb,2.0395e-03_jprb,7.7720e-04_jprb,1.1251e-04_jprb/)
      fracrefao(:, 7) = (/ &
     &  1.1331e-01_jprb,1.3015e-01_jprb,1.5574e-01_jprb,1.5489e-01_jprb,1.2697e-01_jprb,1.0746e-01_jprb, &
     &  8.4777e-02_jprb,6.5145e-02_jprb,4.4293e-02_jprb,4.7426e-03_jprb,3.8383e-03_jprb,2.9065e-03_jprb, &
     &  2.8430e-03_jprb,2.0401e-03_jprb,7.7689e-04_jprb,1.1251e-04_jprb/)
      fracrefao(:, 8) = (/ &
     &  1.0993e-01_jprb,1.2320e-01_jprb,1.4893e-01_jprb,1.4573e-01_jprb,1.3174e-01_jprb,1.1149e-01_jprb, &
     &  9.3326e-02_jprb,6.9942e-02_jprb,4.6762e-02_jprb,4.9309e-03_jprb,3.8583e-03_jprb,4.1889e-03_jprb, &
     &  3.0415e-03_jprb,2.0406e-03_jprb,7.7720e-04_jprb,1.1251e-04_jprb/)
      fracrefao(:, 9) = (/ &
     &  1.2028e-01_jprb,1.2091e-01_jprb,1.3098e-01_jprb,1.3442e-01_jprb,1.3574e-01_jprb,1.1739e-01_jprb, &
     &  9.5343e-02_jprb,7.0224e-02_jprb,5.3456e-02_jprb,6.0206e-03_jprb,5.0758e-03_jprb,4.1906e-03_jprb, &
     &  3.0431e-03_jprb,2.0400e-03_jprb,7.7689e-04_jprb,1.1251e-04_jprb/)


!     ------------------------------------------------------------------

!     the array kao contains absorption coefs at the 16 chosen g-values 
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


!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.4739e-04_jprb,3.1686e-04_jprb,8.5973e-04_jprb,1.9039e-03_jprb,3.1820e-03_jprb,3.6596e-03_jprb, &
     &3.8724e-03_jprb,3.6785e-03_jprb,3.7141e-03_jprb,3.7646e-03_jprb,4.2955e-03_jprb,4.6343e-03_jprb, &
     &5.0612e-03_jprb,4.0227e-03_jprb,4.2966e-03_jprb,4.6622e-03_jprb/)
      forrefo(2,:) = (/ &
     &1.9397e-04_jprb,3.6322e-04_jprb,8.9797e-04_jprb,2.1001e-03_jprb,3.0307e-03_jprb,3.5563e-03_jprb, &
     &3.8498e-03_jprb,3.5741e-03_jprb,3.5914e-03_jprb,3.7658e-03_jprb,3.8895e-03_jprb,4.4072e-03_jprb, &
     &4.7112e-03_jprb,4.2230e-03_jprb,4.2666e-03_jprb,4.6634e-03_jprb/)
      forrefo(3,:) = (/ &
     &3.1506e-04_jprb,7.3687e-04_jprb,1.9678e-03_jprb,2.5531e-03_jprb,2.8345e-03_jprb,2.7809e-03_jprb, &
     &2.9124e-03_jprb,2.7125e-03_jprb,2.6644e-03_jprb,2.4907e-03_jprb,2.7032e-03_jprb,4.0967e-03_jprb, &
     &4.1971e-03_jprb,4.4507e-03_jprb,4.2293e-03_jprb,4.6633e-03_jprb/)
      forrefo(4,:) = (/ &
     &8.8196e-04_jprb,2.1125e-03_jprb,2.8042e-03_jprb,2.8891e-03_jprb,2.4362e-03_jprb,1.8733e-03_jprb, &
     &1.4078e-03_jprb,1.1987e-03_jprb,1.2808e-03_jprb,8.9050e-04_jprb,9.4375e-04_jprb,7.8351e-04_jprb, &
     &1.0756e-03_jprb,1.6586e-03_jprb,1.7511e-03_jprb,4.7803e-03_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.37879e-02_jprb, 2.10719e-02_jprb, 1.86660e-02_jprb, 1.65348e-02_jprb, 1.46469e-02_jprb, &
     & 1.29746e-02_jprb, 1.14932e-02_jprb, 1.01810e-02_jprb, 9.01858e-03_jprb, 7.98888e-03_jprb/)
      selfrefo(:, 2) = (/ &
     & 3.10625e-02_jprb, 2.82664e-02_jprb, 2.57220e-02_jprb, 2.34066e-02_jprb, 2.12997e-02_jprb, &
     & 1.93824e-02_jprb, 1.76377e-02_jprb, 1.60500e-02_jprb, 1.46053e-02_jprb, 1.32906e-02_jprb/)
      selfrefo(:, 3) = (/ &
     & 5.19103e-02_jprb, 4.80004e-02_jprb, 4.43850e-02_jprb, 4.10419e-02_jprb, 3.79506e-02_jprb, &
     & 3.50922e-02_jprb, 3.24491e-02_jprb, 3.00050e-02_jprb, 2.77450e-02_jprb, 2.56553e-02_jprb/)
      selfrefo(:, 4) = (/ &
     & 9.12444e-02_jprb, 8.38675e-02_jprb, 7.70870e-02_jprb, 7.08547e-02_jprb, 6.51263e-02_jprb, &
     & 5.98610e-02_jprb, 5.50214e-02_jprb, 5.05730e-02_jprb, 4.64843e-02_jprb, 4.27262e-02_jprb/)
      selfrefo(:, 5) = (/ &
     & 1.11323e-01_jprb, 1.04217e-01_jprb, 9.75650e-02_jprb, 9.13376e-02_jprb, 8.55076e-02_jprb, &
     & 8.00498e-02_jprb, 7.49403e-02_jprb, 7.01570e-02_jprb, 6.56790e-02_jprb, 6.14868e-02_jprb/)
      selfrefo(:, 6) = (/ &
     & 1.25301e-01_jprb, 1.16877e-01_jprb, 1.09020e-01_jprb, 1.01691e-01_jprb, 9.48543e-02_jprb, &
     & 8.84774e-02_jprb, 8.25293e-02_jprb, 7.69810e-02_jprb, 7.18057e-02_jprb, 6.69784e-02_jprb/)
      selfrefo(:, 7) = (/ &
     & 1.34063e-01_jprb, 1.24662e-01_jprb, 1.15920e-01_jprb, 1.07791e-01_jprb, 1.00232e-01_jprb, &
     & 9.32035e-02_jprb, 8.66676e-02_jprb, 8.05900e-02_jprb, 7.49386e-02_jprb, 6.96836e-02_jprb/)
      selfrefo(:, 8) = (/ &
     & 1.26997e-01_jprb, 1.18306e-01_jprb, 1.10210e-01_jprb, 1.02668e-01_jprb, 9.56417e-02_jprb, &
     & 8.90964e-02_jprb, 8.29991e-02_jprb, 7.73190e-02_jprb, 7.20276e-02_jprb, 6.70984e-02_jprb/)
      selfrefo(:, 9) = (/ &
     & 1.28823e-01_jprb, 1.20235e-01_jprb, 1.12220e-01_jprb, 1.04739e-01_jprb, 9.77569e-02_jprb, &
     & 9.12402e-02_jprb, 8.51579e-02_jprb, 7.94810e-02_jprb, 7.41826e-02_jprb, 6.92374e-02_jprb/)
      selfrefo(:,10) = (/ &
     & 1.35802e-01_jprb, 1.25981e-01_jprb, 1.16870e-01_jprb, 1.08418e-01_jprb, 1.00577e-01_jprb, &
     & 9.33034e-02_jprb, 8.65557e-02_jprb, 8.02960e-02_jprb, 7.44890e-02_jprb, 6.91020e-02_jprb/)
      selfrefo(:,11) = (/ &
     & 1.35475e-01_jprb, 1.27572e-01_jprb, 1.20130e-01_jprb, 1.13122e-01_jprb, 1.06523e-01_jprb, &
     & 1.00309e-01_jprb, 9.44573e-02_jprb, 8.89470e-02_jprb, 8.37582e-02_jprb, 7.88721e-02_jprb/)
      selfrefo(:,12) = (/ &
     & 1.51195e-01_jprb, 1.41159e-01_jprb, 1.31790e-01_jprb, 1.23043e-01_jprb, 1.14876e-01_jprb, &
     & 1.07251e-01_jprb, 1.00132e-01_jprb, 9.34860e-02_jprb, 8.72809e-02_jprb, 8.14877e-02_jprb/)
      selfrefo(:,13) = (/ &
     & 1.57538e-01_jprb, 1.47974e-01_jprb, 1.38990e-01_jprb, 1.30552e-01_jprb, 1.22626e-01_jprb, &
     & 1.15181e-01_jprb, 1.08188e-01_jprb, 1.01620e-01_jprb, 9.54505e-02_jprb, 8.96556e-02_jprb/)
      selfrefo(:,14) = (/ &
     & 1.53567e-01_jprb, 1.41564e-01_jprb, 1.30500e-01_jprb, 1.20300e-01_jprb, 1.10898e-01_jprb, &
     & 1.02231e-01_jprb, 9.42406e-02_jprb, 8.68750e-02_jprb, 8.00851e-02_jprb, 7.38259e-02_jprb/)
      selfrefo(:,15) = (/ &
     & 1.53687e-01_jprb, 1.42981e-01_jprb, 1.33020e-01_jprb, 1.23753e-01_jprb, 1.15132e-01_jprb, &
     & 1.07112e-01_jprb, 9.96500e-02_jprb, 9.27080e-02_jprb, 8.62496e-02_jprb, 8.02412e-02_jprb/)
      selfrefo(:,16) = (/ &
     & 1.65129e-01_jprb, 1.53285e-01_jprb, 1.42290e-01_jprb, 1.32084e-01_jprb, 1.22610e-01_jprb, &
     & 1.13815e-01_jprb, 1.05651e-01_jprb, 9.80730e-02_jprb, 9.10384e-02_jprb, 8.45083e-02_jprb/)

if (lhook) call dr_hook('rrtm_kgb12',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb12:error reading file radrrtm")

end subroutine rrtm_kgb12
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

