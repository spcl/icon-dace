! # 1 "ifsrrtm/rrtm_kgb2.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb2.f90"
subroutine rrtm_kgb2

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 2:  250-500 cm-1 (low - h2o; high - h2o)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     abozzo may 2013 update to rrtmg v4.85
!     band 2:  350-500 cm-1 
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto2 , only : kao     ,kbo     ,selfrefo   ,fracrefao  ,&
 & fracrefbo  ,forrefo  ,kao_d, kbo_d
use yommp0    , only : nproc, myproc

!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 30 "ifsrrtm/rrtm_kgb2.f90" 2

if (lhook) call dr_hook('rrtm_kgb2',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb2:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb2:')
endif


! planck fraction mapping level: p = 1053.630 mbar, t = 294.2 k
      fracrefao(:) = (/ &
      &  1.6388e-01_jprb, 1.5241e-01_jprb, 1.4290e-01_jprb, 1.2864e-01_jprb, &
      &  1.1615e-01_jprb, 1.0047e-01_jprb, 8.0013e-02_jprb, 6.0445e-02_jprb, &
      &  4.0530e-02_jprb, 4.3879e-03_jprb, 3.5726e-03_jprb, 2.7669e-03_jprb, &
      &  2.0078e-03_jprb, 1.2864e-03_jprb, 4.7630e-04_jprb, 6.9109e-05_jprb/)

! planck fraction mapping level: p = 3.206e-2 mb, t = 197.92 k
      fracrefbo(:) = (/ &
      &  1.4697e-01_jprb, 1.4826e-01_jprb, 1.4278e-01_jprb, 1.3320e-01_jprb, &
      &  1.1965e-01_jprb, 1.0297e-01_jprb, 8.4170e-02_jprb, 6.3282e-02_jprb, &
      &  4.2868e-02_jprb, 4.6644e-03_jprb, 3.8619e-03_jprb, 3.0533e-03_jprb, &
      &  2.2359e-03_jprb, 1.4226e-03_jprb, 5.3642e-04_jprb, 7.6316e-05_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     & 2.8549e-03_jprb,4.8281e-03_jprb,6.2570e-03_jprb,8.2731e-03_jprb,7.9056e-03_jprb,7.7840e-03_jprb, &
     & 1.0115e-02_jprb,9.6599e-03_jprb,1.0153e-02_jprb,1.0921e-02_jprb,1.2408e-02_jprb,1.3496e-02_jprb, &
     & 1.5059e-02_jprb,1.4636e-02_jprb,1.6483e-02_jprb,1.2394e-02_jprb/)
      forrefo(2,:) = (/ &
     & 3.0036e-03_jprb,5.1093e-03_jprb,5.7317e-03_jprb,9.2246e-03_jprb,8.9829e-03_jprb,8.6477e-03_jprb, &
     & 1.1448e-02_jprb,1.0391e-02_jprb,1.0211e-02_jprb,1.2921e-02_jprb,1.2726e-02_jprb,1.2426e-02_jprb, &
     & 1.4609e-02_jprb,1.5783e-02_jprb,1.6617e-02_jprb,1.6858e-02_jprb/)
      forrefo(3,:) = (/ &
     & 3.0771e-03_jprb,5.1206e-03_jprb,5.8426e-03_jprb,9.5727e-03_jprb,1.0338e-02_jprb,9.3737e-03_jprb, &
     & 1.2805e-02_jprb,1.1272e-02_jprb,1.1353e-02_jprb,1.1837e-02_jprb,1.1550e-02_jprb,1.3020e-02_jprb, &
     & 1.3536e-02_jprb,1.6226e-02_jprb,1.6039e-02_jprb,2.2578e-02_jprb/)
      forrefo(4,:) = (/ &
     & 3.3072e-03_jprb,5.0240e-03_jprb,6.8474e-03_jprb,8.2736e-03_jprb,8.6151e-03_jprb,8.6762e-03_jprb, &
     & 1.1476e-02_jprb,1.0246e-02_jprb,1.0819e-02_jprb,1.0640e-02_jprb,1.0545e-02_jprb,1.0533e-02_jprb, &
     & 1.0496e-02_jprb,1.0142e-02_jprb,9.7979e-03_jprb,1.5255e-02_jprb/)


!     the following are parameters related to the reference water vapor
!     mixing ratios by refparam(i) = refh2o(i) / (.002+refh2o(i)).
!     these parameters are used for the planck function interpolation.
!refparam( :) = (/&
! & 0.903661_jprb   , 0.859386_jprb   , 0.746542_jprb   , 0.580496_jprb   , 0.412889_jprb   ,&
! & 0.275283_jprb   , 0.162745_jprb   , 7.63929e-02_jprb, 1.82553e-02_jprb, 3.72432e-03_jprb, &
! & 2.14946e-03_jprb, 1.66320e-03_jprb, 1.59940e-03_jprb/)  

!     ------------------------------------------------------------------

!     the array kao contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  the first
!     index in the array, jt, which runs from 1 to 5, corresponds to 
!     different temperatures.  more specifically, jt = 3 means that the 
!     data are for the corresponding tref for this  pressure level, 
!     jt = 2 refers to the temperaturetref-15, jt = 1 is for tref-30, 
!     jt = 4 is for tref+15, and jt = 5 is for tref+30.  the second 
!     index, jp, runs from 1 to 13 and refers to the corresponding 
!     pressure level in pref (e.g. jp = 1 is for a pressure of 1053.63 mb).  
!     the third index, ig, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.



!     the array kbo contains absorption coefs at the 16 chosen g-values 
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



!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

    selfrefo(:, 1) = (/ &
     & 7.25695e-01_jprb, 6.53591e-01_jprb, 5.88650e-01_jprb, 5.30162e-01_jprb, 4.77485e-01_jprb, &
     & 4.30042e-01_jprb, 3.87313e-01_jprb, 3.48830e-01_jprb, 3.14170e-01_jprb, 2.82954e-01_jprb/)
      selfrefo(:, 2) = (/ &
     & 9.61996e-01_jprb, 8.77853e-01_jprb, 8.01070e-01_jprb, 7.31003e-01_jprb, 6.67064e-01_jprb, &
     & 6.08718e-01_jprb, 5.55476e-01_jprb, 5.06890e-01_jprb, 4.62554e-01_jprb, 4.22096e-01_jprb/)
      selfrefo(:, 3) = (/ &
     & 9.72584e-01_jprb, 9.02658e-01_jprb, 8.37760e-01_jprb, 7.77527e-01_jprb, 7.21626e-01_jprb, &
     & 6.69743e-01_jprb, 6.21591e-01_jprb, 5.76900e-01_jprb, 5.35423e-01_jprb, 4.96927e-01_jprb/)
      selfrefo(:, 4) = (/ &
     & 1.24790e+00_jprb, 1.14353e+00_jprb, 1.04790e+00_jprb, 9.60263e-01_jprb, 8.79956e-01_jprb, &
     & 8.06364e-01_jprb, 7.38927e-01_jprb, 6.77130e-01_jprb, 6.20501e-01_jprb, 5.68608e-01_jprb/)
      selfrefo(:, 5) = (/ &
     & 1.23574e+00_jprb, 1.12928e+00_jprb, 1.03200e+00_jprb, 9.43096e-01_jprb, 8.61851e-01_jprb, &
     & 7.87605e-01_jprb, 7.19755e-01_jprb, 6.57750e-01_jprb, 6.01087e-01_jprb, 5.49305e-01_jprb/)
      selfrefo(:, 6) = (/ &
     & 1.20921e+00_jprb, 1.10660e+00_jprb, 1.01270e+00_jprb, 9.26766e-01_jprb, 8.48124e-01_jprb, &
     & 7.76155e-01_jprb, 7.10293e-01_jprb, 6.50020e-01_jprb, 5.94861e-01_jprb, 5.44384e-01_jprb/)
      selfrefo(:, 7) = (/ &
     & 1.38112e+00_jprb, 1.26727e+00_jprb, 1.16280e+00_jprb, 1.06694e+00_jprb, 9.78990e-01_jprb, &
     & 8.98287e-01_jprb, 8.24236e-01_jprb, 7.56290e-01_jprb, 6.93945e-01_jprb, 6.36739e-01_jprb/)
      selfrefo(:, 8) = (/ &
     & 1.30321e+00_jprb, 1.20127e+00_jprb, 1.10730e+00_jprb, 1.02068e+00_jprb, 9.40840e-01_jprb, &
     & 8.67243e-01_jprb, 7.99403e-01_jprb, 7.36870e-01_jprb, 6.79229e-01_jprb, 6.26096e-01_jprb/)
      selfrefo(:, 9) = (/ &
     & 1.26713e+00_jprb, 1.17927e+00_jprb, 1.09750e+00_jprb, 1.02140e+00_jprb, 9.50575e-01_jprb, &
     & 8.84662e-01_jprb, 8.23319e-01_jprb, 7.66230e-01_jprb, 7.13099e-01_jprb, 6.63653e-01_jprb/)
      selfrefo(:,10) = (/ &
     & 1.49824e+00_jprb, 1.37053e+00_jprb, 1.25370e+00_jprb, 1.14683e+00_jprb, 1.04908e+00_jprb, &
     & 9.59651e-01_jprb, 8.77849e-01_jprb, 8.03020e-01_jprb, 7.34569e-01_jprb, 6.71954e-01_jprb/)
      selfrefo(:,11) = (/ &
     & 1.44786e+00_jprb, 1.34594e+00_jprb, 1.25120e+00_jprb, 1.16313e+00_jprb, 1.08125e+00_jprb, &
     & 1.00514e+00_jprb, 9.34392e-01_jprb, 8.68620e-01_jprb, 8.07477e-01_jprb, 7.50639e-01_jprb/)
      selfrefo(:,12) = (/ &
     & 1.38460e+00_jprb, 1.30437e+00_jprb, 1.22880e+00_jprb, 1.15760e+00_jprb, 1.09053e+00_jprb, &
     & 1.02735e+00_jprb, 9.67825e-01_jprb, 9.11750e-01_jprb, 8.58924e-01_jprb, 8.09159e-01_jprb/)
      selfrefo(:,13) = (/ &
     & 1.51953e+00_jprb, 1.42822e+00_jprb, 1.34240e+00_jprb, 1.26173e+00_jprb, 1.18592e+00_jprb, &
     & 1.11465e+00_jprb, 1.04768e+00_jprb, 9.84720e-01_jprb, 9.25548e-01_jprb, 8.69932e-01_jprb/)
      selfrefo(:,14) = (/ &
     & 1.62608e+00_jprb, 1.51021e+00_jprb, 1.40260e+00_jprb, 1.30266e+00_jprb, 1.20983e+00_jprb, &
     & 1.12363e+00_jprb, 1.04356e+00_jprb, 9.69200e-01_jprb, 9.00138e-01_jprb, 8.35998e-01_jprb/)
      selfrefo(:,15) = (/ &
     & 1.65383e+00_jprb, 1.54808e+00_jprb, 1.44910e+00_jprb, 1.35644e+00_jprb, 1.26971e+00_jprb, &
     & 1.18853e+00_jprb, 1.11254e+00_jprb, 1.04140e+00_jprb, 9.74813e-01_jprb, 9.12484e-01_jprb/)
      selfrefo(:,16) = (/ &
     & 1.78105e+00_jprb, 1.61421e+00_jprb, 1.46300e+00_jprb, 1.32595e+00_jprb, 1.20174e+00_jprb, &
     & 1.08917e+00_jprb, 9.87141e-01_jprb, 8.94670e-01_jprb, 8.10861e-01_jprb, 7.34904e-01_jprb/)

if (lhook) call dr_hook('rrtm_kgb2',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb2:error reading file radrrtm")

end subroutine rrtm_kgb2
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

