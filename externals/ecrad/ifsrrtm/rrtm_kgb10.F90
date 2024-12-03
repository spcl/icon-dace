! # 1 "ifsrrtm/rrtm_kgb10.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb10.f90"
subroutine rrtm_kgb10

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     abozzo 101306 updated to rrtmg v4.85
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use yommp0    , only : nproc, myproc
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto10, only : kao, kbo, kao_d, kbo_d, fracrefao, fracrefbo, selfrefo, forrefo

!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle



! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 29 "ifsrrtm/rrtm_kgb10.f90" 2

if (lhook) call dr_hook('rrtm_kgb10',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb10:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb10:')
endif

! planck fraction mapping level : p = 212.7250, t = 223.06 k
      fracrefao(:) = (/ &
     &  1.6909e-01_jprb, 1.5419e-01_jprb, 1.3999e-01_jprb, 1.2637e-01_jprb, &
     &  1.1429e-01_jprb, 9.9676e-02_jprb, 8.0093e-02_jprb, 6.0283e-02_jprb, &
     &  4.1077e-02_jprb, 4.4857e-03_jprb, 3.6545e-03_jprb, 2.9243e-03_jprb, &
     &  2.0407e-03_jprb, 1.2891e-03_jprb, 4.8767e-04_jprb, 6.7748e-05_jprb/)

! planck fraction mapping level : p = 95.58350 mb, t = 215.70 k
      fracrefbo(:) = (/ &
     &  1.7391e-01_jprb, 1.5680e-01_jprb, 1.4419e-01_jprb, 1.2672e-01_jprb, &
     &  1.0708e-01_jprb, 9.7034e-02_jprb, 7.8545e-02_jprb, 5.9784e-02_jprb, &
     &  4.0879e-02_jprb, 4.4704e-03_jprb, 3.7150e-03_jprb, 2.9038e-03_jprb, &
     &  2.1454e-03_jprb, 1.2802e-03_jprb, 4.8328e-04_jprb, 6.7378e-05_jprb/)


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


!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.0515e-02_jprb,1.4860e-02_jprb,1.7181e-02_jprb,1.6642e-02_jprb,1.6644e-02_jprb,1.5649e-02_jprb, &
     &1.7734e-02_jprb,1.7521e-02_jprb,1.7868e-02_jprb,1.8400e-02_jprb,1.9361e-02_jprb,2.1487e-02_jprb, &
     &2.0192e-02_jprb,1.6545e-02_jprb,2.0922e-02_jprb,2.0922e-02_jprb/)
      forrefo(2,:) = (/ &
     &1.0423e-02_jprb,1.4593e-02_jprb,1.6329e-02_jprb,1.7071e-02_jprb,1.7252e-02_jprb,1.6188e-02_jprb, &
     &1.7752e-02_jprb,1.7913e-02_jprb,1.7551e-02_jprb,1.8203e-02_jprb,1.7946e-02_jprb,1.9828e-02_jprb, &
     &2.1566e-02_jprb,1.9707e-02_jprb,2.0944e-02_jprb,2.0944e-02_jprb/)
      forrefo(3,:) = (/ &
     &9.2770e-03_jprb,1.2818e-02_jprb,1.7181e-02_jprb,1.7858e-02_jprb,1.7888e-02_jprb,1.7121e-02_jprb, &
     &1.8116e-02_jprb,1.8230e-02_jprb,1.7719e-02_jprb,1.7833e-02_jprb,1.8438e-02_jprb,1.7995e-02_jprb, &
     &2.0895e-02_jprb,2.1525e-02_jprb,2.0517e-02_jprb,2.0954e-02_jprb/)
      forrefo(4,:) = (/ &
     &8.3290e-03_jprb,1.3483e-02_jprb,1.5432e-02_jprb,2.0793e-02_jprb,1.8404e-02_jprb,1.7470e-02_jprb, &
     &1.7253e-02_jprb,1.7132e-02_jprb,1.7119e-02_jprb,1.7376e-02_jprb,1.7030e-02_jprb,1.6847e-02_jprb, &
     &1.5562e-02_jprb,1.6836e-02_jprb,1.8746e-02_jprb,2.1233e-02_jprb/)

!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.41120e-01_jprb, 2.27071e-01_jprb, 2.13840e-01_jprb, 2.01380e-01_jprb, 1.89646e-01_jprb, &
     & 1.78596e-01_jprb, 1.68190e-01_jprb, 1.58390e-01_jprb, 1.49161e-01_jprb, 1.40470e-01_jprb/)
      selfrefo(:, 2) = (/ &
     & 3.11156e-01_jprb, 2.92249e-01_jprb, 2.74490e-01_jprb, 2.57810e-01_jprb, 2.42144e-01_jprb, &
     & 2.27430e-01_jprb, 2.13610e-01_jprb, 2.00630e-01_jprb, 1.88439e-01_jprb, 1.76988e-01_jprb/)
      selfrefo(:, 3) = (/ &
     & 3.37148e-01_jprb, 3.17767e-01_jprb, 2.99500e-01_jprb, 2.82283e-01_jprb, 2.66056e-01_jprb, &
     & 2.50762e-01_jprb, 2.36347e-01_jprb, 2.22760e-01_jprb, 2.09955e-01_jprb, 1.97885e-01_jprb/)
      selfrefo(:, 4) = (/ &
     & 3.57139e-01_jprb, 3.32763e-01_jprb, 3.10050e-01_jprb, 2.88888e-01_jprb, 2.69170e-01_jprb, &
     & 2.50798e-01_jprb, 2.33680e-01_jprb, 2.17730e-01_jprb, 2.02869e-01_jprb, 1.89022e-01_jprb/)
      selfrefo(:, 5) = (/ &
     & 3.60626e-01_jprb, 3.35433e-01_jprb, 3.12000e-01_jprb, 2.90204e-01_jprb, 2.69931e-01_jprb, &
     & 2.51074e-01_jprb, 2.33534e-01_jprb, 2.17220e-01_jprb, 2.02045e-01_jprb, 1.87931e-01_jprb/)
      selfrefo(:, 6) = (/ &
     & 3.42420e-01_jprb, 3.18795e-01_jprb, 2.96800e-01_jprb, 2.76323e-01_jprb, 2.57258e-01_jprb, &
     & 2.39509e-01_jprb, 2.22985e-01_jprb, 2.07600e-01_jprb, 1.93277e-01_jprb, 1.79942e-01_jprb/)
      selfrefo(:, 7) = (/ &
     & 3.65491e-01_jprb, 3.41599e-01_jprb, 3.19270e-01_jprb, 2.98400e-01_jprb, 2.78895e-01_jprb, &
     & 2.60664e-01_jprb, 2.43625e-01_jprb, 2.27700e-01_jprb, 2.12816e-01_jprb, 1.98905e-01_jprb/)
      selfrefo(:, 8) = (/ &
     & 3.70354e-01_jprb, 3.45005e-01_jprb, 3.21390e-01_jprb, 2.99392e-01_jprb, 2.78899e-01_jprb, &
     & 2.59809e-01_jprb, 2.42026e-01_jprb, 2.25460e-01_jprb, 2.10028e-01_jprb, 1.95652e-01_jprb/)
      selfrefo(:, 9) = (/ &
     & 3.60483e-01_jprb, 3.37846e-01_jprb, 3.16630e-01_jprb, 2.96747e-01_jprb, 2.78112e-01_jprb, &
     & 2.60648e-01_jprb, 2.44280e-01_jprb, 2.28940e-01_jprb, 2.14563e-01_jprb, 2.01090e-01_jprb/)
      selfrefo(:,10) = (/ &
     & 3.71845e-01_jprb, 3.48164e-01_jprb, 3.25990e-01_jprb, 3.05229e-01_jprb, 2.85790e-01_jprb, &
     & 2.67588e-01_jprb, 2.50547e-01_jprb, 2.34590e-01_jprb, 2.19650e-01_jprb, 2.05661e-01_jprb/)
      selfrefo(:,11) = (/ &
     & 3.60606e-01_jprb, 3.40789e-01_jprb, 3.22060e-01_jprb, 3.04361e-01_jprb, 2.87634e-01_jprb, &
     & 2.71826e-01_jprb, 2.56888e-01_jprb, 2.42770e-01_jprb, 2.29428e-01_jprb, 2.16819e-01_jprb/)
      selfrefo(:,12) = (/ &
     & 3.90046e-01_jprb, 3.68879e-01_jprb, 3.48860e-01_jprb, 3.29928e-01_jprb, 3.12023e-01_jprb, &
     & 2.95089e-01_jprb, 2.79075e-01_jprb, 2.63930e-01_jprb, 2.49607e-01_jprb, 2.36061e-01_jprb/)
      selfrefo(:,13) = (/ &
     & 4.38542e-01_jprb, 4.05139e-01_jprb, 3.74280e-01_jprb, 3.45771e-01_jprb, 3.19434e-01_jprb, &
     & 2.95103e-01_jprb, 2.72626e-01_jprb, 2.51860e-01_jprb, 2.32676e-01_jprb, 2.14953e-01_jprb/)
      selfrefo(:,14) = (/ &
     & 4.19448e-01_jprb, 3.81920e-01_jprb, 3.47750e-01_jprb, 3.16637e-01_jprb, 2.88307e-01_jprb, &
     & 2.62513e-01_jprb, 2.39026e-01_jprb, 2.17640e-01_jprb, 1.98168e-01_jprb, 1.80438e-01_jprb/)
      selfrefo(:,15) = (/ &
     & 4.20276e-01_jprb, 3.92281e-01_jprb, 3.66150e-01_jprb, 3.41760e-01_jprb, 3.18995e-01_jprb, &
     & 2.97746e-01_jprb, 2.77912e-01_jprb, 2.59400e-01_jprb, 2.42121e-01_jprb, 2.25993e-01_jprb/)
      selfrefo(:,16) = (/ &
     & 4.20276e-01_jprb, 3.92281e-01_jprb, 3.66150e-01_jprb, 3.41760e-01_jprb, 3.18995e-01_jprb, &
     & 2.97746e-01_jprb, 2.77912e-01_jprb, 2.59400e-01_jprb, 2.42121e-01_jprb, 2.25993e-01_jprb/)

if (lhook) call dr_hook('rrtm_kgb10',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb10:error reading file radrrtm")

end subroutine rrtm_kgb10
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

