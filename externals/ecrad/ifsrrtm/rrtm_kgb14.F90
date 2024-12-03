! # 1 "ifsrrtm/rrtm_kgb14.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb14.f90"
subroutine rrtm_kgb14

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     abozzo 201306 updated to rrtmg v4.85
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use yommp0    , only : nproc, myproc
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto14, only : kao     ,kbo     ,selfrefo, forrefo ,fracrefao  ,fracrefbo, &
  &  kao_d, kbo_d

!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 29 "ifsrrtm/rrtm_kgb14.f90" 2

if (lhook) call dr_hook('rrtm_kgb14',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb14:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb14:')
endif

! planck fraction mapping level : p = 142.5940 mb, t = 215.70 k
      fracrefao(:) = (/ &
     &  1.9360e-01_jprb, 1.7276e-01_jprb, 1.4811e-01_jprb, 1.2238e-01_jprb, &
     &  1.0242e-01_jprb, 8.6830e-02_jprb, 7.1890e-02_jprb, 5.4030e-02_jprb, &
     &  3.5075e-02_jprb, 3.8052e-03_jprb, 3.1458e-03_jprb, 2.4873e-03_jprb, &
     &  1.8182e-03_jprb, 1.1563e-03_jprb, 4.3251e-04_jprb, 5.7744e-05_jprb/)

! planck fraction mapping level : p = 4.758820mb, t = 250.85 k
      fracrefbo(:) = (/ &
     &  1.8599e-01_jprb, 1.6646e-01_jprb, 1.4264e-01_jprb, 1.2231e-01_jprb, &
     &  1.0603e-01_jprb, 9.2014e-02_jprb, 7.5287e-02_jprb, 5.6758e-02_jprb, &
     &  3.8386e-02_jprb, 4.2139e-03_jprb, 3.5399e-03_jprb, 2.7381e-03_jprb, &
     &  1.9202e-03_jprb, 1.2083e-03_jprb, 4.5395e-04_jprb, 6.2699e-05_jprb/)

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
     &2.7075e-06_jprb,2.2609e-06_jprb,1.5633e-06_jprb,8.7484e-07_jprb,5.5470e-07_jprb,4.8456e-07_jprb, &
     &4.7463e-07_jprb,4.6154e-07_jprb,4.4425e-07_jprb,4.2960e-07_jprb,4.2626e-07_jprb,4.1715e-07_jprb, &
     &4.2607e-07_jprb,3.6616e-07_jprb,2.6366e-07_jprb,2.6029e-07_jprb/)
      forrefo(2,:) = (/ &
     &2.6759e-06_jprb,2.2237e-06_jprb,1.4466e-06_jprb,9.3032e-07_jprb,6.4927e-07_jprb,5.4809e-07_jprb, &
     &4.9504e-07_jprb,4.6305e-07_jprb,4.4873e-07_jprb,4.2146e-07_jprb,4.2176e-07_jprb,4.2812e-07_jprb, &
     &4.0529e-07_jprb,4.0969e-07_jprb,2.9442e-07_jprb,2.6821e-07_jprb/)
      forrefo(3,:) = (/ &
     &2.6608e-06_jprb,2.1140e-06_jprb,1.4838e-06_jprb,9.2083e-07_jprb,6.3350e-07_jprb,5.7195e-07_jprb, &
     &6.2253e-07_jprb,5.1783e-07_jprb,4.4749e-07_jprb,4.3261e-07_jprb,4.2553e-07_jprb,4.2175e-07_jprb, &
     &4.1085e-07_jprb,4.0358e-07_jprb,3.5340e-07_jprb,2.7191e-07_jprb/)
      forrefo(4,:) = (/ &
     &2.6412e-06_jprb,1.9814e-06_jprb,1.2672e-06_jprb,8.1129e-07_jprb,7.1447e-07_jprb,7.5026e-07_jprb, &
     &7.4386e-07_jprb,7.2759e-07_jprb,7.3583e-07_jprb,7.6493e-07_jprb,8.8959e-07_jprb,7.5534e-07_jprb, &
     &5.3734e-07_jprb,4.5572e-07_jprb,4.1676e-07_jprb,3.6198e-07_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 4.67262e-03_jprb, 3.95211e-03_jprb, 3.34270e-03_jprb, 2.82726e-03_jprb, 2.39130e-03_jprb, &
     & 2.02256e-03_jprb, 1.71069e-03_jprb, 1.44690e-03_jprb, 1.22379e-03_jprb, 1.03508e-03_jprb/)
      selfrefo(:, 2) = (/ &
     & 4.42593e-03_jprb, 3.73338e-03_jprb, 3.14920e-03_jprb, 2.65643e-03_jprb, 2.24076e-03_jprb, &
     & 1.89014e-03_jprb, 1.59438e-03_jprb, 1.34490e-03_jprb, 1.13446e-03_jprb, 9.56943e-04_jprb/)
      selfrefo(:, 3) = (/ &
     & 3.96072e-03_jprb, 3.33789e-03_jprb, 2.81300e-03_jprb, 2.37065e-03_jprb, 1.99786e-03_jprb, &
     & 1.68369e-03_jprb, 1.41893e-03_jprb, 1.19580e-03_jprb, 1.00776e-03_jprb, 8.49286e-04_jprb/)
      selfrefo(:, 4) = (/ &
     & 3.71833e-03_jprb, 3.10030e-03_jprb, 2.58500e-03_jprb, 2.15535e-03_jprb, 1.79711e-03_jprb, &
     & 1.49841e-03_jprb, 1.24936e-03_jprb, 1.04170e-03_jprb, 8.68558e-04_jprb, 7.24195e-04_jprb/)
      selfrefo(:, 5) = (/ &
     & 3.55755e-03_jprb, 2.95355e-03_jprb, 2.45210e-03_jprb, 2.03578e-03_jprb, 1.69015e-03_jprb, &
     & 1.40320e-03_jprb, 1.16497e-03_jprb, 9.67180e-04_jprb, 8.02973e-04_jprb, 6.66646e-04_jprb/)
      selfrefo(:, 6) = (/ &
     & 3.47601e-03_jprb, 2.88628e-03_jprb, 2.39660e-03_jprb, 1.99000e-03_jprb, 1.65238e-03_jprb, &
     & 1.37204e-03_jprb, 1.13927e-03_jprb, 9.45980e-04_jprb, 7.85487e-04_jprb, 6.52224e-04_jprb/)
      selfrefo(:, 7) = (/ &
     & 3.44479e-03_jprb, 2.86224e-03_jprb, 2.37820e-03_jprb, 1.97602e-03_jprb, 1.64185e-03_jprb, &
     & 1.36420e-03_jprb, 1.13350e-03_jprb, 9.41810e-04_jprb, 7.82539e-04_jprb, 6.50204e-04_jprb/)
      selfrefo(:, 8) = (/ &
     & 3.40154e-03_jprb, 2.82953e-03_jprb, 2.35370e-03_jprb, 1.95789e-03_jprb, 1.62864e-03_jprb, &
     & 1.35476e-03_jprb, 1.12694e-03_jprb, 9.37430e-04_jprb, 7.79788e-04_jprb, 6.48655e-04_jprb/)
      selfrefo(:, 9) = (/ &
     & 3.39380e-03_jprb, 2.82288e-03_jprb, 2.34800e-03_jprb, 1.95301e-03_jprb, 1.62446e-03_jprb, &
     & 1.35119e-03_jprb, 1.12389e-03_jprb, 9.34820e-04_jprb, 7.77560e-04_jprb, 6.46755e-04_jprb/)
      selfrefo(:,10) = (/ &
     & 3.37185e-03_jprb, 2.80654e-03_jprb, 2.33600e-03_jprb, 1.94435e-03_jprb, 1.61837e-03_jprb, &
     & 1.34704e-03_jprb, 1.12120e-03_jprb, 9.33220e-04_jprb, 7.76759e-04_jprb, 6.46530e-04_jprb/)
      selfrefo(:,11) = (/ &
     & 3.37924e-03_jprb, 2.81172e-03_jprb, 2.33950e-03_jprb, 1.94659e-03_jprb, 1.61967e-03_jprb, &
     & 1.34765e-03_jprb, 1.12132e-03_jprb, 9.33000e-04_jprb, 7.76306e-04_jprb, 6.45930e-04_jprb/)
      selfrefo(:,12) = (/ &
     & 3.39658e-03_jprb, 2.82289e-03_jprb, 2.34610e-03_jprb, 1.94984e-03_jprb, 1.62051e-03_jprb, &
     & 1.34680e-03_jprb, 1.11933e-03_jprb, 9.30270e-04_jprb, 7.73146e-04_jprb, 6.42561e-04_jprb/)
      selfrefo(:,13) = (/ &
     & 3.36070e-03_jprb, 2.79913e-03_jprb, 2.33140e-03_jprb, 1.94183e-03_jprb, 1.61735e-03_jprb, &
     & 1.34709e-03_jprb, 1.12199e-03_jprb, 9.34510e-04_jprb, 7.78354e-04_jprb, 6.48292e-04_jprb/)
      selfrefo(:,14) = (/ &
     & 3.40428e-03_jprb, 2.81994e-03_jprb, 2.33590e-03_jprb, 1.93495e-03_jprb, 1.60282e-03_jprb, &
     & 1.32770e-03_jprb, 1.09980e-03_jprb, 9.11020e-04_jprb, 7.54645e-04_jprb, 6.25111e-04_jprb/)
      selfrefo(:,15) = (/ &
     & 3.27075e-03_jprb, 2.70783e-03_jprb, 2.24180e-03_jprb, 1.85597e-03_jprb, 1.53655e-03_jprb, &
     & 1.27210e-03_jprb, 1.05317e-03_jprb, 8.71910e-04_jprb, 7.21849e-04_jprb, 5.97615e-04_jprb/)
      selfrefo(:,16) = (/ &
     & 3.23123e-03_jprb, 2.67891e-03_jprb, 2.22100e-03_jprb, 1.84136e-03_jprb, 1.52661e-03_jprb, &
     & 1.26567e-03_jprb, 1.04932e-03_jprb, 8.69960e-04_jprb, 7.21256e-04_jprb, 5.97970e-04_jprb/)

if (lhook) call dr_hook('rrtm_kgb14',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb14:error reading file radrrtm")

end subroutine rrtm_kgb14
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

