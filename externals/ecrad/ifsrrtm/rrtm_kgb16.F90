! # 1 "ifsrrtm/rrtm_kgb16.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb16.f90"
subroutine rrtm_kgb16

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
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
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto16, only : kao,kbo ,selfrefo,forrefo ,fracrefao,fracrefbo,kao_d,kbo_d
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
! # 29 "ifsrrtm/rrtm_kgb16.f90" 2

if (lhook) call dr_hook('rrtm_kgb16',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d,kbo_d
  close(nulrad,err=1000)
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb16:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb16:')
endif

! planck fraction mapping level: p = 387.6100 mbar, t = 250.17 k
      fracrefao(:, 1) = (/ &
     &  1.1593e-01_jprb,2.3390e-01_jprb,1.9120e-01_jprb,1.3121e-01_jprb,1.0590e-01_jprb,8.4852e-02_jprb, &
     &  6.4168e-02_jprb,4.2537e-02_jprb,2.3220e-02_jprb,2.1767e-03_jprb,1.8203e-03_jprb,1.3724e-03_jprb, &
     &  9.5452e-04_jprb,5.5015e-04_jprb,1.9348e-04_jprb,2.7344e-05_jprb/)
      fracrefao(:, 2) = (/ &
     &  2.8101e-01_jprb,1.9773e-01_jprb,1.4749e-01_jprb,1.1399e-01_jprb,8.8190e-02_jprb,7.0531e-02_jprb, &
     &  4.6356e-02_jprb,3.0774e-02_jprb,1.7332e-02_jprb,2.0054e-03_jprb,1.5950e-03_jprb,1.2760e-03_jprb, &
     &  9.5034e-04_jprb,5.4992e-04_jprb,1.9349e-04_jprb,2.7309e-05_jprb/)
      fracrefao(:, 3) = (/ &
     &  2.9054e-01_jprb,2.1263e-01_jprb,1.4133e-01_jprb,1.1083e-01_jprb,8.5107e-02_jprb,6.5247e-02_jprb, &
     &  4.4542e-02_jprb,2.7205e-02_jprb,1.6495e-02_jprb,1.8453e-03_jprb,1.5222e-03_jprb,1.1884e-03_jprb, &
     &  8.1094e-04_jprb,4.9173e-04_jprb,1.9344e-04_jprb,2.7286e-05_jprb/)
      fracrefao(:, 4) = (/ &
     &  2.9641e-01_jprb,2.1738e-01_jprb,1.4228e-01_jprb,1.0830e-01_jprb,8.2837e-02_jprb,6.1359e-02_jprb, &
     &  4.4683e-02_jprb,2.5027e-02_jprb,1.6057e-02_jprb,1.7558e-03_jprb,1.4193e-03_jprb,1.0970e-03_jprb, &
     &  7.8281e-04_jprb,4.3260e-04_jprb,1.4837e-04_jprb,2.2958e-05_jprb/)
      fracrefao(:, 5) = (/ &
     &  2.9553e-01_jprb,2.2139e-01_jprb,1.4816e-01_jprb,1.0601e-01_jprb,8.0048e-02_jprb,6.0082e-02_jprb, &
     &  4.3952e-02_jprb,2.3788e-02_jprb,1.5734e-02_jprb,1.6586e-03_jprb,1.3434e-03_jprb,1.0281e-03_jprb, &
     &  7.0256e-04_jprb,4.2577e-04_jprb,1.2803e-04_jprb,1.3315e-05_jprb/)
      fracrefao(:, 6) = (/ &
     &  2.9313e-01_jprb,2.2476e-01_jprb,1.5470e-01_jprb,1.0322e-01_jprb,7.8904e-02_jprb,5.8175e-02_jprb, &
     &  4.3097e-02_jprb,2.3618e-02_jprb,1.5385e-02_jprb,1.5942e-03_jprb,1.2702e-03_jprb,9.5566e-04_jprb, &
     &  6.5421e-04_jprb,4.0165e-04_jprb,1.2805e-04_jprb,1.3355e-05_jprb/)
      fracrefao(:, 7) = (/ &
     &  2.9069e-01_jprb,2.2823e-01_jprb,1.5995e-01_jprb,1.0170e-01_jprb,7.7287e-02_jprb,5.6780e-02_jprb, &
     &  4.1752e-02_jprb,2.3899e-02_jprb,1.4937e-02_jprb,1.4916e-03_jprb,1.1909e-03_jprb,9.1307e-04_jprb, &
     &  6.3518e-04_jprb,3.9866e-04_jprb,1.2805e-04_jprb,1.3298e-05_jprb/)
      fracrefao(:, 8) = (/ &
     &  2.8446e-01_jprb,2.2651e-01_jprb,1.7133e-01_jprb,1.0299e-01_jprb,7.4231e-02_jprb,5.6031e-02_jprb, &
     &  4.1368e-02_jprb,2.4318e-02_jprb,1.4135e-02_jprb,1.4216e-03_jprb,1.1465e-03_jprb,8.9800e-04_jprb, &
     &  6.3553e-04_jprb,3.9536e-04_jprb,1.2749e-04_jprb,1.3298e-05_jprb/)
      fracrefao(:, 9) = (/ &
     &  2.0568e-01_jprb,2.5049e-01_jprb,2.0568e-01_jprb,1.1781e-01_jprb,7.5579e-02_jprb,5.8136e-02_jprb, &
     &  4.2397e-02_jprb,2.6544e-02_jprb,1.3067e-02_jprb,1.4061e-03_jprb,1.1455e-03_jprb,8.9408e-04_jprb, &
     &  6.3652e-04_jprb,3.9450e-04_jprb,1.2841e-04_jprb,1.3315e-05_jprb/)

! planck fraction mapping level : p=95.58350 mb, t = 215.70 k
      fracrefbo(:) = (/ &
     &  1.8111e-01_jprb,2.2612e-01_jprb,1.6226e-01_jprb,1.1872e-01_jprb,9.9048e-02_jprb,8.0390e-02_jprb, &
     &  6.1648e-02_jprb,4.1704e-02_jprb,2.2976e-02_jprb,1.9263e-03_jprb,1.4694e-03_jprb,1.1498e-03_jprb, &
     &  7.9906e-04_jprb,4.8310e-04_jprb,1.6188e-04_jprb,2.2651e-05_jprb/)


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
     &5.1629e-06_jprb,7.7578e-06_jprb,1.9043e-05_jprb,1.4802e-04_jprb,2.2980e-04_jprb,2.8057e-04_jprb, &
     &3.2824e-04_jprb,3.4913e-04_jprb,3.6515e-04_jprb,3.8271e-04_jprb,3.7499e-04_jprb,3.6966e-04_jprb, &
     &3.7424e-04_jprb,3.8884e-04_jprb,3.7117e-04_jprb,4.3710e-04_jprb/)
      forrefo(2,:) = (/ &
     &5.0804e-06_jprb,1.3466e-05_jprb,7.2606e-05_jprb,1.6940e-04_jprb,2.1022e-04_jprb,2.5900e-04_jprb, &
     &2.9106e-04_jprb,3.2261e-04_jprb,3.2066e-04_jprb,3.5421e-04_jprb,3.7128e-04_jprb,3.8144e-04_jprb, &
     &3.7854e-04_jprb,3.8347e-04_jprb,3.8921e-04_jprb,3.7339e-04_jprb/)
      forrefo(3,:) = (/ &
     &5.4797e-05_jprb,1.0026e-04_jprb,1.2422e-04_jprb,1.6386e-04_jprb,1.8378e-04_jprb,1.9616e-04_jprb, &
     &2.0711e-04_jprb,2.2492e-04_jprb,2.5240e-04_jprb,2.6187e-04_jprb,2.6058e-04_jprb,2.4892e-04_jprb, &
     &2.6526e-04_jprb,3.2105e-04_jprb,3.6903e-04_jprb,3.7213e-04_jprb/)
      forrefo(4,:) = (/ &
     &4.2782e-05_jprb,1.4775e-04_jprb,1.4588e-04_jprb,1.6964e-04_jprb,1.6667e-04_jprb,1.7192e-04_jprb, &
     &1.9057e-04_jprb,2.0180e-04_jprb,2.1177e-04_jprb,2.2326e-04_jprb,2.3801e-04_jprb,2.9308e-04_jprb, &
     &3.1130e-04_jprb,3.1829e-04_jprb,3.5035e-04_jprb,3.7782e-04_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 1.27793e-03_jprb, 1.05944e-03_jprb, 8.78300e-04_jprb, 7.28133e-04_jprb, 6.03641e-04_jprb, &
     & 5.00434e-04_jprb, 4.14873e-04_jprb, 3.43940e-04_jprb, 2.85135e-04_jprb, 2.36384e-04_jprb/)
      selfrefo(:, 2) = (/ &
     & 1.42785e-03_jprb, 1.17602e-03_jprb, 9.68600e-04_jprb, 7.97765e-04_jprb, 6.57060e-04_jprb, &
     & 5.41172e-04_jprb, 4.45724e-04_jprb, 3.67110e-04_jprb, 3.02361e-04_jprb, 2.49033e-04_jprb/)
      selfrefo(:, 3) = (/ &
     & 2.94095e-03_jprb, 2.27102e-03_jprb, 1.75370e-03_jprb, 1.35422e-03_jprb, 1.04574e-03_jprb, &
     & 8.07525e-04_jprb, 6.23577e-04_jprb, 4.81530e-04_jprb, 3.71841e-04_jprb, 2.87138e-04_jprb/)
      selfrefo(:, 4) = (/ &
     & 3.94894e-03_jprb, 3.48184e-03_jprb, 3.07000e-03_jprb, 2.70687e-03_jprb, 2.38669e-03_jprb, &
     & 2.10439e-03_jprb, 1.85547e-03_jprb, 1.63600e-03_jprb, 1.44249e-03_jprb, 1.27187e-03_jprb/)
      selfrefo(:, 5) = (/ &
     & 4.19971e-03_jprb, 3.86333e-03_jprb, 3.55390e-03_jprb, 3.26925e-03_jprb, 3.00740e-03_jprb, &
     & 2.76652e-03_jprb, 2.54494e-03_jprb, 2.34110e-03_jprb, 2.15359e-03_jprb, 1.98110e-03_jprb/)
      selfrefo(:, 6) = (/ &
     & 4.95922e-03_jprb, 4.57134e-03_jprb, 4.21380e-03_jprb, 3.88422e-03_jprb, 3.58042e-03_jprb, &
     & 3.30038e-03_jprb, 3.04225e-03_jprb, 2.80430e-03_jprb, 2.58496e-03_jprb, 2.38278e-03_jprb/)
      selfrefo(:, 7) = (/ &
     & 5.27379e-03_jprb, 4.91005e-03_jprb, 4.57140e-03_jprb, 4.25611e-03_jprb, 3.96256e-03_jprb, &
     & 3.68925e-03_jprb, 3.43480e-03_jprb, 3.19790e-03_jprb, 2.97734e-03_jprb, 2.77199e-03_jprb/)
      selfrefo(:, 8) = (/ &
     & 5.75341e-03_jprb, 5.31533e-03_jprb, 4.91060e-03_jprb, 4.53669e-03_jprb, 4.19126e-03_jprb, &
     & 3.87212e-03_jprb, 3.57729e-03_jprb, 3.30490e-03_jprb, 3.05325e-03_jprb, 2.82077e-03_jprb/)
      selfrefo(:, 9) = (/ &
     & 5.49849e-03_jprb, 5.14295e-03_jprb, 4.81040e-03_jprb, 4.49935e-03_jprb, 4.20842e-03_jprb, &
     & 3.93629e-03_jprb, 3.68177e-03_jprb, 3.44370e-03_jprb, 3.22102e-03_jprb, 3.01275e-03_jprb/)
      selfrefo(:,10) = (/ &
     & 6.04962e-03_jprb, 5.60945e-03_jprb, 5.20130e-03_jprb, 4.82285e-03_jprb, 4.47194e-03_jprb, &
     & 4.14656e-03_jprb, 3.84485e-03_jprb, 3.56510e-03_jprb, 3.30570e-03_jprb, 3.06518e-03_jprb/)
      selfrefo(:,11) = (/ &
     & 6.40108e-03_jprb, 5.87551e-03_jprb, 5.39310e-03_jprb, 4.95029e-03_jprb, 4.54385e-03_jprb, &
     & 4.17077e-03_jprb, 3.82833e-03_jprb, 3.51400e-03_jprb, 3.22548e-03_jprb, 2.96065e-03_jprb/)
      selfrefo(:,12) = (/ &
     & 6.77938e-03_jprb, 6.15713e-03_jprb, 5.59200e-03_jprb, 5.07874e-03_jprb, 4.61259e-03_jprb, &
     & 4.18922e-03_jprb, 3.80472e-03_jprb, 3.45550e-03_jprb, 3.13834e-03_jprb, 2.85029e-03_jprb/)
      selfrefo(:,13) = (/ &
     & 6.90020e-03_jprb, 6.26766e-03_jprb, 5.69310e-03_jprb, 5.17121e-03_jprb, 4.69717e-03_jprb, &
     & 4.26658e-03_jprb, 3.87546e-03_jprb, 3.52020e-03_jprb, 3.19750e-03_jprb, 2.90439e-03_jprb/)
      selfrefo(:,14) = (/ &
     & 6.92759e-03_jprb, 6.32882e-03_jprb, 5.78180e-03_jprb, 5.28206e-03_jprb, 4.82552e-03_jprb, &
     & 4.40843e-03_jprb, 4.02740e-03_jprb, 3.67930e-03_jprb, 3.36129e-03_jprb, 3.07076e-03_jprb/)
      selfrefo(:,15) = (/ &
     & 7.54539e-03_jprb, 6.81161e-03_jprb, 6.14920e-03_jprb, 5.55120e-03_jprb, 5.01136e-03_jprb, &
     & 4.52402e-03_jprb, 4.08407e-03_jprb, 3.68690e-03_jprb, 3.32836e-03_jprb, 3.00468e-03_jprb/)
      selfrefo(:,16) = (/ &
     & 7.62039e-03_jprb, 7.10834e-03_jprb, 6.63070e-03_jprb, 6.18515e-03_jprb, 5.76955e-03_jprb, &
     & 5.38186e-03_jprb, 5.02023e-03_jprb, 4.68290e-03_jprb, 4.36823e-03_jprb, 4.07471e-03_jprb/)

if (lhook) call dr_hook('rrtm_kgb16',1,zhook_handle)
return

1000 continue
call abor1("rrtm_kgb16:error closing file radrrtm")
1001 continue
call abor1("rrtm_kgb16:error reading file radrrtm")

end subroutine rrtm_kgb16
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

