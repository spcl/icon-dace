! # 1 "ifsrrtm/rrtm_kgb7.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb7.f90"
subroutine rrtm_kgb7

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     abozzo 201306 updated rrtmg v4.85
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto7 , only : kao     ,kbo, kao_mco2     ,kbo_mco2     ,selfrefo   ,forrefo, &
 & fracrefao  , fracrefbo   , kao_d, kbo_d

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
! # 30 "ifsrrtm/rrtm_kgb7.f90" 2

if (lhook) call dr_hook('rrtm_kgb7',0,zhook_handle)

if( myproc==1 )then
  read(nulrad) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb7:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb7:')
endif


! planck fraction mapping level : p = 706.27 mb, t = 278.94 k
      fracrefao(:, 1) = (/ &
      & 1.6312e-01_jprb,1.4949e-01_jprb,1.4305e-01_jprb,1.3161e-01_jprb,1.1684e-01_jprb,9.9900e-02_jprb, &
      & 8.0912e-02_jprb,6.0203e-02_jprb,4.0149e-02_jprb,4.3365e-03_jprb,3.5844e-03_jprb,2.8019e-03_jprb, &
      & 2.0756e-03_jprb,1.3449e-03_jprb,5.0492e-04_jprb,7.1194e-05_jprb/)
      fracrefao(:, 2) = (/ &
      & 1.6329e-01_jprb,1.4989e-01_jprb,1.4328e-01_jprb,1.3101e-01_jprb,1.1691e-01_jprb,9.9754e-02_jprb, &
      & 8.0956e-02_jprb,5.9912e-02_jprb,4.0271e-02_jprb,4.3298e-03_jprb,3.5626e-03_jprb,2.8421e-03_jprb, &
      & 2.1031e-03_jprb,1.3360e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)
      fracrefao(:, 3) = (/ &
      & 1.6236e-01_jprb,1.5081e-01_jprb,1.4341e-01_jprb,1.3083e-01_jprb,1.1684e-01_jprb,9.9701e-02_jprb, &
      & 8.0956e-02_jprb,5.9884e-02_jprb,4.0245e-02_jprb,4.3837e-03_jprb,3.6683e-03_jprb,2.9250e-03_jprb, &
      & 2.0969e-03_jprb,1.3320e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)
      fracrefao(:, 4) = (/ &
      & 1.6096e-01_jprb,1.5183e-01_jprb,1.4354e-01_jprb,1.3081e-01_jprb,1.1687e-01_jprb,9.9619e-02_jprb, &
      & 8.0947e-02_jprb,5.9899e-02_jprb,4.0416e-02_jprb,4.4389e-03_jprb,3.7280e-03_jprb,2.9548e-03_jprb, &
      & 2.0977e-03_jprb,1.3305e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)
      fracrefao(:, 5) = (/ &
      & 1.5661e-01_jprb,1.5478e-01_jprb,1.4414e-01_jprb,1.3097e-01_jprb,1.1695e-01_jprb,9.9823e-02_jprb, &
      & 8.0750e-02_jprb,6.0100e-02_jprb,4.0741e-02_jprb,4.4598e-03_jprb,3.7366e-03_jprb,2.9521e-03_jprb, &
      & 2.0980e-03_jprb,1.3297e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)
      fracrefao(:, 6) = (/ &
      & 1.4879e-01_jprb,1.5853e-01_jprb,1.4586e-01_jprb,1.3162e-01_jprb,1.1729e-01_jprb,1.0031e-01_jprb, &
      & 8.0908e-02_jprb,6.0460e-02_jprb,4.1100e-02_jprb,4.4578e-03_jprb,3.7388e-03_jprb,2.9508e-03_jprb, &
      & 2.0986e-03_jprb,1.3288e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)
      fracrefao(:, 7) = (/ &
      & 1.4117e-01_jprb,1.4838e-01_jprb,1.4807e-01_jprb,1.3759e-01_jprb,1.2218e-01_jprb,1.0228e-01_jprb, &
      & 8.2130e-02_jprb,6.1546e-02_jprb,4.1522e-02_jprb,4.4577e-03_jprb,3.7428e-03_jprb,2.9475e-03_jprb, &
      & 2.0997e-03_jprb,1.3277e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)
      fracrefao(:, 8) = (/ &
      & 1.4018e-01_jprb,1.4207e-01_jprb,1.3919e-01_jprb,1.3332e-01_jprb,1.2325e-01_jprb,1.0915e-01_jprb, &
      & 9.0280e-02_jprb,6.5554e-02_jprb,4.1852e-02_jprb,4.4707e-03_jprb,3.7572e-03_jprb,2.9364e-03_jprb, &
      & 2.1023e-03_jprb,1.3249e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)
      fracrefao(:, 9) = (/ &
      & 1.4863e-01_jprb,1.4926e-01_jprb,1.4740e-01_jprb,1.3558e-01_jprb,1.1999e-01_jprb,1.0044e-01_jprb, &
      & 8.1927e-02_jprb,6.0989e-02_jprb,4.0665e-02_jprb,4.4481e-03_jprb,3.7369e-03_jprb,2.9482e-03_jprb, &
      & 2.0976e-03_jprb,1.3281e-03_jprb,4.8965e-04_jprb,6.8900e-05_jprb/)

! planck fraction mapping level : p=95.58 mbar, t= 215.70 k
fracrefbo(:) = (/ &
 &  1.5872e-01_jprb,1.5443e-01_jprb,1.4413e-01_jprb,1.3147e-01_jprb,1.1634e-01_jprb,9.8914e-02_jprb, &
 &  8.0236e-02_jprb,6.0197e-02_jprb,4.0624e-02_jprb,4.4225e-03_jprb,3.6688e-03_jprb,2.9074e-03_jprb, &
 &  2.0862e-03_jprb,1.3039e-03_jprb,4.8561e-04_jprb,6.8854e-05_jprb/)


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


!     the array kao_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   the first index in the array, js, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference mls column amount value of gas 1 
!     to that of gas2.  the second index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the third index 
!     runs over the g-channel (1 to 16).

      kao_mco2( 1, :, 1) = (/ &
     & 7.38630e-06_jprb, 8.97432e-06_jprb, 1.09037e-05_jprb, 1.32480e-05_jprb, 1.60963e-05_jprb, &
     & 1.95569e-05_jprb, 2.37615e-05_jprb, 2.88701e-05_jprb, 3.50770e-05_jprb, 4.26184e-05_jprb, &
     & 5.17811e-05_jprb, 6.29138e-05_jprb, 7.64400e-05_jprb, 9.28742e-05_jprb, 1.12842e-04_jprb, &
     & 1.37102e-04_jprb, 1.66578e-04_jprb, 2.02392e-04_jprb, 2.45905e-04_jprb/)
      kao_mco2( 2, :, 1) = (/ &
     & 7.03916e-06_jprb, 8.58785e-06_jprb, 1.04773e-05_jprb, 1.27824e-05_jprb, 1.55947e-05_jprb, &
     & 1.90257e-05_jprb, 2.32115e-05_jprb, 2.83183e-05_jprb, 3.45487e-05_jprb, 4.21498e-05_jprb, &
     & 5.14233e-05_jprb, 6.27370e-05_jprb, 7.65398e-05_jprb, 9.33794e-05_jprb, 1.13924e-04_jprb, &
     & 1.38989e-04_jprb, 1.69568e-04_jprb, 2.06874e-04_jprb, 2.52389e-04_jprb/)
      kao_mco2( 3, :, 1) = (/ &
     & 7.80015e-06_jprb, 9.48520e-06_jprb, 1.15343e-05_jprb, 1.40260e-05_jprb, 1.70560e-05_jprb, &
     & 2.07405e-05_jprb, 2.52211e-05_jprb, 3.06695e-05_jprb, 3.72950e-05_jprb, 4.53517e-05_jprb, &
     & 5.51489e-05_jprb, 6.70626e-05_jprb, 8.15499e-05_jprb, 9.91670e-05_jprb, 1.20590e-04_jprb, &
     & 1.46640e-04_jprb, 1.78319e-04_jprb, 2.16841e-04_jprb, 2.63684e-04_jprb/)
      kao_mco2( 4, :, 1) = (/ &
     & 9.24267e-06_jprb, 1.11747e-05_jprb, 1.35105e-05_jprb, 1.63346e-05_jprb, 1.97490e-05_jprb, &
     & 2.38771e-05_jprb, 2.88682e-05_jprb, 3.49025e-05_jprb, 4.21981e-05_jprb, 5.10188e-05_jprb, &
     & 6.16832e-05_jprb, 7.45768e-05_jprb, 9.01656e-05_jprb, 1.09013e-04_jprb, 1.31800e-04_jprb, &
     & 1.59350e-04_jprb, 1.92659e-04_jprb, 2.32930e-04_jprb, 2.81619e-04_jprb/)
      kao_mco2( 5, :, 1) = (/ &
     & 1.59506e-05_jprb, 1.90078e-05_jprb, 2.26509e-05_jprb, 2.69923e-05_jprb, 3.21658e-05_jprb, &
     & 3.83309e-05_jprb, 4.56777e-05_jprb, 5.44325e-05_jprb, 6.48654e-05_jprb, 7.72978e-05_jprb, &
     & 9.21132e-05_jprb, 1.09768e-04_jprb, 1.30807e-04_jprb, 1.55878e-04_jprb, 1.85755e-04_jprb, &
     & 2.21357e-04_jprb, 2.63784e-04_jprb, 3.14342e-04_jprb, 3.74591e-04_jprb/)
      kao_mco2( 6, :, 1) = (/ &
     & 3.53189e-05_jprb, 4.14789e-05_jprb, 4.87131e-05_jprb, 5.72092e-05_jprb, 6.71870e-05_jprb, &
     & 7.89050e-05_jprb, 9.26667e-05_jprb, 1.08829e-04_jprb, 1.27809e-04_jprb, 1.50100e-04_jprb, &
     & 1.76279e-04_jprb, 2.07024e-04_jprb, 2.43131e-04_jprb, 2.85535e-04_jprb, 3.35335e-04_jprb, &
     & 3.93821e-04_jprb, 4.62507e-04_jprb, 5.43172e-04_jprb, 6.37906e-04_jprb/)
      kao_mco2( 7, :, 1) = (/ &
     & 6.63273e-05_jprb, 7.76356e-05_jprb, 9.08718e-05_jprb, 1.06365e-04_jprb, 1.24499e-04_jprb, &
     & 1.45725e-04_jprb, 1.70570e-04_jprb, 1.99651e-04_jprb, 2.33689e-04_jprb, 2.73531e-04_jprb, &
     & 3.20166e-04_jprb, 3.74752e-04_jprb, 4.38644e-04_jprb, 5.13429e-04_jprb, 6.00964e-04_jprb, &
     & 7.03424e-04_jprb, 8.23352e-04_jprb, 9.63726e-04_jprb, 1.12803e-03_jprb/)
      kao_mco2( 8, :, 1) = (/ &
     & 9.01134e-05_jprb, 1.05517e-04_jprb, 1.23553e-04_jprb, 1.44673e-04_jprb, 1.69402e-04_jprb, &
     & 1.98359e-04_jprb, 2.32265e-04_jprb, 2.71967e-04_jprb, 3.18456e-04_jprb, 3.72890e-04_jprb, &
     & 4.36630e-04_jprb, 5.11265e-04_jprb, 5.98657e-04_jprb, 7.00989e-04_jprb, 8.20811e-04_jprb, &
     & 9.61116e-04_jprb, 1.12540e-03_jprb, 1.31777e-03_jprb, 1.54302e-03_jprb/)
      kao_mco2( 9, :, 1) = (/ &
     & 1.14205e-05_jprb, 1.36364e-05_jprb, 1.62823e-05_jprb, 1.94416e-05_jprb, 2.32139e-05_jprb, &
     & 2.77181e-05_jprb, 3.30963e-05_jprb, 3.95181e-05_jprb, 4.71858e-05_jprb, 5.63414e-05_jprb, &
     & 6.72734e-05_jprb, 8.03266e-05_jprb, 9.59124e-05_jprb, 1.14523e-04_jprb, 1.36743e-04_jprb, &
     & 1.63276e-04_jprb, 1.94957e-04_jprb, 2.32784e-04_jprb, 2.77952e-04_jprb/)
      kao_mco2( 1, :, 2) = (/ &
     & 2.01754e-05_jprb, 2.40506e-05_jprb, 2.86701e-05_jprb, 3.41769e-05_jprb, 4.07414e-05_jprb, &
     & 4.85668e-05_jprb, 5.78953e-05_jprb, 6.90155e-05_jprb, 8.22717e-05_jprb, 9.80739e-05_jprb, &
     & 1.16912e-04_jprb, 1.39367e-04_jprb, 1.66136e-04_jprb, 1.98047e-04_jprb, 2.36087e-04_jprb, &
     & 2.81433e-04_jprb, 3.35489e-04_jprb, 3.99928e-04_jprb, 4.76744e-04_jprb/)
      kao_mco2( 2, :, 2) = (/ &
     & 2.08613e-05_jprb, 2.48759e-05_jprb, 2.96631e-05_jprb, 3.53716e-05_jprb, 4.21786e-05_jprb, &
     & 5.02955e-05_jprb, 5.99746e-05_jprb, 7.15163e-05_jprb, 8.52791e-05_jprb, 1.01690e-04_jprb, &
     & 1.21260e-04_jprb, 1.44596e-04_jprb, 1.72422e-04_jprb, 2.05604e-04_jprb, 2.45171e-04_jprb, &
     & 2.92352e-04_jprb, 3.48613e-04_jprb, 4.15702e-04_jprb, 4.95700e-04_jprb/)
      kao_mco2( 3, :, 2) = (/ &
     & 2.06879e-05_jprb, 2.47009e-05_jprb, 2.94924e-05_jprb, 3.52133e-05_jprb, 4.20439e-05_jprb, &
     & 5.01995e-05_jprb, 5.99372e-05_jprb, 7.15637e-05_jprb, 8.54456e-05_jprb, 1.02020e-04_jprb, &
     & 1.21810e-04_jprb, 1.45439e-04_jprb, 1.73651e-04_jprb, 2.07335e-04_jprb, 2.47554e-04_jprb, &
     & 2.95574e-04_jprb, 3.52909e-04_jprb, 4.21366e-04_jprb, 5.03102e-04_jprb/)
      kao_mco2( 4, :, 2) = (/ &
     & 2.12700e-05_jprb, 2.54064e-05_jprb, 3.03472e-05_jprb, 3.62490e-05_jprb, 4.32984e-05_jprb, &
     & 5.17188e-05_jprb, 6.17767e-05_jprb, 7.37906e-05_jprb, 8.81410e-05_jprb, 1.05282e-04_jprb, &
     & 1.25757e-04_jprb, 1.50213e-04_jprb, 1.79425e-04_jprb, 2.14319e-04_jprb, 2.55998e-04_jprb, &
     & 3.05782e-04_jprb, 3.65249e-04_jprb, 4.36280e-04_jprb, 5.21125e-04_jprb/)
      kao_mco2( 5, :, 2) = (/ &
     & 1.88144e-05_jprb, 2.25220e-05_jprb, 2.69602e-05_jprb, 3.22730e-05_jprb, 3.86328e-05_jprb, &
     & 4.62458e-05_jprb, 5.53591e-05_jprb, 6.62682e-05_jprb, 7.93271e-05_jprb, 9.49594e-05_jprb, &
     & 1.13672e-04_jprb, 1.36073e-04_jprb, 1.62887e-04_jprb, 1.94986e-04_jprb, 2.33410e-04_jprb, &
     & 2.79406e-04_jprb, 3.34467e-04_jprb, 4.00377e-04_jprb, 4.79275e-04_jprb/)
      kao_mco2( 6, :, 2) = (/ &
     & 1.20964e-05_jprb, 1.46021e-05_jprb, 1.76268e-05_jprb, 2.12780e-05_jprb, 2.56856e-05_jprb, &
     & 3.10062e-05_jprb, 3.74289e-05_jprb, 4.51820e-05_jprb, 5.45411e-05_jprb, 6.58388e-05_jprb, &
     & 7.94769e-05_jprb, 9.59399e-05_jprb, 1.15813e-04_jprb, 1.39803e-04_jprb, 1.68762e-04_jprb, &
     & 2.03720e-04_jprb, 2.45919e-04_jprb, 2.96859e-04_jprb, 3.58350e-04_jprb/)
      kao_mco2( 7, :, 2) = (/ &
     & 3.07117e-05_jprb, 3.64441e-05_jprb, 4.32465e-05_jprb, 5.13186e-05_jprb, 6.08974e-05_jprb, &
     & 7.22642e-05_jprb, 8.57525e-05_jprb, 1.01758e-04_jprb, 1.20752e-04_jprb, 1.43291e-04_jprb, &
     & 1.70037e-04_jprb, 2.01775e-04_jprb, 2.39436e-04_jprb, 2.84128e-04_jprb, 3.37161e-04_jprb, &
     & 4.00094e-04_jprb, 4.74773e-04_jprb, 5.63391e-04_jprb, 6.68549e-04_jprb/)
      kao_mco2( 8, :, 2) = (/ &
     & 9.34077e-05_jprb, 1.10481e-04_jprb, 1.30675e-04_jprb, 1.54559e-04_jprb, 1.82810e-04_jprb, &
     & 2.16224e-04_jprb, 2.55745e-04_jprb, 3.02491e-04_jprb, 3.57780e-04_jprb, 4.23175e-04_jprb, &
     & 5.00523e-04_jprb, 5.92009e-04_jprb, 7.00217e-04_jprb, 8.28203e-04_jprb, 9.79582e-04_jprb, &
     & 1.15863e-03_jprb, 1.37041e-03_jprb, 1.62089e-03_jprb, 1.91716e-03_jprb/)
      kao_mco2( 9, :, 2) = (/ &
     & 1.15325e-05_jprb, 1.37935e-05_jprb, 1.64978e-05_jprb, 1.97322e-05_jprb, 2.36007e-05_jprb, &
     & 2.82277e-05_jprb, 3.37618e-05_jprb, 4.03808e-05_jprb, 4.82976e-05_jprb, 5.77664e-05_jprb, &
     & 6.90916e-05_jprb, 8.26372e-05_jprb, 9.88384e-05_jprb, 1.18216e-04_jprb, 1.41392e-04_jprb, &
     & 1.69113e-04_jprb, 2.02267e-04_jprb, 2.41922e-04_jprb, 2.89352e-04_jprb/)
      kao_mco2( 1, :, 3) = (/ &
     & 2.56142e-05_jprb, 3.05385e-05_jprb, 3.64096e-05_jprb, 4.34093e-05_jprb, 5.17547e-05_jprb, &
     & 6.17045e-05_jprb, 7.35672e-05_jprb, 8.77104e-05_jprb, 1.04573e-04_jprb, 1.24677e-04_jprb, &
     & 1.48646e-04_jprb, 1.77223e-04_jprb, 2.11294e-04_jprb, 2.51915e-04_jprb, 3.00346e-04_jprb, &
     & 3.58087e-04_jprb, 4.26929e-04_jprb, 5.09006e-04_jprb, 6.06862e-04_jprb/)
      kao_mco2( 2, :, 3) = (/ &
     & 2.49802e-05_jprb, 2.98040e-05_jprb, 3.55593e-05_jprb, 4.24259e-05_jprb, 5.06186e-05_jprb, &
     & 6.03932e-05_jprb, 7.20554e-05_jprb, 8.59696e-05_jprb, 1.02571e-04_jprb, 1.22377e-04_jprb, &
     & 1.46009e-04_jprb, 1.74204e-04_jprb, 2.07844e-04_jprb, 2.47979e-04_jprb, 2.95865e-04_jprb, &
     & 3.52998e-04_jprb, 4.21163e-04_jprb, 5.02491e-04_jprb, 5.99524e-04_jprb/)
      kao_mco2( 3, :, 3) = (/ &
     & 2.54644e-05_jprb, 3.03959e-05_jprb, 3.62825e-05_jprb, 4.33091e-05_jprb, 5.16965e-05_jprb, &
     & 6.17083e-05_jprb, 7.36589e-05_jprb, 8.79240e-05_jprb, 1.04952e-04_jprb, 1.25277e-04_jprb, &
     & 1.49539e-04_jprb, 1.78499e-04_jprb, 2.13068e-04_jprb, 2.54331e-04_jprb, 3.03586e-04_jprb, &
     & 3.62380e-04_jprb, 4.32560e-04_jprb, 5.16331e-04_jprb, 6.16326e-04_jprb/)
      kao_mco2( 4, :, 3) = (/ &
     & 2.55054e-05_jprb, 3.04699e-05_jprb, 3.64007e-05_jprb, 4.34859e-05_jprb, 5.19501e-05_jprb, &
     & 6.20619e-05_jprb, 7.41418e-05_jprb, 8.85731e-05_jprb, 1.05813e-04_jprb, 1.26409e-04_jprb, &
     & 1.51014e-04_jprb, 1.80408e-04_jprb, 2.15523e-04_jprb, 2.57474e-04_jprb, 3.07589e-04_jprb, &
     & 3.67459e-04_jprb, 4.38983e-04_jprb, 5.24428e-04_jprb, 6.26505e-04_jprb/)
      kao_mco2( 5, :, 3) = (/ &
     & 2.48615e-05_jprb, 2.97398e-05_jprb, 3.55754e-05_jprb, 4.25560e-05_jprb, 5.09064e-05_jprb, &
     & 6.08952e-05_jprb, 7.28441e-05_jprb, 8.71375e-05_jprb, 1.04236e-04_jprb, 1.24689e-04_jprb, &
     & 1.49155e-04_jprb, 1.78423e-04_jprb, 2.13433e-04_jprb, 2.55313e-04_jprb, 3.05410e-04_jprb, &
     & 3.65338e-04_jprb, 4.37024e-04_jprb, 5.22777e-04_jprb, 6.25356e-04_jprb/)
      kao_mco2( 6, :, 3) = (/ &
     & 2.09074e-05_jprb, 2.50891e-05_jprb, 3.01072e-05_jprb, 3.61290e-05_jprb, 4.33553e-05_jprb, &
     & 5.20269e-05_jprb, 6.24329e-05_jprb, 7.49202e-05_jprb, 8.99051e-05_jprb, 1.07887e-04_jprb, &
     & 1.29466e-04_jprb, 1.55361e-04_jprb, 1.86435e-04_jprb, 2.23724e-04_jprb, 2.68471e-04_jprb, &
     & 3.22169e-04_jprb, 3.86607e-04_jprb, 4.63933e-04_jprb, 5.56725e-04_jprb/)
      kao_mco2( 7, :, 3) = (/ &
     & 1.25163e-05_jprb, 1.51688e-05_jprb, 1.83835e-05_jprb, 2.22795e-05_jprb, 2.70011e-05_jprb, &
     & 3.27234e-05_jprb, 3.96583e-05_jprb, 4.80630e-05_jprb, 5.82488e-05_jprb, 7.05933e-05_jprb, &
     & 8.55539e-05_jprb, 1.03685e-04_jprb, 1.25659e-04_jprb, 1.52289e-04_jprb, 1.84563e-04_jprb, &
     & 2.23677e-04_jprb, 2.71081e-04_jprb, 3.28530e-04_jprb, 3.98154e-04_jprb/)
      kao_mco2( 8, :, 3) = (/ &
     & 1.00408e-04_jprb, 1.20081e-04_jprb, 1.43608e-04_jprb, 1.71745e-04_jprb, 2.05395e-04_jprb, &
     & 2.45637e-04_jprb, 2.93765e-04_jprb, 3.51322e-04_jprb, 4.20156e-04_jprb, 5.02476e-04_jprb, &
     & 6.00926e-04_jprb, 7.18665e-04_jprb, 8.59472e-04_jprb, 1.02787e-03_jprb, 1.22926e-03_jprb, &
     & 1.47010e-03_jprb, 1.75814e-03_jprb, 2.10261e-03_jprb, 2.51457e-03_jprb/)
      kao_mco2( 9, :, 3) = (/ &
     & 8.50402e-06_jprb, 1.02737e-05_jprb, 1.24116e-05_jprb, 1.49945e-05_jprb, 1.81148e-05_jprb, &
     & 2.18844e-05_jprb, 2.64385e-05_jprb, 3.19403e-05_jprb, 3.85871e-05_jprb, 4.66169e-05_jprb, &
     & 5.63178e-05_jprb, 6.80375e-05_jprb, 8.21959e-05_jprb, 9.93008e-05_jprb, 1.19965e-04_jprb, &
     & 1.44930e-04_jprb, 1.75089e-04_jprb, 2.11525e-04_jprb, 2.55543e-04_jprb/)
      kao_mco2( 1, :, 4) = (/ &
     & 2.68659e-05_jprb, 3.20986e-05_jprb, 3.83506e-05_jprb, 4.58203e-05_jprb, 5.47450e-05_jprb, &
     & 6.54078e-05_jprb, 7.81476e-05_jprb, 9.33687e-05_jprb, 1.11555e-04_jprb, 1.33282e-04_jprb, &
     & 1.59242e-04_jprb, 1.90259e-04_jprb, 2.27316e-04_jprb, 2.71592e-04_jprb, 3.24490e-04_jprb, &
     & 3.87693e-04_jprb, 4.63206e-04_jprb, 5.53426e-04_jprb, 6.61218e-04_jprb/)
      kao_mco2( 2, :, 4) = (/ &
     & 2.74827e-05_jprb, 3.28460e-05_jprb, 3.92560e-05_jprb, 4.69169e-05_jprb, 5.60728e-05_jprb, &
     & 6.70155e-05_jprb, 8.00937e-05_jprb, 9.57241e-05_jprb, 1.14405e-04_jprb, 1.36731e-04_jprb, &
     & 1.63415e-04_jprb, 1.95305e-04_jprb, 2.33419e-04_jprb, 2.78972e-04_jprb, 3.33413e-04_jprb, &
     & 3.98480e-04_jprb, 4.76244e-04_jprb, 5.69184e-04_jprb, 6.80261e-04_jprb/)
      kao_mco2( 3, :, 4) = (/ &
     & 2.84702e-05_jprb, 3.40189e-05_jprb, 4.06490e-05_jprb, 4.85713e-05_jprb, 5.80375e-05_jprb, &
     & 6.93487e-05_jprb, 8.28644e-05_jprb, 9.90142e-05_jprb, 1.18312e-04_jprb, 1.41370e-04_jprb, &
     & 1.68922e-04_jprb, 2.01844e-04_jprb, 2.41182e-04_jprb, 2.88188e-04_jprb, 3.44354e-04_jprb, &
     & 4.11466e-04_jprb, 4.91659e-04_jprb, 5.87481e-04_jprb, 7.01977e-04_jprb/)
      kao_mco2( 4, :, 4) = (/ &
     & 2.92293e-05_jprb, 3.49243e-05_jprb, 4.17289e-05_jprb, 4.98593e-05_jprb, 5.95738e-05_jprb, &
     & 7.11810e-05_jprb, 8.50498e-05_jprb, 1.01621e-04_jprb, 1.21420e-04_jprb, 1.45078e-04_jprb, &
     & 1.73344e-04_jprb, 2.07119e-04_jprb, 2.47473e-04_jprb, 2.95690e-04_jprb, 3.53302e-04_jprb, &
     & 4.22139e-04_jprb, 5.04388e-04_jprb, 6.02662e-04_jprb, 7.20083e-04_jprb/)
      kao_mco2( 5, :, 4) = (/ &
     & 2.88531e-05_jprb, 3.45646e-05_jprb, 4.14067e-05_jprb, 4.96033e-05_jprb, 5.94224e-05_jprb, &
     & 7.11851e-05_jprb, 8.52764e-05_jprb, 1.02157e-04_jprb, 1.22379e-04_jprb, 1.46604e-04_jprb, &
     & 1.75625e-04_jprb, 2.10391e-04_jprb, 2.52038e-04_jprb, 3.01929e-04_jprb, 3.61697e-04_jprb, &
     & 4.33295e-04_jprb, 5.19067e-04_jprb, 6.21818e-04_jprb, 7.44908e-04_jprb/)
      kao_mco2( 6, :, 4) = (/ &
     & 2.79869e-05_jprb, 3.36885e-05_jprb, 4.05516e-05_jprb, 4.88130e-05_jprb, 5.87574e-05_jprb, &
     & 7.07278e-05_jprb, 8.51368e-05_jprb, 1.02481e-04_jprb, 1.23359e-04_jprb, 1.48490e-04_jprb, &
     & 1.78742e-04_jprb, 2.15156e-04_jprb, 2.58988e-04_jprb, 3.11751e-04_jprb, 3.75262e-04_jprb, &
     & 4.51712e-04_jprb, 5.43737e-04_jprb, 6.54510e-04_jprb, 7.87849e-04_jprb/)
      kao_mco2( 7, :, 4) = (/ &
     & 1.45797e-05_jprb, 1.78204e-05_jprb, 2.17815e-05_jprb, 2.66230e-05_jprb, 3.25407e-05_jprb, &
     & 3.97737e-05_jprb, 4.86145e-05_jprb, 5.94203e-05_jprb, 7.26281e-05_jprb, 8.87715e-05_jprb, &
     & 1.08503e-04_jprb, 1.32621e-04_jprb, 1.62100e-04_jprb, 1.98130e-04_jprb, 2.42170e-04_jprb, &
     & 2.95999e-04_jprb, 3.61792e-04_jprb, 4.42210e-04_jprb, 5.40503e-04_jprb/)
      kao_mco2( 8, :, 4) = (/ &
     & 6.32607e-05_jprb, 7.63420e-05_jprb, 9.21282e-05_jprb, 1.11179e-04_jprb, 1.34169e-04_jprb, &
     & 1.61913e-04_jprb, 1.95393e-04_jprb, 2.35797e-04_jprb, 2.84557e-04_jprb, 3.43398e-04_jprb, &
     & 4.14407e-04_jprb, 5.00100e-04_jprb, 6.03512e-04_jprb, 7.28308e-04_jprb, 8.78909e-04_jprb, &
     & 1.06065e-03_jprb, 1.27998e-03_jprb, 1.54466e-03_jprb, 1.86407e-03_jprb/)
      kao_mco2( 9, :, 4) = (/ &
     & 1.52296e-05_jprb, 1.84301e-05_jprb, 2.23032e-05_jprb, 2.69902e-05_jprb, 3.26622e-05_jprb, &
     & 3.95261e-05_jprb, 4.78324e-05_jprb, 5.78844e-05_jprb, 7.00487e-05_jprb, 8.47694e-05_jprb, &
     & 1.02584e-04_jprb, 1.24142e-04_jprb, 1.50230e-04_jprb, 1.81800e-04_jprb, 2.20005e-04_jprb, &
     & 2.66239e-04_jprb, 3.22190e-04_jprb, 3.89897e-04_jprb, 4.71833e-04_jprb/)
      kao_mco2( 1, :, 5) = (/ &
     & 3.43213e-05_jprb, 4.11301e-05_jprb, 4.92896e-05_jprb, 5.90679e-05_jprb, 7.07860e-05_jprb, &
     & 8.48288e-05_jprb, 1.01657e-04_jprb, 1.21825e-04_jprb, 1.45993e-04_jprb, 1.74955e-04_jprb, &
     & 2.09663e-04_jprb, 2.51257e-04_jprb, 3.01103e-04_jprb, 3.60837e-04_jprb, 4.32421e-04_jprb, &
     & 5.18206e-04_jprb, 6.21010e-04_jprb, 7.44208e-04_jprb, 8.91846e-04_jprb/)
      kao_mco2( 2, :, 5) = (/ &
     & 3.14792e-05_jprb, 3.79075e-05_jprb, 4.56485e-05_jprb, 5.49703e-05_jprb, 6.61956e-05_jprb, &
     & 7.97133e-05_jprb, 9.59914e-05_jprb, 1.15594e-04_jprb, 1.39199e-04_jprb, 1.67624e-04_jprb, &
     & 2.01854e-04_jprb, 2.43075e-04_jprb, 2.92712e-04_jprb, 3.52487e-04_jprb, 4.24467e-04_jprb, &
     & 5.11147e-04_jprb, 6.15527e-04_jprb, 7.41222e-04_jprb, 8.92585e-04_jprb/)
      kao_mco2( 3, :, 5) = (/ &
     & 3.21655e-05_jprb, 3.87990e-05_jprb, 4.68006e-05_jprb, 5.64523e-05_jprb, 6.80945e-05_jprb, &
     & 8.21377e-05_jprb, 9.90770e-05_jprb, 1.19510e-04_jprb, 1.44156e-04_jprb, 1.73886e-04_jprb, &
     & 2.09746e-04_jprb, 2.53002e-04_jprb, 3.05179e-04_jprb, 3.68117e-04_jprb, 4.44033e-04_jprb, &
     & 5.35607e-04_jprb, 6.46066e-04_jprb, 7.79304e-04_jprb, 9.40020e-04_jprb/)
      kao_mco2( 4, :, 5) = (/ &
     & 3.22870e-05_jprb, 3.89864e-05_jprb, 4.70759e-05_jprb, 5.68439e-05_jprb, 6.86388e-05_jprb, &
     & 8.28810e-05_jprb, 1.00078e-04_jprb, 1.20844e-04_jprb, 1.45919e-04_jprb, 1.76196e-04_jprb, &
     & 2.12756e-04_jprb, 2.56902e-04_jprb, 3.10207e-04_jprb, 3.74574e-04_jprb, 4.52296e-04_jprb, &
     & 5.46146e-04_jprb, 6.59468e-04_jprb, 7.96304e-04_jprb, 9.61533e-04_jprb/)
      kao_mco2( 5, :, 5) = (/ &
     & 3.31190e-05_jprb, 3.99528e-05_jprb, 4.81967e-05_jprb, 5.81417e-05_jprb, 7.01387e-05_jprb, &
     & 8.46111e-05_jprb, 1.02070e-04_jprb, 1.23131e-04_jprb, 1.48538e-04_jprb, 1.79187e-04_jprb, &
     & 2.16161e-04_jprb, 2.60764e-04_jprb, 3.14570e-04_jprb, 3.79479e-04_jprb, 4.57781e-04_jprb, &
     & 5.52240e-04_jprb, 6.66190e-04_jprb, 8.03652e-04_jprb, 9.69477e-04_jprb/)
      kao_mco2( 6, :, 5) = (/ &
     & 3.31287e-05_jprb, 3.99772e-05_jprb, 4.82413e-05_jprb, 5.82139e-05_jprb, 7.02480e-05_jprb, &
     & 8.47698e-05_jprb, 1.02294e-04_jprb, 1.23440e-04_jprb, 1.48958e-04_jprb, 1.79750e-04_jprb, &
     & 2.16909e-04_jprb, 2.61749e-04_jprb, 3.15858e-04_jprb, 3.81153e-04_jprb, 4.59945e-04_jprb, &
     & 5.55026e-04_jprb, 6.69762e-04_jprb, 8.08216e-04_jprb, 9.75292e-04_jprb/)
      kao_mco2( 7, :, 5) = (/ &
     & 3.35235e-05_jprb, 4.02832e-05_jprb, 4.84061e-05_jprb, 5.81668e-05_jprb, 6.98958e-05_jprb, &
     & 8.39898e-05_jprb, 1.00926e-04_jprb, 1.21277e-04_jprb, 1.45731e-04_jprb, 1.75117e-04_jprb, &
     & 2.10428e-04_jprb, 2.52860e-04_jprb, 3.03847e-04_jprb, 3.65116e-04_jprb, 4.38739e-04_jprb, &
     & 5.27208e-04_jprb, 6.33516e-04_jprb, 7.61260e-04_jprb, 9.14762e-04_jprb/)
      kao_mco2( 8, :, 5) = (/ &
     & 3.57666e-05_jprb, 4.27511e-05_jprb, 5.10995e-05_jprb, 6.10783e-05_jprb, 7.30057e-05_jprb, &
     & 8.72622e-05_jprb, 1.04303e-04_jprb, 1.24671e-04_jprb, 1.49017e-04_jprb, 1.78117e-04_jprb, &
     & 2.12900e-04_jprb, 2.54475e-04_jprb, 3.04169e-04_jprb, 3.63567e-04_jprb, 4.34565e-04_jprb, &
     & 5.19427e-04_jprb, 6.20861e-04_jprb, 7.42103e-04_jprb, 8.87020e-04_jprb/)
      kao_mco2( 9, :, 5) = (/ &
     & 2.96349e-05_jprb, 3.55202e-05_jprb, 4.25743e-05_jprb, 5.10292e-05_jprb, 6.11633e-05_jprb, &
     & 7.33099e-05_jprb, 8.78687e-05_jprb, 1.05319e-04_jprb, 1.26234e-04_jprb, 1.51304e-04_jprb, &
     & 1.81352e-04_jprb, 2.17367e-04_jprb, 2.60534e-04_jprb, 3.12275e-04_jprb, 3.74290e-04_jprb, &
     & 4.48622e-04_jprb, 5.37715e-04_jprb, 6.44502e-04_jprb, 7.72495e-04_jprb/)
      kao_mco2( 1, :, 6) = (/ &
     & 4.14659e-05_jprb, 4.98693e-05_jprb, 5.99757e-05_jprb, 7.21302e-05_jprb, 8.67479e-05_jprb, &
     & 1.04328e-04_jprb, 1.25471e-04_jprb, 1.50899e-04_jprb, 1.81479e-04_jprb, 2.18257e-04_jprb, &
     & 2.62489e-04_jprb, 3.15685e-04_jprb, 3.79660e-04_jprb, 4.56601e-04_jprb, 5.49135e-04_jprb, &
     & 6.60422e-04_jprb, 7.94261e-04_jprb, 9.55224e-04_jprb, 1.14881e-03_jprb/)
      kao_mco2( 2, :, 6) = (/ &
     & 4.25940e-05_jprb, 5.11162e-05_jprb, 6.13434e-05_jprb, 7.36168e-05_jprb, 8.83459e-05_jprb, &
     & 1.06022e-04_jprb, 1.27235e-04_jprb, 1.52691e-04_jprb, 1.83241e-04_jprb, 2.19904e-04_jprb, &
     & 2.63902e-04_jprb, 3.16703e-04_jprb, 3.80068e-04_jprb, 4.56111e-04_jprb, 5.47368e-04_jprb, &
     & 6.56885e-04_jprb, 7.88313e-04_jprb, 9.46036e-04_jprb, 1.13532e-03_jprb/)
      kao_mco2( 3, :, 6) = (/ &
     & 4.44940e-05_jprb, 5.32922e-05_jprb, 6.38303e-05_jprb, 7.64522e-05_jprb, 9.15700e-05_jprb, &
     & 1.09677e-04_jprb, 1.31365e-04_jprb, 1.57341e-04_jprb, 1.88454e-04_jprb, 2.25719e-04_jprb, &
     & 2.70353e-04_jprb, 3.23813e-04_jprb, 3.87844e-04_jprb, 4.64537e-04_jprb, 5.56395e-04_jprb, &
     & 6.66418e-04_jprb, 7.98196e-04_jprb, 9.56032e-04_jprb, 1.14508e-03_jprb/)
      kao_mco2( 4, :, 6) = (/ &
     & 4.83402e-05_jprb, 5.78065e-05_jprb, 6.91265e-05_jprb, 8.26633e-05_jprb, 9.88510e-05_jprb, &
     & 1.18209e-04_jprb, 1.41357e-04_jprb, 1.69038e-04_jprb, 2.02140e-04_jprb, 2.41725e-04_jprb, &
     & 2.89061e-04_jprb, 3.45667e-04_jprb, 4.13357e-04_jprb, 4.94303e-04_jprb, 5.91101e-04_jprb, &
     & 7.06854e-04_jprb, 8.45275e-04_jprb, 1.01080e-03_jprb, 1.20874e-03_jprb/)
      kao_mco2( 5, :, 6) = (/ &
     & 5.14797e-05_jprb, 6.15328e-05_jprb, 7.35491e-05_jprb, 8.79120e-05_jprb, 1.05080e-04_jprb, &
     & 1.25600e-04_jprb, 1.50128e-04_jprb, 1.79445e-04_jprb, 2.14487e-04_jprb, 2.56373e-04_jprb, &
     & 3.06438e-04_jprb, 3.66281e-04_jprb, 4.37809e-04_jprb, 5.23305e-04_jprb, 6.25498e-04_jprb, &
     & 7.47647e-04_jprb, 8.93650e-04_jprb, 1.06816e-03_jprb, 1.27676e-03_jprb/)
      kao_mco2( 6, :, 6) = (/ &
     & 5.71481e-05_jprb, 6.83156e-05_jprb, 8.16652e-05_jprb, 9.76237e-05_jprb, 1.16701e-04_jprb, &
     & 1.39505e-04_jprb, 1.66766e-04_jprb, 1.99354e-04_jprb, 2.38311e-04_jprb, 2.84879e-04_jprb, &
     & 3.40548e-04_jprb, 4.07096e-04_jprb, 4.86647e-04_jprb, 5.81744e-04_jprb, 6.95424e-04_jprb, &
     & 8.31318e-04_jprb, 9.93769e-04_jprb, 1.18796e-03_jprb, 1.42010e-03_jprb/)
      kao_mco2( 7, :, 6) = (/ &
     & 5.69513e-05_jprb, 6.84420e-05_jprb, 8.22512e-05_jprb, 9.88466e-05_jprb, 1.18790e-04_jprb, &
     & 1.42758e-04_jprb, 1.71562e-04_jprb, 2.06177e-04_jprb, 2.47776e-04_jprb, 2.97768e-04_jprb, &
     & 3.57848e-04_jprb, 4.30049e-04_jprb, 5.16817e-04_jprb, 6.21093e-04_jprb, 7.46407e-04_jprb, &
     & 8.97006e-04_jprb, 1.07799e-03_jprb, 1.29549e-03_jprb, 1.55687e-03_jprb/)
      kao_mco2( 8, :, 6) = (/ &
     & 4.39361e-06_jprb, 5.50076e-06_jprb, 6.88690e-06_jprb, 8.62235e-06_jprb, 1.07951e-05_jprb, &
     & 1.35154e-05_jprb, 1.69212e-05_jprb, 2.11851e-05_jprb, 2.65236e-05_jprb, 3.32074e-05_jprb, &
     & 4.15754e-05_jprb, 5.20520e-05_jprb, 6.51687e-05_jprb, 8.15907e-05_jprb, 1.02151e-04_jprb, &
     & 1.27892e-04_jprb, 1.60120e-04_jprb, 2.00469e-04_jprb, 2.50985e-04_jprb/)
      kao_mco2( 9, :, 6) = (/ &
     & 5.75515e-05_jprb, 6.86850e-05_jprb, 8.19722e-05_jprb, 9.78298e-05_jprb, 1.16755e-04_jprb, &
     & 1.39342e-04_jprb, 1.66297e-04_jprb, 1.98468e-04_jprb, 2.36862e-04_jprb, 2.82683e-04_jprb, &
     & 3.37369e-04_jprb, 4.02633e-04_jprb, 4.80523e-04_jprb, 5.73481e-04_jprb, 6.84422e-04_jprb, &
     & 8.16824e-04_jprb, 9.74841e-04_jprb, 1.16342e-03_jprb, 1.38849e-03_jprb/)
      kao_mco2( 1, :, 7) = (/ &
     & 6.84544e-05_jprb, 8.16461e-05_jprb, 9.73799e-05_jprb, 1.16146e-04_jprb, 1.38528e-04_jprb, &
     & 1.65223e-04_jprb, 1.97063e-04_jprb, 2.35039e-04_jprb, 2.80333e-04_jprb, 3.34355e-04_jprb, &
     & 3.98788e-04_jprb, 4.75637e-04_jprb, 5.67296e-04_jprb, 6.76618e-04_jprb, 8.07008e-04_jprb, &
     & 9.62525e-04_jprb, 1.14801e-03_jprb, 1.36924e-03_jprb, 1.63310e-03_jprb/)
      kao_mco2( 2, :, 7) = (/ &
     & 6.88332e-05_jprb, 8.21719e-05_jprb, 9.80955e-05_jprb, 1.17105e-04_jprb, 1.39798e-04_jprb, &
     & 1.66888e-04_jprb, 1.99229e-04_jprb, 2.37836e-04_jprb, 2.83925e-04_jprb, 3.38944e-04_jprb, &
     & 4.04627e-04_jprb, 4.83037e-04_jprb, 5.76641e-04_jprb, 6.88385e-04_jprb, 8.21782e-04_jprb, &
     & 9.81031e-04_jprb, 1.17114e-03_jprb, 1.39809e-03_jprb, 1.66901e-03_jprb/)
      kao_mco2( 3, :, 7) = (/ &
     & 7.49899e-05_jprb, 8.94606e-05_jprb, 1.06724e-04_jprb, 1.27318e-04_jprb, 1.51887e-04_jprb, &
     & 1.81196e-04_jprb, 2.16161e-04_jprb, 2.57873e-04_jprb, 3.07635e-04_jprb, 3.66999e-04_jprb, &
     & 4.37818e-04_jprb, 5.22304e-04_jprb, 6.23092e-04_jprb, 7.43330e-04_jprb, 8.86769e-04_jprb, &
     & 1.05789e-03_jprb, 1.26203e-03_jprb, 1.50556e-03_jprb, 1.79608e-03_jprb/)
      kao_mco2( 4, :, 7) = (/ &
     & 8.26801e-05_jprb, 9.85802e-05_jprb, 1.17538e-04_jprb, 1.40141e-04_jprb, 1.67092e-04_jprb, &
     & 1.99225e-04_jprb, 2.37537e-04_jprb, 2.83217e-04_jprb, 3.37682e-04_jprb, 4.02621e-04_jprb, &
     & 4.80048e-04_jprb, 5.72365e-04_jprb, 6.82435e-04_jprb, 8.13673e-04_jprb, 9.70148e-04_jprb, &
     & 1.15671e-03_jprb, 1.37916e-03_jprb, 1.64438e-03_jprb, 1.96061e-03_jprb/)
      kao_mco2( 5, :, 7) = (/ &
     & 9.29561e-05_jprb, 1.10845e-04_jprb, 1.32176e-04_jprb, 1.57612e-04_jprb, 1.87944e-04_jprb, &
     & 2.24112e-04_jprb, 2.67241e-04_jprb, 3.18669e-04_jprb, 3.79995e-04_jprb, 4.53121e-04_jprb, &
     & 5.40321e-04_jprb, 6.44302e-04_jprb, 7.68293e-04_jprb, 9.16146e-04_jprb, 1.09245e-03_jprb, &
     & 1.30268e-03_jprb, 1.55338e-03_jprb, 1.85231e-03_jprb, 2.20877e-03_jprb/)
      kao_mco2( 6, :, 7) = (/ &
     & 1.09700e-04_jprb, 1.30879e-04_jprb, 1.56148e-04_jprb, 1.86294e-04_jprb, 2.22261e-04_jprb, &
     & 2.65172e-04_jprb, 3.16367e-04_jprb, 3.77446e-04_jprb, 4.50317e-04_jprb, 5.37257e-04_jprb, &
     & 6.40983e-04_jprb, 7.64734e-04_jprb, 9.12376e-04_jprb, 1.08852e-03_jprb, 1.29868e-03_jprb, &
     & 1.54941e-03_jprb, 1.84854e-03_jprb, 2.20543e-03_jprb, 2.63122e-03_jprb/)
      kao_mco2( 7, :, 7) = (/ &
     & 1.43457e-04_jprb, 1.71554e-04_jprb, 2.05153e-04_jprb, 2.45332e-04_jprb, 2.93381e-04_jprb, &
     & 3.50840e-04_jprb, 4.19552e-04_jprb, 5.01722e-04_jprb, 5.99985e-04_jprb, 7.17492e-04_jprb, &
     & 8.58014e-04_jprb, 1.02606e-03_jprb, 1.22701e-03_jprb, 1.46732e-03_jprb, 1.75470e-03_jprb, &
     & 2.09836e-03_jprb, 2.50933e-03_jprb, 3.00078e-03_jprb, 3.58849e-03_jprb/)
      kao_mco2( 8, :, 7) = (/ &
     & 1.52152e-05_jprb, 1.89421e-05_jprb, 2.35819e-05_jprb, 2.93582e-05_jprb, 3.65494e-05_jprb, &
     & 4.55021e-05_jprb, 5.66476e-05_jprb, 7.05233e-05_jprb, 8.77978e-05_jprb, 1.09304e-04_jprb, &
     & 1.36077e-04_jprb, 1.69409e-04_jprb, 2.10905e-04_jprb, 2.62565e-04_jprb, 3.26880e-04_jprb, &
     & 4.06948e-04_jprb, 5.06629e-04_jprb, 6.30726e-04_jprb, 7.85219e-04_jprb/)
      kao_mco2( 9, :, 7) = (/ &
     & 1.15683e-04_jprb, 1.37544e-04_jprb, 1.63535e-04_jprb, 1.94438e-04_jprb, 2.31180e-04_jprb, &
     & 2.74866e-04_jprb, 3.26807e-04_jprb, 3.88563e-04_jprb, 4.61989e-04_jprb, 5.49289e-04_jprb, &
     & 6.53088e-04_jprb, 7.76501e-04_jprb, 9.23234e-04_jprb, 1.09770e-03_jprb, 1.30512e-03_jprb, &
     & 1.55175e-03_jprb, 1.84498e-03_jprb, 2.19362e-03_jprb, 2.60815e-03_jprb/)
      kao_mco2( 1, :, 8) = (/ &
     & 1.18154e-04_jprb, 1.40516e-04_jprb, 1.67111e-04_jprb, 1.98739e-04_jprb, 2.36353e-04_jprb, &
     & 2.81086e-04_jprb, 3.34285e-04_jprb, 3.97553e-04_jprb, 4.72796e-04_jprb, 5.62278e-04_jprb, &
     & 6.68697e-04_jprb, 7.95257e-04_jprb, 9.45770e-04_jprb, 1.12477e-03_jprb, 1.33765e-03_jprb, &
     & 1.59081e-03_jprb, 1.89190e-03_jprb, 2.24996e-03_jprb, 2.67580e-03_jprb/)
      kao_mco2( 2, :, 8) = (/ &
     & 1.40874e-04_jprb, 1.67009e-04_jprb, 1.97993e-04_jprb, 2.34726e-04_jprb, 2.78273e-04_jprb, &
     & 3.29899e-04_jprb, 3.91102e-04_jprb, 4.63661e-04_jprb, 5.49680e-04_jprb, 6.51659e-04_jprb, &
     & 7.72556e-04_jprb, 9.15884e-04_jprb, 1.08580e-03_jprb, 1.28724e-03_jprb, 1.52605e-03_jprb, &
     & 1.80917e-03_jprb, 2.14482e-03_jprb, 2.54273e-03_jprb, 3.01446e-03_jprb/)
      kao_mco2( 3, :, 8) = (/ &
     & 1.55092e-04_jprb, 1.84132e-04_jprb, 2.18609e-04_jprb, 2.59542e-04_jprb, 3.08140e-04_jprb, &
     & 3.65837e-04_jprb, 4.34337e-04_jprb, 5.15664e-04_jprb, 6.12219e-04_jprb, 7.26853e-04_jprb, &
     & 8.62952e-04_jprb, 1.02453e-03_jprb, 1.21637e-03_jprb, 1.44413e-03_jprb, 1.71453e-03_jprb, &
     & 2.03557e-03_jprb, 2.41671e-03_jprb, 2.86923e-03_jprb, 3.40647e-03_jprb/)
      kao_mco2( 4, :, 8) = (/ &
     & 1.80666e-04_jprb, 2.14521e-04_jprb, 2.54721e-04_jprb, 3.02454e-04_jprb, 3.59131e-04_jprb, &
     & 4.26429e-04_jprb, 5.06339e-04_jprb, 6.01223e-04_jprb, 7.13887e-04_jprb, 8.47663e-04_jprb, &
     & 1.00651e-03_jprb, 1.19512e-03_jprb, 1.41908e-03_jprb, 1.68500e-03_jprb, 2.00076e-03_jprb, &
     & 2.37568e-03_jprb, 2.82087e-03_jprb, 3.34947e-03_jprb, 3.97714e-03_jprb/)
      kao_mco2( 5, :, 8) = (/ &
     & 2.21554e-04_jprb, 2.63265e-04_jprb, 3.12829e-04_jprb, 3.71724e-04_jprb, 4.41707e-04_jprb, &
     & 5.24865e-04_jprb, 6.23679e-04_jprb, 7.41096e-04_jprb, 8.80619e-04_jprb, 1.04641e-03_jprb, &
     & 1.24341e-03_jprb, 1.47750e-03_jprb, 1.75567e-03_jprb, 2.08620e-03_jprb, 2.47896e-03_jprb, &
     & 2.94566e-03_jprb, 3.50023e-03_jprb, 4.15920e-03_jprb, 4.94224e-03_jprb/)
      kao_mco2( 6, :, 8) = (/ &
     & 2.78997e-04_jprb, 3.32548e-04_jprb, 3.96378e-04_jprb, 4.72460e-04_jprb, 5.63146e-04_jprb, &
     & 6.71238e-04_jprb, 8.00077e-04_jprb, 9.53647e-04_jprb, 1.13669e-03_jprb, 1.35487e-03_jprb, &
     & 1.61493e-03_jprb, 1.92491e-03_jprb, 2.29438e-03_jprb, 2.73477e-03_jprb, 3.25969e-03_jprb, &
     & 3.88537e-03_jprb, 4.63114e-03_jprb, 5.52005e-03_jprb, 6.57958e-03_jprb/)
      kao_mco2( 7, :, 8) = (/ &
     & 2.84939e-04_jprb, 3.40606e-04_jprb, 4.07149e-04_jprb, 4.86691e-04_jprb, 5.81774e-04_jprb, &
     & 6.95432e-04_jprb, 8.31295e-04_jprb, 9.93700e-04_jprb, 1.18783e-03_jprb, 1.41989e-03_jprb, &
     & 1.69729e-03_jprb, 2.02888e-03_jprb, 2.42526e-03_jprb, 2.89907e-03_jprb, 3.46544e-03_jprb, &
     & 4.14246e-03_jprb, 4.95176e-03_jprb, 5.91915e-03_jprb, 7.07554e-03_jprb/)
      kao_mco2( 8, :, 8) = (/ &
     & 5.30764e-05_jprb, 6.47812e-05_jprb, 7.90673e-05_jprb, 9.65039e-05_jprb, 1.17786e-04_jprb, &
     & 1.43761e-04_jprb, 1.75464e-04_jprb, 2.14159e-04_jprb, 2.61387e-04_jprb, 3.19030e-04_jprb, &
     & 3.89385e-04_jprb, 4.75255e-04_jprb, 5.80062e-04_jprb, 7.07982e-04_jprb, 8.64111e-04_jprb, &
     & 1.05467e-03_jprb, 1.28726e-03_jprb, 1.57113e-03_jprb, 1.91761e-03_jprb/)
      kao_mco2( 9, :, 8) = (/ &
     & 2.76806e-04_jprb, 3.29639e-04_jprb, 3.92556e-04_jprb, 4.67481e-04_jprb, 5.56708e-04_jprb, &
     & 6.62964e-04_jprb, 7.89501e-04_jprb, 9.40190e-04_jprb, 1.11964e-03_jprb, 1.33334e-03_jprb, &
     & 1.58783e-03_jprb, 1.89089e-03_jprb, 2.25180e-03_jprb, 2.68159e-03_jprb, 3.19341e-03_jprb, &
     & 3.80293e-03_jprb, 4.52878e-03_jprb, 5.39316e-03_jprb, 6.42253e-03_jprb/)
      kao_mco2( 1, :, 9) = (/ &
     & 3.30614e-04_jprb, 3.93289e-04_jprb, 4.67844e-04_jprb, 5.56534e-04_jprb, 6.62036e-04_jprb, &
     & 7.87539e-04_jprb, 9.36833e-04_jprb, 1.11443e-03_jprb, 1.32569e-03_jprb, 1.57700e-03_jprb, &
     & 1.87596e-03_jprb, 2.23158e-03_jprb, 2.65463e-03_jprb, 3.15787e-03_jprb, 3.75650e-03_jprb, &
     & 4.46862e-03_jprb, 5.31575e-03_jprb, 6.32345e-03_jprb, 7.52219e-03_jprb/)
      kao_mco2( 2, :, 9) = (/ &
     & 3.78453e-04_jprb, 4.50735e-04_jprb, 5.36824e-04_jprb, 6.39355e-04_jprb, 7.61469e-04_jprb, &
     & 9.06906e-04_jprb, 1.08012e-03_jprb, 1.28642e-03_jprb, 1.53212e-03_jprb, 1.82475e-03_jprb, &
     & 2.17326e-03_jprb, 2.58835e-03_jprb, 3.08271e-03_jprb, 3.67149e-03_jprb, 4.37273e-03_jprb, &
     & 5.20790e-03_jprb, 6.20259e-03_jprb, 7.38725e-03_jprb, 8.79818e-03_jprb/)
      kao_mco2( 3, :, 9) = (/ &
     & 4.57576e-04_jprb, 5.45512e-04_jprb, 6.50348e-04_jprb, 7.75330e-04_jprb, 9.24332e-04_jprb, &
     & 1.10197e-03_jprb, 1.31374e-03_jprb, 1.56621e-03_jprb, 1.86721e-03_jprb, 2.22604e-03_jprb, &
     & 2.65384e-03_jprb, 3.16385e-03_jprb, 3.77187e-03_jprb, 4.49675e-03_jprb, 5.36092e-03_jprb, &
     & 6.39117e-03_jprb, 7.61942e-03_jprb, 9.08370e-03_jprb, 1.08294e-02_jprb/)
      kao_mco2( 4, :, 9) = (/ &
     & 5.18277e-04_jprb, 6.18764e-04_jprb, 7.38735e-04_jprb, 8.81967e-04_jprb, 1.05297e-03_jprb, &
     & 1.25713e-03_jprb, 1.50087e-03_jprb, 1.79187e-03_jprb, 2.13929e-03_jprb, 2.55407e-03_jprb, &
     & 3.04928e-03_jprb, 3.64050e-03_jprb, 4.34635e-03_jprb, 5.18905e-03_jprb, 6.19514e-03_jprb, &
     & 7.39631e-03_jprb, 8.83036e-03_jprb, 1.05425e-02_jprb, 1.25865e-02_jprb/)
      kao_mco2( 5, :, 9) = (/ &
     & 4.45365e-04_jprb, 5.32106e-04_jprb, 6.35742e-04_jprb, 7.59563e-04_jprb, 9.07500e-04_jprb, &
     & 1.08425e-03_jprb, 1.29542e-03_jprb, 1.54773e-03_jprb, 1.84917e-03_jprb, 2.20933e-03_jprb, &
     & 2.63963e-03_jprb, 3.15374e-03_jprb, 3.76797e-03_jprb, 4.50184e-03_jprb, 5.37865e-03_jprb, &
     & 6.42622e-03_jprb, 7.67783e-03_jprb, 9.17320e-03_jprb, 1.09598e-02_jprb/)
      kao_mco2( 6, :, 9) = (/ &
     & 2.87301e-04_jprb, 3.43009e-04_jprb, 4.09519e-04_jprb, 4.88926e-04_jprb, 5.83730e-04_jprb, &
     & 6.96916e-04_jprb, 8.32050e-04_jprb, 9.93386e-04_jprb, 1.18601e-03_jprb, 1.41597e-03_jprb, &
     & 1.69053e-03_jprb, 2.01833e-03_jprb, 2.40969e-03_jprb, 2.87693e-03_jprb, 3.43478e-03_jprb, &
     & 4.10079e-03_jprb, 4.89594e-03_jprb, 5.84527e-03_jprb, 6.97867e-03_jprb/)
      kao_mco2( 7, :, 9) = (/ &
     & 1.10743e-04_jprb, 1.32286e-04_jprb, 1.58020e-04_jprb, 1.88760e-04_jprb, 2.25480e-04_jprb, &
     & 2.69342e-04_jprb, 3.21738e-04_jprb, 3.84326e-04_jprb, 4.59090e-04_jprb, 5.48397e-04_jprb, &
     & 6.55078e-04_jprb, 7.82511e-04_jprb, 9.34734e-04_jprb, 1.11657e-03_jprb, 1.33378e-03_jprb, &
     & 1.59324e-03_jprb, 1.90318e-03_jprb, 2.27340e-03_jprb, 2.71565e-03_jprb/)
      kao_mco2( 8, :, 9) = (/ &
     & 8.63177e-05_jprb, 1.03067e-04_jprb, 1.23066e-04_jprb, 1.46946e-04_jprb, 1.75459e-04_jprb, &
     & 2.09505e-04_jprb, 2.50158e-04_jprb, 2.98698e-04_jprb, 3.56658e-04_jprb, 4.25864e-04_jprb, &
     & 5.08498e-04_jprb, 6.07168e-04_jprb, 7.24982e-04_jprb, 8.65658e-04_jprb, 1.03363e-03_jprb, &
     & 1.23420e-03_jprb, 1.47368e-03_jprb, 1.75963e-03_jprb, 2.10107e-03_jprb/)
      kao_mco2( 9, :, 9) = (/ &
     & 4.52715e-04_jprb, 5.41540e-04_jprb, 6.47792e-04_jprb, 7.74892e-04_jprb, 9.26929e-04_jprb, &
     & 1.10880e-03_jprb, 1.32635e-03_jprb, 1.58658e-03_jprb, 1.89787e-03_jprb, 2.27024e-03_jprb, &
     & 2.71568e-03_jprb, 3.24850e-03_jprb, 3.88587e-03_jprb, 4.64830e-03_jprb, 5.56031e-03_jprb, &
     & 6.65127e-03_jprb, 7.95627e-03_jprb, 9.51732e-03_jprb, 1.13847e-02_jprb/)
      kao_mco2( 1, :,10) = (/ &
     & 9.10418e-04_jprb, 1.08631e-03_jprb, 1.29619e-03_jprb, 1.54662e-03_jprb, 1.84543e-03_jprb, &
     & 2.20198e-03_jprb, 2.62741e-03_jprb, 3.13503e-03_jprb, 3.74073e-03_jprb, 4.46344e-03_jprb, &
     & 5.32580e-03_jprb, 6.35476e-03_jprb, 7.58251e-03_jprb, 9.04748e-03_jprb, 1.07955e-02_jprb, &
     & 1.28812e-02_jprb, 1.53699e-02_jprb, 1.83394e-02_jprb, 2.18826e-02_jprb/)
      kao_mco2( 2, :,10) = (/ &
     & 9.06680e-04_jprb, 1.08622e-03_jprb, 1.30130e-03_jprb, 1.55898e-03_jprb, 1.86768e-03_jprb, &
     & 2.23750e-03_jprb, 2.68056e-03_jprb, 3.21135e-03_jprb, 3.84724e-03_jprb, 4.60905e-03_jprb, &
     & 5.52171e-03_jprb, 6.61508e-03_jprb, 7.92496e-03_jprb, 9.49421e-03_jprb, 1.13742e-02_jprb, &
     & 1.36265e-02_jprb, 1.63247e-02_jprb, 1.95572e-02_jprb, 2.34298e-02_jprb/)
      kao_mco2( 3, :,10) = (/ &
     & 8.17976e-04_jprb, 9.79458e-04_jprb, 1.17282e-03_jprb, 1.40435e-03_jprb, 1.68160e-03_jprb, &
     & 2.01357e-03_jprb, 2.41108e-03_jprb, 2.88707e-03_jprb, 3.45703e-03_jprb, 4.13950e-03_jprb, &
     & 4.95671e-03_jprb, 5.93525e-03_jprb, 7.10696e-03_jprb, 8.51000e-03_jprb, 1.01900e-02_jprb, &
     & 1.22017e-02_jprb, 1.46105e-02_jprb, 1.74949e-02_jprb, 2.09486e-02_jprb/)
      kao_mco2( 4, :,10) = (/ &
     & 3.70314e-04_jprb, 4.41440e-04_jprb, 5.26226e-04_jprb, 6.27298e-04_jprb, 7.47782e-04_jprb, &
     & 8.91407e-04_jprb, 1.06262e-03_jprb, 1.26671e-03_jprb, 1.51001e-03_jprb, 1.80003e-03_jprb, &
     & 2.14576e-03_jprb, 2.55789e-03_jprb, 3.04918e-03_jprb, 3.63483e-03_jprb, 4.33297e-03_jprb, &
     & 5.16520e-03_jprb, 6.15727e-03_jprb, 7.33988e-03_jprb, 8.74963e-03_jprb/)
      kao_mco2( 5, :,10) = (/ &
     & 1.00859e-04_jprb, 1.19692e-04_jprb, 1.42041e-04_jprb, 1.68563e-04_jprb, 2.00038e-04_jprb, &
     & 2.37389e-04_jprb, 2.81715e-04_jprb, 3.34318e-04_jprb, 3.96742e-04_jprb, 4.70823e-04_jprb, &
     & 5.58736e-04_jprb, 6.63065e-04_jprb, 7.86874e-04_jprb, 9.33801e-04_jprb, 1.10816e-03_jprb, &
     & 1.31508e-03_jprb, 1.56064e-03_jprb, 1.85204e-03_jprb, 2.19786e-03_jprb/)
      kao_mco2( 6, :,10) = (/ &
     & 9.24477e-05_jprb, 1.09659e-04_jprb, 1.30074e-04_jprb, 1.54290e-04_jprb, 1.83015e-04_jprb, &
     & 2.17087e-04_jprb, 2.57503e-04_jprb, 3.05442e-04_jprb, 3.62307e-04_jprb, 4.29759e-04_jprb, &
     & 5.09768e-04_jprb, 6.04672e-04_jprb, 7.17245e-04_jprb, 8.50776e-04_jprb, 1.00917e-03_jprb, &
     & 1.19704e-03_jprb, 1.41990e-03_jprb, 1.68425e-03_jprb, 1.99780e-03_jprb/)
      kao_mco2( 7, :,10) = (/ &
     & 8.42943e-05_jprb, 1.00044e-04_jprb, 1.18735e-04_jprb, 1.40919e-04_jprb, 1.67248e-04_jprb, &
     & 1.98496e-04_jprb, 2.35582e-04_jprb, 2.79597e-04_jprb, 3.31836e-04_jprb, 3.93835e-04_jprb, &
     & 4.67418e-04_jprb, 5.54748e-04_jprb, 6.58395e-04_jprb, 7.81407e-04_jprb, 9.27402e-04_jprb, &
     & 1.10067e-03_jprb, 1.30632e-03_jprb, 1.55039e-03_jprb, 1.84005e-03_jprb/)
      kao_mco2( 8, :,10) = (/ &
     & 6.86464e-05_jprb, 8.18163e-05_jprb, 9.75129e-05_jprb, 1.16221e-04_jprb, 1.38518e-04_jprb, &
     & 1.65093e-04_jprb, 1.96767e-04_jprb, 2.34517e-04_jprb, 2.79509e-04_jprb, 3.33133e-04_jprb, &
     & 3.97046e-04_jprb, 4.73220e-04_jprb, 5.64008e-04_jprb, 6.72214e-04_jprb, 8.01179e-04_jprb, &
     & 9.54887e-04_jprb, 1.13808e-03_jprb, 1.35643e-03_jprb, 1.61666e-03_jprb/)
      kao_mco2( 9, :,10) = (/ &
     & 1.03095e-04_jprb, 1.21985e-04_jprb, 1.44335e-04_jprb, 1.70781e-04_jprb, 2.02072e-04_jprb, &
     & 2.39096e-04_jprb, 2.82904e-04_jprb, 3.34739e-04_jprb, 3.96070e-04_jprb, 4.68639e-04_jprb, &
     & 5.54505e-04_jprb, 6.56103e-04_jprb, 7.76316e-04_jprb, 9.18556e-04_jprb, 1.08686e-03_jprb, &
     & 1.28599e-03_jprb, 1.52162e-03_jprb, 1.80041e-03_jprb, 2.13029e-03_jprb/)
      kao_mco2( 1, :,11) = (/ &
     & 1.01275e-03_jprb, 1.21433e-03_jprb, 1.45605e-03_jprb, 1.74587e-03_jprb, 2.09339e-03_jprb, &
     & 2.51007e-03_jprb, 3.00970e-03_jprb, 3.60878e-03_jprb, 4.32711e-03_jprb, 5.18842e-03_jprb, &
     & 6.22117e-03_jprb, 7.45950e-03_jprb, 8.94430e-03_jprb, 1.07247e-02_jprb, 1.28594e-02_jprb, &
     & 1.54191e-02_jprb, 1.84882e-02_jprb, 2.21683e-02_jprb, 2.65809e-02_jprb/)
      kao_mco2( 2, :,11) = (/ &
     & 1.06856e-03_jprb, 1.27885e-03_jprb, 1.53052e-03_jprb, 1.83171e-03_jprb, 2.19218e-03_jprb, &
     & 2.62359e-03_jprb, 3.13990e-03_jprb, 3.75781e-03_jprb, 4.49732e-03_jprb, 5.38236e-03_jprb, &
     & 6.44158e-03_jprb, 7.70924e-03_jprb, 9.22637e-03_jprb, 1.10421e-02_jprb, 1.32151e-02_jprb, &
     & 1.58157e-02_jprb, 1.89281e-02_jprb, 2.26531e-02_jprb, 2.71110e-02_jprb/)
      kao_mco2( 3, :,11) = (/ &
     & 7.34896e-04_jprb, 8.77863e-04_jprb, 1.04864e-03_jprb, 1.25265e-03_jprb, 1.49634e-03_jprb, &
     & 1.78744e-03_jprb, 2.13516e-03_jprb, 2.55054e-03_jprb, 3.04672e-03_jprb, 3.63943e-03_jprb, &
     & 4.34745e-03_jprb, 5.19321e-03_jprb, 6.20349e-03_jprb, 7.41032e-03_jprb, 8.85192e-03_jprb, &
     & 1.05740e-02_jprb, 1.26311e-02_jprb, 1.50883e-02_jprb, 1.80236e-02_jprb/)
      kao_mco2( 4, :,11) = (/ &
     & 5.89491e-05_jprb, 7.12560e-05_jprb, 8.61322e-05_jprb, 1.04114e-04_jprb, 1.25850e-04_jprb, &
     & 1.52124e-04_jprb, 1.83883e-04_jprb, 2.22272e-04_jprb, 2.68676e-04_jprb, 3.24768e-04_jprb, &
     & 3.92571e-04_jprb, 4.74528e-04_jprb, 5.73595e-04_jprb, 6.93346e-04_jprb, 8.38096e-04_jprb, &
     & 1.01307e-03_jprb, 1.22457e-03_jprb, 1.48022e-03_jprb, 1.78924e-03_jprb/)
      kao_mco2( 5, :,11) = (/ &
     & 5.32400e-05_jprb, 6.45465e-05_jprb, 7.82542e-05_jprb, 9.48731e-05_jprb, 1.15021e-04_jprb, &
     & 1.39448e-04_jprb, 1.69063e-04_jprb, 2.04966e-04_jprb, 2.48495e-04_jprb, 3.01268e-04_jprb, &
     & 3.65248e-04_jprb, 4.42816e-04_jprb, 5.36856e-04_jprb, 6.50868e-04_jprb, 7.89092e-04_jprb, &
     & 9.56672e-04_jprb, 1.15984e-03_jprb, 1.40615e-03_jprb, 1.70478e-03_jprb/)
      kao_mco2( 6, :,11) = (/ &
     & 5.31408e-05_jprb, 6.42409e-05_jprb, 7.76597e-05_jprb, 9.38814e-05_jprb, 1.13491e-04_jprb, &
     & 1.37198e-04_jprb, 1.65856e-04_jprb, 2.00500e-04_jprb, 2.42381e-04_jprb, 2.93010e-04_jprb, &
     & 3.54214e-04_jprb, 4.28203e-04_jprb, 5.17647e-04_jprb, 6.25774e-04_jprb, 7.56486e-04_jprb, &
     & 9.14503e-04_jprb, 1.10553e-03_jprb, 1.33645e-03_jprb, 1.61561e-03_jprb/)
      kao_mco2( 7, :,11) = (/ &
     & 5.24517e-05_jprb, 6.32485e-05_jprb, 7.62676e-05_jprb, 9.19667e-05_jprb, 1.10897e-04_jprb, &
     & 1.33725e-04_jprb, 1.61251e-04_jprb, 1.94443e-04_jprb, 2.34467e-04_jprb, 2.82730e-04_jprb, &
     & 3.40928e-04_jprb, 4.11106e-04_jprb, 4.95728e-04_jprb, 5.97770e-04_jprb, 7.20816e-04_jprb, &
     & 8.69190e-04_jprb, 1.04811e-03_jprb, 1.26385e-03_jprb, 1.52400e-03_jprb/)
      kao_mco2( 8, :,11) = (/ &
     & 5.01768e-05_jprb, 6.02217e-05_jprb, 7.22774e-05_jprb, 8.67466e-05_jprb, 1.04112e-04_jprb, &
     & 1.24955e-04_jprb, 1.49969e-04_jprb, 1.79991e-04_jprb, 2.16024e-04_jprb, 2.59270e-04_jprb, &
     & 3.11173e-04_jprb, 3.73467e-04_jprb, 4.48231e-04_jprb, 5.37962e-04_jprb, 6.45656e-04_jprb, &
     & 7.74910e-04_jprb, 9.30039e-04_jprb, 1.11622e-03_jprb, 1.33968e-03_jprb/)
      kao_mco2( 9, :,11) = (/ &
     & 5.46391e-05_jprb, 6.58765e-05_jprb, 7.94252e-05_jprb, 9.57603e-05_jprb, 1.15455e-04_jprb, &
     & 1.39200e-04_jprb, 1.67829e-04_jprb, 2.02346e-04_jprb, 2.43962e-04_jprb, 2.94137e-04_jprb, &
     & 3.54632e-04_jprb, 4.27568e-04_jprb, 5.15504e-04_jprb, 6.21526e-04_jprb, 7.49353e-04_jprb, &
     & 9.03471e-04_jprb, 1.08929e-03_jprb, 1.31331e-03_jprb, 1.58342e-03_jprb/)
      kao_mco2( 1, :,12) = (/ &
     & 1.18469e-03_jprb, 1.41755e-03_jprb, 1.69619e-03_jprb, 2.02959e-03_jprb, 2.42854e-03_jprb, &
     & 2.90589e-03_jprb, 3.47708e-03_jprb, 4.16055e-03_jprb, 4.97836e-03_jprb, 5.95691e-03_jprb, &
     & 7.12782e-03_jprb, 8.52889e-03_jprb, 1.02053e-02_jprb, 1.22113e-02_jprb, 1.46116e-02_jprb, &
     & 1.74837e-02_jprb, 2.09204e-02_jprb, 2.50325e-02_jprb, 2.99530e-02_jprb/)
      kao_mco2( 2, :,12) = (/ &
     & 1.09092e-03_jprb, 1.30288e-03_jprb, 1.55602e-03_jprb, 1.85834e-03_jprb, 2.21940e-03_jprb, &
     & 2.65061e-03_jprb, 3.16560e-03_jprb, 3.78064e-03_jprb, 4.51519e-03_jprb, 5.39245e-03_jprb, &
     & 6.44016e-03_jprb, 7.69143e-03_jprb, 9.18580e-03_jprb, 1.09705e-02_jprb, 1.31020e-02_jprb, &
     & 1.56476e-02_jprb, 1.86878e-02_jprb, 2.23187e-02_jprb, 2.66550e-02_jprb/)
      kao_mco2( 3, :,12) = (/ &
     & 3.97521e-04_jprb, 4.74103e-04_jprb, 5.65438e-04_jprb, 6.74369e-04_jprb, 8.04285e-04_jprb, &
     & 9.59228e-04_jprb, 1.14402e-03_jprb, 1.36442e-03_jprb, 1.62727e-03_jprb, 1.94076e-03_jprb, &
     & 2.31464e-03_jprb, 2.76055e-03_jprb, 3.29237e-03_jprb, 3.92663e-03_jprb, 4.68309e-03_jprb, &
     & 5.58528e-03_jprb, 6.66128e-03_jprb, 7.94456e-03_jprb, 9.47505e-03_jprb/)
      kao_mco2( 4, :,12) = (/ &
     & 7.18557e-05_jprb, 8.56230e-05_jprb, 1.02028e-04_jprb, 1.21576e-04_jprb, 1.44870e-04_jprb, &
     & 1.72626e-04_jprb, 2.05701e-04_jprb, 2.45112e-04_jprb, 2.92075e-04_jprb, 3.48035e-04_jprb, &
     & 4.14718e-04_jprb, 4.94176e-04_jprb, 5.88858e-04_jprb, 7.01682e-04_jprb, 8.36121e-04_jprb, &
     & 9.96319e-04_jprb, 1.18721e-03_jprb, 1.41467e-03_jprb, 1.68572e-03_jprb/)
      kao_mco2( 5, :,12) = (/ &
     & 7.33026e-05_jprb, 8.69077e-05_jprb, 1.03038e-04_jprb, 1.22162e-04_jprb, 1.44836e-04_jprb, &
     & 1.71717e-04_jprb, 2.03588e-04_jprb, 2.41375e-04_jprb, 2.86175e-04_jprb, 3.39289e-04_jprb, &
     & 4.02262e-04_jprb, 4.76923e-04_jprb, 5.65440e-04_jprb, 6.70387e-04_jprb, 7.94812e-04_jprb, &
     & 9.42331e-04_jprb, 1.11723e-03_jprb, 1.32459e-03_jprb, 1.57044e-03_jprb/)
      kao_mco2( 6, :,12) = (/ &
     & 7.44053e-05_jprb, 8.82167e-05_jprb, 1.04592e-04_jprb, 1.24007e-04_jprb, 1.47025e-04_jprb, &
     & 1.74317e-04_jprb, 2.06674e-04_jprb, 2.45038e-04_jprb, 2.90523e-04_jprb, 3.44451e-04_jprb, &
     & 4.08389e-04_jprb, 4.84196e-04_jprb, 5.74074e-04_jprb, 6.80637e-04_jprb, 8.06979e-04_jprb, &
     & 9.56774e-04_jprb, 1.13437e-03_jprb, 1.34494e-03_jprb, 1.59459e-03_jprb/)
      kao_mco2( 7, :,12) = (/ &
     & 7.68762e-05_jprb, 9.11305e-05_jprb, 1.08028e-04_jprb, 1.28058e-04_jprb, 1.51802e-04_jprb, &
     & 1.79949e-04_jprb, 2.13315e-04_jprb, 2.52868e-04_jprb, 2.99754e-04_jprb, 3.55334e-04_jprb, &
     & 4.21220e-04_jprb, 4.99322e-04_jprb, 5.91905e-04_jprb, 7.01656e-04_jprb, 8.31756e-04_jprb, &
     & 9.85979e-04_jprb, 1.16880e-03_jprb, 1.38551e-03_jprb, 1.64241e-03_jprb/)
      kao_mco2( 8, :,12) = (/ &
     & 8.45996e-05_jprb, 1.00214e-04_jprb, 1.18711e-04_jprb, 1.40622e-04_jprb, 1.66577e-04_jprb, &
     & 1.97323e-04_jprb, 2.33743e-04_jprb, 2.76885e-04_jprb, 3.27991e-04_jprb, 3.88529e-04_jprb, &
     & 4.60241e-04_jprb, 5.45189e-04_jprb, 6.45816e-04_jprb, 7.65016e-04_jprb, 9.06216e-04_jprb, &
     & 1.07348e-03_jprb, 1.27161e-03_jprb, 1.50632e-03_jprb, 1.78434e-03_jprb/)
      kao_mco2( 9, :,12) = (/ &
     & 7.73583e-05_jprb, 9.16767e-05_jprb, 1.08645e-04_jprb, 1.28755e-04_jprb, 1.52586e-04_jprb, &
     & 1.80829e-04_jprb, 2.14299e-04_jprb, 2.53964e-04_jprb, 3.00970e-04_jprb, 3.56678e-04_jprb, &
     & 4.22696e-04_jprb, 5.00934e-04_jprb, 5.93652e-04_jprb, 7.03533e-04_jprb, 8.33751e-04_jprb, &
     & 9.88072e-04_jprb, 1.17096e-03_jprb, 1.38769e-03_jprb, 1.64454e-03_jprb/)
      kao_mco2( 1, :,13) = (/ &
     & 1.20952e-03_jprb, 1.44504e-03_jprb, 1.72642e-03_jprb, 2.06260e-03_jprb, 2.46423e-03_jprb, &
     & 2.94407e-03_jprb, 3.51735e-03_jprb, 4.20226e-03_jprb, 5.02053e-03_jprb, 5.99814e-03_jprb, &
     & 7.16612e-03_jprb, 8.56153e-03_jprb, 1.02287e-02_jprb, 1.22204e-02_jprb, 1.46000e-02_jprb, &
     & 1.74430e-02_jprb, 2.08395e-02_jprb, 2.48974e-02_jprb, 2.97455e-02_jprb/)
      kao_mco2( 2, :,13) = (/ &
     & 8.47667e-04_jprb, 1.01027e-03_jprb, 1.20407e-03_jprb, 1.43505e-03_jprb, 1.71034e-03_jprb, &
     & 2.03843e-03_jprb, 2.42946e-03_jprb, 2.89550e-03_jprb, 3.45094e-03_jprb, 4.11293e-03_jprb, &
     & 4.90192e-03_jprb, 5.84225e-03_jprb, 6.96296e-03_jprb, 8.29866e-03_jprb, 9.89058e-03_jprb, &
     & 1.17879e-02_jprb, 1.40492e-02_jprb, 1.67442e-02_jprb, 1.99562e-02_jprb/)
      kao_mco2( 3, :,13) = (/ &
     & 1.45612e-04_jprb, 1.71739e-04_jprb, 2.02554e-04_jprb, 2.38897e-04_jprb, 2.81762e-04_jprb, &
     & 3.32318e-04_jprb, 3.91945e-04_jprb, 4.62271e-04_jprb, 5.45215e-04_jprb, 6.43041e-04_jprb, &
     & 7.58421e-04_jprb, 8.94503e-04_jprb, 1.05500e-03_jprb, 1.24430e-03_jprb, 1.46756e-03_jprb, &
     & 1.73088e-03_jprb, 2.04145e-03_jprb, 2.40774e-03_jprb, 2.83975e-03_jprb/)
      kao_mco2( 4, :,13) = (/ &
     & 1.40167e-04_jprb, 1.65266e-04_jprb, 1.94858e-04_jprb, 2.29750e-04_jprb, 2.70889e-04_jprb, &
     & 3.19394e-04_jprb, 3.76585e-04_jprb, 4.44016e-04_jprb, 5.23522e-04_jprb, 6.17264e-04_jprb, &
     & 7.27791e-04_jprb, 8.58110e-04_jprb, 1.01176e-03_jprb, 1.19293e-03_jprb, 1.40654e-03_jprb, &
     & 1.65839e-03_jprb, 1.95534e-03_jprb, 2.30547e-03_jprb, 2.71828e-03_jprb/)
      kao_mco2( 5, :,13) = (/ &
     & 1.37406e-04_jprb, 1.61990e-04_jprb, 1.90973e-04_jprb, 2.25141e-04_jprb, 2.65423e-04_jprb, &
     & 3.12911e-04_jprb, 3.68896e-04_jprb, 4.34898e-04_jprb, 5.12709e-04_jprb, 6.04442e-04_jprb, &
     & 7.12587e-04_jprb, 8.40082e-04_jprb, 9.90387e-04_jprb, 1.16758e-03_jprb, 1.37648e-03_jprb, &
     & 1.62276e-03_jprb, 1.91310e-03_jprb, 2.25539e-03_jprb, 2.65892e-03_jprb/)
      kao_mco2( 6, :,13) = (/ &
     & 1.35356e-04_jprb, 1.59577e-04_jprb, 1.88132e-04_jprb, 2.21797e-04_jprb, 2.61485e-04_jprb, &
     & 3.08276e-04_jprb, 3.63440e-04_jprb, 4.28475e-04_jprb, 5.05147e-04_jprb, 5.95539e-04_jprb, &
     & 7.02106e-04_jprb, 8.27743e-04_jprb, 9.75861e-04_jprb, 1.15048e-03_jprb, 1.35635e-03_jprb, &
     & 1.59906e-03_jprb, 1.88520e-03_jprb, 2.22255e-03_jprb, 2.62025e-03_jprb/)
      kao_mco2( 7, :,13) = (/ &
     & 1.33359e-04_jprb, 1.57252e-04_jprb, 1.85424e-04_jprb, 2.18645e-04_jprb, 2.57817e-04_jprb, &
     & 3.04007e-04_jprb, 3.58472e-04_jprb, 4.22695e-04_jprb, 4.98425e-04_jprb, 5.87722e-04_jprb, &
     & 6.93017e-04_jprb, 8.17177e-04_jprb, 9.63581e-04_jprb, 1.13621e-03_jprb, 1.33978e-03_jprb, &
     & 1.57981e-03_jprb, 1.86284e-03_jprb, 2.19659e-03_jprb, 2.59012e-03_jprb/)
      kao_mco2( 8, :,13) = (/ &
     & 1.29667e-04_jprb, 1.53001e-04_jprb, 1.80534e-04_jprb, 2.13022e-04_jprb, 2.51356e-04_jprb, &
     & 2.96589e-04_jprb, 3.49961e-04_jprb, 4.12938e-04_jprb, 4.87249e-04_jprb, 5.74931e-04_jprb, &
     & 6.78393e-04_jprb, 8.00473e-04_jprb, 9.44521e-04_jprb, 1.11449e-03_jprb, 1.31505e-03_jprb, &
     & 1.55170e-03_jprb, 1.83094e-03_jprb, 2.16042e-03_jprb, 2.54920e-03_jprb/)
      kao_mco2( 9, :,13) = (/ &
     & 1.37892e-04_jprb, 1.62557e-04_jprb, 1.91635e-04_jprb, 2.25914e-04_jprb, 2.66324e-04_jprb, &
     & 3.13963e-04_jprb, 3.70124e-04_jprb, 4.36330e-04_jprb, 5.14379e-04_jprb, 6.06389e-04_jprb, &
     & 7.14858e-04_jprb, 8.42730e-04_jprb, 9.93473e-04_jprb, 1.17118e-03_jprb, 1.38068e-03_jprb, &
     & 1.62765e-03_jprb, 1.91880e-03_jprb, 2.26202e-03_jprb, 2.66665e-03_jprb/)
      kao_mco2( 1, :,14) = (/ &
     & 1.28098e-03_jprb, 1.52939e-03_jprb, 1.82597e-03_jprb, 2.18007e-03_jprb, 2.60284e-03_jprb, &
     & 3.10759e-03_jprb, 3.71022e-03_jprb, 4.42972e-03_jprb, 5.28874e-03_jprb, 6.31435e-03_jprb, &
     & 7.53885e-03_jprb, 9.00081e-03_jprb, 1.07463e-02_jprb, 1.28302e-02_jprb, 1.53183e-02_jprb, &
     & 1.82889e-02_jprb, 2.18355e-02_jprb, 2.60699e-02_jprb, 3.11255e-02_jprb/)
      kao_mco2( 2, :,14) = (/ &
     & 1.27275e-04_jprb, 1.48842e-04_jprb, 1.74064e-04_jprb, 2.03561e-04_jprb, 2.38055e-04_jprb, &
     & 2.78395e-04_jprb, 3.25570e-04_jprb, 3.80740e-04_jprb, 4.45259e-04_jprb, 5.20710e-04_jprb, &
     & 6.08947e-04_jprb, 7.12137e-04_jprb, 8.32812e-04_jprb, 9.73937e-04_jprb, 1.13898e-03_jprb, &
     & 1.33198e-03_jprb, 1.55769e-03_jprb, 1.82165e-03_jprb, 2.13034e-03_jprb/)
      kao_mco2( 3, :,14) = (/ &
     & 1.27744e-04_jprb, 1.49255e-04_jprb, 1.74389e-04_jprb, 2.03755e-04_jprb, 2.38066e-04_jprb, &
     & 2.78155e-04_jprb, 3.24995e-04_jprb, 3.79722e-04_jprb, 4.43666e-04_jprb, 5.18376e-04_jprb, &
     & 6.05668e-04_jprb, 7.07660e-04_jprb, 8.26826e-04_jprb, 9.66059e-04_jprb, 1.12874e-03_jprb, &
     & 1.31881e-03_jprb, 1.54089e-03_jprb, 1.80037e-03_jprb, 2.10354e-03_jprb/)
      kao_mco2( 4, :,14) = (/ &
     & 1.28543e-04_jprb, 1.50136e-04_jprb, 1.75357e-04_jprb, 2.04814e-04_jprb, 2.39219e-04_jprb, &
     & 2.79404e-04_jprb, 3.26339e-04_jprb, 3.81159e-04_jprb, 4.45188e-04_jprb, 5.19972e-04_jprb, &
     & 6.07319e-04_jprb, 7.09339e-04_jprb, 8.28496e-04_jprb, 9.67670e-04_jprb, 1.13022e-03_jprb, &
     & 1.32008e-03_jprb, 1.54184e-03_jprb, 1.80084e-03_jprb, 2.10335e-03_jprb/)
      kao_mco2( 5, :,14) = (/ &
     & 1.29218e-04_jprb, 1.50897e-04_jprb, 1.76214e-04_jprb, 2.05778e-04_jprb, 2.40302e-04_jprb, &
     & 2.80618e-04_jprb, 3.27698e-04_jprb, 3.82678e-04_jprb, 4.46881e-04_jprb, 5.21855e-04_jprb, &
     & 6.09409e-04_jprb, 7.11652e-04_jprb, 8.31048e-04_jprb, 9.70475e-04_jprb, 1.13330e-03_jprb, &
     & 1.32343e-03_jprb, 1.54547e-03_jprb, 1.80476e-03_jprb, 2.10755e-03_jprb/)
      kao_mco2( 6, :,14) = (/ &
     & 1.30502e-04_jprb, 1.52368e-04_jprb, 1.77898e-04_jprb, 2.07706e-04_jprb, 2.42508e-04_jprb, &
     & 2.83141e-04_jprb, 3.30583e-04_jprb, 3.85974e-04_jprb, 4.50646e-04_jprb, 5.26153e-04_jprb, &
     & 6.14313e-04_jprb, 7.17244e-04_jprb, 8.37422e-04_jprb, 9.77736e-04_jprb, 1.14156e-03_jprb, &
     & 1.33283e-03_jprb, 1.55616e-03_jprb, 1.81690e-03_jprb, 2.12133e-03_jprb/)
      kao_mco2( 7, :,14) = (/ &
     & 1.32820e-04_jprb, 1.55041e-04_jprb, 1.80980e-04_jprb, 2.11259e-04_jprb, 2.46604e-04_jprb, &
     & 2.87862e-04_jprb, 3.36022e-04_jprb, 3.92240e-04_jprb, 4.57864e-04_jprb, 5.34467e-04_jprb, &
     & 6.23886e-04_jprb, 7.28265e-04_jprb, 8.50107e-04_jprb, 9.92334e-04_jprb, 1.15836e-03_jprb, &
     & 1.35215e-03_jprb, 1.57838e-03_jprb, 1.84244e-03_jprb, 2.15069e-03_jprb/)
      kao_mco2( 8, :,14) = (/ &
     & 1.40203e-04_jprb, 1.63590e-04_jprb, 1.90879e-04_jprb, 2.22720e-04_jprb, 2.59872e-04_jprb, &
     & 3.03221e-04_jprb, 3.53801e-04_jprb, 4.12819e-04_jprb, 4.81681e-04_jprb, 5.62031e-04_jprb, &
     & 6.55783e-04_jprb, 7.65175e-04_jprb, 8.92814e-04_jprb, 1.04174e-03_jprb, 1.21552e-03_jprb, &
     & 1.41828e-03_jprb, 1.65486e-03_jprb, 1.93091e-03_jprb, 2.25301e-03_jprb/)
      kao_mco2( 9, :,14) = (/ &
     & 1.30642e-04_jprb, 1.52513e-04_jprb, 1.78046e-04_jprb, 2.07853e-04_jprb, 2.42651e-04_jprb, &
     & 2.83275e-04_jprb, 3.30699e-04_jprb, 3.86063e-04_jprb, 4.50696e-04_jprb, 5.26149e-04_jprb, &
     & 6.14234e-04_jprb, 7.17066e-04_jprb, 8.37113e-04_jprb, 9.77259e-04_jprb, 1.14087e-03_jprb, &
     & 1.33186e-03_jprb, 1.55484e-03_jprb, 1.81514e-03_jprb, 2.11902e-03_jprb/)
      kao_mco2( 1, :,15) = (/ &
     & 1.37603e-03_jprb, 1.64035e-03_jprb, 1.95543e-03_jprb, 2.33105e-03_jprb, 2.77881e-03_jprb, &
     & 3.31257e-03_jprb, 3.94887e-03_jprb, 4.70739e-03_jprb, 5.61162e-03_jprb, 6.68952e-03_jprb, &
     & 7.97449e-03_jprb, 9.50627e-03_jprb, 1.13323e-02_jprb, 1.35091e-02_jprb, 1.61039e-02_jprb, &
     & 1.91973e-02_jprb, 2.28848e-02_jprb, 2.72806e-02_jprb, 3.25208e-02_jprb/)
      kao_mco2( 2, :,15) = (/ &
     & 1.67843e-04_jprb, 1.93707e-04_jprb, 2.23557e-04_jprb, 2.58007e-04_jprb, 2.97765e-04_jprb, &
     & 3.43650e-04_jprb, 3.96606e-04_jprb, 4.57722e-04_jprb, 5.28256e-04_jprb, 6.09659e-04_jprb, &
     & 7.03606e-04_jprb, 8.12031e-04_jprb, 9.37163e-04_jprb, 1.08158e-03_jprb, 1.24825e-03_jprb, &
     & 1.44060e-03_jprb, 1.66259e-03_jprb, 1.91880e-03_jprb, 2.21448e-03_jprb/)
      kao_mco2( 3, :,15) = (/ &
     & 1.67595e-04_jprb, 1.93410e-04_jprb, 2.23200e-04_jprb, 2.57579e-04_jprb, 2.97253e-04_jprb, &
     & 3.43039e-04_jprb, 3.95876e-04_jprb, 4.56852e-04_jprb, 5.27220e-04_jprb, 6.08426e-04_jprb, &
     & 7.02141e-04_jprb, 8.10291e-04_jprb, 9.35098e-04_jprb, 1.07913e-03_jprb, 1.24534e-03_jprb, &
     & 1.43716e-03_jprb, 1.65853e-03_jprb, 1.91398e-03_jprb, 2.20879e-03_jprb/)
      kao_mco2( 4, :,15) = (/ &
     & 1.67354e-04_jprb, 1.93130e-04_jprb, 2.22877e-04_jprb, 2.57206e-04_jprb, 2.96823e-04_jprb, &
     & 3.42541e-04_jprb, 3.95301e-04_jprb, 4.56187e-04_jprb, 5.26452e-04_jprb, 6.07539e-04_jprb, &
     & 7.01116e-04_jprb, 8.09106e-04_jprb, 9.33728e-04_jprb, 1.07755e-03_jprb, 1.24352e-03_jprb, &
     & 1.43505e-03_jprb, 1.65608e-03_jprb, 1.91116e-03_jprb, 2.20553e-03_jprb/)
      kao_mco2( 5, :,15) = (/ &
     & 1.67437e-04_jprb, 1.93232e-04_jprb, 2.23002e-04_jprb, 2.57358e-04_jprb, 2.97006e-04_jprb, &
     & 3.42763e-04_jprb, 3.95570e-04_jprb, 4.56511e-04_jprb, 5.26842e-04_jprb, 6.08007e-04_jprb, &
     & 7.01677e-04_jprb, 8.09778e-04_jprb, 9.34533e-04_jprb, 1.07851e-03_jprb, 1.24466e-03_jprb, &
     & 1.43642e-03_jprb, 1.65771e-03_jprb, 1.91310e-03_jprb, 2.20783e-03_jprb/)
      kao_mco2( 6, :,15) = (/ &
     & 1.67267e-04_jprb, 1.93027e-04_jprb, 2.22753e-04_jprb, 2.57057e-04_jprb, 2.96645e-04_jprb, &
     & 3.42328e-04_jprb, 3.95047e-04_jprb, 4.55885e-04_jprb, 5.26092e-04_jprb, 6.07110e-04_jprb, &
     & 7.00606e-04_jprb, 8.08500e-04_jprb, 9.33010e-04_jprb, 1.07669e-03_jprb, 1.24251e-03_jprb, &
     & 1.43385e-03_jprb, 1.65467e-03_jprb, 1.90949e-03_jprb, 2.20355e-03_jprb/)
      kao_mco2( 7, :,15) = (/ &
     & 1.67354e-04_jprb, 1.93130e-04_jprb, 2.22877e-04_jprb, 2.57206e-04_jprb, 2.96823e-04_jprb, &
     & 3.42541e-04_jprb, 3.95301e-04_jprb, 4.56187e-04_jprb, 5.26452e-04_jprb, 6.07539e-04_jprb, &
     & 7.01116e-04_jprb, 8.09106e-04_jprb, 9.33728e-04_jprb, 1.07755e-03_jprb, 1.24352e-03_jprb, &
     & 1.43505e-03_jprb, 1.65608e-03_jprb, 1.91116e-03_jprb, 2.20553e-03_jprb/)
      kao_mco2( 8, :,15) = (/ &
     & 1.67276e-04_jprb, 1.93038e-04_jprb, 2.22769e-04_jprb, 2.57079e-04_jprb, 2.96673e-04_jprb, &
     & 3.42365e-04_jprb, 3.95094e-04_jprb, 4.55944e-04_jprb, 5.26166e-04_jprb, 6.07203e-04_jprb, &
     & 7.00722e-04_jprb, 8.08643e-04_jprb, 9.33186e-04_jprb, 1.07691e-03_jprb, 1.24277e-03_jprb, &
     & 1.43417e-03_jprb, 1.65506e-03_jprb, 1.90996e-03_jprb, 2.20412e-03_jprb/)
      kao_mco2( 9, :,15) = (/ &
     & 1.67437e-04_jprb, 1.93232e-04_jprb, 2.23002e-04_jprb, 2.57358e-04_jprb, 2.97006e-04_jprb, &
     & 3.42763e-04_jprb, 3.95570e-04_jprb, 4.56511e-04_jprb, 5.26842e-04_jprb, 6.08007e-04_jprb, &
     & 7.01677e-04_jprb, 8.09778e-04_jprb, 9.34533e-04_jprb, 1.07851e-03_jprb, 1.24466e-03_jprb, &
     & 1.43642e-03_jprb, 1.65771e-03_jprb, 1.91310e-03_jprb, 2.20783e-03_jprb/)
      kao_mco2( 1, :,16) = (/ &
     & 1.42104e-03_jprb, 1.69791e-03_jprb, 2.02872e-03_jprb, 2.42399e-03_jprb, 2.89626e-03_jprb, &
     & 3.46055e-03_jprb, 4.13478e-03_jprb, 4.94038e-03_jprb, 5.90294e-03_jprb, 7.05303e-03_jprb, &
     & 8.42720e-03_jprb, 1.00691e-02_jprb, 1.20309e-02_jprb, 1.43749e-02_jprb, 1.71757e-02_jprb, &
     & 2.05221e-02_jprb, 2.45205e-02_jprb, 2.92979e-02_jprb, 3.50061e-02_jprb/)
      kao_mco2( 2, :,16) = (/ &
     & 1.63777e-04_jprb, 1.88736e-04_jprb, 2.17498e-04_jprb, 2.50643e-04_jprb, 2.88839e-04_jprb, &
     & 3.32857e-04_jprb, 3.83582e-04_jprb, 4.42037e-04_jprb, 5.09401e-04_jprb, 5.87030e-04_jprb, &
     & 6.76490e-04_jprb, 7.79583e-04_jprb, 8.98386e-04_jprb, 1.03530e-03_jprb, 1.19307e-03_jprb, &
     & 1.37488e-03_jprb, 1.58441e-03_jprb, 1.82586e-03_jprb, 2.10411e-03_jprb/)
      kao_mco2( 3, :,16) = (/ &
     & 1.63679e-04_jprb, 1.88621e-04_jprb, 2.17365e-04_jprb, 2.50489e-04_jprb, 2.88661e-04_jprb, &
     & 3.32650e-04_jprb, 3.83342e-04_jprb, 4.41759e-04_jprb, 5.09079e-04_jprb, 5.86657e-04_jprb, &
     & 6.76057e-04_jprb, 7.79080e-04_jprb, 8.97804e-04_jprb, 1.03462e-03_jprb, 1.19228e-03_jprb, &
     & 1.37397e-03_jprb, 1.58335e-03_jprb, 1.82464e-03_jprb, 2.10269e-03_jprb/)
      kao_mco2( 4, :,16) = (/ &
     & 1.63679e-04_jprb, 1.88621e-04_jprb, 2.17365e-04_jprb, 2.50489e-04_jprb, 2.88661e-04_jprb, &
     & 3.32650e-04_jprb, 3.83342e-04_jprb, 4.41759e-04_jprb, 5.09079e-04_jprb, 5.86657e-04_jprb, &
     & 6.76057e-04_jprb, 7.79080e-04_jprb, 8.97804e-04_jprb, 1.03462e-03_jprb, 1.19228e-03_jprb, &
     & 1.37397e-03_jprb, 1.58335e-03_jprb, 1.82464e-03_jprb, 2.10269e-03_jprb/)
      kao_mco2( 5, :,16) = (/ &
     & 1.63586e-04_jprb, 1.88513e-04_jprb, 2.17239e-04_jprb, 2.50343e-04_jprb, 2.88490e-04_jprb, &
     & 3.32451e-04_jprb, 3.83111e-04_jprb, 4.41490e-04_jprb, 5.08766e-04_jprb, 5.86292e-04_jprb, &
     & 6.75633e-04_jprb, 7.78588e-04_jprb, 8.97231e-04_jprb, 1.03395e-03_jprb, 1.19151e-03_jprb, &
     & 1.37307e-03_jprb, 1.58231e-03_jprb, 1.82342e-03_jprb, 2.10128e-03_jprb/)
      kao_mco2( 6, :,16) = (/ &
     & 1.63679e-04_jprb, 1.88621e-04_jprb, 2.17365e-04_jprb, 2.50489e-04_jprb, 2.88661e-04_jprb, &
     & 3.32650e-04_jprb, 3.83342e-04_jprb, 4.41759e-04_jprb, 5.09079e-04_jprb, 5.86657e-04_jprb, &
     & 6.76057e-04_jprb, 7.79080e-04_jprb, 8.97804e-04_jprb, 1.03462e-03_jprb, 1.19228e-03_jprb, &
     & 1.37397e-03_jprb, 1.58335e-03_jprb, 1.82464e-03_jprb, 2.10269e-03_jprb/)
      kao_mco2( 7, :,16) = (/ &
     & 1.63679e-04_jprb, 1.88621e-04_jprb, 2.17365e-04_jprb, 2.50489e-04_jprb, 2.88661e-04_jprb, &
     & 3.32650e-04_jprb, 3.83342e-04_jprb, 4.41759e-04_jprb, 5.09079e-04_jprb, 5.86657e-04_jprb, &
     & 6.76057e-04_jprb, 7.79080e-04_jprb, 8.97804e-04_jprb, 1.03462e-03_jprb, 1.19228e-03_jprb, &
     & 1.37397e-03_jprb, 1.58335e-03_jprb, 1.82464e-03_jprb, 2.10269e-03_jprb/)
      kao_mco2( 8, :,16) = (/ &
     & 1.63479e-04_jprb, 1.88391e-04_jprb, 2.17098e-04_jprb, 2.50180e-04_jprb, 2.88303e-04_jprb, &
     & 3.32236e-04_jprb, 3.82863e-04_jprb, 4.41205e-04_jprb, 5.08437e-04_jprb, 5.85914e-04_jprb, &
     & 6.75198e-04_jprb, 7.78087e-04_jprb, 8.96654e-04_jprb, 1.03329e-03_jprb, 1.19074e-03_jprb, &
     & 1.37219e-03_jprb, 1.58129e-03_jprb, 1.82226e-03_jprb, 2.09994e-03_jprb/)
      kao_mco2( 9, :,16) = (/ &
     & 1.63586e-04_jprb, 1.88513e-04_jprb, 2.17239e-04_jprb, 2.50343e-04_jprb, 2.88490e-04_jprb, &
     & 3.32451e-04_jprb, 3.83111e-04_jprb, 4.41490e-04_jprb, 5.08766e-04_jprb, 5.86292e-04_jprb, &
     & 6.75633e-04_jprb, 7.78588e-04_jprb, 8.97231e-04_jprb, 1.03395e-03_jprb, 1.19151e-03_jprb, &
     & 1.37307e-03_jprb, 1.58231e-03_jprb, 1.82342e-03_jprb, 2.10128e-03_jprb/)

!     the array kbo_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kbo_mco2(:, 1) = (/ &
     & 3.72069e-06_jprb, 4.81866e-06_jprb, 6.24064e-06_jprb, 8.08226e-06_jprb, 1.04673e-05_jprb, &
     & 1.35562e-05_jprb, 1.75567e-05_jprb, 2.27376e-05_jprb, 2.94475e-05_jprb, 3.81375e-05_jprb, &
     & 4.93918e-05_jprb, 6.39674e-05_jprb, 8.28441e-05_jprb, 1.07291e-04_jprb, 1.38953e-04_jprb, &
     & 1.79958e-04_jprb, 2.33064e-04_jprb, 3.01840e-04_jprb, 3.90913e-04_jprb/)
      kbo_mco2(:, 2) = (/ &
     & 8.14357e-06_jprb, 1.06031e-05_jprb, 1.38056e-05_jprb, 1.79752e-05_jprb, 2.34041e-05_jprb, &
     & 3.04728e-05_jprb, 3.96763e-05_jprb, 5.16596e-05_jprb, 6.72622e-05_jprb, 8.75770e-05_jprb, &
     & 1.14027e-04_jprb, 1.48467e-04_jprb, 1.93307e-04_jprb, 2.51691e-04_jprb, 3.27708e-04_jprb, &
     & 4.26685e-04_jprb, 5.55555e-04_jprb, 7.23346e-04_jprb, 9.41814e-04_jprb/)
      kbo_mco2(:, 3) = (/ &
     & 1.09367e-05_jprb, 1.42063e-05_jprb, 1.84533e-05_jprb, 2.39701e-05_jprb, 3.11362e-05_jprb, &
     & 4.04446e-05_jprb, 5.25358e-05_jprb, 6.82417e-05_jprb, 8.86432e-05_jprb, 1.15144e-04_jprb, &
     & 1.49567e-04_jprb, 1.94281e-04_jprb, 2.52363e-04_jprb, 3.27809e-04_jprb, 4.25810e-04_jprb, &
     & 5.53109e-04_jprb, 7.18466e-04_jprb, 9.33256e-04_jprb, 1.21226e-03_jprb/)
      kbo_mco2(:, 4) = (/ &
     & 1.76192e-05_jprb, 2.27752e-05_jprb, 2.94401e-05_jprb, 3.80553e-05_jprb, 4.91916e-05_jprb, &
     & 6.35867e-05_jprb, 8.21944e-05_jprb, 1.06247e-04_jprb, 1.37339e-04_jprb, 1.77529e-04_jprb, &
     & 2.29480e-04_jprb, 2.96635e-04_jprb, 3.83440e-04_jprb, 4.95648e-04_jprb, 6.40691e-04_jprb, &
     & 8.28180e-04_jprb, 1.07054e-03_jprb, 1.38381e-03_jprb, 1.78876e-03_jprb/)
      kbo_mco2(:, 5) = (/ &
     & 3.72142e-05_jprb, 4.78603e-05_jprb, 6.15520e-05_jprb, 7.91605e-05_jprb, 1.01806e-04_jprb, &
     & 1.30931e-04_jprb, 1.68387e-04_jprb, 2.16558e-04_jprb, 2.78510e-04_jprb, 3.58185e-04_jprb, &
     & 4.60653e-04_jprb, 5.92435e-04_jprb, 7.61915e-04_jprb, 9.79881e-04_jprb, 1.26020e-03_jprb, &
     & 1.62071e-03_jprb, 2.08436e-03_jprb, 2.68064e-03_jprb, 3.44751e-03_jprb/)
      kbo_mco2(:, 6) = (/ &
     & 7.74131e-05_jprb, 9.98876e-05_jprb, 1.28887e-04_jprb, 1.66305e-04_jprb, 2.14587e-04_jprb, &
     & 2.76886e-04_jprb, 3.57272e-04_jprb, 4.60994e-04_jprb, 5.94831e-04_jprb, 7.67521e-04_jprb, &
     & 9.90348e-04_jprb, 1.27787e-03_jprb, 1.64886e-03_jprb, 2.12755e-03_jprb, 2.74522e-03_jprb, &
     & 3.54221e-03_jprb, 4.57059e-03_jprb, 5.89752e-03_jprb, 7.60968e-03_jprb/)
      kbo_mco2(:, 7) = (/ &
     & 1.32294e-04_jprb, 1.70977e-04_jprb, 2.20973e-04_jprb, 2.85587e-04_jprb, 3.69095e-04_jprb, &
     & 4.77022e-04_jprb, 6.16507e-04_jprb, 7.96779e-04_jprb, 1.02976e-03_jprb, 1.33088e-03_jprb, &
     & 1.72004e-03_jprb, 2.22299e-03_jprb, 2.87301e-03_jprb, 3.71310e-03_jprb, 4.79884e-03_jprb, &
     & 6.20207e-03_jprb, 8.01561e-03_jprb, 1.03594e-02_jprb, 1.33886e-02_jprb/)
      kbo_mco2(:, 8) = (/ &
     & 3.59868e-05_jprb, 4.63611e-05_jprb, 5.97261e-05_jprb, 7.69439e-05_jprb, 9.91253e-05_jprb, &
     & 1.27701e-04_jprb, 1.64515e-04_jprb, 2.11941e-04_jprb, 2.73040e-04_jprb, 3.51752e-04_jprb, &
     & 4.53155e-04_jprb, 5.83790e-04_jprb, 7.52085e-04_jprb, 9.68897e-04_jprb, 1.24821e-03_jprb, &
     & 1.60804e-03_jprb, 2.07161e-03_jprb, 2.66882e-03_jprb, 3.43818e-03_jprb/)
      kbo_mco2(:, 9) = (/ &
     & 5.09543e-05_jprb, 6.60510e-05_jprb, 8.56205e-05_jprb, 1.10988e-04_jprb, 1.43872e-04_jprb, &
     & 1.86498e-04_jprb, 2.41753e-04_jprb, 3.13380e-04_jprb, 4.06228e-04_jprb, 5.26585e-04_jprb, &
     & 6.82601e-04_jprb, 8.84842e-04_jprb, 1.14700e-03_jprb, 1.48684e-03_jprb, 1.92735e-03_jprb, &
     & 2.49839e-03_jprb, 3.23861e-03_jprb, 4.19814e-03_jprb, 5.44196e-03_jprb/)
      kbo_mco2(:,10) = (/ &
     & 2.08253e-05_jprb, 2.64900e-05_jprb, 3.36954e-05_jprb, 4.28609e-05_jprb, 5.45194e-05_jprb, &
     & 6.93491e-05_jprb, 8.82125e-05_jprb, 1.12207e-04_jprb, 1.42728e-04_jprb, 1.81551e-04_jprb, &
     & 2.30935e-04_jprb, 2.93751e-04_jprb, 3.73653e-04_jprb, 4.75290e-04_jprb, 6.04572e-04_jprb, &
     & 7.69021e-04_jprb, 9.78201e-04_jprb, 1.24428e-03_jprb, 1.58273e-03_jprb/)
      kbo_mco2(:,11) = (/ &
     & 2.08953e-05_jprb, 2.65543e-05_jprb, 3.37459e-05_jprb, 4.28852e-05_jprb, 5.44996e-05_jprb, &
     & 6.92595e-05_jprb, 8.80169e-05_jprb, 1.11854e-04_jprb, 1.42147e-04_jprb, 1.80644e-04_jprb, &
     & 2.29568e-04_jprb, 2.91741e-04_jprb, 3.70752e-04_jprb, 4.71161e-04_jprb, 5.98764e-04_jprb, &
     & 7.60925e-04_jprb, 9.67005e-04_jprb, 1.22889e-03_jprb, 1.56171e-03_jprb/)
      kbo_mco2(:,12) = (/ &
     & 2.65295e-05_jprb, 3.36318e-05_jprb, 4.26356e-05_jprb, 5.40498e-05_jprb, 6.85198e-05_jprb, &
     & 8.68636e-05_jprb, 1.10118e-04_jprb, 1.39599e-04_jprb, 1.76972e-04_jprb, 2.24350e-04_jprb, &
     & 2.84412e-04_jprb, 3.60553e-04_jprb, 4.57079e-04_jprb, 5.79446e-04_jprb, 7.34572e-04_jprb, &
     & 9.31230e-04_jprb, 1.18053e-03_jprb, 1.49658e-03_jprb, 1.89724e-03_jprb/)
      kbo_mco2(:,13) = (/ &
     & 3.45358e-05_jprb, 4.36743e-05_jprb, 5.52309e-05_jprb, 6.98455e-05_jprb, 8.83273e-05_jprb, &
     & 1.11700e-04_jprb, 1.41256e-04_jprb, 1.78634e-04_jprb, 2.25902e-04_jprb, 2.85678e-04_jprb, &
     & 3.61271e-04_jprb, 4.56867e-04_jprb, 5.77758e-04_jprb, 7.30639e-04_jprb, 9.23973e-04_jprb, &
     & 1.16847e-03_jprb, 1.47765e-03_jprb, 1.86865e-03_jprb, 2.36311e-03_jprb/)
      kbo_mco2(:,14) = (/ &
     & 3.99721e-05_jprb, 5.12343e-05_jprb, 6.56698e-05_jprb, 8.41725e-05_jprb, 1.07888e-04_jprb, &
     & 1.38286e-04_jprb, 1.77249e-04_jprb, 2.27190e-04_jprb, 2.91201e-04_jprb, 3.73248e-04_jprb, &
     & 4.78412e-04_jprb, 6.13207e-04_jprb, 7.85980e-04_jprb, 1.00743e-03_jprb, 1.29128e-03_jprb, &
     & 1.65510e-03_jprb, 2.12144e-03_jprb, 2.71916e-03_jprb, 3.48529e-03_jprb/)
      kbo_mco2(:,15) = (/ &
     & 8.51533e-06_jprb, 1.23021e-05_jprb, 1.77730e-05_jprb, 2.56767e-05_jprb, 3.70953e-05_jprb, &
     & 5.35918e-05_jprb, 7.74243e-05_jprb, 1.11855e-04_jprb, 1.61598e-04_jprb, 2.33461e-04_jprb, &
     & 3.37283e-04_jprb, 4.87275e-04_jprb, 7.03968e-04_jprb, 1.01703e-03_jprb, 1.46930e-03_jprb, &
     & 2.12271e-03_jprb, 3.06670e-03_jprb, 4.43047e-03_jprb, 6.40072e-03_jprb/)
      kbo_mco2(:,16) = (/ &
     & 2.93050e-06_jprb, 3.65298e-06_jprb, 4.55358e-06_jprb, 5.67622e-06_jprb, 7.07564e-06_jprb, &
     & 8.82006e-06_jprb, 1.09945e-05_jprb, 1.37051e-05_jprb, 1.70840e-05_jprb, 2.12959e-05_jprb, &
     & 2.65461e-05_jprb, 3.30908e-05_jprb, 4.12490e-05_jprb, 5.14185e-05_jprb, 6.40952e-05_jprb, &
     & 7.98972e-05_jprb, 9.95951e-05_jprb, 1.24149e-04_jprb, 1.54757e-04_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296_rb,260_rb,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &2.0677e-07_jprb,2.0363e-07_jprb,2.0583e-07_jprb,2.0547e-07_jprb,2.0267e-07_jprb,2.0154e-07_jprb, &
     &2.0190e-07_jprb,2.0103e-07_jprb,1.9869e-07_jprb,1.9663e-07_jprb,1.9701e-07_jprb,2.0103e-07_jprb, &
     &2.0527e-07_jprb,2.0206e-07_jprb,2.0364e-07_jprb,2.0364e-07_jprb/)
      forrefo(2,:) = (/ &
     &2.2427e-07_jprb,2.1489e-07_jprb,2.0453e-07_jprb,1.9710e-07_jprb,1.9650e-07_jprb,1.9738e-07_jprb, &
     &1.9767e-07_jprb,1.9769e-07_jprb,1.9940e-07_jprb,1.9846e-07_jprb,1.9898e-07_jprb,1.9853e-07_jprb, &
     &2.0000e-07_jprb,2.0517e-07_jprb,2.0482e-07_jprb,2.0482e-07_jprb/)
      forrefo(3,:) = (/ &
     &2.2672e-07_jprb,2.1706e-07_jprb,2.0571e-07_jprb,1.9747e-07_jprb,1.9706e-07_jprb,1.9698e-07_jprb, &
     &1.9781e-07_jprb,1.9774e-07_jprb,1.9724e-07_jprb,1.9714e-07_jprb,1.9751e-07_jprb,1.9758e-07_jprb, &
     &1.9840e-07_jprb,1.9968e-07_jprb,1.9931e-07_jprb,1.9880e-07_jprb/)
      forrefo(4,:) = (/ &
     &2.2191e-07_jprb,2.0899e-07_jprb,2.0265e-07_jprb,2.0101e-07_jprb,2.0034e-07_jprb,2.0021e-07_jprb, &
     &1.9987e-07_jprb,1.9978e-07_jprb,1.9902e-07_jprb,1.9742e-07_jprb,1.9672e-07_jprb,1.9615e-07_jprb, &
     &1.9576e-07_jprb,1.9540e-07_jprb,1.9588e-07_jprb,1.9590e-07_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 5.18832e-02_jprb, 4.28690e-02_jprb, 3.54210e-02_jprb, 2.92670e-02_jprb, 2.41822e-02_jprb, &
     & 1.99808e-02_jprb, 1.65093e-02_jprb, 1.36410e-02_jprb, 1.12710e-02_jprb, 9.31280e-03_jprb/)
      selfrefo(:, 2) = (/ &
     & 4.36030e-02_jprb, 3.78379e-02_jprb, 3.28350e-02_jprb, 2.84936e-02_jprb, 2.47262e-02_jprb, &
     & 2.14569e-02_jprb, 1.86199e-02_jprb, 1.61580e-02_jprb, 1.40216e-02_jprb, 1.21677e-02_jprb/)
      selfrefo(:, 3) = (/ &
     & 4.26492e-02_jprb, 3.71443e-02_jprb, 3.23500e-02_jprb, 2.81745e-02_jprb, 2.45379e-02_jprb, &
     & 2.13707e-02_jprb, 1.86124e-02_jprb, 1.62100e-02_jprb, 1.41177e-02_jprb, 1.22955e-02_jprb/)
      selfrefo(:, 4) = (/ &
     & 4.03591e-02_jprb, 3.54614e-02_jprb, 3.11580e-02_jprb, 2.73769e-02_jprb, 2.40546e-02_jprb, &
     & 2.11355e-02_jprb, 1.85706e-02_jprb, 1.63170e-02_jprb, 1.43369e-02_jprb, 1.25970e-02_jprb/)
      selfrefo(:, 5) = (/ &
     & 3.94512e-02_jprb, 3.46232e-02_jprb, 3.03860e-02_jprb, 2.66674e-02_jprb, 2.34038e-02_jprb, &
     & 2.05397e-02_jprb, 1.80260e-02_jprb, 1.58200e-02_jprb, 1.38839e-02_jprb, 1.21848e-02_jprb/)
      selfrefo(:, 6) = (/ &
     & 3.90567e-02_jprb, 3.40694e-02_jprb, 2.97190e-02_jprb, 2.59241e-02_jprb, 2.26138e-02_jprb, &
     & 1.97261e-02_jprb, 1.72072e-02_jprb, 1.50100e-02_jprb, 1.30933e-02_jprb, 1.14214e-02_jprb/)
      selfrefo(:, 7) = (/ &
     & 3.85397e-02_jprb, 3.36462e-02_jprb, 2.93740e-02_jprb, 2.56443e-02_jprb, 2.23881e-02_jprb, &
     & 1.95454e-02_jprb, 1.70636e-02_jprb, 1.48970e-02_jprb, 1.30055e-02_jprb, 1.13541e-02_jprb/)
      selfrefo(:, 8) = (/ &
     & 3.79692e-02_jprb, 3.31360e-02_jprb, 2.89180e-02_jprb, 2.52369e-02_jprb, 2.20245e-02_jprb, &
     & 1.92209e-02_jprb, 1.67742e-02_jprb, 1.46390e-02_jprb, 1.27756e-02_jprb, 1.11493e-02_jprb/)
      selfrefo(:, 9) = (/ &
     & 3.68819e-02_jprb, 3.22827e-02_jprb, 2.82570e-02_jprb, 2.47333e-02_jprb, 2.16490e-02_jprb, &
     & 1.89494e-02_jprb, 1.65863e-02_jprb, 1.45180e-02_jprb, 1.27076e-02_jprb, 1.11229e-02_jprb/)
      selfrefo(:,10) = (/ &
     & 3.65157e-02_jprb, 3.20121e-02_jprb, 2.80640e-02_jprb, 2.46028e-02_jprb, 2.15685e-02_jprb, &
     & 1.89084e-02_jprb, 1.65764e-02_jprb, 1.45320e-02_jprb, 1.27397e-02_jprb, 1.11685e-02_jprb/)
      selfrefo(:,11) = (/ &
     & 3.59917e-02_jprb, 3.16727e-02_jprb, 2.78720e-02_jprb, 2.45274e-02_jprb, 2.15841e-02_jprb, &
     & 1.89940e-02_jprb, 1.67148e-02_jprb, 1.47090e-02_jprb, 1.29439e-02_jprb, 1.13907e-02_jprb/)
      selfrefo(:,12) = (/ &
     & 3.66963e-02_jprb, 3.20483e-02_jprb, 2.79890e-02_jprb, 2.44439e-02_jprb, 2.13478e-02_jprb, &
     & 1.86438e-02_jprb, 1.62824e-02_jprb, 1.42200e-02_jprb, 1.24189e-02_jprb, 1.08459e-02_jprb/)
      selfrefo(:,13) = (/ &
     & 3.66422e-02_jprb, 3.19026e-02_jprb, 2.77760e-02_jprb, 2.41832e-02_jprb, 2.10551e-02_jprb, &
     & 1.83317e-02_jprb, 1.59605e-02_jprb, 1.38960e-02_jprb, 1.20986e-02_jprb, 1.05336e-02_jprb/)
      selfrefo(:,14) = (/ &
     & 3.81260e-02_jprb, 3.29322e-02_jprb, 2.84460e-02_jprb, 2.45709e-02_jprb, 2.12237e-02_jprb, &
     & 1.83325e-02_jprb, 1.58352e-02_jprb, 1.36780e-02_jprb, 1.18147e-02_jprb, 1.02052e-02_jprb/)
      selfrefo(:,15) = (/ &
     & 3.51264e-02_jprb, 3.05081e-02_jprb, 2.64970e-02_jprb, 2.30133e-02_jprb, 1.99876e-02_jprb, &
     & 1.73597e-02_jprb, 1.50773e-02_jprb, 1.30950e-02_jprb, 1.13733e-02_jprb, 9.87800e-03_jprb/)
      selfrefo(:,16) = (/ &
     & 3.51264e-02_jprb, 3.05081e-02_jprb, 2.64970e-02_jprb, 2.30133e-02_jprb, 1.99876e-02_jprb, &
     & 1.73597e-02_jprb, 1.50773e-02_jprb, 1.30950e-02_jprb, 1.13733e-02_jprb, 9.87800e-03_jprb/)


if (lhook) call dr_hook('rrtm_kgb7',1,zhook_handle)
return

1001 continue

call abor1("rrtm_kgb7:error reading file radrrtm")

end subroutine rrtm_kgb7
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

