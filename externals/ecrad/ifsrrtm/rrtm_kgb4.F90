! # 1 "ifsrrtm/rrtm_kgb4.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb4.f90"
subroutine rrtm_kgb4

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
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
use yommp0    , only : nproc, myproc

use yoerrto4 , only :  kao     ,kbo     ,selfrefo   ,forrefo, fracrefao  ,fracrefbo, &
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
! # 30 "ifsrrtm/rrtm_kgb4.f90" 2

if (lhook) call dr_hook('rrtm_kgb4',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001)  kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb4:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb4:')
endif


! planck fraction mapping level : p = 142.5940 mbar, t = 215.70 k
      fracrefao(:, 1) = (/ &
     &   1.5572e-01_jprb,1.4925e-01_jprb,1.4107e-01_jprb,1.3126e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2393e-02_jprb,4.2146e-02_jprb,4.5907e-03_jprb,3.7965e-03_jprb,2.9744e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 2) = (/ &
     &   1.5572e-01_jprb,1.4925e-01_jprb,1.4107e-01_jprb,1.3126e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2392e-02_jprb,4.2146e-02_jprb,4.5906e-03_jprb,3.7965e-03_jprb,2.9745e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 3) = (/ &
     &   1.5572e-01_jprb,1.4925e-01_jprb,1.4107e-01_jprb,1.3126e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2393e-02_jprb,4.2146e-02_jprb,4.5907e-03_jprb,3.7965e-03_jprb,2.9745e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 4) = (/ &
     &   1.5572e-01_jprb,1.4925e-01_jprb,1.4107e-01_jprb,1.3126e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2393e-02_jprb,4.2146e-02_jprb,4.5907e-03_jprb,3.7964e-03_jprb,2.9744e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 5) = (/ &
     &   1.5572e-01_jprb,1.4925e-01_jprb,1.4107e-01_jprb,1.3126e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2393e-02_jprb,4.2146e-02_jprb,4.5907e-03_jprb,3.7965e-03_jprb,2.9744e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 6) = (/ &
     &   1.5572e-01_jprb,1.4925e-01_jprb,1.4107e-01_jprb,1.3126e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2393e-02_jprb,4.2146e-02_jprb,4.5907e-03_jprb,3.7965e-03_jprb,2.9744e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 7) = (/ &
     &   1.5572e-01_jprb,1.4926e-01_jprb,1.4107e-01_jprb,1.3126e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2393e-02_jprb,4.2146e-02_jprb,4.5908e-03_jprb,3.7964e-03_jprb,2.9745e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 8) = (/ &
     &   1.5571e-01_jprb,1.4926e-01_jprb,1.4107e-01_jprb,1.3125e-01_jprb,1.1791e-01_jprb,1.0173e-01_jprb, &
     &   8.2949e-02_jprb,6.2393e-02_jprb,4.2146e-02_jprb,4.5907e-03_jprb,3.7964e-03_jprb,2.9744e-03_jprb, &
     &   2.2074e-03_jprb,1.4063e-03_jprb,5.3012e-04_jprb,7.4595e-05_jprb/)
      fracrefao(:, 9) = (/ &
     &   1.5952e-01_jprb,1.5155e-01_jprb,1.4217e-01_jprb,1.3077e-01_jprb,1.1667e-01_jprb,1.0048e-01_jprb, &
     &   8.1511e-02_jprb,6.1076e-02_jprb,4.1111e-02_jprb,4.4432e-03_jprb,3.6910e-03_jprb,2.9076e-03_jprb, &
     &   2.1329e-03_jprb,1.3566e-03_jprb,5.2235e-04_jprb,7.9935e-05_jprb/)

! planck fraction mapping level : p = 95.58350 mb, t = 215.70 k
      fracrefbo(:, 1) = (/ &
     &   1.5558e-01_jprb,1.4931e-01_jprb,1.4104e-01_jprb,1.3124e-01_jprb,1.1793e-01_jprb,1.0160e-01_jprb, &
     &   8.3142e-02_jprb,6.2403e-02_jprb,4.2170e-02_jprb,4.5935e-03_jprb,3.7976e-03_jprb,2.9986e-03_jprb, &
     &   2.1890e-03_jprb,1.4061e-03_jprb,5.3005e-04_jprb,7.4587e-05_jprb/)
      fracrefbo(:, 2) = (/ &
     &   1.5558e-01_jprb,1.4932e-01_jprb,1.4104e-01_jprb,1.3124e-01_jprb,1.1792e-01_jprb,1.0159e-01_jprb, &
     &   8.3142e-02_jprb,6.2403e-02_jprb,4.2170e-02_jprb,4.5935e-03_jprb,3.7976e-03_jprb,2.9986e-03_jprb, &
     &   2.1890e-03_jprb,1.4061e-03_jprb,5.3005e-04_jprb,7.4587e-05_jprb/)
      fracrefbo(:, 3) = (/ &
     &   1.5558e-01_jprb,1.4933e-01_jprb,1.4103e-01_jprb,1.3124e-01_jprb,1.1792e-01_jprb,1.0159e-01_jprb, &
     &   8.3142e-02_jprb,6.2403e-02_jprb,4.2170e-02_jprb,4.5935e-03_jprb,3.7976e-03_jprb,2.9986e-03_jprb, &
     &   2.1890e-03_jprb,1.4061e-03_jprb,5.3005e-04_jprb,7.4587e-05_jprb/)
      fracrefbo(:, 4) = (/ &
     &   1.5569e-01_jprb,1.4926e-01_jprb,1.4102e-01_jprb,1.3122e-01_jprb,1.1791e-01_jprb,1.0159e-01_jprb, &
     &   8.3141e-02_jprb,6.2403e-02_jprb,4.2170e-02_jprb,4.5935e-03_jprb,3.7976e-03_jprb,2.9986e-03_jprb, &
     &   2.1890e-03_jprb,1.4061e-03_jprb,5.3005e-04_jprb,7.4587e-05_jprb/)
      fracrefbo(:, 5) = (/ &
     &   1.5947e-01_jprb,1.5132e-01_jprb,1.4195e-01_jprb,1.3061e-01_jprb,1.1680e-01_jprb,1.0054e-01_jprb, &
     &   8.1785e-02_jprb,6.1212e-02_jprb,4.1276e-02_jprb,4.4424e-03_jprb,3.6628e-03_jprb,2.8943e-03_jprb, &
     &   2.1134e-03_jprb,1.3457e-03_jprb,5.1024e-04_jprb,7.3998e-05_jprb/)

!     ------------------------------------------------------------------

!     the array kao contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to co2.  the first index in the array, js, runs
!     from 1 to 9 and corresponds to different water vapor to co2 ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = h2o/(h20 + (rat) * co2), where rat is the ratio of the integrated
!     line strength in the band of co2 to that of h2o.  for instance,
!     js=1 refers to dry air (eta = 0), js = 9 corresponds to eta = 1.0.
!     the 2nd index in the array, jt, which runs from 1 to 5, corresponds 
!     to different temperatures.  more specifically, jt = 3 means that the 
!     data are for the reference temperature tref for this pressure 
!     level, jt = 2 refers to the temperature tref-15,
!     jt = 1 is for tref-30, jt = 4 is for tref+15, and jt = 5
!     is for tref+30.  the third index, jp, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. jp = 1 is for a
!     pressure of 1053.63 mb).  the fourth index, ig, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.



!     the array kbo contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels  < ~100mb, temperatures, and ratios
!     of o3 to co2.  the first index in the array, js, runs from 1 to 6, 
!     and corresponds to different o3 to co2 ratios, as expressed through 
!     the binary species parameter eta, defined as eta = o3/(o3+rat*h2o), 
!     where rat is the ratio of the integrated line strength in the band 
!     of co2 to that of o3.  for instance, js=1 refers to no o3 (eta = 0) 
!     and js = 5 corresponds to eta = 1.0.  the second index, jt, which
!     runs from 1 to 5, corresponds to different temperatures.  more 
!     specifically, jt = 3 means that the data are for the corresponding 
!     reference temperature tref for this  pressure level, jt = 2 refers 
!     to the tref-15, jt = 1 is for tref-30, jt = 4 is for tref+15, and
!     jt = 5 is for tref+30.  the third index, jp, runs from 13 to 59 and
!     refers to the corresponding pressure level in pref (e.g. jp = 13 is
!     for a pressure of 95.5835 mb).  the fourth index, ig, goes from 1 to
!     16, and tells us which g-interval the absorption coefficients are for.



!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &3.3839e-04_jprb,2.4739e-04_jprb,2.2846e-04_jprb,2.3376e-04_jprb,2.2622e-04_jprb,2.3188e-04_jprb, &
     &2.2990e-04_jprb,2.2532e-04_jprb,2.1233e-04_jprb,2.0593e-04_jprb,2.0716e-04_jprb,2.0809e-04_jprb, &
     &2.0889e-04_jprb,2.0932e-04_jprb,2.0944e-04_jprb,2.0945e-04_jprb/)
      forrefo(2,:) = (/ &
     &3.4391e-04_jprb,2.6022e-04_jprb,2.3449e-04_jprb,2.4544e-04_jprb,2.3831e-04_jprb,2.3014e-04_jprb, &
     &2.3729e-04_jprb,2.2726e-04_jprb,2.1892e-04_jprb,1.9223e-04_jprb,2.1291e-04_jprb,2.1406e-04_jprb, &
     &2.1491e-04_jprb,2.1548e-04_jprb,2.1562e-04_jprb,2.1567e-04_jprb/)
      forrefo(3,:) = (/ &
     &3.4219e-04_jprb,2.7334e-04_jprb,2.3727e-04_jprb,2.4515e-04_jprb,2.5272e-04_jprb,2.4212e-04_jprb, &
     &2.3824e-04_jprb,2.3615e-04_jprb,2.2724e-04_jprb,2.2381e-04_jprb,1.9634e-04_jprb,2.1625e-04_jprb, &
     &2.1963e-04_jprb,2.2032e-04_jprb,2.2057e-04_jprb,2.2058e-04_jprb/)
      forrefo(4,:) = (/ &
     &3.1684e-04_jprb,2.4823e-04_jprb,2.4890e-04_jprb,2.4577e-04_jprb,2.4106e-04_jprb,2.4353e-04_jprb, &
     &2.4038e-04_jprb,2.3932e-04_jprb,2.3604e-04_jprb,2.3773e-04_jprb,2.4243e-04_jprb,2.2597e-04_jprb, &
     &2.2879e-04_jprb,2.2440e-04_jprb,2.1104e-04_jprb,2.1460e-04_jprb/)

!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.62922e-01_jprb, 2.29106e-01_jprb, 1.99640e-01_jprb, 1.73964e-01_jprb, 1.51589e-01_jprb, &
     & 1.32093e-01_jprb, 1.15104e-01_jprb, 1.00300e-01_jprb, 8.74000e-02_jprb, 7.61592e-02_jprb/)
      selfrefo(:, 2) = (/ &
     & 2.45448e-01_jprb, 2.13212e-01_jprb, 1.85210e-01_jprb, 1.60886e-01_jprb, 1.39756e-01_jprb, &
     & 1.21401e-01_jprb, 1.05457e-01_jprb, 9.16070e-02_jprb, 7.95759e-02_jprb, 6.91249e-02_jprb/)
      selfrefo(:, 3) = (/ &
     & 2.41595e-01_jprb, 2.09697e-01_jprb, 1.82010e-01_jprb, 1.57979e-01_jprb, 1.37121e-01_jprb, &
     & 1.19016e-01_jprb, 1.03302e-01_jprb, 8.96630e-02_jprb, 7.78246e-02_jprb, 6.75492e-02_jprb/)
      selfrefo(:, 4) = (/ &
     & 2.44818e-01_jprb, 2.12172e-01_jprb, 1.83880e-01_jprb, 1.59360e-01_jprb, 1.38110e-01_jprb, &
     & 1.19694e-01_jprb, 1.03733e-01_jprb, 8.99010e-02_jprb, 7.79131e-02_jprb, 6.75238e-02_jprb/)
      selfrefo(:, 5) = (/ &
     & 2.43458e-01_jprb, 2.10983e-01_jprb, 1.82840e-01_jprb, 1.58451e-01_jprb, 1.37315e-01_jprb, &
     & 1.18998e-01_jprb, 1.03125e-01_jprb, 8.93690e-02_jprb, 7.74480e-02_jprb, 6.71171e-02_jprb/)
      selfrefo(:, 6) = (/ &
     & 2.40186e-01_jprb, 2.08745e-01_jprb, 1.81420e-01_jprb, 1.57672e-01_jprb, 1.37032e-01_jprb, &
     & 1.19095e-01_jprb, 1.03505e-01_jprb, 8.99560e-02_jprb, 7.81806e-02_jprb, 6.79467e-02_jprb/)
      selfrefo(:, 7) = (/ &
     & 2.42752e-01_jprb, 2.10579e-01_jprb, 1.82670e-01_jprb, 1.58460e-01_jprb, 1.37459e-01_jprb, &
     & 1.19240e-01_jprb, 1.03437e-01_jprb, 8.97280e-02_jprb, 7.78359e-02_jprb, 6.75200e-02_jprb/)
      selfrefo(:, 8) = (/ &
     & 2.39620e-01_jprb, 2.08166e-01_jprb, 1.80840e-01_jprb, 1.57101e-01_jprb, 1.36479e-01_jprb, &
     & 1.18563e-01_jprb, 1.03000e-01_jprb, 8.94790e-02_jprb, 7.77332e-02_jprb, 6.75292e-02_jprb/)
      selfrefo(:, 9) = (/ &
     & 2.38856e-01_jprb, 2.07166e-01_jprb, 1.79680e-01_jprb, 1.55841e-01_jprb, 1.35165e-01_jprb, &
     & 1.17232e-01_jprb, 1.01678e-01_jprb, 8.81880e-02_jprb, 7.64877e-02_jprb, 6.63397e-02_jprb/)
      selfrefo(:,10) = (/ &
     & 2.29821e-01_jprb, 2.00586e-01_jprb, 1.75070e-01_jprb, 1.52800e-01_jprb, 1.33363e-01_jprb, &
     & 1.16398e-01_jprb, 1.01591e-01_jprb, 8.86680e-02_jprb, 7.73887e-02_jprb, 6.75443e-02_jprb/)
      selfrefo(:,11) = (/ &
     & 2.39945e-01_jprb, 2.08186e-01_jprb, 1.80630e-01_jprb, 1.56722e-01_jprb, 1.35978e-01_jprb, &
     & 1.17980e-01_jprb, 1.02364e-01_jprb, 8.88150e-02_jprb, 7.70594e-02_jprb, 6.68598e-02_jprb/)
      selfrefo(:,12) = (/ &
     & 2.40271e-01_jprb, 2.08465e-01_jprb, 1.80870e-01_jprb, 1.56927e-01_jprb, 1.36154e-01_jprb, &
     & 1.18131e-01_jprb, 1.02494e-01_jprb, 8.89260e-02_jprb, 7.71545e-02_jprb, 6.69412e-02_jprb/)
      selfrefo(:,13) = (/ &
     & 2.40503e-01_jprb, 2.08670e-01_jprb, 1.81050e-01_jprb, 1.57086e-01_jprb, 1.36294e-01_jprb, &
     & 1.18254e-01_jprb, 1.02602e-01_jprb, 8.90210e-02_jprb, 7.72380e-02_jprb, 6.70147e-02_jprb/)
      selfrefo(:,14) = (/ &
     & 2.40670e-01_jprb, 2.08811e-01_jprb, 1.81170e-01_jprb, 1.57188e-01_jprb, 1.36380e-01_jprb, &
     & 1.18327e-01_jprb, 1.02663e-01_jprb, 8.90730e-02_jprb, 7.72819e-02_jprb, 6.70517e-02_jprb/)
      selfrefo(:,15) = (/ &
     & 2.40711e-01_jprb, 2.08846e-01_jprb, 1.81200e-01_jprb, 1.57213e-01_jprb, 1.36402e-01_jprb, &
     & 1.18346e-01_jprb, 1.02679e-01_jprb, 8.90870e-02_jprb, 7.72939e-02_jprb, 6.70621e-02_jprb/)
      selfrefo(:,16) = (/ &
     & 2.40727e-01_jprb, 2.08859e-01_jprb, 1.81210e-01_jprb, 1.57221e-01_jprb, 1.36408e-01_jprb, &
     & 1.18350e-01_jprb, 1.02682e-01_jprb, 8.90890e-02_jprb, 7.72952e-02_jprb, 6.70627e-02_jprb/)

if (lhook) call dr_hook('rrtm_kgb4',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb4:error reading file radrrtm")

end subroutine rrtm_kgb4
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

